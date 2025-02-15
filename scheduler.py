import streamlit as st
import numpy as np

from models import Job
from utils import get_num_from_0
from operation import suspend_hpc_jobs, resched_start, resched_jobs


def sort_key_fcfs(job: Job):
    type_order = 0 if job.type.startswith('QC') else 1
    return (type_order, job.timestamp)


def sort_key_sjf(job: Job):
    type_order = 0 if job.type.startswith('QC') else 1
    return (type_order, job.nnodes*job.elapsed, job.timestamp)


def sort_key_priority(job: Job):
    type_order = 0 if job.type.startswith('QC') else 1
    return (type_order, job.priority, job.timestamp)


def sort_key_qc(job: Job):
    qc = get_num_from_0(job.type)
    return (qc, job.priority, job.timestamp)


def do_mapping(job: Job, col: int, row: int):
    # fill in sched map with job id
    st.session_state[f'job_manager_{job.src}'].sched_map[row:(row+job.nnodes), col:(col+job.elapsed)] = int(job.id)
    # change semaphore status for qc
    if job.type.startswith('QC'):
        qc = get_num_from_0(job.type)
        st.session_state['semaphore'].qc_flag[qc][col:col+job.elapsed] = job.vid
    # change job status   
    if col == 0:
        job.status = 'RUNNING'
    else:
        job.status = 'QUEUED'
    # record map location
    job.map = (col, row)


def map(job: Job, algo: str) -> bool:
    if job.status == 'ACCEPT' or job.status == 'SUSPEND' or job.status == 'HOLD':
        if job.start + job.elapsed < st.session_state['SCHED_MAP_TIME'] + 1:
            qc_start = -1
            nsuspend = job.nnodes
            ids = []
            col_row = None 

            for col in range(job.start, st.session_state['SCHED_MAP_TIME']-job.elapsed+1): # x-axis
                if job.type.startswith('QC'):
                    qc = get_num_from_0(job.type)
                    if np.any(st.session_state['semaphore'].qc_flag[qc][col:col+job.elapsed] > 0):
                        continue     
                    if qc_start < 0:
                        qc_start = col # qc job should start at time within (qc_start, qc_start+st.session_state['SUSPEND_TOLERANCE'])

                for row in range(st.session_state['HPC_NODES']-job.nnodes+1): # y-axis
                    if np.all(st.session_state[f'job_manager_{job.src}'].sched_map[row:(row+job.nnodes), col:(col+job.elapsed)] == 0): 
                        do_mapping(job, col, row)
                        return True
                    elif job.type.startswith('QC') and algo == "QPriority":
                        include_qc = False
                        m = st.session_state[f'job_manager_{job.src}'].sched_map[row:(row+job.nnodes), col:(col+job.elapsed)]
                        unique_occupied = list(np.unique(m[m>0]))
                        for id in unique_occupied:
                            if (id//1000)%100 > 90:
                                include_qc = True
                                break
                        if not include_qc and len(unique_occupied) < nsuspend:
                            # number of hpc jobs to be suspended
                            nsuspend = len(unique_occupied)
                            # hpc job ids to be suspended
                            ids = unique_occupied
                            # qc job location to be mapped
                            col_row = (col, row)

                if qc_start >= 0 and col == min(qc_start+st.session_state['SUSPEND_TOLERANCE'], st.session_state['SCHED_MAP_TIME']-job.elapsed) and col_row and job.type.startswith('QC') and algo == 'QPriority': # not suspend hpc job if it will finish in short time (1)
                    # suspend running hpc jobs
                    suspend_hpc_jobs(job.src, ids)    
                    # map qc job            
                    do_mapping(job, col_row[0], col_row[1])

                    return True
                
        # cannot be queued within scheduling period
        job.status = 'HOLD'
    return False


# consider inter-hpc priority
def schedule(algo: str, resched: str):
    qc_sched = [] # all qc jobs

    start = -1 # reschedule starting point, -1 if not needed, 0 if needed for all queued jobs
    if resched == 'Yes':
        start = resched_start()

    for src in range(st.session_state['NUM_HPC']):
        st.session_state[f'job_manager_{src}'].jobs_scheduled = st.session_state[f'job_manager_{src}'].jobs_submitted.copy() # separate individuals, but the same child attribute, ex jobs can be separately added or removed but job attribute modification will reflect to both
        if start >= 0:
            resched_jobs(src, start) # suspend job, release nodes, release semaphor (qc job)

        if algo == 'FCFS':
            st.session_state[f'job_manager_{src}'].jobs_scheduled.sort(key=sort_key_fcfs)
        elif algo == 'SJF':
            st.session_state[f'job_manager_{src}'].jobs_scheduled.sort(key=sort_key_sjf)
        elif algo.endswith('Priority'):
            st.session_state[f'job_manager_{src}'].jobs_scheduled.sort(key=sort_key_priority)
        
        for job in st.session_state[f'job_manager_{src}'].jobs_scheduled: # qc jobs always are in front
            if job.type.startswith('QC'):
                qc_sched.append(job)
            else:
                break

    qc_sched.sort(key=sort_key_qc)

    ############ mapping begins ############

    # first schdule qc jobs
    for qc_job in qc_sched:
        sched = False
        for src in range(st.session_state['NUM_HPC']):
            for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
                if qc_job.id == job.id:
                    map(job, algo)
                    sched = True
                    break
            if sched:
                break

    # then schedule hpc jobs
    for src in range(st.session_state['NUM_HPC']):
        for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
            if not job.type.startswith('QC'):              
                map(job, algo)  