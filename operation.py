import streamlit as st
import pandas as pd

from models import Job, Semaphore, JobManager
from utils import get_num_from_0, get_id, get_vid


def init(mode):
    if 'semaphore' not in st.session_state:
        st.session_state['HPC_NODES'] = 10 if mode == 'Demo' else 256
        st.session_state['SCHED_MAP_TIME'] = 20 if mode == 'Demo' else 48
        st.session_state['NUM_HPC'] = 3
        st.session_state['NUM_QC'] = 2
        st.session_state['SUSPEND_TOLERANCE'] = 1

        # co-scheduler
        st.session_state['semaphore'] = Semaphore() 
        for src in range(st.session_state['NUM_HPC']): 
            # per hpc
            st.session_state[f'job_manager_{src}'] = JobManager() 


def append_job(src: int, type: str, nnodes: int, elapsed: int, start: int, priority: int):
    id = get_id(src, type)
    vid = get_vid(src)

    # record all jobs at src hpc
    st.session_state[f'job_manager_{src}'].jobs_submitted.append(Job(src=src, 
                                                                        id=id, 
                                                                        vid=vid, 
                                                                        type=type, 
                                                                        nnodes=nnodes, 
                                                                        elapsed=elapsed, 
                                                                        start=start, 
                                                                        priority=priority))


@st.cache_data
def import_file(uploaded_file):
    df = pd.read_csv(uploaded_file, comment='#', header=None)
    for idx, row in df.iterrows():
        src = get_num_from_0(row[0])
        type = row[1]
        nnodes = int(row[2])
        elapsed = int(row[3])
        start = int(row[4])
        priority = int(row[5]) 

        append_job(src=src, 
                type=type, 
                nnodes=nnodes, 
                elapsed=elapsed,
                start=start,
                priority=priority)


def submit(hpc: str, type: str, nnodes: int, elapsed: int, start: int, priority: int):
    src = get_num_from_0(hpc) 

    append_job(src=src, 
               type=type, 
               nnodes=nnodes, 
               elapsed=elapsed,
               start=start,
               priority=priority)
    

def suspend(job: Job):
    # change job status
    job.status = 'SUSPEND'
    # restore sched map with 0
    release_nodes(job) 

    if job.type.startswith('QC'):
        release_semaphor(job)


def suspend_hpc_jobs(src: int, ids: list):
    for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
        if int(job.id) in ids:
            suspend(job) 


def release_nodes(job: Job):
    col = job.map[0]
    row = job.map[1]
    # restore sched map with 0
    st.session_state[f'job_manager_{job.src}'].sched_map[row:(row+job.nnodes), col:(col+job.elapsed)] = 0


def release_semaphor(job: Job):
    qc = get_num_from_0(job.type)
    st.session_state['semaphore'].qc_flag[qc][job.map[0]:job.map[0]+job.elapsed] = 0


def del_job(vid: str) -> bool:
    if not vid.isdigit():
        return False
    src = get_num_from_0(vid[0])
    for job in st.session_state[f'job_manager_{src}'].jobs_submitted:
        if job.vid == vid:
            if job.status == 'ACCEPT' or job.status == 'HOLD' or job.status == 'QUEUED':
                if job.status == 'QUEUED':
                    release_nodes(job)
                    if job.type.startswith('QC'):
                        release_semaphor(job)
                job.status = 'SUSPEND' # reflected in st.session_state[f'job_manager_{src}'].jobs_scheduled
                st.session_state[f'job_manager_{src}'].jobs_submitted.remove(job) # retained in st.session_state[f'job_manager_{src}'].jobs_scheduled
                return True
            break
    return False 


def resched_jobs(src: int, start: int):
    for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
        if job.status == 'QUEUED' and job.map[0] > start:
            suspend(job)  


def resched_start() -> int:
    cols = []
    for src in range(st.session_state['NUM_HPC']):
        if not st.session_state[f'job_manager_{src}'].jobs_scheduled:
            continue
        
        scheduled_vids = {job.vid for job in st.session_state[f'job_manager_{src}'].jobs_scheduled}
        submitted_vids = {job.vid for job in st.session_state[f'job_manager_{src}'].jobs_submitted}

        added = submitted_vids - scheduled_vids
        if added:
            starts = [job.start for job in added]
            return min(starts)
        
        deleted = scheduled_vids - submitted_vids
        for vid in deleted:
            for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
                if vid == job.vid:
                    if job.status == 'SUSPEND':
                        cols.append(job.map[0])
                    break

    if cols:
        return min(cols)
    else:
        return -1    
    

def finish(job: Job):
    # restore sched map with 0
    release_nodes(job) 

    if job.type.startswith('QC'):
        release_semaphor(job)


def update_mapping_ideal(nsteps: int):
    # update sched map
    for src in range(st.session_state['NUM_HPC']):
        # move to left
        st.session_state[f'job_manager_{src}'].sched_map[:, :-nsteps] = st.session_state[f'job_manager_{src}'].sched_map[:, nsteps:]
        # fill rightmost with 0
        st.session_state[f'job_manager_{src}'].sched_map[:, -nsteps:] = 0

        # update job.map
        for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
            if job.status == 'RUNNING':
                job.elapsed -= nsteps
                if job.elapsed <= 0:
                    job.status = 'FINISH'
            elif job.status == 'QUEUED':
                start = job.map[0] - nsteps
                end = start + job.elapsed
                if end <= 0: # start < 0
                    job.status = 'FINISH'
                elif start <= 0: # end > 0
                    job.status = 'RUNNING'        
                    job.map = (0, job.map[1])
                    job.elapsed = end
                else: # start > 0, end > 0
                    job.map = (job.map[0]-nsteps, job.map[1])   

            # update job.start
            job.start = max(0, job.start-nsteps)

    # update semaphore status for qc
    for i in range(st.session_state['NUM_QC']):
        # move to left
        st.session_state['semaphore'].qc_flag[i][:-nsteps] = st.session_state['semaphore'].qc_flag[i][nsteps:]
        # fill rightmost with 0
        st.session_state['semaphore'].qc_flag[i][-nsteps:] = 0


def update_mapping(nsteps: int):
    # update semaphore status for qc
    for i in range(st.session_state['NUM_QC']):
        # move to left
        st.session_state['semaphore'].qc_flag[i][:-nsteps] = st.session_state['semaphore'].qc_flag[i][nsteps:]
        # fill rightmost with 0
        st.session_state['semaphore'].qc_flag[i][-nsteps:] = 0

    # update sched map
    for src in range(st.session_state['NUM_HPC']):
        # move to left
        st.session_state[f'job_manager_{src}'].sched_map[:, :-nsteps] = st.session_state[f'job_manager_{src}'].sched_map[:, nsteps:]
        # fill rightmost with 0
        st.session_state[f'job_manager_{src}'].sched_map[:, -nsteps:] = 0

        fin_jobs = []
        for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
            # update job.start
            job.start = max(0, job.start-nsteps) 

            # update job.map
            if job.status == 'RUNNING':
                job.elapsed = max(0, job.elapsed-nsteps)
                job.relapsed = max(0, job.relapsed-nsteps)

                if job.relapsed == 0:
                    job.status = 'FINISH'
                    if job.elapsed > 0:
                        fin_jobs.append(job)

            elif job.status == 'QUEUED':
                start = job.map[0] - nsteps
                end = start + job.elapsed
                rend = start + job.relapsed

                job.map = (max(0, start), job.map[1])

                if rend <= 0: # start < 0
                    job.status = 'FINISH'
                    job.elapsed = max(0, end)
                    job.relapsed = 0
                    if end > 0:
                        fin_jobs.append(job)
                elif start <= 0: # rend > 0
                    job.status = 'RUNNING'        
                    job.elapsed = end
                    job.relapsed = rend  

        if fin_jobs:
            for job in fin_jobs:
                finish(job)
            resched_jobs(src, 0)
