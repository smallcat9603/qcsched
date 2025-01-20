import streamlit as st
import pandas as pd

from constants import SCHED_MAP_TIME, NUM_HPC, NUM_QC
from models import Job, Semaphore, JobManager


def init():
    if 'semaphore' not in st.session_state:
        # co-scheduler
        st.session_state['semaphore'] = Semaphore() 
        for src in range(NUM_HPC): 
            # per hpc
            st.session_state[f'job_manager_{src}'] = JobManager() 


def get_num_from_0(text: str) -> int:
    return int(text[-1]) - 1


def get_id(src: int, type: str) -> str:
    if type.startswith('QC'):
        qc = get_num_from_0(type)
        # record qc job count at src hpc
        st.session_state[f"job_manager_{src}"].qc_cnt[qc] += 1
        # allocate qc job id, ex 191000
        id = f'{src+1}9{qc+1}{st.session_state[f"job_manager_{src}"].qc_cnt[qc]:03d}'            
    else:
        # record hpc job count at src hpc
        st.session_state[f'job_manager_{src}'].hpc_cnt += 1 
        # allocate hpc job id, ex 100000
        id = f'{src+1}{st.session_state[f"job_manager_{src}"].hpc_cnt:05d}' 
    return id


def get_vid(src: int) -> str:
    vid = st.session_state[f"job_manager_{src}"].hpc_cnt
    for i in range(NUM_QC):
        vid += st.session_state[f"job_manager_{src}"].qc_cnt[i]
    vid = f'{src+1}{vid:05d}'            
    return vid


def color(text: str) -> str:
    colors = {
        'hpc': 'green',
        'qc1': 'orange',
        'qc2': 'violet',
    }
    return colors[text]


def show_submitted_jobs():
    st.header('Jobs Submitted')

    tabs = st.tabs([f'HPC{src+1}' for src in range(NUM_HPC)])

    for src in range(NUM_HPC):
        if st.session_state[f'job_manager_{src}'].jobs_submitted:
            df = pd.DataFrame([{
                'Job ID': job.vid,
                'Status': job.status,
                'Type': job.type,
                'Nodes': job.nnodes,
                'Elapsed Time': job.elapsed,
                'Start Time': job.start,
                'Priority': job.priority,
                } for job in st.session_state[f'job_manager_{src}'].jobs_submitted]
            )
            tabs[src].table(df)
        else:
            tabs[src].info(f'No job submitted in HPC{src+1}')


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


def qc_util(arr_qc_flag) -> list: 
    """
    Provide list of (start, len) according to non-zero elements in qc_flag 
    """
    result = []
    cur = 0
    for t in range(1, SCHED_MAP_TIME):
        if arr_qc_flag[t] != arr_qc_flag[t-1]:
            if arr_qc_flag[t-1] > 0:
                result.append((cur, t-cur))
            if arr_qc_flag[t] > 0:
                cur = t
        if t == SCHED_MAP_TIME - 1 and arr_qc_flag[t] > 0:
            result.append((cur, t-cur+1))
    return result
