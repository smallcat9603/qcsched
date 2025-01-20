import streamlit as st
import pandas as pd

from constants import NUM_HPC, NUM_QC
from models import Job
from utils import get_num_from_0, get_id, get_vid

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


def update_mapping(nsteps: int):
    # update sched map
    for src in range(NUM_HPC):
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

    # update semaphore status for qc
    for i in range(NUM_QC):
        # move to left
        st.session_state['semaphore'].qc_flag[i][:-nsteps] = st.session_state['semaphore'].qc_flag[i][nsteps:]
        # fill rightmost with 0
        st.session_state['semaphore'].qc_flag[i][-nsteps:] = 0