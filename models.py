import streamlit as st
import numpy as np
import time

from utils import rtime

class Job:
    def __init__(self, src: int, id: str, vid: str, type: str, nnodes: int, elapsed: int, start: int, priority: int):
        self.src = src # 0 (HPC1), 1 (HPC2), 2 (HPC3)
        self.id = id # for node mapping use, can identify qc job or not
        self.vid = vid # for display, cannot identify qc job or not
        self.type = type # HPC, QC1, QC2
        self.nnodes = nnodes # number of HPC nodes
        self.elapsed = elapsed # expected elasped time
        self.start = start # starting time
        self.priority = priority # 1 for highest priority
        self.relapsed = rtime(elapsed) # real elapsed time
        self.wait = 0 # add 1 until running starts
        self.status = 'ACCEPT'
        self.map = None # (col, row)
        self.timestamp = time.time()
        
    def __repr__(self):
        return f'Job(src={self.src}, id={self.id}, vid={self.vid}, type={self.type}, nnodes={self.nnodes}, elapsed={self.elapsed}, start={self.start}, priority={self.priority}, relapsed={self.relapsed}, wait={self.wait}, status={self.status}, map={self.map}, timestamp={self.timestamp})'
    

class Semaphore:
    def __init__(self): 
        self.qc_flag = {i: np.zeros(st.session_state['SCHED_MAP_TIME'], dtype=int) for i in range(st.session_state['NUM_QC'])} # t represents (t, t+1) occupation, 0: available, 1: occupied

    def __repr__(self):
        return f'qc_flag={self.qc_flag}'  


class JobManager:  
    def __init__(self):
        self.hpc_cnt = 0 # add 1 if submitting an hpc job
        self.qc_cnt = {i: 0 for i in range(st.session_state['NUM_QC'])} # add 1 if submitting a qc job
        self.jobs_submitted = [] # all submitted jobs
        self.jobs_scheduled = [] # all scheduled jobs
        self.sched_map = np.zeros((st.session_state['HPC_NODES'], st.session_state['SCHED_MAP_TIME']), dtype=int) # (n, t) represents (t, t+1) occupation at n-th node

    def __repr__(self):
        return self.sched_map[::-1] # np.array index align with axis
