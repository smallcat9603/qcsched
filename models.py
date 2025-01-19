import numpy as np
import time

from constants import HPC_NODES, SCHED_MAP_TIME, NUM_QC

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
        self.rstart = 0
        self.end = 0
        self.status = 'ACCEPT'
        self.map = None
        self.timestamp = time.time()
        
    def __repr__(self):
        return f'Job(src={self.src}, id={self.id}, vid={self.vid}, type={self.type}, nnodes={self.nnodes}, elapsed={self.elapsed}, start={self.start}, priority={self.priority}, status={self.status}, map={self.map}, timestamp={self.timestamp})'
    

class Semaphore:
    def __init__(self): 
        self.qc_flag = {i: np.zeros(SCHED_MAP_TIME, dtype=int) for i in range(NUM_QC)} # t represents (t, t+1) occupation, 0: available, 1: occupied

    def __repr__(self):
        return f'qc_flag={self.qc_flag}'  


class JobManager:  
    def __init__(self):
        self.hpc_cnt = 0 # add 1 if submitting an hpc job
        self.qc_cnt = {i: 0 for i in range(NUM_QC)} # add 1 if submitting a qc job
        self.jobs_submitted = [] # all submitted jobs
        self.jobs_scheduled = [] # all scheduled jobs
        self.sched_map = np.zeros((HPC_NODES, SCHED_MAP_TIME), dtype=int) # (n, t) represents (t, t+1) occupation at n-th node
