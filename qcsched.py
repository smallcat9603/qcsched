import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches

HPC_NODES = 10
SCHED_MAP_TIME = 20
NUM_HPC = 3

class Job:
    def __init__(self, id, type, nnodes, time, start, priority):
        self.id = id
        self.type = type
        self.nnodes = nnodes
        self.time = time 
        self.start = start
        self.priority = priority
        self.rstart = 0
        self.end = 0
        self.status = 'ACCEPT'
    def __repr__(self):
        return f"Job(id={self.id}, type={self.type}, nnodes={self.nnodes}, time={self.time}, start={self.start}, priority={self.priority})"
    
class Semaphore:
    def __init__(self, qc1_flag, qc2_flag): # 0: available, 1: occupied
        self.qc1_flag = qc1_flag
        self.qc2_flag = qc2_flag
    def __repr__(self):
        return f"qc1_flag={self.qc1_flag}, qc2_flag={self.qc2_flag}"    

# def turn_around_time(arrival,finish):
#     """
#     Provide arrival time and finish time as array inputs returns 
#     """
#     turn_around_time=[]
#     for x in range(len(arrival)):
#         turn_around_time.append(finish[x]-arrival[x])
#     return turn_around_time

# def wait_time(turn_around,burst):
#     """
#     provide array input and will produce array as output
#     """
#     wait=[]
#     for x in range(len(burst)):
#         wait.append(turn_around[x]-burst[x])
#     return wait

# def make_dataframe(process,start,burst,finish,turn_around,wait):
#     """
#     provide the arrays for all attributes and converted dataframe will be returned
#     """
#     df=pd.DataFrame({"Process":process,'Arrival Time':start,"Burst Time":burst,'Completion Time':finish,'Turn Around Time':turn_around,'Waiting Time':wait}).sort_values(by ='Process')
#     return df

def map(job):
    if job.status == "ACCEPT":
        if job.start + job.time < SCHED_MAP_TIME+1:
            for col in range(job.start, SCHED_MAP_TIME-job.time+1): # x-axis
                if job.type == 'QC1':
                    if set(list(range(col, col+job.time))) & set(st.session_state["semaphore"].qc1_flag):
                        continue
                elif job.type == 'QC2':
                    if set(list(range(col, col+job.time))) & set(st.session_state["semaphore"].qc2_flag):
                        continue                    
                for row in range(HPC_NODES-job.nnodes+1): # y-axis
                    idx = job.id[0]
                    if np.all(st.session_state[f"sched_map_{idx}"][row:(row+job.nnodes), col:(col+job.time)] == 0):
                        if job.type == 'QC1':
                            st.session_state[f"sched_map_{idx}"][row:(row+job.nnodes), col:(col+job.time)] = 91
                            st.session_state["semaphore"].qc1_flag += list(range(col, col+job.time))
                        elif job.type == 'QC2':
                            st.session_state[f"sched_map_{idx}"][row:(row+job.nnodes), col:(col+job.time)] = 92
                            st.session_state["semaphore"].qc2_flag += list(range(col, col+job.time))
                        else:
                            st.session_state[f"sched_map_{idx}"][row:(row+job.nnodes), col:(col+job.time)] = 1
                        job.status = "QUEUED"
                        return (col, row)
    if job.status != "QUEUED": 
        job.status = "HOLD"
    return None


def schedule(algo, i):
    st.session_state[f"jobs_scheduled_{i}"] = st.session_state[f"jobs_submitted_{i}"].copy() # separate individuals, but the same child attribute
    if algo == 'FCFS':
        st.session_state[f"jobs_scheduled_{i}"].sort(key=sort_key_fcfs)
    elif algo == 'SJF':
        st.session_state[f"jobs_scheduled_{i}"].sort(key=sort_key_sjf)
    elif algo == 'Priority':
        st.session_state[f"jobs_scheduled_{i}"].sort(key=sort_key_priority)
    for job in st.session_state[f"jobs_scheduled_{i}"]:
        ret = map(job)
        if ret:
            col = ret[0]
            row = ret[1]
            fc = 'green'
            if job.type == 'QC1':
                fc = 'orange'
            elif job.type == 'QC2':
                fc = 'violet'
            rect = patches.Rectangle((col,row), job.time, job.nnodes, edgecolor='black', facecolor=fc)
            st.session_state["axes"][i].add_patch(rect)
            st.session_state["axes"][i].text(col+job.time/2, row+job.nnodes/2, f"{job.id}-{job.type}", size=10, horizontalalignment='center', verticalalignment='center')
    # st.write(sched_map[::-1]) # np.array index align with axis


def show_jobs(str):
    df = pd.DataFrame([{
        'Job ID': job.id,
        'Status': job.status,
        'Type': job.type,
        'Nodes': job.nnodes,
        'Elapsed Time': job.time,
        'Start Time': job.start,
        'Priority': job.priority,
        } for job in st.session_state[str]]
    )
    st.table(df)

def sort_key_fcfs(job):
    type_order = 0 if job.type == 'QC' else 1
    return type_order

def sort_key_sjf(job):
    type_order = 0 if job.type == 'QC' else 1
    return (type_order, job.nnodes*job.time)

def sort_key_priority(job):
    type_order = 0 if job.type == 'QC' else 1
    return (type_order, job.priority)

def app_layout():
    st.title("QCSchduler Simulation")
    st.divider()

    st.header('Jobs Submitted')
    col1, col2 = st.columns([1,1])
    hpc = col1.segmented_control(
        'HPC Selection', ['HPC1', 'HPC2', 'HPC3'], default='HPC1'
    )
    type = col2.segmented_control(
        'Job Type', ['HPC Only', 'QC1', 'QC2'], default='HPC Only'
    )
    col1, col2, col3, col4 = st.columns([1,1,1,1])
    nnodes = col1.number_input('HPC Nodes', min_value=1, max_value=96, value=1, step=1)
    time = col2.number_input('Elapsed Time', min_value=1, max_value=60, value=5, step=1)
    start = col3.number_input('Start Time', min_value=0, max_value=120, value=0, step=1)
    priority = col4.number_input('Priority', min_value=1, max_value=20, value=1, step=1)

    if "semaphore" not in st.session_state:
        st.session_state["semaphore"] = Semaphore([], [])
        st.session_state["fig"], st.session_state["axes"] = plt.subplots(NUM_HPC,1,figsize=(8, 6*NUM_HPC))
        for i in range(NUM_HPC): 
            st.session_state[f"jobs_submitted_{i}"] = []
            st.session_state[f"jobs_scheduled_{i}"] = []
            st.session_state[f"sched_map_{i}"] = np.zeros((HPC_NODES, SCHED_MAP_TIME))
            st.session_state["axes"][i].set_title(f"HPC{i+1}")
            st.session_state["axes"][i].set_ylim(0, HPC_NODES)
            st.session_state["axes"][i].set_xlim(0, SCHED_MAP_TIME)
            st.session_state["axes"][i].set_xlabel('Time')
            st.session_state["axes"][i].set_ylabel('Nodes')
            st.session_state["axes"][i].xaxis.set_major_locator(ticker.MultipleLocator(2))
            st.session_state["axes"][i].xaxis.set_minor_locator(ticker.MultipleLocator(1))
            st.session_state["axes"][i].spines['right'].set_color('none')
            st.session_state["axes"][i].spines['top'].set_color('none')

    if st.button(label='Submit'):
        str = f"jobs_submitted_{int(hpc[-1])-1}"
        st.session_state[str].append(Job(id=f'{int(hpc[-1])-1}{len(st.session_state[str]):05d}', 
                                        type=type, 
                                        nnodes=nnodes,
                                        time=time,
                                        start=start,
                                        priority= priority
                                        )
                                    )
    for i in range(NUM_HPC):
        st.write(f"HPC{i+1}")
        show_jobs(f'jobs_submitted_{i}')
    if st.button(label='Clear'):
        for key in st.session_state.keys():
            del st.session_state[key] 
        st.rerun()

    st.header('Jobs Scheduled')
    algo = st.segmented_control(
        'Select Scheduling Algorithm', ['FCFS', 'SJF', 'Priority'], default='Priority'
    )
    if st.button(label='Schedule'):
        for i in range(NUM_HPC):
            schedule(algo, i)
        st.rerun()

    st.pyplot(st.session_state["fig"]) 
      
if __name__=='__main__':
    app_layout()