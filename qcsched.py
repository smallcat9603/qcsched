import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches

HPC_NODES = 10
SCHED_MAP_TIME = 20
NUM_HPC = 3
NUM_QC = 2

class Job:
    def __init__(self, src, id, type, nnodes, time, start, priority):
        self.src = src
        self.id = id
        self.type = type
        self.nnodes = nnodes
        self.time = time 
        self.start = start
        self.priority = priority
        self.rstart = 0
        self.end = 0
        self.status = 'ACCEPT'
        self.map = None
        
    def __repr__(self):
        return f'Job(id={self.id}, type={self.type}, nnodes={self.nnodes}, time={self.time}, start={self.start}, priority={self.priority}, status={self.status}, map={self.map}, src={self.src})'
    
class Semaphore:
    def __init__(self): # 0: available, 1: occupied
        self.qc_flag = {}
        for i in range(NUM_QC):
            self.qc_flag[i] = []
    def __repr__(self):
        return f'qc_flag={self.qc_flag}'  

class JobManager:  
    def __init__(self):
        self.hpc_count = 0
        self.qc_count = {}
        for i in range(NUM_QC):
            self.qc_count[i] = 0
        self.jobs_submitted = []
        self.jobs_scheduled = []
        self.sched_map = np.zeros((HPC_NODES, SCHED_MAP_TIME))

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

def init():
    if 'semaphore' not in st.session_state:
        # co-scheduler
        st.session_state['semaphore'] = Semaphore() 
        for src in range(NUM_HPC): 
            # per hpc
            st.session_state[f'job_manager_{src}'] = JobManager() 

def get_id(src, type):
    if type.startswith('QC'):
        qc = int(type[-1]) - 1
        # record qc job count at src hpc
        st.session_state[f"job_manager_{src}"].qc_count[qc] += 1
        # allocate qc job id, ex 191000
        id = f'{src+1}9{type[-1]}{st.session_state[f"job_manager_{src}"].qc_count[qc]:03d}'            
    else:
        # record hpc job count at src hpc
        st.session_state[f'job_manager_{src}'].hpc_count += 1 
        # allocate hpc job id, ex 100000
        id = f'{src+1}{st.session_state[f"job_manager_{src}"].hpc_count:05d}' 
    return id

def submit(hpc, type, nnodes, time, start, priority):
    src = int(hpc[-1]) - 1  
    id = get_id(src, type)
    # record all jobs at src hpc
    st.session_state[f'job_manager_{src}'].jobs_submitted.append(Job(src=src, id=id, type=type, nnodes=nnodes, time=time, start=start, priority=priority))

@st.cache_data
def import_file(uploaded_file):
    df = pd.read_csv(uploaded_file, comment='#', header=None)
    for idx, row in df.iterrows():
        src = int(row[0][-1]) - 1
        type = row[1]
        nnodes = int(row[2])
        time = int(row[3])
        start = int(row[4])
        priority = int(row[5]) 
        id = get_id(src, type)
        # record all jobs at src hpc
        st.session_state[f'job_manager_{src}'].jobs_submitted.append(Job(src=src, id=id, type=type, nnodes=nnodes, time=time, start=start, priority=priority))

def do_mapping(job, col, row):
    # fill in sched map with job id
    st.session_state[f'job_manager_{job.src}'].sched_map[row:(row+job.nnodes), col:(col+job.time)] = int(job.id)
    # change semaphore status for qc
    if job.type.startswith('QC'):
        st.session_state['semaphore'].qc_flag[int(job.type[-1])-1] += list(range(col, col+job.time))
    # change job status   
    if col == 0:
        job.status = 'RUNNING'
    else:
        if job.status == 'STOP':
            job.status = 'REQUEUED'
        else:
            job.status = 'QUEUED'
    # record map location
    job.map = (col, row)

def release_nodes(job):
    col = job.map[0]
    row = job.map[1]
    # restore sched map with 0
    st.session_state[f'job_manager_{job.src}'].sched_map[row:(row+job.nnodes), col:(col+job.time)] = 0

def stop_jobs(src, ids):
    for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
        if int(job.id) in ids:
            # change job status
            job.status = 'STOP'
            # restore sched map with 0
            release_nodes(job)    

def map(job, algo):
    if job.status == 'ACCEPT' or job.status == 'STOP' or job.status == 'HOLD':
        if job.start + job.time < SCHED_MAP_TIME + 1:
            for col in range(job.start, SCHED_MAP_TIME-job.time+1): # x-axis
                if job.type.startswith('QC'):
                    if set(list(range(col, col+job.time))) & set(st.session_state['semaphore'].qc_flag[int(job.type[-1])-1]):
                        continue                    

                nstop = HPC_NODES-job.nnodes+1
                ids = []
                col_row = None 
                for row in range(HPC_NODES-job.nnodes+1): # y-axis
                    if np.all(st.session_state[f'job_manager_{job.src}'].sched_map[row:(row+job.nnodes), col:(col+job.time)] == 0):
                        do_mapping(job, col, row)
                        return True
                    elif job.type.startswith('QC') and algo == "QPriority":
                        include_qc = False
                        m = st.session_state[f'job_manager_{job.src}'].sched_map[row:(row+job.nnodes), col:(col+job.time)]
                        unique_occupied = list(np.unique(m[m>0]))
                        for id in unique_occupied:
                            if (id//1000)%100 > 90:
                                include_qc = True
                                break
                        if not include_qc and len(unique_occupied) < nstop:
                            # number of hpc jobs to be stopped
                            nstop = len(unique_occupied)
                            # hpc job ids to be stopped
                            ids = unique_occupied
                            # qc job location to be mapped
                            col_row = (col, row)

                if job.type.startswith('QC') and algo == 'QPriority':
                    # stop running hpc jobs
                    stop_jobs(job.src, ids)    
                    # map qc job            
                    do_mapping(job, col_row[0], col_row[1])
                    return True
        # cannot be queued within scheduling period
        job.status = 'HOLD'
    return False

def color(str):
    colors = {
        'hpc': 'green',
        'qc1': 'orange',
        'qc2': 'violet',
    }
    return colors[str]

def schedule(algo):
    for src in range(NUM_HPC):
        st.session_state[f'job_manager_{src}'].jobs_scheduled = st.session_state[f'job_manager_{src}'].jobs_submitted.copy() # separate individuals, but the same child attribute, ex jobs can be separately added or removed but job attribute modification will reflect to both
        if algo == 'FCFS':
            st.session_state[f'job_manager_{src}'].jobs_scheduled.sort(key=sort_key_fcfs)
        elif algo == 'SJF':
            st.session_state[f'job_manager_{src}'].jobs_scheduled.sort(key=sort_key_sjf)
        elif algo.endswith('Priority'):
            st.session_state[f'job_manager_{src}'].jobs_scheduled.sort(key=sort_key_priority)
        for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
            map(job, algo)
        # st.write(sched_map[::-1]) # np.array index align with axis

def show_submitted_jobs():
    for src in range(NUM_HPC):
        st.write(f'HPC{src+1}')
        df = pd.DataFrame([{
            'Job ID': job.id,
            'Status': job.status,
            'Type': job.type,
            'Nodes': job.nnodes,
            'Elapsed Time': job.time,
            'Start Time': job.start,
            'Priority': job.priority,
            } for job in st.session_state[f'job_manager_{src}'].jobs_submitted]
        )
        st.table(df)

def sort_key_fcfs(job):
    type_order = 0 if job.type.startswith('QC') else 1
    return type_order

def sort_key_sjf(job):
    type_order = 0 if job.type.startswith('QC') else 1
    return (type_order, job.nnodes*job.time)

def sort_key_priority(job):
    type_order = 0 if job.type.startswith('QC') else 1
    return (type_order, job.priority)

def plot():
    for src in range(NUM_HPC):
        if len(st.session_state[f'job_manager_{src}'].jobs_scheduled) > 0:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.set_title(f'HPC{src+1}')
            ax.set_ylim(0, HPC_NODES)
            ax.set_xlim(0, SCHED_MAP_TIME)
            ax.set_xlabel('Time')
            ax.set_ylabel('Nodes')
            ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')

            for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
                # st.write(job)
                if job.status == 'RUNNING' or job.status == 'QUEUED' or job.status == 'REQUEUED':
                    fc = color('hpc')
                    if job.type.startswith('QC'):
                        fc = color(f'qc{job.type[-1]}')
                
                    rect = patches.Rectangle(job.map, job.time, job.nnodes, edgecolor='black', facecolor=fc)
                    ax.add_patch(rect)
                    ax.text(job.map[0]+job.time/2, job.map[1]+job.nnodes/2, f'{job.id}-{job.type}', size=10, horizontalalignment='center', verticalalignment='center')

            st.pyplot(fig) 

def update_mapping(nsteps):
    # update sched map
    for src in range(NUM_HPC):
        # move to left
        st.session_state[f'job_manager_{src}'].sched_map[:, :-nsteps] = st.session_state[f'job_manager_{src}'].sched_map[:, nsteps:]
        # fill rightmost with 0
        st.session_state[f'job_manager_{src}'].sched_map[:, -nsteps:] = 0

        # update map location
        for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
            if job.status == 'RUNNING':
                job.time -= nsteps
                if job.time <= 0:
                    job.status = 'FINISH'
            elif job.status == 'QUEUED' or job.status == 'REQUEUED':
                start = job.map[0] - nsteps
                end = start + job.time
                if end <= 0: # start < 0
                    job.status = 'FINISH'
                elif start <= 0: # end > 0
                    job.status = 'RUNNING'        
                    job.map = (0, job.map[1])
                    job.time = end
                else: # start > 0, end > 0
                    job.map = (job.map[0]-nsteps, job.map[1])   
    # update semaphore status for qc
    for qc in range(NUM_QC):
        st.session_state['semaphore'].qc_flag[qc] = [t - nsteps for t in st.session_state['semaphore'].qc_flag[qc] if t - nsteps >= 0] 

def app_layout():

    st.title('QCSched Simulation')

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

    init()

    if st.button(label='Submit'):
        submit(hpc, type, nnodes, time, start, priority)

    expander = st.expander('Import')
    uploaded_file = expander.file_uploader("Choose a file")
    if uploaded_file is not None:
        import_file(uploaded_file)
        
    show_submitted_jobs()

    if st.button(label='Clear'):
        for key in st.session_state.keys():
            del st.session_state[key] 
        st.rerun()

    st.header('Jobs Scheduled')

    algo = st.segmented_control(
        'Select Scheduling Algorithm', ['FCFS', 'SJF', 'Priority', 'QPriority'], default='QPriority'
    )

    if st.button(label='Schedule'):
        schedule(algo)
        st.rerun()

    expander = st.expander('Time Flies')
    nsteps = expander.number_input('Number of Steps Forward', min_value=1, max_value=10, value=1, step=1)
    if expander.button(label='Proceed'):
        update_mapping(nsteps)
        st.rerun()

    plot()
      
if __name__=='__main__':
    app_layout()