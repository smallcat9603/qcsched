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
SUSPEND_TOLERANCE = 1

class Job:
    def __init__(self, src, id, vid, type, nnodes, time, start, priority):
        self.src = src
        self.id = id # for node mapping use, can identify qc job or not
        self.vid = vid # for display, cannot identify qc job or not
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
        return f'Job(src={self.src}, id={self.id}, vid={self.vid}, type={self.type}, nnodes={self.nnodes}, time={self.time}, start={self.start}, priority={self.priority}, status={self.status}, map={self.map})'
    
class Semaphore:
    def __init__(self): # 0: available, 1: occupied
        self.qc_flag = {}
        for i in range(NUM_QC): 
            self.qc_flag[i] = np.zeros(SCHED_MAP_TIME, dtype=int) # t represents (t, t+1) occupation

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
        self.sched_map = np.zeros((HPC_NODES, SCHED_MAP_TIME), dtype=int) # (n, t) represents (t, t+1) occupation at n-th node

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

def wait_time(scheduled, scheduled_qc):
    wtime = 0.0
    wtime_qc = 0.0

    for src in range(NUM_HPC):
        for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
            if job.status == 'QUEUED':
                wtime += job.map[0]
                if job.type.startswith('QC'):
                    wtime_qc += job.map[0]

    avg_wtime = 'na'
    avg_wtime_qc = 'na'
    if scheduled > 0:
        avg_wtime = str(wtime/scheduled)
    if scheduled_qc > 0:
        avg_wtime_qc = str(wtime_qc/scheduled_qc)

    return avg_wtime, avg_wtime_qc

def show_statistics():
    total = 0
    running = 0
    queued = 0
    hold = 0
    finish = 0

    total_qc = 0
    running_qc = 0
    queued_qc = 0
    hold_qc = 0
    finish_qc = 0

    for src in range(NUM_HPC):
        total += len(st.session_state[f'job_manager_{src}'].jobs_submitted)
        for job in st.session_state[f'job_manager_{src}'].jobs_submitted:
            if job.type.startswith('QC'):
                total_qc += 1
            if job.status == 'RUNNING':
                running += 1
                if job.type.startswith('QC'):
                    running_qc += 1
            elif job.status == 'QUEUED':
                queued += 1
                if job.type.startswith('QC'):
                    queued_qc += 1
            elif job.status == 'HOLD':
                hold += 1
                if job.type.startswith('QC'):
                    hold_qc += 1
            elif job.status == 'FINISH':
                finish += 1
                if job.type.startswith('QC'):
                    finish_qc += 1

    avg_wtime, avg_wtime_qc = wait_time(running+queued, running_qc+queued_qc)

    st.header('Statistics')
    cols = st.columns(6)
    cols[0].metric('Total', f'{total} ({total_qc})')
    cols[1].metric('Running', f'{running} ({running_qc})')
    cols[2].metric('Queued', f'{queued} ({queued_qc})')
    cols[3].metric('Hold', f'{hold} ({hold_qc})')
    cols[4].metric('Finish', f'{finish} ({finish_qc})')
    cols[5].metric('Wait', f'{avg_wtime} ({avg_wtime_qc})')

def init():
    if 'semaphore' not in st.session_state:
        # co-scheduler
        st.session_state['semaphore'] = Semaphore() 
        for src in range(NUM_HPC): 
            # per hpc
            st.session_state[f'job_manager_{src}'] = JobManager() 

def get_num_from_0(str):
    return int(str[-1]) - 1

def get_id(src, type):
    if type.startswith('QC'):
        qc = get_num_from_0(type)
        # record qc job count at src hpc
        st.session_state[f"job_manager_{src}"].qc_count[qc] += 1
        # allocate qc job id, ex 191000
        id = f'{src+1}9{qc+1}{st.session_state[f"job_manager_{src}"].qc_count[qc]:03d}'            
    else:
        # record hpc job count at src hpc
        st.session_state[f'job_manager_{src}'].hpc_count += 1 
        # allocate hpc job id, ex 100000
        id = f'{src+1}{st.session_state[f"job_manager_{src}"].hpc_count:05d}' 
    return id

def get_vid(src):
    vid = st.session_state[f"job_manager_{src}"].hpc_count
    for i in range(NUM_QC):
        vid += st.session_state[f"job_manager_{src}"].qc_count[i]
    vid = f'{src+1}{vid:05d}'            
    return vid

def submit(hpc, type, nnodes, time, start, priority):
    src = get_num_from_0(hpc) 
    id = get_id(src, type)
    vid = get_vid(src)
    # record all jobs at src hpc
    st.session_state[f'job_manager_{src}'].jobs_submitted.append(Job(src=src, id=id, vid=vid, type=type, nnodes=nnodes, time=time, start=start, priority=priority))

@st.cache_data
def import_file(uploaded_file):
    df = pd.read_csv(uploaded_file, comment='#', header=None)
    for idx, row in df.iterrows():
        src = get_num_from_0(row[0])
        type = row[1]
        nnodes = int(row[2])
        time = int(row[3])
        start = int(row[4])
        priority = int(row[5]) 
        id = get_id(src, type)
        vid = get_vid(src)
        # record all jobs at src hpc
        st.session_state[f'job_manager_{src}'].jobs_submitted.append(Job(src=src, id=id, vid=vid, type=type, nnodes=nnodes, time=time, start=start, priority=priority))

def do_mapping(job, col, row):
    # fill in sched map with job id
    st.session_state[f'job_manager_{job.src}'].sched_map[row:(row+job.nnodes), col:(col+job.time)] = int(job.id)
    # change semaphore status for qc
    if job.type.startswith('QC'):
        qc = get_num_from_0(job.type)
        st.session_state['semaphore'].qc_flag[qc][col:col+job.time] = job.vid
    # change job status   
    if col == 0:
        job.status = 'RUNNING'
    else:
        job.status = 'QUEUED'
    # record map location
    job.map = (col, row)

def release_nodes(job):
    col = job.map[0]
    row = job.map[1]
    # restore sched map with 0
    st.session_state[f'job_manager_{job.src}'].sched_map[row:(row+job.nnodes), col:(col+job.time)] = 0

def suspend(job):
    # change job status
    job.status = 'SUSPEND'
    # restore sched map with 0
    release_nodes(job)      

def suspend_jobs(src, ids):
    for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
        if int(job.id) in ids:
            suspend(job) 

def check_mapping(src):
    st.write(st.session_state[f'job_manager_{src}'].sched_map[::-1]) # np.array index align with axis  

def resched_jobs(src, qc_start):
    for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
        if not job.type.startswith('QC') and job.status == 'QUEUED' and job.map[0] > qc_start:
            suspend(job)    

def map(job, algo, resched):
    if job.status == 'ACCEPT' or job.status == 'SUSPEND' or job.status == 'HOLD':
        if job.start + job.time < SCHED_MAP_TIME + 1:
            qc_start = -1
            nsuspend = job.nnodes
            ids = []
            col_row = None 

            for col in range(job.start, SCHED_MAP_TIME-job.time+1): # x-axis
                if job.type.startswith('QC'):
                    qc = get_num_from_0(job.type)
                    if np.any(st.session_state['semaphore'].qc_flag[qc][col:col+job.time] > 0):
                        continue     
                    if qc_start < 0:
                        qc_start = col # qc job should start at time within (qc_start, qc_start+SUSPEND_TOLERANCE)

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
                        if not include_qc and len(unique_occupied) < nsuspend:
                            # number of hpc jobs to be suspended
                            nsuspend = len(unique_occupied)
                            # hpc job ids to be suspended
                            ids = unique_occupied
                            # qc job location to be mapped
                            col_row = (col, row)

                if qc_start >= 0 and col == min(qc_start+SUSPEND_TOLERANCE, SCHED_MAP_TIME-job.time) and col_row and job.type.startswith('QC') and algo == 'QPriority': # not suspend hpc job if it will finish in short time (1)
                    # suspend running hpc jobs
                    suspend_jobs(job.src, ids)    
                    # map qc job            
                    do_mapping(job, col_row[0], col_row[1])

                    # reschedule hpc jobs starting after qc_start
                    if resched == "Yes":
                        resched_jobs(job.src, qc_start)  

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

# consider inter-hpc priority
def schedule(algo, resched):
    qc_scheduled = [] # all qc jobs
    for src in range(NUM_HPC):
        st.session_state[f'job_manager_{src}'].jobs_scheduled = st.session_state[f'job_manager_{src}'].jobs_submitted.copy() # separate individuals, but the same child attribute, ex jobs can be separately added or removed but job attribute modification will reflect to both
        if algo == 'FCFS':
            st.session_state[f'job_manager_{src}'].jobs_scheduled.sort(key=sort_key_fcfs)
        elif algo == 'SJF':
            st.session_state[f'job_manager_{src}'].jobs_scheduled.sort(key=sort_key_sjf)
        elif algo.endswith('Priority'):
            st.session_state[f'job_manager_{src}'].jobs_scheduled.sort(key=sort_key_priority)
        
        for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
            if job.type.startswith('QC'):
                qc_scheduled.append(job)
            else:
                break

    qc_scheduled.sort(key=sort_key_qc)

    # first schdule qc jobs
    for qc_job in qc_scheduled:
        for src in range(NUM_HPC):
            for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
                if qc_job.id == job.id:
                    map(job, algo, resched)

    # then schedule hpc jobs
    for src in range(NUM_HPC):
        for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
            if not job.type.startswith('QC'):              
                map(job, algo, resched)  

def show_submitted_jobs():
    st.header('Jobs Submitted')

    tabs = st.tabs([f'HPC{src+1}' for src in range(NUM_HPC)])

    for src in range(NUM_HPC):
        if len(st.session_state[f'job_manager_{src}'].jobs_submitted) > 0:
            df = pd.DataFrame([{
                'Job ID': job.vid,
                'Status': job.status,
                'Type': job.type,
                'Nodes': job.nnodes,
                'Elapsed Time': job.time,
                'Start Time': job.start,
                'Priority': job.priority,
                } for job in st.session_state[f'job_manager_{src}'].jobs_submitted]
            )
            tabs[src].table(df)
        else:
            tabs[src].info(f'No job submitted in HPC{src+1}')

def sort_key_fcfs(job):
    type_order = 0 if job.type.startswith('QC') else 1
    return type_order

def sort_key_sjf(job):
    type_order = 0 if job.type.startswith('QC') else 1
    return (type_order, job.nnodes*job.time)

def sort_key_priority(job):
    type_order = 0 if job.type.startswith('QC') else 1
    return (type_order, job.priority)

def sort_key_qc(job):
    qc = get_num_from_0(job.type)
    return (qc, job.priority)

def qc_util(arr_qc_flag): 
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

def plot_qc():
    for i in range(NUM_QC):
        # st.write(st.session_state['semaphore'].qc_flag[i])
        if np.any(st.session_state['semaphore'].qc_flag[i] > 0):
            fig, ax = plt.subplots(figsize=(8, 1))
            # ax.set_title(f'QC{i+1}')
            ax.set_ylim(0, 1)
            ax.set_xlim(0, SCHED_MAP_TIME)
            ax.set_xlabel('Time')
            ax.set_ylabel('Utilization')
            ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')

            fc = color(f'qc{i+1}')
            utils = qc_util(st.session_state['semaphore'].qc_flag[i])
            for util in utils:
                rect = patches.Rectangle((util[0], 0), util[1], 1, edgecolor='black', facecolor=fc)
                ax.add_patch(rect)
                ax.text(util[0]+util[1]/2, 1/2, f'{st.session_state['semaphore'].qc_flag[i][util[0]]}', size=10, horizontalalignment='center', verticalalignment='center')   

            st.header(f'QC{i+1} Schedule')
            st.pyplot(fig) 

def plot_hpc():
    for src in range(NUM_HPC):
        if len(st.session_state[f'job_manager_{src}'].jobs_scheduled) > 0:
            fig, ax = plt.subplots(figsize=(8, 4))
            # ax.set_title(f'HPC{src+1}')
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
                if job.status == 'RUNNING' or job.status == 'QUEUED':
                    fc = color('hpc')
                    if job.type.startswith('QC'):
                        qc = get_num_from_0(job.type)
                        fc = color(f'qc{qc+1}')
                
                    rect = patches.Rectangle(job.map, job.time, job.nnodes, edgecolor='black', facecolor=fc)
                    ax.add_patch(rect)
                    ax.text(job.map[0]+job.time/2, job.map[1]+job.nnodes/2, f'{job.vid}-{job.priority}-{job.type}', size=10, horizontalalignment='center', verticalalignment='center')

            st.header(f'HPC{src+1} Schedule')
            st.pyplot(fig)

def plot():
    plot_qc()
    plot_hpc()

def update_mapping(nsteps):
    # update sched map
    for src in range(NUM_HPC):
        # move to left
        st.session_state[f'job_manager_{src}'].sched_map[:, :-nsteps] = st.session_state[f'job_manager_{src}'].sched_map[:, nsteps:]
        # fill rightmost with 0
        st.session_state[f'job_manager_{src}'].sched_map[:, -nsteps:] = 0

        # update job.map
        for job in st.session_state[f'job_manager_{src}'].jobs_scheduled:
            if job.status == 'RUNNING':
                job.time -= nsteps
                if job.time <= 0:
                    job.status = 'FINISH'
            elif job.status == 'QUEUED':
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
    for i in range(NUM_QC):
        # move to left
        st.session_state['semaphore'].qc_flag[i][:-nsteps] = st.session_state['semaphore'].qc_flag[i][nsteps:]
        # fill rightmost with 0
        st.session_state['semaphore'].qc_flag[i][-nsteps:] = 0

def app_layout():

    st.set_page_config(layout="wide")

    st.sidebar.title('QCSched Simulation')

    ################################## Jobs Submitted ##################################

    hpc = st.sidebar.segmented_control(
        'HPC Selection', ['HPC1', 'HPC2', 'HPC3'], default='HPC1'
    )
    type = st.sidebar.segmented_control(
        'Job Type', ['HPC Only', 'QC1', 'QC2'], default='HPC Only'
    )

    col1, col2, col3, col4 = st.sidebar.columns([1,1,1,1])
    nnodes = col1.number_input('HPC Nodes', min_value=1, max_value=96, value=1, step=1)
    time = col2.number_input('Elapsed Time', min_value=1, max_value=60, value=5, step=1)
    start = col3.number_input('Start Time', min_value=0, max_value=120, value=0, step=1)
    priority = col4.number_input('Job Priority', min_value=1, max_value=20, value=1, step=1)

    init()

    col1, col2 = st.sidebar.columns([1,1])
    if col1.button(label='Submit', type='primary'):
        submit(hpc, type, nnodes, time, start, priority)
    if col2.button(label='Clear'):
        st.cache_data.clear() # clear cache data via @st.cache_data, not including st.session_state
        for key in st.session_state.keys():
            del st.session_state[key] 
        st.rerun()

    expander = st.sidebar.expander('Import')
    uploaded_file = expander.file_uploader("Choose a file")
    if uploaded_file is not None:
        import_file(uploaded_file)
        
    show_submitted_jobs()

    ################################## Jobs Scheduled ##################################

    algo = st.sidebar.segmented_control(
        'Select Scheduling Algorithm', ['FCFS', 'SJF', 'Priority', 'QPriority'], default='QPriority'
    )
    resched = st.sidebar.segmented_control(
        'Allow Re-scheduling after HPC Jobs are suspended', ['No', 'Yes'], default='Yes'
    )

    if st.sidebar.button(label='Schedule', type='primary'):
        schedule(algo, resched)
        st.rerun()

    expander = st.sidebar.expander('Time Flies')
    nsteps = expander.number_input('Number of Steps Forward', min_value=1, max_value=10, value=1, step=1)
    if expander.button(label='Proceed'):
        update_mapping(nsteps)
        schedule(algo, resched)
        st.rerun()

    show_statistics()

    plot()

if __name__=='__main__':
    app_layout()