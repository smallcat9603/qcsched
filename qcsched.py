import streamlit as st
import plotly.express as px
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

HPC_NODES = 10
SCHED_MAP_TIME = 20

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
        self.status = ""
    def __repr__(self):
        return f"Job(id={self.id}, type={self.type}, nnodes={self.nnodes}, time={self.time}, start={self.start}, priority={self.priority})"

def turn_around_time(arrival,finish):
    """
    Provide arrival time and finish time as array inputs returns 
    """
    turn_around_time=[]
    for x in range(len(arrival)):
        turn_around_time.append(finish[x]-arrival[x])
    return turn_around_time

def wait_time(turn_around,burst):
    """
    provide array input and will produce array as output
    """
    wait=[]
    for x in range(len(burst)):
        wait.append(turn_around[x]-burst[x])
    return wait

def make_dataframe(process,start,burst,finish,turn_around,wait):
    """
    provide the arrays for all attributes and converted dataframe will be returned
    """
    df=pd.DataFrame({"Process":process,'Arrival Time':start,"Burst Time":burst,'Completion Time':finish,'Turn Around Time':turn_around,'Waiting Time':wait}).sort_values(by ='Process')
    return df



def app_layout():

    # st.set_page_config(layout="wide")
    st.title("QCSchduler Simulation")

    st.divider()

    st.header('Jobs Submitted')

    col1, col2, col3, col4, col5 = st.columns(5)
    type = col1.selectbox(
        "Job Type",
        ("QC", "HPC"),
    )
    nnodes = col2.number_input('HPC Nodes', min_value=1, max_value=96, value=1, step=1)
    time = col3.number_input('Elapsed Time', min_value=1, max_value=60, value=5, step=1)
    start = col4.number_input('Start Time', min_value=-1, max_value=120, value=1, step=1)
    priority = col5.number_input('Priority', min_value=1, max_value=20, value=1, step=1)

    if "jobs" not in st.session_state:
        st.session_state["jobs"] = []

    submit = st.button(label='Submit')
    if submit:
        st.session_state["jobs"].append(Job(id=f'{len(st.session_state["jobs"]):05d}', 
                                            type=type, 
                                            nnodes=nnodes,
                                            time=time,
                                            start=start,
                                            priority= priority
                                            )
                                        )

    df = pd.DataFrame([{
                    'Job ID': job.id,
                    'Type': job.type,
                    'Nodes': job.nnodes,
                    'Elapsed Time': job.time,
                    'Start Time': job.start,
                    'Priority': job.priority,
                    } for job in st.session_state["jobs"]]
    )

    st.table(df)

    clear = st.button(label='Clear')
    if clear:
        for key in st.session_state.keys():
            del st.session_state[key] 
        st.rerun()

    st.divider()

    st.header('Jobs Scheduled')

    algo = st.selectbox(
        'Select Scheduling Algorithm', 
        ('FCFS', 'SJF', 'Priority'),
    )

    if "sched_map" not in st.session_state:
        st.session_state["sched_map"] = [[0] * HPC_NODES for _ in range(SCHED_MAP_TIME)]

    schedule = st.button(label='Schedule')
    if schedule:
        fig, ax = plt.subplots()
        ax.set_ylim(0, HPC_NODES)
        ax.set_xlim(0, SCHED_MAP_TIME)
        ax.set_xlabel('Time')
        ax.set_ylabel('Nodes')
        # ax.grid(True)
        # ax.set_xticklabels([])
        # ax.set_yticklabels([])
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        if len(st.session_state["jobs"]) > 0:
            for job in st.session_state["jobs"]:
                rect = patches.Rectangle((job.start,0), job.time, job.nnodes, edgecolor="black", facecolor="orange")
                ax.add_patch(rect)
                ax.text(job.start+job.time/2, job.nnodes/2, f"{job.id}-{job.type}", size=20, horizontalalignment='center', verticalalignment='center')

        st.pyplot(fig) 
      

if __name__=='__main__':
    app_layout()