import streamlit as st
import plotly.express as px
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

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



def show_df():
    if "jtype" not in st.session_state:
        st.session_state["jtype"] = []
        st.session_state["jnodes"] = []
        st.session_state["jtime"] = []
        st.session_state["jstart"] = []
        st.session_state["jprior"] = []
    df = pd.DataFrame(data = {
        "Job Type": st.session_state["jtype"],
        "Nodes": st.session_state["jnodes"],
        "Elapsed Time": st.session_state["jtime"],
        "Start Time": st.session_state["jstart"],
        "Priority": st.session_state["jprior"],
    })

    st.table(df)

def app_layout():

    # st.set_page_config(layout="wide")
    st.title("QCSchduler Simulation")

    st.divider()

    st.header('Jobs Submitted')

    col1, col2, col3, col4, col5 = st.columns(5)
    jtype = col1.selectbox(
        "Job Type",
        ("QC", "HPC"),
    )
    jnodes = col2.number_input('HPC Nodes', min_value=1, max_value=96, value=1, step=1)
    jtime = col3.number_input('Elapsed Time', min_value=1, max_value=60, value=5, step=1)
    jstart = col4.number_input('Start Time', min_value=-1, max_value=120, value=1, step=1)
    jprior = col5.number_input('Priority', min_value=1, max_value=20, value=1, step=1)

    submit = st.button(label='Submit')
    if submit:
        st.session_state["jtype"].append(jtype)
        st.session_state["jnodes"].append(jnodes)
        st.session_state["jtime"].append(jtime)
        st.session_state["jstart"].append(jstart)
        st.session_state["jprior"].append(jprior)

    show_df()

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

    schedule = st.button(label='Schedule')
    if schedule:
        fig, ax = plt.subplots()
        njobs = len(st.session_state["jtype"])
        for i in range(njobs):
            type = st.session_state["jtype"][i]
            time = st.session_state["jtime"][i]
            nnodes = st.session_state["jnodes"][i]
            start = st.session_state["jstart"][i]
            rect = patches.Rectangle((start,1), time, nnodes, edgecolor="black", facecolor="orange")
            ax.add_patch(rect)
            ax.text(start+time/2, 1+nnodes/2, type, size=50, horizontalalignment='center', verticalalignment='center')
        ax.set_ylim(0, 5)
        ax.set_xlim(0, 10)
        ax.set_xlabel('Time')
        ax.set_ylabel('Nodes')
        # ax.set_yticks([15, 25], labels=['Bill', 'Jim'])     # Modify y-axis tick labels
        # ax.grid(True)                                       # Make grid lines visible
        # ax.set_xticklabels([])
        # ax.set_yticklabels([])
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        st.pyplot(fig) 
      

if __name__=='__main__':
    app_layout()