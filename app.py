import streamlit as st

from utils import show_submitted_jobs
from stats import show_statistics, show_results
from visualizer import plot
from operation import init, import_file, submit, update_mapping, del_job
from scheduler import schedule
from simulator import simulate


def run():

    st.set_page_config(layout="wide")

    st.sidebar.title('QCSched Simulation')

    mode = st.sidebar.radio('Mode', ['Demo', 'Simulation'], horizontal=True, label_visibility='collapsed')

    init(mode) # initialize constants, st.session_state['semaphore'] and st.session_state[f'job_manager_{src}']

    if st.sidebar.button(label='Reset'):
        st.cache_data.clear() # clear cache data via @st.cache_data, not including st.session_state
        for key in st.session_state.keys():
            del st.session_state[key] 
        st.rerun()

    ################################## Jobs Submitted ##################################

    if mode == 'Demo':
        hpc = st.sidebar.segmented_control(
            'HPC Selection', ['HPC1', 'HPC2', 'HPC3'], default='HPC1'
        )
        type = st.sidebar.segmented_control(
            'Job Type', ['HPC Only', 'QC1', 'QC2'], default='HPC Only'
        )

        col1, col2, col3, col4 = st.sidebar.columns([1,1,1,1])
        nnodes = col1.number_input('HPC Nodes', min_value=1, max_value=96, value=5, step=1)
        elapsed = col2.number_input('Elapsed Time', min_value=1, max_value=60, value=5, step=1)
        start = col3.number_input('Start Time', min_value=0, max_value=120, value=0, step=1)
        priority = col4.number_input('Job Priority', min_value=1, max_value=20, value=1, step=1)

        if st.sidebar.button(label='Submit', type='primary'):
            submit(hpc, type, nnodes, elapsed, start, priority)

    expander = st.sidebar.expander('Import')
    uploaded_file = expander.file_uploader("Choose a file")
    if uploaded_file is not None:
        import_file(uploaded_file)

    expander = st.sidebar.expander('Delete')
    vid = expander.text_input('Job ID', placeholder='100001')
    if expander.button(label='Delete'):
        if del_job(vid):
            expander.success(f'Job {vid} deleted!')
        else:
            expander.error(f'Job {vid} deletion failed!')
        
    show_submitted_jobs()

    ################################## Jobs Scheduled ##################################

    algo = st.sidebar.segmented_control(
        'Select Scheduling Algorithm', ['FCFS', 'SJF', 'Priority', 'QPriority'], default='QPriority'
    )
    resched = st.sidebar.segmented_control(
        'Allow Re-scheduling (e.g., after HPC Jobs are suspended)', ['No', 'Yes'], default='Yes'
    )

    if st.sidebar.button(label='Schedule', type='primary'):
        schedule(algo, resched)
        st.rerun()   

    plot(mode) # include qc and hpc   

    if mode == 'Demo':

        expander = st.sidebar.expander('Time Flies')
        nsteps = expander.number_input('Number of Steps Forward', min_value=1, max_value=10, value=1, step=1)
        if expander.button(label='Proceed'):
            update_mapping(nsteps)
            schedule(algo, resched)
            st.rerun()

        show_statistics()

    ################################## Simulation ##################################

    else: # mode == 'simulation'      

        if st.sidebar.button(label='Simulate', type='primary'):
            simulate(algo, resched)
            st.rerun()

        show_results()

if __name__=='__main__':
    run()