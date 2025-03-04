import streamlit as st
from scheduler import schedule
from operation import update_mapping
import time


def simulate(algo: str, resched: str):
    njobs = sum(len(st.session_state[f'job_manager_{src}'].jobs_submitted) for src in range(st.session_state['NUM_HPC']))

    progress_text = "Operation in progress. Please wait."
    my_bar = st.progress(0.0, text=progress_text)

    start_time = time.time()

    while True:
        njobs_fin = sum(1 for src in range(st.session_state['NUM_HPC'])
                        for job in st.session_state[f'job_manager_{src}'].jobs_submitted
                        if job.status == 'FINISH')

        if njobs_fin < njobs:
            update_mapping(1)
            schedule(algo, resched)           
            my_bar.progress(njobs_fin/njobs, text=progress_text)
        else:
            my_bar.empty()
            break       

    end_time = time.time()

    st.info(f'Simulation time: {round(end_time-start_time, 1)}s')
