import streamlit as st
from scheduler import schedule
from operation import update_mapping


def simulate(algo: str, resched: str):
    while any(job.status != 'FINISH' for src in range(st.session_state['NUM_HPC']) for job in st.session_state[f'job_manager_{src}'].jobs_submitted):
        update_mapping(1)
        schedule(algo, resched)
