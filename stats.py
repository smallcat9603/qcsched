import streamlit as st
import numpy as np

from constants import NUM_HPC


def show_statistics():
    total = 0
    running = 0
    queued = 0
    hold = 0
    finish = 0
    wtime = []

    total_qc = 0
    running_qc = 0
    queued_qc = 0
    hold_qc = 0
    finish_qc = 0
    wtime_qc = []

    for src in range(NUM_HPC):
        total += len(st.session_state[f'job_manager_{src}'].jobs_submitted)
        for job in st.session_state[f'job_manager_{src}'].jobs_submitted:
            if job.type.startswith('QC'):
                total_qc += 1
            if job.status == 'RUNNING':
                running += 1
                wtime.append(0)
                if job.type.startswith('QC'):
                    running_qc += 1
                    wtime_qc.append(0)
            elif job.status == 'QUEUED':
                queued += 1
                wtime.append(job.map[0])
                if job.type.startswith('QC'):
                    queued_qc += 1
                    wtime_qc.append(job.map[0])
            elif job.status == 'HOLD':
                hold += 1
                if job.type.startswith('QC'):
                    hold_qc += 1
            elif job.status == 'FINISH':
                finish += 1
                if job.type.startswith('QC'):
                    finish_qc += 1

    avg_wtime = round(np.mean(wtime), 1) if wtime else 0.0
    avg_wtime_qc = round(np.mean(wtime_qc), 1) if wtime_qc else 0.0

    st.header('Statistics (QC)')
    cols = st.columns(6)
    cols[0].metric('Total', f'{total} ({total_qc})')
    cols[1].metric('Running', f'{running} ({running_qc})')
    cols[2].metric('Queued', f'{queued} ({queued_qc})')
    cols[3].metric('Hold', f'{hold} ({hold_qc})')
    cols[4].metric('Finish', f'{finish} ({finish_qc})')
    cols[5].metric('Wait', f'{avg_wtime} ({avg_wtime_qc})')
