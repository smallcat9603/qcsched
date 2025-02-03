import streamlit as st
import numpy as np
import pandas as pd


def show_statistics():
    total = running = queued = hold = finish = 0
    wtime = []

    total_qc = running_qc = queued_qc = hold_qc = finish_qc = 0
    wtime_qc = []

    for src in range(st.session_state['NUM_HPC']):
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


def show_results():
    njobs = njobs_qc = njobs_qc_priority_1 = njobs_qc_priority_2 = njobs_qc_priority_3 = 0
    wtime = []
    wtime_qc = []
    wtime_qc_priority_1 = []
    wtime_qc_priority_2 = []
    wtime_qc_priority_3 = []
    rtime = []
    rtime_qc = []
    rtime_qc_priority_1 = []
    rtime_qc_priority_2 = []
    rtime_qc_priority_3 = []

    for src in range(st.session_state['NUM_HPC']):
        njobs += len(st.session_state[f'job_manager_{src}'].jobs_submitted)
        for job in st.session_state[f'job_manager_{src}'].jobs_submitted:
            wtime.append(job.wait)
            rtime.append(job.wait+job.info[2])
            if job.type.startswith('QC'):
                njobs_qc += 1
                wtime_qc.append(job.wait)
                rtime_qc.append(job.wait+job.info[2])
                if job.priority == 1:
                    njobs_qc_priority_1 += 1
                    wtime_qc_priority_1.append(job.wait)
                    rtime_qc_priority_1.append(job.wait+job.info[2])
                if job.priority <= 2:
                    njobs_qc_priority_2 += 1
                    wtime_qc_priority_2.append(job.wait)
                    rtime_qc_priority_2.append(job.wait+job.info[2])        
                if job.priority <= 3:
                    njobs_qc_priority_3 += 1
                    wtime_qc_priority_3.append(job.wait)
                    rtime_qc_priority_3.append(job.wait+job.info[2])
    
    avg_wtime = round(np.mean(wtime), 1) if wtime else 0.0
    avg_wtime_qc = round(np.mean(wtime_qc), 1) if wtime_qc else 0.0  
    avg_wtime_qc_priority_1 = round(np.mean(wtime_qc_priority_1), 1) if wtime_qc_priority_1 else 0.0
    avg_wtime_qc_priority_2 = round(np.mean(wtime_qc_priority_2), 1) if wtime_qc_priority_2 else 0.0
    avg_wtime_qc_priority_3 = round(np.mean(wtime_qc_priority_3), 1) if wtime_qc_priority_3 else 0.0
    avg_rtime = round(np.mean(rtime), 1) if rtime else 0.0
    avg_rtime_qc = round(np.mean(rtime_qc), 1) if rtime_qc else 0.0
    avg_rtime_qc_priority_1 = round(np.mean(rtime_qc_priority_1), 1) if rtime_qc_priority_1 else 0.0
    avg_rtime_qc_priority_2 = round(np.mean(rtime_qc_priority_2), 1) if rtime_qc_priority_2 else 0.0
    avg_rtime_qc_priority_3 = round(np.mean(rtime_qc_priority_3), 1) if rtime_qc_priority_3 else 0.0

    st.header('Simulation Result')
    df = pd.DataFrame({
        'Num of all jobs': njobs,
        'Num of QC jobs': njobs_qc,
        'Num of HPC jobs': njobs - njobs_qc,
        'Avg wait time of all jobs': avg_wtime,
        'Avg wait time of QC jobs': avg_wtime_qc,
        'Avg wait time of QC (Priority=1) jobs': avg_wtime_qc_priority_1,
        'Avg wait time of QC (Priority<=2) jobs': avg_wtime_qc_priority_2,
        'Avg wait time of QC (Priority<=3) jobs': avg_wtime_qc_priority_3,
        'Avg response time of all jobs': avg_rtime,
        'Avg response time of QC jobs': avg_rtime_qc, 
        'Avg response time of QC (Priority=1) jobs': avg_rtime_qc_priority_1,
        'Avg response time of QC (Priority<=2) jobs': avg_rtime_qc_priority_2,
        'Avg response time of QC (Priority<=3) jobs': avg_rtime_qc_priority_3,
    }.items(), columns=['Metric','Value'])
    st.dataframe(df, hide_index=True)
