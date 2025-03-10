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
    wtime_hpc = [[] for _ in range(5)]
    wtime_qc = [[] for _ in range(5)]
    rtime_hpc = [[] for _ in range(5)]
    rtime_qc = [[] for _ in range(5)]

    for src in range(st.session_state['NUM_HPC']):
        for job in st.session_state[f'job_manager_{src}'].jobs_submitted:
            if job.type.startswith('QC'):
                wtime_qc[job.priority-1].append(job.wait)
                rtime_qc[job.priority-1].append(job.wait+job.info[2])
            else:
                wtime_hpc[job.priority-1].append(job.wait)
                rtime_hpc[job.priority-1].append(job.wait+job.info[2])

    wtime_hpc_list = [element for sublist in wtime_hpc for element in sublist]
    wtime_qc_list = [element for sublist in wtime_qc for element in sublist]
    rtime_hpc_list = [element for sublist in rtime_hpc for element in sublist]
    rtime_qc_list = [element for sublist in rtime_qc for element in sublist]
    
    avg_wtime_hpc = float(round(np.mean(wtime_hpc_list), 1)) if wtime_hpc_list else 0.0
    avg_wtime_qc = float(round(np.mean(wtime_qc_list), 1)) if wtime_qc_list else 0.0 

    avg_rtime_hpc = float(round(np.mean(rtime_hpc_list), 1)) if rtime_hpc_list else 0.0
    avg_rtime_qc = float(round(np.mean(rtime_qc_list), 1)) if rtime_qc_list else 0.0


    st.header('Simulation Result')
    df = pd.DataFrame({
        'Num of All jobs': len(wtime_hpc_list) + len(wtime_qc_list),
        'Num of HPC jobs': len(wtime_hpc_list),
        'Num of QC jobs': len(wtime_qc_list),

        'Wait time of HPC jobs': (min(wtime_hpc_list), avg_wtime_hpc, max(wtime_hpc_list)),
        'Wait time of HPC (Priority=1) jobs': (min(wtime_hpc[0]), round(sum(wtime_hpc[0])/len(wtime_hpc[0]), 1), max(wtime_hpc[0]), wtime_hpc[0]),
        'Wait time of HPC (Priority=2) jobs': (min(wtime_hpc[1]), round(sum(wtime_hpc[1])/len(wtime_hpc[1]), 1), max(wtime_hpc[1]), wtime_hpc[1]),
        'Wait time of HPC (Priority=3) jobs': (min(wtime_hpc[2]), round(sum(wtime_hpc[2])/len(wtime_hpc[2]), 1), max(wtime_hpc[2]), wtime_hpc[2]),
        'Wait time of HPC (Priority=4) jobs': (min(wtime_hpc[3]), round(sum(wtime_hpc[3])/len(wtime_hpc[3]), 1), max(wtime_hpc[3]), wtime_hpc[3]),
        'Wait time of HPC (Priority=5) jobs': (min(wtime_hpc[4]), round(sum(wtime_hpc[4])/len(wtime_hpc[4]), 1), max(wtime_hpc[4]), wtime_hpc[4]),

        'Wait time of QC jobs': (min(wtime_qc_list), avg_wtime_qc, max(wtime_qc_list)),
        'Wait time of QC (Priority=1) jobs': (min(wtime_qc[0]), round(sum(wtime_qc[0])/len(wtime_qc[0]), 1), max(wtime_qc[0]), wtime_qc[0]),
        'Wait time of QC (Priority=2) jobs': (min(wtime_qc[1]), round(sum(wtime_qc[1])/len(wtime_qc[1]), 1), max(wtime_qc[1]), wtime_qc[1]),
        'Wait time of QC (Priority=3) jobs': (min(wtime_qc[2]), round(sum(wtime_qc[2])/len(wtime_qc[2]), 1), max(wtime_qc[2]), wtime_qc[2]),
        'Wait time of QC (Priority=4) jobs': (min(wtime_qc[3]), round(sum(wtime_qc[3])/len(wtime_qc[3]), 1), max(wtime_qc[3]), wtime_qc[3]),
        'Wait time of QC (Priority=5) jobs': (min(wtime_qc[4]), round(sum(wtime_qc[4])/len(wtime_qc[4]), 1), max(wtime_qc[4]), wtime_qc[4]),

        'Response time of HPC jobs': (min(rtime_hpc_list), avg_rtime_hpc, max(rtime_hpc_list)),
        'Response time of HPC (Priority=1) jobs': (min(rtime_hpc[0]), round(sum(rtime_hpc[0])/len(rtime_hpc[0]), 1), max(rtime_hpc[0]), rtime_hpc[0]),
        'Response time of HPC (Priority=2) jobs': (min(rtime_hpc[1]), round(sum(rtime_hpc[1])/len(rtime_hpc[1]), 1), max(rtime_hpc[1]), rtime_hpc[1]),
        'Response time of HPC (Priority=3) jobs': (min(rtime_hpc[2]), round(sum(rtime_hpc[2])/len(rtime_hpc[2]), 1), max(rtime_hpc[2]), rtime_hpc[2]),
        'Response time of HPC (Priority=4) jobs': (min(rtime_hpc[3]), round(sum(rtime_hpc[3])/len(rtime_hpc[3]), 1), max(rtime_hpc[3]), rtime_hpc[3]),
        'Response time of HPC (Priority=5) jobs': (min(rtime_hpc[4]), round(sum(rtime_hpc[4])/len(rtime_hpc[4]), 1), max(rtime_hpc[4]), rtime_hpc[4]),

        'Response time of QC jobs': (min(rtime_qc_list), avg_rtime_qc, max(rtime_qc_list)),
        'Response time of QC (Priority=1) jobs': (min(rtime_qc[0]), round(sum(rtime_qc[0])/len(rtime_qc[0]), 1), max(rtime_qc[0]), rtime_qc[0]),
        'Response time of QC (Priority=2) jobs': (min(rtime_qc[1]), round(sum(rtime_qc[1])/len(rtime_qc[1]), 1), max(rtime_qc[1]), rtime_qc[1]),
        'Response time of QC (Priority=3) jobs': (min(rtime_qc[2]), round(sum(rtime_qc[2])/len(rtime_qc[2]), 1), max(rtime_qc[2]), rtime_qc[2]),
        'Response time of QC (Priority=4) jobs': (min(rtime_qc[3]), round(sum(rtime_qc[3])/len(rtime_qc[3]), 1), max(rtime_qc[3]), rtime_qc[3]),
        'Response time of QC (Priority=5) jobs': (min(rtime_qc[4]), round(sum(rtime_qc[4])/len(rtime_qc[4]), 1), max(rtime_qc[4]), rtime_qc[4]),
    }.items(), columns=['Metric','Value'])
    st.dataframe(df, hide_index=True)
