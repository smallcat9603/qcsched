import streamlit as st
import pandas as pd
import random


def get_num_from_0(text: str) -> int:
    return int(text[-1]) - 1


def get_id(src: int, type: str) -> str:
    if type.startswith('QC'):
        qc = get_num_from_0(type)
        # record qc job count at src hpc
        st.session_state[f"job_manager_{src}"].qc_cnt[qc] += 1
        # allocate qc job id, ex 191000
        id = f'{src+1}9{qc+1}{st.session_state[f"job_manager_{src}"].qc_cnt[qc]:03d}'            
    else:
        # record hpc job count at src hpc
        st.session_state[f'job_manager_{src}'].hpc_cnt += 1 
        # allocate hpc job id, ex 100000
        id = f'{src+1}{st.session_state[f"job_manager_{src}"].hpc_cnt:05d}' 
    return id


def get_vid(src: int) -> str:
    vid = st.session_state[f"job_manager_{src}"].hpc_cnt
    for i in range(st.session_state['NUM_QC']):
        vid += st.session_state[f"job_manager_{src}"].qc_cnt[i]
    vid = f'{src+1}{vid:05d}'            
    return vid


def color(text: str) -> str:
    colors = {
        'hpc': 'green',
        'qc1': 'orange',
        'qc2': 'violet',
    }
    return colors[text]


def show_submitted_jobs():
    st.header('Jobs Submitted')

    tabs = st.tabs([f'HPC{src+1}' for src in range(st.session_state['NUM_HPC'])])

    for src in range(st.session_state['NUM_HPC']):
        if st.session_state[f'job_manager_{src}'].jobs_submitted:
            df = pd.DataFrame([{
                'Job ID': job.vid,
                'Status': job.status,
                'Type': job.type,
                'Nodes': job.nnodes,
                'Elapsed Time': job.elapsed,
                # 'Real Elapsed Time': job.relapsed,
                'Start Time': job.start,
                'Priority': job.priority,
                } for job in st.session_state[f'job_manager_{src}'].jobs_submitted]
            )
            tabs[src].dataframe(df, hide_index=True)
        else:
            tabs[src].info(f'No job submitted in HPC{src+1}')


def qc_util(arr_qc_flag) -> list: 
    """
    Provide list of (start, len) according to non-zero elements in qc_flag 
    """
    result = []
    cur = 0
    for t in range(1, st.session_state['SCHED_MAP_TIME']):
        if arr_qc_flag[t] != arr_qc_flag[t-1]:
            if arr_qc_flag[t-1] > 0:
                result.append((cur, t-cur))
            if arr_qc_flag[t] > 0:
                cur = t
        if t == st.session_state['SCHED_MAP_TIME'] - 1 and arr_qc_flag[t] > 0:
            result.append((cur, t-cur+1))
    return result


def cdf(data_list: list):
    all_sum = sum(data_list)
    cur_sum = 0
    for prob in data_list:
        cur_sum += prob
        yield cur_sum/all_sum


def rtime(elapsed: int) -> int:
    # data from AXIES 2024 'スーパーコンピュータ WisteriaBDEC-01 における利用状況を考慮した運用の再検討'
    exp_to_real_1, exp_to_real_1_5, exp_to_real_2, exp_to_real_5, exp_to_real_10, exp_to_real_50, exp_to_real_100, exp_to_real_500, exp_to_real_1000 = 26_942, 31_422, 68_526, 249_048, 328_095, 272_384, 21_913, 165_838, 16_621 
    data_list = [exp_to_real_1000, exp_to_real_500, exp_to_real_100, exp_to_real_50, exp_to_real_10, exp_to_real_5, exp_to_real_2, exp_to_real_1_5, exp_to_real_1]
    fin_1000, fin_500, fin_100, fin_50, fin_10, fin_5, fin_2, fin_1_5, fin_1 = cdf(data_list)

    r = random.random()
    if r < fin_1000:
        return elapsed//1000 + 1
    elif r < fin_500:
        return elapsed//500 + 1
    elif r < fin_100:
        return elapsed//100 + 1
    elif r < fin_50:
        return elapsed//50 + 1
    elif r < fin_10:
        return elapsed//10 + 1
    elif r < fin_5:
        return elapsed//5 + 1
    elif r < fin_2:
        return elapsed//2 + 1
    elif r < fin_1_5:
        return elapsed*2//3 + 1
    elif r < fin_1:
        return elapsed

