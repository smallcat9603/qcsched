import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches

from utils import color, qc_util, get_num_from_0


def plot_qc(mode):
    for qc in range(st.session_state['NUM_QC']):
        # st.write(st.session_state['semaphore'].qc_flag[i])
        if np.any(st.session_state['semaphore'].qc_flag[qc] > 0):
            fig, ax = plt.subplots(figsize=(8, 1))
            # ax.set_title(f'QC{i+1}')
            ax.set_ylim(0, 1)
            ax.set_xlim(0, st.session_state['SCHED_MAP_TIME'])
            ax.set_xlabel('Time')
            ax.set_ylabel('Utilization')
            ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')

            fc = color(f'qc{qc+1}')
            utils = qc_util(st.session_state['semaphore'].qc_flag[qc])
            for util in utils:
                rect = patches.Rectangle((util[0], 0), util[1], 1, edgecolor='black', facecolor=fc)
                ax.add_patch(rect)
                if mode == 'Demo':
                    ax.text(util[0]+util[1]/2, 1/2, f'{st.session_state['semaphore'].qc_flag[qc][util[0]]}', size=10, horizontalalignment='center', verticalalignment='center')   

            st.header(f'QC{qc+1} Schedule')
            st.pyplot(fig) 


def plot_hpc(mode):
    for src in range(st.session_state['NUM_HPC']):
        if st.session_state[f'job_manager_{src}'].jobs_scheduled and any(job.status == 'RUNNING' for job in st.session_state[f'job_manager_{src}'].jobs_scheduled):
            fig, ax = plt.subplots(figsize=(8, 4))
            # ax.set_title(f'HPC{src+1}')
            ax.set_ylim(0, st.session_state['HPC_NODES'])
            ax.set_xlim(0, st.session_state['SCHED_MAP_TIME'])
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
                
                    rect = patches.Rectangle(job.map, job.elapsed, job.nnodes, edgecolor='black', facecolor=fc)
                    ax.add_patch(rect)
                    if mode == 'Demo':
                        ax.text(job.map[0]+job.elapsed/2, job.map[1]+job.nnodes/2, f'{job.vid}-{job.priority}-{job.type}', size=10, horizontalalignment='center', verticalalignment='center')

            st.header(f'HPC{src+1} Schedule')
            st.pyplot(fig)


def plot(mode):
    plot_qc(mode)
    plot_hpc(mode)
