#######################
####
#### python3 gen_wo_rtime.py --> generate csv file
####
#######################

import random
import pandas as pd
import datetime


# data from AXIES 2024 'スーパーコンピュータ WisteriaBDEC-01 における利用状況を考慮した運用の再検討'

n1, n2, n4, n8, n16, n32, n64, n128, n256, n512, n1024 = 488_476, 79_250, 134_038, 124_586, 167_627, 49_502, 42_705, 54_301, 38_885, 3_046, 848
njobs_n = sum([n1, n2, n4, n8, n16, n32, n64, n128, n256])

prob_n1 = n1/njobs_n
prob_n2 = n2/njobs_n
prob_n4 = n4/njobs_n
prob_n8 = n8/njobs_n
prob_n16 = n16/njobs_n
prob_n32 = n32/njobs_n
prob_n64 = n64/njobs_n
prob_n128 = n128/njobs_n
prob_n256 = n256/njobs_n

t1, t6, t12, t18, t24, t30, t36, t42, t48 = 1_009_000, 107_205, 30_421, 9_559, 8_993, 3_955, 2_784, 2_103, 9_244
njobs_t = sum([t1, t6, t12, t18, t24, t30, t36, t42, t48])

prob_t1 = t1/njobs_t
prob_t6 = t6/njobs_t
prob_t12 = t12/njobs_t
prob_t18 = t18/njobs_t
prob_t24 = t24/njobs_t
prob_t30 = t30/njobs_t
prob_t36 = t36/njobs_t
prob_t42 = t42/njobs_t
prob_t48 = t48/njobs_t

# parameters

qc_ratio = 0.1
nrows = 1000


def generate_job():
    hpc = random.choice(['HPC1', 'HPC2', 'HPC3'])
    type = random.choices(['HPC', 'QC1', 'QC2'], [1-qc_ratio, qc_ratio/2, qc_ratio/2])[0]
    nnodes = random.choices([1, 2, 4, 8, 16, 32, 64, 128, 256], 
                            [prob_n1, prob_n2, prob_n4, prob_n8, prob_n16, prob_n32, prob_n64, prob_n128, prob_n256])[0]
    elapsed = random.choices([1, 6, 12, 18, 24, 30, 36, 42, 48], 
                            [prob_t1, prob_t6, prob_t12, prob_t18, prob_t24, prob_t30, prob_t36, prob_t42, prob_t48])[0]
    start = random.randint(0, 24)
    priority = random.randint(1, 5)
    
    return [hpc, type, nnodes, elapsed, start, priority]


def generate_csv():
    data = [generate_job() for _ in range(nrows)]
    df = pd.DataFrame(data, columns=['HPC', 'Type', 'Nodes', 'Time', 'Start', 'Priority'])
    now = datetime.datetime.now()
    timestamp = now.strftime("%Y%m%d%H%M")
    df.to_csv(f'wor_{timestamp}.txt', index=False, header=False)


generate_csv()