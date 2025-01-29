#######################
####
#### python3 gen_jobs.py --> generate csv file
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

exp_to_real_1, exp_to_real_1_5, exp_to_real_2, exp_to_real_5, exp_to_real_10, exp_to_real_50, exp_to_real_100, exp_to_real_500, exp_to_real_1000 = 26_942, 31_422, 68_526, 249_048, 328_095, 272_384, 21_913, 165_838, 16_621 
total = sum([exp_to_real_1000, exp_to_real_500, exp_to_real_100, exp_to_real_50, exp_to_real_10, exp_to_real_5, exp_to_real_2, exp_to_real_1_5, exp_to_real_1])


prob_1 = exp_to_real_1/total
prob_1_5 = exp_to_real_1_5/total
prob_2 = exp_to_real_2/total
prob_5 = exp_to_real_5/total
prob_10 = exp_to_real_10/total
prob_50 = exp_to_real_50/total
prob_100 = exp_to_real_100/total
prob_500 = exp_to_real_500/total
prob_1000 = exp_to_real_1000/total


# parameters

qc_ratio = 0.1 # ratio of qc jobs to all jobs
nrows = 1000 # num of all jobs

def generate_job():
    hpc = random.choice(['HPC1', 'HPC2', 'HPC3'])
    type = random.choices(['HPC', 'QC1', 'QC2'], [1-qc_ratio, qc_ratio/2, qc_ratio/2])[0]
    nnodes = random.choices([1, 2, 4, 8, 16, 32, 64, 128, 256], 
                            [prob_n1, prob_n2, prob_n4, prob_n8, prob_n16, prob_n32, prob_n64, prob_n128, prob_n256])[0]
    relapsed = random.choices([1, 6, 12, 18, 24, 30, 36, 42, 48], 
                            [prob_t1, prob_t6, prob_t12, prob_t18, prob_t24, prob_t30, prob_t36, prob_t42, prob_t48])[0]
    start = random.randint(0, 24)
    priority = random.randint(1, 5)
    elapsed = relapsed * random.choices([1, 1.5, 2, 5, 10, 50, 100, 500, 1000], 
                            [prob_1, prob_1_5, prob_2, prob_5, prob_10, prob_50, prob_100, prob_500, prob_1000])[0]
    if elapsed > 48:
        elapsed = 48
    
    return [hpc, type, nnodes, elapsed, start, priority, relapsed]


def generate_csv():
    data = [generate_job() for _ in range(nrows)]
    df = pd.DataFrame(data, columns=['HPC', 'Type', 'Nodes', 'Time', 'Start', 'Priority', 'rTime'])
    now = datetime.datetime.now()
    timestamp = now.strftime("%Y%m%d%H%M")
    df.to_csv(f'JOBS_{timestamp}.txt', index=False, header=False)


generate_csv()