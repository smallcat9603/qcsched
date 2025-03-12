import Pyro5.api
import platform
import pandas as pd
import time

class Job:
    def __init__(self, rscgroup: str, id: int, vid: int, type: int, nnodes: int, elapsed: int, priority: int, relapsed: int):
        self.rscgroup = rscgroup # qc-hpc-a-m-o, qc-hpc-b-m-o, qc-hpc-c-m-o, qc-hpc-ibm-a, qc-hpc-quan-a
        self.id = id # hybrid id
        self.vid = vid # for display, cannot identify qc job or not
        self.type = type # "QC" == 0 or "HPC" == 1
        self.nnodes = nnodes # number of HPC nodes
        self.elapsed = elapsed # expected elasped time
        self.priority = priority # 1 for highest priority
        self.relapsed = relapsed # real elapsed time
        self.wait = 0 # add 1 until running starts
        self.status = 'ACCEPT'
        self.timestamp = time.time()
        
    def __repr__(self):
        return f'Job(rscgroup={self.rscgroup}, id={self.id}, vid={self.vid}, type={self.type}, nnodes={self.nnodes}, elapsed={self.elapsed}, priority={self.priority}, relapsed={self.relapsed}, wait={self.wait}, status={self.status}, timestamp={self.timestamp})'


def import_joblist(filename):
    joblist = []
    id = 0
    vid = 100000
    with open(filename, 'r') as f: # pd.read_csv(filename, comment='#', header=None) cannot read data with different columns in different lines
        lines = f.readlines()
        for line in lines:
            row = line.split(',')
            priority = int(row[0])
            id += 1
            subjoblist = []
            type = 0
            for i in range(1, len(row)):
                if i == 1 and 'hpc-' in row[i]:
                    type = 1
                vid += 1
                rscgroup = ''
                if 'hpc-a-' in row[i]:
                    rscgroup = 'qc-hpc-a-m-o'
                elif 'hpc-b-' in row[i]:
                    rscgroup = 'qc-hpc-b-m-o'
                elif 'hpc-c-' in row[i]:
                    rscgroup = 'qc-hpc-c-m-o'
                elif 'ibm-' in row[i]:
                    rscgroup = 'qc-hpc-ibm-a'
                elif 'quan-' in row[i]:
                    rscgroup = 'qc-hpc-quan-a'
                nnodes = row[i].split('-')[-2]
                elapsed = row[i].split('-')[-1]
                relapsed = elapsed
                subjoblist.append(Job(rscgroup=rscgroup, id=id, vid=vid, type=type, nnodes=nnodes, elapsed=elapsed, priority=priority, relapsed=relapsed))
            joblist.append(subjoblist)
    return joblist

def sort_key_qc(subjoblist: list[Job]):
    return (subjoblist[0].type, subjoblist[0].priority, subjoblist[0].id)

joblist = import_joblist('test.txt')
joblist.sort(key=sort_key_qc)

uname = platform.uname()
system = uname[0]
node = uname[1]
release = uname[2]

ibm_semaphor = 1
quan_semaphor = 1

for subjoblist in joblist:
    if 'ibm-' in subjoblist[0].rscgroup:
        if ibm_semaphor:
            submit(subjoblist)
            ibm_semaphor = 0
        else:
            set_queued(subjoblist)
    elif 'quan' in subjoblist[0].rscgroup:
        if quan_semaphor:
            submit(subjoblist)
            quan_semaphor = 0
        else:
            set_queued(subjoblist)
    else:
        submit(subjoblist)



# host = 'localhost'
# port = 9091
# if 'raspberry' in release:
#     host = '192.168.3.69'

# uri = f"PYRO:HPC-Job@{host}:{port}"
# server = Pyro5.client.Proxy(uri)
# t = int(input("The time of the HPC job? ").strip())
# print(server.run(t))

# host = 'localhost'
# port = 9092
# if 'raspberry' in release:
#     host = '192.168.3.80'

# uri = f"PYRO:QC-Job@{host}:{port}"
# server = Pyro5.client.Proxy(uri)
# t = int(input("The time of the QC job? ").strip())
# print(server.run(t))
