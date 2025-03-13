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
    

class Sched:
    def __init__(self):
        uname = platform.uname()
        system = uname[0]
        node = uname[1]
        release = uname[2]

        host_hpc = '192.168.3.69' if 'raspberrypi' in node else 'localhost'
        host_qc = '192.168.3.80' if 'raspberrypi' in node else 'localhost'
        port_hpc = 9091
        port_qc = 9092

        uri_qc = f"PYRO:QC-Job@{host_qc}:{port_qc}"
        uri_hpc = f"PYRO:HPC-Job@{host_hpc}:{port_hpc}"

        self.server_qc = Pyro5.client.Proxy(uri_qc)
        self.server_hpc = Pyro5.client.Proxy(uri_hpc)

        self.ibm_semaphor = 1
        self.quan_semaphor = 1

        self.joblist = []

    def sort_key_qc(self, subjoblist: list[Job]):
        return (subjoblist[0].type, subjoblist[0].priority, subjoblist[0].id)   

    def run(self, subjoblist: list[Job]):
        for job in subjoblist:
            job.status = 'RUNNING'
            if 'ibm-' or 'quan-' in job.rscgroup:
                self.server_qc.run(job.vid, job.relapsed) 
            else:
                self.server_hpc.run(job.vid, job.relapsed)  

    def queue(self, subjoblist: list[Job]):
        for job in subjoblist:     
            job.status = 'HOLD'     

    @Pyro5.server.expose
    def submit_joblist(self, filename: str):
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
                self.joblist.append(subjoblist)

        self.joblist.sort(key=self.sort_key_qc)

        for subjoblist in self.joblist:
            job_vids = [job.vid for job in subjoblist]
            vids = ', '.join(map(str, job_vids))
            print(f'The jobs were submitted: {vids}')
            if 'ibm-' in subjoblist[0].rscgroup:
                if ibm_semaphor:
                    self.run(subjoblist)
                    ibm_semaphor = 0
                else:
                    self.queue(subjoblist)
            elif 'quan-' in subjoblist[0].rscgroup:
                if quan_semaphor:
                    self.run(subjoblist)
                    quan_semaphor = 0
                else:
                    self.queue(subjoblist)
            else:
                self.run(subjoblist)

    @Pyro5.server.expose
    def finish(self, vid):    
        for subjoblist in self.joblist:
            for job in subjoblist:
                if vid == job.vid:
                    job.status = 'FINISH'  

def main():   

    uname = platform.uname()
    system = uname[0]
    node = uname[1]
    release = uname[2]   

    host = '192.168.3.13' if 'raspberrypi' in node else None
    port = 9093

    daemon = Pyro5.api.Daemon(host=host, port=port)
    uri = daemon.register(Sched, objectId="Sched")
    print("Scheduler running ...")
    daemon.requestLoop()

if __name__ == "__main__":
    main()
