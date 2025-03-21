import Pyro5.api
import time
from setting import *


class Job:
    def __init__(self, rscgroup: str, id: int, vid: int, type: int, nnodes: int, elapsed: int, priority: int, relapsed: int):
        self.rscgroup = rscgroup # qc-hpc-a-m-o, qc-hpc-b-m-o, qc-hpc-c-m-o, qc-hpc-ibm-a, qc-hpc-quan-a
        self.id = id # job id
        self.vid = vid # subjob id
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
        self.ibm_semaphor = 1
        self.quan_semaphor = 1

        self.joblist_unsorted = []
        self.joblist_sorted = []

        self.id = 0
        self.vid = 100000


    def sort_key_qc(self, subjoblist: list[Job]):
        return (subjoblist[0].type, subjoblist[0].priority, subjoblist[0].id)   
    

    def run(self, subjoblist: list[Job]):
        for job in subjoblist:
            job.status = 'RUNNING'
            
            if 'ibm-' in job.rscgroup or 'quan-' in job.rscgroup:
                server_qc = Pyro5.client.Proxy(URI_QC)
                server_qc.run(job.vid, job.relapsed)
            else:
                server_hpc = Pyro5.client.Proxy(URI_HPC)
                server_hpc.run(job.vid, job.relapsed)  


    def hold(self, subjoblist: list[Job]):
        for job in subjoblist:     
            job.status = 'HOLD'    


    @Pyro5.server.expose
    def submit_joblist(self, filename: str) -> str:
        msg = []

        with open(filename, 'r') as f: # pd.read_csv(filename, comment='#', header=None) cannot read data with different columns in different lines
            lines = f.readlines()
            for line in lines: # 2,ibm-4,hpc-a-3-4
                row = line.split(',')
                subjoblist = []

                self.id += 1
                priority = int(row[0])
                type = 0 # qc
                for i in range(1, len(row)):
                    self.vid += 1
                    if i == 1 and 'hpc-' in row[i]:
                        type = 1 # hpc
                    if 'hpc-a-' in row[i]:
                        rscgroup = 'qc-hpc-a-m-o'
                        nnodes = int(row[i].split('-')[-2])
                    elif 'hpc-b-' in row[i]:
                        rscgroup = 'qc-hpc-b-m-o'
                        nnodes = int(row[i].split('-')[-2])
                    elif 'hpc-c-' in row[i]:
                        rscgroup = 'qc-hpc-c-m-o'
                        nnodes = int(row[i].split('-')[-2])
                    elif 'ibm-' in row[i]:
                        rscgroup = 'qc-hpc-ibm-a'
                        nnodes = 1
                    elif 'quan-' in row[i]:
                        rscgroup = 'qc-hpc-quan-a'
                        nnodes = 1
                    elapsed = int(row[i].split('-')[-1])

                    subjoblist.append(Job(rscgroup=rscgroup, id=self.id, vid=self.vid, type=type, nnodes=nnodes, elapsed=elapsed, priority=priority, relapsed=elapsed))
                self.joblist_unsorted.append(subjoblist)

                job_vids = [job.vid for job in subjoblist]
                vids = ', '.join(map(str, job_vids))
                msg.append(f'The jobs were submitted: {vids}')

        return '\n'.join(msg)
    

    @Pyro5.server.expose
    @Pyro5.server.oneway
    def sched_joblist(self):    

        self.joblist_sorted = self.joblist_unsorted.copy()
        self.joblist_sorted.sort(key=self.sort_key_qc)

        while True:
            has_unfinished = False
            for subjoblist in self.joblist_sorted:
                first_sub_job = subjoblist[0] # a qc sub job must be subjoblist[0]
                if first_sub_job.status in ['ACCEPT', 'HOLD']:
                    has_unfinished = True
                    if 'ibm-' in first_sub_job.rscgroup: 
                        if self.ibm_semaphor:
                            self.run(subjoblist)
                            self.ibm_semaphor = 0
                        else:
                            self.hold(subjoblist)
                    elif 'quan-' in first_sub_job.rscgroup:
                        if self.quan_semaphor:
                            self.run(subjoblist)
                            self.quan_semaphor = 0
                        else:
                            self.hold(subjoblist)
                    else:
                        self.run(subjoblist)
            if not has_unfinished:
                break


    @Pyro5.server.expose
    def finish(self, vid):    
        for subjoblist in self.joblist_sorted:
            for job in subjoblist:
                if vid == job.vid:
                    job.status = 'FINISH'
                    if 'ibm-' in job.rscgroup:
                        self.ibm_semaphor = 1
                    elif 'quan-' in job.rscgroup:
                        self.quan_semaphor = 1


    @Pyro5.server.expose
    def stat_job_status(self):   
        job_status = {}
        for subjoblist in self.joblist_unsorted:
            for job in subjoblist:
                job_status[job.vid] = job.status
        return job_status   


    @Pyro5.server.expose
    def stat_qc_semaphor(self):   
        return f'ibm_semaphor: {self.ibm_semaphor}, quan_semaphor: {self.quan_semaphor}'     


def main():   

    daemon = Pyro5.api.Daemon(host=HOST_SCHED, port=PORT_SCHED)
    sched = Sched()
    daemon.register(sched, objectId="sched")
    print("Scheduler running ...")
    daemon.requestLoop()


if __name__ == "__main__":
    main()
