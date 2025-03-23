import Pyro5.api
import threading
import time
from datetime import datetime
from setting import *


class Job:
    def __init__(self, name: str, rscgroup: str, id: int, vid: int, type: int, nnodes: int, elapsed: int, priority: int, relapsed: int):
        self.name = name
        self.rscgroup = rscgroup # qc-hpc-a-m-o, qc-hpc-b-m-o, qc-hpc-c-m-o, qc-hpc-ibm-a, qc-hpc-quan-a
        self.id = id # job id
        self.vid = vid # subjob id
        self.type = type # 0 for QC, 1 for HPC, all sub jobs have the same type
        self.nnodes = nnodes # number of HPC nodes
        self.elapsed = elapsed # expected elasped time
        self.priority = priority # 1 for highest priority
        self.relapsed = relapsed # real elapsed time
        self.wait = 0
        self.status = 'ACCEPT'
        self.timestamp = time.time()
        self.starttime = '--/-- --:--:--'
        self.endtime = '--/-- --:--:--'

        
    def __repr__(self):
        return f'Job(name={self.name}, rscgroup={self.rscgroup}, id={self.id}, vid={self.vid}, type={self.type}, nnodes={self.nnodes}, elapsed={self.elapsed}, priority={self.priority}, relapsed={self.relapsed}, wait={self.wait}, status={self.status}, timestamp={self.timestamp}, starttime={self.starttime}, endtime={self.endtime})'
    

class Sched:
    def __init__(self):
        self.ibm_semaphor = 1
        self.quan_semaphor = 1

        self.joblist = [] # list[subjoblist]

        self.id = 0
        self.vid = 100000


    def sort_key_qc(self, subjoblist: list[Job]):
        return (subjoblist[0].type, subjoblist[0].priority, subjoblist[0].id)   
    

    def run(self, subjoblist: list[Job]):
        for job in subjoblist:
            job.status = 'RUNNING'
            job.wait = int(time.time() - job.timestamp)
            job.starttime = datetime.now().strftime('%m/%d %H:%M:%S')
            
            if job.rscgroup in ['qc-hpc-ibm-a', 'qc-hpc-quan-a']:
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
        '''
        put into joblist
        '''
        msg = []

        with open(filename, 'r') as f: # pd.read_csv(filename, comment='#', header=None) cannot read data with different columns in different lines
            lines = f.readlines()
            for line in lines: # 2,qc-ibm-4,hpc-a-3-4
                row = line.split(',')
                subjoblist = []

                self.id += 1
                priority = int(row[0])
                type = 0 # qc
                for i in range(1, len(row)):
                    name = row[i].strip() # remove \n
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
                    elif 'qc-ibm-' in row[i]:
                        rscgroup = 'qc-hpc-ibm-a'
                        nnodes = 1
                    elif 'qc-quan-' in row[i]:
                        rscgroup = 'qc-hpc-quan-a'
                        nnodes = 1
                    elapsed = int(row[i].split('-')[-1])

                    subjoblist.append(Job(name=name, rscgroup=rscgroup, id=self.id, vid=self.vid, type=type, nnodes=nnodes, elapsed=elapsed, priority=priority, relapsed=elapsed))
                self.joblist.append(subjoblist)

                job_vids = [job.vid for job in subjoblist]
                vids = ', '.join(map(str, job_vids))
                msg.append(f'The jobs were submitted: {vids}')

        return '\n'.join(msg)
    

    def sched_joblist(self):
        '''
        accept -> running
        accept -> hold
        hold -> running
        '''
        while True:
            time.sleep(SCHED_INTERVAL)

            self.joblist.sort(key=self.sort_key_qc)

            for subjoblist in self.joblist:
                first_sub_job = subjoblist[0] # a qc sub job must be subjoblist[0]
                if first_sub_job.status in ['ACCEPT', 'HOLD']:
                    if first_sub_job.rscgroup == 'qc-hpc-ibm-a': 
                        if self.ibm_semaphor:
                            self.run(subjoblist)
                            self.ibm_semaphor = 0
                        else:
                            self.hold(subjoblist)
                    elif first_sub_job.rscgroup == 'qc-hpc-quan-a':
                        if self.quan_semaphor:
                            self.run(subjoblist)
                            self.quan_semaphor = 0
                        else:
                            self.hold(subjoblist)
                    else:
                        self.run(subjoblist)


    @Pyro5.server.expose
    def finish(self, vid): 
        found = False
        for subjoblist in self.joblist:
            for job in subjoblist:
                if vid == job.vid:
                    found = True
                    job.status = 'FINISH'
                    job.endtime = datetime.now().strftime('%m/%d %H:%M:%S')
                    if job.rscgroup == 'qc-hpc-ibm-a':
                        self.ibm_semaphor = 1
                    elif job.rscgroup == 'qc-hpc-quan-a':
                        self.quan_semaphor = 1
                    break
            if found:
                break
        

    @Pyro5.server.expose
    def stat_job_status(self):
        job_status = {}
        joblist_ = sorted(self.joblist, key=lambda subjoblist: subjoblist[0].id)
        for subjoblist in joblist_:
            for job in subjoblist:
                elapse = '00:00:00'
                token = f'({job.elapsed:.1f})'
                if job.status in ['RUNNING', 'FINISH']:
                    curtime = datetime.now()
                    starttime = datetime.strptime(job.starttime, "%m/%d %H:%M:%S").replace(year=curtime.year)
                    difftime = curtime - starttime
                    if job.status == 'FINISH':
                        endtime = datetime.strptime(job.endtime, "%m/%d %H:%M:%S").replace(year=curtime.year)
                        difftime = endtime - starttime
                    hours, remainder = divmod(difftime.seconds, 3600)
                    minutes, seconds = divmod(remainder, 60)
                    elapse = f'{hours:02}:{minutes:02}:{seconds:02}'
                    token = f'{hours*3600+minutes*60+seconds:.1f}'
                if job.rscgroup in ['qc-hpc-ibm-a', 'qc-hpc-quan-a']:
                    job_status[job.vid] = [job.name, job.status, 'gz00', job.rscgroup, job.starttime, elapse, token, '-', str(job.nnodes)]
                else:
                    job_status[job.vid] = [job.name, job.status, 'gz00', job.rscgroup, job.starttime, elapse, token, str(job.nnodes), '-']
        return job_status   


    @Pyro5.server.expose
    def stat_qc_semaphor(self):   
        return f'ibm_semaphor: {self.ibm_semaphor}, quan_semaphor: {self.quan_semaphor}'     


def main():   

    daemon = Pyro5.api.Daemon(host=HOST_SCHED, port=PORT_SCHED)
    sched = Sched()

    thread = threading.Thread(target=sched.sched_joblist)
    thread.start()

    daemon.register(sched, objectId=NAME_SCHED)
    print("Scheduler running ...")
    daemon.requestLoop()


if __name__ == "__main__":
    main()
