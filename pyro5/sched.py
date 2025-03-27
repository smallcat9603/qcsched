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
        self.status = 'HOLD'
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
            job.starttime = datetime.now().strftime(DATE_FORMAT)
            
            if job.rscgroup in [RSCGROUP_QC_IBM, RSCGROUP_QC_QUAN]:
                server_qc = Pyro5.client.Proxy(URI_QC)
                server_qc.run(job.vid, job.relapsed)
            else:
                server_hpc = Pyro5.client.Proxy(URI_HPC)
                server_hpc.run(job.vid, job.relapsed)  


    def hold(self, subjoblist: list[Job]):
        for job in subjoblist:     
            job.status = 'HOLD'    


    def get_nnodes(self, name):
        return int(name.removesuffix(".sh").split('-')[-2])
    

    def get_elapsed(self, name):
        return int(name.removesuffix(".sh").split('-')[-1])


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
                    if i == 1 and 'hpc-' in name:
                        type = 1 # hpc
                    if 'hpc-a-' in name:
                        rscgroup = RSCGROUP_HPC_A
                        nnodes = self.get_nnodes(name)
                    elif 'hpc-b-' in name:
                        rscgroup = RSCGROUP_HPC_B
                        nnodes = self.get_nnodes(name)
                    elif 'hpc-c-' in name:
                        rscgroup = RSCGROUP_HPC_C
                        nnodes = self.get_nnodes(name)
                    elif 'qc-ibm-' in name:
                        rscgroup = RSCGROUP_QC_IBM
                        nnodes = 1
                    elif 'qc-quan-' in name:
                        rscgroup = RSCGROUP_QC_QUAN
                        nnodes = 1
                    elapsed = self.get_elapsed(name)

                    subjoblist.append(Job(name=name, rscgroup=rscgroup, id=self.id, vid=self.vid, type=type, nnodes=nnodes, elapsed=elapsed, priority=priority, relapsed=elapsed))
                self.joblist.append(subjoblist)

                job_vids = [job.vid for job in subjoblist]
                vids = ', '.join(map(str, job_vids))
                msg.append(f'The following qcsubjobs were submitted: {vids}')

        return '\n'.join(msg)
    

    @Pyro5.server.expose
    def add_subjob(self, vid: int, filename: str) -> str:
        '''
        add subjob to existing joblist
        '''
        msg = []

        with open(filename, 'r') as f: # pd.read_csv(filename, comment='#', header=None) cannot read data with different columns in different lines
            lines = f.readlines()
            for line in lines: # 2,hpc-a-3-4
                row = line.split(',')

                for i in range(1, len(row)):
                    name = row[i].strip() # remove \n
                    self.vid += 1
                    if i == 1 and 'qc-' in name:
                        return 'QC jobs cannot be added.'
                    if 'hpc-a-' in name:
                        rscgroup = RSCGROUP_HPC_A
                        nnodes = self.get_nnodes(name)
                    elif 'hpc-b-' in name:
                        rscgroup = RSCGROUP_HPC_B
                        nnodes = self.get_nnodes(name)
                    elif 'hpc-c-' in name:
                        rscgroup = RSCGROUP_HPC_C
                        nnodes = self.get_nnodes(name)
                    elapsed = self.get_elapsed(name)

                    found = False
                    for subjoblist in self.joblist:
                        for subjob in subjoblist:
                            if subjob.vid == vid:   
                                if subjob.status == 'HOLD':
                                    subjoblist.append(Job(name=name, rscgroup=rscgroup, id=subjob.id, vid=self.vid, type=1, nnodes=nnodes, elapsed=elapsed, priority=subjob.priority, relapsed=elapsed))
                                    msg.append(f'The following qcsubjob has been added: {self.vid}')
                                    
                                    found = True
                                    break
                                else:
                                    return f'Subjob cannot be added to {subjob.status} job.'      
                        if found:
                            break

                    if not found:
                        return f'Subjob {vid} not found.'

        return '\n'.join(msg)


    def sched_joblist(self):
        '''
        hold -> running
        '''
        while True:
            time.sleep(SCHED_INTERVAL)

            self.joblist.sort(key=self.sort_key_qc)

            for subjoblist in self.joblist:
                first_sub_job = subjoblist[0] # a qc sub job must be subjoblist[0]
                if first_sub_job.status == 'HOLD':
                    if first_sub_job.rscgroup == RSCGROUP_QC_IBM: 
                        if self.ibm_semaphor:
                            self.run(subjoblist)
                            self.ibm_semaphor = 0
                        else:
                            self.hold(subjoblist)
                    elif first_sub_job.rscgroup == RSCGROUP_QC_QUAN:
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
                    job.endtime = datetime.now().strftime(DATE_FORMAT)
                    if job.rscgroup == RSCGROUP_QC_IBM:
                        self.ibm_semaphor = 1
                    elif job.rscgroup == RSCGROUP_QC_QUAN:
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
                    starttime = datetime.strptime(job.starttime, DATE_FORMAT).replace(year=curtime.year)
                    difftime = curtime - starttime
                    if job.status == 'FINISH':
                        endtime = datetime.strptime(job.endtime, DATE_FORMAT).replace(year=curtime.year)
                        difftime = endtime - starttime
                    hours, remainder = divmod(difftime.seconds, 3600)
                    minutes, seconds = divmod(remainder, 60)
                    elapse = f'{hours:02}:{minutes:02}:{seconds:02}'
                    token = f'{hours*3600+minutes*60+seconds:.1f}'
                if job.rscgroup in [RSCGROUP_QC_IBM, RSCGROUP_QC_QUAN]:
                    job_status[job.vid] = [job.id, job.name, job.status, 'gz00', job.rscgroup, job.starttime, elapse, token, '-', str(job.nnodes)]
                else:
                    job_status[job.vid] = [job.id, job.name, job.status, 'gz00', job.rscgroup, job.starttime, elapse, token, str(job.nnodes), '-']
        return job_status   


    @Pyro5.server.expose
    def stat_qc_semaphor(self):   
        return f'ibm_semaphor: {self.ibm_semaphor}, quan_semaphor: {self.quan_semaphor}'     


def main():   

    daemon = Pyro5.api.Daemon(host=HOST_SCHED, port=PORT_SCHED)
    sched = Sched()

    thread = threading.Thread(target=sched.sched_joblist)
    thread.daemon = True # exits too if main thread exits
    thread.start()

    daemon.register(sched, objectId=NAME_SCHED)
    print("Scheduler running ...")
    daemon.requestLoop()


if __name__ == "__main__":
    main()
