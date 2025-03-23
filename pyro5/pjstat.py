import Pyro5.api
import sys
from setting import *


def main():
    if len(sys.argv) != 1:
        print("Usage: python pjstat.py")
        sys.exit(1)

    server_sched = Pyro5.client.Proxy(URI_SCHED)

    stat = []
    job_status = server_sched.stat_job_status()
    for key, values in job_status.items():
        if values[1] != 'FINISH':
            status = '\t'.join(values)
            stat.append(f'{key}\t{status}') 
    
    if stat:
        print('JOB_ID\tJOB_NAME\tSTATUS\tPROJECT\tRSCGROUP\tSTART_DATE\tELAPSE\t\tTOKEN\tNODE\tGPU')
        print('\n'.join(stat))
    else:
        print('No jobs were submitted.')

    # print(server_sched.stat_qc_semaphor())


if __name__ == "__main__":
    main()
