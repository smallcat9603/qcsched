import Pyro5.api
import sys
import time
from setting import *


def main():
    if len(sys.argv) != 2:
        print("Usage: python qcsub.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]  

    server_sched = Pyro5.client.Proxy(URI_SCHED)

    print(server_sched.submit_joblist(filename))
    server_sched.sched_joblist()
    
    for i in range(10):
        time.sleep(2)
        print(server_sched.stat_job_status())
        print(server_sched.stat_qc_semaphor())


if __name__ == "__main__":
    main()
