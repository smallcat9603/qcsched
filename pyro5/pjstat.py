import Pyro5.api
import sys
import time
from setting import *


def main():
    if len(sys.argv) != 1:
        print("Usage: python pjstat.py")
        sys.exit(1)

    server_sched = Pyro5.client.Proxy(URI_SCHED)

    print(server_sched.stat_job_status())
    print(server_sched.stat_qc_semaphor())


if __name__ == "__main__":
    main()
