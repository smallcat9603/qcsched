import Pyro5.api
import platform
import sys
import time

def main():
    if len(sys.argv) != 2:
        print("Usage: python qcsub.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]  

    uname = platform.uname()
    system = uname[0]
    node = uname[1]
    release = uname[2]

    host_sched = '192.168.3.13' if 'raspberrypi' in node else 'localhost'
    port_sched = 9093

    uri_sched = f"PYRO:sched@{host_sched}:{port_sched}"
    server_sched = Pyro5.client.Proxy(uri_sched)

    server_sched.submit_joblist(filename)
    
    for i in range(10):
        time.sleep(2)
        print(server_sched.stat_job_status())
        print(server_sched.stat_qc_semaphor())

if __name__ == "__main__":
    main()
