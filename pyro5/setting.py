import platform

def get_node_info():

    uname = platform.uname()
    system = uname[0]
    node = uname[1]
    release = uname[2]  

    return node


NAME_HPC = 'hpc'
NAME_QC = 'qc'
NAME_SCHED = 'sched'

HOST_HPC = '192.168.3.69' if 'raspberrypi' in get_node_info() else 'localhost'
HOST_QC = '192.168.3.80' if 'raspberrypi' in get_node_info() else 'localhost'
HOST_SCHED = '192.168.3.13' if 'raspberrypi' in get_node_info() else 'localhost'

PORT_HPC = 9092
PORT_QC = 9093
PORT_SCHED = 9091

URI_HPC = f"PYRO:{NAME_HPC}@{HOST_HPC}:{PORT_HPC}"
URI_QC = f"PYRO:{NAME_QC}@{HOST_QC}:{PORT_QC}"
URI_SCHED = f"PYRO:{NAME_SCHED}@{HOST_SCHED}:{PORT_SCHED}"

SCHED_INTERVAL = 10
