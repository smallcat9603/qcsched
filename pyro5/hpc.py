import Pyro5.api
import time
import platform

@Pyro5.api.expose
class Job(object):
    def run(self, t):
        time.sleep(t)
        return "HPC-Job finished."
    
uname = platform.uname()
print(uname)
system = uname[0]
node = uname[1]
release = uname[2]

host = None
port = 9091
if 'raspberry' in release:
    host = '192.168.3.69'

daemon = Pyro5.api.Daemon(host=host, port=port)
uri = daemon.register(Job, objectId="HPC-Job")
print("HPC running ...")
daemon.requestLoop()
