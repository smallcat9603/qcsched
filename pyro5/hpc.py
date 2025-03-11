import Pyro5.api
import time

@Pyro5.api.expose
class Job(object):
    def run(self, t):
        time.sleep(t)
        return "HPC-Job finished."

daemon = Pyro5.api.Daemon(host='192.168.3.69', port=9091)
uri = daemon.register(Job, objectId="HPC-Job")
print("HPC starting ...")
daemon.requestLoop()
