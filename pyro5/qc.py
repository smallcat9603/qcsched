import Pyro5.api
import time

@Pyro5.api.expose
class Job(object):
    def run(self, t):
        time.sleep(t)
        return "QC-Job finished."

daemon = Pyro5.api.Daemon(host='192.168.3.80', port=9092)
uri = daemon.register(Job, objectId="QC-Job")
print("QC starting ...")
daemon.requestLoop()
