import Pyro5.api
import time
import threading
from setting import *


class HPC():
    def __init__(self, uri_sched: str):
        self.uri_sched = uri_sched

    def run_threaded(self, vid: int, t: int):
        time.sleep(t)
        server_sched = Pyro5.client.Proxy(self.uri_sched)
        server_sched.finish(vid)

    @Pyro5.server.expose
    @Pyro5.server.oneway
    def run(self, vid: int, t: int):
        thread = threading.Thread(target=self.run_threaded, args=[vid, t])
        thread.start()
        
    
def main():

    daemon = Pyro5.api.Daemon(host=HOST_HPC, port=PORT_HPC)
    hpc = HPC(URI_SCHED)
    uri = daemon.register(hpc, objectId="hpc")
    print("HPC running ...")
    daemon.requestLoop()

if __name__ == "__main__":
    main()
