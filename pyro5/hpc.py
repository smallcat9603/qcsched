import Pyro5.api
import time
import threading
from setting import *


class HPC():
    def __init__(self, name: str):
        self.name = name

    def run_threaded(self, vid: int, t: int):
        time.sleep(t)
        server_sched = Pyro5.client.Proxy(URI_SCHED)
        server_sched.finish(vid)

    @Pyro5.server.expose
    @Pyro5.server.oneway
    def run(self, vid: int, t: int):
        thread = threading.Thread(target=self.run_threaded, args=[vid, t])
        thread.start()
        
    
def main():

    daemon = Pyro5.api.Daemon(host=HOST_HPC, port=PORT_HPC)
    hpc = HPC(NAME_HPC)
    daemon.register(hpc, objectId=NAME_HPC)
    print("HPC running ...")
    daemon.requestLoop()


if __name__ == "__main__":
    main()
