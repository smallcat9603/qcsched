import Pyro5.api
import time
import threading
from setting import *


class HPC():
    def __init__(self, name: str):
        self.name = name
        self.events = {}


    def run_threaded(self, vid: int, t: int, event):
        for i in range(t):
            time.sleep(1)
            event.wait()
        server_sched = Pyro5.client.Proxy(URI_SCHED)
        server_sched.finish(vid)


    @Pyro5.server.expose
    @Pyro5.server.oneway
    def run(self, vid: int, t: int):  
        self.events[vid] = threading.Event()      
        thread = threading.Thread(target=self.run_threaded, args=[vid, t, self.events[vid]])
        thread.start()
        self.events[vid].set()


    @Pyro5.server.expose
    @Pyro5.server.oneway
    def pause(self, vid: int):
        self.events[vid].clear()


    @Pyro5.server.expose
    @Pyro5.server.oneway
    def resume(self, vid: int):
        self.events[vid].set()

    
def main():

    daemon = Pyro5.api.Daemon(host=HOST_HPC, port=PORT_HPC)
    hpc = HPC(NAME_HPC)
    daemon.register(hpc, objectId=NAME_HPC)
    print("HPC running ...")
    daemon.requestLoop()


if __name__ == "__main__":
    main()
