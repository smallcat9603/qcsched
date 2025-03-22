import Pyro5.api
import time
import threading
from setting import *


class QC():
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

    daemon = Pyro5.api.Daemon(host=HOST_QC, port=PORT_QC)
    qc = QC(NAME_QC)
    daemon.register(qc, objectId=NAME_QC)
    print("QC running ...")
    daemon.requestLoop()


if __name__ == "__main__":
    main()
