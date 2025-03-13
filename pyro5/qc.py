import Pyro5.api
import time
import platform
import threading

class QC():
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
    uname = platform.uname()
    print(uname)
    system = uname[0]
    node = uname[1]
    release = uname[2]

    host_sched = '192.168.3.13' if 'raspberrypi' in node else 'localhost'
    port_sched = 9093
    uri_sched = f"PYRO:sched@{host_sched}:{port_sched}"

    host = '192.168.3.80' if 'raspberrypi' in node else None
    port = 9092

    daemon = Pyro5.api.Daemon(host=host, port=port)
    qc = QC(uri_sched)
    uri = daemon.register(qc, objectId="qc")
    print("QC running ...")
    daemon.requestLoop()

if __name__ == "__main__":
    main()
