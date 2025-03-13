import Pyro5.api
import time
import platform

class Jobs(object):

    @Pyro5.server.expose
    @Pyro5.server.oneway
    def run(self, vid, t):
        time.sleep(t)
        
        # return "QC-Job finished."

    
def main():
    uname = platform.uname()
    print(uname)
    system = uname[0]
    node = uname[1]
    release = uname[2]

    host = '192.168.3.80' if 'raspberrypi' in node else None
    port = 9092

    daemon = Pyro5.api.Daemon(host=host, port=port)
    uri = daemon.register(Jobs, objectId="QC-Job")
    print("QC running ...")
    daemon.requestLoop()

if __name__ == "__main__":
    main()
