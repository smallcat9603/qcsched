import Pyro5.api
import time
import platform

class Jobs(object):

    @Pyro5.server.expose
    @Pyro5.server.oneway
    def run(self, vid, t):
        time.sleep(t)
        
        # return "HPC-Job finished."

    
def main():
    uname = platform.uname()
    print(uname)
    system = uname[0]
    node = uname[1]
    release = uname[2]

    host = '192.168.3.69' if 'raspberrypi' in node else None
    port = 9091

    daemon = Pyro5.api.Daemon(host=host, port=port)
    uri = daemon.register(Jobs, objectId="HPC-Job")
    print("HPC running ...")
    daemon.requestLoop()

if __name__ == "__main__":
    main()
