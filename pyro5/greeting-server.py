# saved as greeting-server.py
import Pyro5.api

@Pyro5.api.expose
class GreetingMaker(object):
    def get_fortune(self, name):
        return "Hello, {0}. Here is your fortune message:\n" \
               "Behold the warranty -- the bold print giveth and the fine print taketh away.".format(name)

daemon = Pyro5.api.Daemon(host='192.168.3.13')             # make a Pyro daemon
uri = daemon.register(GreetingMaker)    # register the greeting maker as a Pyro object

print("Ready. Object uri =", uri)       # print the uri so we can use it in the client later
daemon.requestLoop()                    # start the event loop of the server to wait for calls
