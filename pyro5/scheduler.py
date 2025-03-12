import Pyro5.api
import platform

uname = platform.uname()
print(uname)
system = uname[0]
node = uname[1]
release = uname[2]

host = 'localhost'
port = 9091
if 'raspberry' in release:
    host = '192.168.3.69'

uri = f"PYRO:HPC-Job@{host}:{port}"
server = Pyro5.client.Proxy(uri)
t = int(input("The time of the HPC job? ").strip())
print(server.run(t))

host = 'localhost'
port = 9092
if 'raspberry' in release:
    host = '192.168.3.80'

uri = f"PYRO:QC-Job@{host}:{port}"
server = Pyro5.client.Proxy(uri)
t = int(input("The time of the QC job? ").strip())
print(server.run(t))
