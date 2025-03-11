# saved as greeting-client.py
import Pyro5.api

uri = "PYRO:HPC-Job@192.168.3.69:9091"
server = Pyro5.client.Proxy(uri)
t = int(input("The time of the HPC job? ").strip())
print(server.run(t))

uri = "PYRO:QC-Job@192.168.3.80:9092"
server = Pyro5.client.Proxy(uri)
t = int(input("The time of the QC job? ").strip())
print(server.run(t))
