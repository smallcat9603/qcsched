import Pyro5.api
import sys
from setting import *


def main():
    if len(sys.argv) != 2:
        print("Usage: python qcsub.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]  

    server_sched = Pyro5.client.Proxy(URI_SCHED)

    print(server_sched.submit_joblist(filename))


if __name__ == "__main__":
    main()
