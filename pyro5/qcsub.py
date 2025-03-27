import Pyro5.api
import argparse
from setting import *


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='job list file')
    parser.add_argument('-add', type=int, help='added job vid')
    args = parser.parse_args()

    server_sched = Pyro5.client.Proxy(URI_SCHED)

    if args.add:
        print(server_sched.add_subjob(args.add, args.file))
    else:
        print(server_sched.submit_joblist(args.file))


if __name__ == "__main__":
    main()
