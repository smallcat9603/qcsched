import Pyro5.api
import argparse
from setting import *


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('vids', nargs='+', help='subjob vids to be deleted')
    args = parser.parse_args()

    server_sched = Pyro5.client.Proxy(URI_SCHED)

    print(server_sched.del_subjobs(args.vids))


if __name__ == "__main__":
    main()
