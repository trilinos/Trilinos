import argparse
import sys

try:
    parser = argparse.ArgumentParser()
    parser.add_argument("--merged", type=int, required=True, help="first argument")
    parser.add_argument("--failed", type=int, required=True, help="second argument")
    parser.add_argument("--wip", type=int, required=True, help="third argument")
    parser.add_argument("--reviewed", type=int, required=True, help="fourth argument")
    parser.add_argument("--review", type=int, required=True, help="fifth argument")
    parser.add_argument("--waiting", type=int, required=True, help="sixth argument")
    parser.add_argument("--open", type=int, required=True, help="seventh argument")
    parser.add_argument("--master", type=int, required=True, help="eighth argument")
    parser.add_argument("--jira", type=int, required=True, help="ninth argument")

    args = parser.parse_args()

    print("Argument 1:", args.merged)
    print("Argument 2:", args.failed)
    print("Argument 3:", args.wip)
    print("Argument 4:", args.reviewed)
    print("Argument 5:", args.review)
    print("Argument 6:", args.waiting)
    print("Argument 7:", args.open)
    print("Argument 8:", args.master)
    print("Argument 9:", args.jira)

except:
    print(sys.stderr)


