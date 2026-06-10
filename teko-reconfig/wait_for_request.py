"""
wait_for_request.py — waits for a reconfiguration request to appear at
request.json (written by Teko::KrylovSurrogate::writeRequestJson) and reads
it, then writes reconfiguration.json (read by
Teko::KrylovSurrogate::waitForOrdering).

Run with:
    source /workspace/mypy/bin/activate
    python /workspace/Trilinos/teko-reconfig/wait_for_request.py

Set TEKO_RECONFIG_REQUESTS_DIR to this directory before running the solve.
"""

import json
import os
import time

REQUESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REQUEST_PATH = os.path.join(REQUESTS_DIR, "request.json")
RECONFIG_PATH = os.path.join(REQUESTS_DIR, "reconfiguration.json")
POLL_INTERVAL_SEC = 0.5


def wait_for_request(path=REQUEST_PATH, poll_interval=POLL_INTERVAL_SEC):
    print(f"Waiting for request at {path} ...")
    while not os.path.exists(path):
        time.sleep(poll_interval)

    with open(path) as f:
        request = json.load(f)

    return request


def make_ordering(n_blocks):
    if n_blocks == 0:
        raise ValueError("n_blocks must be greater than zero")
    if n_blocks == 1:
        return [0]
    return [0] * (n_blocks - 1) + [1]


def write_reconfig(ordering, path=RECONFIG_PATH):
    with open(path, "w") as f:
        json.dump({"ordering": ordering}, f, indent=2)


if __name__ == "__main__":
    request = wait_for_request()

    print("Received request:")
    print(f"  n_blocks   = {request['n_blocks']}")
    print(f"  block_sizes = {request['block_sizes']}")
    print(f"  krylov_dim = {request['krylov_dim']}")
    print(f"  ranks      = {request['ranks']}")
    print(f"  C_hat keys = {sorted(request['C_hat'].keys())}")
    print(f"  b_hat keys = {sorted(request['b_hat'].keys())}")

    ordering = make_ordering(request["n_blocks"])
    write_reconfig(ordering)
    print(f"Wrote ordering {ordering} to {RECONFIG_PATH}")
