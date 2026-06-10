"""
wait_for_request.py — continuously watches for reconfiguration requests
(request_<N>.json, written by Teko::KrylovSurrogate::writeRequestJson) and
answers each one by writing reconfiguration_<N>.json (read by
Teko::KrylovSurrogate::waitForOrdering).

Runs until Ctrl-C or until RUN_SECONDS (default 300 = 5 minutes) have
elapsed, whichever comes first.

Run with:
    source /workspace/mypy/bin/activate
    python /workspace/Trilinos/teko-reconfig/wait_for_request.py

Defaults to watching this directory; override with TEKO_RECONFIG_REQUESTS_DIR
to match the C++ side's default (kDefaultRequestsDir in
Teko_KrylovSurrogate.hpp).
"""

import glob
import json
import os
import re
import time

REQUESTS_DIR = os.environ.get(
    "TEKO_RECONFIG_REQUESTS_DIR",
    os.path.dirname(os.path.abspath(__file__)),
)
POLL_INTERVAL_SEC = 0.5
RUN_SECONDS = 300

REQUEST_RE = re.compile(r"request_(\d+)\.json$")


def make_ordering(n_blocks):
    if n_blocks == 0:
        raise ValueError("n_blocks must be greater than zero")
    if n_blocks == 1:
        return [0]
    return [0] * (n_blocks - 1) + [1]


def write_reconfig(ordering, path):
    # Write to a temp file and atomically rename into place, so the C++ side's
    # fs::exists(path) poll never observes a partially-written file.
    tmp_path = path + ".tmp"
    with open(tmp_path, "w") as f:
        json.dump({"ordering": ordering}, f, indent=2)
    os.replace(tmp_path, path)


def handle_request(request_id):
    request_path = os.path.join(REQUESTS_DIR, f"request_{request_id}.json")
    reconfig_path = os.path.join(REQUESTS_DIR, f"reconfiguration_{request_id}.json")

    with open(request_path) as f:
        request = json.load(f)

    print(f"Received request_{request_id}.json:")
    print(f"  n_blocks   = {request['n_blocks']}")
    print(f"  block_sizes = {request['block_sizes']}")
    print(f"  krylov_dim = {request['krylov_dim']}")
    print(f"  ranks      = {request['ranks']}")
    print(f"  C_hat keys = {sorted(request['C_hat'].keys())}")
    print(f"  b_hat keys = {sorted(request['b_hat'].keys())}")

    ordering = make_ordering(request["n_blocks"])
    write_reconfig(ordering, reconfig_path)
    print(f"Wrote ordering {ordering} to {reconfig_path}")


def main():
    print(f"Watching {REQUESTS_DIR} for request_<N>.json "
          f"(until Ctrl-C or {RUN_SECONDS}s)...")

    answered = set()
    start = time.monotonic()
    try:
        while time.monotonic() - start < RUN_SECONDS:
            for path in glob.glob(os.path.join(REQUESTS_DIR, "request_*.json")):
                m = REQUEST_RE.match(os.path.basename(path))
                if not m:
                    continue
                request_id = int(m.group(1))
                if request_id in answered:
                    continue
                handle_request(request_id)
                answered.add(request_id)
            time.sleep(POLL_INTERVAL_SEC)
    except KeyboardInterrupt:
        print("Interrupted.")

    elapsed = time.monotonic() - start
    print(f"Stopping after {elapsed:.1f}s; answered requests: {sorted(answered)}")


if __name__ == "__main__":
    main()
