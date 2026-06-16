"""
wait_for_request.py — continuously watches for reconfiguration requests
(s<N>_request.json, written by Teko::KrylovSurrogate::writeRequestJson) and
answers each one by writing s<N>_reconfig.json (read by
Teko::KrylovSurrogate::waitForOrdering).

The C++ side leaves s<N>_reconfig.json in place after reading it (it does not
delete it), and writes s<N>_conv.json once it has consumed the ordering. This
watcher simply keeps looping: it answers each new s<N>_request.json once and
goes back to watching for the next request.

Runs until Ctrl-C or until RUN_SECONDS (default 300 = 5 minutes) have
elapsed, whichever comes first.

Run with:
    source /home/node/codespace/mypy/bin/activate
    python /home/node/codespace/Trilinos/teko-reconfig/wait_for_request.py

All three file types (s<N>_request.json, s<N>_reconfig.json, s<N>_conv.json)
live together in one directory. Defaults to the "requests" sibling of this
script; override with TEKO_RECONFIG_REQUESTS_DIR to match the C++ side's
default (kDefaultRequestsDir in Teko_KrylovSurrogate.hpp).
"""

import glob
import json
import os
import re
import time

REQUESTS_DIR = os.environ.get(
    "TEKO_RECONFIG_REQUESTS_DIR",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "requests"),
)
POLL_INTERVAL_SEC = 0.5
RUN_SECONDS = 300

REQUEST_RE = re.compile(r"s(\d+)_request\.json$")


def already_handled(request_id):
    # s<N>_conv.json is written (atomically) by the C++ side once it has
    # consumed s<N>_reconfig.json (or given up waiting for it). Its presence
    # means s<N>_request.json is a stale leftover from a previous run's
    # pyTeko() call, not a request this watcher should answer.
    return os.path.exists(os.path.join(REQUESTS_DIR, f"s{request_id}_conv.json"))


def make_ordering(n_blocks):
    if n_blocks == 0:
        raise ValueError("n_blocks must be greater than zero")
    if n_blocks == 1:
        return [0]
    return [0] * (n_blocks - 1) + [1]


# selection_mode written into the reconfig response, switching what the C++
# side returns after timing the test_orderings sweep:
#   "chosen"    — apply use_ordering (the sweep is still timed & recorded).
#   "best_conv" — apply the swept ordering with the best convergence.
#   "best_time" — apply the swept ordering with the best (factor+iterate) time.
SELECTION_MODE = "chosen"

# When True, emit a test_orderings sweep for the C++ side to build/solve/time
# (results land in s<N>_solved.json). For small block counts this is the full
# set of partitions (every way to group the blocks); above the cap only a few
# representative candidates are sent to keep the sweep bounded.
EMIT_TEST_ORDERINGS = True
MAX_BLOCKS_FOR_EXHAUSTIVE = 6


def restricted_growth_strings(n):
    """Every set partition of n blocks, as a restricted-growth string. These
    are exactly the valid orderings (contiguous group ids starting at 0). The
    count is the Bell number B(n)."""
    if n <= 0:
        return [[]]
    result = []

    def rec(prefix, max_used):
        if len(prefix) == n:
            result.append(list(prefix))
            return
        for v in range(max_used + 2):
            prefix.append(v)
            rec(prefix, max(max_used, v))
            prefix.pop()

    rec([], -1)
    return result


def make_test_orderings(n_blocks, use_ordering):
    if not EMIT_TEST_ORDERINGS or n_blocks == 0:
        return []
    if n_blocks <= MAX_BLOCKS_FOR_EXHAUSTIVE:
        return restricted_growth_strings(n_blocks)
    # Too many partitions to sweep exhaustively: send a few representatives.
    candidates = [list(use_ordering), list(range(n_blocks)), [0] * n_blocks]
    seen, out = set(), []
    for c in candidates:
        key = tuple(c)
        if key not in seen:
            seen.add(key)
            out.append(c)
    return out


def write_reconfig(ordering, path, n_blocks):
    # Write to a temp file and atomically rename into place, so the C++ side's
    # fs::exists(path) poll never observes a partially-written file.
    #
    # Fields:
    #   selection_mode   — what the C++ side returns after the sweep (see
    #                      SELECTION_MODE). Written first.
    #   use_ordering     — the ordering applied in "chosen" mode (the field C++
    #                      always reads).
    #   opt_ordering     — the "optimal" ordering from the cheap surrogate
    #                      search. For now a copy of use_ordering.
    #   exh_opt_ordering — the optimal ordering from an exhaustive search.
    #                      Left blank for now.
    #   test_orderings   — candidate orderings to build/solve/time; results are
    #                      written by C++ to s<N>_solved.json. May be empty.
    tmp_path = path + ".tmp"
    payload = {
        "selection_mode": SELECTION_MODE,
        "use_ordering": ordering,
        "opt_ordering": list(ordering),
        "exh_opt_ordering": [],
        "test_orderings": make_test_orderings(n_blocks, ordering),
    }
    with open(tmp_path, "w") as f:
        json.dump(payload, f, indent=2)
    os.replace(tmp_path, path)


def handle_request(request_id):
    request_path = os.path.join(REQUESTS_DIR, f"s{request_id}_request.json")
    reconfig_path = os.path.join(REQUESTS_DIR, f"s{request_id}_reconfig.json")

    with open(request_path) as f:
        request = json.load(f)

    c_hat = request["C_hat"]            # single R x R matrix (list of rows)
    b_hat = request["b_hat"]            # single length-R vector
    n_rows = len(c_hat)
    n_cols = len(c_hat[0]) if c_hat else 0

    print(f"Received s{request_id}_request.json:")
    print(f"  n_blocks      = {request['n_blocks']}")
    print(f"  block_sizes   = {request['block_sizes']}")
    print(f"  krylov_dim    = {request['krylov_dim']}")
    print(f"  ranks         = {request['ranks']}")
    print(f"  equation_ends = {request['equation_ends']}")
    print(f"  C_hat shape   = {n_rows} x {n_cols}")
    print(f"  b_hat length  = {len(b_hat)}")

    ordering = make_ordering(request["n_blocks"])
    write_reconfig(ordering, reconfig_path, request["n_blocks"])
    n_tests = len(make_test_orderings(request["n_blocks"], ordering))
    print(f"Wrote reconfig (mode={SELECTION_MODE}, use_ordering={ordering}, "
          f"{n_tests} test orderings) to {reconfig_path}")


def main():
    print(f"Watching {REQUESTS_DIR} for s<N>_request.json "
          f"(until Ctrl-C or {RUN_SECONDS}s)...")

    answered = set()
    start = time.monotonic()
    try:
        while time.monotonic() - start < RUN_SECONDS:
            for path in glob.glob(os.path.join(REQUESTS_DIR, "s*_request.json")):
                m = REQUEST_RE.match(os.path.basename(path))
                if not m:
                    continue
                request_id = int(m.group(1))
                if request_id in answered:
                    continue
                if already_handled(request_id):
                    answered.add(request_id)
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
