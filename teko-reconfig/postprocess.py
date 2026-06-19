#!/usr/bin/env python3
"""Post-process Teko adaptive reconfiguration convergence files."""

import argparse
import json
import math
import re
from pathlib import Path
import os
import sys


def _ensure_mypy_venv():
    script_dir = Path(__file__).resolve().parent
    venv_dir = script_dir.parents[1] / "trilinos-teko-pyfront" / "mypy"
    venv_python = venv_dir / "bin" / "python3"

    if Path(sys.prefix).resolve() == venv_dir.resolve():
        return
    if not venv_python.exists():
        raise RuntimeError(f"Expected Python venv at {venv_python}")

    os.environ["VIRTUAL_ENV"] = str(venv_dir)
    os.environ["PATH"] = str(venv_dir / "bin") + os.pathsep + os.environ.get("PATH", "")
    os.environ.pop("PYTHONHOME", None)
    os.execv(str(venv_python), [str(venv_python), *sys.argv])


_ensure_mypy_venv()

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


CONV_RE = re.compile(r"s(\d+)_conv\.json$")


def load_convergence(requests_dir):
    rows = []
    for path in sorted(requests_dir.glob("s*_conv.json")):
        match = CONV_RE.match(path.name)
        if match is None:
            continue

        with path.open() as f:
            data = json.load(f)

        request_id = int(data.get("request_id", match.group(1)))
        solve1 = data.get("solve1") or {}
        solve2 = data.get("solve2") or {}
        rows.append(
            {
                "request_id": request_id,
                "solve1_iterations": solve1.get("iterations", math.nan),
                "solve2_iterations": solve2.get("iterations", math.nan),
                "solve1_wall_time_sec": solve1.get(
                    "total_wall_time_sec", solve1.get("wall_time_sec", math.nan)
                ),
                "solve2_wall_time_sec": solve2.get(
                    "total_wall_time_sec", solve2.get("wall_time_sec", math.nan)
                ),
            }
        )

    return sorted(rows, key=lambda row: row["request_id"])


def plot_convergence(rows, output_path):
    if not rows:
        raise RuntimeError("No s<N>_conv.json files found")

    request_ids = [row["request_id"] for row in rows]
    solve1_iters = [row["solve1_iterations"] for row in rows]
    solve2_iters = [row["solve2_iterations"] for row in rows]
    solve1_times = [row["solve1_wall_time_sec"] for row in rows]
    solve2_times = [row["solve2_wall_time_sec"] for row in rows]

    fig, (ax_iters, ax_time) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    ax_iters.plot(request_ids, solve1_iters, marker="o", label="Solve 1")
    ax_iters.plot(request_ids, solve2_iters, marker="s", label="Solve 2")
    ax_iters.set_ylabel("Iterations")
    ax_iters.set_title("Teko Adaptive Reconfiguration")
    ax_iters.grid(True, alpha=0.3)
    ax_iters.legend()

    ax_time.plot(request_ids, solve1_times, marker="o", label="Solve 1")
    ax_time.plot(request_ids, solve2_times, marker="s", label="Solve 2")
    ax_time.set_xlabel("Request number")
    ax_time.set_ylabel("Total wall time (s)")
    ax_time.grid(True, alpha=0.3)
    ax_time.legend()

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main():
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--requests-dir",
        type=Path,
        default=script_dir / "requests",
        help="Directory containing s<N>_conv.json files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=script_dir / "figs" / "reconfig_convergence.png",
        help="Output figure path",
    )
    args = parser.parse_args()

    rows = load_convergence(args.requests_dir)
    plot_convergence(rows, args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
