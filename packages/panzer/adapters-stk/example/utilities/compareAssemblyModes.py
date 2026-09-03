#!/usr/bin/env python3
"""Compare scalar values from two Panzer runs that differ only in assembly path.

The FE assembly path (Tpetra::FECrsMatrix, one object with owned/owned+shared views
migrated by endAssembly) and the classic path (separate owned/ghosted matrices joined
by an explicit export) solve the identical problem, so the reported scalars must agree.

Values are parsed as doubles and compared with a RELATIVE tolerance rather than by
string equality: identical text would be a stricter check, but it makes the test
sensitive to last-digit roundoff differences across compilers, MPI rank counts and
machines, which is instability unrelated to what is being verified.

Usage: compareAssemblyModes.py <fe_output> <classic_output> <spec> [<spec> ...]

Each spec is either

    LABEL             compare the two runs against each other only
    LABEL=EXPECTED    additionally require BOTH runs to match EXPECTED

The second form is the stronger check: agreeing with each other only shows the two
paths are consistent, not that either is right, so where an analytic answer exists it
is compared against directly.

A label selects the lines containing it that carry a number after the label; the last
such number is the value compared. A label must select the same number of values in
both files, and at least one.
"""

import re
import sys

RELATIVE_TOLERANCE = 1.0e-12
ANALYTIC_RELATIVE_TOLERANCE = 1.0e-8

# matches C/C++ style floating point output, including exponents
NUMBER = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")

# e.g. "p=0 | Response Value ..." -> "Response Value ..."
RANK_PREFIX = re.compile(r"^p=\d+\s*\|\s*")


def values_for_label(path, label):
    """Last number appearing AFTER `label` on each line containing it.

    Searching only the remainder of the line matters: labels themselves may contain
    digits ("L2 Error", "H1 Error"), so scanning the whole line would pick the "2" out
    of a heading like "This is the L2 Error" and report a spurious value. Lines with
    nothing numeric after the label are headings and are skipped.
    """
    found = []
    with open(path) as f:
        for line in f:
            line = RANK_PREFIX.sub("", line.strip())
            pos = line.find(label)
            if pos < 0:
                continue
            numbers = NUMBER.findall(line[pos + len(label):])
            if numbers:
                found.append(float(numbers[-1]))
    return found


def relative_difference(a, b):
    denom = max(abs(a), abs(b))
    if denom == 0.0:
        return 0.0
    return abs(a - b) / denom


def main():
    if len(sys.argv) < 4:
        print("Usage: compareAssemblyModes.py <fe_output> <classic_output> "
              "<label> [<label> ...]")
        return 1

    fe_path, classic_path = sys.argv[1], sys.argv[2]
    specs = sys.argv[3:]

    failed = False
    for spec in specs:
        if "=" in spec:
            label, expected_text = spec.rsplit("=", 1)
            expected = float(expected_text)
        else:
            label, expected = spec, None

        fe = values_for_label(fe_path, label)
        classic = values_for_label(classic_path, label)

        if not fe or not classic:
            print("Test Failed: label '%s' selected no values "
                  "(fe: %d, classic: %d) -- did both runs complete?"
                  % (label, len(fe), len(classic)))
            failed = True
            continue

        if len(fe) != len(classic):
            print("Test Failed: label '%s' selected %d value(s) from the FE run but "
                  "%d from the classic run" % (label, len(fe), len(classic)))
            failed = True
            continue

        for i, (a, b) in enumerate(zip(fe, classic)):
            suffix = "" if len(fe) == 1 else " [%d]" % i

            rel = relative_difference(a, b)
            status = "ok" if rel <= RELATIVE_TOLERANCE else "FAILED"
            if rel > RELATIVE_TOLERANCE:
                failed = True
            print("  %-28s fe=%.17g classic=%.17g rel_diff=%.3g  %s"
                  % (label + suffix, a, b, rel, status))

            if expected is not None:
                for name, value in (("fe", a), ("classic", b)):
                    rel_e = relative_difference(value, expected)
                    status_e = "ok" if rel_e <= ANALYTIC_RELATIVE_TOLERANCE else "FAILED"
                    if rel_e > ANALYTIC_RELATIVE_TOLERANCE:
                        failed = True
                    print("  %-28s %s=%.17g expected=%.17g rel_diff=%.3g  %s"
                          % (label + suffix + " vs analytic", name, value,
                             expected, rel_e, status_e))

    if failed:
        print("Test Failed: checks did not pass (run-to-run tolerance %g, "
              "analytic tolerance %g)" % (RELATIVE_TOLERANCE,
                                          ANALYTIC_RELATIVE_TOLERANCE))
        return 1

    print("All checks passed (run-to-run tolerance %g, analytic tolerance %g)"
          % (RELATIVE_TOLERANCE, ANALYTIC_RELATIVE_TOLERANCE))
    print("Test Passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
