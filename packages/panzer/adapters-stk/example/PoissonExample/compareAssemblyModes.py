#!/usr/bin/env python3
"""Compare PoissonExample output from two assembly modes.

The FE assembly path (Tpetra::FECrsMatrix, one object with owned/owned+shared
views migrated by endAssembly) and the classic path (separate owned/ghosted
matrices joined by an explicit export) solve the identical problem, so every
reported error norm must match exactly -- not merely to a tolerance. The norms
are printed to full precision, so comparing the text of the "... Error = ..."
lines is a strict check that the two paths agree numerically.

Usage: compareAssemblyModes.py <fe_output> <classic_output>
"""

import sys


def error_lines(path):
    with open(path) as f:
        return [line.strip() for line in f if "Error =" in line]


def main():
    if len(sys.argv) != 3:
        print("Usage: compareAssemblyModes.py <fe_output> <classic_output>")
        return 1

    fe_path, classic_path = sys.argv[1], sys.argv[2]

    fe = error_lines(fe_path)
    classic = error_lines(classic_path)

    if not fe or not classic:
        print("Test Failed: no 'Error =' lines found "
              "(fe: %d, classic: %d) -- did both runs complete?" % (len(fe), len(classic)))
        return 1

    if fe != classic:
        print("Test Failed: FE assembly does not match classic assembly")
        print("  FE      (%s):" % fe_path)
        for line in fe:
            print("    " + line)
        print("  classic (%s):" % classic_path)
        for line in classic:
            print("    " + line)
        return 1

    print("FE assembly matches classic assembly on %d error norms:" % len(fe))
    for line in fe:
        print("  " + line)
    print("Test Passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
