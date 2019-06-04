#!/usr/bin/env python
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
#
# Change shebang line to '/usr/bin/python -3' for python 3.x porting warnings
"""
This script drives a PR testing build.  It assume that Trilinos is already
cloned under $PWD/Trilinos and that the 'origin' remote points to
$TRILINOS_TARGET_REPO (but that is not checked here).

As long as the ${PWD}/Trilinos git repo has the correct 'origin', this
script will automatically set it up to do the merge correctly, no matter
what its state before this script is called (i.e. from a past PR
attempt). Unless the Trilinos/.git directory becomes corrupted, there should
*NEVER* be any need to delete and reclone this Trilinos git repo.

This script can be run in a mode where the driver scripts are run from one
Trilinos git repo and operate on another Trilinos git repo that gets
manipulated in order to merge the "source" topic branch into the "target"
branch.  This makes it easy to test changes to the PR scripts.  But if this
script is run from ${PWD}/Trilinos, then these repos are one and the same
and we get the correct behavior for PR testing.
"""
from __future__ import print_function

# turn off generation of the .pyc/.pyo files.
import sys
sys.dont_write_bytecode = True


def run():
    return_value = True

    return return_value


if __name__ == '__main__':  # pragma nocover
    returnValue = run()
    if returnValue:
        exit(0)
    else:
        exit(1)
