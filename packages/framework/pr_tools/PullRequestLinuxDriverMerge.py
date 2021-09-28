#!/usr/bin/env python3 -u
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
# Change above to '/usr/bin/python -3' for python 3.x porting warnings
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
import sys
sys.dont_write_bytecode = True

import argparse
import os
import subprocess
import sys
from textwrap import dedent


def print_wrapper(text: str, prefix="PRLinuxDriverMerge> ", end="\n"):
    """
    """
    rval = print(f"{prefix}{text}", end=end)
    #sys.stdout.flush()
    return rval


def write_header():
    """
    """
    print_wrapper("-"*80)
    print_wrapper("-")
    print_wrapper("- Begin: PullRequestLinuxDriverMerge.py")
    print("-")
    print_wrapper("-"*80)


def echoJenkinsVars(workspace):
    """
    """
    print_wrapper(80*"-")
    print_wrapper(f"Jenkins Environment Variables:")
    print_wrapper(f"WORKSPACE: {workspace}")
    print_wrapper(f"Environment:")
    print_wrapper(f"-- pwd = {os.getcwd()}")
    print_wrapper(f"Environment Variables:")
    for key in os.environ:
        print_wrapper(f"-- {key} = {os.environ[key]}")
    print_wrapper("")
    print_wrapper(80*"-")


def parseArgs():
    """
    Parse the arguments - no  options are available at this time
    """
    parser = argparse.ArgumentParser(description='Parse the repo and merge information')
    parser.add_argument('sourceRepo',
                        help='Repo with the new changes',
                        action='store')
    parser.add_argument('sourceBranch',
                        help='Branch with the new changes',
                        action='store')
    parser.add_argument('targetRepo',
                        help='Repo to merge into',
                        action='store')
    parser.add_argument('targetBranch',
                        help='Branch to merge to',
                        action='store')
    parser.add_argument('sourceSHA',
                        help='SHA1 of the commit to use from the source branch',
                        action='store')
    parser.add_argument('workspaceDir',
                        help='The local workspace directory jenkins set up')

    return parser.parse_args()


def check_call_wrapper(args):
    """
    A simple wrapper for subprocess.check_call() that prints out
    a more verbose bit of output to stdout for console logging, etc.

    Args:
        args (list): A list of arguments to be executed.

    Returns:
        int: Returns a 0.
    """
    # print("PRLinuxDriverMerge> {}".format(" ".join(args)))
    print_wrapper("Checked Call:")
    print_wrapper(" ".join( [str(x) for x in args] ))
    subprocess.check_call(args)
    sys.stdout.flush()
    sys.stderr.flush()
    print_wrapper("OK")
    print_wrapper("")
    return None


def check_output_wrapper(args):
    """
    A simple wrapper for ``subprocess.check_output()`` that prints out
    a more verbose bit of output to stdout for console logging, etc.

    Args:
        args (list): A list of arguments to be executed.

    Returns:
        int: Passes along the value returned by ``subprocess.check_output()``
    """
    print_wrapper(" ".join(args))
    output = subprocess.check_output(args)
    print_wrapper("")
    sys.stdout.flush()
    sys.stderr.flush()
    return output


def merge_branch(source_url, source_branch, target_branch, sourceSHA):
    """
    TODO: add docstring.
    """
    source_url    = source_url.strip()
    source_branch = source_branch.strip()
    target_branch = target_branch.strip()
    sourceSHA     = sourceSHA.strip()

    remote_list = check_output_wrapper(['git', 'remote', '-v'])

    if isinstance(remote_list, bytes):
        remote_list = remote_list.decode('utf-8')

    if 'source_remote' in remote_list:
        print_wrapper("git remote exists, removing it")
        check_call_wrapper(['git', 'remote', 'rm', 'source_remote'])

    check_call_wrapper(['git', 'remote', 'add', 'source_remote', source_url])

    fetch_succeeded = False
    for i in range(3):
        try:
            check_call_wrapper(['git', 'fetch', 'source_remote', source_branch])
            fetch_succeeded = True
            break
        except subprocess.CalledProcessError:
            pass

    if not fetch_succeeded:
        print_wrapper("ERROR: Fetch did not succeed.")
        raise SystemExit(12)

    check_call_wrapper(['git', 'fetch', 'origin', target_branch])
    check_call_wrapper(['git', 'reset', '--hard', 'HEAD'])
    check_call_wrapper(['git', 'checkout', '-B', target_branch, 'origin/' + target_branch])
    check_call_wrapper(['git', 'merge', '--no-edit', 'source_remote/' + source_branch]),

    actual_source_SHA = check_output_wrapper(['git', 'rev-parse', 'source_remote/' + source_branch])

    if isinstance(actual_source_SHA, bytes):
        actual_source_SHA = actual_source_SHA.decode('utf-8')

    if isinstance(actual_source_SHA, bytes):
        actual_source_SHA = actual_source_SHA.decode('utf-8')

    actual_source_SHA = actual_source_SHA.strip()

    if actual_source_SHA != sourceSHA:
        print_wrapper(f"The SHA ({actual_source_SHA}) for the last commit on branch {source_branch}")
        print_wrapper(f"in repo {source_url} is different from the expected SHA,")
        print_wrapper(f"which is: {sourceSHA}.")
        raise SystemExit(-1)

    return 0



def run():
    return_value = True
    try:
        arguments = parseArgs()
    except SystemExit:
        return_value = False
    if return_value:
        os.chdir(os.path.join(arguments.workspaceDir, 'Trilinos'))
        print_wrapper(f"Set CWD = {os.path.join(arguments.workspaceDir, 'Trilinos')}")
        write_header()
        echoJenkinsVars(arguments.workspaceDir)
        try:
            merge_branch(arguments.sourceRepo,
                         arguments.sourceBranch,
                         arguments.targetBranch,
                         arguments.sourceSHA)
        except SystemExit:
            return_value = False
        except subprocess.CalledProcessError as cpe:
            return_value = False

            print_wrapper(f"Received subprocess.CalledProcessError - returned {cpe.returncode}")
            print_wrapper(f"from command {cpe.cmd}")
            print_wrapper(f"output {cpe.output}")

            try:
                print_wrapper(f"  stdout {cpe.stdout}")
            except AttributeError:
                pass
            try:
                print_wrapper(f"  stderr {cpe.stderr}")
            except AttributeError:
                pass

    return return_value


if __name__ == '__main__':  # pragma nocover
    returnValue = run()
    if returnValue:
        print_wrapper(f" Finished Normally")
        exit(0)
    else:
        print_wrapper(f" Finished with error(s)")
        exit(1)


