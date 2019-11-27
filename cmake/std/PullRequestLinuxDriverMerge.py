#!/usr/bin/env python
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

import os
import argparse
import subprocess


def write_header():
    print('''--------------------------------------------------------------------------------
-
- Begin: PullRequestLinuxDriver-Merge.py
-
--------------------------------------------------------------------------------''',
          file=sys.stdout)

def echoJenkinsVars(workspace):
    print('''
================================================================================
Jenkins Environment Variables:
- WORKSPACE    : {workspace}

================================================================================
Environment:

  pwd = {cwd}
'''.format(workspace=workspace,
           cwd=os.getcwd()), file=sys.stdout)

    for key in os.environ:
        print(key +' = ' + os.environ[key])
    print('''
================================================================================''',
          file=sys.stdout)

def parseArgs():
    '''Parse the arguments - no  options are available at this time'''
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


def merge_branch(source_url, source_branch, target_branch, sourceSHA):

    source_url    = source_url.strip()
    source_branch = source_branch.strip()
    target_branch = target_branch.strip()
    sourceSHA     = sourceSHA.strip()

    remote_list = subprocess.check_output(['git', 'remote', '-v'])
    if 'source_remote' in remote_list:
        print('git remote exists, removing it', file=sys.stdout)
        subprocess.check_call(['git', 'remote', 'rm', 'source_remote'])
    subprocess.check_call(['git', 'remote', 'add', 'source_remote',
                           source_url])

    fetch_succeeded = False
    for i in range(3):
        try:
            subprocess.check_call(['git', 'fetch', 'source_remote',
                                   source_branch])
            fetch_succeeded = True
            break
        except subprocess.CalledProcessError:
            pass
    if not fetch_succeeded:
        raise SystemExit(12)

    subprocess.check_call(['git', 'fetch', 'origin',
                           target_branch])
    subprocess.check_call(['git', 'reset', '--hard',
                           'HEAD'])
    subprocess.check_call(['git', 'checkout', '-B',
                           target_branch, 'origin/' + target_branch])
    subprocess.check_call(['git', 'merge', '--no-edit',
                           'source_remote/' + source_branch]),

    actual_source_SHA = subprocess.check_output(['git', 'rev-parse',
                                                 'source_remote/' + source_branch])

    actual_source_SHA = actual_source_SHA.strip()

    if actual_source_SHA != sourceSHA:
        print('The SHA ({source_sha}) for the last commit on branch {source_branch}'.format(source_sha=actual_source_SHA,
                                                                                            source_branch=source_branch),
              file=sys.stdout)
        print('  in repo {source_repo} is different than the expected SHA,'.format(source_repo=source_url),
              file=sys.stdout)
        print('  which is: {source_sha}.'.format(source_sha=sourceSHA),
              file=sys.stdout)
        raise SystemExit(-1)


def run():
    return_value = True
    try:
        arguments = parseArgs()
    except SystemExit:
        return_value = False
    if return_value:
        os.chdir(os.path.join(arguments.workspaceDir, 'Trilinos'))
        print("Set CWD = {dirName}".format(dirName=os.path.join(
                                                 arguments.workspaceDir,
                                                 'Trilinos')),
              file=sys.stdout)
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
            print('Recieved subprocess.CalledProcessError - returned {error_num}'.format(error_num=cpe.returncode),
                  file=sys.stdout)
            print('  from command {cmd}'.format(cmd=cpe.cmd),
                  file=sys.stdout)
            print('  output {out}'.format(out=cpe.output),
                  file=sys.stdout)
            try:
                print('  stdout {out}'.format(out=cpe.stdout),
                      file=sys.stdout)
            except AttributeError:
                pass
            try:
                print('  stderr {eout}'.format(eout=cpe.stderr),
                      file=sys.stdout)
            except AttributeError:
                pass

    return return_value


if __name__ == '__main__':  # pragma nocover
    returnValue = run()
    if returnValue:
        exit(0)
    else:
        exit(1)
