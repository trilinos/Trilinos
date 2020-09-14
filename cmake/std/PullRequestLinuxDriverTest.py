#!/usr/bin/env python
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
#
# Change shebang line to '/usr/bin/python -3' for python 3.x porting warnings
"""
This script drives a PR testing build.  It assume that Trilinos is already
cloned under $WORKSPACE/Trilinos and that the 'origin' remote points to
$TRILINOS_TARGET_REPO (but that is not checked here).

As long as the ${WORKSPACE}/Trilinos git repo has the correct 'origin', this
script will automatically set it up to do the merge correctly, no matter
what its state before this script is called (i.e. from a past PR
attempt). Unless the Trilinos/.git directory becomes corrupted, there should
*NEVER* be any need to delete and reclone this Trilinos git repo.

This script can be run in a mode where the driver scripts are run from one
Trilinos git repo and operate on another Trilinos git repo that gets
manipulated in order to merge the "source" topic branch into the "target"
branch.  This makes it easy to test changes to the PR scripts.  But if this
script is run from ${WORKSPACE}/Trilinos, then these repos are one and the same
and we get the correct behavior for PR testing.

Expectations
------------

### Required Environment Variables
- MODULESHOME : Path to the location where modulefiles are.
- CC : C Compiler
- FC : Fortran Compiler
- PULLREQUEST_CDASH_TRACK : Which CDash track should this result be published to?

### Other Expectations?

"""
from __future__ import print_function

# turn off generation of the .pyc/.pyo files.
import sys
sys.dont_write_bytecode = True

import argparse
import os
import sys

import trilinosprhelpers



def parse_args():
    """
    This function parses the command line arguments for the script.

    Returns:
        Namespace object containing the arguments from the command line.
    """
    parser    = argparse.ArgumentParser(description='Parse the repo and build information', add_help=False)
    required  = parser.add_argument_group('Required Arguments')
    optional  = parser.add_argument_group('Optional Arguments')
    cwd       = os.getcwd()

    default_workspace = cwd
    if "WORKSPACE" in os.environ.keys():
        default_workspace = os.environ["WORKSPACE"]

    default_package_enables = os.path.join("..", "packageEnables.cmake")
    default_subprojects_file = os.path.join("..", "package_subproject_list.cmake")
    #default_package_enables = os.path.join(default_workspace, "packageEnables.cmake")


    required.add_argument('--sourceRepo',
                          dest="sourceRepo",
                          action='store',
                          help='Repo with the new changes',
                          required=True)

    required.add_argument('--sourceBranch',
                          dest="sourceBranch",
                          action='store',
                          help='Branch with the new changes',
                          required=True)

    required.add_argument('--targetRepo',
                          dest="targetRepo",
                          action='store',
                          help='Repo to merge into',
                          required=True)

    required.add_argument('--targetBranch',
                          dest="targetBranch",
                          action='store',
                          help='Branch to merge into',
                          required=True)

    required.add_argument('--job_base_name',
                          dest="job_base_name",
                          action='store',
                          help='The Jenkins job base name',
                          required=True)

    required.add_argument('--github_pr_number',
                          dest="github_pr_number",
                          action='store',
                          help='The github PR number',
                          required=True)

    required.add_argument('--job_number',
                          dest="job_number",
                          action='store',
                          help='The Jenkins build number',
                          required=True)

    optional.add_argument('--config',
                          dest='configfile',
                          action='store',
                          default=os.path.join(cwd, "pr_config/pullrequest.ini"),
                          help="Configuration file. Default=%(default)s",
                          required=False)

    optional.add_argument('--workspaceDir',
                          dest="workspaceDir",
                          action='store',
                          default=default_workspace,
                          help="The local workspace directory that Jenkins set up." +
                               " Default={}".format(default_workspace),
                          required=False)

    optional.add_argument('--packageEnables',
                          dest="package_enables",
                          action="store",
                          default=default_package_enables,
                          help="Custom packageEnables.cmake override." +
                               " Default={}".format(default_package_enables))

    optional.add_argument('--subprojects_file',
                          dest="subprojects_file",
                          action="store",
                          default=default_subprojects_file,
                          help="Custom subprojects file." +
                               " Default={}".format(default_subprojects_file))

    optional.add_argument('--mode',
                          dest='mode',
                          action='store',
                          default='standard',
                          help="PR testing mode. 'standard' for normal PR tests, 'installation'" +
                               " for installation testing. Default = %(default)s")

    optional.add_argument('--req-mem-per-core',
                          dest='req_mem_per_core',
                          action='store',
                          default=3.0,
                          help="Minimum required memory per core (Gigabytes)." +
                               "Default = %(default)s")

    optional.add_argument('--max-cores-allowed',
                          dest='max_cores_allowed',
                          action='store',
                          default=12,
                          help="Max cores allowed, if >= 0 we will use the # of " +
                               "detected cores on the system. Default = %(default)s")

    optional.add_argument('--num-concurrent-tests',
                          dest='num_concurrent_tests',
                          action='store',
                          default=-1,
                          help="Number of concurrent tests to run. Default = %(default)s")

    optional.add_argument("--dry-run",
                          dest="dry_run",
                          action="store_true",
                          default=False,
                          help="Enable dry-run mode. Runs script but don't execute the build(s) Default = %(default)s")

    optional.add_argument('-h', '--help',
                          action='help',
                          default=argparse.SUPPRESS,
                          help='show this help message and exit')

    arguments = parser.parse_args()

    # Type conversions
    arguments.max_cores_allowed    = int(arguments.max_cores_allowed)
    arguments.num_concurrent_tests = int(arguments.num_concurrent_tests)
    arguments.req_mem_per_core     = float(arguments.req_mem_per_core)

    #if arguments.workspaceDir is None:
        #arguments.workspaceDir = default_workspace
        #if arguments.package_enables is None:
            #arguments.package_enables = default_package_enables
    #else:
        #if arguments.package_enables is None:
            #arguments.package_enables = os.path.join(arguments.workspaceDir, "packageEnables.cmake")


    # Print the arguments to the console
    print("\n" + "="*90)
    print("Parameters")
    print("----------")
    print("- CONFIGURATION_FILE     : {configfile}".format(**vars(arguments)))
    print("- MODE                   : {mode}".format(**vars(arguments)))
    print("- REQ_MEM_PER_CORE       : {req_mem_per_core}".format(**vars(arguments)))
    print("- MAX_CORES_ALLOWED      : {max_cores_allowed}".format(**vars(arguments)))
    print("- NUM_CONCURRENT_TESTS   : {num_concurrent_tests}".format(**vars(arguments)))
    print("")
    print("- JOB_BASE_NAME          : {job_base_name}".format(**vars(arguments)))
    print("- WORKSPACE              : {workspaceDir}".format(**vars(arguments)))
    print("")
    print("- TRILINOS_SOURCE_REPO   : {sourceRepo}".format(**vars(arguments)))
    print("- TRILINOS_SOURCE_BRANCH : {sourceBranch}".format(**vars(arguments)))
    print("")
    print("- TRILINOS_TARGET_REPO   : {targetRepo}".format(**vars(arguments)))
    print("- TRILINOS_TARGET_BRANCH : {targetBranch}".format(**vars(arguments)))
    print("")
    print("- PULLREQUESTNUM         : {github_pr_number}".format(**vars(arguments)))
    print("- BUILD_NUMBER           : {job_number}".format(**vars(arguments)))
    print("- PACKAGE_ENABLES        : {package_enables}".format(**vars(arguments)))
    print("-")
    print("- DRY_RUN                : {dry_run}".format(**vars(arguments)))
    print("" + "="*90 + "\n")

    return arguments



def main(args):

    pr_config = None

    # Banner
    print("="*80)
    print("=")
    print("=   T R I L I N O S   P U L L R E Q U E S T   D R I V E R   S T A R T")
    print("=")
    if args.dry_run:
        print("=                D R Y   R U N   M O D E   E N A B L E D")
        print("=")
    print("="*80)
    print("")

    if 'standard' == args.mode:
        pr_config = trilinosprhelpers.TrilinosPRConfigurationStandard(args)
    elif 'installation' == args.mode:
        pr_config = trilinosprhelpers.TrilinosPRConfigurationInstallation(args)
    else:
        raise KeyError("ERROR: Unknown test mode, {}, was provided.".format(args.mode))

    pr_config.prepare_test()
    pr_config.execute_test()
    return 0



if __name__ == "__main__":
    args = parse_args()
    status = main(args)
    print("Done.")
    sys.exit(status)


