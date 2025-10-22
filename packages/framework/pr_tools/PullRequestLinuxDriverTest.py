#!/usr/bin/env python3 -u
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
#
# Change shebang line to '/usr/bin/python -3' for python 3.x porting warnings
"""
This script drives a PR testing build (configure, build, test, report) given
an existing directory of Trilinos source code.
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
    parser    = argparse.ArgumentParser(description='Parse the repo and build information')
    required  = parser.add_argument_group('Required Arguments')
    optional  = parser.add_argument_group('Optional Arguments')
    cwd       = os.getcwd()

    default_workspace = cwd
    if "WORKSPACE" in os.environ.keys():
        default_workspace = os.environ["WORKSPACE"]

    # Set up path to the packageenables aand subprojects files.
    default_filename_packageenables = os.path.join("..", "packageEnables.cmake")
    default_filename_subprojects = os.path.join("..", "package_subproject_list.cmake")

    required.add_argument('--target-branch-name',
                          dest="target_branch_name",
                          action='store',
                          help='Branch to merge into',
                          required=True)

    required.add_argument('--genconfig-build-name',
                          dest="genconfig_build_name",
                          action='store',
                          help='The job base name for the cmake configuration',
                          required=True)

    required.add_argument('--pullrequest-number',
                          dest="pullrequest_number",
                          action='store',
                          help='The github PR number',
                          required=True)

    required.add_argument('--source-dir',
                          dest="source_dir",
                          action='store',
                          help="Directory containing the source code to compile/test.",
                          required=True)

    required.add_argument('--build-dir',
                          dest="build_dir",
                          action='store',
                          help="Path to the build directory.",
                          required=True)

    optional.add_argument('--pullrequest-build-name',
                          dest="pullrequest_build_name",
                          action='store',
                          help='The Jenkins job base name',
                          required=False)

    optional.add_argument('--jenkins-job-number',
                          dest="jenkins_job_number",
                          action='store',
                          help='The Jenkins build number',
                          required=False)

    optional.add_argument('--dashboard-build-name',
                          dest="dashboard_build_name",
                          action='store',
                          help='The build name posted by ctest to a dashboard',
                          required=False)

    optional.add_argument('--use-explicit-cachefile',
                          dest='use_explicit_cachefile',
                          action='store_true',
                          default=False,
                          help="Use -DTrilinos_CONFIGURE_OPTIONS_FILE instead of -C.",
                          required=False)

    optional.add_argument('--ctest-driver',
                          dest="ctest_driver",
                          action='store',
                          help="Location of the CTest driver script to load via `-S`.",
                          required=False)

    optional.add_argument('--ctest-drop-site',
                          dest="ctest_drop_site",
                          action='store',
                          default="testing.sandia.gov",
                          help="URL of the cdash server to post to.",
                          required=False)

    optional.add_argument('--pullrequest-cdash-track',
                          dest='pullrequest_cdash_track',
                          action='store',
                          default="Pull Request",
                          help="The CDash Track to add results to. Default=%(default)s",
                          required=False)

    optional.add_argument('--pullrequest-env-config-file',
                          dest='pullrequest_env_config_file',
                          action='store',
                          default=os.path.join(cwd, "pr_config/pullrequest.ini"),
                          help="The Trilinos PR driver configuration file " + \
                               "containing job mappings to environment specifications. Default=%(default)s",
                          required=False)

    optional.add_argument('--pullrequest-gen-config-file',
                          dest='pullrequest_gen_config_file',
                          action='store',
                          default=os.path.join(cwd, "../GenConfig/src/gen-config.ini"),
                          help="The Trilinos PR driver configuration file " + \
                               "containing job mappings to cmake specifications. Default=%(default)s",
                          required=False)

    optional.add_argument('--workspace-dir',
                          dest="workspace_dir",
                          action='store',
                          default=default_workspace,
                          help="The local workspace directory that Jenkins set up." +
                               " Default={}".format(default_workspace),
                          required=False)

    desc_package_enables = "The packageEnables.cmake is usually generated by TriBiTS infrastructure " + \
                           "based on which packages contain the changes between the source and target " + \
                           "branches."

    optional.add_argument('--filename-packageenables',
                          dest="filename_packageenables",
                          action="store",
                          default=default_filename_packageenables,
                          help="{} Default={}".format(desc_package_enables, default_filename_packageenables))

    optional.add_argument('--skip-create-packageenables',
                          dest="skip_create_packageenables",
                          action="store_true",
                          default=False,
                          help="Skip the creation of the packageEnables.cmake fragment file generated by " + \
                          "the TriBITS infrastructure indicating which packages are to be enabled based on file " + \
                          "changes between a source and target branch. Default=False")

    desc_subprojects_file = "The subprojects_file is used by the testing infrastructure. This parameter " + \
                            "allows the default, generated file, to be overridden. Generally this should " + \
                            "not be changed from the defaults."

    optional.add_argument('--filename-subprojects',
                          dest="filename_subprojects",
                          action="store",
                          default=default_filename_subprojects,
                          help="{}. Default={}".format(desc_subprojects_file, default_filename_subprojects))

    optional.add_argument('--test-mode',
                          dest='test_mode',
                          action='store',
                          default='standard',
                          help="PR testing mode. Use 'standard' for normal PR tests, 'installation'" +
                               " for installation testing. Default = %(default)s")

    optional.add_argument('--req-mem-per-core',
                          dest='req_mem_per_core',
                          action='store',
                          default=3.0,
                          help="Minimum required memory per core (GB) to build Trilinos." +
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
                          help="Set the number of concurrent tests allowed in CTest. " + \
                               "This is equivalent to `ctest -j <num-concurrent-tests>`. "
                               "If > 0 then this value is used, otherwise the value is calculated " + \
                               "based on number_of_available_cores / max_test_parallelism" + \
                               " Default = %(default)s")

    optional.add_argument('--slots-per-gpu',
                          dest='slots_per_gpu',
                          action='store',
                          default=2,
                          help="If the machine has GPUs, this value will be the number " + \
                               "of resource slots allowed per GPU (e.g. 4 would allow 4 MPI " + \
                               "ranks to talk to each GPU)" + \
                               " Default = %(default)s")

    optional.add_argument("--enable-ccache",
                          dest="ccache_enable",
                          action="store_true",
                          default=False,
                          help="Enable ccache object caching to improve build times. Default = %(default)s")

    optional.add_argument("--dry-run",
                          dest="dry_run",
                          action="store_true",
                          default=False,
                          help="Enable dry-run mode. Script will run but not execute the build steps. Default = %(default)s")

    optional.add_argument("--skip-run-tests",
                          dest="skip_run_tests",
                          action="store_true",
                          default=False,
                          help="Skip running tests in SimpleTesting test stage. Default = %(default)s")

    optional.add_argument("--extra-configure-args",
                          dest="extra_configure_args",
                          action="store",
                          default="",
                          help="Extra arguments that will be passed to CMake for configuring Trilinos.")

    arguments = parser.parse_args()

    if not arguments.ctest_driver:
        arguments.ctest_driver = os.path.join(arguments.source_dir, 'cmake', 'SimpleTesting', 'cmake', 'ctest-driver.cmake')

    # Type conversions
    arguments.max_cores_allowed    = int(arguments.max_cores_allowed)
    arguments.num_concurrent_tests = int(arguments.num_concurrent_tests)
    arguments.req_mem_per_core     = float(arguments.req_mem_per_core)

    # Print the arguments to the console
    print("\n")
    print("+" + "="*78 + "+")
    print("| PullRequestLinuxDriverTest Parameters")
    print("+" + "="*78 + "+")
    print("| - [R] target_branch_name          : {target_branch_name}".format(**vars(arguments)))
    print("| - [R] pullrequest-build-name      : {pullrequest_build_name}".format(**vars(arguments)))
    print("| - [R] genconfig-build-name        : {genconfig_build_name}".format(**vars(arguments)))
    print("| - [R] pullrequest-number          : {pullrequest_number}".format(**vars(arguments)))
    print("| - [R] jenkins-job-number          : {jenkins_job_number}".format(**vars(arguments)))
    print("| - [R] source-dir                  : {source_dir}".format(**vars(arguments)))
    print("| - [R] build-dir                   : {build_dir}".format(**vars(arguments)))
    print("| - [R] ctest-driver                : {ctest_driver}".format(**vars(arguments)))
    print("| - [R] ctest-drop-site             : {ctest_drop_site}".format(**vars(arguments)))
    print("|")
    print("| - [O] dry-run                     : {dry_run}".format(**vars(arguments)))
    print("| - [O] enable-ccache               : {ccache_enable}".format(**vars(arguments)))
    print("| - [O] filename-packageenables     : {filename_packageenables}".format(**vars(arguments)))
    print("| - [O] max-cores-allowed           : {max_cores_allowed}".format(**vars(arguments)))
    print("| - [O] num-concurrent-tests        : {num_concurrent_tests}".format(**vars(arguments)))
    print("| - [O] pullrequest-cdash-track     : {pullrequest_cdash_track}".format(**vars(arguments)))
    print("| - [O] pullrequest-env-config-file : {pullrequest_env_config_file}".format(**vars(arguments)))
    print("| - [O] pullrequest-gen-config-file : {pullrequest_gen_config_file}".format(**vars(arguments)))
    print("| - [O] req-mem-per-core            : {req_mem_per_core}".format(**vars(arguments)))
    print("| - [O] test-mode                   : {test_mode}".format(**vars(arguments)))
    print("| - [O] workspace-dir               : {workspace_dir}".format(**vars(arguments)))
    print("| - [O] skip-run-tests              : {skip_run_tests}".format(**vars(arguments)))
    print("| - [O] extra_configure_args        : {extra_configure_args}".format(**vars(arguments)))
    print("| - [O] dashboard_build_name        : {dashboard_build_name}".format(**vars(arguments)))
    print("| - [O] use_explicit_cachefile       : {use_explicit_cachefile}".format(**vars(arguments)))
    #print("| - [O] : {}".format(**vars(arguments)))
    print("+" + "="*78 + "+")

    return arguments



def main(args):

    pr_config = None

    # Banner
    print("+" + "="*78 + "+")
    print("|")
    print("|   T R I L I N O S   P U L L R E Q U E S T   D R I V E R   S T A R T")
    print("|")
    if args.dry_run:
        print("|                D R Y   R U N   M O D E   E N A B L E D")
        print("|")
    print("+" + "="*78 + "+")
    print("")

    if 'standard' == args.test_mode:
        pr_config = trilinosprhelpers.TrilinosPRConfigurationStandard(args)
    elif 'installation' == args.test_mode:
        pr_config = trilinosprhelpers.TrilinosPRConfigurationInstallation(args)
    else:
        raise KeyError("ERROR: Unknown test mode, {}, was provided.".format(args.test_mode))

    pr_config.prepare_test()

    status = pr_config.execute_test()

    print("PullRequestLinuxDriverTest.py main()> Done.")

    return status



if __name__ == "__main__":
    args = parse_args()
    status = main(args)
    sys.exit(status)


