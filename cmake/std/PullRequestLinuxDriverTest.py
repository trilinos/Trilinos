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
"""
from __future__ import print_function

# turn off generation of the .pyc/.pyo files.
import sys
sys.dont_write_bytecode = True

import argparse
import os
import re
import sys

import subprocess
from multiprocessing import cpu_count

sys.path.insert(1, os.path.join(os.environ['MODULESHOME'], 'init'))

try:
    from env_modules_python import module
except ImportError:
    def module(*args):
        if type(args[0]) == type([]):
            args = args[0]
        else:
            args = list(args)
        (output, error) = subprocess.Popen(['/usr/bin/modulecmd', 'python'] +
                                           args,
                                           stdout=subprocess.PIPE).communicate()
        exec(output)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Parse the repo and build information')
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
    parser.add_argument('job_base_name',
                        help='The jenkins job base name')
    parser.add_argument('github_pr_number',
                        help='The github PR number')
    parser.add_argument('job_number',
                        help='The jenkins build number')
    parser.add_argument('workspaceDir',
                        help='The local workspace directory jenkins set up')
    arguments = parser.parse_args()
    return arguments


def print_input_variables(arguments):
    print(
        "\n==========================================================================================",
        file=sys.stdout)
    print("Jenkins Input Variables:", file=sys.stdout)
    print("- JOB_BASE_NAME: {job_base_name}".format(
        job_base_name=arguments.job_base_name),
          file=sys.stdout)
    print("- WORKSPACE    : {workspace}".format(
        workspace=arguments.workspaceDir),
          file=sys.stdout)
    print(
        "\n==========================================================================================",
        file=sys.stdout)
    print("Parameters:", file=sys.stdout)
    print("- TRILINOS_SOURCE_REPO  : {source_repo}".format(
        source_repo=arguments.sourceRepo),
          file=sys.stdout)
    print("- TRILINOS_SOURCE_BRANCH: {source_branch}\n".format(
        source_branch=arguments.sourceBranch),
          file=sys.stdout)
    print("- TRILINOS_TARGET_REPO  : {target_repo}".format(
        target_repo=arguments.targetRepo),
          file=sys.stdout)
    print("- TRILINOS_TARGET_BRANCH: {target_branch}\n".format(
        target_branch=arguments.targetBranch),
          file=sys.stdout)
    print("- PULLREQUESTNUM        : {num}".format(
        num=arguments.github_pr_number),
          file=sys.stdout)
    print("- BUILD_NUMBER          : {num}".format(
        num=arguments.job_number),
          file=sys.stdout)
    print(
        "\n==========================================================================================",
        file=sys.stdout)


def confirmGitVersion():
    git_version_string = subprocess.check_output(['git', '--version'])
    git_version_number_string = git_version_string[git_version_string.rfind(' '):]
    major_git_version = int(git_version_number_string[:git_version_number_string.find('.')])
    minor_git_version = int(git_version_number_string[git_version_number_string.find('.')+1:
                                                      git_version_number_string.rfind('.')])

    if major_git_version  <  2 or (major_git_version == 2 and
                                   minor_git_version < 10):
        raise SystemExit("Git version  should be 2.10 or better - Exiting!")
    else:
        print(git_version_string)


def setBuildEnviron(arguments):
    moduleMap = {'Trilinos_pullrequest_gcc_4.8.4':
                     ['sems-env',
                     'sems-git/2.10.1',
                     'sems-gcc/4.8.4',
                     'sems-openmpi/1.10.1',
                     'sems-python/2.7.9',
                     'sems-boost/1.63.0/base',
                     'sems-zlib/1.2.8/base',
                     'sems-hdf5/1.8.12/parallel',
                     'sems-netcdf/4.4.1/exo_parallel',
                     'sems-parmetis/4.0.3/parallel',
                     'sems-scotch/6.0.3/nopthread_64bit_parallel',
                     'sems-superlu/4.3/base',
                     'sems-cmake/3.10.3',
                     'atdm-env',
                     'atdm-ninja_fortran/1.7.2'],
                 'Trilinos_pullrequest_gcc_4.9.3_SERIAL':
                     ['sems-env',
                      'sems-git/2.10.1',
                      'sems-gcc/4.9.3',
                      'sems-python/2.7.9',
                      'sems-boost/1.63.0/base',
                      'sems-zlib/1.2.8/base',
                      'sems-hdf5/1.8.12/base',
                      'sems-netcdf/4.4.1/exo',
                      'sems-metis/5.1.0/base',
                      'sems-superlu/4.3/base',
                      'sems-cmake/3.10.3',
                      'atdm-env',
                      'atdm-ninja_fortran/1.7.2'],
                 'Trilinos_pullrequest_python_2':
                     ['sems-git/2.10.1',
                      'sems-gcc/7.2.0',
                      ('', 'sems-python/2.7.9'),
                      'sems-cmake/3.10.3',
                      'atdm-env',
                      'atdm-ninja_fortran/1.7.2'],
                'Trilinos_pullrequest_python_3':
                     ['sems-git/2.10.1',
                      'sems-gcc/7.2.0',
                      ('', 'sems-python/2.7.9'),
                      'sems-cmake/3.10.3',
                      'atdm-env',
                      'atdm-ninja_fortran/1.7.2'],
                'Trilinos_pullrequest_gcc_7.2.0':
                     ['sems-env',
                     'sems-git/2.10.1',
                     'sems-gcc/7.2.0',
                     'sems-openmpi/1.10.1',
                     'sems-python/2.7.9',
                     'sems-boost/1.63.0/base',
                     'sems-zlib/1.2.8/base',
                     'sems-hdf5/1.8.12/parallel',
                     'sems-netcdf/4.4.1/exo_parallel',
                     'sems-parmetis/4.0.3/parallel',
                     'sems-scotch/6.0.3/nopthread_64bit_parallel',
                     'sems-superlu/4.3/base',
                     'sems-cmake/3.10.3',
                     'atdm-env',
                     'atdm-ninja_fortran/1.7.2'],
                'Trilinos_pullrequest_intel_17.0.1':
                     ['sems-env',
                     'sems-git/2.10.1',
                     'sems-gcc/4.9.3',
                     'sems-intel/17.0.1',
                     'sems-mpich/3.2',
                     'sems-python/2.7.9',
                     'sems-boost/1.63.0/base',
                     'sems-zlib/1.2.8/base',
                     'sems-hdf5/1.8.12/parallel',
                     'sems-netcdf/4.4.1/exo_parallel',
                     'sems-parmetis/4.0.3/parallel',
                     'sems-scotch/6.0.3/nopthread_64bit_parallel',
                     'sems-superlu/4.3/base',
                     'sems-cmake/3.10.3',
                     'atdm-env',
                     'atdm-ninja_fortran/1.7.2'],
                'Trilinos_pullrequest_cuda_9.2':
                     ['git/2.10.1',
                     'devpack/20180521/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88',
                      ('openblas/0.2.20/gcc/7.2.0', 'netlib/3.8.0/gcc/7.2.0')]}

    environMap = {'Trilinos_pullrequest_gcc_4.8.4':
                      {'OMP_NUM_THREADS': '2'},
                  'Trilinos_pullrequest_gcc_4.9.3_SERIAL':
                      {'OMP_NUM_THREADS': '2'},
                 'Trilinos_pullrequest_python_2':
                      {'PYTHONPATH':
                           os.path.join(os.path.sep,
                                        'projects',
                                        'sierra',
                                        'linux_rh7',
                                        'install',
                                        'Python',
                                        'extras',
                                        'lib',
                                        'python2.7',
                                        'site-packages'),
                       'MANPATH':
                           os.path.join(os.path.sep,
                                        'projects',
                                        'sierra',
                                        'linux_rh7',
                                        'install',
                                        'Python',
                                        '2.7.15',
                                        'share',
                                        'man'),
                       'PATH': os.path.join(os.path.sep,
                                            'projects',
                                            'sierra',
                                            'linux_rh7',
                                            'install',
                                            'Python',
                                            '2.7.15',
                                            'bin') + os.pathsep +
                               os.path.join(os.path.sep,
                                            'projects',
                                            'sierra',
                                            'linux_rh7',
                                            'install',
                                            'Python',
                                            'extras'
                                            'bin')}
        ,
                 'Trilinos_pullrequest_python_3':
                      {'PYTHONPATH':
                           os.path.join(os.path.sep,
                                        'projects',
                                        'sierra',
                                        'linux_rh7',
                                        'install',
                                        'Python',
                                        'extras',
                                        'lib',
                                        'python3.6',
                                        'site-packages'),
                       'MANPATH':
                           os.path.join(os.path.sep,
                                        'projects',
                                        'sierra',
                                        'linux_rh7',
                                        'install',
                                        'Python',
                                        '3.6.3',
                                        'share',
                                        'man'),
                       'PATH': os.path.join(os.path.sep,
                                            'projects',
                                            'sierra',
                                            'linux_rh7',
                                            'install',
                                            'Python',
                                            '3.6.3',
                                            'bin') + os.pathsep +
                               os.path.join(os.path.sep,
                                            'projects',
                                            'sierra',
                                            'linux_rh7',
                                            'install',
                                            'Python',
                                            'extras'
                                            'bin')},
                 'Trilinos_pullrequest_gcc_7.2.0':
                      {'SEMS_FORCE_LOCAL_COMPILER_VERSION': '4.9.3',
                       'OMP_NUM_THREADS': '2'},
                 'Trilinos_pullrequest_intel_17.0.1':
                      {'SEMS_FORCE_LOCAL_COMPILER_VERSION': '4.9.3',
                       'OMP_NUM_THREADS': '2'},
                 'Trilinos_pullrequest_cuda_9.2':
                      {'OMPI_CXX':
                       os.path.join(os.environ['WORKSPACE'],
                                    'Trilinos',
                                    'packages',
                                    'kokkos',
                                    'bin',
                                    'nvcc_wrapper'),
                       'OMPI_CC': os.environ.get('CC', ''),
                       'OMPI_FC': os.environ.get('FC', ''),
                       'CUDA_LAUNCH_BLOCKING': '1',
                       'CUDA_MANAGED_FORCE_DEVICE_ALLOC': '1',
                       'PATH': os.path.join(os.path.sep,
                                            'ascldap',
                                            'users',
                                            'rabartl',
                                            'install',
                                            'white-ride',
                                            'cmake-3.11.2',
                                            'bin') + os.pathsep +
                               os.path.join(os.path.sep,
                                            'ascldap',
                                            'users',
                                            'rabartl',
                                            'install',
                                            'white-ride',
                                            'ninja-1.8.2',
                                            'bin')} }

    try:
        moduleList = moduleMap[arguments.job_base_name]
        l_environMap = environMap[arguments.job_base_name]
    except KeyError:
        sys.exit("""ERROR: Unable to find matching environment for job: Trilinos_pullrequest_UNKOWN
       Error code was: 42""")

    if 'sems-env' == moduleList[0]:
        module('use', '/projects/sems/modulefiles/projects')
    for mod in moduleList:
        if isinstance(mod, str):
            module('load', mod)
        else:
            unl_mod, load_mod = mod
            if '' == unl_mod:
                module('unload', load_mod)
            else:
                module('swap', unl_mod, load_mod)

    if 'OMPI_CC' in l_environMap:
        l_environMap['OMPI_CC'] = os.environ.get('CC', '')
    if 'OMPI_FC' in l_environMap:
        l_environMap['OMPI_FC'] = os.environ.get('FC', '')

    for key, value in l_environMap.items():
        if key in os.environ:
            # we are assuming these are paths to be prepended
            os.environ[key] = value + os.pathsep + os.environ[key]
        else:
            os.environ[key] = value
    confirmGitVersion()

    print ("Environment:\n", file=sys.stdout)
    print("  pwd = {cwd}".format(cwd=os.getcwd()), file=sys.stdout)
    print("", file=sys.stdout)
    for key in os.environ:
        print(key + ' = ' + os.environ[key],
              file=sys.stdout)
    print(
        "\n==========================================================================================",
        file=sys.stdout)
    print(module('list'))


def getCDashTrack():
    returnValue = 'Pull Request'
    if 'PULLREQUEST_CDASH_TRACK' in os.environ:
        returnValue = os.environ['PULLREQUEST_CDASH_TRACK']
        print('PULLREQUEST_CDASH_TRACK is set. Setting CDASH_TRACK={}'.format(
            returnValue) )
        returnValue = os.environ['PULLREQUEST_CDASH_TRACK']
    else:
        print('PULLREQUEST_CDASH_TRACK isn\'t set, using default value')

    return returnValue

def compute_n():
    '''given the default and the hardware environment determine the
     number of processors  to use'''
    try:
        environment_weight = int(os.environ['JENKINS_JOB_WEIGHT'])
    except KeyError:
        environment_weight = 29
    requested_mem_per_core =  3.0

    n_cpu = cpu_count()
    # this assumes unix/linux systems. A more general
    # solution is available in psutil at the cost of
    # using anaconda - which is not  bad.
    with open('/proc/meminfo') as f_ptr:
        meminfo = dict((i.split()[0].rstrip(':'), int(i.split()[1])) for i in
                        f_ptr.readlines())
    mem_kib = meminfo['MemTotal']
    mem_G = mem_kib/(1024*1024)
    number_possible_jobs = max(1, int(n_cpu/environment_weight))
    parallel_level = int(mem_G/(requested_mem_per_core*number_possible_jobs))

    if parallel_level > environment_weight:
        parallel_level = environment_weight
    return parallel_level


config_map = {'Trilinos_pullrequest_gcc_4.8.4': 'PullRequestLinuxGCC4.8.4TestingSettings.cmake',
              'Trilinos_pullrequest_intel_17.0.1': 'PullRequestLinuxIntelTestingSettings.cmake',
              'Trilinos_pullrequest_gcc_4.9.3_SERIAL': 'PullRequestLinuxGCC4.9.3TestingSettingsSERIAL.cmake',
              'Trilinos_pullrequest_gcc_7.2.0': 'PullRequestLinuxGCC7.2.0TestingSettings.cmake',
              'Trilinos_pullrequest_cuda_9.2': 'PullRequestLinuxCuda9.2TestingSettings.cmake',
              'Trilinos_pullrequest_python_2': 'PullRequestLinuxPython2.cmake',
              'Trilinos_pullrequest_python_3': 'PullRequestLinuxPython3.cmake'}


def createPackageEnables(arguments):
    enable_map = {'Trilinos_pullrequest_python_2': 'TrilinosFrameworkTests',
                  'Trilinos_pullrequest_python_3': 'TrilinosFrameworkTests'}

    try:
        if arguments.job_base_name not in enable_map:
            subprocess.check_call([os.path.join(arguments.workspaceDir,
                                                'Trilinos',
                                                'commonTools',
                                                'framework',
                                                'get-changed-trilinos-packages.sh'),
                                   os.path.join('origin', arguments.targetBranch),
                                   'HEAD',
                                   'packageEnables.cmake'])
        else:
            with open('packageEnables.cmake',  'w') as f_out:
                f_out.write('''
MACRO(PR_ENABLE_BOOL  VAR_NAME  VAR_VAL)
  MESSAGE("-- Setting ${VAR_NAME} = ${VAR_VAL}")
  SET(${VAR_NAME} ${VAR_VAL} CACHE BOOL "Set in $CMAKE_PACKAGE_ENABLES_OUT")
ENDMACRO()

PR_ENABLE_BOOL(Trilinos_ENABLE_''' + enable_map[arguments.job_base_name] + ''' ON)
''')
        print('Enabled packages:')
        cmake_rstring = subprocess.check_output(['cmake',
                                                 '-P',
                                                 'packageEnables.cmake'],
                                                stderr=subprocess.STDOUT)
        if sys.version_info.major is not 3:
            print(cmake_rstring)
        else:
            print(str(cmake_rstring, 'ASCII'))
    except subprocess.CalledProcessError as cpe:
        print('There was an issue generating packageEnables.cmake.  '
              'The error code was: {}'.format(cpe.returncode))

def run():
    return_value = True
    arguments = parse_args()

    os.chdir(arguments.workspaceDir)

    print_input_variables(arguments)

    re_src_branchname = "master_merge_[0-9]{8}_[0-9]{6}"
    if "master" == arguments.targetBranch:
        if re.match(re_src_branchname, arguments.sourceBranch) is None:
            sys.exit("""------------------------------------------------------------------------------------------
NOTICE: Destination branch is trilinos/Trilnos::master
ERROR : Source branch is NOT trilinos/Trilinos::master_merge_YYYYMMDD_HHMMSS
      : This violates Trilinos policy, pull requests into the master branch are restricted.
      : Perhaps you forgot to specify the develop branch as the target in your PR?
------------------------------------------------------------------------------------------
""")
        else:
            print("""NOTICE: Source branch IS trilinos/Trilinos::{}
        : This is allowed, proceeding with testing.""".format(arguments.sourceBranch))
    else:
        print("""NOTICE: Destination branch is NOT trilinos/Trilinos::master"
      : PR testing will proceed.""")

    setBuildEnviron(arguments)

    CDash_Track = getCDashTrack()

    createPackageEnables(arguments)

    parallel_level = compute_n()

    build_name = "PR-{PULLREQUESTNUM}-test-{JOB_BASE_NAME}-{BUILD_NUMBER}". \
        format(PULLREQUESTNUM=arguments.github_pr_number,
               JOB_BASE_NAME=arguments.job_base_name,
               BUILD_NUMBER=arguments.job_number)

    config_script = config_map[arguments.job_base_name]

    os.chdir('TFW_testing_single_configure_prototype')
    print('Set CWD = {}'.format(os.getcwd()))

    subprocess.check_call(['ctest', '-S', 'simple_testing.cmake',
                           '-Dbuild_name={}'.format(build_name),
                           '-Dskip_by_parts_submit=OFF',
                           '-Dskip_update_step=ON',
                           '-Ddashboard_model=Experimental',
                           '-Ddashboard_track={}'.format(CDash_Track),
                           '-DPARALLEL_LEVEL={}'.format(parallel_level),
                           '-Dbuild_dir={}/pull_request_test'.format(arguments.workspaceDir),
                           '-Dconfigure_script=' +
                               os.path.join(arguments.workspaceDir,
                                            'Trilinos',
                                            'cmake',
                                            'std',
                                            config_script),
                           '-Dpackage_enables=../packageEnables.cmake',
                           '-Dsubprojects_file=../TFW_single_configure_support_scripts/'+
                           'package_subproject_list.cmake'])

    return return_value


if __name__ == '__main__':  # pragma nocover
    returnValue = run()
    if returnValue:
        exit(0)
    else:
        exit(1)
