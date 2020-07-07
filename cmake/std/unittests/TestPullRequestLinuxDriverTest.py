#!/usr/bin/env python
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
'''
Tests for the Test chunk of the Driver script
'''
from __future__ import print_function
import sys
sys.dont_write_bytecode = True

import os
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from textwrap import dedent

import unittest

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

try:
    import mock
except ImportError:  # pragma nocover
    import unittest.mock as mock


from argparse import Namespace
from subprocess import CalledProcessError

if 'MODULESHOME' not in os.environ: # for things like our macs
    os.environ['MODULESHOME'] = os.getcwd()

import PullRequestLinuxDriverTest



class Test_run(unittest.TestCase):
    """Does the run script exist?"""

    def setUp(self):
        self.source_branch = 'incoming_branch'
        self.source_url = '/dev/null/source/Trilinos.git'
        self.target_branch = 'base_branch'
        self.target_url = '/dev/null/target/Trilinos.git'
        self.job_base_name = 'JenkinsBaseName'
        self.github_pr_number = '8888'
        self.jenkins_build_number = '7777'
        self.jenkins_workspace='/dev/null/workspace'

        self.IOredirect = mock.patch('sys.stdout', new_callable=StringIO)
        self.m_chdir = mock.patch('os.chdir')
        self.m_getcwd = mock.patch('os.getcwd',
                                   side_effect=['/dev/null/workspace',
                                                '/dev/null/workspace/TFW_testing_single_configure_prototype'])
        self.m_check_out = mock.patch('subprocess.check_output',
                                      return_value='git version 2.10.1')
        self.m_check_call = mock.patch('subprocess.check_call')
        self.m_config_map = mock.patch.dict(PullRequestLinuxDriverTest.config_map,
                                            {'JenkinsBaseName': 'dummyConfig.cmake'},
                                            clear=False)
        self.m_argv = mock.patch.object(sys, 'argv', ['programName',
                                        self.source_branch,
                                        self.source_url,
                                        self.target_branch,
                                        self.target_url,
                                        self.job_base_name,
                                        self.github_pr_number,
                                        self.jenkins_build_number,
                                        self.jenkins_workspace])
        self.m_environ = mock.patch.dict(os.environ, {'JOB_BASE_NAME': self.job_base_name,
                                                     'JOB_NAME': 'TEST_JOB_NAME',
                                                     'WORKSPACE': self.jenkins_workspace,
                                                     'NODE_NAME': 'TEST_NODE_NAME'},
                                         clear=True)


    def test_verifyGit_fails_with_old_version(self):
        """Check to see that git is in path"""
        with self.m_check_out as m_check_out:
            m_check_out.return_value='git version 1.10.1'

            bad_git_string = 'Git version  should be 2.10 or better - Exiting!'
            if(sys.version_info.major != 3):
                with self.assertRaisesRegexp(SystemExit, bad_git_string):
                    PullRequestLinuxDriverTest.confirmGitVersion()
            else:
                with self.assertRaisesRegex(SystemExit, bad_git_string):
                    PullRequestLinuxDriverTest.confirmGitVersion()

            m_check_out.assert_called_once_with(['git', '--version'])


    def test_verifyGit_passes_with_2_10(self):
        """Check to see that git is in path"""
        with self.m_check_out as m_check_out:
            m_check_out.return_value='git version 2.10.1'

            with self.IOredirect:
                PullRequestLinuxDriverTest.confirmGitVersion()
            m_check_out.assert_called_once_with(['git', '--version'])


    def test_verifyGit_passes_with_2_12(self):
        """Check to see that git is in path"""
        with self.m_check_out as m_check_out:
            m_check_out.return_value='git version 2.12.4'

            with self.IOredirect:
                PullRequestLinuxDriverTest.confirmGitVersion()
        m_check_out.assert_called_once_with(['git', '--version'])


    def test_verifyGit_passes_with_3_x(self):
        """Check to see that git is in path"""
        with self.m_check_out as m_check_out:
            m_check_out.return_value='git version 3.6.1'

            with self.IOredirect:
                PullRequestLinuxDriverTest.confirmGitVersion()
        m_check_out.assert_called_once_with(['git', '--version'])


    def test_verifyTargetBranch_fails_with_master_target_non_mm_source(self):
        """Check to see that git is in path"""
        l_argv = mock.patch.object(sys, 'argv', ['programName',
                                   self.source_url,
                                   self.source_branch,
                                   self.target_url,
                                   'master',
                                   self.job_base_name,
                                   self.github_pr_number,
                                   self.jenkins_build_number,
                                   self.jenkins_workspace])
        with self.IOredirect, \
                self.m_chdir, \
                self.m_check_out, \
                l_argv, \
                self.m_environ:

            bad_branch_string = """------------------------------------------------------------------------------------------
NOTICE: Destination branch is trilinos/Trilnos::master
ERROR : Source branch is NOT trilinos/Trilinos::master_merge_YYYYMMDD_HHMMSS
      : This violates Trilinos policy, pull requests into the master branch are restricted.
      : Perhaps you forgot to specify the develop branch as the target in your PR?
*"""
            if(sys.version_info.major != 3):
                with self.assertRaisesRegexp(SystemExit, bad_branch_string):
                    PullRequestLinuxDriverTest.run()
            else:
                with self.assertRaisesRegex(SystemExit, bad_branch_string):
                    PullRequestLinuxDriverTest.run()


    def test_verifyTargetBranch_passes_with_master_target_mm_source(self):
        """Check to see that git is in path"""
        l_argv = mock.patch.object(sys, 'argv', ['programName',
                                   self.source_url,
                                   'master_merge_20200130_120155',
                                   self.target_url,
                                   'master',
                                   self.job_base_name,
                                   self.github_pr_number,
                                   self.jenkins_build_number,
                                   self.jenkins_workspace])
        l_environ = mock.patch.dict(os.environ, {'JOB_BASE_NAME': self.job_base_name,
                                                 'JOB_NAME': 'TEST_JOB_NAME',
                                                 'WORKSPACE': self.jenkins_workspace,
                                                 'NODE_NAME': 'TEST_NODE_NAME',
                                                 'JENKINS_TEST_WEIGHT': '8'},
                                         clear=True)
        with l_environ:
            env_string_io = StringIO()
            for key in os.environ:
                print(key + ' = ' + os.environ[key],
                      file=env_string_io)
        expected_out = """
==========================================================================================
Jenkins Input Variables:
- JOB_BASE_NAME: JenkinsBaseName
- WORKSPACE    : /dev/null/workspace

==========================================================================================
Parameters:
- TRILINOS_SOURCE_REPO  : /dev/null/source/Trilinos.git
- TRILINOS_SOURCE_BRANCH: master_merge_20200130_120155

- TRILINOS_TARGET_REPO  : /dev/null/target/Trilinos.git
- TRILINOS_TARGET_BRANCH: master

- PULLREQUESTNUM        : 8888
- BUILD_NUMBER          : 7777

==========================================================================================
NOTICE: Source branch IS trilinos/Trilinos::master_merge_20200130_120155
        : This is allowed, proceeding with testing.
Set CWD = /dev/null/workspace
""".format(environ=env_string_io.getvalue())

        with self.IOredirect as m_output, \
                self.m_getcwd, \
                self.m_chdir, \
                self.m_check_out, \
                self.m_config_map, \
                self.m_check_call as m_call, \
                l_argv, \
                l_environ, \
                mock.patch('PullRequestLinuxDriverTest.createPackageEnables'), \
                mock.patch('PullRequestLinuxDriverTest.setBuildEnviron'), \
                mock.patch('PullRequestLinuxDriverTest.compute_n', return_value=20), \
                mock.patch('PullRequestLinuxDriverTest.getCDashTrack') as m_cdtr:
            PullRequestLinuxDriverTest.run()

        self.assertEqual(expected_out, m_output.getvalue())
        m_call.assert_called_once_with(['ctest', '-S', 'simple_testing.cmake',
                                        '-Dbuild_name=PR-8888-test-JenkinsBaseName-7777',
                                        '-Dskip_by_parts_submit=OFF',
                                        '-Dskip_update_step=ON',
                                        '-Ddashboard_model=Experimental',
                                        '-Ddashboard_track={}'.format(m_cdtr.return_value),
                                        '-DPARALLEL_LEVEL=20',
                                        '-DTEST_PARALLEL_LEVEL=8',
                                        '-Dbuild_dir=/dev/null/workspace/pull_request_test',
                                        '-Dconfigure_script=/dev/null/workspace/Trilinos/cmake/std/dummyConfig.cmake',
                                        '-Dpackage_enables=../packageEnables.cmake',
                                        '-Dsubprojects_file=../package_subproject_list.cmake'])


    def test_verifyTargetBranch_passes_with_develop_target(self):
        """Check to see that git is in path"""
        l_argv = mock.patch.object(sys, 'argv', ['programName',
                                   self.source_url,
                                   self.source_branch,
                                   self.target_url,
                                   'develop',
                                   self.job_base_name,
                                   self.github_pr_number,
                                   self.jenkins_build_number,
                                   self.jenkins_workspace])

        with self.m_environ:
            env_string_io = StringIO()
            for key in os.environ:
                print(key + ' = ' + os.environ[key], file=env_string_io)

        expected_out = """
==========================================================================================
Jenkins Input Variables:
- JOB_BASE_NAME: JenkinsBaseName
- WORKSPACE    : /dev/null/workspace

==========================================================================================
Parameters:
- TRILINOS_SOURCE_REPO  : /dev/null/source/Trilinos.git
- TRILINOS_SOURCE_BRANCH: incoming_branch

- TRILINOS_TARGET_REPO  : /dev/null/target/Trilinos.git
- TRILINOS_TARGET_BRANCH: develop

- PULLREQUESTNUM        : 8888
- BUILD_NUMBER          : 7777

==========================================================================================
NOTICE: Destination branch is NOT trilinos/Trilinos::master"
      : PR testing will proceed.
Set CWD = /dev/null/workspace
""".format(environ=env_string_io.getvalue())

        with self.IOredirect as m_output, \
                self.m_chdir, \
                self.m_getcwd, \
                self.m_config_map, \
                self.m_check_out, \
                self.m_check_call as m_call, \
                l_argv, \
                self.m_environ, \
                mock.patch('PullRequestLinuxDriverTest.createPackageEnables'), \
                mock.patch('PullRequestLinuxDriverTest.setBuildEnviron'), \
                mock.patch('PullRequestLinuxDriverTest.compute_n', return_value=20), \
                mock.patch('PullRequestLinuxDriverTest.getCDashTrack', return_value='testTrack'):
            PullRequestLinuxDriverTest.run()

        self.assertEqual(expected_out, m_output.getvalue())
        m_call.assert_called_once_with(['ctest', '-S', 'simple_testing.cmake',
                                        '-Dbuild_name=PR-8888-test-JenkinsBaseName-7777',
                                        '-Dskip_by_parts_submit=OFF',
                                        '-Dskip_update_step=ON',
                                        '-Ddashboard_model=Experimental',
                                        '-Ddashboard_track=testTrack',
                                        '-DPARALLEL_LEVEL=20',
                                        '-DTEST_PARALLEL_LEVEL=20',
                                        '-Dbuild_dir=/dev/null/workspace/pull_request_test',
                                        '-Dconfigure_script=/dev/null/workspace/Trilinos/cmake/std/dummyConfig.cmake',
                                        '-Dpackage_enables=../packageEnables.cmake',
                                        '-Dsubprojects_file=../package_subproject_list.cmake'])



class Test_createPackageEnables(unittest.TestCase):

    def setUp(self):
        self.source_branch = 'incoming_branch'
        self.source_url = '/dev/null/source/Trilinos.git'
        self.target_branch = 'base_branch'
        self.target_url = '/dev/null/target/Trilinos.git'
        self.job_base_name = 'JenkinsBaseName'
        self.github_pr_number = '8888'
        self.jenkins_build_number = '7777'
        self.jenkins_workspace='/dev/null/workspace'

        self.arguments = Namespace()
        setattr(self.arguments, 'sourceBranch', self.source_branch)
        setattr(self.arguments, 'sourceRepo', self.source_url)
        setattr(self.arguments, 'targetBranch', self.target_branch)
        setattr(self.arguments, 'targetRepo', self.target_url)
        setattr(self.arguments, 'job_base_name', self.job_base_name)
        setattr(self.arguments, 'github_pr_number', self.github_pr_number)
        setattr(self.arguments, 'workspaceDir', self.jenkins_workspace)


    def success_side_effect(self):
        with open('packageEnables.cmake',  'w') as f_out:
            f_out.write(dedent('''\
            MACRO(PR_ENABLE_BOOL  VAR_NAME  VAR_VAL)
              MESSAGE("-- Setting ${VAR_NAME} = ${VAR_VAL}")
              SET(${VAR_NAME} ${VAR_VAL} CACHE BOOL "Set in $CMAKE_PACKAGE_ENABLES_OUT")
            ENDMACRO()
            '''))
            f_out.write("PR_ENABLE_BOOL(Trilinos_ENABLE_FooPackageBar ON)")
        with open ('package_subproject_list.cmake', 'w') as f_out:
            f_out.write(dedent('''\
            set(CTEST_LABELS_FOR_SUBPROJECTS TrilinosFrameworkTests '''))

    def test_call_success(self):
        expected_output = '''Enabled packages:
-- Setting Trilinos_ENABLE_FooPackageBar = ON

'''
        with mock.patch('subprocess.check_call',
                        side_effect=self.success_side_effect()) as m_out, \
                mock.patch('sys.stdout',
                           new_callable=StringIO) as m_stdout:
            PullRequestLinuxDriverTest.createPackageEnables(self.arguments)
        m_out.assert_called_once_with([os.path.join(self.jenkins_workspace,
                                                    'Trilinos',
                                                    'commonTools',
                                                    'framework',
                                                    'get-changed-trilinos-packages.sh'),
                                       os.path.join('origin',
                                                    self.target_branch),
                                       'HEAD', 'packageEnables.cmake',
                                       'package_subproject_list.cmake'])
        self.assertEqual(expected_output, m_stdout.getvalue())
        os.unlink('packageEnables.cmake')
        os.unlink('package_subproject_list.cmake')


    def test_call_python2(self):
        expected_output = '''Enabled packages:
-- Setting Trilinos_ENABLE_TrilinosFrameworkTests = ON

'''

        l_arguments = self.arguments
        l_arguments.job_base_name = 'Trilinos_pullrequest_python_2'
        with mock.patch('subprocess.check_call',
                        side_effect=self.success_side_effect()) as m_out, \
            mock.patch('sys.stdout', new_callable=StringIO) as m_stdout:
            PullRequestLinuxDriverTest.createPackageEnables(l_arguments)
        m_out.assert_not_called()
        self.assertEqual(expected_output, m_stdout.getvalue())
        os.unlink('packageEnables.cmake')
        os.unlink('package_subproject_list.cmake')


    def test_call_failure(self):
        expected_output = '''There was an issue generating packageEnables.cmake.  The error code was: 39
'''
        with mock.patch('subprocess.check_call',
                        side_effect=CalledProcessError(cmd='cmake',
                                                       returncode=39)) as m_out, \
                 mock.patch('sys.stdout',
                            new_callable=StringIO) as m_stdout:
             PullRequestLinuxDriverTest.createPackageEnables(self.arguments)
        m_out.assert_called_once_with([os.path.join(self.jenkins_workspace,
                                                    'Trilinos',
                                                    'commonTools',
                                                    'framework',
                                                    'get-changed-trilinos-packages.sh'),
                                       os.path.join('origin',
                                                    self.target_branch),
                                       'HEAD', 'packageEnables.cmake',
                                       'package_subproject_list.cmake'])
        self.assertEqual(expected_output, m_stdout.getvalue())



class Test_setEnviron(unittest.TestCase):
    """Does the script exist?"""

    def setUp(self):
        self.source_branch = 'incoming_branch'
        self.source_url = '/dev/null/source/Trilinos.git'
        self.target_branch = 'base_branch'
        self.target_url = '/dev/null/target/Trilinos.git'
        self.job_base_name='Trilinos_pullrequest_UNKOWN'
        self.github_pr_number= '8888'
        self.jenkins_workspace='/dev/null/workspace'

        self.IOredirect = mock.patch('sys.stdout', new_callable=StringIO)
        self.m_chdir = mock.patch('os.chdir')
        self.m_check_out = mock.patch('subprocess.check_output',
                                      return_value='git version 2.10.1')

        self.m_environ = mock.patch.dict(os.environ, {'JOB_BASE_NAME': self.job_base_name,
                                                     'JOB_NAME': 'TEST_JOB_NAME',
                                                     'WORKSPACE': self.jenkins_workspace,
                                                     'NODE_NAME': 'TEST_NODE_NAME',
                                                     'PATH': '/fake/path',},
                                         clear=True)
        self.arguments = Namespace()
        setattr(self.arguments, 'sourceBranch', self.source_branch)
        setattr(self.arguments, 'sourceRepo', self.source_url)
        setattr(self.arguments, 'targetBranch', self.target_branch)
        setattr(self.arguments, 'targetRepo', self.target_url)
        setattr(self.arguments, 'job_base_name', self.job_base_name)
        setattr(self.arguments, 'github_pr_number', self.github_pr_number)
        setattr(self.arguments, 'workspaceDir', self.jenkins_workspace)


    def test_buildEnv_fails_with_unknown(self):
        """Find the function"""
        expected_output = """ERROR: Unable to find matching environment for job: Trilinos_pullrequest_UNKOWN
       Error code was: 42"""
        with self.IOredirect, \
                self.m_chdir, \
                self.m_check_out, \
                self.m_environ:
            if(sys.version_info.major != 3):
                with self.assertRaisesRegexp(SystemExit, expected_output):
                    PullRequestLinuxDriverTest.setBuildEnviron(self.arguments)
            else:
                with self.assertRaisesRegex(SystemExit, expected_output):
                    PullRequestLinuxDriverTest.setBuildEnviron(self.arguments)


    def buildEnv_passes(self, PR_name, expected_list,
                        test_ENV={'PATH': '/add/CC/path'}):
        setattr(self.arguments,
                'job_base_name',
                PR_name)

        def add_CC(*args, **kwargs):
            os.environ['CC'] = '/fake/gcc/path/bin/gcc'
            os.environ['FC'] = '/fake/gcc/path/bin/gfortran'
            os.environ['PATH'] = '/add/CC/path' + os.environ['PATH']

        with self.IOredirect, \
             self.m_chdir, \
             self.m_check_out, \
             self.m_environ, \
             mock.patch('PullRequestLinuxDriverTest.module',
                        side_effect=add_CC) as m_mod:
            PullRequestLinuxDriverTest.setBuildEnviron(self.arguments)
            for key, value in test_ENV.items():
                if isinstance(value, str):
                    self.assertEqual(os.environ[key], value)
                else:
                    for l_str in value:
                        self.assertNotEqual(os.environ[key].find(l_str), -1)

        m_mod.assert_has_calls(expected_list, any_order=True)


    def test_buildEnv_passes_with_python2(self):
        """Find the function"""
        PR_name = 'Trilinos_pullrequest_python_2'
        expected_list = [mock.call('load', 'sems-git/2.10.1'),
                         mock.call('load', 'sems-gcc/7.2.0'),
                         mock.call('unload', 'sems-python/2.7.9'),
                         mock.call('load', 'sems-cmake/3.10.3'),
                         mock.call('load', 'atdm-env'),
                         mock.call('load', 'atdm-ninja_fortran/1.7.2'),
                         ]
        expected_env = {'PYTHONPATH':
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
                        'PATH': [os.path.join(os.path.sep,
                                              'projects',
                                              'sierra',
                                              'linux_rh7',
                                              'install',
                                              'Python',
                                              '2.7.15',
                                              'bin'),
                                 os.path.join(os.path.sep,
                                              'projects',
                                              'sierra',
                                              'linux_rh7',
                                              'install',
                                              'Python',
                                              'extras'
                                              'bin'),
                                 '/fake/path']}

        self.buildEnv_passes(PR_name, expected_list, expected_env)


    def test_buildEnv_passes_with_python3(self):
        """Find the function"""
        PR_name = 'Trilinos_pullrequest_python_3'
        expected_list = [mock.call('load', 'sems-git/2.10.1'),
                         mock.call('load', 'sems-gcc/7.2.0'),
                         mock.call('unload', 'sems-python/2.7.9'),
                         mock.call('load', 'sems-cmake/3.10.3'),
                         mock.call('load', 'atdm-env'),
                         mock.call('load', 'atdm-ninja_fortran/1.7.2'),
                         ]
        expected_env = {'PYTHONPATH':
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
                        'PATH': [os.path.join(os.path.sep,
                                              'projects',
                                              'sierra',
                                              'linux_rh7',
                                              'install',
                                              'Python',
                                              '3.6.3',
                                              'bin'),
                                 os.path.join(os.path.sep,
                                              'projects',
                                              'sierra',
                                              'linux_rh7',
                                              'install',
                                              'Python',
                                              'extras'
                                              'bin'),
                                 '/fake/path']}

        self.buildEnv_passes(PR_name, expected_list, expected_env)


    def test_buildEnv_passes_with_gcc_484(self):
        """Find the function"""
        PR_name = 'Trilinos_pullrequest_gcc_4.8.4'
        expected_list = [mock.call('use', '/projects/sems/modulefiles/projects'),
                         mock.call('load', 'sems-env'),
                         mock.call('load', 'sems-git/2.10.1'),
                         mock.call('load', 'sems-gcc/4.8.4'),
                         mock.call('load', 'sems-openmpi/1.10.1'),
                         mock.call('load', 'sems-python/2.7.9'),
                         mock.call('load', 'sems-boost/1.63.0/base'),
                         mock.call('load', 'sems-zlib/1.2.8/base'),
                         mock.call('load', 'sems-hdf5/1.10.6/parallel'),
                         mock.call('load', 'sems-netcdf/4.7.3/parallel'),
                         mock.call('load', 'sems-parmetis/4.0.3/parallel'),
                         mock.call('load', 'sems-scotch/6.0.3/nopthread_64bit_parallel'),
                         mock.call('load', 'sems-superlu/4.3/base'),
                         mock.call('load', 'sems-cmake/3.10.3'),
                         mock.call('load', 'atdm-env'),
                         mock.call('load', 'atdm-ninja_fortran/1.7.2'),
                         ]
        self.buildEnv_passes(PR_name, expected_list,
                             test_ENV={'OMP_NUM_THREADS': '2'})


    def test_buildEnv_passes_with_gcc_493_Serial(self):
        """Find the function"""
        PR_name = 'Trilinos_pullrequest_gcc_4.9.3_SERIAL'
        expected_list = [mock.call('use', '/projects/sems/modulefiles/projects'),
                         mock.call('load', 'sems-env'),
                         mock.call('load', 'sems-git/2.10.1'),
                         mock.call('load', 'sems-gcc/4.9.3'),
                         mock.call('load', 'sems-python/2.7.9'),
                         mock.call('load', 'sems-boost/1.63.0/base'),
                         mock.call('load', 'sems-zlib/1.2.8/base'),
                         mock.call('load', 'sems-hdf5/1.10.6/base'),
                         mock.call('load', 'sems-netcdf/4.7.3/base'),
                         mock.call('load', 'sems-metis/5.1.0/base'),
                         mock.call('load', 'sems-superlu/4.3/base'),
                         mock.call('load', 'sems-cmake/3.10.3'),
                         mock.call('load', 'atdm-env'),
                         mock.call('load', 'atdm-ninja_fortran/1.7.2'),
                         ]
        self.buildEnv_passes(PR_name, expected_list,
                             test_ENV={'OMP_NUM_THREADS': '2'})


    def test_buildEnv_passes_with_gcc_720(self):
        """Find the function"""
        PR_name = 'Trilinos_pullrequest_gcc_7.2.0'
        expected_list = [mock.call('use', '/projects/sems/modulefiles/projects'),
                         mock.call('load', 'sems-env'),
                         mock.call('load', 'sems-git/2.10.1'),
                         mock.call('load', 'sems-gcc/7.2.0'),
                         mock.call('load', 'sems-openmpi/1.10.1'),
                         mock.call('load', 'sems-python/2.7.9'),
                         mock.call('load', 'sems-boost/1.63.0/base'),
                         mock.call('load', 'sems-zlib/1.2.8/base'),
                         mock.call('load', 'sems-hdf5/1.10.6/parallel'),
                         mock.call('load', 'sems-netcdf/4.7.3/parallel'),
                         mock.call('load', 'sems-parmetis/4.0.3/parallel'),
                         mock.call('load', 'sems-scotch/6.0.3/nopthread_64bit_parallel'),
                         mock.call('load', 'sems-superlu/4.3/base'),
                         mock.call('load', 'sems-cmake/3.10.3'),
                         mock.call('load', 'atdm-env'),
                         mock.call('load', 'atdm-ninja_fortran/1.7.2'),
                         ]
        self.buildEnv_passes(PR_name, expected_list,
                             test_ENV={'OMP_NUM_THREADS': '2'})


    def test_buildEnv_passes_with_gcc_830(self):
        """Find the function"""
        PR_name = 'Trilinos_pullrequest_gcc_8.3.0'
        expected_list = [mock.call('use', '/projects/sems/modulefiles/projects'),
                         mock.call('load', 'sems-env'),
                         mock.call('load', 'sems-git/2.10.1'),
                         mock.call('load', 'sems-gcc/8.3.0'),
                         mock.call('load', 'sems-openmpi/1.10.1'),
                         mock.call('load', 'sems-python/2.7.9'),
                         mock.call('load', 'sems-boost/1.66.0/base'),
                         mock.call('load', 'sems-zlib/1.2.8/base'),
                         mock.call('load', 'sems-hdf5/1.10.6/parallel'),
                         mock.call('load', 'sems-netcdf/4.7.3/parallel'),
                         mock.call('load', 'sems-parmetis/4.0.3/parallel'),
                         mock.call('load', 'sems-scotch/6.0.3/nopthread_64bit_parallel'),
                         mock.call('load', 'sems-superlu/4.3/base'),
                         mock.call('load', 'sems-cmake/3.17.1'),
                         mock.call('load', 'sems-ninja_fortran/1.10.0'),
                         mock.call('load', 'atdm-env'),
                         ]
        self.buildEnv_passes(PR_name, expected_list,
                             test_ENV={'OMP_NUM_THREADS': '2'})


    def test_buildEnv_passes_with_intel_1701(self):
        """Find the function"""
        PR_name = 'Trilinos_pullrequest_intel_17.0.1'
        expected_list = [mock.call('use', '/projects/sems/modulefiles/projects'),
                         mock.call('load', 'sems-env'),
                         mock.call('load', 'sems-git/2.10.1'),
                         mock.call('load', 'sems-gcc/4.9.3'),
                         mock.call('load', 'sems-intel/17.0.1'),
                         mock.call('load', 'sems-mpich/3.2'),
                         mock.call('load', 'sems-python/2.7.9'),
                         mock.call('load', 'sems-boost/1.63.0/base'),
                         mock.call('load', 'sems-zlib/1.2.8/base'),
                         mock.call('load', 'sems-hdf5/1.10.6/parallel'),
                         mock.call('load', 'sems-netcdf/4.7.3/parallel'),
                         mock.call('load', 'sems-parmetis/4.0.3/parallel'),
                         mock.call('load', 'sems-scotch/6.0.3/nopthread_64bit_parallel'),
                         mock.call('load', 'sems-superlu/4.3/base'),
                         mock.call('load', 'sems-cmake/3.10.3'),
                         mock.call('load', 'atdm-env'),
                         mock.call('load', 'atdm-ninja_fortran/1.7.2'),
                         ]

        self.buildEnv_passes(PR_name, expected_list,
                             test_ENV={'OMP_NUM_THREADS': '2'})


    def test_buildEnv_passes_with_clang_701(self):
        """Find the function"""
        PR_name = 'Trilinos_pullrequest_clang_7.0.1'
        expected_list = [mock.call('use', '/projects/sems/modulefiles/projects'),
                         mock.call('load', 'sems-env'),
                         mock.call('load', 'sems-git/2.10.1'),
                         mock.call('load', 'sems-gcc/5.3.0'),
                         mock.call('load', 'sems-clang/7.0.1'),
                         mock.call('load', 'sems-openmpi/1.10.1'),
                         mock.call('load', 'sems-python/2.7.9'),
                         mock.call('load', 'sems-boost/1.63.0/base'),
                         mock.call('load', 'sems-zlib/1.2.8/base'),
                         mock.call('load', 'sems-hdf5/1.10.6/parallel'),
                         mock.call('load', 'sems-netcdf/4.7.3/parallel'),
                         mock.call('load', 'sems-parmetis/4.0.3/parallel'),
                         mock.call('load', 'sems-scotch/6.0.3/nopthread_64bit_parallel'),
                         mock.call('load', 'sems-superlu/4.3/base'),
                         mock.call('load', 'sems-cmake/3.12.2'),
                         mock.call('load', 'atdm-env'),
                         mock.call('load', 'atdm-ninja_fortran/1.7.2'),
                         ]

        self.buildEnv_passes(PR_name, expected_list,
                             test_ENV={'OMP_NUM_THREADS': '2'})


    def test_buildEnv_passes_with_clang_900(self):
        """Find the function"""
        PR_name = 'Trilinos_pullrequest_clang_9.0.0'
        expected_list = [mock.call('use', '/projects/sems/modulefiles/projects'),
                         mock.call('load', 'sems-env'),
                         mock.call('load', 'sems-git/2.10.1'),
                         mock.call('load', 'sems-gcc/5.3.0'),
                         mock.call('load', 'sems-clang/9.0.0'),
                         mock.call('load', 'sems-openmpi/1.10.1'),
                         mock.call('load', 'sems-python/2.7.9'),
                         mock.call('load', 'sems-boost/1.63.0/base'),
                         mock.call('load', 'sems-zlib/1.2.8/base'),
                         mock.call('load', 'sems-hdf5/1.10.6/parallel'),
                         mock.call('load', 'sems-netcdf/4.7.3/parallel'),
                         mock.call('load', 'sems-parmetis/4.0.3/parallel'),
                         mock.call('load', 'sems-scotch/6.0.3/nopthread_64bit_parallel'),
                         mock.call('load', 'sems-superlu/4.3/base'),
                         mock.call('load', 'sems-cmake/3.12.2'),
                         mock.call('load', 'atdm-env'),
                         mock.call('load', 'atdm-ninja_fortran/1.7.2'),
                         ]

        self.buildEnv_passes(PR_name, expected_list,
                             test_ENV={'OMP_NUM_THREADS': '2'})


    def test_buildEnv_passes_with_clang_1000(self):
        """Find the function"""
        PR_name = 'Trilinos_pullrequest_clang_10.0.0'
        expected_list = [mock.call('use', '/projects/sems/modulefiles/projects'),
                         mock.call('load', 'sems-env'),
                         mock.call('load', 'sems-git/2.10.1'),
                         mock.call('load', 'sems-gcc/5.3.0'),
                         mock.call('load', 'sems-clang/10.0.0'),
                         mock.call('load', 'sems-openmpi/1.10.1'),
                         mock.call('load', 'sems-python/2.7.9'),
                         mock.call('load', 'sems-boost/1.69.0/base'),
                         mock.call('load', 'sems-zlib/1.2.8/base'),
                         mock.call('load', 'sems-hdf5/1.10.6/parallel'),
                         mock.call('load', 'sems-netcdf/4.7.3/parallel'),
                         mock.call('load', 'sems-parmetis/4.0.3/parallel'),
                         mock.call('load', 'sems-scotch/6.0.3/nopthread_64bit_parallel'),
                         mock.call('load', 'sems-superlu/4.3/base'),
                         mock.call('load', 'sems-cmake/3.17.1'),
                         mock.call('load', 'sems-ninja_fortran/1.10.0'),
                         ]

        self.buildEnv_passes(PR_name, expected_list,
                             test_ENV={'OMP_NUM_THREADS': '2'})


    def test_buildEnv_passes_with_cuda_92(self):
        """Find the function"""
        PR_name = 'Trilinos_pullrequest_cuda_9.2'
        expected_list = [mock.call('load', 'git/2.10.1'),
                         mock.call('load', 'devpack/20180521/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88'),
                         mock.call('swap', 'openblas/0.2.20/gcc/7.2.0', 'netlib/3.8.0/gcc/7.2.0'),]
        expected_env = {'OMPI_CXX':
                       os.path.join(self.jenkins_workspace,
                                    'Trilinos',
                                    'packages',
                                    'kokkos',
                                    'bin',
                                    'nvcc_wrapper'),
                       'OMPI_CC': '/fake/gcc/path/bin/gcc',
                       'OMPI_FC': '/fake/gcc/path/bin/gfortran',
                       'CUDA_LAUNCH_BLOCKING': '1',
                       'CUDA_MANAGED_FORCE_DEVICE_ALLOC': '1',
                       'PATH': [os.path.join(os.path.sep,
                                            'ascldap',
                                            'users',
                                            'rabartl',
                                            'install',
                                            'white-ride',
                                            'cmake-3.11.2',
                                            'bin'),
                               os.path.join(os.path.sep,
                                            'ascldap',
                                            'users',
                                            'rabartl',
                                            'install',
                                            'white-ride',
                                            'ninja-1.8.2',
                                            'bin'),
                                '/fake/path']}

        self.buildEnv_passes(PR_name, expected_list, test_ENV=expected_env)



class Test_GetCDashTrack(unittest.TestCase):
    """use the default or override from environment"""

    def setUp(self):
        self.test_track = 'Experimental Test N'
        self.IOredirect = mock.patch('sys.stdout', new_callable=StringIO)


    def test_default(self):
        expected_output = 'PULLREQUEST_CDASH_TRACK isn\'t set, using default value\n'
        with self.IOredirect as m_output, \
            mock.patch.dict(os.environ,
                            {'NOT_PULLREQUEST_CDASH_TRACK': self.test_track},
                            clear=True):
            self.assertEqual('Pull Request',
                             PullRequestLinuxDriverTest.getCDashTrack())
        self.assertEqual(expected_output, m_output.getvalue())


    def test_FromEnvironment(self):
        expected_output = 'PULLREQUEST_CDASH_TRACK is set. Setting CDASH_TRACK={}\n'.format(
            self.test_track)
        with self.IOredirect as m_output, \
            mock.patch.dict(os.environ,
                            {'PULLREQUEST_CDASH_TRACK': self.test_track},
                            clear=True):
            self.assertEqual(self.test_track,
                             PullRequestLinuxDriverTest.getCDashTrack())
        self.assertEqual(expected_output, m_output.getvalue())



class testCompute_n(unittest.TestCase):
    '''How many processors will we use based on memory limits'''


    def setUp(self):
        self.m_environ = mock.patch.dict(os.environ,
                                         {'JENKINS_JOB_WEIGHT': '29'},
                                         clear=True)


    def test_over_weight(self):
        '''if the jenkins job weight is set higher than the  machines nprocs
        we should still  run at least one job'''
        m_open =  mock.mock_open(read_data='''MemTotal 32965846 kB''')
        with self.m_environ, \
            mock.patch('PullRequestLinuxDriverTest.open', m_open, create=True), \
            mock.patch('PullRequestLinuxDriverTest.cpu_count', return_value=18):
            parallel_level = PullRequestLinuxDriverTest.compute_n()
        self.assertEqual(10, parallel_level)


    def helper_call_compute_n(self, num_cpu, mem_kb, expected_cpu_count):
        """
        """
        mem_gb   = mem_kb / (1024**2)

        mem_info = {"mem_kb": mem_kb, "mem_gb": mem_gb}
        m_open   = mock.mock_open(read_data="MemTotal %s kB"%(mem_kb))

        with self.m_environ, \
            mock.patch('PullRequestLinuxDriverTest.open', m_open, create=True), \
            mock.patch('PullRequestLinuxDriverTest.cpu_count', return_value=num_cpu), \
            mock.patch('PullRequestLinuxDriverTest.get_memory_info', return_value=mem_info):
            parallel_level = PullRequestLinuxDriverTest.compute_n()
        self.assertEqual(expected_cpu_count, parallel_level)


    def test_32p_64g(self):
        '''check we match whats on 113-115 and the cloud'''
        num_cpu  = 32
        mem_kb   = 65805212
        expected_cpu_count = 20

        self.helper_call_compute_n(num_cpu, mem_kb, expected_cpu_count)


    def test_72p_64g(self):
        '''match whats on the 14x series'''
        num_cpu  = 72
        mem_kb   = 65805212
        expected_cpu_count = 10

        self.helper_call_compute_n(num_cpu, mem_kb, expected_cpu_count)


    def test_88p_128g(self):
        '''match ascic158'''
        num_cpu  = 88
        mem_kb   = 131610424
        expected_cpu_count = 13

        self.helper_call_compute_n(num_cpu, mem_kb, expected_cpu_count)


    def test_80p_128g(self):
        '''this matches ascic166'''
        num_cpu  = 80
        mem_kb   = 131610424
        expected_cpu_count = 20

        self.helper_call_compute_n(num_cpu, mem_kb, expected_cpu_count)



if __name__ == '__main__':
    unittest.main()  # pragma nocover
