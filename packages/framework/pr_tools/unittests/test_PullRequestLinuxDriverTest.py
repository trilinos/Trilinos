#!/usr/bin/env python3
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

import PullRequestLinuxDriverTest



class Test_parse_args(unittest.TestCase):
    """the argument parser"""

    def setUp(self):
        """
        """
        self.maxDiff = None

        self.stdoutRedirect = mock.patch('sys.stdout', new_callable=StringIO)
        self.stderrRedirect = mock.patch('sys.stderr', new_callable=StringIO)

        self.m_argv = mock.patch.object(sys, 'argv', ['programName',
                                                      '--source-repo-url',
                                                      os.path.join(os.path.sep,
                                                                   'dev',
                                                                   'null',
                                                                   'source_repo'),
                                                      '--source-branch-name', 'foobar',
                                                      '--target-repo-url',
                                                      os.path.join(os.path.sep,
                                                                   'dev',
                                                                   'null',
                                                                   'target_repo'),
                                                      '--target-branch-name', 'real_trash',
                                                      '--pullrequest-build-name',
                                                      'Some_odd_compiler',
                                                      '--genconfig-build-name',
                                                      'Some_odd_compiler_and_options',
                                                      '--pullrequest-number', '4242',
                                                      '--jenkins-job-number', '2424'])

        self.default_options = Namespace(source_repo_url='/dev/null/source_repo',
                                         source_branch_name='foobar',
                                         target_repo_url='/dev/null/target_repo',
                                         target_branch_name='real_trash',
                                         pullrequest_build_name='Some_odd_compiler',
                                         genconfig_build_name='Some_odd_compiler_and_options',
                                         pullrequest_number='4242',
                                         jenkins_job_number='2424',
                                         source_dir='UNKNOWN',
                                         build_dir='UNKNOWN',
                                         ctest_driver='UNKNOWN',
                                         ctest_drop_site='testing.sandia.gov',
                                         pullrequest_cdash_track='Pull Request',
                                         pullrequest_env_config_file='/dev/null/Trilinos_clone/pr_config/pullrequest.ini',
                                         pullrequest_gen_config_file='/dev/null/Trilinos_clone/../GenConfig/src/gen-config.ini',
                                         workspace_dir='/dev/null/Trilinos_clone',
                                         filename_packageenables='../packageEnables.cmake',
                                         filename_subprojects='../package_subproject_list.cmake',
                                         test_mode='standard',
                                         req_mem_per_core=3.0,
                                         max_cores_allowed=12,
                                         num_concurrent_tests=-1,
                                         ccache_enable=False,
                                         dry_run=False)

        self.default_stdout = dedent('''\
                | - [R] source-repo-url             : /dev/null/source_repo
                | - [R] source-branch-name          : foobar
                | - [R] target_repo_url             : /dev/null/target_repo
                | - [R] target_branch_name          : real_trash
                | - [R] pullrequest-build-name      : Some_odd_compiler
                | - [R] genconfig-build-name        : Some_odd_compiler_and_options
                | - [R] pullrequest-number          : 4242
                | - [R] jenkins-job-number          : 2424
                | - [R] source-dir                  : UNKNOWN
                | - [R] build-dir                   : UNKNOWN
                | - [R] ctest-driver                : UNKNOWN
                | - [R] ctest-drop-site             : testing.sandia.gov
                |
                | - [O] dry-run                     : False
                | - [O] enable-ccache               : False
                | - [O] filename-packageenables     : ../packageEnables.cmake
                | - [O] max-cores-allowed           : 12
                | - [O] num-concurrent-tests        : -1
                | - [O] pullrequest-cdash-track     : Pull Request
                | - [O] pullrequest-env-config-file : /dev/null/Trilinos_clone/pr_config/pullrequest.ini
                | - [O] pullrequest-gen-config-file : /dev/null/Trilinos_clone/../GenConfig/src/gen-config.ini
                | - [O] req-mem-per-core            : 3.0
                | - [O] test-mode                   : standard
                | - [O] workspace-dir               : /dev/null/Trilinos_clone
                ''')

        self.help_output = dedent('''\
                usage: programName [-h] --source-repo-url SOURCE_REPO_URL --source-branch-name
                                   SOURCE_BRANCH_NAME --target-repo-url TARGET_REPO_URL
                                   --target-branch-name TARGET_BRANCH_NAME
                                   --pullrequest-build-name PULLREQUEST_BUILD_NAME
                                   --genconfig-build-name GENCONFIG_BUILD_NAME
                                   --pullrequest-number PULLREQUEST_NUMBER
                                   --jenkins-job-number JENKINS_JOB_NUMBER
                                   [--source-dir SOURCE_DIR] [--build-dir BUILD_DIR]
                                   [--ctest-driver CTEST_DRIVER]
                                   [--ctest-drop-site CTEST_DROP_SITE]
                                   [--pullrequest-cdash-track PULLREQUEST_CDASH_TRACK]
                                   [--pullrequest-env-config-file PULLREQUEST_ENV_CONFIG_FILE]
                                   [--pullrequest-gen-config-file PULLREQUEST_GEN_CONFIG_FILE]
                                   [--workspace-dir WORKSPACE_DIR]
                                   [--filename-packageenables FILENAME_PACKAGEENABLES]
                                   [--filename-subprojects FILENAME_SUBPROJECTS]
                                   [--test-mode TEST_MODE]
                                   [--req-mem-per-core REQ_MEM_PER_CORE]
                                   [--max-cores-allowed MAX_CORES_ALLOWED]
                                   [--num-concurrent-tests NUM_CONCURRENT_TESTS]
                                   [--enable-ccache] [--dry-run]

                Parse the repo and build information

                optional arguments:
                  -h, --help            show this help message and exit

                Required Arguments:
                  --source-repo-url SOURCE_REPO_URL
                                        Repo with the new changes
                  --source-branch-name SOURCE_BRANCH_NAME
                                        Branch with the new changes
                  --target-repo-url TARGET_REPO_URL
                                        Repo to merge into
                  --target-branch-name TARGET_BRANCH_NAME
                                        Branch to merge into
                  --pullrequest-build-name PULLREQUEST_BUILD_NAME
                                        The Jenkins job base name
                  --genconfig-build-name GENCONFIG_BUILD_NAME
                                        The job base name for the cmake configuration
                  --pullrequest-number PULLREQUEST_NUMBER
                                        The github PR number
                  --jenkins-job-number JENKINS_JOB_NUMBER
                                        The Jenkins build number

                Optional Arguments:
                  --source-dir SOURCE_DIR
                                        Directory containing the source code to compile/test.
                  --build-dir BUILD_DIR
                                        Path to the build directory.
                  --ctest-driver CTEST_DRIVER
                                        Location of the CTest driver script to load via `-S`.
                  --ctest-drop-site CTEST_DROP_SITE
                                        URL of the cdash server to post to.
                  --pullrequest-cdash-track PULLREQUEST_CDASH_TRACK
                                        The CDash Track to add results to. Default=Pull
                                        Request
                  --pullrequest-env-config-file PULLREQUEST_ENV_CONFIG_FILE
                                        The Trilinos PR driver configuration file containing
                                        job mappings to environment specifications. Default=/d
                                        ev/null/Trilinos_clone/pr_config/pullrequest.ini
                  --pullrequest-gen-config-file PULLREQUEST_GEN_CONFIG_FILE
                                        The Trilinos PR driver configuration file containing
                                        job mappings to cmake specifications.
                                        Default=/dev/null/Trilinos_clone/../GenConfig/src/gen-
                                        config.ini
                  --workspace-dir WORKSPACE_DIR
                                        The local workspace directory that Jenkins set up.
                                        Default=/dev/null/Trilinos_clone
                  --filename-packageenables FILENAME_PACKAGEENABLES
                                        The packageEnables.cmake is usually generated by
                                        TriBiTS infrastructure based on which packages contain
                                        the changes between the source and target branches.
                                        Default=../packageEnables.cmake
                  --filename-subprojects FILENAME_SUBPROJECTS
                                        The subprojects_file is used by the testing
                                        infrastructure. This parameter allows the default,
                                        generated file, to be overridden. Generally this
                                        should not be changed from the defaults..
                                        Default=../package_subproject_list.cmake
                  --test-mode TEST_MODE
                                        PR testing mode. Use 'standard' for normal PR tests,
                                        'installation' for installation testing. Default =
                                        standard
                  --req-mem-per-core REQ_MEM_PER_CORE
                                        Minimum required memory per core (GB) to build
                                        Trilinos.Default = 3.0
                  --max-cores-allowed MAX_CORES_ALLOWED
                                        Max cores allowed, if >= 0 we will use the # of
                                        detected cores on the system. Default = 12
                  --num-concurrent-tests NUM_CONCURRENT_TESTS
                                        Set the number of concurrent tests allowd in CTest.
                                        This is equivalent to `ctest -j <num-concurrent-
                                        tests>`. If > 0 then this value is used, otherwise the
                                        value is calculated based on number_of_available_cores
                                        / max_test_parallelism Default = -1
                  --enable-ccache       Enable ccache object caching to improve build times.
                                        Default = False
                  --dry-run             Enable dry-run mode. Script will run but not execute
                                        the build steps. Default = False
                ''')

        self.usage_output = dedent('''\
                usage: programName [-h] --source-repo-url SOURCE_REPO_URL --source-branch-name
                                   SOURCE_BRANCH_NAME --target-repo-url TARGET_REPO_URL
                                   --target-branch-name TARGET_BRANCH_NAME
                                   --pullrequest-build-name PULLREQUEST_BUILD_NAME
                                   --genconfig-build-name GENCONFIG_BUILD_NAME
                                   --pullrequest-number PULLREQUEST_NUMBER
                                   --jenkins-job-number JENKINS_JOB_NUMBER
                                   [--source-dir SOURCE_DIR] [--build-dir BUILD_DIR]
                                   [--ctest-driver CTEST_DRIVER]
                                   [--ctest-drop-site CTEST_DROP_SITE]
                                   [--pullrequest-cdash-track PULLREQUEST_CDASH_TRACK]
                                   [--pullrequest-env-config-file PULLREQUEST_ENV_CONFIG_FILE]
                                   [--pullrequest-gen-config-file PULLREQUEST_GEN_CONFIG_FILE]
                                   [--workspace-dir WORKSPACE_DIR]
                                   [--filename-packageenables FILENAME_PACKAGEENABLES]
                                   [--filename-subprojects FILENAME_SUBPROJECTS]
                                   [--test-mode TEST_MODE]
                                   [--req-mem-per-core REQ_MEM_PER_CORE]
                                   [--max-cores-allowed MAX_CORES_ALLOWED]
                                   [--num-concurrent-tests NUM_CONCURRENT_TESTS]
                                   [--enable-ccache] [--dry-run]
                programName: error: the following arguments are required: --source-repo-url, --source-branch-name, --target-repo-url, --target-branch-name, --pullrequest-build-name, --genconfig-build-name, --pullrequest-number, --jenkins-job-number
                ''')

        self.m_cwd = mock.patch('PullRequestLinuxDriverTest.os.getcwd',
                                return_value = os.path.join(os.path.sep,
                                                           'dev',
                                                           'null',
                                                           'Trilinos_clone'))
        self.m_cwd.start()
        return


    def tearDown(self):
        self.m_cwd.stop()


    def test_parse_args_returns_defaults(self):
        '''
        No inputs
        '''
        import difflib
        with self.m_argv, self.stdoutRedirect as m_stdout:
            returned_default = PullRequestLinuxDriverTest.parse_args()

        self.assertEqual(self.default_options, returned_default)

        self.assertIn(self.default_stdout, m_stdout.getvalue())
        return


    def test_parse_args_uses_workspace_environ(self):
        '''No inputs'''
        tmp_workspace_dir = os.path.join(os.path.sep, 'dev', 'null', 'Trilinos_workspace')

        l_options               = self.default_options
        l_options.workspace_dir = tmp_workspace_dir

        l_stdout = self.default_stdout
        l_stdout = l_stdout.replace('workspace-dir               : /dev/null/Trilinos_clone',
                                    'workspace-dir               : /dev/null/Trilinos_workspace')

        with self.m_argv, mock.patch.dict(os.environ, {'WORKSPACE': tmp_workspace_dir}, clear=True), self.stdoutRedirect as m_stdout:
            returned_default = PullRequestLinuxDriverTest.parse_args()

        self.assertEqual(l_options, returned_default)
        self.assertIn(l_stdout, m_stdout.getvalue())
        return


    def test_help(self):
        """
        Compare the help message to the expected
        """

        with mock.patch.object(sys, 'argv', ['programName', '--help']), self.assertRaises(SystemExit), self.stdoutRedirect as m_stdout:
            PullRequestLinuxDriverTest.parse_args()

        self.assertEqual(self.help_output, m_stdout.getvalue())
        return


    def test_usage(self):
        '''
        Compare the usage message to the expected
        '''

        with mock.patch.object(sys, 'argv', ['programName', '--usage']), self.assertRaises(SystemExit), self.stderrRedirect as m_stderr:
            PullRequestLinuxDriverTest.parse_args()

        self.assertIn(self.usage_output, m_stderr.getvalue())
        return


class Test_main(unittest.TestCase):
    """the run script"""

    def setUp(self):
        self.prstd_mock = \
            mock.patch('PullRequestLinuxDriverTest.trilinosprhelpers.'
                       'TrilinosPRConfigurationStandard')
        self.prinstall_mock = \
            mock.patch('PullRequestLinuxDriverTest.trilinosprhelpers.'
                       'TrilinosPRConfigurationInstallation')

        self.IOredirect = mock.patch('sys.stdout', new_callable=StringIO)


    def test_calls_standard(self):
        """
        If the test_mode is standard
        """
        m_args = Namespace()
        setattr(m_args, 'dry_run', False)
        setattr(m_args, 'test_mode', 'standard')

        with self.prstd_mock as m_pr_helper, self.IOredirect as m_io:
            self.assertTrue(PullRequestLinuxDriverTest.main(m_args))

        expected_calls_list = [mock.call.prepare_test(), mock.call.execute_test()]

        m_pr_helper.assert_called_once_with(m_args)
        m_pr_helper.return_value.assert_has_calls(expected_calls_list)

        self.assertIn("P U L L R E Q U E S T   D R I V E R   S T A R T", m_io.getvalue())
        return


    def test_calls_standard_dry_run(self):
        """
        If the test_mode is standard
        """
        m_args = Namespace()
        setattr(m_args, 'dry_run', True)
        setattr(m_args, 'test_mode', 'standard')

        with self.prstd_mock as m_pr_helper, self.IOredirect as m_io:
            self.assertTrue(PullRequestLinuxDriverTest.main(m_args))

        expected_calls_list = [mock.call.prepare_test(), mock.call.execute_test()]

        m_pr_helper.assert_called_once_with(m_args)
        m_pr_helper.return_value.assert_has_calls(expected_calls_list)

        self.assertIn("P U L L R E Q U E S T   D R I V E R   S T A R T", m_io.getvalue())
        self.assertIn("D R Y   R U N   M O D E   E N A B L E D", m_io.getvalue())

        return


    def test_calls_install(self):
        """
        If the test_mode is installation
        """
        m_args = Namespace()
        setattr(m_args, 'dry_run', False)
        setattr(m_args, 'test_mode', 'installation')

        with self.prinstall_mock as m_pr_helper, self.IOredirect as m_io:
            self.assertTrue(PullRequestLinuxDriverTest.main(m_args))

        expected_calls_list = [ mock.call.prepare_test(), mock.call.execute_test() ]

        m_pr_helper.assert_called_once_with(m_args)
        m_pr_helper.return_value.assert_has_calls(expected_calls_list)

        self.assertIn("P U L L R E Q U E S T   D R I V E R   S T A R T", m_io.getvalue())

        return


    def test_calls_install_dry_run(self):
        """If the test_mode is installation"""
        m_args = Namespace()
        setattr(m_args, 'dry_run', True)
        setattr(m_args, 'test_mode', 'installation')

        with self.prinstall_mock as m_pr_helper, self.IOredirect as m_io:
            self.assertTrue(PullRequestLinuxDriverTest.main(m_args))

        expected_calls_list = [ mock.call.prepare_test(), mock.call.execute_test() ]

        m_pr_helper.assert_called_once_with(m_args)
        m_pr_helper.return_value.assert_has_calls( expected_calls_list )

        self.assertIn("P U L L R E Q U E S T   D R I V E R   S T A R T", m_io.getvalue())
        self.assertIn("D R Y   R U N   M O D E   E N A B L E D", m_io.getvalue())

        return


    def test_raises_on_unkown(self):
        """If the test_mode is installation"""
        m_args = Namespace()
        setattr(m_args, 'dry_run', False)
        setattr(m_args, 'test_mode', 'foobar')

        with self.prinstall_mock as m_pr_helper, \
             self.IOredirect as m_io, \
             self.assertRaisesRegex(KeyError,
                                    'ERROR: Unknown test mode, foobar, was provided.'):
            PullRequestLinuxDriverTest.main(m_args)

        m_pr_helper.assert_not_called()
        self.assertEqual(dedent('''\
                                +==============================================================================+
                                |
                                |   T R I L I N O S   P U L L R E Q U E S T   D R I V E R   S T A R T
                                |
                                +==============================================================================+

'''),
                                m_io.getvalue())
        return


    def test_always_fail(self):
        """
        This test exists to make an easy 'always fail' test for debugging changes
        to our testing framework when we want to verify that the framework will catch
        a failing python unit test.
        """
        force_test_to_fail = True

        # uncomment this line to make this test pass
        # This should never be commented out in an actual commit
        force_test_to_fail = False

        # if force_test_to_fail is not False then we should fail this test.
        self.assertEqual(force_test_to_fail, False)
        return


if __name__ == '__main__':
    unittest.main()  # pragma nocover
