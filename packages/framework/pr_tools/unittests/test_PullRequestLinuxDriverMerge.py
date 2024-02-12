#!/usr/bin/env python3
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
'''
Tests for the Merge chunk of the Driver script
'''
from __future__ import print_function
import sys
sys.dont_write_bytecode = True

import os
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


import unittest
from unittest.mock import patch

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

try:
    import unittest.mock as mock
except ImportError:  # pragma nocover
    import mock

from argparse import Namespace
import contextlib
from subprocess import CalledProcessError

import PullRequestLinuxDriverMerge as PRMerge



class Test_header(unittest.TestCase):
    '''Test that we can properly echo the header information'''

    def test_writeHeader(self):
        print("BEGIN:  Test PullRequestLinuxDriverMerge.py :: `write_header()`:")

        #m_stdout = StringIO()
        #with contextlib.redirect_stdout(m_stdout):
        with mock.patch('sys.stdout', new_callable=StringIO) as m_stdout:
            PRMerge.write_header()
            self.assertIn("Begin: PullRequestLinuxDriverMerge.py", m_stdout.getvalue())

        print("FINISH: Test PullRequestLinuxDriverMerge.py :: `write_header()`:")
        return



class Test_EchoJenkinsVars(unittest.TestCase):
    '''Test that the Jenkins environment is echoed properly'''

    def setUp(self):
        tmp_environ = {}
        tmp_environ['JOB_BASE_NAME'] = 'TEST_JOB_BASE_NAME'
        tmp_environ['JOB_NAME']      = 'TEST_JOB_NAME'
        tmp_environ['WORKSPACE']     = os.path.join(os.sep, 'dev', 'null', 'TEST_WORKSPACE')
        tmp_environ['NODE_NAME']     = 'TEST_NODE_NAME'
        self.m_environ = mock.patch.dict(os.environ, tmp_environ, clear=True)
        return


    def test_echoJenkinsVars(self):
        print("BEGIN:  Test PullRequestLinuxDriverMerge.py :: `echoJenkinsVars()`:")
        with self.m_environ:
            env_string_io = StringIO()
            for key in os.environ:
                print(key + ' = ' + os.environ[key],
                      file=env_string_io)

        tmp_path = os.path.join(os.sep, 'dev', 'null', 'TEST_WORKSPACE')

        with mock.patch('sys.stdout', new_callable=StringIO) as m_stdout, self.m_environ:
            PRMerge.echoJenkinsVars( tmp_path )

        stdout_actual = m_stdout.getvalue()
        self.assertIn("JOB_BASE_NAME = TEST_JOB_BASE_NAME", stdout_actual)
        self.assertIn("JOB_NAME = TEST_JOB_NAME", stdout_actual)
        self.assertIn("WORKSPACE = /dev/null/TEST_WORKSPACE", stdout_actual)
        self.assertIn("NODE_NAME = TEST_NODE_NAME", stdout_actual)
        print("FINISH: Test PullRequestLinuxDriverMerge.py :: `echoJenkinsVars()`:")
        return



class Test_parsing(unittest.TestCase):
    '''Finally given a decent parser I want to now pass
       all parameters and check them'''

    def test_parseArgs(self):
        test_namespace = Namespace()
        setattr(test_namespace, 'sourceRepo', '/dev/null/source/Trilinos.git')
        setattr(test_namespace, 'targetRepo', '/dev/null/target/Trilinos.git')
        setattr(test_namespace, 'targetBranch', 'fake_develop')
        setattr(test_namespace, 'sourceSHA', '0123456789abcdef')
        setattr(test_namespace, 'workspaceDir', '/dev/null/workspace')

        with mock.patch.object(sys, 'argv', ['programName',
                                             os.path.join(os.path.sep,
                                                          'dev',
                                                          'null',
                                                          'source',
                                                          'Trilinos.git'),
                                             os.path.join(os.path.sep,
                                                          'dev',
                                                          'null',
                                                          'target',
                                                          'Trilinos.git'),
                                             'fake_develop',
                                             '0123456789abcdef',
                                             os.path.join(os.path.sep,
                                                          'dev',
                                                          'null',
                                                          'workspace')]):

            args = PRMerge.parseArgs()
        self.assertEqual(test_namespace, args)
        return


    def test_parseInsufficientArgs_fails(self):
        test_namespace = Namespace()
        setattr(test_namespace, 'sourceRepo', '/dev/null/source/Trilinos.git')
        expected_output = '''usage: programName [-h]
                   sourceRepo targetRepo targetBranch sourceSHA workspaceDir
programName: error: the following arguments are required: sourceRepo, \
targetRepo, targetBranch, sourceSHA, workspaceDir
'''
        if sys.version_info.major != 3:
            expected_output = '''usage: programName [-h]
                   sourceRepo targetRepo targetBranch sourceSHA
                   workspaceDir\nprogramName: error: too few arguments
'''
        with mock.patch.object(sys, 'argv', ['programName']), \
                mock.patch('sys.stderr', new_callable=StringIO) as m_stderr:
            if sys.version_info.major != 3:
                with self.assertRaisesRegexp(SystemExit, '2'):
                    PRMerge.parseArgs()
            else:
                with self.assertRaisesRegex(SystemExit, '2'):
                    PRMerge.parseArgs()
        self.assertEqual(expected_output, m_stderr.getvalue())
        return



class Test_mergeBranch(unittest.TestCase):
    '''Verify that we call the correct sequence to merge the source branch/SHA
       into the target branch'''

    def test_mergeBranch_without_source_remote(self):
        with mock.patch('subprocess.check_output', side_effect=['origin /dev/null/target/Trilinos', '']) as m_check_out, \
            mock.patch('subprocess.check_call') as m_check_call:
            PRMerge.merge_branch(os.path.join(os.path.sep,
                                                                  'dev',
                                                                  'null',
                                                                  'source',
                                                                  'Trilinos.git'),
                                                    'fake_develop',
                                                    'df324ae')
        m_check_out.assert_has_calls([mock.call(['git', 'remote', '-v'])])

        m_check_call.assert_has_calls([mock.call(['git', 'remote', 'add', 'source_remote', '/dev/null/source/Trilinos.git']),
                                       mock.call(['git', 'prune']),
                                       mock.call(['git', 'gc']),
                                       mock.call(['git', 'fetch', 'source_remote', 'df324ae']),
                                       mock.call(['git', 'fetch', 'origin', 'fake_develop']),
                                       mock.call(['git', 'reset', '--hard', 'HEAD']),
                                       mock.call(['git', 'checkout', '-B', 'fake_develop', 'origin/fake_develop']),
                                       mock.call(['git', 'merge', '--no-ff', '--no-edit', 'df324ae']),
                                       ])
        return


    @patch('subprocess.check_call')
    def test_mergeBranch_with_source_remote(self, m_check_call):
        """
        """
        tmp_path = os.path.join(os.path.sep, 'dev', 'null', 'source', 'Trilinos.git')

        side_effect_list = ["origin /dev/null/target/Trilinos\nsource_remote /dev/null/source12/Trilinos.git", '']

        with mock.patch('subprocess.check_output', side_effect=side_effect_list) as m_check_out:
            with mock.patch('sys.stdout', new_callable=StringIO) as m_stdout:
                PRMerge.merge_branch(tmp_path, 'fake_develop', 'df324ae')

        expected_calls = []
        expected_calls.append(mock.call(['git', 'remote', '-v']))
        expected_calls.append(mock.call('git rev-parse --verify --quiet source_remote/df324ae || true', shell=True))
        m_check_out.assert_has_calls(expected_calls)

        m_check_call.assert_has_calls([mock.call(['git', 'remote', 'rm', 'source_remote']),
                                       mock.call(['git', 'remote', 'add', 'source_remote', '/dev/null/source/Trilinos.git']),
                                       mock.call(['git', 'prune']),
                                       mock.call(['git', 'gc']),
                                       mock.call(['git', 'fetch', 'source_remote', 'df324ae']),
                                       mock.call(['git', 'fetch', 'origin', 'fake_develop']),
                                       mock.call(['git', 'reset', '--hard', 'HEAD']),
                                       mock.call(['git', 'checkout', '-B', 'fake_develop', 'origin/fake_develop']),
                                       mock.call(['git', 'merge', '--no-ff', '--no-edit', 'df324ae']),
                                       ])
        self.assertIn("git remote exists, removing it", m_stdout.getvalue())
        return


    @patch('subprocess.check_call')
    def test_mergeBranch_ref_is_remote_branch(self, m_check_call):
        """
        """
        tmp_path = os.path.join(os.path.sep, 'dev', 'null', 'source', 'Trilinos.git')

        side_effect_list = ["origin /dev/null/target/Trilinos\nsource_remote /dev/null/source12/Trilinos.git", 'df324ae']

        with mock.patch('subprocess.check_output', side_effect=side_effect_list) as m_check_out:
            with mock.patch('sys.stdout', new_callable=StringIO) as m_stdout:
                PRMerge.merge_branch(tmp_path, 'fake_develop', 'some_ref')

        expected_calls = []
        expected_calls.append( mock.call(['git', 'remote', '-v']) )
        m_check_out.assert_has_calls(expected_calls)

        m_check_call.assert_has_calls([mock.call(['git', 'remote', 'rm', 'source_remote']),
                                       mock.call(['git', 'remote', 'add', 'source_remote', '/dev/null/source/Trilinos.git']),
                                       mock.call(['git', 'prune']),
                                       mock.call(['git', 'gc']),
                                       mock.call(['git', 'fetch', 'source_remote', 'some_ref']),
                                       mock.call(['git', 'fetch', 'origin', 'fake_develop']),
                                       mock.call(['git', 'reset', '--hard', 'HEAD']),
                                       mock.call(['git', 'checkout', '-B', 'fake_develop', 'origin/fake_develop']),
                                       mock.call(['git', 'merge', '--no-ff', '--no-edit', 'source_remote/some_ref']),
                                       ])
        self.assertIn("git remote exists, removing it", m_stdout.getvalue())
        return


    @patch('subprocess.check_output', side_effect=['origin /dev/null/target/Trilinos', 'df324ae'])
    def test_mergeBranch_fails_on_source_fetch(self, m_check_out):
        """
        Tests that ``merge_branch`` should fail on source fetch if a failure
        occurrs in the *fetch* operation
        """
        tmp_path = os.path.join(os.path.sep, 'dev', 'null', 'source', 'Trilinos.git')

        side_effect_list = [
            None,
            None,
            None,
            CalledProcessError(-1, 'cmd'),
            CalledProcessError(-2, 'cmd'),
            CalledProcessError(-3, 'cmd')
        ]

        with self.assertRaises(SystemExit):
            with mock.patch('subprocess.check_call', side_effect=side_effect_list) as m_check_call:
                PRMerge.merge_branch(tmp_path, 'fake_develop', 'df324ae')

        expected_calls = []
        expected_calls.append(mock.call(['git', 'remote', '-v']))
        m_check_out.assert_has_calls(expected_calls)

        expected_calls = []
        expected_calls.append(mock.call(['git', 'remote', 'add', 'source_remote', '/dev/null/source/Trilinos.git']))
        expected_calls.append(mock.call(['git', 'prune']))
        expected_calls.append(mock.call(['git', 'gc']))
        expected_calls.append(mock.call(['git', 'fetch', 'source_remote', 'df324ae']))
        expected_calls.append(mock.call(['git', 'fetch', 'source_remote', 'df324ae']))
        expected_calls.append(mock.call(['git', 'fetch', 'source_remote', 'df324ae']))
        m_check_call.assert_has_calls(expected_calls)

        return



class Test_run(unittest.TestCase):
    '''This is the main function that ties everything together in order'''

    @patch('PullRequestLinuxDriverMerge.parseArgs')
    @patch('os.chdir')
    @patch('PullRequestLinuxDriverMerge.write_header')
    @patch('PullRequestLinuxDriverMerge.echoJenkinsVars')
    @patch('PullRequestLinuxDriverMerge.merge_branch')
    @patch('os.path.join')
    def test_run(self, m_join, m_mergeBranch, m_echoJenkins, m_writeHeader, m_chdir, m_parser):
        """
        TODO: Add docstring that explains this test's purpose
        """
        with mock.patch('sys.stdout', new_callable=StringIO) as m_stdout:
            self.assertTrue(PRMerge.run())

        m_parser.assert_called_once_with()
        m_chdir.assert_called_once_with(m_join(m_parser().workspaceDir, 'Trilinos'))
        m_writeHeader.assert_called_once_with()
        m_echoJenkins.assert_called_once_with(m_parser().workspaceDir)
        m_mergeBranch.assert_called_once_with(m_parser().sourceRepo,
                                              m_parser().targetBranch,
                                              m_parser().sourceSHA)

        self.assertIn('Set CWD = ' + str(m_join()) + '\n', m_stdout.getvalue())
        return


    def test_run_fails_on_bad_parse(self):
        with mock.patch('PullRequestLinuxDriverMerge.parseArgs', side_effect=SystemExit(2)):
            self.assertFalse(PRMerge.run())
        return


    @patch('os.path.join')
    @patch('os.chdir')
    @patch('PullRequestLinuxDriverMerge.parseArgs')
    def test_run_fails_on_bad_fetch(self, m_parseArgs, m_chdir, m_join):
        """
        TODO: Add a docstring explaining this test and its purpose.
        """
        side_effect_check_output = ['origin /dev/null/target/Trilinos', 'df324ae']

        side_effect_check_call = []
        side_effect_check_call.append( None )
        side_effect_check_call.append( CalledProcessError(-1, 'cmd') )
        side_effect_check_call.append( CalledProcessError(-2, 'cmd') )
        side_effect_check_call.append( CalledProcessError(-3, 'cmd') )

        with mock.patch('subprocess.check_output', side_effect=side_effect_check_output):
            with mock.patch('subprocess.check_call', side_effect=side_effect_check_call):
                # I'm not sure if we really needed to mock out stdio, we don't check it
                # so would it be more useful to provide that to the output if the test
                # should fail? Otherwise, are we flying blin?
                with mock.patch('sys.stdout', new_callable=StringIO):
                    self.assertFalse(PRMerge.run())

        return


    @patch('os.path.join')
    @patch('os.chdir')
    @patch('PullRequestLinuxDriverMerge.parseArgs')
    @patch('sys.stdout', new_callable=StringIO)
    def test_run_fails_on_bad_remote_add(self, m_stdout, m_parseArgs, m_chdir, m_path_join):

        expected_string = '''Recieved subprocess.CalledProcessError - returned -1
  from command test_cmd
  output None
  stdout None
  stderr None\n'''

        if sys.version_info.major != 3:
            expected_string = '''Recieved subprocess.CalledProcessError - returned -1
  from command test_cmd
  output None\n'''

        side_effect_check_output = ['origin /dev/null/target/Trilinos', 'df324ae']
        side_effect_check_call   = CalledProcessError(-1, 'test_cmd')

        with mock.patch('subprocess.check_output', side_effect=side_effect_check_output), \
            mock.patch('subprocess.check_call', side_effect=side_effect_check_call):
            self.assertFalse(PRMerge.run())

        expected_string_list = []
        expected_string_list.append("Received subprocess.CalledProcessError - returned -1")
        expected_string_list.append("from command test_cmd")
        expected_string_list.append("output None")

        if sys.version_info.major == 3:
            expected_string_list.append("output None")
            expected_string_list.append("stdout None")
            expected_string_list.append("stderr None")

        for expected_string_i in expected_string_list:
            self.assertIn(expected_string_i, m_stdout.getvalue())

        return


if __name__ == '__main__':
    unittest.main()  # pragma nocover
