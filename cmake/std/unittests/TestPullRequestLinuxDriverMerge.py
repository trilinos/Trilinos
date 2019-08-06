#!/usr/bin/env python
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
'''
Tests for the Merge chunk of the Driver script
'''
from __future__ import print_function
import sys
sys.dont_write_bytecode = True

import os
sys.path.insert(1, os.path.dirname(os.path.dirname(__file__)))


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

import PullRequestLinuxDriverMerge


class Test_header(unittest.TestCase):
    '''Test that we can properly echo the header information'''
    def test_writeHeader(self):
        with mock.patch('sys.stdout', new_callable=StringIO) as m_stdout:
            PullRequestLinuxDriverMerge.write_header()
        self.assertEqual('''--------------------------------------------------------------------------------
-
- Begin: PullRequestLinuxDriver-Merge.py
-
--------------------------------------------------------------------------------
''',
                         m_stdout.getvalue())


class Test_EchoJenkinsVars(unittest.TestCase):
    '''Test that the Jenkins environment is echoed properly'''


    def setUp(self):
        self.m_environ = mock.patch.dict(os.environ, {'JOB_BASE_NAME':'TEST_JOB_BASE_NAME',
                                         'JOB_NAME':'TEST_JOB_NAME',
                                         'WORKSPACE':os.path.join(os.sep,
                                                                  'dev',
                                                                  'null',
                                                                  'TEST_WORKSPACE'),
                                         'NODE_NAME':'TEST_NODE_NAME'},
                            clear=True)


    def test_echoJenkinsVars(self):
        with self.m_environ:
            env_string_io = StringIO()
            for key in os.environ:
                print(key + ' = ' + os.environ[key],
                      file=env_string_io)

        expected_string = '''
================================================================================
Jenkins Environment Variables:
- WORKSPACE    : /dev/null/TEST_WORKSPACE

================================================================================
Environment:

  pwd = {cwd}

{environ}
================================================================================
'''.format(cwd=os.getcwd(),
           environ=env_string_io.getvalue())

        with mock.patch('sys.stdout', new_callable=StringIO) as m_stdout, \
            self.m_environ:
            PullRequestLinuxDriverMerge.echoJenkinsVars(os.path.join(os.sep,
                                                                     'dev',
                                                                     'null',
                                                                     'TEST_WORKSPACE'))
        self.assertEqual(expected_string, m_stdout.getvalue())


class Test_parsing(unittest.TestCase):
    '''Finally given a decent parser I want to now pass
       all parameters and check them'''

    def test_parseArgs(self):
        test_namespace = Namespace()
        setattr(test_namespace, 'sourceRepo', '/dev/null/source/Trilinos.git')
        setattr(test_namespace, 'sourceBranch', 'fake_test_branch_fixing_issue_neverland')
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
                                             'fake_test_branch_fixing_issue_neverland',
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

            args = PullRequestLinuxDriverMerge.parseArgs()
        self.assertEqual(test_namespace, args)

    def test_parseInsufficientArgs_fails(self):
        test_namespace = Namespace()
        setattr(test_namespace, 'sourceRepo', '/dev/null/source/Trilinos.git')
        expected_output = '''usage: programName [-h]
                   sourceRepo sourceBranch targetRepo targetBranch sourceSHA
                   workspaceDir
programName: error: the following arguments are required: sourceRepo, \
sourceBranch, targetRepo, targetBranch, sourceSHA, workspaceDir
'''
        if sys.version_info.major is not 3:
            expected_output = '''usage: programName [-h]
                   sourceRepo sourceBranch targetRepo targetBranch sourceSHA
                   workspaceDir\nprogramName: error: too few arguments
'''
        with mock.patch.object(sys, 'argv', ['programName']), \
                mock.patch('sys.stderr', new_callable=StringIO) as m_stderr:
            if sys.version_info.major is not 3:
                with self.assertRaisesRegexp(SystemExit, '2'):
                    PullRequestLinuxDriverMerge.parseArgs()
            else:
                with self.assertRaisesRegex(SystemExit, '2'):
                    PullRequestLinuxDriverMerge.parseArgs()
        self.assertEqual(expected_output, m_stderr.getvalue())


class Test_mergeBranch(unittest.TestCase):
    '''Verify that we call the correct sequence to merge the source branch/SHA
       into the target branch'''
    def test_mergeBranch_without_source_remote(self):
        with mock.patch('subprocess.check_output',
                        side_effect=['origin /dev/null/target/Trilinos',
                                      'df324ae']) as m_check_out, \
            mock.patch('subprocess.check_call') as m_check_call:
            PullRequestLinuxDriverMerge.merge_branch(os.path.join(os.path.sep,
                                                                  'dev',
                                                                  'null',
                                                                  'source',
                                                                  'Trilinos.git'),
                                                    'neverland',
                                                    'fake_develop',
                                                    'df324ae')
        m_check_out.assert_has_calls([mock.call(['git', 'remote', '-v']),
                                      mock.call(['git', 'rev-parse',
                                                 'source_remote/neverland'])])

        m_check_call.assert_has_calls([mock.call(['git', 'remote', 'add',
                                                 'source_remote',
                                                 '/dev/null/source/Trilinos.git']),
                                       mock.call(['git', 'fetch', 'source_remote',
                                                  'neverland']),
                                       mock.call(['git', 'fetch', 'origin',
                                                  'fake_develop']),
                                       mock.call(['git', 'reset', '--hard',
                                                  'HEAD']),
                                       mock.call(['git', 'checkout',
                                                  '-B', 'fake_develop',
                                                  'origin/fake_develop']),
                                       mock.call(['git', 'merge',
                                                  '--no-edit',
                                                  'source_remote/neverland']),
                                       ])


    def test_mergeBranch_with_source_remote(self):
        with mock.patch('subprocess.check_output',
                        side_effect=['''origin /dev/null/target/Trilinos
    source_remote /dev/null/source12/Trilinos.git''', 'df324ae']) as m_check_out, \
            mock.patch('subprocess.check_call') as m_check_call, \
            mock.patch('sys.stdout', new_callable=StringIO) as m_stdout:
            PullRequestLinuxDriverMerge.merge_branch(os.path.join(os.path.sep,
                                                                  'dev',
                                                                  'null',
                                                                  'source',
                                                                  'Trilinos.git'),
                                                    'neverland',
                                                    'fake_develop',
                                                    'df324ae')
        m_check_out.assert_has_calls([mock.call(['git', 'remote', '-v']),
                                      mock.call(['git', 'rev-parse',
                                                'source_remote/neverland'])])

        m_check_call.assert_has_calls([mock.call(['git', 'remote', 'rm',
                                                 'source_remote']),
                                       mock.call(['git', 'remote', 'add',
                                                  'source_remote',
                                                  '/dev/null/source/Trilinos.git']),
                                       mock.call(['git', 'fetch', 'source_remote',
                                                  'neverland']),
                                       mock.call(['git', 'fetch', 'origin',
                                                  'fake_develop']),
                                       mock.call(['git', 'reset', '--hard',
                                                  'HEAD']),
                                       mock.call(['git', 'checkout',
                                                  '-B', 'fake_develop',
                                                  'origin/fake_develop']),
                                       mock.call(['git', 'merge',
                                                  '--no-edit',
                                                  'source_remote/neverland']),
                                       ])
        self.assertEqual("git remote exists, removing it\n", m_stdout.getvalue())


    def test_mergeBranch_fails_on_source_fetch(self):
        with mock.patch('subprocess.check_output',
                        side_effect=['origin /dev/null/target/Trilinos',
                                      'df324ae']) as m_check_out, \
            mock.patch('subprocess.check_call',
                       side_effect=[None,
                                    CalledProcessError(-1, 'cmd'),
                                    CalledProcessError(-2, 'cmd'),
                                    CalledProcessError(-3, 'cmd')]) as m_check_call:
            if sys.version_info.major is not 3:
                with self.assertRaisesRegexp(SystemExit, '12'):
                    PullRequestLinuxDriverMerge.merge_branch(os.path.join(os.path.sep,
                                                                          'dev',
                                                                          'null',
                                                                          'source',
                                                                          'Trilinos.git'),
                                                            'neverland',
                                                            'fake_develop',
                                                            'df324ae')
            else:
                with self.assertRaisesRegex(SystemExit, '12'):
                    PullRequestLinuxDriverMerge.merge_branch(os.path.join(os.path.sep,
                                                                          'dev',
                                                                          'null',
                                                                          'source',
                                                                          'Trilinos.git'),
                                                             'neverland',
                                                             'fake_develop',
                                                             'df324ae')
        m_check_out.assert_has_calls([mock.call(['git', 'remote', '-v'])])

        m_check_call.assert_has_calls([mock.call(['git', 'remote', 'add',
                                                 'source_remote',
                                                 '/dev/null/source/Trilinos.git']),
                                       mock.call(['git', 'fetch', 'source_remote',
                                                  'neverland']),
                                       mock.call(['git', 'fetch', 'source_remote',
                                                  'neverland']),
                                       mock.call(['git', 'fetch', 'source_remote',
                                                  'neverland']),
                                       ])


    def test_mergeBranch_fails_on_SHA_mismatch(self):
        with mock.patch('subprocess.check_output',
                        side_effect=['origin /dev/null/target/Trilinos',
                                      'df324ae']) as m_check_out, \
            mock.patch('subprocess.check_call') as m_check_call, \
            mock.patch('sys.stdout', new_callable=StringIO) as m_stdout:
            if sys.version_info.major is not 3:
                with self.assertRaisesRegexp(SystemExit, '-1'):
                    PullRequestLinuxDriverMerge.merge_branch(os.path.join(os.path.sep,
                                                                          'dev',
                                                                          'null',
                                                                          'source',
                                                                          'Trilinos.git'),
                                                            'neverland',
                                                            'fake_develop',
                                                            'foobar')
            else:
                with self.assertRaisesRegex(SystemExit, '-1'):
                    PullRequestLinuxDriverMerge.merge_branch(os.path.join(os.path.sep,
                                                                          'dev',
                                                                          'null',
                                                                          'source',
                                                                          'Trilinos.git'),
                                                             'neverland',
                                                             'fake_develop',
                                                             'foobar')
        m_check_out.assert_has_calls([mock.call(['git', 'remote', '-v']),
                                      mock.call(['git', 'rev-parse',
                                                 'source_remote/neverland'])])

        m_check_call.assert_has_calls([mock.call(['git', 'remote', 'add',
                                                 'source_remote',
                                                 '/dev/null/source/Trilinos.git']),
                                       mock.call(['git', 'fetch', 'source_remote',
                                                  'neverland']),
                                       ])
        self.assertEqual('''The SHA (df324ae) for the last commit on branch neverland
  in repo /dev/null/source/Trilinos.git is different than the expected SHA,
  which is: foobar.\n''',
                         m_stdout.getvalue())


class Test_run(unittest.TestCase):
    '''This is the main function that ties everything together in order'''
    def test_run(self):
        with mock.patch('PullRequestLinuxDriverMerge.parseArgs') as m_parser, \
            mock.patch('os.chdir') as m_chdir, \
            mock.patch('PullRequestLinuxDriverMerge.write_header') as m_writeHeader, \
            mock.patch('PullRequestLinuxDriverMerge.echoJenkinsVars') as m_echoJenkins, \
            mock.patch('PullRequestLinuxDriverMerge.merge_branch') as m_mergeBranch, \
            mock.patch('os.path.join') as m_join, \
            mock.patch('sys.stdout', new_callable=StringIO) as m_stdout:
            self.assertTrue(PullRequestLinuxDriverMerge.run())
        m_parser.assert_called_once_with()
        m_chdir.assert_called_once_with(m_join(m_parser().workspaceDir, 'Trilinos'))
        m_writeHeader.assert_called_once_with()
        m_echoJenkins.assert_called_once_with(m_parser().workspaceDir)
        m_mergeBranch.assert_called_once_with(m_parser().sourceRepo,
                                              m_parser().sourceBranch,
                                              m_parser().targetBranch,
                                              m_parser().sourceSHA)
        self.assertEqual('Set CWD = ' + str(m_join()) + '\n',
                         m_stdout.getvalue())


    def test_run_fails_on_bad_parse(self):
        with mock.patch('PullRequestLinuxDriverMerge.parseArgs',
                        side_effect=SystemExit(2)):
            self.assertFalse(PullRequestLinuxDriverMerge.run())

    def test_run_fails_on_bad_fetch(self):
        with mock.patch('subprocess.check_output',
                        side_effect=['origin /dev/null/target/Trilinos',
                                      'df324ae']), \
            mock.patch('subprocess.check_call',
                       side_effect=[None,
                                    CalledProcessError(-1, 'cmd'),
                                    CalledProcessError(-2, 'cmd'),
                                    CalledProcessError(-3, 'cmd')]), \
            mock.patch('PullRequestLinuxDriverMerge.parseArgs'), \
            mock.patch('sys.stdout', new_callable=StringIO), \
            mock.patch('os.path.join'), \
            mock.patch('os.chdir'):
            self.assertFalse(PullRequestLinuxDriverMerge.run())


    def test_run_fails_on_bad_remote_add(self):
        expected_string = '''Recieved subprocess.CalledProcessError - returned -1
  from command test_cmd
  output None
  stdout None
  stderr None\n'''
        if sys.version_info.major is not 3:
            expected_string = '''Recieved subprocess.CalledProcessError - returned -1
  from command test_cmd
  output None\n'''
        with mock.patch('subprocess.check_output',
                        side_effect=['origin /dev/null/target/Trilinos',
                                      'df324ae']), \
            mock.patch('subprocess.check_call',
                       side_effect=CalledProcessError(-1, 'test_cmd')), \
            mock.patch('PullRequestLinuxDriverMerge.parseArgs'), \
            mock.patch('sys.stdout', new_callable=StringIO) as m_stdout, \
            mock.patch('os.path.join'), \
            mock.patch('os.chdir'):
            self.assertFalse(PullRequestLinuxDriverMerge.run())
        self.assertTrue(m_stdout.getvalue().endswith(expected_string))


if __name__ == '__main__':
    unittest.main()  # pragma nocover
