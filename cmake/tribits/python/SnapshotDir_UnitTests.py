# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER

########################################
# Unit testing code for SnapshotDir.py #
########################################


from SnapshotDir import *
import unittest
import re


#
# Unit test support code
#


scriptsDir = getScriptBaseDir()


class WriteToString:
  def __init__(self):
    self.str = ""
  def write(self, s):
    self.str += s;
  def flush(self):
    None
  def getStr(self):
    return self.str


def getDummyDefaultOptions():
  dummyDefaultOptions = DefaultOptions()
  dummyDefaultOptions.setDefaultOrigDir("dummy/orig/dir/")
  dummyDefaultOptions.setDefaultDestDir("dummy/dest/dir/")
  return dummyDefaultOptions


def runSnapshotDirTestCase(testObject, cmndLineArgsList, cmndInterceptList,
  passRegexExpressionsList, defaultOptions=None \
  ):

  # Set up default options
  if not defaultOptions:
    defaultOptions = getDummyDefaultOptions()

  # Set up the command intercepts
  g_sysCmndInterceptor.readCommandsFromStr("".join(cmndInterceptList))
  g_sysCmndInterceptor.setAllowExtraCmnds(False)

  # Run the command, intercept the output, and test it
  sout = WriteToString()
  try:
    rtn = snapshotDirMainDriver(cmndLineArgsList, defaultOptions, sout)
    ostr = sout.getStr()
    #print ostr
    for passRegexExpr in passRegexExpressionsList:
      testObject.assert_(re.search(passRegexExpr, ostr))
  except Exception, e:
    print "\n\nGenerated output:\n\n" + sout.getStr() + "\n\n"
    printStackTrace()
    raise e



#
# Standard commands used in snapshot-dir.py
#

g_gitDiffHead = "IT: git diff --name-status HEAD -- \.; 0;''\n"

g_gitRemote = "IT: git remote -v; 0; 'origin\tsome-url-location (fetch)'\n"

g_gitLog = "IT: git log  --pretty=.*; 0; 'one commit msg'\n"

g_rsync = "IT: rsync -av --delete --exclude=.* dummy/orig/dir/ dummy/dest/dir/; 0; 'sync passed'\n"

g_gitAdd = "IT: git add \.; 0; 'added some files'\n"

g_gitCommit = "IT: git commit .+; 0; 'did a commit'\n"


#
# Unit test SnapshotDir
#

class test_snapshot_dir(unittest.TestCase):


  def test_show_defaults(self):
    runSnapshotDirTestCase(
      self,
      ["--show-defaults"],
      [],
      [
        "Script: snapshot-dir\.py",
        "--orig-dir='dummy/orig/dir/'",
        "--dest-dir='dummy/dest/dir/'"
        ]
      )


  def test_override_orig_dest_dirs(self):
    runSnapshotDirTestCase(
      self,
      ["--show-defaults", "--orig-dir=new/orig/dir", "--dest-dir=new/dest/dir"],
      [],
      [
        "Script: snapshot-dir\.py",
        "--orig-dir='new/orig/dir'",
        "--dest-dir='new/dest/dir'"
        ]
     )


  def test_full_snapshot(self):
    runSnapshotDirTestCase(
      self,
      [""],
      [
        g_gitDiffHead,
        g_gitDiffHead,
        g_gitRemote,
        g_gitLog,
        g_rsync,
        g_gitAdd,
        g_gitCommit,
        ],
      [
        "Script: snapshot-dir\.py",
        "--orig-dir='dummy/orig/dir/'",
        "--dest-dir='dummy/dest/dir/'"
        ]
     )

  # ToDo: Test assert failure of clean origDir ...

  # ToDo: Test skipping test of clean origDir ...

  # ToDo: Test assert failure of clean destDir ...

  # ToDo: Test skipping test of clean destDir ...

  # ToDo: Test failing to acquire origin remote URL ...

  # ToDo: Test failure to acquire origin commit ...

  # ToDo: Test skipping getting origin info ...

  # ToDo: Test failing rysnc ...

  # ToDo: Test failing to create commit in dest repo ...

  # ToDo: Test skipping creation of comit in dest repo ...


if __name__ == '__main__':
  unittest.main()
