# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
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
# ************************************************************************
# @HEADER

##################################################
# Unit testing code for extract_rst_cmake_doc.py #
##################################################


from GeneralScriptSupport import *
from extract_rst_cmake_doc import *

import unittest
import filecmp


scriptBaseDir = getScriptBaseDir()


#
# Test getBaseEntityName()
#


class test_getBaseEntityName(unittest.TestCase):

  def test_noparenth(self):
    self.assertEqual(getBaseEntityName("SomeName"), "SomeName")

  def test_withparenth_closed(self):
    self.assertEqual(getBaseEntityName("SomeName()"), "SomeName")

  def test_withparenth_open(self):
    self.assertEqual(getBaseEntityName("SomeName( )"), "SomeName")


#
# Test data
#


# This extracts the RST comment block
simpleDocText1 = """
blah blah blah

something

#
# @MACRO: SOME_MACRO_NAME1()
#
# Good documenation
#
MACRO(SOME_MACRO_NAME1 ...)
some other stuff 
...
"""

simpleDocText1_rstDocBlocks_expected = {
  "SOME_MACRO_NAME1()" : {
    "type" : "MACRO",
    "body" : "\nGood documenation\n\n"
    }
  }


# This extracts the same RST comment block as simpleDocText1
simpleDocWithSpacesText1 = """
blah blah blah

something

#
# @MACRO :   SOME_MACRO_NAME1()
#
# Good documenation
#
MACRO(  SOME_MACRO_NAME1 ...)
some other stuff 
...
"""


# This extracts the same RST comment block as simpleDocText1
simpleDocNoArgsText1 = """
blah blah blah

something

#
# @MACRO: SOME_MACRO_NAME1()
#
# Good documenation
#
MACRO(SOME_MACRO_NAME1)
some other stuff 
...
"""

# This extracts the RST comment block
simpleDocText2 = """
blah blah blah

something

#
# @FUNCTION: SOME_FUNC_NAME2()
# Better documenation
#
# Usage::
#
#   SOME_FUNC_NAME2(blah
#     goat
#     )
#
FUNCTION(SOME_FUNC_NAME2 ...)
some other stuff 
...
"""

simpleDocText2_rstDocBlocks_expected = {
  "SOME_FUNC_NAME2()" : {
    "type" : "FUNCTION",
    "body" : "Better documenation\n\nUsage::\n\n  SOME_FUNC_NAME2(blah\n    goat\n    )\n\n"
    }
  }


simpleDocText1_and_2 = simpleDocText1 + "\n\n\n" + simpleDocText2


simpleDocText1_and_2_rstDocBlocks_expected = {}
simpleDocText1_and_2_rstDocBlocks_expected.update(simpleDocText1_rstDocBlocks_expected)
simpleDocText1_and_2_rstDocBlocks_expected.update(simpleDocText2_rstDocBlocks_expected)


# This results in an error where the comment block is not extracted
funcMissingColon = """
#
# @FUNCTION SOME_MACRO_NAME1()
#
# Blah blah
#
FUNCTION(SOME_MACRO_NAME1 ...)
"""


# This results in an error where the comment block is not extracted
funcTerminateOnMacroText = """
#
# @FUNCTION: SOME_MACRO_NAME1()
#
# Blah blah
#
MACRO(SOME_MACRO_NAME1 ...)
"""


# This results in an error where the comment block is not extracted
funcNameMistmatchText = """
#
# @FUNCTION: SOME_FUNC_NAME1()
# Blah blah
#
FUNCTION(SOME_FUNC_NAME2 ...)
"""


# This results in an error where the comment block is not extracted
badVerticalSpaceInCommentBlockText = """
#
# @FUNCTION: SOME_FUNC_NAME1()

# Blah blah
#
FUNCTION(SOME_FUNC_NAME1 ...)
"""


# This results in an error where the comment block is not extracted
msisingHorizontalSpaceInCommentBlockText = """
#
# @FUNCTION: SOME_FUNC_NAME1()
#
#Blah blah
#
FUNCTION(SOME_FUNC_NAME1 ...)
"""


# This results in the comment block not being extracted
tooManySpacesBeforeBlockIndentifierText = """
#
#  @FUNCTION: SOME_FUNC_NAME1()
#
# Blah blah
#
FUNCTION(SOME_FUNC_NAME1 ...)
"""


#
# Test extractRstDocBlocksFromText()
#


rstBlockTypes = ["MACRO", "FUNCTION"]


class test_extractRstDocBlocksFromText(unittest.TestCase):


  def test_extract_1_block_simple_1(self):
    rstDocBlocks = extractRstDocBlocksFromText(simpleDocText1, rstBlockTypes, "")
    self.assertEqual(rstDocBlocks, simpleDocText1_rstDocBlocks_expected)


  def test_extract_1_block_simple_2(self):
    rstDocBlocks = extractRstDocBlocksFromText(simpleDocText2, rstBlockTypes, "")
    self.assertEqual(rstDocBlocks, simpleDocText2_rstDocBlocks_expected)


  def test_extract_2_blocks_simle_1_2(self):
    rstDocBlocks = extractRstDocBlocksFromText(simpleDocText1_and_2, rstBlockTypes, "")
    self.assertEqual(rstDocBlocks, simpleDocText1_and_2_rstDocBlocks_expected)


  def test_extract_1_block_simple_with_spaces_1(self):
    rstDocBlocks = extractRstDocBlocksFromText(simpleDocWithSpacesText1, rstBlockTypes, "")
    self.assertEqual(rstDocBlocks, simpleDocText1_rstDocBlocks_expected)


  def test_extract_1_block_simple_no_args_1(self):
    rstDocBlocks = extractRstDocBlocksFromText(simpleDocNoArgsText1, rstBlockTypes, "")
    self.assertEqual(rstDocBlocks, simpleDocText1_rstDocBlocks_expected)


  def test_func_mussing_colon(self):
    exceptMessage = "NO EXCEPTION WAS THOWN"
    try:
      rstDocBlocks = extractRstDocBlocksFromText(funcMissingColon, rstBlockTypes,
        "someFile1.cmake")
    except Exception, e:
      exceptMessage = e.args[0]
    exceptMessage_expected = "someFile1.cmake:3: error: '# @FUNCTION SOME_MACRO_NAME1()' is missing the colon ':' separator!"
    self.assertEqual(exceptMessage, exceptMessage_expected)


  def test_func_terminate_on_macro(self):
    exceptMessage = "NO EXCEPTION WAS THOWN"
    try:
      rstDocBlocks = extractRstDocBlocksFromText(funcTerminateOnMacroText, rstBlockTypes,
        "someFile1.cmake")
    except Exception, e:
      exceptMessage = e.args[0]
    exceptMessage_expected = "someFile1.cmake:7: error: expecting 'FUNCTION' but got type 'MACRO''"
    self.assertEqual(exceptMessage, exceptMessage_expected)


  def test_func_name_mismatch(self):
    exceptMessage = "NO EXCEPTION WAS THOWN"
    try:
      rstDocBlocks = extractRstDocBlocksFromText(funcNameMistmatchText, rstBlockTypes,
        "someFile2.cmake")
    except Exception, e:
      exceptMessage = e.args[0]
    exceptMessage_expected = "someFile2.cmake:6: error: expecting FUNCTION SOME_FUNC_NAME1() but got wrong FUNCTION name 'SOME_FUNC_NAME2'"
    self.assertEqual(exceptMessage, exceptMessage_expected)


  def test_bad_vertical_space_in_comment_block(self):
    exceptMessage = "NO EXCEPTION WAS THOWN"
    try:
      rstDocBlocks = extractRstDocBlocksFromText(badVerticalSpaceInCommentBlockText, rstBlockTypes,
        "someFile3.cmake")
    except Exception, e:
      exceptMessage = e.args[0]
    exceptMessage_expected = "someFile3.cmake:4: error: expecting FUNCTION(SOME_FUNC_NAME1 ...) on this line.  RST comment block must terminate in the stated entity!"
    self.assertEqual(exceptMessage, exceptMessage_expected)


  def test_missing_horizontal_space_in_comment_block(self):
    exceptMessage = "NO EXCEPTION WAS THOWN"
    try:
      rstDocBlocks = extractRstDocBlocksFromText(msisingHorizontalSpaceInCommentBlockText,
         rstBlockTypes, "someFile4.cmake")
    except Exception, e:
      exceptMessage = e.args[0]
    exceptMessage_expected = "someFile4.cmake:5: error: Comment blocks must have at least one space after '#' in line = '#Blah blah'"
    self.assertEqual(exceptMessage, exceptMessage_expected)


  def test_too_many_spaces_before_block_identifier(self):
    rstDocBlocks = extractRstDocBlocksFromText(tooManySpacesBeforeBlockIndentifierText,
    rstBlockTypes, "")
    self.assertEqual(rstDocBlocks, {})


#
# Test replaceWithRstDocBlocksInText()
#



textNoReplacement = """

something 1

other
"""

textToReplace1 = """

something 1

other

@FUNCTION: SOME_FUNC_NAME2() -

something 2

@MACRO: SOME_MACRO_NAME1() +

something else

"""


# Gives the same substitutions as textToReplace1
textToReplaceWithSpaces1 = """

something 1

other

@FUNCTION :  SOME_FUNC_NAME2()   -

something 2

@MACRO  :   SOME_MACRO_NAME1()  +

something else

"""


# Gives the same substitutions as textToReplace1
textToReplaceNoSpaceAfterColon1 = """

something 1

other

@FUNCTION:SOME_FUNC_NAME2() -

something 2

@MACRO:SOME_MACRO_NAME1() +

something else

"""


# NOTE that this adds what looks like an extra line after each replacement!
replacedText1_expected = """

something 1

other

SOME_FUNC_NAME2()
-----------------
Better documenation

Usage::

  SOME_FUNC_NAME2(blah
    goat
    )


something 2

SOME_MACRO_NAME1()
++++++++++++++++++

Good documenation


something else

"""


textMissingColon = """
a
@FUNCTION SOME_FUNC_NAME2() -
b
"""


textMissingSecChar = """
a
@FUNCTION: SOME_FUNC_NAME2()
b
"""


textSepCharTooLong = """
a
@FUNCTION: SOME_FUNC_NAME2() --
b
"""


textWrongBlockType = """
a
@MACRO: SOME_FUNC_NAME2() -
b
"""

textMisspelledOrMissingBlockName = """
a
@MACRO: MissingBlockName() -
b
"""


def lineByLineCompareAssert(testObj, text, textExpected):
  textArray = text.split("\n")
  textExpectedArray = textExpected.split("\n")
  
  for i in range(0,min(len(textArray), len(textExpectedArray))):
    testObj.assertEqual(textArray[i], textExpectedArray[i],
      "AssertionError: textArray["+str(i)+"] = "+textArray[i]+" != "+\
      "textExpectedArray["+str(i)+"] = "+textExpectedArray[i]+"\n\n"+\
      "text = {\n"+text+"}\n\n"+\
      "textExpected = {\n"+textExpected+"}\n")
  testObj.assertEqual(len(textArray), len(textExpectedArray))


class test_replaceWithRstDocBlocksInText(unittest.TestCase):
 

  def test_noreplace(self):
    replacedText = replaceWithRstDocBlocksInText(textNoReplacement,
      rstBlockTypes, simpleDocText1_and_2_rstDocBlocks_expected, "")
    self.assertEqual(replacedText, textNoReplacement)
 

  def test_replace_1(self):
    replacedText = replaceWithRstDocBlocksInText(textToReplace1,
      rstBlockTypes, simpleDocText1_and_2_rstDocBlocks_expected, "")
    lineByLineCompareAssert(self, replacedText, replacedText1_expected)
    self.assertEqual(replacedText, replacedText1_expected)
 

  def test_replace_with_spaces_1(self):
    replacedText = replaceWithRstDocBlocksInText(textToReplaceWithSpaces1,
      rstBlockTypes, simpleDocText1_and_2_rstDocBlocks_expected, "")
    lineByLineCompareAssert(self, replacedText, replacedText1_expected)
    self.assertEqual(replacedText, replacedText1_expected)
 

  def test_replace_no_space_after_colon_1(self):
    replacedText = replaceWithRstDocBlocksInText(textToReplaceNoSpaceAfterColon1,
      rstBlockTypes, simpleDocText1_and_2_rstDocBlocks_expected, "")
    lineByLineCompareAssert(self, replacedText, replacedText1_expected)
    self.assertEqual(replacedText, replacedText1_expected)
 

  def test_missing_colon(self):
    replacedText = replaceWithRstDocBlocksInText(textMissingColon,
      rstBlockTypes, simpleDocText1_and_2_rstDocBlocks_expected, "")
    self.assertEqual(replacedText, textMissingColon)
 

  def test_missing_sec_char(self):
    exceptMessage = "NO EXCEPTION WAS THOWN"
    try:
      replacedText = replaceWithRstDocBlocksInText(textMissingSecChar,
       rstBlockTypes, simpleDocText1_and_2_rstDocBlocks_expected, "someFile1.cmake")
    except Exception, e:
      exceptMessage = e.args[0]
    exceptMessage_expected = "someFile1.cmake:3: Error, line does not match format '@FUNCTION: <blockName> <sepChar>' for line = '@FUNCTION: SOME_FUNC_NAME2()'"
    self.assertEqual(exceptMessage, exceptMessage_expected)


  def test_sec_char_too_long(self):
    exceptMessage = "NO EXCEPTION WAS THOWN"
    try:
      replacedText = replaceWithRstDocBlocksInText(textSepCharTooLong,
       rstBlockTypes, simpleDocText1_and_2_rstDocBlocks_expected, "someFile2.cmake")
    except Exception, e:
      exceptMessage = e.args[0]
    exceptMessage_expected = "someFile2.cmake:3: Error, separation char must be on char in line = '@FUNCTION: SOME_FUNC_NAME2() --'"
    self.assertEqual(exceptMessage, exceptMessage_expected)


  def test_wrong_block_type(self):
    exceptMessage = "NO EXCEPTION WAS THOWN"
    try:
      replacedText = replaceWithRstDocBlocksInText(textWrongBlockType,
       rstBlockTypes, simpleDocText1_and_2_rstDocBlocks_expected, "someFile3.cmake")
    except Exception, e:
      exceptMessage = e.args[0]
    exceptMessage_expected = "someFile3.cmake:3: Error, looked up block type '@FUNCTION' does not match block type for line = '@MACRO: SOME_FUNC_NAME2() -'"
    self.assertEqual(exceptMessage, exceptMessage_expected)


  def test_wrong_block_type(self):
    exceptMessage = "NO EXCEPTION WAS THOWN"
    try:
      replacedText = replaceWithRstDocBlocksInText(textMisspelledOrMissingBlockName,
       rstBlockTypes, simpleDocText1_and_2_rstDocBlocks_expected, "someFile4.cmake")
    except Exception, e:
      exceptMessage = e.args[0]
    exceptMessage_expected = "someFile4.cmake:3: Error, block name 'MissingBlockName()' does not exist in rstDocBlocks for line = '@MACRO: MissingBlockName() -'"
    self.assertEqual(exceptMessage, exceptMessage_expected)


#
# Mock comamndline options
#



class MockOptions:
  def __init__(self):
    self.extractFrom = ""
    self.fileExtensions = ".cmake"
    self.rstFilePairs = ""


#
# Test getExtractFilesList()
#


class test_getExtractFilesList(unittest.TestCase):

  def test_1_file(self):
    opts = MockOptions()
    opts.extractFrom = "file1.cmake"
    self.assertEqual(getExtractFilesList(opts), ["file1.cmake"] )

  def test_2_files(self):
    opts = MockOptions()
    opts.extractFrom = "file2.cmake,file1.cmake"
    self.assertEqual(getExtractFilesList(opts), ["file2.cmake","file1.cmake"] )

  def test_glob_data_files(self):
    baseDir = scriptBaseDir+"/UnitTests/extract_rst_cmake_doc"
    opts = MockOptions()
    opts.extractFrom = baseDir+"/"
    filesList = getExtractFilesList(opts)
    filesList.sort()
    filesList_expected = [
      baseDir+"/simpleDocText1.cmake",
      baseDir+"/simpleDocText2.cmake"
      ]
    self.assertEqual(filesList, filesList_expected )

  def test_glob_data_files_2_dirs(self):
    baseDir = scriptBaseDir+"/UnitTests/extract_rst_cmake_doc"
    opts = MockOptions()
    opts.extractFrom = baseDir+"/,"+baseDir+"/"
    filesList = getExtractFilesList(opts)
    filesList.sort()
    filesList_expected = [
      baseDir+"/simpleDocText1.cmake",
      baseDir+"/simpleDocText1.cmake",
      baseDir+"/simpleDocText2.cmake",
      baseDir+"/simpleDocText2.cmake"
      ]
    self.assertEqual(filesList, filesList_expected )


#
# Test getExtractFilesList()
#


class test_getRstFilesList(unittest.TestCase):

  def test_empty(self):
    opts = MockOptions()
    opts.rstFilePairs = ""
    self.assertEqual(getRstFilesList(opts), [])

  def test_1_pair(self):
    opts = MockOptions()
    opts.rstFilePairs = "somebase/someTemplate.rst:otherbase/some.rst"
    self.assertEqual(getRstFilesList(opts),
      [ ["somebase/someTemplate.rst" ,"otherbase/some.rst"] ] )

  def test_2_pairs(self):
    opts = MockOptions()
    opts.rstFilePairs = "atmp.rst:a.rst,bt.rst.rst:b.rst"
    self.assertEqual(getRstFilesList(opts),
      [ ["atmp.rst", "a.rst" ], [ "bt.rst.rst", "b.rst"] ] )


#
# Test extractRstDocBlocksFromFileList()
#


class test_extractRstDocBlocksFromFileList(unittest.TestCase):

  def test_extract_1_block_simple_1(self):
    fileList = ["test_extract_1_block_simple_1.cmake"]
    open(fileList[0], 'w').write(simpleDocText1)
    rstDocBlocks = extractRstDocBlocksFromFileList(fileList, rstBlockTypes)
    self.assertEqual(rstDocBlocks, simpleDocText1_rstDocBlocks_expected)



#
# Test replaceWithRstDocBlocksInTemplateFile()
#


class test_replaceWithRstDocBlocksInTemplateFile(unittest.TestCase):

  def test_replace_1_block_1_file(self):
    baseDir = scriptBaseDir+"/UnitTests/extract_rst_cmake_doc"
    templateFileName = baseDir+"/simpleTemplate1.rst"
    fileName = "test_replace_1_block_1_file.rst"
    rstFileList = [ [templateFileName , fileName] ] 
    if os.path.exists(fileName): os.remove(fileName)
    replaceWithRstDocBlocksInTemplateFileList(rstFileList, rstBlockTypes,
      simpleDocText1_and_2_rstDocBlocks_expected)
    self.assertTrue(filecmp.cmp(fileName, baseDir+"/"+fileName+".gold"))


#
# Run the run tests!
#


if __name__ == '__main__':
  unittest.main()
