#!/usr/bin/env python

"""
Python module/program that counts source code stats for a series of directories
for files of different types.
"""

# ToDo: I want to add the following statistics as well:
# 0) Count the number of semi-colons
# 1) Count the number of tokens in all of the source code lines
# 2) Count the number of total non-space chars in all of the source code lines
# 3) Count the number of total non-space chars in all of the documentation lines
# 4) Separate out the copyright header comment lines from the other comment
#    lines and print the line counts for these separately

import sys
import os
import re

reBeginCStyleCommentLine = re.compile(r"^\s*/\*.*")
reEndCStyleCommentLine = re.compile(r".*\*/.*")
reCppStyleCommentLine = re.compile(r"^\s*//.*")
reBlankLine = re.compile(r"^\s*$")

def countLines(fileName):
  #
  file = open(fileName)
  #
  numBlankLines = 0
  numCopyrightHeaderLines = 0
  numCommentLines = 0
  numCodeLines = 0
  numTotalLines = 0
  numSrcTokens = 0
  numSrcChars = 0
  numDocChars = 0
  #
  activeCommentBlock = 0
  #
  for line in file:
    line = line.rstrip()
    #print "\nline = \"" + line + "\"\n"
    numTotalLines += 1
    if not activeCommentBlock:
      #print "\nComment block is not active!\n";
      if re.match(reBeginCStyleCommentLine,line):
        #print "\nMatches the beginning of a C-style comment line!\n";
        numCommentLines += 1
        if not re.match(reEndCStyleCommentLine,line):
          #print "\nDoes not match the end of a C-stlye comment line!\n";
          activeCommentBlock = 1
        #else: print "\nMatches the end of a C-stlye comment line!\n";
      elif re.match(reCppStyleCommentLine,line):
        #print "\nMatches a C++-style comment line!\n";
        numCommentLines += 1
      elif re.match(reBlankLine,line):
        #print "\nMatches a blank line!\n";
        numBlankLines += 1
      else:
        #print "\nMust be a code line!\n";
        numCodeLines += 1
      # end if
    else:
      #print "\nComment block is active!\n";
      numCommentLines += 1
      if re.match(reEndCStyleCommentLine,line):
        #print "\nMatches end of C-sytle comment line!\n";
        activeCommentBlock = 0
    # end if
  # end for
  return (numBlankLines,numCommentLines,numCodeLines,numTotalLines)

if __name__ == '__main__':
  (numBlankLines,numCommentLines,numCodeLines,numTotalLines) \
    = countLines(sys.argv[1])
  print "File: ", os.path.basename(sys.argv[1])
  print "Number of code lines     =", numCodeLines
  print "Number of comment lines  =", numCommentLines
  print "Number of blank lines    =", numBlankLines
  print "numTotalLines (", numTotalLines, ")"

