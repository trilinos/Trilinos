#!/usr/bin/env python
#ML timelines for matrix matrix multiply in RAP
LABELS    = ['i&x', 'right serialCore', 'right fc', 'left serialCore', 'left fc']
TIMELINES    = ['RAP right: pre-multiply', 'RAP right: multiply time', 'RAP right: post-multiply', 'RAP left:  multiply', 'RAP left:  post-multiply']
def PARSEFUNC(fileName,stringToFind):
  return "grep -i \"" + stringToFind + "\" " + fileName + " | cut -f2 -d'='"
