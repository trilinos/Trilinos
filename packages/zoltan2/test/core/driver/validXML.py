#!/usr/bin/python
#
# Figure out if input file is valid XML
#
#  Type "validXML {xmlFileName}
#

import elementtree.ElementTree as ET
import sys

##
## Begin
##

if len(sys.argv) < 2:
  print "usage: ", sys.argv[0], " {xmlFileName}"

else:
  fname = sys.argv[1]
  print "Checking ",fname

  tree = ET.parse(fname)

  root = tree.getroot()

  if root.tag == "Tests":
    print "format is OK"
  else:
    print "Invalid, \"Tests\" is not defined in the XML file."
