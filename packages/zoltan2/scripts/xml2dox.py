#!/usr/bin/python
#
# Element is available in python 1.5.2 and higher.
# It can be downloaded to use with older python.
#
# This script runs in zoltan2/doc.  It reads an XML description
# of the Zoltan2 parameters, and writes a doxygen page with
# this information.
#
import elementtree.ElementTree as ET

tree = ET.parse("../data/parameters.xml")

root = tree.getroot()
