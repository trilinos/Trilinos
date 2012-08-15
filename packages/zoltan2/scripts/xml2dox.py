#!/usr/bin/python
#
# Element is available in python 1.5.2 and higher.
# It can be downloaded to use with older python.
#
# This script reads an XML description
# of the Zoltan2 parameters, and writes a doxygen page with
# this information.
#
#
#  The Doxygen page we want to create is formatted this way:
#
# /*! \page z2_parameters Zoltan2 Parameters
# 
# <DL>
# 
# <DT> parameter_1 \anchor parameter_1 </DT>
# <DD> explanation of what it means </DD>
# 
# <DT> parameter_2 \anchor parameter_2 </DT>
# <DD> explanation of what it means </DD>
# 
# </DL>
# 
# */

import elementtree.ElementTree as ET

tree = ET.parse("../data/parameters.xml")

root = tree.getroot()
