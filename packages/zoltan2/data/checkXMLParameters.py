#!/usr/bin/python
#
# Element is available in python 1.5.2 and higher.
# It can be downloaded to use with older python.
#
# This script does some basic validity checks on parameters.xml.

import elementtree.ElementTree as ET
import sys

##
## Begin
##

tree = ET.parse("./parameters.xml")

root = tree.getroot()

if root.tag != "ParameterList":
  print "Error: Root tag is not ParameterList"
  sys.exit(1)

validators = []
for node in root:
  if node.tag == "Validators":
    validators = node
    break

if len(validators) == 0:
  print "Error: This is not a valid Zoltan2 parameter list."
  sys.exit(1)
  

# Create a dictionary of Validators

vals={}

for node in validators:
  id = node.get("validatorId")
  if id in vals:
    print "Error: Validator ID ",id," is repeated."
    sys.exit(1)

  vals[id] = node 

##
# Match up parameters to validators
##

idList = []
fail = False

for node in root:
  if node.tag != "Parameter":
    continue
  paramName = node.get("name")
  paramId = node.get("id")
  validatorId = node.get("validatorId")

  if paramId in idList:
    print "Error: Parameter id ",paramId," is reused for ",paramName
    fail = True
    break

  idList.append(paramId)

  if validatorId not in vals.keys():
    print "Error: Parameter ",paramName," has invalid validator ID ",validatorId
    fail = True
    break

if fail == True:
  print "Looks OK"
  
  
