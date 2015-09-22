#!/usr/bin/python
#
# Element is available in python 1.5.2 and higher.
# It can be downloaded to use with older python.
#
# This script reads an XML description of the Zoltan2 parameters, 
# and writes a doxygen page with this information.

import elementtree.ElementTree as ET

outfile = open("parameters.dox", "w")

def parameterDocumentationHeader():
  outfile.write("/*! \\page z2_parameters Zoltan2 Parameters\n\n")
  outfile.write("This page lists each Zoltan2 parameter and how to use it.  The validators\n")
  outfile.write("are classes of interest to Zoltan2 developers.  They are used to evaluate\n")
  outfile.write("the validity of parameter values at run time.\n\n")

def parameterDocumentationFooter():
  outfile.write("*/\n")

def parameterEnhancedNumber(pname, pinfo, pval):
  desc = pinfo.get("docString")
  min = pval.get("min")
  max = pval.get("max")
  step = pval.get("step")

  outfile.write("- \\b "+pname+" \\anchor "+pname+"\n")
  outfile.write("  - Description: "+desc+"\n")
  outfile.write("  - Valid values:\n")
  outfile.write("     - minimum is "+min+"\n")
  outfile.write("     - maximum is "+max+"\n")
  outfile.write("     - step is "+step+"\n")
  outfile.write("  - Validator type: Teuchos::EnhancedNumberValidator\n")
  outfile.write("\n")

def parameterIntegerRangeList(pname, pinfo, pval):
  desc = pinfo.get("docString")
  unsorted = pval.get("unsorted")
  min="unset"
  max="unset"
  if "min" in pinfo.keys():
    min = pval.get("min")
  if "max" in pinfo.keys():
    max = pval.get("max")

  outfile.write("- \\b "+pname+" \\anchor "+pname+"\n")
  outfile.write("  - Description: "+desc+"\n")
  outfile.write("  - Valid values: a comma-separated list of\n")
  outfile.write("     - numbers\n")
  outfile.write("     - number ranges separated by a dash (\"1-5\")\n")
  outfile.write("     - the word \"all\" to indicate all possible values\n")
  if min != "unset":
    outfile.write("     - minimum is: "+min+"\n")
  if max != "unset":
    outfile.write("     - maximum is: "+max+"\n")

  outfile.write("  - Examples: \"1,2,7\", \"all\", \"5,1-15,80-82,99\"\n")
  outfile.write("  - Validator type: Zoltan2::IntegerRangeListValidator\n")
  if unsorted == "true":
    outfile.write(      "(list is not changed during processing)\n")
  else:
    outfile.write(      "(list will be sorted, and duplicates removed, during processing)\n")

  outfile.write("\n")

def parameterFileName(pname, pinfo, pval):
  desc = pinfo.get("docString")
  outfile.write("- \\b "+pname+" \\anchor "+pname+"\n")
  outfile.write("  - Description: "+desc+"\n")
  outfile.write("  - Validator type: Teuchos::FileNameValidator\n")
  outfile.write("\n")

def parameterAnyNumber(pname, pinfo, pval):
  desc = pinfo.get("docString")
  validTypes = []
  if pval.get("allowDouble") == "true":
    validTypes.append("double")
  if pval.get("allowInt") == "true":
    validTypes.append("int")
  if pval.get("allowString") == "true":
    validTypes.append("string")

  outfile.write("- \\b "+pname+" \\anchor "+pname+"\n")
  outfile.write("  - Description: "+desc+"\n")
  outfile.write("  - Valid values are any values of type:\n")
  for vtype in validTypes:
    outfile.write("      \\li "+vtype+"\n")
  outfile.write("  - Validator type: Teuchos::AnyNumberParameterEntryValidator\n")
  outfile.write("\n")

def parameterStringToIntegral(pname, pinfo, pval):
  desc = pinfo.get("docString")
  outfile.write("- \\b "+pname+" \\anchor "+pname+"\n")
  outfile.write("  - Description: "+desc+"\n")
  outfile.write("  - Valid values:\n")
  for node in pval:
    if node.tag == "String":
      svalue = node.get("stringValue")
      sdoc = "unset"
      if "stringDoc" in node.keys():
        sdoc = node.get("stringDoc")
      if sdoc == "unset":
        outfile.write("    \\li \\e "+svalue+"\n")
      else:
        outfile.write("    \\li \\e "+svalue+" "+sdoc+"\n")
  outfile.write("  - Validator type: Teuchos::StringToIntegralParameterEntryValidator\n")
  outfile.write("\n")

def parameterFileName(pname, pinfo, pval):
  desc = pinfo.get("docString")
  mustExist = pinfo.get("fileMustExist")
  outfile.write("- \\b "+pname+" \\anchor "+pname+"\n")
  outfile.write("  - Description: "+desc+"\n")
  if mustExist == "true":
    outfile.write("  File must exist.\n")
  else:
    outfile.write("  File need not already exist.\n")
  outfile.write("  - Validator type: Teuchos::FileNameValidator\n")
  outfile.write("\n")

def parameterString(pname, pinfo, pval):
  desc = pinfo.get("docString")
  outfile.write("- \\b "+pname+" \\anchor "+pname+"\n")
  outfile.write("  - Description: "+desc+"\n")
  outfile.write("  - Valid values:\n")
  for node in pval:
    if node.tag == "String":
      outfile.write("    \\li \\e "+node.get("value")+"\n")
  outfile.write("  - Validator type: Teuchos::StringValidator\n")
  outfile.write("\n")

def writeInfo(param):
  pname = param[0]
  pinfo = param[1]
  pval = param[2]

  pvalidatorType = pval.get("type")
  
  if pvalidatorType == "anynumberValidator":
    parameterAnyNumber(pname, pinfo, pval)

  elif pvalidatorType == "FilenameValidator":
    parameterFileName(pname, pinfo, pval)

  elif pvalidatorType == "StringValidator":
    parameterString(pname, pinfo, pval)

  elif "StringIntegralValidator" in pvalidatorType:
    parameterStringToIntegral(pname, pinfo, pval)

  elif "IntegerRangeListValidator" in pvalidatorType:
    parameterIntegerRangeList(pname, pinfo, pval)

  elif "EnhancedNumberValidator" in pvalidatorType:
    parameterEnhancedNumber(pname, pinfo, pval)

  else:
    print "Error 4: This is not a valid Zoltan2 parameter list."
    exit

##
## Begin
##

tree = ET.parse("../data/parameters.xml")

root = tree.getroot()

if root.tag != "ParameterList":
  print "Error 1: This is not a valid Zoltan2 parameter list."
  exit

validators = []
for node in root:
  if node.tag == "Validators":
    validators = node
    break

if len(validators) == 0:
  print "Error 1: This is not a valid Zoltan2 parameter list."
  exit

# Create a dictionary of Validators

vals={}

for node in validators:
  id = node.get("validatorId")
  vals[id] = node 

##
# Create list of a 3-tuples for each parameter, including
# the parameter name, its data, and its validator.
##

parameterInfo = []

for node in root:
  if node.tag != "Parameter":
    continue
  id = node.get("validatorId")
  if id not in vals.keys():
    print "Error 3: This is not a valid Zoltan2 parameter list."
    exit

  paramName = node.get("name")

  parameterInfo.append((paramName, node, vals[id]))
  
##
# Write the doxygen documentation for these parameters
##

parameterDocumentationHeader()

for info in sorted(set(parameterInfo)):
  print "Parameter: ",info[0]
  writeInfo(info)

parameterDocumentationFooter()

outfile.close()
  
  
