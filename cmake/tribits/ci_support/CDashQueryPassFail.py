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

try:
  # Python 2
  from urllib2 import urlopen
except ImportError:
  # Python 3
  from urllib.request import urlopen

import json
import datetime
import pprint

from FindGeneralScriptSupport import *

pp = pprint.PrettyPrinter()


# Validate a date format
def validateYYYYMMDD(dateText):
  try:
    return datetime.datetime.strptime(dateText, '%Y-%m-%d')
  except ValueError:
    raise ValueError("Incorrect data format for '"+dateText+"', should be YYYY-MM-DD")


# Construct the full query URL given the pieces
def getCDashIndexQueryUrl(cdashUrl, projectName, date, filterFields):
  return cdashUrl+"/api/v1/index.php?project="+projectName+"&date="+date \
    + "&"+filterFields


# Given a CDash query URL, return the full Python CDash data-structure
def extractCDashApiQueryData(cdashApiQueryUrl):
  response = urlopen(cdashApiQueryUrl)
  return json.load(response)


# Collect CDash index.php build summary fields
def collectCDashIndexBuildSummaryFields(fullCDashIndexBuild):
  summaryBuild = {
    u('buildname') : fullCDashIndexBuild.get('buildname', 'missing_build_name'),
    u('update') : \
      fullCDashIndexBuild.get('update', {'errors':9999,'this_field_was_missing':1}),
    u('configure') : \
      fullCDashIndexBuild.get('configure', {'error':9999,'this_field_was_missing':1}),
    u('compilation') : \
      fullCDashIndexBuild.get('compilation', {'error':9999,'this_field_was_missing':1}),
    u('test') : \
     fullCDashIndexBuild.get('test', {'fail':9999, 'notrun':9999,'this_field_was_missing':1} ),
    }
  return summaryBuild


# Given the full Python CDash API builds data-structure returned from the
# CDash index.php page and query, return an reduced data-structure to be used
# for pass/fail examination.
#
# This function takes in the data-structre directly returned from:
#
#   <cdash-url>/api/v1/index.php?project=<project>&date=<YYYY-MM-DD>&<filter-fields>
#
# The input full CDash API collapsed builds data-structure that has the
# following structure and fields of interest:
#
#  fullCDashIndexBuilds =
#  {
#    'all_buildgroups': [ {'id':1,'name:"Nightly"}, ...],
#    'buildgroups': [
#      {
#        'builds":[
#          {
#            'buildname':"???",
#            'update': {'errors':???, ...},
#            'configure':{'error': ???, ...},
#            'compilation':{'error': ???, ...},
#            'test': {'fail':???, 'notrun':???, 'pass':???, ...},
#            ...
#            },
#            ...
#          ]
#        },
#        ...
#      ...
#      ]
#      },
#      ...
#    }
#
# This function gets the data from *all* of the collapsed builds and returns
# the reduced data-structure:
#
#   [
#     {
#       'buildname':"???",
#       'update': {'errors':???, ...},
#       'configure':{'error': ???, ...},
#       'compilation':{'error': ???, ...},
#       'test': {'fail':???, 'notrun':???, 'pass':???, ...},
#       ...
#       },
#       ...
#       }
#
# This collects *all* of the builds from all of the build groups, not just the
# 'Nighlty' build group.  Therefore, if you want to only consider on set of
# build groups, you need to add that to the CDash query URL
# (e.g. group='Nighlty').
#
def getCDashIndexBuildsSummary(fullCDashIndexBuilds):
  summaryCDashIndexBuilds = []
  for buildgroup in fullCDashIndexBuilds["buildgroups"]:
    for build in buildgroup["builds"]:
      summaryBuild = collectCDashIndexBuildSummaryFields(build)
      summaryCDashIndexBuilds.append(summaryBuild)
  return summaryCDashIndexBuilds
  

# Return if a CDash Index build passes
def cdashIndexBuildPasses(cdashIndexBuild):
  if cdashIndexBuild['update']['errors'] > 0:
    return False
  if cdashIndexBuild['configure']['error'] > 0:
    return False
  if cdashIndexBuild['compilation']['error'] > 0:
    return False
  if (cdashIndexBuild['test']['fail'] + cdashIndexBuild['test']['notrun'])  > 0:
    return False
  return True
  

# Return if a list of CDash builds pass or fail and return error string if
# they fail.
def cdashIndexBuildsPass(summaryCDashIndexBuilds):
  buildsPass = True
  buildFailedMsg = ""
  for build in summaryCDashIndexBuilds:
    if not cdashIndexBuildPasses(build):
      buildsPass = False
      buildFailedMsg = "Error, the build " + sorted_dict_str(build) + " failed!"
      break
  return (buildsPass, buildFailedMsg)


# Extract the set of build names from a list of build names
def getCDashIndexBuildNames(summaryCDashIndexBuilds):
  buildNames = []
  for build in summaryCDashIndexBuilds:
    buildNames.append(build['buildname'])
  return buildNames


# Return if all of the expected builds exist and an error message if they
# don't.
def doAllExpectedBuildsExist(buildNames, expectedBuildNames):
  allExpectedBuildsExist = True
  errMsg = ""
  for expectedBuildName in expectedBuildNames:
    if findInSequence(buildNames, expectedBuildName) == -1:
      allExpectedBuildsExist = False
      errMsg = "Error, the expected build '"+expectedBuildName+"'" \
        +" does not exist in the list of builds "+str(buildNames) 
      break
  return (allExpectedBuildsExist, errMsg)    


# Return if a list of summary CDash index.php builds pass and has all of the
# expected builds.
def cdashIndexBuildsPassAndExpectedExist(summaryCDashIndexBuilds, 
  expectedBuildNames \
  ):
  cdashIndexBuildsPassAndExpectedExist_pass = True
  errMsg = ""
  # Check that all of the builds pass!
  if cdashIndexBuildsPassAndExpectedExist_pass:
    (buildsPass, buildFailedMsg) = cdashIndexBuildsPass(summaryCDashIndexBuilds)
    if not buildsPass:
      cdashIndexBuildsPassAndExpectedExist_pass = False
      errMsg = buildFailedMsg
  # Check that all of the expected builds are listed
  if cdashIndexBuildsPassAndExpectedExist_pass:
    buildNames = getCDashIndexBuildNames(summaryCDashIndexBuilds)
    (allExpectedBuildsExist, errMsg) = \
      doAllExpectedBuildsExist(buildNames, expectedBuildNames)
    if not allExpectedBuildsExist:
      cdashIndexBuildsPassAndExpectedExist_pass = False
      errMsg = errMsg
  return (cdashIndexBuildsPassAndExpectedExist_pass, errMsg)


# Determine if CDash index.php query builds all pass and has all expected
# builds.
def queryCDashAndDeterminePassFail(cdashUrl, projectName, date, filterFields,
  expectedBuildNames, printCDashUrl=True,
  extractCDashApiQueryData_in=extractCDashApiQueryData \
  ):
  # Get the query data
  cdashQueryUrl = getCDashIndexQueryUrl(cdashUrl, projectName, date, filterFields)
  if printCDashUrl:
    print("Getting data from:\n\n  " + cdashQueryUrl )
  fullCDashIndexBuilds = extractCDashApiQueryData_in(cdashQueryUrl)
  summaryCDashIndexBuilds = getCDashIndexBuildsSummary(fullCDashIndexBuilds)
  # Determine pass/fail
  (cdashIndexBuildsPassAndExpectedExist_pass, errMsg) = \
    cdashIndexBuildsPassAndExpectedExist(summaryCDashIndexBuilds, expectedBuildNames)
  if not cdashIndexBuildsPassAndExpectedExist_pass:
    return (False, errMsg)
  return (True, "")
    



