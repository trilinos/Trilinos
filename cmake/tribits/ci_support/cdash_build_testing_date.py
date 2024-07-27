#!/usr/bin/env python

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

import datetime

#
# Help message
#


usageHelp = r"""cdash_testing_date.py --cdash-project-start-time="hh:mm" [other options]

Returns to STDOUT the date string YYYY-MM-DD that corresponds the the matching
CDash test day.  Is should be compatible with the 'date=YYYY-MM-DD' argument
to many CDash PHP pages for a build that has a starting time stamp that
matches the requested date.

This tool is meant to match the CDash logic for what testing day a given CDash
build matches for the CDash 'date' field in various PHP pages for a given
CDash build with a give build start time (as reported to CDash from CTest).

For example, if the CDash project start time is 04:00 UTC and the CDash
build start time is 2018-01-22T00:00 UTC, then the relative CDash build date
time would be 2018-01-21T20:00 UTC which is the previous testing day
2018-01-21.  But if the CDash build start time was 2018-01-22T04:10 UTC,
then the relative CDash build date time would be 2018-01-22:00:10 UTC and
the testing day would be 2018-01-22.
"""

#
# Helper functions
#


def injectCmndLineOptionsInParser(clp, gitoliteRootDefault=""):
  
  clp.add_option(
    "--cdash-project-start-time", dest="cdashProjectStartTimeUtcStr", type="string", default="",
    help="Starting time for the CDash testing day in 'hh:mm' in UTC."\
      + " Check the CDash project settings for the testing day start time." )
  
  clp.add_option(
    "--day-incr", dest="dayIncrInt", type="int", default="0",
    help="Increment for the testing date (default '0') [optional]" )
  
  clp.add_option(
    "--cdash-build-start-time", dest="cdashBuildStartTime", type="string", default="",
    help="CDash build start time in format 'YYYY-MM-DDThh:mm UTC'.  If empty ''"\
      +" then the current date/time in UTC is used.  (default '') [optional]" )
  
  clp.add_option(
    "--debug-level", dest="debugLevel", type="int", default="0",
    help="Debug level.  An integer >= 0 (default '0')" )


def getCmndLineOptions():
  from optparse import OptionParser
  clp = OptionParser(usage=usageHelp)
  injectCmndLineOptionsInParser(clp)
  (options, args) = clp.parse_args()
  if options.cdashProjectStartTimeUtcStr == "":
    raise Exception("Error, input argument --cdash-project-start-time must be non"\
      +"-empty and must be of the format hh:mm UTC")
  return options


# Return the current time in UTC as a datetime object
def getCurrentDateTimeUtc():
  return datetime.datetime.utcnow()


# Get the timezone offset as a timedelta object w.r.t to UTC
#
# The supported timezones for timeZoneStr are the strings:
#
# * UTC: 0
# * EDT: 4
# * EST: 5
# * CDT: 5
# * CST: 6
# * MDT: 6
# * MST: 7
#
# NOTE: Any timezone that CDash returns for the 'buildstarttime' field must be
# added below.
#
def getTimeZoneOffset(timeZoneStr):
  if timeZoneStr == "UTC": timezoneOffsetInt = 0
  elif timeZoneStr == "EDT": timezoneOffsetInt = 4
  elif timeZoneStr == "EST": timezoneOffsetInt = 5
  elif timeZoneStr == "CDT": timezoneOffsetInt = 5
  elif timeZoneStr == "CST": timezoneOffsetInt = 6
  elif timeZoneStr == "MDT": timezoneOffsetInt = 6
  elif timeZoneStr == "MST": timezoneOffsetInt = 7
  else: raise Exception("Error, unrecognized timezone '"+timeZoneStr+"'!")
  return datetime.timedelta(hours=timezoneOffsetInt)


# Return a timezone aware datetime object given an input date and time given
# in the format "<YYYY>-<MM>-<DD>T<hh>:<mm>:<ss> <TZ>".
#
# Note. the timzone <TZ> can be any of those supported by the function
# getTimeZoneOffset()
def getBuildStartTimeUtcFromStr(buildStartTimeStr):
  buildStartTimeStrArray = buildStartTimeStr.split(" ")
  if len(buildStartTimeStrArray) == 2:
    timezoneOffset = getTimeZoneOffset(buildStartTimeStrArray[1])
  else:
    timezoneOffset = 0
  localDateTime = datetime.datetime.strptime(buildStartTimeStrArray[0], "%Y-%m-%dT%H:%M:%S")
  return localDateTime + timezoneOffset 


# Return a timedelta object for the CDash Project start time passed in as a
# string in the format "<hh>:<mm>" in UTC.
def getProjectTestingDayStartTimeDeltaFromStr(cdashProjectStartTimeUtcStr):
  t = datetime.datetime.strptime(cdashProjectStartTimeUtcStr, '%H:%M')
  return datetime.timedelta(hours=t.hour, minutes=t.minute, seconds=t.second)


# Return a timedelta object for a day increment pass in as an signed integer.
def getDayIncrTimeDeltaFromInt(dayIncrInt):
  return datetime.timedelta(days=dayIncrInt)


# Return the string "YYYY-MM-DDThh:mm:ss UTC" corresponding to the input
# datetime object.
def getBuildStartTimeUtcStrFromUtcDT(datetimeUtcDT, removeSpace=False):
  dateTimeStr = datetimeUtcDT.strftime("%Y-%m-%dT%H:%M:%S UTC")
  if removeSpace:
    return dateTimeStr.replace(" ", "")
  return dateTimeStr


# Return the string "<YYYY>-<MM>-<DD>" for an input datetime object.
def getDateStrFromDateTime(dateTime):
  return dateTime.strftime("%Y-%m-%d")


# Return the datetime object for just <YYYY>-<MM>-<DD> for an input datetime
# object.
def getDateOnlyFromDateTime(dateTime):
  return datetime.datetime(year=dateTime.year, month=dateTime.month, day=dateTime.day)


# Compute the CDash testing day for a given build using its 'buildstarttime'
# field and the testing day start time and return it as a datetime object.
#
# buildStartTimeStr [in]: The 'buildstarttime' field string as returned from
# CDash in the format "<YYYY>-<MM>-<DD>T<hh>:<mm>:<ss> <TZ>".  Note. the
# timzone <TZ> can be any of those supported by the function
# getTimeZoneOffset()
#
# testingDayStartTimeUtcTD [in]: The testing day start time as the a timedelta
# object.
#
# ToDo: Take into account the CDash server timezone which determine what noon
# means!
# 
def getTestingDayDateFromBuildStartTimeDT(
  buildStartTimeStr, testingDayStartTimeUtcTD,
  ):

  #print()
  #print("testingDayStartTimeUtcTD = "+str(testingDayStartTimeUtcTD))
  #print("testingDayStartTimeUtcTD.seconds = "+str(testingDayStartTimeUtcTD.seconds))

  # Constants
  oneDayTD = datetime.timedelta(days=1)
  noonTD = getProjectTestingDayStartTimeDeltaFromStr("12:00")

  # Convert input 'buildstartime' to datetime object in UtC
  buildStartTimeUtcDT = getBuildStartTimeUtcFromStr(buildStartTimeStr)

  # Get the date and time separately for the input 'buildstartime'
  buildStartTimeDateDT = datetime.datetime(
    year=buildStartTimeUtcDT.year,
    month=buildStartTimeUtcDT.month,
    day=buildStartTimeUtcDT.day,
    )
  #print("buildStartTimeDateDT = "+str(buildStartTimeDateDT))
  buildStartTimeTimeTD = buildStartTimeUtcDT - buildStartTimeDateDT
  #print("buildStartTimeTimeTD = "+str(buildStartTimeTimeTD))
  #print("buildStartTimeTimeTD.seconds = "+str(buildStartTimeTimeTD.seconds))

  if buildStartTimeTimeTD.seconds < testingDayStartTimeUtcTD.seconds:
    buildStartTimeDateDT -= oneDayTD
    #print("buildStartTimeDateDT = "+str(buildStartTimeDateDT))

  if testingDayStartTimeUtcTD > noonTD:
    buildStartTimeDateDT += oneDayTD
    #print("buildStartTimeDateDT = "+str(buildStartTimeDateDT))

  return buildStartTimeDateDT


# Compute the CDash testing day for a given build using its 'buildstarttime'
# field and the testing day start time and return as a string "YYYY-MM-DD".
#
# See getTestingDayDateFromBuildStartTimeDT()
#
def getTestingDayDateFromBuildStartTimeStr(
  buildStartTimeStr, testingDayStartTimeUtcTD,
  ):
  return getDateStrFromDateTime(
    getTestingDayDateFromBuildStartTimeDT(buildStartTimeStr, testingDayStartTimeUtcTD) )


# Return the shifted CDash build start relative to the given CDash project
# testing day start time (as configured on CDash).
#
# This function is meant to match the CDash logic for what testing day a given
# CDash build matches for the CDash 'date' field in various PHP pages for a
# given CDash build with a give build start time (as reported to CDash from
# CTest).
#
# For example, if the CDash project start time is 04:00 UTC and the CDash
# build start time is 2018-01-22T00:00 UTC, then the relative CDash build date
# time would be 2018-01-21T20:00 UTC which is the previous testing day
# 2018-01-21.  But if the CDash build start time was 2018-01-22T04:10 UTC,
# then the relative CDash build date time would be 2018-01-22:00:10 UTC and
# the testing day would be 2018-01-22.
#
# To extract just the year, month and day from the returned datetime object as
# a datetime object, use the function getDateOnlyFromDateTime().
#
# To extract the "<YYYY>-<MM>-<DD>" string from the returned datetime object,
# use the function getDateStrFromDateTime().
# 
def getRelativeCDashBuildStartTime(
  cdashBuildStartTimeUtc,   # CDash build start time in UTC (datetime object)
  cdashProjectStartTimeUtc, # CDash start time in UTC (timedelta object)
  dayIncr,                  # Day increment in past or future (timedelta object)
  ):
  relativeCDashBuildDateTime = \
    cdashBuildStartTimeUtc - cdashProjectStartTimeUtc + dayIncr
  return relativeCDashBuildDateTime


# Return the shifted CDash build start relative to the given CDash project
# testing day start time given input from the command-line arguments.
#
# This function just converts the input command-line args and then calls
# getRelativeCDashBuildStartTime().
#
# Arguments:
#
# cdashBuildStartTimeStr [in] Reference build start time in
# "YYY-MM-DDThh:mm:ss TZ".  If 'None', then is taken from
# getCurrentDateTimeUtc().
#
# cdashProjectStartTimeUtcStr [in] CDash project build start time in
# UTC in the format "hh:mm".
#
# dayIncrInt [in] Integer for the day to return relative to
# cdashProjectStartTimeUtcStr.  The current testing day would be 0.  Yesterday
# would be -1 and so on.
#
# debugLevel [in] If 0, then no debugging.  If > 0, then debug output will be
# printed.
#
def getRelativeCDashBuildStartTimeFromCmndLineArgs(
  cdashBuildStartTimeStr,
  cdashProjectStartTimeUtcStr,
  dayIncrInt,
  debugLevel=0,
  ):

  if cdashBuildStartTimeStr:
    buildStartTimeUtc = getBuildStartTimeUtcFromStr(cdashBuildStartTimeStr)
  else:
    buildStartTimeUtc = getCurrentDateTimeUtc()
  if debugLevel: print("buildStartTimeUtc = "+str(buildStartTimeUtc))

  cdashStartTime = \
    getProjectTestingDayStartTimeDeltaFromStr(cdashProjectStartTimeUtcStr)
  if debugLevel: print("cdashStartTime = "+str(cdashStartTime))

  dayIncr = getDayIncrTimeDeltaFromInt(dayIncrInt)
  if debugLevel: print("dayIncr = "+str(dayIncr))

  return getRelativeCDashBuildStartTime(buildStartTimeUtc, cdashStartTime, dayIncr)


# Class to encapsulate handling of CDash testing day logic
#
class CDashProjectTestingDay(object):

  # Construct
  #
  # currentTestingDayDateStr [in]: The current project CDash testing date in
  # string format "YYYY-MM-DD".
  #
  # projectTestingDayStartTimeUtcStr [in]: The CDash projects's testing day start
  # time in UTC.  This is a string in the format "hh:mm".
  #
  def __init__(self, currentTestingDayDateStr, projectTestingDayStartTimeUtcStr):
    # Store input args
    self.__currentTestingDayDateStr = currentTestingDayDateStr
    self.__projectTestingDayStartTimeUtcStr = projectTestingDayStartTimeUtcStr
    # Set the start of the testing day in UTC
    self.__currentTestingDayDateDT = \
      datetime.datetime.strptime(currentTestingDayDateStr, "%Y-%m-%d")
    self.__projectTestingDayStartTimeUtcTD = \
      getProjectTestingDayStartTimeDeltaFromStr(projectTestingDayStartTimeUtcStr)
    noonTD = getProjectTestingDayStartTimeDeltaFromStr("12:00")
    oneDayTD = datetime.timedelta(days=1)
    if self.__projectTestingDayStartTimeUtcTD >= noonTD:
      self.__testingDayStartDateTimeUtcDT = \
         self.__currentTestingDayDateDT - oneDayTD + self.__projectTestingDayStartTimeUtcTD
    else:
      self.__testingDayStartDateTimeUtcDT = \
        self.__currentTestingDayDateDT + self.__projectTestingDayStartTimeUtcTD

  # Return the testing day start in UTC as a datetime object
  def getTestingDayStartUtcDT(self):
    return self.__testingDayStartDateTimeUtcDT

  # Return the testing day date (not time) as a datetime object
  def getCurrentTestingDayDateDT(self):
    return self.__currentTestingDayDateDT

  # Return the testing day as a datetime object for the input 'buildstartime'
  # string.
  def getTestingDayDateFromBuildStartTimeDT(self, buildStartTimeStr):
    return getTestingDayDateFromBuildStartTimeDT(
      buildStartTimeStr, self.__projectTestingDayStartTimeUtcTD)

  # Return the testing day string "YYYY-MM-DD" for the input 'buildstartime'
  # string.
  def getTestingDayDateFromBuildStartTimeStr(self, buildStartTimeStr):
    return getDateStrFromDateTime(
       self.getTestingDayDateFromBuildStartTimeDT(buildStartTimeStr) )


#
# Run the script
#

if __name__ == '__main__':

  inOptions = getCmndLineOptions()

  relativeCDashBuildBuildDateTime = getRelativeCDashBuildStartTimeFromCmndLineArgs(
    inOptions.cdashBuildStartTime, inOptions.cdashProjectStartTimeUtcStr,
    inOptions.dayIncrInt, inOptions.debugLevel,
    )

  print(getDateStrFromDateTime(relativeCDashBuildBuildDateTime))
