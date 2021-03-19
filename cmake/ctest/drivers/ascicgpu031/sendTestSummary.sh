#!/bin/sh

# Parse command line options.
DEBUGMODE=0
USAGE="sendTestSummary.sh [-d] <logfile>"
while getopts d OPT; do
    case "$OPT" in
        d)
            # debug mode, send email summary to me only
            DEBUGMODE=1
            ;;
        \?)
            # getopts issues an error message
            echo ${USAGE}
            echo
            exit 1
            ;;
    esac
done

# Remove the options we parsed above.
shift `expr $OPTIND - 1`
# The logfile is required. Error out if it's not provided.
if [ $# -eq 0 ]; then
  echo $USAGE >&2
  exit 1
fi
### end parsing ###

#########################################################################
# Variables you might want to modify.
#########################################################################

#Perl script to produce prettified HTML
HTMLPERLSCRIPT="/ascldap/users/csiefer/Trilinos/ascicgpu-testing/Trilinos/cmake/ctest/drivers/ascicgpu031/drakify-email.pl"

#root of cdash testing directory
TESTLOCATION="/home/csiefer/Trilinos/ascicgpu-testing"
LOGBACKUPDIRECTORY="/home/csiefer/Trilinos/ascicgpu-testing/logs"

#packages to be summarized
PATTERN="(Amesos2|Ifpack2|Xpetra|MueLu|Tpetra|Zoltan2)"

#variables to be passed to the perl script
MACHINENAME=`hostname -s`
USER=`whoami`

#who gets the email summary
if [[ $DEBUGMODE == 1 ]]; then
  RECIPIENTS=(
    "${USER}@sandia.gov"
  )
else
  RECIPIENTS=(
    "tpetra-regression@software.sandia.gov"
  )
fi
#suffix for all the log files
timeStamp="$(date +%F_%R)"

#cron driver log file
INFILE=$1
#root of file to be emailed.  The correct suffix must be appended whenever you use this.
OUTFILE="test-summary-${timeStamp}"
MAILCOMMAND="/usr/sbin/sendmail"
#########################################################################

cd ${TESTLOCATION}

backupFile="cron_driver.log.$timeStamp"
cp cron_driver.log $backupFile

testStartString=`egrep "Starting nightly Trilinos development" cron_driver.log`
testStartDate=`echo $testStartString | sed "s/:/#/" | cut -f 2 -d#`
ttt=`echo $testStartString | cut -f 1 -d:`
testMachine=${ttt##* }
testEndString=`egrep "Ending nightly Trilinos development" cron_driver.log`
testEndDate=`echo $testEndString | sed "s/:/#/" | cut -f 2 -d#`

awk -v packagesToMatch="$PATTERN" -v summaryFile="${OUTFILE}.txt" -v machine="$testMachine" -v startTime="$testStartDate" -v endTime="$testEndDate" '

###################################################
# Commands to run before the file is processed
###################################################
BEGIN {
  print "Machine    :  " machine   > summaryFile
  print "Start time : " startTime > summaryFile
  print "End time   : " endTime   > summaryFile
  testctr=0
  gitUpdateFailed=0
  dashboardErrors=0
}

###################################################
# Commands to run while processing the file
###################################################
{

  if ($0 ~ "Update command failed")
  {
    gitUpdateFailed=1
  }

  if ($0 ~ "^test [0-9]*$")
  {
    #start of test found, e.g., test 4
    FOUND=2
    testNum=$0
    sub(/test /,"",testNum)
    testNum=testNum":"
    #\x27 is hex code for single quote
    packageLibBuild=testNum" Building target: \x27" packagesToMatch "_libs"
    packageTestBuild=testNum" Build ALL target for \x27" packagesToMatch "\x27"
    runTestPattern=testNum" Running test for package \x27" packagesToMatch "\x27"
    next #skip any more processing, go on to next line
  }

  if (FOUND==2)
  {
    FOUND--
    dashboardName=$0
    sub(/Start [ ]*[0-9]*: /,"",dashboardName)
    dashboardName=RemoveWhiteSpace(dashboardName)
    listOfDashboardNames[testctr] = dashboardName
    testctr++
    dashBoardPattern="Test [ ]*#[0-9]*: " dashboardName
  }

# Record the "track" for this dashboard, which could be "Nightly", "Experimental", or "Specialized"
  if (FOUND && $0 ~ "-- Trilinos_TRACK=")
  {
    thisLine=$0 
    sub(/^[0-9]*: -- Trilinos_TRACK=\x27/,"",thisLine)
    sub(/\x27/,"",thisLine)
    if (length(thisLine) > 0) {
      dashboardTrack[dashboardName] = thisLine
      trackTypes[thisLine]++
    }
  }

  if (FOUND && $0 ~ dashBoardPattern)
  {
    thisLine=$0
    thisLine=RemoveWhiteSpace(thisLine)
    if (dashboardErrors == 0)
      dashBoardSummary[dashboardName] = "passed"
    else
      dashBoardSummary[dashboardName] = "FAILED"
    match(thisLine,"[0-9]*\\.[0-9]* sec$")
    timeSummary[dashboardName] = substr(thisLine,RSTART,RLENGTH)
    #done with this dashboard, reset error flag
    dashboardErrors=0
  }

  # library build
  if (FOUND && $0 ~ packageLibBuild)
  {
    getCompilerSummary=2
    thisLine = $0
    pat = "\x27" packagesToMatch "_libs\x27"
    match(thisLine,pat)
    currentPackage = substr(thisLine,RSTART+1,RLENGTH-2)
    listOfPackages[currentPackage] = currentPackage
  }

  # tests build
  if (FOUND && $0 ~ packageTestBuild)
  {
    getCompilerSummary=2
    thisLine = $0
    pat = "\x27" packagesToMatch "\x27"
    match(thisLine,pat)
    currentPackage = substr(thisLine,RSTART+1,RLENGTH-2)
    listOfPackages[currentPackage] = currentPackage
  }

  if (getCompilerSummary>0 && $0 ~ " Compiler errors")
  {
    thisLine=$0
    sub(testNum,"",thisLine)
    pat="[0-9]*"
    thisLine=RemoveWhiteSpace(thisLine)

    match(thisLine,pat)
    numErrors = substr(thisLine,RSTART,RLENGTH)
    errorSummary[dashboardName,currentPackage] = numErrors

    getCompilerSummary--
  }

  if (getCompilerSummary>0 && $0 ~ " Compiler warnings")
  {
    thisLine=$0
    sub(testNum,"",thisLine)
    thisLine=RemoveWhiteSpace(thisLine)
    pat="[0-9]*"
    match(thisLine,pat)
    numWarnings = substr(thisLine,RSTART,RLENGTH);
    warningSummary[dashboardName,currentPackage] = numWarnings
    getCompilerSummary--
  }

  #Look for pattern indicating that the tests of interest have in fact run.
  if (FOUND && match($0,runTestPattern))
  {
    packageTested = substr($0,RSTART,RLENGTH);
    sub(testNum,"",packageTested)
    packageTested=RemoveWhiteSpace(packageTested)
    getTestSummary=1
  }

  if (getTestSummary && $0 ~ "No tests were found!!!")
  {
    getTestSummary=0
  }

  #Calculate the number of failing, passing, and total tests.
  if (getTestSummary && $0 ~ "tests failed out of")
  {
    thisLine=$0
    sub(testNum,"",thisLine)
    thisLine=RemoveWhiteSpace(thisLine)
    getTestSummary=0
    pat = "[0-9]* tests failed out of [0-9]*"
    match(thisLine,pat)
    ttt = substr(thisLine,RSTART,RLENGTH);
    pat = "^[0-9]*"
    match(ttt,pat)
    numFailed = substr(ttt,RSTART,RLENGTH);
    pat = "[0-9]*$"
    match(ttt,pat)
    numTotal = substr(ttt,RSTART,RLENGTH);
    failSummary[dashboardName,currentPackage] = numFailed
    passSummary[dashboardName,currentPackage] = numTotal+0-numFailed
    totalSummary[dashboardName,currentPackage] = numTotal
    if (numFailed != 0)
      dashboardErrors=1
  }
}

###################################################
# helper functions
###################################################
function RemoveWhiteSpace(theString)
{
  sub(/^[ ]*/,"",theString); sub(/[ ]*$/,"",theString);
  return (theString)
}

###################################################
# Commands to run after the file is processed
###################################################
END {

  if (gitUpdateFailed == 1) {
    print "\n *** git update FAILED ***\n" > summaryFile
  }

  # do some nice formatting
  numPlusses=73
  thePluses="  "
  while (jj++<numPlusses)
    thePluses = thePluses "+"

  printf("\n---------------------------------- Summary ----------------------------------\n") > summaryFile
  for (track in trackTypes) {
    printf("%s\n",thePluses) > summaryFile
    trackNameLength = length(track)
    numPlussesToTheRight = numPlusses - trackNameLength - 4
    plussesToTheRight=""
    kk = 0
    while (kk++<numPlussesToTheRight)
      plussesToTheRight = plussesToTheRight "+"
    printf("  ++ %*-s %s\n",trackNameLength,track,plussesToTheRight) > summaryFile
    printf("%s\n",thePluses) > summaryFile
    for (i in listOfDashboardNames) {
      db=listOfDashboardNames[i]
      if (dashboardTrack[db] == track)
        printf("  %61-s  ... %s\n",db,dashBoardSummary[db]) > summaryFile;
    }
  }
  printf("-----------------------------------------------------------------------------\n\n") > summaryFile

  for (i in listOfDashboardNames) {
    db=listOfDashboardNames[i]
    spaces="     "
    printf("%55-s\n%s%8-s, %5.1f seconds\n",db, spaces, dashBoardSummary[db], timeSummary[db]) > summaryFile;
    for (k in listOfPackages) {
      pat = "_lib"
      if (match(k,pat)) isLib = 1;
      else              isLib = 0;
      if ((db,k) in warningSummary) nwarn = warningSummary[db,k]
      else                                nwarn = "-";
      if ((db,k) in errorSummary) nerr = errorSummary[db,k]
      else                              nerr = "-";
      if ((db,k) in failSummary) nfail = failSummary[db,k]
      else                             nfail = "-";
      if ((db,k) in passSummary) npass = passSummary[db,k]
      else                             npass = "-";
      if ((db,k) in totalSummary) ntotal = totalSummary[db,k]
      else                              ntotal = "-";
      if (isLib) {
        summaryString = sprintf("%15s | %3d warnings | %3d errors",k,nwarn,nerr);
      }
      else {
        summaryString = sprintf("%15s | %3d warnings | %3d errors | %d/%d passed",k,nwarn,nerr,npass,ntotal);
      }
      print spaces summaryString > summaryFile
    }
  }

}
' $INFILE

date2=`echo $(date) | sed "s/ /_/g"`
cdashDate="$(date +%F)"
cat ${OUTFILE}.txt | perl ${HTMLPERLSCRIPT} ${date2} ${cdashDate} ${MACHINENAME} ${USER} > ${OUTFILE}.html

${MAILCOMMAND} -it <<END_MESSAGE
To: ${RECIPIENTS[@]}
$(cat ${OUTFILE}.html | grep -v UNMATCHED)
END_MESSAGE

#clean up
bzip2 --best $backupFile
mv ${backupFile}.bz2 ${OUTFILE}.txt ${LOGBACKUPDIRECTORY}
rm -f ${OUTFILE}.html
#mv ${OUTFILE}.html  ${LOGBACKUPDIRECTORY}
