# awk script to provide context for MueLu timers
#
# Don't run this script directly, instead use compareTE.sh.
#
# Note: This requires GNU awk and won't run in compatibility mode.
#
# Note: You might need to define AWKPATH if mueprof.sh can't find this script.
#
# Thanks to http://www.grymoire.com/Unix/Awk.html and CMS for help.

###############################################################################
#stuff that happens before processing any files
BEGIN {
    whichBlockToAnalyze=blockNumber; #passed in from shell script
    foundTimerBlock=0;
    foundTimersToReport=0;
    possibleTotalLabels[1] = "ScalingTest: 2 - MueLu Setup";
    possibleTotalLabels[2] = "MueLu: Hierarchy: Setup [(]total[)]";
    possibleTotalLabels[3] = "nalu MueLu preconditioner setup";
}

###############################################################################
#stuff that happens when processing all files
{
    if (match($0,"Linear algebra library: ")) {
      after = substr($0,RSTART+RLENGTH);
      linalg[FILENAME] = after;
    }

    #indicates start of MueLu timer output
    if (match($0,"^Timer Name")) {
      #logic to ensure only the ith solver block is analyzed
      whichBlockToAnalyze--;
      if (whichBlockToAnalyze>0)
        foundTimerBlock=0;
      if (whichBlockToAnalyze==0)
        foundTimerBlock=1;
      if (whichBlockToAnalyze<0)
        nextfile; #we've moved past the solver block of interest,
                  #stop analyzing this file
    }

    if (foundTimerBlock) {
      if (match($0,"^MueLu: ")) {
        # matched a MueLu timer
        if (match($0,"[(]level=[0-9][)]")) {
          # timer is level-specific (and by its nature excludes calls to child factories)
          foundTimersToReport=1;
          factAndLevel = substr($0,1,RSTART-1+RLENGTH);
          alltimes = substr($0,RSTART+RLENGTH);
          cutCmd="cut -f3 -d')' | cut -f1 -d'('"
          maxtime = ExtractTime(alltimes,cutCmd);
          if (match(factAndLevel,"MueLu: Hierarchy: Solve")) {
            #TODO figure out which solve labels to pull out
            solveLabels[factAndLevel] = factAndLevel;
            solveTimes[factAndLevel,linalg[FILENAME]] = maxtime;
          } else {
            setupLabels[factAndLevel] = factAndLevel;
            setupTimes[factAndLevel,linalg[FILENAME]] = maxtime;
          }
        }
      }

      # Check for any reported total setup time.  This is printed as a sanity check
      # against the running total.
      for (i in possibleTotalLabels) {
        if (match($0,possibleTotalLabels[i])) {
          pattern = substr($0,RSTART,RLENGTH);
          alltimes = substr($0,RSTART+RLENGTH);
          cutCmd="cut -f3 -d')' | cut -f1 -d'('"
          TotalSetup[pattern,linalg[FILENAME]] = ExtractTime(alltimes,cutCmd);
        }
      }
    } #if (foundTimerBlock)

}

###############################################################################
function PrintHeader(description,linalg)
{
  space = " ";
  printf("%60s      ",toupper(description));
  for (j in linalg) {
    printf("%10s  ",linalg[j]);
    printf(" (total)");
  }
  printf("\n%60s          ------------------\n",space);
}

###############################################################################
# Reorder all the timings by sorting either the Epetra or Tpetra timing values
# in ascending order.
function SortAndReportTimings(libToSortOn,mylabels,tallies,linalg)
{
  numInds=1;
  for (i in mylabels) {
    # the plus 0 is just a trick to make asort treat "talliesToSort" as numeric
    talliesToSort[numInds] = tallies[i,libToSortOn]+0;
    newlabels[talliesToSort[numInds]] = mylabels[i];
    for (j in linalg) {
      newtallies[talliesToSort[numInds],linalg[j]] = tallies[i,linalg[j]];
    }
    numInds++;
  }

  # Sort timings for library "libToSortOn" in ascending order.
  # These sorted timings are used as the indices into newlabels and newtallies.
  SortAndPrint(talliesToSort,newlabels,newtallies,linalg);

  #awk has weird scoping rules
  delete newlabels
  delete newtallies
  delete talliesToSort
}

###############################################################################
# helper function
function SortAndPrint(arrayToSort,labels,tallies,linalg)
{
  asort(arrayToSort);
  for (j in linalg) {
    runningTotals[j] = 0;
  }
  arrayLeng=0
  for (i in arrayToSort) arrayLeng++; #length(arrayToSort) doesn't work here for some reason
  for (i=1; i<arrayLeng+1; i++) {
    ind = arrayToSort[i];
    printf("%60s  ==> ",labels[ind]);
    for (j in linalg) {
      runningTotals[j] += tallies[ind,linalg[j]];
      printf("%10.4f   (%5.2f)",tallies[ind,linalg[j]],runningTotals[j]);
    }
    if (tallies[ind,"Epetra"] != 0) {
      printf("%8.2f   ",tallies[ind,"Tpetra"] / tallies[ind,"Epetra"]);
    }
    printf("\n");
  }
}

###############################################################################
# For sanity purposes, print the total that MueLu prints
function PrintTotalSetupTime()
{
  printf("%60s          ----------------\n"," ");
  for (i in TotalSetup) {
    split(i,sep,SUBSEP); #breaks multiarray index i up into its constituent parts.
                         #we only want the first one, sep[1]
    printf("%60s  ==>   %6.3f seconds       \n",sep[1],TotalSetup[i]);
  }
}
###############################################################################

function ExtractTime(alltimes,cutCmd)
{
  alltimes = RemoveWhiteSpace(alltimes);
  # system call to cutCmd
  # see GNU Awk manual section 12.3, Two-Way Communications with Another Process
  print alltimes |& cutCmd
  close(cutCmd,"to")
  cutCmd |& getline themax
  close(cutCmd)
  themax = RemoveWhiteSpace(themax);
  return (themax);
}

###############################################################################

function RemoveWhiteSpace(theString)
{
  sub(/^[ ]*/,"",theString); sub(/[ ]*$/,"",theString);
  return (theString)
}
###############################################################################

###############################################################################
#stuff that happens after processing all files
###############################################################################
END {

  if (foundTimersToReport) {
    PrintHeader("Setup times (level specific) excluding child calls ",linalg);
    SortAndReportTimings(linalg[FILENAME],setupLabels,setupTimes,linalg);
    PrintTotalSetupTime();
    #PrintHeader("Solve times ",linalg);
    #SortAndReportTimings(linalg[FILENAME],solveLabels,solveTimes,linalg);
  } else {
    printf("\nFound only %d solver blocks, you requested block %d.\n",blockNumber-whichBlockToAnalyze,blockNumber);
  }
  printf("\n");

}
