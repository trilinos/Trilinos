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
    # regular expression for timing output
    #regex="[0-9]*[.][0-9e-]*";
    regex="[0-9]+[.]?[0-9]*[e]?[-]?[0-9]* [(][0-9][)]";
    startParsingTimers=0;
    startParsingPC=0;
}

###############################################################################
#stuff that happens when processing all files
{
    #fix minor difference between Epetra and Tpetra smoother tags
    if (match($0,"Linear algebra library: ")) {
      after = substr($0,RSTART+RLENGTH);
      linalg[FILENAME] = after;
    }

    #indicates start of MueLu timer output
    if (match($0,"^Timer Name")) {
      startParsingTimers=1;
    }

    if (startParsingTimers && match($0,"^MueLu: ")) {

      # level-specific timers that do not include calls to child factories
      if (match($0,"[(]level=[0-9][)]")) {
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

    # Pull out the reported total setup time.  This is printed as a sanity check.
    foundTotal=0;
    if (startParsingTimers && match($0,"^ScalingTest: 2 - MueLu Setup")) {
      foundTotal=1;
      #match($0,regex);
      #TotalSetup[linalg[FILENAME]] = substr($0,RSTART,RLENGTH);
      alltimes = substr($0,RSTART+RLENGTH);
      cutCmd="cut -f3 -d')' | cut -f1 -d'('"
      TotalSetup[linalg[FILENAME]] = ExtractTime(alltimes,cutCmd);
    }
    if (foundTotal==0) {
      #This duplicates the previous if-block.  Only the regex in the "if" is different.
      if (startParsingTimers && match($0,"^MueLu: Hierarchy: Setup [(]total[)]")) {
        alltimes = substr($0,RSTART+RLENGTH);
        cutCmd="cut -f3 -d')' | cut -f1 -d'('"
        TotalSetup[linalg[FILENAME]] = ExtractTime(alltimes,cutCmd);
      }
    }
}

###############################################################################
function PrintHeader(description,linalg)
{
  space = " ";
  printf("%80s      ",toupper(description));
  for (j in linalg) {
    printf("%10s  ",linalg[j]);
    printf(" (total)");
  }
  printf("\n%80s          ------------------\n",space);
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
  foo = 0;
  arrayLeng=0
  for (i in arrayToSort) arrayLeng++; #length(arrayToSort) doesn't work here for some reason
  for (i=1; i<arrayLeng+1; i++) {
    ind = arrayToSort[i];
    printf("%80s  ==> ",labels[ind]);
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
  printf("%80s          ----------------\n%83s"," "," ");
  for (i in linalg)
    printf("Hierarchy Setup: %5.2f seconds       ",TotalSetup[linalg[i]]);
    #printf("%77s%-17s%-17s%-17s%-11s\n"," "," "," ","MueLu reported",TotalSetup[linalg[j]] " sec.");
}
###############################################################################

function ExtractTime(alltimes,cutCmd)
{
  alltimes = RemoveWhiteSpace(alltimes);
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

  printLevelStats = 1;
  if (printLevelStats) {
    PrintHeader("Setup times (level specific) excluding child calls ",linalg);
    SortAndReportTimings(linalg[FILENAME],setupLabels,setupTimes,linalg);
    PrintTotalSetupTime();
    #PrintHeader("Solve times ",linalg);
    #SortAndReportTimings(linalg[FILENAME],solveLabels,solveTimes,linalg);
    printf("\n\n");
  }

}
