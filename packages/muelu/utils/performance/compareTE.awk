# awk script to compare Epetra and Tpetra timings from MueLu_ScalingTest.exe side-by-side.
#
# Don't run this script directly, instead use compareTE.sh.
#
# Note: This requires GNU awk and won't run in compatibility mode.
#
# Note: You might need to define AWKPATH if compareTE.sh can't find this script.
#
# Thanks to http://www.grymoire.com/Unix/Awk.html and CMS for help.

###############################################################################
#stuff that happens before processing any files
BEGIN {
    # regular expression for timing output
    regex="[0-9]*[.][0-9e-]*";
    startParsingTimers=0;
}

###############################################################################
#stuff that happens when processing all files
{
    #fix minor difference between Epetra and Tpetra smoother tags
    sub(/Ifpack2Smoother/,"IfpackSmoother");
    sub(/Amesos2Smoother/,"AmesosSmoother");
    if (match($0,"Linear algebra library: ")) {
      after = substr($0,RSTART+RLENGTH);
      linalg[FILENAME] = after;
    }

    #indicates start of MueLu timer output
    if (match($0,"^Timer Name")) {
      startParsingTimers=1;
    }

    if (startParsingTimers && match($0,"^MueLu: ")) {

        #match time
        match($0,regex);
        #RSTART is where the pattern starts
        #RLENGTH is the length of the pattern
        before = substr($0,1,RSTART-1);
        #remove trailing white space in variable "before"
        sub(/[ ]*$/,"",before);
        pattern = substr($0,RSTART,RLENGTH);
        after = substr($0,RSTART+RLENGTH);

        # totals, not level-specific
        if (match($0,"[(]total[)]")) {
          tlabels[before] = before;
          ttallies[before,linalg[FILENAME]] = pattern;
        }

        # totals by level
        if (match($0,"[(]total, level=")) {
          tllabels[before] = before;
          tltallies[before,linalg[FILENAME]] = pattern;
        }

        # no totals, not level-specific
        if (!match($0,"level") && !match($0,"total")) {
          ntlabels[before] = before;
          nttallies[before,linalg[FILENAME]] = pattern;
        }

        # no totals, level-specific
        if (match($0,"[(]level=")) {
          #print linalg[FILENAME] " | " before " | " pattern " | " after   #debugging
          ntllabels[before] = before;
          ntltallies[before,linalg[FILENAME]] = pattern;
        }

        #subfactory timings, summed  over all levels
        if (match($0,"[(]sub, total[)]")) {
          stlabels[before] = before;
          sttallies[before,linalg[FILENAME]] = pattern;
        }

        #subfactory timings, level-specific
        if (match($0,"[(]sub, total, level")) {
          stllabels[before] = before;
          stltallies[before,linalg[FILENAME]] = pattern;
        }

    }

    # Pull out the reported total setup time.  This is printed as a sanity check.
    if (startParsingTimers && match($0,"^ScalingTest: 2 - MueLu Setup")) {
      match($0,regex);
      TotalSetup[linalg[FILENAME]] = substr($0,RSTART,RLENGTH);
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
  printf("%10s   ","T/E ratio");
  printf("\n%80s          ----------------------------------------------\n",space);
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
# Reorder all the timings by sorting either the Epetra or Tpetra timer numbers
# in ascending order.
function SortAndReportByTimerNumber(libToSortOn,mylabels,tallies,linalg)
{
  numInds=1;
  for (i in mylabels) {
    if (match(i," [0-9]* :")) {
      # pull out the timer number
      #RSTART is where the pattern starts
      #RLENGTH is the length of the pattern
      timerNumber = substr(i,RSTART,RLENGTH);
      sub(/[ ]*/,"",timerNumber); #remove spaces
      sub(/:/,"",timerNumber); #remove colon
      # the plus 0 is just a trick to make asort treat "timerNumbersToSort" as numeric
      timerNumbersToSort[numInds] = timerNumber+0;
      newlabels[timerNumbersToSort[numInds]] = mylabels[i];
      for (j in linalg) {
        newtallies[timerNumbersToSort[numInds],linalg[j]] = tallies[i,linalg[j]];
      }
      numInds++;
    }
  }

  SortAndPrint(timerNumbersToSort,newlabels,newtallies,linalg);

  #awk has weird scoping rules
  delete newlabels
  delete newtallies
  delete timerNumbersToSort
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
  for (i=1; i<arrayLeng; i++) {
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
  printf("%80s          ----------------------------------------------\n%100s"," "," ");
  for (i in linalg)
    printf("%5.2f               ",TotalSetup[linalg[i]]);
}
###############################################################################
#stuff that happens after processing all files
END {

  if (printSummedStats) {
    PrintHeader("Timings including child calls, summed over all levels",linalg);
    SortAndReportTimings("Tpetra",tlabels,ttallies,linalg);
    printf("\n\n");
  }

  if (printLevelStats) {
    PrintHeader("Timings including child calls, level specific",linalg);
    SortAndReportTimings("Tpetra",tllabels,tltallies,linalg);
    printf("\n\n");
  }

  if (printSummedStats) {
    PrintHeader("Timings without child calls, summed over all levels",linalg);
    SortAndReportTimings("Tpetra",ntlabels,nttallies,linalg);
    PrintTotalSetupTime();
    printf("\n\n");
  }

  if (printLevelStats) {
    PrintHeader("Timings without child calls, level specific",linalg);
    SortAndReportByTimerNumber("Tpetra",ntllabels,ntltallies,linalg);
    PrintTotalSetupTime();
    printf("\n\n");
  }

  if (printSummedStats) {
    PrintHeader("Timings for subfactories, summed over all levels",linalg);
    SortAndReportTimings("Tpetra",stlabels,sttallies,linalg);
    printf("\n\n");
  }

  if (printLevelStats) {
    PrintHeader("Timings for subfactories, level specific",linalg);
    SortAndReportTimings("Tpetra",stllabels,stltallies,linalg);
  }
}
