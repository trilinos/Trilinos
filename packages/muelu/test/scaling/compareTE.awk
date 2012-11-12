# awk script to compare Epetra and Tpetra timings from MueLu_ScalingTest.exe side-by-side.
#
# syntax:
#    awk -f compareTE.awk file1 file2
#
# where file1 and file2 are epetra/tpetra screen dumps.  file1 and file2 can have any name you like, and order doesn't matter.
#
# You should see something like the following:
#                                                                                          Tpetra       Epetra    T/E ratio
#                                                                                          --------------------------------
#                          MueLu: IfpackSmoother: Setup Smoother (total, level=0)  ==>     0.0325       0.0119       2.72
#                MueLu: AmalgamationFactory: AmalgamationFactory (total, level=2)  ==>     0.0173       0.0130       1.34
#                              MueLu: TransPFactory: Transpose P (total, level=4)  ==>     0.0021       0.0006       3.70
#                         MueLu: AmalgamationFactory: AmalgamationFactory (total)  ==>     0.3636       0.3431       1.06
#                                   MueLu: AmesosSmoother: Setup Smoother (total)  ==>     0.0001       0.0008       0.10
#
# The first two columns are times in seconds, the last column is the ratio of the first two columns.
#
# Thanks to http://www.grymoire.com/Unix/Awk.html and CMS for help.

#stuff that happens before processing any files
BEGIN {
    # regular expression for timing output
    regex="[0-9]*[.][0-9e-]*";
}

#stuff that happens when processing all files
{
    #fix minor difference between Epetra and Tpetra smoother tags
    sub(/Ifpack2Smoother/,"IfpackSmoother");
    sub(/Amesos2Smoother/,"AmesosSmoother");
    if (match($0,"Linear algebra library: ")) {
      after = substr($0,RSTART+RLENGTH);
      linalg[FILENAME] = after;
    }
    if (match($0,"[(]total") && match($0,regex)) {
      #RSTART is where the pattern starts
      #RLENGTH is the length of the pattern
      before = substr($0,1,RSTART-1);
      #remove trailing white space in variable "before"
      sub(/[ ]*$/,"",before);
      pattern = substr($0,RSTART,RLENGTH);
      after = substr($0,RSTART+RLENGTH);
      labels[before] = before;
      tallies[before,linalg[FILENAME]] = pattern;
      foundATotal=0;
    }
}

#stuff that happens after processing all files
END {
  space = " ";
  printf("%80s      ",space);
  for (j in linalg) {
    printf("%10s   ",linalg[j]);
  }
  printf("%10s   ","T/E ratio");
  printf("\n%80s          --------------------------------\n",space);
  for (i in labels) {
    printf("%80s  ==> ",labels[i]);
    for (j in linalg) {
      printf("%10.4f   ",tallies[i,linalg[j]]);
    }
    if (tallies[i,"Epetra"] != 0) {
      printf("%8.2f   ",tallies[i,"Tpetra"] / tallies[i,"Epetra"]);
    }
    printf("\n");
  }
}
