#!/usr/bin/perl

use File::Copy;
use File::Compare;

##############################################################################
##############################################################################
### Remove white-space from a file.
sub nowhite($) {
  my $file = $_[0];
  for ($file) {
    s/ //g;
    s/\t//g;
  }
  return uc($file);
}

##############################################################################
##############################################################################

### Text string needed by CTEST to provide unlimited output.
### Otherwise, output is limited to 50K characters.
print "CTEST_FULL_OUTPUT\n";

### Check command line arguments.
### Usage:  ctest_zoltan.pl #processors package [debug] [mpiexec_path].
foreach $argnum (0 .. $#ARGV) {
   print "$ARGV[$argnum]";
}
print "\n";

$numArgs = $#ARGV + 1;
if ($numArgs < 2) {
  print "Usage:  ctest_zoltan.pl #processors package [debug] [mpiexec_path]\n";
  exit -1;
}

### Get debug level for ctest_zoltan.pl.
$debug = 0;
if ($numArgs > 2) {$debug = $ARGV[2];}

### Get number of processors
$np = $ARGV[0];

### Get the package number indicating which tests to run.
$package = $ARGV[1];
if ($debug) {print "DEBUG:  package $package\n";}

### Get the path to mpiexec if it is specified.
$mpiexec = "mpiexec";
if ($numArgs > 3) {$mpiexec = $ARGV[3];}
$mpiexecargs = "--mca mpi_yield_when_idle 1";
if ($debug) {print "DEBUG:  mpiexec $mpiexec $mpiexecargs\n";}


### Assign the executable.
$zdrive = "../../src/driver/zdrive.exe";

### Get current working directory name
use Cwd;
$dir = getcwd();

### Compute base name for zdrive output files.
### Remove the ch_, hg_ or nem_ prefix.
@tmparr = split('/', $dir);
$dirname = $tmparr[-1];
$dirname =~ s/ch_//g;
$dirname =~ s/hg_//g;
$dirname =~ s/nem_//g;
$zoutfilebase = sprintf("%s.out.%d.", $dirname, $np);
$zoutdropbase = sprintf("%s.drops.%d.", $dirname, $np);
if ($debug) {
  print "DEBUG  Dir $dir dirname $dirname\n";
  print "DEBUG  Outfilebase: $zoutfilebase;  Dropbase: $zoutdropbase\n";
}

### Open a logfile
$zoutlogfile = sprintf("%s.%s.logfile", $dirname, $package);
open(LOG, "> $zoutlogfile");
$time = localtime;
print LOG "$time\n";


### If output subdirectory does not exist, create it.
mkdir "output" unless -d "output";

### Get list of input files
if ($package eq "Zoltan") {
  ###  Zoltan native algorithms
  @inpfiles = glob("zdrive.inp.rcb* 
                    zdrive.inp.rib*
                    zdrive.inp.reftree* 
                    zdrive.inp.hsfc*
                    zdrive.inp.oct*
                    zdrive.inp.phg* 
                    zdrive.inp.color*
                    zdrive.inp.block*
                    zdrive.inp.cyclic*
                    zdrive.inp.random* 
                    zdrive.inp.graph-repartit*
                    zdrive.inp.graph-partit*");
} 
if ($package eq "ParMETIS") {
  ### ParMETIS tests
  if ($np > 1) {
    @inpfiles = glob("zdrive.inp.adp*
                      zdrive.inp.part* 
                      zdrive.inp.*metis*");
  } else {
    ### ParMETIS adaptive methods do not work on one processor.
    @inpfiles = glob("zdrive.inp.part* 
                      zdrive.inp.*metis*");
  }
} 
if ($package eq "Scotch") {
  ### Scotch tests
  @inpfiles = glob("zdrive.inp.*scotch*");
} 
if ($package eq "PaToH") {
  ### PaToH tests
  @inpfiles = glob("zdrive.inp.*patoh*");
}

### Set test counters
$testcnt = 0;
$failcnt = 0;
$passcnt = 0;

### For each zdrive.inp file, run the test.
TEST:  foreach $file (@inpfiles) {
  if ($debug) {print "DEBUG  Running test $testcnt on $file\n";}

  ### Remove zdrive output files from previous runs.
  unlink glob("$zoutfilebase*");
  unlink glob("$zoutdropbase*");

  ### Create filenames for soon-to-be-created zdrive output files.
  $testname = $file;
  $testname =~ s/zdrive\.inp\.//g;

  $archfilebase = sprintf("%s.%s.%s.", $dirname, $testname, $np);
  $archdropbase = sprintf("%s.%s.drops.%s.", $dirname, $testname, $np);
  if ($debug) {
    print "DEBUG  Test name:  $testname\n";
    print "DEBUG  Archfilebase: $archfilebase; Dropbase: $archdropbase\n";
  }

  ### For serial tests only...if answer file doesn't exist, skip the test.
  if ($np == 1) {
    $zoutfile = sprintf("%s%d", $archfilebase, 0);
    if (!(-e "answers/$zoutfile")) {
      print LOG "Test $dirname:$testname SKIPPED (no answer file)\n";
      print "Test $dirname:$testname SKIPPED (no answer file)\n";
      next TEST;
    }
  }

  ### Execute zdrive.exe.
  $zouterrfile = sprintf("%s.%s.%s.outerr", $dirname, $testname, $np);
  if ($np > 1) {
    $cmd = sprintf("$mpiexec $mpiexecargs -np %d %s %s 2>&1 | tee %s\n", 
                    $np, $zdrive, $file, $zouterrfile);
  }
  else {
    $cmd = sprintf("%s %s 2>&1 | tee %s\n", $zdrive, $file, $zouterrfile);
  }
  if ($debug) {print "DEBUG Executing now:  $cmd\n";}
  $result = system($cmd);
  if ($debug) {print "DEBUG system results $result\n";}

  ### Copy zdrive output files to output directory.
  $failed = 0;
  move($zouterrfile, "output/$zouterrfile");
  foreach $ii (0..$np-1) {
    if ($np < 10) {$format = "%s%d";}
    else {$format = "%s%02d";}
    $zoutfile = sprintf("$format", $zoutfilebase, $ii);
    $zoutdrop = sprintf("$format", $zoutdropbase, $ii);
    $archfile = sprintf("output/$format", $archfilebase, $ii);
    $archdrop = sprintf("output/$format", $archdropbase, $ii);
    $answfile = sprintf("answers/$format", $archfilebase, $ii);
    $answdrop = sprintf("answers/$format", $archdropbase, $ii);
    if ($debug) {print "DEBUG moving files:  $zoutfile $archfile\n";}
    if (-e "$zoutfile") {
      move($zoutfile, $archfile);

      ### Diff the zdrive output files with the accepted answer.
      ### File comparison, ignoring whitespace.
      if ($debug) {print "DEBUG comparing files:  $answfile $archfile\n";}
      $result = compare($archfile,$answfile,sub{nowhite($_[0]) ne nowhite($_[1])});
      if ($result != 0) {$failed++;}
    }
    else {
      ### Failure if no output files.
      $failed = -1;
      last;
    }  
    if (!(-e "$answfile")) {
      ### Failure if no answer files.
      $failed = -2;
      last;
    }

    ### Diff the drop test output files (if any) with the accepted answer.
    ### File comparison, ignoring whitespace.
    if (-e "$zoutdrop") {
      move($zoutdrop, $archdrop);
      if ($debug) {print "DEBUG comparing files:  $answdrop $archdrop\n";}
      $result = compare($archdrop,$answdrop,sub{nowhite($_[0]) ne nowhite($_[1])});
      if ($result != 0) {$failed = 1;}
    }
    if ($debug) {print "DEBUG COMPARISON $result   $failed\n";}
  }

  if ($failed ne 0) {
    $reason = "(Diff failed on $failed files)";
    if ($failed == -1) {
      $reason = "(Missing output files)";
    }
    if ($failed == -2) {
      $reason = "(Missing answer files)";
    }
    print LOG "Test $dirname:$testname FAILED $reason\n";
    print "Test $dirname:$testname FAILED $reason\n";
    $failcnt++;
  }
  else {
    print LOG "Test $dirname:$testname PASSED\n";
    print "Test $dirname:$testname PASSED\n";
    $passcnt++;
  }

  $testcnt++;
}

print LOG "Test $dirname:  $passcnt out of $testcnt tests PASSED.\n";
print "Test $dirname:  $passcnt out of $testcnt tests PASSED.\n";
if ($failcnt > 0) {
  print LOG "Test $dirname:  $failcnt out of $testcnt tests FAILED.\n";
  print "Test $dirname:  $failcnt out of $testcnt tests FAILED.\n";
}
$time = localtime;
print LOG "$time\n";
close(LOG);
