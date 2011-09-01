#!/usr/bin/perl

use File::Copy;
use File::Compare;
use Getopt::Long;

##############################################################################
##############################################################################
### Remove white-space from a file.
sub nowhite($) {
  my $file = $_[0];
  for ($file) {
    s/ //g;
    s/\t//g;
    s/\r//g;
  }
  return uc($file);
}

##############################################################################
##############################################################################

### Text string needed by CTEST to provide unlimited output.
### Otherwise, output is limited to 50K characters.
print "CTEST_FULL_OUTPUT\n";

### Check command line arguments.
### Usage:  ctest_zoltan.pl --np #processors --pkg package [--debug 0/1] [--mpiexec mpiexec_path] [--mpiexecarg mpiexec_arguments] [--zdrive zdrive_path]
foreach $argnum (0 .. $#ARGV) {
   print "$ARGV[$argnum]";
}
print "\n";

$numArgs = $#ARGV + 1;
if ($numArgs < 2) {
  print "Usage:  ctest_zoltan.pl --np #processors --pkg package [--debug 0/1] [--mpiexec mpiexec_path] [--mpiexecarg mpiexec_arguments] [--zdrive zdrive_path]\n";
  exit -1;
}

### Set default parameter values.
$debug = 0;
$np = -1;
$package = "zzzzz";
$mpiexec = "mpiexec";
$mpiexecarg = "-np";

### Assign the executable.
$zdrive = "../zdrive.exe";
unless (-x $zdrive) {
# Maybe the copy didn't work; look harder.
  $zdrive = "../../src/driver/zdrive.exe";
}
unless (-x $zdrive) {
# Maybe you're on a Windows platform; look harder still.
  $zdrive = "../../src/driver/$ENV{TEST_CONFIG}/zdrive.exe";
}

### Get user-supplied parameter values.
GetOptions('np=i' => \$np,
           'pkg=s' => \$package,
           'debug' => \$debug,
           'mpiexec=s' => \$mpiexec,
           'mpiexecarg=s' => \$mpiexecarg,
           'zdrive=s' => \$zdrive);

if ($debug) {print "DEBUG:  package $package\n";}

### Error test the inputs
die "zdrive executable not found\n" unless (-x $zdrive);
die "number of processors (--np) must be specified\n" unless ($np > 0);
die "package (--pkg) must be specified\n" unless ($package ne "zzzzz");

### Test if MPI implementation supports --mca option
$mpiexecextraargs = "--mca mpi_yield_when_idle 1";
$result = system("$mpiexec $mpiexecextraargs $mpiexecarg 1 uptime");
if ($result) {
  $mpiexecextraargs = "";
}
if ($debug) {print "DEBUG:  mpiexec $mpiexec $mpiexecextraargs $mpiexecarg\n";}


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
                    zdrive.inp.phg* 
                    zdrive.inp.block*
                    zdrive.inp.color*
                    zdrive.inp.cyclic*
                    zdrive.inp.random* 
                    zdrive.inp.graph-repartit*
                    zdrive.inp.graph-partit*");
} 
if ($package eq "ParMETIS") {
  ### ParMETIS tests
  @oneprocfiles = glob("zdrive.inp.*order-metis-*");
  if ($np > 1) {
    @inpfiles = glob("zdrive.inp.adap*
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
  @oneprocfiles = glob("zdrive.inp.*order-scotch*");
} 
if ($package eq "PaToH") {
  ### PaToH tests
  @inpfiles = glob("zdrive.inp.*patoh*");
}

### Set test counters
$testcnt = 0;
$failcnt = 0;
$passcnt = 0;
$purifyerrcnt = 0;
$zmemerrcnt = 0;

### For each zdrive.inp file, run the test.
TEST:  foreach $file (@inpfiles) {
  ### By default run on $np cpus
  $loop_np = $np;

  if ($debug) {print "DEBUG  Running test $testcnt on $file\n";}

  ### Create filenames for soon-to-be-created zdrive output files.
  $testname = $file;
  $testname =~ s/zdrive\.inp\.//g;

  ### Skip serial ordering tests run on more than 1 cpu
  if ($loop_np > 1) {
    if (grep /$file/, @oneprocfiles) {
      $loop_np = 1;
      print LOG "Test $dirname:$testname run on 1 cpu as serial only\n";
      print "Test $dirname:$testname run on 1 cpu as serial only\n";
    }
  }
  ### Define template names here as number of processors may change.
  $zoutfilebase = sprintf("%s.out.%d.", $dirname, $loop_np);
  $zoutdropbase = sprintf("%s.drops.%d.", $dirname, $loop_np);

  ### Remove zdrive output files from previous runs.
  unlink glob("$zoutfilebase*");
  unlink glob("$zoutdropbase*");

  $archfilebase = sprintf("%s.%s.%s.", $dirname, $testname, $loop_np);
  $archdropbase = sprintf("%s.%s.drops.%s.", $dirname, $testname, $loop_np);
  if ($debug) {
    print "DEBUG  Test name:  $testname\n";
    print "DEBUG  Archfilebase: $archfilebase; Dropbase: $archdropbase\n";
  }

  ### For serial tests only...if answer file doesn't exist, skip the test.
  if ($loop_np == 1) {
    $zoutfile = sprintf("%s%d", $archfilebase, 0);
    if (!(-e "answers/$zoutfile")) {
      print LOG "Test $dirname:$testname SKIPPED (no answer file)\n";
      print "Test $dirname:$testname SKIPPED (no answer file)\n";
      next TEST;
    }
  }

  ### Execute zdrive.exe.
  $zouterrfile = sprintf("%s.%s.%s.outerr", $dirname, $testname, $loop_np);
  if ($np > 1) {  # Test on $np because we want to know is the binary needs mpiexec
    $cmd = sprintf("$mpiexec $mpiexecextraargs $mpiexecarg %d %s %s 2>&1 | tee %s\n", 
                    $loop_np, $zdrive, $file, $zouterrfile);
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
  foreach $ii (0..$loop_np-1) {
    if ($loop_np < 10) {$format = "%s%d";}
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
  ### look for purify errors in outerr file.
  ### look for Zoltan memory tool errors, too.
  ### OK if didn't run under purify; it just won't find any purify errors.
  open (ZOUTERR, "output/$zouterrfile") || print "cannot open output/$zouterrfile\n";
  while ($text = <ZOUTERR>) {
    chomp($text);
    @foo = grep(/Possible memory error/, $text);
    $zmemerrcnt += @foo;
    @foo = grep(/ABR:|ABW:|BRK:|BSR:|BSW:|COR:|FMM:|FMR:|FMW:|FNH:|FUM:|IPR:|IPW:|MAF:|MIU:|MLK:|MRE:|MSE:|NPR:|NPW:|PLK:|SBR:|SBW:|SIG:|SOF:|UMC:|UMR:|WPF:|WPM:|WPN:|WPR:|WPW:|WPX:|ZPR:|ZPW:/, $text);
    $nerr = @foo;
    if ($nerr) {
      # Check for errors from the third-party libraries.
      # This check is really a hack, as it is very specific to certain errors
      # and makes assumptions about what the output file looks like.
      # A more general solution would be better.
      $nparmerr = 0;
      $nscotcherr = 0;
      # Look ahead a few lines in the file to find ParMETIS or Scotch errors.
      foreach $ii (1..4) {
        $text = <ZOUTERR>;
        chomp($text);
        # ParMETIS has MLK in CheckInputs; we'll ignore it.
        @parmerr = grep(/CheckInputs/, $text);
        $nparmerr += @parmerr;
        # Scotch has UMR in _SCOTCHhdgraphGather, called by hdgraphOrderNdFold2;
        # we'll ignore it
        @scotcherr = grep(/hdgraphOrderNdFold2/, $text);
        $nscotcherr += @scotcherr;
      }
      if (!$nparmerr && !$nscotcherr) {
        # There are purify errors that don't come from ParMETIS or Scotch.
        print LOG "     Purify error in $zouterrfile:  $_\n" foreach @foo;
        print "     Purify error in $zouterrfile:  $_\n" foreach @foo;
        $purifyerrcnt += $nerr;
      }
    }
  }
  close(ZOUTERR);

  $testcnt++;
}

print LOG "Test $dirname:  $passcnt out of $testcnt tests PASSED.\n";
print "Test $dirname:  $passcnt out of $testcnt tests PASSED.\n";
if ($failcnt > 0) {
  print LOG "Test $dirname:  $failcnt out of $testcnt tests FAILED.\n";
  print "Test $dirname:  $failcnt out of $testcnt tests FAILED.\n";
}
if ($purifyerrcnt) {
  print LOG "Test $dirname:  $purifyerrcnt purify errors; test FAILED.\n";
  print "Test $dirname:  $purifyerrcnt purify errors; test FAILED.\n";
}
if ($zmemerrcnt) {
  print LOG "Test $dirname:  $zmemerrcnt Zoltan memory errors; test FAILED.\n";
  print "Test $dirname:  $zmemerrcnt Zoltan memory errors; test FAILED.\n";
}
$time = localtime;
print LOG "$time\n";
close(LOG);
