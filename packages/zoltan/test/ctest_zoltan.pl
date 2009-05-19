#!/usr/bin/perl

use File::Copy;
use File::Compare;

##############################################################################
sub print_time {
    my $fp = shift;
    my @abbr = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $year += 1900;
    $tt = sprintf("%02d:%02d:%02d", $hour, $min, $sec);
    print $fp $mday, $abbr[$mon], $year, " ", $tt, "\n";
}

##############################################################################
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

### Check command line arguments.
### Usage:  ctest_zoltan.pl #processors [debug].
foreach $argnum (0 .. $#ARGV) {
   print "$ARGV[$argnum]";
}
print "\n";

$numArgs = $#ARGV + 1;
if ($numArgs < 2) {
  print "Usage:  ctest_zoltan.pl #processors package [debug]\n";
  exit -1;
}

### Get debug level for ctest_zoltan.pl.
$debug = 0;
if ($numArgs > 1) {$debug = $ARGV[2];}

### Get number of processors
$np = $ARGV[0];

### Get the package number indicating which tests to run.
$package = $ARGV[1];
if ($debug) {print "DEBUG:  package $package\n";}

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
print "KDDKDD $package\n";
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
                    zdrive.inp.random* 
                    zdrive.inp.graph
                    zdrive.inp.graph-re*
                    zdrive.inp.graph-partition");
} 
if ($package eq "ParMETIS") {
  ### ParMETIS tests
  @inpfiles = glob("zdrive.inp.ada* 
                    zdrive.inp.part* 
                    zdrive.inp.*metis*");
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
foreach $file (@inpfiles) {
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

  ### Execute zdrive.exe.
  $zouterrfile = sprintf("%s.%s.%s.outerr", $dirname, $testname, $np);
  if ($np > 1) {
    $cmd = sprintf("mpiexec -np %d %s %s | tee %s\n", $np, $zdrive, $file, 
                   $zouterrfile);
  }
  else {
    $cmd = ("%s %s | tee %s\n", $zdrive, $file, $zouterrfile);
  }
  if ($debug) {print "DEBUG Executing now:  $cmd\n";}
  $result = system($cmd);
  if ($debug) {print "DEBUG system results $result\n";}

  ### Copy zdrive output files to output directory.
  $failed = 0;
  copy($zouterrfile, "output/$zouterrfile");
  foreach $ii (0..$np-1) {
    $zoutfile = sprintf("%s%d", $zoutfilebase, $ii);
    $zoutdrop = sprintf("%s%d", $zoutdropbase, $ii);
    $archfile = sprintf("output/%s%d", $archfilebase, $ii);
    $archdrop = sprintf("output/%s%d", $archdropbase, $ii);
    $answfile = sprintf("answers/%s%d", $archfilebase, $ii);
    $answdrop = sprintf("answers/%s%d", $archdropbase, $ii);
    if ($debug) {print "DEBUG copying files:  $zoutfile $archfile\n";}
    if (-e "$zoutfile") {
      copy($zoutfile, $archfile);

      ### Diff the zdrive output files with the accepted answer.
      ### File comparison, ignoring whitespace.
      if ($debug) {print "DEBUG comparing files:  $answfile $archfile\n";}
      $result = compare($archfile,$answfile,sub{nowhite($_[0]) ne nowhite($_[1])});
      if ($result != 0) {$failed = 1;}
    }
    else {
      ### Failure if no output files.
      $failed = 1;
    }  

    ### Diff the drop test output files (if any) with the accepted answer.
    ### File comparison, ignoring whitespace.
    if (-e "$zoutdrop") {
      copy($zoutdrop, $archdrop);
      if ($debug) {print "DEBUG comparing files:  $answdrop $archdrop\n";}
      $result = compare($archdrop,$answdrop,sub{nowhite($_[0]) ne nowhite($_[1])});
      if ($result != 0) {$failed = 1;}
    }
    if ($debug) {print "DEBUG COMPARISON $result   $failed\n";}
  }

  if ($failed) {
    print LOG "Test $dirname:$testname FAILED\n";
    print "Test $dirname:$testname FAILED\n";
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
