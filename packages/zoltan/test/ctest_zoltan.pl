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

### Check command line arguments.
### Usage:  ctest_zoltan.pl #processors.
foreach $argnum (0 .. $#ARGV) {
   print "$ARGV[$argnum]";
}
print "\n";

$numArgs = $#ARGV + 1;
if ($numArgs < 1) {
  print "Usage:  ctest_zoltan.pl #processors [debug]\n";
  exit -1;
}

### Get number of processors
$np = $ARGV[0];

### Get debug level for ctest_zoltan.pl.
$debug = 0;
if ($numArgs > 1) {$debug = $ARGV[1];}

### Assign the executable.
$zdrive = "../../src/driver/zdrive";

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
$zoutlogfile = sprintf("%s.logfile", $dirname);
open(LOG, "> $zoutlogfile");

### Get list of input files
@inpfiles = glob("zdrive.inp.*");

### Set test counters
$testcnt = 0;
$failcnt = 0;
$passcnt = 0;

### For each zdrive.inp file, run the test.
foreach $file (@inpfiles) {
  if ($debug) {print "DEBUG  Running test $testcnt on $file\n";}

  ### Remove zdrive output files from previous runs.
  system("/bin/rm -f $zoutfilebase*");
  system("/bin/rm -f $zoutdropbase*");

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
  ### Perhaps should collect the result from zdrive.
  system(@cmd);

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
    copy($zoutfile, $archfile);
    copy($zoutdrop, $archdrop);

    ### Diff the zdrive output files with the accepted answer.
    ### File comparison, ignoring whitespace.
    if ($debug) {print "DEBUG comparing files:  $answfile $archfile\n";}
    $result = compare($archfile,$answfile,sub{nowhite($_[0]) ne nowhite($_[1])});
    if ($result != 0) {
      $failed = 1;
    }
    if ($debug) {print "DEBUG COMPARISON $result   $failed\n";}

    ### Need to add drop-test comparison here, too.
  }
  if ($failed) {
    print LOG "Test $dirname:$testname FAILED\n";
    $failcnt++;
  }
  else {
    print LOG "Test $dirname:$testname PASSED\n";
    $passcnt++;
  }
  $testcnt++;
}

print LOG "Test $dirname:  $passcnt out of $testcnt tests PASSED.\n";
print LOG "Test $dirname:  $failcnt out of $testcnt tests FAILED.\n";
close(LOG);
