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
print "DEBUG  ";
foreach $argnum (0 .. $#ARGV) {
   print "$ARGV[$argnum]";
}
print "\n";

$numArgs = $#ARGV + 1;
if ($numArgs != 1) {
  print "Usage:  ctest_zoltan.pl #processors\n";
  exit -1;
}

### Get number of processors
$np = $ARGV[0];

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
print "DEBUG  Dir $dir dirname $dirname\n";
$zoutfilebase = sprintf("%s.out.%d.", $dirname, $np);
$zoutdropbase = sprintf("%s.drops.%d.", $dirname, $np);
print "DEBUG  Outfilebase: $zoutfilebase;  Dropbase: $zoutdropbase\n";

### Get list of input files
@inpfiles = glob("zdrive.inp.*");

### Set test counters
$testcnt = 0;
$failcnt = 0;
$passcnt = 0;

### For each zdrive.inp file, run the test.
foreach $file (@inpfiles) {
  print "DEBUG  Running test $testcnt on $file\n";

  # Remove zdrive output files from previous runs.
  system("/bin/rm -f $zoutfilebase*");
  system("/bin/rm -f $zoutdropbase*");

  # Create filenames for soon-to-be-created zdrive output files.
  $testname = $file;
  $testname =~ s/zdrive\.inp\.//g;
  print "DEBUG  Test name:  $testname\n";

  $archfilebase = sprintf("%s.%s.%s.", $dirname, $testname, $np);
  $archdropbase = sprintf("%s.%s.drops.%s.", $dirname, $testname, $np);
  print "DEBUG  Archfilebase: $archfilebase; Dropbase: $archdropbase\n";

  # Execute zdrive.exe.
  if ($np > 1) {
    @cmd = ("mpiexec", "-np", $np, $zdrive, $file);
  }
  else {
    @cmd = ($zdrive, $file);
  }
  print "DEBUG Executing now:  @cmd\n";
  system(@cmd);

  ### Copy zdrive output files to output directory.
  $failed = 0;
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
    print "DEBUG comparing files:  $answfile $archfile\n";
    $result = compare($archfile,$answfile,sub{nowhite($_[0]) ne nowhite($_[1])});
    if ($result != 0) {
      $failed = 1;
    }
    print "DEBUG COMPARISON $result   $failed\n";

    ### Need to add drop-test comparison here, too.
  }
  if ($failed) {
    print "Test $dirname:$testname FAILED\n";
    $failcnt++;
  }
  else {
    print "Test $dirname:$testname PASSED\n";
    $passcnt++;
  }
  $testcnt++;
}

print "Test $dirname:  $passcnt out of $testcnt tests PASSED.\n";
print "Test $dirname:  $failcnt out of $testcnt tests FAILED.\n";
