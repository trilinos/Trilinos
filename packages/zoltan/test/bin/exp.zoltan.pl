#!/usr/bin/perl
#

##
## BEGIN COPY HERE
##
use Cwd;

#
# See if the first argument is '-I' which specifies an include path
# for the TestLib module
#
$ctr=0;
if ($ARGV[0] eq "-I") {
   unshift @INC, $ARGV[1];
   $ENV{'PATH'} = $ENV{'PATH'} . ":.:" . $ARGV[1];
   $ctr = 2;
   }
#
# See if the TESTLIBDIR environmental variable is set, which provides
# an include path for the TestLib module
#
$pwd = `pwd`;
chomp($pwd);
if ($ENV{TESTLIBDIR} ne "") {
   $tmp = $ENV{TESTLIBDIR};
   @words = split /[ ]+/, $tmp, $#tmp;
   foreach $word (@words) {
     unshift @INC, $word;
     }
   }
require TestLib;
##
## END COPY HERE
##

my (%bal, %cutl, %cutn);
my $xml;
my $iters;

#
# Statistics collected by the analyzeData call-back routine
#

sub analyzeData {
  $iters = 0;

  ($bal{MAX},$bal{MIN},$bal{AVG}) = ("ERROR","ERROR","ERROR");
  ($cutl{MAX},$cutl{MIN},$cutl{AVG}) = ("ERROR","ERROR","ERROR");
  ($cutn{MAX},$cutn{MIN},$cutn{AVG}) = ("ERROR","ERROR","ERROR");
  $timeavg = -1;

  my $filename = @_[0];
  open ZOUT, "$filename" or die "can't open zdrive output $filename";

  while (<ZOUT>) {
    if (/.*STATS Runs ([0-9]+)/) {
      my $vname="";
      $iters = "$1";
      if (/ bal (.*)/) {
	%bal = (split ' ', $1);
      } elsif (/ cutl (.*)/) {
	%cutl = (split ' ', $1);
      } elsif (/ cutn (.*)/) {
	%cutn = (split ' ', $1);
      }
    } elsif (/^FILE .*Average:  ([\d.e+]+) seconds per Iteration/) {
      $timeavg = $1 + 0;
    }
  }
  close ZOUT;

#  print "iters: $iters\n";
#  print "bal: $bal{MAX} $bal{MIN} $bal{AVG}\n";
#  print "cutl: $cutl{MAX} $cutl{MIN} $cutl{AVG}\n";
#  print "cutn: $cutn{MAX} $cutn{MIN} $cutn{AVG}\n";
#  print "time: $timeavg\n";

  if ($timeavg == -1) {
     @_[1]  = "Fail";
  }
}

#
# Default description value
#
$description="  <Description>ERROR</Description>\n";

#
# A callback that is used to print measurement information
#

sub printMeasurements {
  print "  <Trial id=\"0\" seed=\"0\" data=\"?\">\n";
  print "     <Value name=\"BalMax\" type=\"numeric/double\">$bal{MAX}</Value>\n";
  print "     <Value name=\"BalMin\" type=\"numeric/double\">$bal{MIN}</Value>\n";
  print "     <Value name=\"BalAvg\" type=\"numeric/double\">$bal{AVG}</Value>\n";
  print "     <Value name=\"CutlMax\" type=\"numeric/double\">$cutl{MAX}</Value>\n";
  print "     <Value name=\"CutlMin\" type=\"numeric/double\">$cutl{MIN}</Value>\n";
  print "     <Value name=\"CutlAvg\" type=\"numeric/double\">$cutl{AVG}</Value>\n";
  print "     <Value name=\"CutnMax\" type=\"numeric/double\">$cutn{MAX}</Value>\n";
  print "     <Value name=\"CutnMin\" type=\"numeric/double\">$cutn{MIN}</Value>\n";
  print "     <Value name=\"CutnAvg\" type=\"numeric/double\">$cutn{AVG}</Value>\n";
  print "     <Value name=\"SolverTime\" type=\"numeric/double\">$timeavg</Value>\n";
  print "     <Value name=\"Iterations\" type=\"numeric/int\">$iters</Value>\n";

  print "  </Trial>\n";

}

##
## MAIN ROUTINE
##
if (! @ARGV) {
   print "\n";
   print "$0 - Zoltan big experiment driver\n";
   print "\n";
   print "usage:\n";
   print "\n";
   print "   $0 [options] <expname> <testname> <factor-id> <factor-value> [...]\n";
   print "\n";
   print "options:\n";
   print "\n";
   print "   --debug\t\tPrint debugging information.\n";
   print "\n";
   print "notes:\n";
   print "\n";
   print "   This script is launched by 'runexp' to execute and process a\n";
   print "   single test from an experiment file.  The output of this test\n";
   print "   is stored in the file <expname>-<testname>.out, which can be\n";
   print "   reviewed after the test is performed.  This script generates an\n";
   print "   XML summary that is included in the <expname>.results.xml file\n";
   print "   by 'runexp'.\n";
   print "\n";
   print "   This script is not intended for interactive use.  Consequently,\n";
   print "   the processing of the command-line options is somewhat fragile.\n";
   print "   For example, it depends on the order in which the options are\n";
   print "   processed in this script!\n";
   print "\n";
   print "   This script assumes that the first factor is the test information,\n";
   print "   and the remaining factors specify command-line options.\n";
   exit;
}
#
# Process command line arguments
#
$debug = 0;
$testoptions = "";
if ($ARGV[$ctr] eq "--debug") {
   $ctr += 1;
   $debug = 1;
}

if ($debug == 1) {
   print "<debuginfo>\n"
}

$description = TestLib::generate_description($ctr,@ARGV);
$expname = $ARGV[$ctr];
$ctr +=1;
$testname = $ARGV[$ctr];
$ctr +=1;
$ntrials = $ARGV[$ctr];
$ctr+=3;
while ($ctr <= $#ARGV) {
  $testoptions .= $ARGV[$ctr] . " ";
  $ctr+=3;
}


#
# Process experiment tag/value pairs
#
# This is a convention that is not explicitly supported by the 
# 'runexp' utility, but this is a convenient way for
# structuring experiment testing.
#
$solnvalue = "unknown";
$tolerance = 1e-12;
$files = "";
%testparams = ();
$auxoptions = TestLib::process_testoptions($testoptions, \%testparams);

if (defined($testparams{"_optimum"})) {
   $solnvalue = $testparams{"_optimum"};
} elsif (defined($testparams{"_value"})) {
   $solnvalue = $testparams{"_value"};
}

$iters = 10; # default
foreach $word (TestLib::get_tokens($testparams{"_iterations"})) {
  $iters = $word;
}

foreach $word (TestLib::get_tokens($testparams{"_inp"})) {
  $files = $files . $word . " ";
}

# can configure a short script to create a custom zdrive here, if desired
# but it uses the commandline to specify the name
my $solver = "/Net/local/mpi/build/solaris/ch_p4/bin/mpirun -np $testparams{_np} ../../Obj_solaris/zdrive";

# $solver = "/bin/cat test.out"; #debug

if ($debug == 1) {
   print "options: $testoptions\n";
   print "exp:       $expname\n";
   print "value:     $solnvalue\n";
   print "input:     $files\n";
   print "numproc:   $testparams{_np}\n";
   print "iters:     $testparams{_iterations}\n";
}
   #
   # Setup command line
   #
   $cmdline = "$solver ";
   @options = split(/ /,$auxoptions,$#auxoptions);
   while (@options) {
     local $tmp = shift @options;
     if ($tmp ne "") {
        $cmdline .= " --" . $tmp;
     }
   }
   $cmdline .= " " . $files;
   print "CMDLINE: $cmdline\n" if ($debug == 1);

# open Z
# open ZINP ">zdrive.inp"

if ($debug == 1) {
   print "</debuginfo>\n"
}

#
# Launch the test with the TestLib driver.
#
TestLib::run_experiments(\*STDOUT, $expname, $testname, $cmdline, $ntrials, \&analyzeData, \&printMeasurements, $auxoptions, \%testparams);
