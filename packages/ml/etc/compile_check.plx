#!/usr/bin/perl -w
#
# compile_check - a dumb script to try out different configure options by 
#                 continuously compiling with a random set of things enabled. 
#                 To add more configuration options, just keep appending new
#                 strings to the variable @config_choices.
#
# Ray Tuminaro, 1/1/05
#
@config_choices = ("--enable-epetra", "--enable-anasazi", "--enable-aztecoo","--enable-ifpack", "--with-ml_superlu",
            "--enable-amesos", "--enable-ifpack-metis", "--with-ml_metis", "--with-ml_parmetis3x","--with-ml_zoltan",
            "--enable-teuchos","--with-ml_arpack","--enable-triutils"); 
$config_always= "--disable-default-packages --disable-examples --disable-tests --enable-ml-examples --enable-ml-tests --enable-ml\n";
$lib_success_file="packages/ml/src/libml.a";
$TRI_HOME="/home/tuminaro/Trilinos";
$SCRIPT_HOME= $TRI_HOME . "/packages/ml/test/testharness/marzio-laptop";
$TESTHARNESS_INPUTFILE="compile-check-input";
#
# Shouldn't have to change anything below this line!!!!
#
#
sub trimwhitespace($);
$TESTHARNESS_INPUTFILE=trimwhitespace($TESTHARNESS_INPUTFILE);
$SCRIPT_HOME          =trimwhitespace($SCRIPT_HOME);
$TRI_HOME             =trimwhitespace($TRI_HOME);
#
#
#  Parse command line options
#
#    Usage: compile_check.plx [--trihome TrilinosHome] [--scripthome ScriptHome] [--inputfile TestHarnessInputFile]
use Getopt::Long;
GetOptions("trihome=s"=>\$TRI_HOME,
	   "scripthome=s"=>\$SCRIPT_HOME,
	   "help"=>\$HELP_REQUESTED,
	   "inputfile=s"=>\$TESTHARNESS_INPUTFILE);

#          "verbose!"=>\$verboseornoverbose,"optional:s",\$optionalstring,"int=i"=>\$mandatoryinteger,"float=f"=>\$mandatoryfloat,
#
if (@ARGV || $HELP_REQUESTED) {
    if ($HELP_REQUESTED) {
        print  "\n * A dumb script to try out different configure options by continuously compiling with a random set\n";
	print  " * of things enabled. To add more configuration options, just edit this script and append new\n";
        print  " * strings to the variable \@config_choices. The current set of config choices include:\n";
        print  " *     @config_choices\n\n";
        print  " *** The Detailed Description: This script creates the TRILINOS_CONFIG_FILE expect by test_harness.plx.\n";
        print  " ***                           It is ASSUMED that if all of the above options were put into a file\n";
        print  " ***                           called 'compile-check-config', test_harness.plx would properly configure\n";
        print  " ***                           and compile Trilinos. compile_check.plx simply tries different subsets of\n";
        print  " ***                           these choices and checks for success by seeing if the file given by\n";
        print  " ***                           \$lib_success_file (currently = $lib_success_file)\n";
        print  " ***                           exists. If this is not found, it is assumed that something failed and\n";
        print  " ***                           the infinite loop is terminated.\n";
        print  " ***                           Note: \$config_always contains configure options which are always on.\n\n";




    }

    print "\nUsage: compile_check.plx [--trihome TrilinosHome] [--scripthome ScriptHome] [--inputfile TestHarnessInputFile] [--help]\n";
    print "\nwhere\n";
    print "     TrilinosHome:        The main home directory where Trilinos lives.\n\n";
    print "     ScriptHome:          The directory where the inputfiles for the Trilinos script testharness.plx live.\n";
    print "                          These input files are used to guide the configuration process (for more details \n";
    print "                          see http://software.sandia.gov/Trilinos/developer/test_harness.html ).\n";
    print "                          Note: compile_check.plx adds a file to this directory 'compile-check-config'\n\n";
    print "     TestHarnessInputFile:The main input file used by test_harness.plx. This file must contain a line\n";
    print "                          of the form: TRILINOS_CONFIG_FILE  = \$SCRIPT_HOME/compile-check-config\n";
    print "                          where \$SCRIPT_HOME is set by the --scripthome option described above\n";
    print "                          (see  http://software.sandia.gov/Trilinos/developer/test_harness.html).\n";
    exit;
}

#

$FLAG=1;

$Nconfig_choices = @config_choices;

chdir($SCRIPT_HOME);

#
# seed the random number generator
#
srand(time() ^($$ + ($$ <<15))) ;


while($FLAG) {

    # This will generate a random integer between (1 and $Nconfig_choices)
    $Npackages=int($Nconfig_choices*rand() + .9999);
    print "Number of packages = $Npackages \n";

    # This creates the list of configure choices with a random number preceeding
    # each choice. This list will then be sorted and then the random number will
    # be stripped off when writing the choices to a file.
    #
    for ( $i=0; $i < $Nconfig_choices; $i++) {
       $tlist[$i] = rand() . "\#" . $config_choices[$i];
    }
    @sortedlist = sort(@tlist);
    open (CONFIG_FILE, ">compile-check-config")
            or die "can't open compile-check-config";
    print CONFIG_FILE $config_always;
    for ( $i=0; $i < $Npackages; $i++) {
	@temp= split("#",$sortedlist[$i]);
	print CONFIG_FILE "$temp[1] \n";
    }
    close CONFIG_FILE;
    system "cat compile-check-config";
#
#   run the Trilinos perl script
#
    chdir($TRI_HOME . "/testharness");
    system "perl test-harness.plx -f ${SCRIPT_HOME}/$TESTHARNESS_INPUTFILE";
    chdir($SCRIPT_HOME);
    if ( !(-e $TESTHARNESS_INPUTFILE)) { $FLAG=0; }

#
#   check if the library was created in MPI_DIR or SERIAL_DIR (depending on what was
#   requested in the $TESTHARNESS_INPUTFILE. If the library is not there set the flag
#   so that the script stops.
#
    $MPI_DIR=`grep MPI_DIR  $TESTHARNESS_INPUTFILE | sed "s/^.*=//" | sed "s/ //g"` ;
    $SERIAL_DIR=`grep SERIAL_DIR  $TESTHARNESS_INPUTFILE | sed "s/^.*=//" | sed "s/ //g"`;
    chop($MPI_DIR);
    chop($SERIAL_DIR);
    if ( $MPI_DIR ) {
       $LIBFILE=$TRI_HOME . "/" . $MPI_DIR . "/" . $lib_success_file;

       system "ls -l $LIBFILE";
       if (!(-e $LIBFILE)) { $FLAG=0; }
    }
    if ( $SERIAL_DIR && ($FLAG) ) {
       $LIBFILE=$TRI_HOME . "/" . $SERIAL_DIR . "/" . $lib_success_file;
       system "ls -l $LIBFILE";
       if (!(-e $LIBFILE)) { $FLAG=0; }
    }

}

# Remove whitespace from the start and end of the string
sub trimwhitespace($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
