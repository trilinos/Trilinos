#!/usr/bin/perl -w
#
# compile_check - a dumb script to try out different configure options by 
#                 continuously compiling with a random set of things enabled. 
#                 To add more configuration options, just keep appending new
#                 strings to the variable @config_choices.
#
# Ray Tuminaro, 1/1/05
#

#
# Options that we want to randomly turn on and off
#
@config_choices = ("--enable-epetra", "--enable-anasazi", "--enable-aztecoo",
                   "--enable-ifpack", "--with-ml_superlu",
                   "--enable-amesos", "--enable-ifpack-metis", 
                   "--with-ml_metis", "--with-ml_parmetis3x",
                   "--with-ml_zoltan","--enable-teuchos",
                   "--with-ml_arpack","--enable-triutils"); 
#
# Options that we always want on 
#
$config_always= "--disable-default-packages --disable-examples --disable-tests --enable-ml-examples --enable-ml-tests --enable-ml\n";
#
# Library file that we want to check to make sure that it was in fact
# created and the build was successful.
#
$lib_success_file="packages/ml/src/libml.a";
#
# Location of Trilinos
#
$TRI_HOME="/home/tuminaro/Trilinos";
#
# Directory where the machine file and basic input file exist
# for the Trilinos test-harness.plx script.
#
$SCRIPT_HOME= $TRI_HOME . "/packages/ml/test/testharness/marzio-laptop";
#
# Name of the input data given to test-harness.plx
#
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
       print "\n * A dumb script to try out different configure options by\n";
       print " * repeatedly compiling with a random set of things enabled.\n";
       print " * To add more configuration options, just edit this script\n"; 
       print " * and append strings to the variable \@config_choices. The\n";
       print " * current set of config choices include:\n";
       print " *     @config_choices\n\n";
       print " *** The Detailed Description: This script creates the";
       print " TRILINOS_CONFIG_FILE\n";
       print " ***                          ";
       print " expect by test_harness.plx. It is ASSUMED that if\n";
       print " ***                          ";
       print " all of the above options were put into a file\n";
       print " ***                          ";
       print " called 'compile-check-config', test_harness.plx\n";
       print " ***                          ";
       print " would properly configure and compile Trilinos.\n";
       print " ***                          ";
       print " compile_check.plx simply tries different subsets\n";
       print " ***                          ";
       print " of these choices and checks for success by seeing\n";
       print " ***                          ";
       print " if the file given by \$lib_success_file\n";
       print " ***                          ";
       print " (currently = $lib_success_file) exits. If\n";
       print " ***                          ";
       print " this is not found, it is assumed that\n";
       print " ***                          ";
       print " something failed and the infinite loop is\n";
       print " ***                          ";
       print " terminated.\n";
       print " ***                          ";
       print " Note: \$config_always contains configure\n";
       print " ***                          ";
       print " options which are always on.\n\n";




    }

    print "\nUsage: compile_check.plx [--trihome TrilinosHome] [--scripthome ScriptHome] [--inputfile TestHarnessInputFile] [--help]\n";
    print "\nwhere\n";
    print "     TrilinosHome:        The main home directory where Trilinos lives.\n\n";
    print "     ScriptHome:          The directory where the inputfiles for the\n";
    print "                          Trilinos script testharness.plx live. These\n";
    print "                          input files are used to guide the configuration\n";
    print "                          process (for more details see \n";
    print "                          http://software.sandia.gov/Trilinos/developer/test_harness.html).\n";
    print "                          Note: compile_check.plx adds a file to this\n";
    print "                          directory called 'compile-check-config'\n\n";
    print "     TestHarnessInputFile:The main input file used by test_harness.plx. This\n";
    print "                          file must contain a line of the form:\n";
    print "                            TRILINOS_CONFIG_FILE  = \$SCRIPT_HOME/compile-check-config\n";
    print "                          where \$SCRIPT_HOME is set by the --scripthome\n";
    print "                          option described above (see \n";
    print "                          http://software.sandia.gov/Trilinos/developer/test_harness.html).\n";
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

    # This creates the list of configure choices with a random number
    # preceeding each choice. This list will then be sorted and then the random
    # number will be stripped off when writing the choices to a file.
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
    system "cat compile-check-config";  # Used for printing, can be removed
#
#   run the Trilinos perl script
#
    chdir($TRI_HOME . "/testharness");
    system "perl test-harness.plx -f ${SCRIPT_HOME}/$TESTHARNESS_INPUTFILE";
    chdir($SCRIPT_HOME);
    if ( !(-e $TESTHARNESS_INPUTFILE)) { $FLAG=0; }

#
#   check if the library was created in MPI_DIR or SERIAL_DIR (depending on 
#   what was requested in the $TESTHARNESS_INPUTFILE. If the library is not
#   there set the flag so that the script stops.
#
#    $MPI_DIR=`grep MPI_DIR  $TESTHARNESS_INPUTFILE | sed "s/^.*=//" | sed "s/ //g"` ;
#    $SERIAL_DIR=`grep SERIAL_DIR  $TESTHARNESS_INPUTFILE | sed "s/^.*=//" | sed "s/ //g"`;
#    chop($MPI_DIR);
#    chop($SERIAL_DIR);
#
# More portable version of the above
#
    $MPI_DIR    = myfirstgrep($TESTHARNESS_INPUTFILE,"MPI_DIR");
    $MPI_DIR    =~ s/^.+\s+=\s+//;
    $SERIAL_DIR = myfirstgrep($TESTHARNESS_INPUTFILE,"SERIAL_DIR");
    $SERIAL_DIR =~ s/^.+\s+=\s+//;
    $MPI_DIR    =trimwhitespace($MPI_DIR);
    $SERIAL_DIR =trimwhitespace($SERIAL_DIR);

    if ( $MPI_DIR ) {
       $LIBFILE=$TRI_HOME . "/" . $MPI_DIR . "/" . $lib_success_file;

       system "ls -l $LIBFILE"; # Only used for printing, can be removed
       if (!(-e $LIBFILE)) { $FLAG=0; }
    }
    if ( $SERIAL_DIR && ($FLAG) ) {
       $LIBFILE=$TRI_HOME . "/" . $SERIAL_DIR . "/" . $lib_success_file;
       system "ls -l $LIBFILE"; # Only used for printing, can be removed
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

#
# special grep that first strips out
# any characters after a '#' and then
# does a grep. Additionally, this function
# only returns the first line that matches
# the grep
#
sub myfirstgrep {
  my $filename = $_[0];
  my $pattern  = $_[1];

  open (THEFILE, "<$filename") or die "can't open $filename";

  while (<THEFILE>) {
      $theline = $_;
      $theline =~ s/#.+//;
      $_ = $theline;
      if (  (m/=/i ) && (m/$pattern/i ))  {
	  close THEFILE;
	  return $theline;
      }
  }
  close THEFILE;
  return "";
}


#
# To set up compile_check on a new machine:
#   1) make sure that any old executables are wiped out of
#      a previously configured directory structure:
#          cd Trilinos
#          find . -name '*.o' -exec rm {} \;
#          find . -name '*.a' -exec rm {} \;
#          find . -name 'config.log' -exec rm {} \;
#          find . -name 'config.status' -exec rm {} \;
#          find . -name 'Makefile.export' -exec rm {} \;
#          find . -name 'Makefile' -exec rm {} \;
#  2) cd ml/etc; compile_check.plx --help
#  3) take the list of packages and put them into a file called
#     compile-check-config
#  4) you may need to build some of third party libraries
#     manually.
#  5) edit the machine and machine-mpi directories to see if they
#     are ok.
#  6) cd testharness and do: 
#          perl test-harness.plx -f /home/tuminaro/Trilinos/packages/ml/test/testharness/marzio-laptop/compile-check-input
#     Note: if things don't work, you might try going to the build
#           directory and running 'invoke-configure' manually.
#           Then, check the config.log files of the package
#           that is bombing.
#
#  Here is some fun stuff to look at the warnings generated
#  by compiling.
#
#  cd ml/src; rm -f *.o warnings; script warnings
#  make; exit
#  vi warnings
#  :g/^if mpic/d 
#  :g/^then mv /d 
#  :/^rm -f libml.a/,$d 
#  :1,/^make\[1\]: Entering directory/d 

