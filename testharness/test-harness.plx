#!/usr/bin/perl -w
# /Trilinos/testharness/test-harness.plx

################################################################################
# The Trilinos Project - Test Harness
#
# Jim Willenbring, Mike Phenow, Ken Stanley
#
################################################################################

use strict;

use Getopt::Std;    # for command line options parsing
use lib "lib";      # ended up needing to use a relative address - don't need?

# Variable Declarations ========================================================

# Options variables
my %flags;                  # command line flags (boolean or string)
my %options;                # config-file hash of arrays (keys, lists of values)

# Host variables
my $hostOS;                 # host operating system
my $hostName;               # host name
my $hostArch;               # host architecture

# Error variables
my $errorCount;             # tally of error log appends
my $warningCount;           # tally of warning log appends
my $updateError;            # update error boolean
my $trilinosCompileError;   # trilinos compile error boolean
my $testCompileError;       # test compile error boolean

# File variables        
# ERRORS;                   # error log         # localize? pass to sendMail()?
# WARNINGS;                 # warning log
# BODY;                     # email body

# Compile variables   
my $comm;   # unclear use...
my $mpi;    # unclear use...

# Test variables
my @potentialScripts;   # used by test and then sendMail
my $potentialScript;    # can these be made local?
        
################################################################################
# Execution ####################################################################
################################################################################

# Preparation ==================================================================

# Get command line options
getopts("gshf:", \%flags);
if ($flags{h}) { printHelp(); }
if ($flags{g}) { genConfigTemp($flags{s}); }
if ($flags{s} && !$flags{g}) {
    print "Error: -s option must be accompanied by the -g option.\n";
    exit;
}
if (!$flags{f}) { 
    print "Error: must specify a config file using the -f option.\n"; 
    print "Run test-harness.plx with the -h option for a complete list of options.\n"; 
    exit;
} 

# two-tiered config system...?
# one and only one main test-harness config (excluded packages, default to and from emails, ...)
# any number of possible config files with things more likely to change

# Trilinos/testharness directory structure...
#
# test-harness(.plx)
# dependencies?
# README
#
# config/
# elements-machine/
# elements-trilinos/

# Parse test-harness-config
# >>> # for 'checkin'-test-harness functionality, maybe check in
# >>> # a standard config file to use for that purpose?
# >>> # what were all the differences between this and the 'checkin' version?
%options = parseConfig($flags{f});

# Determine host OS
# >>> # should this be a config-file option?
# >>> # is it only used for sendMail()?
# >>> # if so, do we care what about the OS, or just the mail program?
# >>> # should the additional libraries be supplied in the config?
getHostOs();

# Update Trilinos from CVS Repository
cvsUpdate();

# Prepare output files
# >>> # are there more errors/warnings we could track?
# >>> # or is there a way to track them in more detail?
# >>> # do we need to use files, or could this be stored in memory?
# >>> # if they do have to be stored in files, is there a universal place to put them?
prepareFiles();

# Boot LAM if any tests are parallel
#   (If mpich is used, comment out this section and the lamhalt section)
# >>> # should we allow for other MPI implementations to be specified in the config file?
# >>> # this way, we wouldn't have to hard-code the implementation at all  
lamboot();

# Main Execution ===============================================================

# Build and test from bottom up, recording and reporting results as necessary
# >>> # how can we have this logic driven by the config?
# >>> # have sets of packages in a specified order...?
# >>> # maintain a dependency table?
run();          # core test-harness logic

# >>> # should each of these run without global variables...?
# >>> # pass directory, etc...
# >>> # return (pass/fail, output, etc...)
# >>> # would there be much benefit to this?

# compile();    # called by run()
# build();      # called by run()
# test();       # called by run()
# sendMail();   # called by run()


# Clean Up =====================================================================

# Halt LAM if any tests are parallel
#   (If mpich is used, comment out this section and the lamboot section)
lamhalt();

# Close filehandles and delete temp files
cleanupFiles();

################################################################################
# Subroutines ##################################################################
################################################################################

    ############################################################################
    # getHostOs()
    #
    # Discovers the host's operating system.
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: 
    #   - returns: 

    sub getHostOs {
        chomp ($hostOS=`uname -a`);
        chomp ($hostName=`uname -n`);
        chomp ($hostArch=`uname`); 
# ???   # $hostArch must match the middle portion of the name of the log file 
# ???   # written to by the testAll file - for example on a Linux platform, 
# ???   # the log file is called logLinux.txt, so $hostArch=Linux
            
        # List of current supported platforms: 
        #   AIX, CPLANT, DEC, IBMSP, LINUX, SGI32, SGI64, SMOS, SOLARIS, TFLOP, sun
# ???   # If MIME::Lite is used in the end, include the file somewhere in Trilinos
        SWITCH: {
            if ($hostOS =~ m/^Linux/) {
                use lib "lib"; # ended up needing to use a relative address
                last SWITCH; 
            };
            if ($hostOS =~ m/^SunOS/) {
                use lib "/local/homes/jmwille/lib";
                last SWITCH; 
            };
            if ($hostOS =~ m/^CYGWIN/) {                                
                # Cygwin doesn't have sendmail, so configure MIME::Lite to use smtp instead
                MIME::Lite->send('smtp',"mailgate.sandia.gov");    
                                             
# >>> # export to config?
                use lib "/home/jmwille/lib";
                last SWITCH; 
            };
            # Fix the rest of the switch statement, currently functional for LINUX only
            # if ($hostOS=~/^Linux.+cluster/||$hostOS=~/^Linux.+node/) 
            #   {$TRILINOS_hostArch='linux-lam-cluster'; last SWITCH; };
            # if ($hostOS=~/^Linux.+cluster/) {$TRILINOS_hostArch='linux-lam-cluster'; last SWITCH; };
            die "hostOS does not match any of the OS choices.\n";
        }
    } # getHostOs()

    ############################################################################
    # cvsUpdate()
    #
    # Updates Trilinos from the CVS repository.
    #   - global variables used: yes
    #   - sends mail: on error
    #   - args: 
    #   - returns: 

    sub cvsUpdate {    
        chdir "$options{'HOMEDIR'}[0]";
        
        my $result;
        print "$options{'CVS_CMD'}[0] update -dP > $options{'HOMEDIR'}[0]/updatelog.txt 2>&1\n";
        $result=system "$options{'CVS_CMD'}[0] update -dP > $options{'HOMEDIR'}[0]/updatelog.txt 2>&1";
        if ($result) {
            $errorCount++;
            $updateError++;
            print ERRORS "*** Error updating Trilinos ***\n";
            print ERRORS "Aborting tests.\n";
            print ERRORS "See updatelog.txt for further information";
            
# >>> #     # sendMail(); # COMMENTING OUT EMAIL
            die " *** ERROR updating TRILINOS! ***\n";
        }
    } # cvsUpdate()

    ############################################################################
    # prepareFiles()
    #
    # Prepares the output files.
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: 
    #   - returns: 

    sub prepareFiles {
        chdir "$options{'HOMEDIR'}[0]";
    
        # Doesn't count the actual number of errors, but rather
        # the number of times ERRORS is written to.
        $errorCount=0; 
        
        # Analogous to errorCount
        $warningCount=0; 
        
        # If an error occurs during the cvs update, this 
        # indicates that updatelog.txt should be attached
        $updateError=0; 
        
        # If an error occurs during the compilation of 
        # Trilinos, this indicates that trilinosCompileLog.txt
        # should be attached        
        $trilinosCompileError=0; 
        
        # If an error occurs during the compilation of one of
        # the tests, this indicates that testCompileLog.txt
        # should be attached        
        $testCompileError=0; 
        
        open (ERRORS, ">>errors.txt")       or die "$! error trying to open file";        
        open (WARNINGS, ">>warnings.txt")   or die "$! error trying to open file";   
        open (BODY, ">>EmailBody.txt")      or die "$! error trying to open file";
    } # prepareFiles()

    ############################################################################
    # lamboot()
    #
    # Check if any tests are parallel tests--if so, do a lamboot.
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: 
    #   - returns: 

    sub lamboot {  
        chdir "$options{'HOMEDIR'}[0]";
        
        if ($mpi) {
            unless (! -f $options{'HOST_FILE'}[0]) {
                system "lamboot $options{'HOST_FILE'}[0] -v";
            } else {
                system "lamboot";
            }
        }
    } # lamboot()

    ############################################################################
    # run()
    #
    # Core test-harness logic. Calls compile(), build(), test(), and sendMail()
    #   - global variables used: yes
    #   - sends mail: yes
    #   - args: 
    #   - returns: 

    sub run {  
        chdir "$options{'HOMEDIR'}[0]";
        
        # run...
        
    } # run()

    ############################################################################
    # compile()
    #
    # Compile the required Trilinos builds for tests determined above
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: 
    #   - returns: 

    sub compile {   
    
# >>> # why is this explicity set?
        $comm = "serial";
        my $compileFail = 0;
        
        # for each build dir (as specified in config)
        for (my $i=0; $i <= $#{$options{'BUILD_DIRS'}}; $i++) {
        
            # descend into build dir
            chdir"$options{'HOMEDIR'}[0]/$options{'BUILD_DIRS'}[$i]";

# >>> # could these all be centrally located under testharness?            
            # test for rx permissions on invoke configure
            if (!-r "invoke-configure") {
                my $error = "$options{'HOMEDIR'}[0]/$options{'BUILD_DIRS'}[$i]/";
                $error += "invoke-configure must be readable\n";
                die "$! $error";
            } if (!-x "invoke-configure") {
                my $error = "$options{'HOMEDIR'}[0]/$options{'BUILD_DIRS'}[$i]/";
                $error += "invoke-configure must be executable\n";
                die "$! $error";
            }

            # open invoke-configure
            open (SCRIPT, "invoke-configure") 
                or die "$! error trying to open file\n$options{'BUILD_DIRS'}[$i]/invoke-configure\n";
                
# >>> # isn't this only reading the first line?
            my $configScript=<SCRIPT>;
            close SCRIPT;
            
# >>> # can this be generalized more and extracted to config?           
            if ($configScript =~ m/--enable-mpi/i) {
        	    $comm="mpi";
        	} else {
        	    $comm="serial";
        	}

# >>> # before this is all preliminary setup
# >>> # this is where we need to start to implement our cascading scheme        	
        	system "make clean";
        	
# >>> # does this have to be stored in a file?
            my $command = "./invoke-configure >> $options{'HOMEDIR'}[0]";
            $command += "/trilinosCompileLog.txt 2>&1";
        	$compileFail += system "$command";
              
# >>> # terminology: compile or configure?
# >>> # what's the difference, what are we doing here?
# >>> # and for that matter: make, make clean, make dist, build...
  
            # compilation failed, skip associated tests, report failure
            if ($compileFail) {
                $errorCount++;
                $trilinosCompileError++;
                print ERRORS "Trilinos configure process failed for the";
                print ERRORS "$options{'BUILD_DIRS'}[$i] test.\n";
                print ERRORS "Will skip associated tests.\n";
# >>> #         # sendMail(); # COMMENTING OUT EMAIL
                $compileFail=0; # reset 
            } 
            
            # compilation succeeded
            else {
        	    
    	        # build failed, skip associated tests, report failure
        		if (!build($i)) {
        	  		$errorCount++;
        	  		$trilinosCompileError++;
        	  		print ERRORS "Trilinos compilation process failed for the ";
        	  		print ERRORS "$options{'BUILD_DIRS'}[$i] test.\nWill skip associated tests.\n";
# >>> #             # sendMail(); # COMMENTING OUT EMAIL
        	  		$compileFail=0; # reset         		        
        		}
        		
        		# compiled and built successfully - test
        		else {
        	        test($i);
        	    } # else
        	    
            } # else (compile success)
        } # for (buildDirs)
    } # compile()

    ############################################################################
    # build()
    #
    # Build the given Trilinos build.
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: $i (number of build directory in $options{BUILD_DIRS})
    #   - returns: 

    sub build { 
        my $i = $_[0];  
        
        return system "make >> $options{'HOMEDIR'}[0]/trilinosCompileLog.txt 2>&1";
        
    } # build()
    
    ############################################################################
    # test()
    #
    # Run tests
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: $i (number of build directory in $options{BUILD_DIRS})
    #   - returns: 

    sub test {    
        my $i = $_[0];        
        my $compileFail = 0;
        
        # locate all test dirs under Trilinos/packages        
        chdir "$options{'HOMEDIR'}[0]";        
# >>> # can we rely on the system having 'find'?
# >>> # we should compile a list of required programs, modules, permissions, etc.
        my @testDirs = `find packages/ -name test -print`;
	          
	    # run all tests 
	    foreach my $testDir (@testDirs) {	        
	        
# ??? #     # for packages that have not been ported to a particular platform,
# ??? #     # the tests associated with such a package could be skipped below
# ??? #     # use unless (uname =.. && $testDir=..) to skip
# ??? #     # Packages excluded below (jpetra, tpetra, etc) do not yet build
# ??? #     # with autotools
	  
	        # exclude unsupported packages and invalid directories
	        unless ($testDir =~ m/^\s+$/ || $testDir =~ m/^$/ || 
	                $testDir =~ m/tpetra/ || $testDir =~ m/jpetra/ ) {
	                
	            # descend into test directory
	            chdir "$options{'HOMEDIR'}[0]/$options{'BUILD_DIRS'}[$i]/$testDir";

                # find potential scripts in test/scripts/<frequency>/<comm>
                my $command = "find $options{'HOMEDIR'}[0]/$options{'BUILD_DIRS'}[$i]/";
                $command += "$testDir/scripts/$options{'FREQUENCY'}[0]/$comm -type f > list_of_files";
	            system "$command";
	            
	            # extract list of potential scripts
	            open FILELIST, "list_of_files";
	            chomp(@potentialScripts=<FILELIST>);
	            close FILELIST;
	            
	            # run each test
# >>> # $potentialScript should be made local and passed to sendMail()
# >>> # all such variables, for that matter, should be made local and passed to sendMail()
	            foreach $potentialScript (@potentialScripts) {
	            
	                chdir "$options{'HOMEDIR'}[0]/$options{'BUILD_DIRS'}[$i]/$testDir";
	                
	                # if potential script file is executable...
		            if (-x $potentialScript ) {
		            
		                # run test script
		                my $command = "$potentialScript $options{'BUILD_DIRS'}[$i] ";
		                $command += "True >> $options{'HOMEDIR'}[0]/testCompileLog.txt 2>&1";
		                $compileFail += system "$command";
		                
		                # script failed
# >>> # not really a 'compile' fail, is it?
# >>> # can this be localized?
		                if ($compileFail) {
		                    $errorCount++;
		                    $testCompileError++;
		                    print ERRORS "Trilinos test suite failed for ";
		                    print ERRORS "$options{'BUILD_DIRS'}[$i] tests.\n\n";
		                }
		                
		                # Create and send email
		                system "rm -f list_of_files";
# >>> #                 # sendMail(); # COMMENTING OUT EMAIL
		                $compileFail=0; # reset 
		            } # if (executable)
		            
	            } # foreach (potentialScript)
	        } # unless (unsupported package)
	    } # foreach (testDir)
    } # test()
    
    ############################################################################
    # lamhalt()
    #
    # Check if any tests were parallel tests--if so, do a lamhalt.
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: 
    #   - returns: 

    sub lamhalt {    
        chdir"$options{'HOMEDIR'}[0]";        
        if ($mpi) { system "lamhalt"; }
    } # lamhalt()

    ############################################################################
    # cleanupFiles()
    #
    # Close filehandles and delete temp files
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: 
    #   - returns: 

    sub cleanupFiles {    
        chdir"$options{'HOMEDIR'}[0]";        
        close BODY;
        close ERRORS;
        close WARNINGS;
    } # cleanupFiles()
    
    ############################################################################
    # sendMail()
    #
    # This subroutine is called when it is time to send the email - either
    # when the tests are complete or when an error occurs from which the
    # script cannot recover
    #   - global variables used: yes
    #   - sends mail: yes
    #   - args: 
    #   - returns: 

    sub sendMail {
        chdir"$options{'HOMEDIR'}[0]";
        
        use MIME::Lite;
        
# >>>   # Warn me that not all calls to SENDEMAIL were commented out
        print "I'M SENDING EMAIL TO SOMEBODY!"; 

        # Compile list of mail recipients ######################################

# ???   # The following line constructs the name of the log that is produced 
# vvv   # during tests whether errors occurred or not.  (-v not used.)
        my $scriptowner;
        if (-f "log$hostArch.txt") {
            system "mv $options{'HOMEDIR'}[0]/log$hostArch.txt EmailBody.txt";
            open (OWNER, "EmailBody.txt") or die "$! error trying to open file";
            $scriptowner=<OWNER>;
            close OWNER;
        }
# ^^^
# ???   
# ???
# vvv     
        # List addresses for all emails to be sent to below
        chomp $scriptowner;
        my $mailTo;
        # A valid email address must contain an '@'
        if ($scriptowner =~ m/@/i) {
            $mailTo = join ", ", 'jmwille@sandia.gov',$scriptowner 
                or die "$! error forming 'To:' field";
        } else {
            # List where emails without a script owner should be sent
            # By sending them to the trilinos-regresstion list, we assure that
            # results from tests without script owners and cvs and compilation 
            # errors are recorded in Mailman

            $mailTo = join ", ", 'trilinos-regression@software.sandia.gov';
        }
# ^^^
# ???

        # Create email #########################################################

        print "\n**$mailTo**\n";
        open (SUMMARY, ">>summary.txt") or die "$! error trying to open file";
        
        # At least one problem occured =========================================
        if ($errorCount || -f "$options{'HOMEDIR'}[0]/logErrors.txt" ||  
            -f "$options{'HOMEDIR'}[0]/logMpiErrors.txt" || $warningCount || 
            $updateError || $testCompileError || $trilinosCompileError) {

# >>> # why the initialization and the else?
            # create subject line
            my $mailSubject = "$hostArch - $hostName - At least one error occurred";
            if ($updateError) {
                $mailSubject = "$hostArch - $hostName - CVS update error";
            } elsif ($trilinosCompileError) {
                $mailSubject = "$hostArch - $hostName - $comm - Trilinos compile error";
            } else {
                $mailSubject = "$hostArch - $hostName - $comm - At least one test failed";
            }
    
            # construct message
            my $msg;
            $msg=MIME::Lite->new(
    			From =>     'trilinos-regression@software.sandia.gov',
# >>> #			# To =>     'jmwille@sandia.gov',
    			To =>       $mailTo,
    			Subject =>  $mailSubject,
    			Type =>     'multipart/mixed'
    			);
            $msg->attach(Type=>'TEXT', Path=>'summary.txt', Disposition=>'inline');
    	        
    	    # append script name
            if ($potentialScript) {
                print SUMMARY "Test script: ";
                print SUMMARY $potentialScript;
                print SUMMARY "\n";
            }
        
            # append frequency
            if ($options{'FREQUENCY'}[0] =~ m/weekly/) {
                print SUMMARY "Results of weekly test:\n";
            } elsif ($options{'FREQUENCY'}[0] =~ m/daily/) {
                print SUMMARY "Results of daily test:\n";
            } else {
                print SUMMARY "Warning: Unknown test frequency, results below:\n";
            }
  
            # parallel errors
            if ( -f "$options{'HOMEDIR'}[0]/logMpiErrors.txt") {
                print SUMMARY "--> Parallel testing did not complete successfully.\n";

                if ($testCompileError) {
                    print SUMMARY "This failure is probably due to the compile ";
                    print SUMMARY "time errors listed in the attachment \"testCompileLog.txt\".\n";
                    print SUMMARY "If that does not appear to be the case, ";
                }

                print SUMMARY "See attachment \"logMpiErrors.txt\".\n\n";
                $msg->attach(Type=>'TEXT', Path=>'logMpiErrors.txt', Disposition=>'attachment');
            }

            # serial errors
            if ( -f "$options{'HOMEDIR'}[0]/logErrors.txt") {
                print SUMMARY "--> Serial testing did not complete successfully.\n";

                if ($testCompileError) {
                    print SUMMARY "This failure is probably due to the compile ";
                    print SUMMARY "time errors listed in the attachment \"testCompileLog.txt\".\n";
                    print SUMMARY "If that does not appear to be the case, ";
                }
                print SUMMARY "See attachment \"logErrors.txt\".\n\n";
                $msg->attach(Type=>'TEXT', Path=>'logErrors.txt', Disposition=>'attachment');
            }

            # test error
            if ($testCompileError && -f "testCompileLog.txt") {
                $msg->attach(Type=>'TEXT', Path=>'testCompileLog.txt', Disposition=>'attachment');
            }

            # errors            
            if ($errorCount && -f "errors.txt") {
                print SUMMARY "--> For additional information about errors that occured,\n";
                print SUMMARY "See attachment \"errors.txt\".\n\n";
                $msg->attach(Type=>'TEXT', Path=>'errors.txt', Disposition=>'attachment');
            }
        
            # warnings
            if ($warningCount && -f "warnings.txt") {
                print SUMMARY "--> At least one non-testing warning occurred ";
                print SUMMARY "(i.e. script/paramter related).\n";
                print SUMMARY "See attachment \"warnings.txt\".\n\n";
                $msg->attach(Type=>'TEXT', Path=>'warnings.txt', Disposition =>'attachment');
            }
  
            # update error
            if ($updateError && -f "updatelog.txt") {
                print SUMMARY "--> At least one error occurred during CVS update.\n";
                print SUMMARY "See attachment \"updatelog.txt\".\n\n";
                $msg->attach(Type=>'TEXT', Path=>'updatelog.txt', Disposition=>'attachment');
            }
  
            # trilinos compile error
            if ($trilinosCompileError && -f "trilinosCompileLog.txt") {
                print SUMMARY "--> At least one error occurred while compiling Trilinos.\n";
                print SUMMARY "See attachment \"trilinosCompileLog.txt\".\n\n";
                $msg->attach(Type=>'TEXT', Path=>'trilinosCompileLog.txt', Disposition=>'attachment');
            }
        
            # preface script output
            print SUMMARY "************************************************\n";
            print SUMMARY "The following is the output from the test script listed ";
            print SUMMARY "above (or failed compile attempt).\n";
            print SUMMARY "Please note that the -v option was not selected for this log.\n";
            print SUMMARY "NOTE:Depending on your Mail User Agent (MUA), the test ";
            print SUMMARY "summary will either appear below, or it will be attached ";
            print SUMMARY "as a file called \"EmailBody.txt\".\n";
            print SUMMARY "See any attachments listed above for more details.\n";
            print SUMMARY "*************************************************\n";
        
            # Send message =====================================================
            if (-f "EmailBody.txt") {
                $msg->attach(Type =>'TEXT', Path=>'EmailBody.txt', Disposition=>'inline');
            }        
            $msg->send;
            
            system "rm -f logErrors.txt logMpiErrors.txt";

        # Everything successful ================================================
        } else {
            # No problems occurred of any kind
            my $mailSubject = "$hostArch - $hostName - $comm - All tests passed";
            my $msg=MIME::Lite->new(
			    From =>'trilinos-regression@software.sandia.gov',
			    To => $mailTo,
			    Subject => $mailSubject,
			    Type =>'multipart/mixed'
			    );
            $msg->attach(Type=>'TEXT',
	            Path=>'summary.txt',
	            Disposition=>'inline'
	            );
	     
            if ($potentialScript) {
                print SUMMARY "Test script: ";
                print SUMMARY $potentialScript;
                print SUMMARY "\n";
            }
    
            if ($options{'FREQUENCY'}[0] =~/weekly/) {
                print SUMMARY "\nResults of weekly tests:\n";
            } elsif ($options{'FREQUENCY'}[0] =~/daily/) {
                print SUMMARY "\nResults of daily tests:\n";
            } else {
                print SUMMARY "\nWarning: Unknown test frequency, results below:\n";
            }

            print SUMMARY "*****************************************************\n";
            print SUMMARY "The following is the output from the test script listed above.\n";
            print SUMMARY "Please note that the -v option was not selected for this log.\n";
            print SUMMARY "While no errors occurred during this test, this log can ";
            print SUMMARY "still be examined to see which tests were run.\n";
            print SUMMARY "NOTE:Depending on your Mail User Agent (MUA), the test ";
            print SUMMARY "summary will either appear below, or it will be attached ";
            print SUMMARY "as a file called 'EmailBody.txt'\n";
            print SUMMARY "******************************************************";

            # Compile and send message =========================================
            $msg->attach(Type =>'TEXT',
	            Path=>'EmailBody.txt',
		        Disposition=>'inline'
	            );
            $msg->send;
            
        } # else

        # Clean up #############################################################
                
        close SUMMARY;    
        
        # Reset error/warning indicators 
        $errorCount=0;
        $warningCount=0;
        $updateError=0;
        $trilinosCompileError=0;
        $testCompileError=0;
        
        # remove temporary files
        system "rm -f errors.txt warnings.txt updatelog.txt EmailBody.txt";
        system "rm -f trilinosCompileLog.txt testCompileLog.txt summary.txt";    
        
    } # sendMail()

    ############################################################################
    # printUsage()
    #
    # Prints Test-Harness usage to standart output and exits.
    #   - global variables used: no
    #   - sends mail: no
    #   - args: 
    #   - returns: 

    sub printHelp {
        print "Test-Harness\n";
        print "\n";
        print "options:\n";
        print "  -f <file> : (required) test-harness-config file to use\n";
        print "  -g        : generate template configuration file (with default values) and exit\n";
        print "  -s        : omit comments from generated configuration file\n";
        print "  -h        : print this help page and exit\n";
        exit;
    } # printHelp()

    ############################################################################
    # parseConfig()
    #
    # Parses test-harness-config and returns hash of arrays of the form
    # ({VARIABLE_A, [valueA1, valueA2, ...]}, {VARIABLE_B, [valueB1, ...]})
    #   - global variables used: no
    #   - sends mail: no
    #   - args: String test-harness-config filename
    #   - returns: Hash of Arrays options

    sub parseConfig {
        my $filename = $_[0];
        my %options; # Hash of Arrays
        my $line;
        my $name;
        my $value;
        my $continue;
        
        open (CONFIG, "<$filename")
            or die "can't open $filename";
            
        while ($line = <CONFIG>) {
            chomp($line);       # trim newline
            $line =~ s/^\s+//;  # trim leading spaces
            $line =~ s/\s+$//;  # trim trailing spaces
            
            if (!$line || $line =~ m/^[# ]/) {
                # skip comments and blank lines
            } else {
                if ($continue) {    # previous line ended with a (\) continuation char
                    $line =~ m/^((\S*?\s*?|".*?"\s*?)*)\s*?(\\?)$/;
                    $value = $1;
                    $continue = $3;                    
                    # print "\$value: $value\n\n"; # debugging   
                                     
                } else {            # expecting a $name [+]= $value [$value ...] pair
                    $line =~ m/^(\S+)\s*?\+?=\s*((\S*?\s*?|".*?"\s*?)*)\s*?(\\?)$/;
                    $name = $1;
                    $value = $2;
                    $continue = $4;                    
                    # print "\$name: $name, \$value: $value\n\n"; # debugging
                }
                
                if (!exists $options{$name}) {  # if there isn't an option with this $name...
                    $options{$name} = ();       # add a new one
                }
                
                while ($value) {
                    $value =~ s/^\s+//;  # trim leading spaces
                    $value =~ s/\s+$//;  # trim trailing spaces
                    # print "\$value: $value\n"; # debugging
                    
                    $value =~ m/^(".*?"|\S+)/;      # grab leftmost value
                    my $v = $1;                     # store temporarily in $v
                    $value =~ s/^$v//;              # remove $v from remaining list 
                    $v =~ s/"//g;                   # remove any quotes from $v
                    push (@{$options{$name}}, $v);  # add $v to %options{$name}     
                }
            } # else (non-comment/blank)
        } # while ($line)
        
        close CONFIG;
        
        # convert psuedo-variables =============================================
        
        # convert <USER_HOME_DIR> psuedo-variable
        if (defined $options{'HOMEDIR'} && defined $options{'HOMEDIR'}[0]) { 
            $options{'HOMEDIR'}[0] =~ s/<USER_HOME_DIR>/$ENV{'HOME'}/;
        } else {
            print "Error: no HOMEDIR defined.\n";
            exit;
        }
        
        # convert <HOMEDIR> psuedo-variable
        for my $name (keys %options) {
            for my $i (0 .. $#{$options{$name}}) {
                $options{$name}[$i] =~ s/<HOMEDIR>/$options{'HOMEDIR'}[0]/;
            }         
        }
        
        # convert <BASE_HOST_FILE> psuedo-variable
        if (defined $options{'HOST_FILE'} && defined $options{'HOST_FILE'}[0]) {
            if (defined $options{'BASE_HOST_FILE'} && defined $options{'BASE_HOST_FILE'}[0]) { 
                $options{'HOST_FILE'}[0] =~ s/<BASE_HOST_FILE>/$options{'BASE_HOST_FILE'}[0]/;
            }
        } 
        
        # Validate options =====================================================
# >>> # should we reset bad values to default values or just bail?
        
        # check for valid build directories and check for --enable-mpi
        if (defined $options{'BUILD_DIRS'} && defined $options{'BUILD_DIRS'}[0]) {
            for my $i (0 .. $#{$options{'BUILD_DIRS'}}) {
                chdir "$options{'HOMEDIR'}[0]/$options{'BUILD_DIRS'}[$i]";            
                print "$options{'HOMEDIR'}[0]/$options{'BUILD_DIRS'}[$i]\n";
                            
    	        open (INVOKE_CONFIGURE, "invoke-configure") or die "$! error trying to open file";
    	        my $file=<INVOKE_CONFIGURE>;
        	    close INVOKE_CONFIGURE;
        	    
        	    # Check if the configure script indicates a parallel or serial build
        	    if ($file =~ m/--enable-mpi/i) { $mpi = 1; }
            }
        } else {
            $warningCount++;
            print WARNINGS "WARNING: No valid BUILD_DIRS found in config file.\n";
# >>> # set default or just bail?
        }        
        
        # print %options # debugging
        #
        # print "\n\%options:\n\n";
        # for my $name (keys %options) {
        #     my $numElements = $#{$options{$name}};
        #     print "  $name (".($numElements+1)."): \n";
        #     for my $i (0 .. $numElements) {
        #         print "    $i = $options{$name}[$i]\n";
        #     }         
        # }
        
        return %options;        
    } # parseConfig()
    
    ############################################################################
    # genConfigTemp()
    #
    # Generates Test-Harness template config file named in current directory
    # named "test-harness-config" and exits.
    #   - global variables used: no
    #   - sends mail: no
    #   - args: boolean isShort;
    #   - returns: 

    sub genConfigTemp {
        my $short = $_[0];
        my $outFile;
                
        open (outFile, "> test-harness-config-template")
            or die "can't open test-harness-config-template";
        
        print outFile "# test-harness-config\n";  
        print outFile "\n";

        if (!$short) {
            print outFile "################################################################################\n";
            print outFile "# This file contains the settings to be used by the test-harness\n";
            print outFile "#\n";
            print outFile "# All text after a hash (#) is considered a comment and will be ignored\n";
            print outFile "# The format is:\n";
            print outFile "#     TAG = value [value ...]\n";
            print outFile "# For lists items can also be appended using:\n";
            print outFile "#     TAG += value [value ...]\n";
            print outFile "# Values that contain spaces should be placed between double quotes (\"\")\n";
            print outFile "#     (not in single quotes and double quotes can't be escaped)\n";
            print outFile "# Line continuation characters (\\) are allowed as an alternative to (+=)\n";
            print outFile "#     (for adding more values, not for extending a value)\n";
            print outFile "# Adding duplicate tags with (=) and not (+=) will not block older values, but\n";
            print outFile "#     will instead simply append the new values, synonymous to (+=)\n";
            print outFile "################################################################################\n";
            print outFile "\n";
        }
        
        print outFile "#===============================================================================\n";
        print outFile "# General configuration options\n";
        print outFile "#===============================================================================\n";
        print outFile "\n";
        
        if (!$short) {        
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "#\n";
            print outFile "# === can't we just use ../? === \n";
            print outFile "#\n";
            print outFile "# Location of the directory in which this file is located. Default location is\n";
            print outFile "# just below user's home directory. The value <USER_HOME_DIR> can be used to \n";
            print outFile "# indicate the current user's home directory.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - the HOME system environment variable can be referred to with\n";
            print outFile "#   the value <USER_HOME_DIR>\n";
            print outFile "\n";
        }
        
        print outFile "HOMEDIR                = <USER_HOME_DIR>/Trilinos \n";
        
        if (!$short) {        
            print outFile "\n";
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# Indicate how often this test should be run.\n";
            print outFile "# Recognized values include:\n";
            print outFile "#   DAILY WEEKDAYS WEEKENDS \n";
            print outFile "#   MONDAY TUESDAY WEDNESDAY THURSDAY FRIDAY SATURDAY SUNDAY.\n";
            print outFile "# (Use one value from the first row OR any combination of values from the\n";
            print outFile "# second row)\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: YES\n";
            print outFile "\n";
        }
        
        print outFile "FREQUENCY              = DAILY \n";
        
        if (!$short) {      
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# Set location of CVS command on this system.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "\n";
        }
        
        print outFile "CVS_CMD                = cvs \n";
        
        if (!$short) {       
            print outFile "\n"; 
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "#\n";
            print outFile "# === Unsure about this option ===\n";
            print outFile "#\n";
            print outFile "# (LAM only) name of the file containing the names of the machines to be used\n"; 
            print outFile "# for parallel jobs. If this file doesn't exist, parallel jobs will be run on\n";
            print outFile "# the local machine only.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - HOMEDIR can be referred to with the value <HOMEDIR>\n";
            print outFile "\n";
        }
        
        print outFile "BASE_HOST_FILE         = hostfile \n";
        
        if (!$short) {      
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "#\n";
            print outFile "# === Unsure about this option ===\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - HOMEDIR can be referred to with the value <HOMEDIR>\n";
            print outFile "# - BASE_HOST_FILE can be referred to with the value <BASE_HOST_FILE>\n";            
            print outFile "\n";
        }
        
        print outFile "HOST_FILE              = <HOMEDIR>/<BASE_HOST_FILE> \n";
        
        if (!$short) {        
            print outFile "\n";
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# List the names of all of the build directories that should be configured,\n";
            print outFile "# compiled and tested for the test harness.  Each build directory must be a\n"; 
            print outFile "# subdirectory of 'Trilinos/'.  One subdirectory per line.  Each directory\n"; 
            print outFile "# listed must contain an \"invoke-configure\" file that contains the options\n";
            print outFile "# that should be passed to configure.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: YES\n";
            print outFile "\n";
        }
        
        print outFile "BUILD_DIRS             = MPI SERIAL\n";
        
        if (!$short) {        
            print outFile "\n";
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# List the email addresses to which summaries and error messages will be sent\n";
            print outFile "# by default (i.e. when no script owner exists for a particular test).\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: YES\n";
            print outFile "\n";
        }
        
        print outFile "DEFAULT_EMAILS         = \n";
        
        print outFile "\n";
        print outFile "# end test-harness-config";
            
        close outFile;
        exit;
    } # genConfigTemp()