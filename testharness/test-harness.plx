#!/usr/bin/perl -w
# /Trilinos/testharness/test-harness.plx

################################################################################
# The Trilinos Project - Test Harness
# 
# Jim Willenbring, Mike Phenow, Ken Stanley
#
# For more information, see Trilinos/testharness/README.
################################################################################

use strict;

# Variable Declarations ========================================================

# Options variables
my %flags;                  # command line flags (boolean or string)
my %options;                # config-file hash of arrays (keys, lists of values)
my %summary;                # overall test-harness summary

# Error constants
use constant FILE_SYSTEM_ERROR => 0;            # files absent or have incorrect permissions
use constant SYSTEM_COMMAND_ERROR => 1;         # system command failed or doesn't exist
use constant CONFIG_ERROR => 2;                 # test-harness configure file isn't valid
use constant UPDATE_ERROR => 3;                 # cvs update failed
use constant TRILINOS_CONFIGURE_ERROR => 4;     # Trilinos configure failed
use constant TRILINOS_BUILD_ERROR => 5;         # Trilinos build failed
use constant TEST_COMPILE_ERROR => 6;           # test compile failed
use constant TEST_FAILED => 7;                  # test failed
use constant TEST_PASSED => 8;                  # test passed
use constant SUMMARY => 9;                      # test-harness summary

# Error labels
my @codes;
$codes[FILE_SYSTEM_ERROR] = "file system error";           
$codes[SYSTEM_COMMAND_ERROR] = "system command error";        
$codes[CONFIG_ERROR] = "test-harness config error";
$codes[UPDATE_ERROR] = "cvs udate error";
$codes[TRILINOS_CONFIGURE_ERROR] = "Trilinos configure error";
$codes[TRILINOS_BUILD_ERROR] = "Trilinos build error";
$codes[TEST_COMPILE_ERROR] = "test compile error";
$codes[TEST_FAILED] = "test failed";
$codes[TEST_PASSED] = "test passed";
$codes[SUMMARY] = "test-harness summary";
        
################################################################################
# Execution ####################################################################
################################################################################

# Preparation ==================================================================

# Parse command line flags
parseFlags();

# Parse test-harness-config
parseConfig($flags{f});

# Confirm that all options are acceptable        
validateOptions();

# Update Trilinos from CVS Repository
cvsUpdate();

# Start up MPI implementation if any tests are parallel
mpiStartup();

# Delete and recreate build dirs
prepareBuildDirs();

# Main Execution ===============================================================

# Configure, build, and test--mailing as necessary
run();

# Send summary email
sendMail(SUMMARY);

# Clean Up =====================================================================

# Shut down MPI implementation if any tests are parallel
mpiShutdown();

################################################################################
# Subroutines ##################################################################
################################################################################

    ############################################################################
    # parseFlags()
    #
    # Parse command line flags.
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: 
    #   - returns: 

    sub parseFlags {    

        # grab flabs
        use Getopt::Std;
        getopts("f:p:g:sh", \%flags);
        
        # parse config file and exit
        if ($flags{p}) {                         
            parseConfig($flags{p});                  
            validateOptions();        
            exit;            
        }
        
        # generate config file and exit
        if ($flags{g}) { 
            genConfigTemp($flags{g}, $flags{s}); 
            exit;            
        }
        
        # print help and exit
        if ($flags{h}) { 
            printHelp();
            exit; 
        }
        
        # nonsensical flags passed ---------------------------------------------
        
        # -s flag without -g flag
        if ($flags{s} && !$flags{g}) {
            print "Error: -s option must be accompanied by the -g option.\n";
            exit;
        }
        
        # no -f flag passed (no flags at all)
        if (!$flags{f}) { 
            print "Error: must specify a config file using the -f option.\n"; 
            print "Run test-harness.plx with the -h option for a complete list of options.\n"; 
            exit;
        } 
    } # parseFlags()

    ############################################################################
    # prepareSummary()
    #
    # Prepares summary hash.
    #   - global variables used: yes
    #   - sends mail: on error
    #   - args: 
    #   - returns: 

    sub prepareSummary {    
            
        $summary{FILE_SYSTEM_ERROR} = ();
        $summary{SYSTEM_COMMAND_ERROR} = ();
        $summary{CONFIG_ERROR} = ();
        $summary{UPDATE_ERROR} = ();
        $summary{TRILINOS_CONFIGURE_ERROR} = ();
        $summary{TRILINOS_BUILD_ERROR} = ();
        $summary{TEST_COMPILE_ERROR} = ();
        $summary{TEST_FAILED} = ();
        $summary{TEST_PASSED} = ();
        
    } # prepareSummary()

    ############################################################################
    # cvsUpdate()
    #
    # Updates Trilinos from the CVS repository.
    #   - global variables used: yes
    #   - sends mail: on error
    #   - args: 
    #   - returns: 

    sub cvsUpdate {    
    
        if ($options{'CVS_UPDATE'}[0] =~ m/YES/i) {
            chdir "$options{'TRILINOS_DIR'}[0]";
            
            my $command = "";
            $command .= "$options{'CVS_CMD'}[0] update -dP > $options{'TRILINOS_DIR'}[0]";
            $command .= "/testharness/temp/update_log.txt 2>&1";
            my $result = system $command;
            if ($result) {
                sendMail(UPDATE_ERROR);
                die " *** Error updating Trilinos - aborting test-harness ***\n";
            }
            system "rm -f $options{'TRILINOS_DIR'}[0]/testharness/temp/update_log.txt";
        }
    } # cvsUpdate()

    ############################################################################
    # mpiStartup()
    #
    # Run specified mpi startup command if an mpi build directory is given
    #   - global variables used: yes
    #   - sends mail: on error
    #   - args: 
    #   - returns: 

    sub mpiStartup {  
        chdir "$options{'TRILINOS_DIR'}[0]";                
        if (defined $options{'MPI_DIR'} && defined $options{'MPI_DIR'}[0]) { 
            my $command = "$options{'MPI_STARTUP_CMD'}[0]";
            my $commandFailed = system $command; 
            if ($commandFailed) {
                sendMail(SYSTEM_COMMAND_ERROR, $command);
                print $command;          
                die " *** Error running system command - aborting test-harness ***\n";
            }
        }
    } # mpiStartup()

    ############################################################################
    # prepareBuildDirs()
    #
    # Delete and recreate build dirs
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: 
    #   - returns: 

    sub prepareBuildDirs {  
        chdir "$options{'TRILINOS_DIR'}[0]";  
        
        # delete and recreate mpi dir
        system "rm -rf $options{'MPI_DIR'}[0]";
        system "mkdir $options{'MPI_DIR'}[0]";
        
        # delete and recreate serial dir
        system "rm -rf $options{'SERIAL_DIR'}[0]";
        system "mkdir $options{'SERIAL_DIR'}[0]";
        
    } # prepareBuildDirs()

    ############################################################################
    # run()
    #
    # Core test-harness logic. Calls configure(), build(), test(), and sendMail()
    #   - global variables used: yes
    #   - sends mail: yes
    #   - args: 
    #   - returns: 

    sub run {  
               
        # Configure, build and test each configuration (as given by the list of
        # machine-independent invoke-configure elements) 
        for (my $i=0; $i <= $#{$options{'TRILINOS_CONFIG_FILES'}}; $i++) {        
            
            prepareConfigure($i);
            
            # create array of build dirs =======================================            
            my @buildDir = ();
            if (defined $options{'MPI_DIR'} && defined $options{'MPI_DIR'}[0]) {
                push (@buildDir, $options{'MPI_DIR'}[0]);
            }
            if (defined $options{'SERIAL_DIR'} && defined $options{'SERIAL_DIR'}[0]) {
                push (@buildDir, $options{'SERIAL_DIR'}[0]);
            }
            
            # for each build directory =========================================
            for (my $j=0; $j <= $#buildDir; $j++) {

                my $comm; 
                if ($buildDir[$j] eq $options{'MPI_DIR'}[0]) {
                    $comm = "mpi";
                } elsif ($buildDir[$j] eq $options{'SERIAL_DIR'}[0]) {
                    $comm = "serial";
                } 
                
                # descend into build dir
                chdir"$options{'TRILINOS_DIR'}[0]/$buildDir[$j]";
                
                # test for existence of invoke-configure
                if (!-f "invoke-configure") {
                    my $message = "";
                    $message .= "$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/";
                    $message .= "invoke-configure must be present\n";
                    sendMail(FILE_SYSTEM_ERROR, $message, $comm);
                    print $message;
                    die " *** File missing - aborting test-harness ***\n";
                } 
                
                # test for read permissions on invoke-configure
                if (!-r "invoke-configure") {
                    my $message = "";
                    $message .= "$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/";
                    $message .= "invoke-configure must be readable\n";
                    sendMail(FILE_SYSTEM_ERROR, $message, $comm);
                    print $message;
                    die " *** File permission wrong - aborting test-harness ***\n";
                } 
                
                # test for executable permission on invoke-configure
                if (!-x "invoke-configure") {
                    my $message = "";
                    $message .= "$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/";
                    $message .= "invoke-configure must be executable\n";
                    sendMail(FILE_SYSTEM_ERROR, $message, $comm);
                    print $message;
                    die " *** File permission wrong - aborting test-harness ***\n";
                }
                
                # configure ----------------------------------------------------
                
                my $configureFailed = configure($buildDir[$j]);
                if ($configureFailed) {    
                    sendMail(TRILINOS_CONFIGURE_ERROR, $configureFailed, $comm); 
                } else {     
                    print "$comm - Trilinos configured successfully.\n";                      
                }
                
                # build --------------------------------------------------------
                
                my $buildFailed = build($buildDir[$j]);
                if ($buildFailed) {
                    sendMail(TRILINOS_BUILD_ERROR, $buildFailed, $comm); 
                } else {    
                    print "$comm - Trilinos built successfully.\n";                  
                }
                
                # test ---------------------------------------------------------
                
                # locate all test dirs under Trilinos/packages        
                chdir "$options{'TRILINOS_DIR'}[0]/$buildDir[$j]";        
                my @testDirs = `find packages/ -name test -print`; 
                
        	    # run all tests 
        	    foreach my $testDir (@testDirs) {
                    $testDir =~ s/\s*$//;  # trim trailing whitespace
                                    
        	        # exclude unsupported packages and invalid directories
        	        unless ($testDir =~ m/^\s+$/ || $testDir =~ m/^$/ || 
        	                $testDir =~ m/tpetra/ || $testDir =~ m/jpetra/ ) {
        	                
        	            # descend into test directory
                        chdir "$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/$testDir";
                
                        # find potential scripts in test/scripts/<frequency>/<comm>
                        my $command = "";
                        $command .= "find $options{'TRILINOS_DIR'}[0]/$buildDir[$j]/";
                        $command .= "$testDir/scripts/$options{'FREQUENCY'}[0]/$comm -type f";
                        my $output = `$command`;
                        my @potentialScripts = split (/\s+/, $output);
                                
                        # run each test
                        foreach my $potentialScript (@potentialScripts) {
                            $potentialScript =~ s/\s*$//;  # trim trailing whitespace
                            
                            # if potential script file is executable...
                            if (-x $potentialScript) {
                
                                # test
                                my $testFailed = test($buildDir[$j], $potentialScript);
                                if ($testFailed) {                
                                    print "Test failed.\n";  
                                    sendMail(TEST_FAILED, $testFailed, $comm, $testDir, $potentialScript); 
                                } else {                                    
                                    print "Test passed.\n";  
                                    sendMail(TEST_PASSED, $testFailed, $comm, $testDir, $potentialScript);
                                }
                                
                            } # if (executable)
                        } # foreach ($potentialScript)                    
                    } # unless (unsupported)
                } # foreach $testDir
                
                # email --------------------------------------------------------
                
            } # for (buildDirs)                       
        } # for (TRILINOS_CONFIG_FILES)               
    } # run()

    ############################################################################
    # prepareConfigure()
    #
    # Compile the invoke configure files and move them into place
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: $i      (which machine-independent invoke-configure element)
    #   - returns: 

    sub prepareConfigure {   
        my $i = $_[0];    
        	
        # create complete invoke-configure files from elements =================              
        my $rawInvokeConfigure = "";  
        
        # prepare machine-dependent portion of invoke-configure
        if (-f "$options{'TRILINOS_DIR'}[0]/testharness/elements-machine/$options{'MACHINE_CONFIG_FILE'}[0]") {
            open (MACHINE, "<$options{'TRILINOS_DIR'}[0]/testharness/elements-machine/$options{'MACHINE_CONFIG_FILE'}[0]") 
                or die "$! error trying to open file";  
            undef $/;                   # undefine input record separator
	        my $file=<MACHINE>;         # copy entire file
    	    close MACHINE;
    	    $rawInvokeConfigure .= $file."\n";
        } else {
            my $message = "";
            $message .= "$options{'TRILINOS_DIR'}[0]/testharness/elements-machine";
            $message .= "/$options{'MACHINE_CONFIG_FILE'}[0] does not exist\n";
            sendMail(FILE_SYSTEM_ERROR, $message);
            print $message;
            die " *** File missing - aborting test-harness ***\n";
        }   
        
        # prepare machine-independent portion of invoke-configure
        if (-f "$options{'TRILINOS_DIR'}[0]/testharness/elements-trilinos/$options{'TRILINOS_CONFIG_FILES'}[$i]") {
            open (COMPONENT, "<$options{'TRILINOS_DIR'}[0]/testharness/elements-trilinos/$options{'TRILINOS_CONFIG_FILES'}[$i]") 
                or die "$! error trying to open file";  
            undef $/;                   # undefine input record separator
	        my $file=<COMPONENT>;       # copy entire file
    	    close COMPONENT;
    	    $rawInvokeConfigure .= $file;
        } else {
            my $message = "";
            $message .= "$options{'TRILINOS_DIR'}[0]/testharness/elements-machine";
            $message .= "/$options{'TRILINOS_CONFIG_FILES'}[$i] does not exist\n";
            sendMail(FILE_SYSTEM_ERROR, $message);
            print $message;
            die " *** File missing - aborting test-harness ***\n";
        }        
	    
	    $rawInvokeConfigure =~ s/$/ \\/mg;      # append line-continuation to each line
	    $rawInvokeConfigure =~ s/\\$//;         # remove last line-continuation
	                
	    # create and copy MPI invoke configure
	    if (defined $options{'MPI_DIR'} && defined $options{'MPI_DIR'}[0]) {	    
        
            # prepare machine-dependent-mpi portion of invoke-configure
            my $mpiFlags = "";
            if (-f "$options{'TRILINOS_DIR'}[0]/testharness/elements-machine/$options{'MACHINE_MPI_CONFIG_FILE'}[0]") {
                open (MACHINE_MPI, "<$options{'TRILINOS_DIR'}[0]/testharness/elements-machine/$options{'MACHINE_MPI_CONFIG_FILE'}[0]") 
                    or die "$! error trying to open file";  
                undef $/;                   # undefine input record separator
    	        my $file=<MACHINE_MPI>;     # copy entire file
        	    close MACHINE_MPI;
        	    $mpiFlags .= $file."\n";
        	    $mpiFlags =~ s/$/ \\/mg;    # append line-continuation to each line
	            $mpiFlags =~ s/\\$//;       # remove last line-continuation
            } else {
                my $message = "";
                $message .= "$options{'TRILINOS_DIR'}[0]/testharness/elements-machine";
                $message .= "/$options{'MACHINE_MPI_CONFIG_FILE'}[0] does not exist\n";
                sendMail(FILE_SYSTEM_ERROR, $message);
                print $message;
                die " *** File missing - aborting test-harness ***\n";
            }   	    
	    
            # open INVOKE_CONFIGURE_MPI for writing
            open (INVOKE_CONFIGURE_MPI, ">$options{'TRILINOS_DIR'}[0]/testharness/temp/invoke-configure-mpi")
                or die "$! error trying to open file";
            print INVOKE_CONFIGURE_MPI ".././configure ";
            print INVOKE_CONFIGURE_MPI $mpiFlags;
            print INVOKE_CONFIGURE_MPI $rawInvokeConfigure;
            close INVOKE_CONFIGURE_MPI;
            
            # move invoke-configure file into place
            my $command;                     
            $command = "mv $options{'TRILINOS_DIR'}[0]/testharness/temp/invoke-configure-mpi ";
            $command .= "$options{'TRILINOS_DIR'}[0]/$options{'MPI_DIR'}[0]/invoke-configure";
            system $command;
            
            # set invoke-configure permissions
            system "chmod a+rx $options{'TRILINOS_DIR'}[0]/$options{'MPI_DIR'}[0]/invoke-configure";
        }
        
        # create and copy SERIAL invoke configure
	    if (defined $options{'SERIAL_DIR'} && defined $options{'SERIAL_DIR'}[0]) {
            # open INVOKE_CONFIGURE_SERIAL for writing
            open (INVOKE_CONFIGURE_SERIAL, ">$options{'TRILINOS_DIR'}[0]/testharness/temp/invoke-configure-serial")
                or die "$! error trying to open file";
            print INVOKE_CONFIGURE_SERIAL ".././configure \\\n";
            print INVOKE_CONFIGURE_SERIAL $rawInvokeConfigure;
            close INVOKE_CONFIGURE_SERIAL;
            
            # move invoke-configure file into place
            my $command;                     
            $command = "mv $options{'TRILINOS_DIR'}[0]/testharness/temp/invoke-configure-serial ";
            $command .= "$options{'TRILINOS_DIR'}[0]/$options{'SERIAL_DIR'}[0]/invoke-configure";
            system $command;
            
            # set invoke-configure permissions
            system "chmod a+rx $options{'TRILINOS_DIR'}[0]/$options{'SERIAL_DIR'}[0]/invoke-configure";
        }
        
    } # prepareConfigure()

    ############################################################################
    # configure()
    #
    # Compile the required Trilinos builds for tests determined above
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: $buildDir (name of build directory)
    #   - returns: 

    sub configure {   
        my $buildDir = $_[0];    
        	
        my $command = "";
        $command .= "./invoke-configure >> $options{'TRILINOS_DIR'}[0]";
        $command .= "/testharness/temp/trilinos_configure_log.txt 2>&1";
    	return system $command;
        
    } # configure()

    ############################################################################
    # build()
    #
    # Build the given Trilinos build.
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: $buildDir (name of build directory)
    #   - returns: 

    sub build { 
        my $buildDir = $_[0];    
    
        my $command = "";
        $command .= "make $options{'MAKE_FLAGS'}[0] >> $options{'TRILINOS_DIR'}[0]";
        $command .= "/testharness/temp/trilinos_build_log.txt 2>&1";
        return system $command;
        
    } # build()
    
    ############################################################################
    # test()
    #
    # Run tests
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: $buildDir (name of build directory) 
    #           $script (test script to run)
    #   - returns: 

    sub test {    
        my $buildDir = $_[0];   
        my $script = $_[1];
        
        # run test script
        my $command = "";
        $command .= "$script $buildDir True >> $options{'TRILINOS_DIR'}[0]";
        $command .= "/testharness/temp/test_compile_log.txt 2>&1";
        return system $command;
        
    } # test()
    
    ############################################################################
    # mpiShutdown()
    #
    # Run the specified mpi shutdown command if an mpi build directory is given
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: 
    #   - returns: 
        
    sub mpiShutdown {    
        chdir"$options{'TRILINOS_DIR'}[0]";        
        if (defined $options{'MPI_DIR'} && defined $options{'MPI_DIR'}[0]) { 
            system "$options{'MPI_SHUTDOWN_CMD'}[0]"; 
        }
    } # mpiShutdown()
    
    ############################################################################
    # sendMail()
    #
    # This subroutine is called when it is time to send the email - either
    # when the tests are complete or when an error occurs from which the
    # script cannot recover
    #   - global variables used: yes
    #   - sends mail: yes
    #   - args: $code       (global constant describing error)
    #           $message    (message accompanying error)
    #           $comm       (current comm, if relevant)
    #           $testDir    (current package directory, if relevant)
    #           $testName   (name of test, if code is test-related)
    #   - returns: 

    sub sendMail {
        my $code = $_[0];     
        my $message = $_[1];      
        my $comm = $_[2];  
        my $testDir = $_[3];
        my $testName = $_[4];
        
        print "\n=== SENDING EMAIL ===\n\n"; 
        
        # preparation ##########################################################
        
        chdir "$options{'TRILINOS_DIR'}[0]/testharness/";
        
        use lib "./lib/MIME-Lite-3.01/lib";     # cvs checkin MIME::Lite
        use MIME::Lite;                         # might need Net::SMTP
        
        chdir "$options{'TRILINOS_DIR'}[0]/testharness/temp";
            
        # host information
        chomp (my $hostOS=`uname`);             # host operating system
        chomp (my $hostOSRelease=`uname -r`);   # host operating system release
        chomp (my $hostOSVersion=`uname -v`);   # host operating system version 
        chomp (my $hostHardware=`uname -m`);    # host hardware 
        chomp (my $hostName=`uname -n`);        # host name    
        
        # remove extra newlines
        $hostOS =~ s/\s*$//; 
        $hostOSRelease =~ s/\s*$//; 
        $hostOSVersion =~ s/\s*$//;  
        $hostHardware =~ s/\s*$//; 
        $hostName =~ s/\s*$//; 
        
        # remove path from test name
        if (defined $testName) {
            $testName =~ s/.*\///;
        }
        
        # extract summary data -------------------------------------------------        
        if ($code == FILE_SYSTEM_ERROR || $code == SYSTEM_COMMAND_ERROR ||
               $code == CONFIG_ERROR || $code == UPDATE_ERROR) {
            push (@{$summary{$code}}, $message);
        }         
        if ($code == TRILINOS_CONFIGURE_ERROR || $code == TRILINOS_BUILD_ERROR) {
            push (@{$summary{$code}}, $comm);
        }        
        if ($code == TEST_COMPILE_ERROR || $code == TEST_FAILED || $code == TEST_PASSED) {
            push (@{$summary{$code}}, $testName);
        }

        # compile list of mail recipients --------------------------------------

        # extract $scriptOwner email from first line of log$hostOS.txt
        # (note: the name "log$hostOS.txt" is functional--it is written-to 
        # by the test scripts.
        my $scriptOwner = "";
        if (-f "$options{'TRILINOS_DIR'}[0]/log$hostOS.txt") {
            open (OWNER, "$options{'TRILINOS_DIR'}[0]/log$hostOS.txt") or die "$! error trying to open file";
            $scriptOwner=<OWNER>;   # read first line (email of script owner is on first line)
            chomp $scriptOwner;     # trim newline
            close OWNER;
        }
        
        my $mailTo;
        if ($scriptOwner =~ m/\S+?\@\S+?/i) {    # look for ...@...
            $mailTo = join ", ", $scriptOwner, "$options{'ALL_SCRIPTS_EMAIL'}";
        } else {
            $mailTo = join ", ", "$options{'NO_SCRIPT_OWNER_EMAIL'}", "$options{'ALL_SCRIPTS_EMAIL'}";
        }
        
        # OVERRIDE OTHER EMAILS
        $mailTo = 'mnphenow@software.sandia.gov';

        # Create email #########################################################

        # subject ==============================================================
        
        my $subject = "$hostOS - $hostName - ";
        if (defined $comm) {$subject .= "$comm - ";}
        if (defined $testName) {$subject .= "$testName - ";}
        $subject .= $codes[$code];
        
        # construct email ======================================================    
        
        my $email;
        $email=MIME::Lite->new(
			From =>     'Trilinos test-harness <trilinos-regression@software.sandia.gov>',
			To =>       $mailTo,
			Subject =>  $subject,
			Type =>     'multipart/mixed'
		);
            
        # body =================================================================
        
        my $body = "";
                
        # header section
        $body .= "Host OS:          $hostOS\n";
        $body .= "Host OS Release:  $hostOSRelease\n";
        $body .= "Host OS Version:  $hostOSVersion\n";
        $body .= "Host Hardware:    $hostHardware\n";
        $body .= "Host Name:        $hostName\n";
        if (defined $comm) {
            $body .= "\n";
            $body .= "Comm:             $comm\n";}
        if (defined $testDir) 
            {$body .= "Test Directory:   ...$testDir...\n";}
        if (defined $testName) 
            {$body .= "Test Name:        $testName\n";}
        if (defined $testName) 
            {$body .= "Frequency:        $options{'FREQUENCY'}[0]\n";}
        $body .= "\n";        
        $body .= "Result:           $codes[$code]\n";
        
        $body .= "\n";        
        
        # fatal error ----------------------------------------------------------
        if ($code == FILE_SYSTEM_ERROR || $code == SYSTEM_COMMAND_ERROR ||
            $code == CONFIG_ERROR || $code == UPDATE_ERROR) {
        
            $body .= "------------------------------------------------------------\n";
            $body .= "FATAL ERROR\n";
            $body .= "\n";       
            $body .= "This error caused the test-harness to quit prematurely. It is\n";
            $body .= "the sort of errorthat the test-harness is not designed to\n";
            $body .= "recover from. It is probably either asimple human oversight,\n";
            $body .= "or once fixed, should not occur again on this machine.\n";    
            $body .= "\n";       
            $body .= "Message/exit status:      $message\n";           
            $body .= "\n";        
        }
        
        # attachments ----------------------------------------------------------
        
        my $attachmentText = "";
        my $attachment = 0;
        
        # update failed
        if ($code == UPDATE_ERROR && -f "update_log.txt") {
            $attachment = 1;
            $attachmentText .= "    update_log.txt\n";
            $email->attach(Type=>'TEXT', Path=>'update_log.txt', Disposition=>'attachment');
        }
        
        # trilinos configure failed
        if ($code == TRILINOS_CONFIGURE_ERROR && -f "trilinos_configure_log.txt") {
            $attachment = 1;
            $attachmentText .= "    trilinos_configure_log.txt\n";
            $email->attach(Type=>'TEXT', Path=>'trilinos_configure_log.txt', Disposition=>'attachment');
        }       
        
        # trilinos build failed
        if ($code == TRILINOS_BUILD_ERROR && -f "trilinos_build_log.txt") {
            $attachment = 1;
            $attachmentText .= "    trilinos_build_log.txt\n";
            $email->attach(Type=>'TEXT', Path=>'trilinos_build_log.txt', Disposition=>'attachment');
        }
        
        # test compile failed
        if ($code == TEST_COMPILE_ERROR && -f "test_compile_log.txt") {
            $attachment = 1;
            $attachmentText .= "    test_compile_log.txt\n";
            $email->attach(Type=>'TEXT', Path=>'test_compile_log.txt', Disposition=>'attachment');
        }
        
        # parallel test failed
        if (-f "$options{'TRILINOS_DIR'}[0]/logMpiErrors.txt") {   
            $attachment = 1;
            $attachmentText .= "    logMpiErrors.txt\n";
            $email->attach(Type=>'TEXT', Path=>"$options{'TRILINOS_DIR'}[0]/logMpiErrors.txt", Disposition=>'attachment');
        }
        
        # serial test failed
        if (-f "$options{'TRILINOS_DIR'}[0]/logErrors.txt") {
            $attachment = 1;
            $attachmentText .= "    logErrors.txt\n";
            $email->attach(Type=>'TEXT', Path=>"$options{'TRILINOS_DIR'}[0]/logErrors.txt", Disposition=>'attachment');
        }
        
        if (-f "$options{'TRILINOS_DIR'}[0]/log$hostOS.txt") {
            $attachment = 1;
            $attachmentText .= "    log$hostOS.txt\n";
            $email->attach(Type =>'TEXT', Path=>"$options{'TRILINOS_DIR'}[0]/log$hostOS.txt", Disposition=>'attachment');
        }             
        
        if ($attachment) {
            $body .= "------------------------------------------------------------\n";
            $body .= "Attachments: \n";
            $body .= "\n";            
            $body .= "$attachmentText";
            $body .= "\n";        
        }
        
        # notes ----------------------------------------------------------------          
        
        if ($code == TEST_FAILED || $code == TEST_PASSED) {
            $body .= "------------------------------------------------------------\n";
            $body .= "Notes: \n";
            $body .= "\n";        
            
            if ($code == TEST_FAILED) {
                $body .= "log$hostOS.txt is the output from the test script listed\n";
                $body .= "above (or failed compile attempt). Please note that the -v\n";
                $body .= "option was not selected for this log.\n"; 
            }  
            
            if ($code == TEST_PASSED) {
                $body .= "log$hostOS.txt is the output from the test script listed\n";
                $body .= "above. Please note that the -v option was not selected for\n";
                $body .= "this log. While no errors occurred during this test, this\n";
                $body .= "log can still be examined to see which tests were run.\n";
            }   
            
            $body .= "\n";        
        }
        
        # print the summary hash of arrays -------------------------------------
        
        if ($code == SUMMARY) {
            $body .= "------------------------------------------------------------\n";
            $body .= "Summary: \n";
            
            for my $id (keys %summary) { 
                my $numElements = $#{$summary{$id}}; 
                $body .= "\n- $codes[$id] (".($numElements+1)."): \n\n";
                for my $i (0 .. $numElements) {
                    $body .= "    $summary{$id}[$i]\n";
                }
            } 
        }           

        # send email ===========================================================
        
		$email->attach(Type=>'TEXT', Data=>$body);
        $email->send();
        
        system "rm -f update_log.txt";
        system "rm -f trilinos_configure_log.txt";
        system "rm -f trilinos_build_log.txt";
        system "rm -f test_compile_log.txt";
        system "rm -f $options{'TRILINOS_DIR'}[0]/logErrors.txt";
        system "rm -f $options{'TRILINOS_DIR'}[0]/logMpiErrors.txt";
        system "rm -f $options{'TRILINOS_DIR'}[0]/log$hostOS.txt";
                
    } # sendMail()

    ############################################################################
    # printHelp()
    #
    # Prints Test-Harness usage to standart output and exits.
    #   - global variables used: no
    #   - sends mail: no
    #   - args: 
    #   - returns: 

    sub printHelp {
        print "Trilinos Test-Harness\n";
        print "\n";
        print "Usage:  ./testharness.plx -f FILENAME\n";
        print "\n";
        print "Options:\n";
        print "  -f FILE : Run test harness normally with given test-harness-config file\n";
        print "\n";
        print "  -p FILE : Parse given test-harness-config file and exit. This is useful\n";
        print "            for catching errors and inconsistencies without running the\n";
        print "            entire test-harness. (no ouput indicates a valid config file)\n";
        print "\n";
        print "  -g FILE : Generate template configuration file (with defaults) named \n";
        print "            FILE and exit\n";
        print "\n";
        print "  -s      : Omit comments from generated configuration file\n";
        print "            (must be of the form: -sg FILE) (do not use -gs)\n";
        print "\n";
        print "  -h      : Print this help page and exit\n";
        print "\n";
        print "Notes:\n";
        print "  - Options with FILE require a filename--absolute or relative to\n";
        print "    Trilinos/testharness.\n";
        print "  - On some systems, \"./testharness.plx\" will not work; try\n";
        print "    \"perl testharness.plx\" instead.\n";
        print "  - See README in Trilinos/testharness for more information.\n";
        print "\n";
    } # printHelp()

    ############################################################################
    # parseConfig()
    #
    # Parses test-harness-config and fills global hash of arrays of the form
    # ({VARIABLE_A, [valueA1, valueA2, ...]}, {VARIABLE_B, [valueB1, ...]})
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: filename of config file to parse
    #   - returns: 

    sub parseConfig {
        my $filename = $_[0];
        my $line;
        my $name;
        my $value;
        my $continue;
        
        open (CONFIG, "<$filename")
            or die "can't open $filename";
            
        while ($line = <CONFIG>) {
            $line =~ s/^\s*//;  # trim leading spaces
            $line =~ s/\s*$//;  # trim trailing spaces
            
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
                    $value =~ s/^\s*//;  # trim leading spaces
                    $value =~ s/\s*$//;  # trim trailing spaces
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
        
        # detect current directory =============================================
        #   this should still be ?/Trilinos/testharness
        use Cwd;
        use File::Basename;
        
        my $cwd = getcwd();             # ?/Trilinos/testharness
        $cwd = dirname($cwd);           # ?/Trilinos
        $options{'TRILINOS_DIR'} = ();
        push (@{$options{'TRILINOS_DIR'}}, $cwd);
        
        # detect frequency =====================================================
        $options{'FREQUENCY'} = ();
        my $date = `date`;
        if ($date =~ m/sun/i) {
            push (@{$options{'FREQUENCY'}}, "weekly");
        } else {
            push (@{$options{'FREQUENCY'}}, "daily");
        }            
        
    } # parseConfig()
    
    ############################################################################
    # validateOptions()
    #
    # Validate options given by test-harness-config
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: 
    #   - returns: 

    sub validateOptions {
             
        # convert <TRILINOS_DIR> psuedo-variable
        for my $name (keys %options) {
            for my $i (0 .. $#{$options{$name}}) {
                $options{$name}[$i] =~ s/<TRILINOS_DIR>/$options{'TRILINOS_DIR'}[0]/;
            }         
        }
        
        # if HOST_FILE isn't specified, assign empty string
        if (!defined $options{'HOST_FILE'} || !defined $options{'HOST_FILE'}[0]) {
            $options{'HOST_FILE'} = ();
            push (@{$options{'HOST_FILE'}}, "");
        } 
        
        # convert <HOST_FILE> psuedo-variable
        for my $name (keys %options) {
            for my $i (0 .. $#{$options{$name}}) {
                $options{$name}[$i] =~ s/<HOST_FILE>/$options{'HOST_FILE'}[0]/;
            }         
        }      
        
        # if MAKE_FLAGS are specified, but don't begin with a '-', add one
        if (defined $options{'MAKE_FLAGS'} && defined $options{'MAKE_FLAGS'}[0]) {
            if (! $options{'MAKE_FLAGS'}[0] =~ m/^-/) {
                $options{'MAKE_FLAGS'}[0] =~ s/^/-/;
            }
        } 
        
        # if MAKE_FLAGS weren't specified, assign the empty string
        else {
            $options{'MAKE_FLAGS'} = ();
            push (@{$options{'MAKE_FLAGS'}}, "");
        }
        
        # if SUMMARY_EMAIL are not specified, assign the empty string
        if (!defined $options{'SUMMARY_EMAIL'} || !defined $options{'SUMMARY_EMAIL'}[0]) {
            $options{'SUMMARY_EMAIL'} = ();
            push (@{$options{'SUMMARY_EMAIL'}}, "");
        } 
        
        # if NO_SCRIPT_OWNER_EMAIL are not specified, assign the empty string
        if (!defined $options{'NO_SCRIPT_OWNER_EMAIL'} || !defined $options{'NO_SCRIPT_OWNER_EMAIL'}[0]) {
            $options{'NO_SCRIPT_OWNER_EMAIL'} = ();
            push (@{$options{'NO_SCRIPT_OWNER_EMAIL'}}, "");
        } 
        
        # if ALL_SCRIPTS_EMAIL are not specified, assign the empty string
        if (!defined $options{'ALL_SCRIPTS_EMAIL'} || !defined $options{'ALL_SCRIPTS_EMAIL'}[0]) {
            $options{'ALL_SCRIPTS_EMAIL'} = ();
            push (@{$options{'ALL_SCRIPTS_EMAIL'}}, "");
        } 
        
        # print the options hash of arrays -------------------------------------
        #print "\n\%options:\n\n";                       # debugging
        #for my $name (keys %options) {                  # debugging
        #    my $numElements = $#{$options{$name}};      # debugging
        #    print "  $name (".($numElements+1)."): \n"; # debugging
        #    for my $i (0 .. $numElements) {             # debugging
        #        print "    $i = $options{$name}[$i]\n"; # debugging
        #    }                                           # debugging
        #}                
    } # validateOptions()
    
    ############################################################################
    # genConfigTemp()
    #
    # Generates Test-Harness template config file named in current directory
    # named "test-harness-config" and exits.
    #   - global variables used: no
    #   - sends mail: no
    #   - args: string/boolean filename, boolean short
    #   - returns: 

    sub genConfigTemp {
        my $filename = $_[0];
        my $short = $_[1];
        my $outFile;
                
        open (outFile, "> $filename")
            or die "can't open $filename";
        
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
        print outFile "# Frequently changed configuration options\n";
        print outFile "#===============================================================================\n";
        print outFile "\n";
        
        if (!$short) {        
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# The name of the file in \"<TRILINOS_DIR>/testharness/elements-machine\"\n";
            print outFile "# containing the machine-dependent configure options for this machine.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES\n";
            print outFile "\n";
        }
        
        print outFile "MACHINE_CONFIG_FILE             = \n";
        
        if (!$short) {        
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# The name of the file in \"<TRILINOS_DIR>/testharness/elements-machine\"\n";
            print outFile "# containing the machine-independent mpi configure options for this\n";
            print outFile "# configuration of the test-harness.\n";
            print outFile "\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES (if MPI_DIR supplied)\n";
            print outFile "\n";
        }
        
        print outFile "MACHINE_MPI_CONFIG_FILE         = \n";
        
        if (!$short) {        
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# The name of the files in \"<TRILINOS_DIR>/testharness/elements-trilinos\"\n";
            print outFile "# containing the machine-independent configure options for this configuration\n";
            print outFile "# of the test-harness. List multiple files in order to have the test-harness\n";
            print outFile "# attempt to configure, build, and test each configuration in order.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: YES\n";
            print outFile "# - value required: YES\n";
            print outFile "\n";
        }
        
        print outFile "TRILINOS_CONFIG_FILES           = \n";
        
        if (!$short) {      
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";         
            print outFile "# Indicate how the test-harness should report results. If the value is EMAIL,\n";              
            print outFile "# then emails will be sent to ________. If the value is anything but EMAIL,\n";               
            print outFile "# then no emails will be sent at all.\n";      
            print outFile "#\n";    
            print outFile "# - recognized values: EMAIL <outputFile>\n";          
            print outFile "#   (where <outputFile> is the name of the local file to which the test-harness\n"; 
            print outFile "#   should write the results)\n"; 
            print outFile "# - multiple values recognized: NO\n"; 
            print outFile "# - value required: YES\n";
            print outFile "# - the absolute path of the Trilinos directly above 'testharness' can be\n";
            print outFile "#   referred to with the value <TRILINOS_DIR>\n";
            print outFile "# - relative paths are relative to <TRILINOS_DIR>/testharness\n";   
            print outFile "\n";
        }
        
        print outFile "REPORT_METHOD                   = EMAIL\n";
        
        print outFile "\n";
        print outFile "#===============================================================================\n";
        print outFile "# Less frequently changed configuration options\n";
        print outFile "#===============================================================================\n";
        print outFile "\n";
        
        if (!$short) {        
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# Provide the name of the serial build directory that should be configured,\n";
            print outFile "# compiled and tested by the test harness, or leave blank to indicate that there\n";
            print outFile "# should be no serial build. Directory must be a subdirectory of 'Trilinos/'.\n"; 
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: NO (unless MPI_DIR is omitted)\n";
            print outFile "\n";
        }
        
        print outFile "SERIAL_DIR                      = SERIAL\n";
        
        if (!$short) {    
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n"; 
            print outFile "# Provide the name of the mpi build directory that should be configured,\n";
            print outFile "# compiled and tested by the test harness, or leave blank to indicate that there\n";
            print outFile "# should be no mpi build. Directory must be a subdirectory of 'Trilinos/'.\n"; 
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: NO (unless SERIAL_DIR is omitted\n";        
            print outFile "\n";
        }
        
        print outFile "MPI_DIR                         = MPI\n";
        
        if (!$short) {    
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";     
            print outFile "# (LAM only) path and filename of the file containing the names of the machines\n"; 
            print outFile "# to be used for parallel jobs. If this file doesn't exist, parallel jobs will\n";
            print outFile "# be run on the local machine only. (If this is the case, make sure\n";  
            print outFile "# MPI_STARTUP_CMD reflects this.)\n";           
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: NO\n";
            print outFile "# - the absolute path of the Trilinos directly above 'testharness' can be\n";
            print outFile "#   referred to with the value <TRILINOS_DIR>\n";          
            print outFile "\n";
        }
        
        print outFile "HOST_FILE                       = <TRILINOS_DIR>/hostfile \n";
        
        if (!$short) {    
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";     
            print outFile "# Indicate how mail should be sent. The unix program sendmail is the default.\n"; 
            print outFile "#\n";
            print outFile "# - recognized values: sendmail smtp::<mail_server> \n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES\n";      
            print outFile "\n";
        }
        
        print outFile "MAIL_METHOD                     = sendmail \n";
        
        if (!$short) {      
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";         
            print outFile "# Give the flags you would like passed to 'make'. These should be given as one\n";  
            print outFile "# string, as it will be passed verbatim.\n";         
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";  
            print outFile "# - value required: NO\n";  
            print outFile "\n";
        }
        
        print outFile "MAKE_FLAGS                      = \n";
        
        if (!$short) {      
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# Indicate whether or not the test-harness should attempt to update the CVS\n";
            print outFile "# working copy.\n";
            print outFile "#\n";
            print outFile "# - recognized values: YES NO \n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES\n";
            print outFile "\n";
        }
        
        print outFile "CVS_UPDATE                      = YES \n";
        
        if (!$short) {      
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# CVS command on this system. Note that CVS access must not require a password.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES\n";
            print outFile "\n";
        }
        
        print outFile "CVS_CMD                         = cvs \n";
        
        if (!$short) {      
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# Specify the command to start up the MPI implementation on this machine.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES (if there are parallel tests)\n";
            print outFile "# - the value of the HOST_FILE option can be referred to with the value\n";
            print outFile "#   <HOST_FILE>\n";   
            print outFile "\n";
        }
        
        print outFile "MPI_STARTUP_CMD                 = \"lamboot <HOST_FILE> -v\"\n";
        
        if (!$short) {      
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# Specify the command (if any) to shut down the MPI implementation on this\n";
            print outFile "# machine.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: NO\n";
            print outFile "\n";
        }
        
        print outFile "MPI_SHUTDOWN_CMD                = lamhalt \n";
        
        if (!$short) {        
            print outFile "\n";
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# List the email addresses to which an overall summary of this run of the\n";
            print outFile "# testharness should be sent.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: YES\n";
            print outFile "# - value required: NO\n";
            print outFile "\n";
        }
        
        print outFile "SUMMARY_EMAIL                   = \n";
        
        if (!$short) {        
            print outFile "\n";
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# List the email addresses to which test summaries and error messages will be\n";
            print outFile "# sent for test where no script owner exists.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: YES\n";
            print outFile "# - value required: NO\n";
            print outFile "\n";
        }
        
        print outFile "NO_SCRIPT_OWNER_EMAIL           = \n";
        
        if (!$short) {        
            print outFile "\n";
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# List the email addresses to which summaries and error messages will be sent\n";
            print outFile "# regardless of whether or not the script has an owner.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: YES\n";
            print outFile "# - value required: NO\n";
            print outFile "\n";
        }
        
        print outFile "ALL_SCRIPTS_EMAIL               = \n";
        
        print outFile "\n";
        print outFile "# end test-harness-config";
            
        close outFile;
    } # genConfigTemp()