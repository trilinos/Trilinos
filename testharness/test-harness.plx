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
my @optionsOrder;           # order of options (for printing in summary)

# Report variables
my %summary;                # overall test-harness summary
my $reportCount;            # count of how many reports have been generated
my %dependencies;           # package dependencies
my %emails;                 # packageDirName/regressionEmailPrefix pairs
my @codes;                  # strings corresponding to error code constants

# Error constants
my $FILE_SYSTEM_ERROR = 0;            # files absent or have incorrect permissions
my $SYSTEM_COMMAND_ERROR = 1;         # system command failed or doesn't exist
my $CONFIG_ERROR = 2;                 # test-harness configure file isn't valid
my $UPDATE_ERROR = 3;                 # cvs update failed
my $TRILINOS_CONFIGURE_ERROR = 4;     # Trilinos configure failed
my $TRILINOS_BUILD_ERROR = 5;         # Trilinos build failed
my $TEST_FAILED = 6;                  # test failed
my $TEST_PASSED = 7;                  # test passed
my $SUMMARY = 8;                      # test-harness summary
        
################################################################################
# Execution ####################################################################
################################################################################

# Preparation ==================================================================

# Prepare global variables      
prepareVariables();

# Parse command line flags
parseFlags();

# Parse dependencies
parseDependencies();

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

# Clean Up =====================================================================

# Shut down MPI implementation if any tests are parallel
mpiShutdown();

# Send summary email
report($SUMMARY);

################################################################################
# Subroutines ##################################################################
################################################################################

    ############################################################################
    # prepareVariables()
    #
    # Prepares global varibles.
    #   - global variables used: yes
    #   - sends mail: on error
    #   - args: 
    #   - returns: 

    sub prepareVariables {    
                
        $reportCount = "000";          
        @optionsOrder = ();
        
        $codes[$FILE_SYSTEM_ERROR] = "file system error";           
        $codes[$SYSTEM_COMMAND_ERROR] = "system command error";        
        $codes[$CONFIG_ERROR] = "test-harness config error";
        $codes[$UPDATE_ERROR] = "cvs update error";
        $codes[$TRILINOS_CONFIGURE_ERROR] = "Trilinos configure error";
        $codes[$TRILINOS_BUILD_ERROR] = "Trilinos build error";
        $codes[$TEST_FAILED] = "test failed";
        $codes[$TEST_PASSED] = "test passed";
        $codes[$SUMMARY] = "summary";
            
        $summary{$TRILINOS_CONFIGURE_ERROR} = ();
        $summary{$TRILINOS_BUILD_ERROR} = ();
        $summary{$TEST_FAILED} = ();
        $summary{$TEST_PASSED} = ();
        
        $emails{"amesos"} = "amesos";
        $emails{"anasazi"} = "anasazi";
        $emails{"aztecoo"} = "aztecoo";
        $emails{"belos"} = "belos";
        $emails{"epetra"} = "epetra";
        $emails{"epetraext"} = "epetra";
        $emails{"ifpack"} = "ifpack";
        $emails{"jpetra"} = "jpetra";
        $emails{"kokkos"} = "kokkos";
        $emails{"komplex"} = "komplex";
        $emails{"meros"} = "meros";
        $emails{"ml"} = "ml";
        $emails{"new_package"} = "trilinos";
        $emails{"nox"} = "nox";
        $emails{"superludist"} = "trilinos";
        $emails{"teuchos"} = "teuchos";
        $emails{"tpetra"} = "tpetra";
        $emails{"triutils"} = "triutils";
        $emails{"TSF"} = "tsf";
        $emails{"TSFCore"} = "tsf";
        $emails{"TSFCoreUtils"} = "tsf";
        $emails{"TSFExtended"} = "tsf";
        $emails{"trilinos"} = "trilinos";
        $emails{"error"} = "trilinos";
        
    } # prepareVariables()

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
        getopts("f:p:g:snth", \%flags);
        
        # parse config file and exit
        if ($flags{p}) {                         
            parseConfig($flags{p});                  
            validateOptions();        
            exit;            
        }
        
        # generate config file and exit
        if ($flags{g}) { 
            generateConfig($flags{g}, $flags{s}); 
            exit;            
        } 
        
        # fill optionsOrder array for correctly ordered output of options
        else {
            generateConfig(0, 1, 1);
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
            printEvent("cvs update\n\n");
            my $result = system $command;
            if ($result) {
                report($UPDATE_ERROR);
                die " *** error updating Trilinos - aborting test-harness ***\n";
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
        if (defined $options{'MPI_DIR'} && defined $options{'MPI_DIR'}[0]
            && defined $options{'MPI_STARTUP_CMD'} && defined $options{'MPI_STARTUP_CMD'}[0]) { 
            printEvent("mpi startup\n", "\n");
            my $command = "$options{'MPI_STARTUP_CMD'}[0]";
            my $commandFailed = system $command; 
            if ($commandFailed) {
                report($SYSTEM_COMMAND_ERROR, $command);
                printEvent($command);          
                die " *** error running system command - aborting test-harness ***\n";
            }
        }
    } # mpiStartup()
    
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
        if (defined $options{'MPI_DIR'} && defined $options{'MPI_DIR'}[0]
            && defined $options{'MPI_SHUTDOWN_CMD'} && defined $options{'MPI_SHUTDOWN_CMD'}[0]) {  
            printEvent("mpi shutdown\n", "\n");            
            my $command = "$options{'MPI_SHUTDOWN_CMD'}[0]";
            my $commandFailed = system $command; 
            if ($commandFailed) {
                report($SYSTEM_COMMAND_ERROR, $command);
                printEvent($command);          
                die " *** error running system command - aborting test-harness ***\n";
            }
        }
    } # mpiShutdown()

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
        
        if (!$flags{t}) {
            # delete and recreate mpi dir
            if (defined $options{'MPI_DIR'} && defined $options{'MPI_DIR'}[0]) {
                system "rm -rf $options{'MPI_DIR'}[0]";
                system "mkdir $options{'MPI_DIR'}[0]";
            }
            
            # delete and recreate serial dir
            if (defined $options{'SERIAL_DIR'} && defined $options{'SERIAL_DIR'}[0]) {
                system "rm -rf $options{'SERIAL_DIR'}[0]";
                system "mkdir $options{'SERIAL_DIR'}[0]";
            }
        }
        
    } # prepareBuildDirs()

    ############################################################################
    # run()
    #
    # Core test-harness logic. Calls configure(), build(), test(), and report()
    #   - global variables used: yes
    #   - sends mail: yes
    #   - args: 
    #   - returns: 

    sub run {  
            
        prepareInvokeConfigure();
        
        # create array of build dirs to iterate upon ===========================            
        my @buildDir = ();
        if (defined $options{'MPI_DIR'} && defined $options{'MPI_DIR'}[0]) {
            push (@buildDir, $options{'MPI_DIR'}[0]);
        }
        if (defined $options{'SERIAL_DIR'} && defined $options{'SERIAL_DIR'}[0]) {
            push (@buildDir, $options{'SERIAL_DIR'}[0]);
        }
        
        # for each build directory =============================================
        for (my $j=0; $j<=$#buildDir; $j++) {
        
            my $comm; 
            if (defined $options{'MPI_DIR'}[0] && $buildDir[$j] eq $options{'MPI_DIR'}[0]) {
                $comm = "mpi";
            } elsif (defined $options{'SERIAL_DIR'}[0] && $buildDir[$j] eq $options{'SERIAL_DIR'}[0]) {
                $comm = "serial";
            } 

            # -t (test only) flag is absent...
            if (!$flags{t}) {
                
                # descend into build dir
                chdir"$options{'TRILINOS_DIR'}[0]/$buildDir[$j]";
                
                # test for existence of invoke-configure
                if (!-f "invoke-configure") {
                    my $message = "";
                    $message .= "$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/";
                    $message .= "invoke-configure must be present\n";
                    report($FILE_SYSTEM_ERROR, $message, $comm);
                    printEvent($message);
                    die " *** file missing - aborting test-harness ***\n";
                } 
                
                # test for read permissions on invoke-configure
                if (!-r "invoke-configure") {
                    my $message = "";
                    $message .= "$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/";
                    $message .= "invoke-configure must be readable\n";
                    report($FILE_SYSTEM_ERROR, $message, $comm);
                    printEvent($message);
                    die " *** file permission wrong - aborting test-harness ***\n";
                } 
                
                # test for executable permission on invoke-configure
                if (!-x "invoke-configure") {
                    my $message = "";
                    $message .= "$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/";
                    $message .= "invoke-configure must be executable\n";
                    report($FILE_SYSTEM_ERROR, $message, $comm);
                    printEvent($message);
                    die " *** file permission wrong - aborting test-harness ***\n";
                }
                
                # should we quit trying? (no invoke-configure left)
                my $quitTrying = 0;
                
                # configure --------------------------------------------------------
                
                # configure passed?
                my $configurePassed = 0;            
                
                while (!$configurePassed && !$quitTrying) {
                
                    # attempt to configure
                    my $configureOutput = configure($buildDir[$j]);
                    
                    # configure failed
                    if ($configureOutput) {    
                        
                        # fix invoke configure
                        my $log = "$options{'TRILINOS_DIR'}[0]/testharness/temp/trilinos_configure_log.txt";
                        my $invokeConfigure ="$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/invoke-configure";                    
                        (my $brokenPackage, my $brokenPackageDir) = fixInvokeConfigure($log, $invokeConfigure);      
                            
                        # quit if error fixing invoke-configure
                        if ($brokenPackage eq "error") {
                            $quitTrying = 1; 
                            report($TRILINOS_BUILD_ERROR, $brokenPackage, $comm); 
                            last; # equivalent to break                          
                        }       
                        
                        # remove broken package
                        system "rm -rf $options{'TRILINOS_DIR'}[0]/$buildDir[$j]/packages/$brokenPackageDir";
                        
                        printEvent("$comm - Trilinos configuration failed--$brokenPackage broke\n");  
                        
                        # there is no invoke-configure left or it's empty
                        if (!-f $invokeConfigure || -z $invokeConfigure) { 
                            printEvent("$comm - invoke-configure missing or empty\n");   
                            $quitTrying = 1;
                        }
                        
                        # send email
                        report($TRILINOS_CONFIGURE_ERROR, $brokenPackage, $comm);
                    } 
                    
                    # configure succeeded
                    else {     
                        $configurePassed = 1;
                        printEvent("$comm - Trilinos configured successfully\n");                      
                    }
                
                    # build --------------------------------------------------------
                    
                    # build passed?            
                    my $buildPassed = 0;            
                    
                    while ($configurePassed && !$buildPassed && !$quitTrying) {
                    
                        # attempt to build
                        my $buildOutput = build($buildDir[$j]);
                        
                        # build failed
                        if ($buildOutput) {
                                       
                            # fix invoke configure         
                            my $log = "$options{'TRILINOS_DIR'}[0]/testharness/temp/trilinos_build_log.txt";
                            my $invokeConfigure ="$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/invoke-configure";                    
                            (my $brokenPackage, my $brokenPackageDir) = fixInvokeConfigure($log, $invokeConfigure);
                            
                            # quit if error fixing invoke-configure
                            if ($brokenPackage eq "error") {
                                $quitTrying = 1; 
                                report($TRILINOS_BUILD_ERROR, $brokenPackage, $comm); 
                                last; # equivalent to break                               
                            }
                            
                            # remove broken package
                            system "rm -rf $options{'TRILINOS_DIR'}[0]/$buildDir[$j]/packages/$brokenPackageDir";
                            
                            printEvent("$comm - Trilinos build failed--$brokenPackage broke\n");  
                        
                            # there is no invoke-configure left or it's empty
                            if (!-f $invokeConfigure || -z $invokeConfigure) {  
                                printEvent("$comm - invoke-configure missing or empty\n"); 
                                $quitTrying = 1;
                            }
                            
                            # force reconfigure
                            $configurePassed = 0;
                            
                            # send email
                            report($TRILINOS_BUILD_ERROR, $brokenPackage, $comm); 
                        } 
                        
                        # build succeeded
                        else {    
                            $buildPassed = 1;
                            printEvent("$comm - Trilinos built successfully\n");                  
                        }
                        
                    } # while ($configurePassed && !$buildPassed && !$quitTrying)
                    
                } # while (!$configurePassed && !$quitTrying)                       
            } # if (!-t)
            
            # test -------------------------------------------------------------
            
            # if -n flag absent -- run tests
            if (!$flags{n}) {
            
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
                        my $potentialTestDir = "";
                        $potentialTestDir .= "$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/";
                        $potentialTestDir .= "$testDir/scripts/$options{'FREQUENCY'}[0]/$comm";
                        
                        if (-d $potentialTestDir) {
                            my $command = "find $potentialTestDir -type f";
                            my $output = `$command`;
                            my @potentialScripts = split (/\s+/, $output);
                                
                            # run each test
                            foreach my $potentialScript (@potentialScripts) {
                                $potentialScript =~ s/\s*$//;  # trim trailing whitespace
                                
                                # if potential script file is executable...
                                if (-x $potentialScript) {
                    
                                    # test
                                    my $testFailed = test($buildDir[$j], $potentialScript);
                                    
                                    # extract test name from path (for printing)
                                    my $testNameOnly = $potentialScript;
                                    $testNameOnly =~ s/.*\///; 
                                    
                                    if ($testFailed) {                                               
                                        printEvent("$testNameOnly - Test failed\n");  
                                        report($TEST_FAILED, $testFailed, $comm, $testDir, $potentialScript); 
                                    } else {                                    
                                        printEvent("$testNameOnly - Test passed\n");  
                                        report($TEST_PASSED, $testFailed, $comm, $testDir, $potentialScript);
                                    }
                                    
                                } # if (executable)
                            } # foreach ($potentialScript)  
                        } # if (-f $potentialTestDir)
                                          
                    } # unless (unsupported)
                } # foreach $testDir
            } # if (!-n)
            
            printEvent("\n");
        } # for (buildDirs)                   
    } # run()

    ############################################################################
    # prepareInvokeConfigure()
    #
    # Compile the invoke configure files and move them into place
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: 
    #   - returns: 

    sub prepareInvokeConfigure {   
        	
        # create complete invoke-configure files from elements =================              
        my $rawInvokeConfigureMPI = "";               
        my $rawInvokeConfigureSERIAL = "";  
        
        if (defined $options{'MPI_DIR'} && defined $options{'MPI_DIR'}[0]) {	    
	        # prepare machine-dependent-mpi portion of invoke-configure
	        if (-f "$options{'TRILINOS_DIR'}[0]/testharness/elements-machine/$options{'MACHINE_MPI_CONFIG_FILE'}[0]") {
                open (MACHINE_MPI, "<$options{'TRILINOS_DIR'}[0]/testharness/elements-machine/$options{'MACHINE_MPI_CONFIG_FILE'}[0]") 
                    or die "$! error trying to open file";  
                while (my $line = <MACHINE_MPI>) {
                    $line =~ s/^\s*//;  # trim leading spaces
                    $line =~ s/\s*$//;  # trim trailing spaces                
                    $rawInvokeConfigureMPI .= "$line \\\n";              
                }
        	    close MACHINE_MPI;  
            } else {
                my $message = "";
                $message .= "$options{'TRILINOS_DIR'}[0]/testharness/elements-machine";
                $message .= "/$options{'MACHINE_MPI_CONFIG_FILE'}[0] does not exist\n";
                report($FILE_SYSTEM_ERROR, $message);
                printEvent($message);
                die " *** file missing - aborting test-harness ***\n";
            }   
        }
        
        # prepare machine-dependent portion of invoke-configure
        if (-f "$options{'TRILINOS_DIR'}[0]/testharness/elements-machine/$options{'MACHINE_CONFIG_FILE'}[0]") {
            open (MACHINE, "<$options{'TRILINOS_DIR'}[0]/testharness/elements-machine/$options{'MACHINE_CONFIG_FILE'}[0]") 
                or die "$! error trying to open file";                  
            while (my $line = <MACHINE>) {
                $line =~ s/^\s*//;  # trim leading spaces
                $line =~ s/\s*$//;  # trim trailing spaces                
                $rawInvokeConfigureMPI .= "$line \\\n";
                $rawInvokeConfigureSERIAL .= "$line \\\n";               
            }
    	    close MACHINE;                
        } else {
            my $message = "";
            $message .= "$options{'TRILINOS_DIR'}[0]/testharness/elements-machine";
            $message .= "/$options{'MACHINE_CONFIG_FILE'}[0] does not exist\n";
            report($FILE_SYSTEM_ERROR, $message);
            printEvent($message);
            die " *** file missing - aborting test-harness ***\n";
        }   
        
        # prepare machine-independent portion of invoke-configure
        if (-f "$options{'TRILINOS_DIR'}[0]/testharness/elements-trilinos/$options{'TRILINOS_CONFIG_FILE'}[0]") {
            open (COMPONENT, "<$options{'TRILINOS_DIR'}[0]/testharness/elements-trilinos/$options{'TRILINOS_CONFIG_FILE'}[0]") 
                or die "$! error trying to open file";  
            while (my $line = <COMPONENT>) {
                $line =~ s/^\s*//;  # trim leading spaces
                $line =~ s/\s*$//;  # trim trailing spaces                
                $rawInvokeConfigureMPI .= "$line \\\n";
                $rawInvokeConfigureSERIAL .= "$line \\\n";               
            }
    	    close COMPONENT;      
        } else {
            my $message = "";
            $message .= "$options{'TRILINOS_DIR'}[0]/testharness/elements-machine";
            $message .= "/$options{'TRILINOS_CONFIG_FILE'}[0] does not exist\n";
            report($FILE_SYSTEM_ERROR, $message);
            printEvent($message);
            die " *** file missing - aborting test-harness ***\n";
        }        
	    
	    # clean up raw invoke-configures
	    $rawInvokeConfigureMPI =~ s/\\\n$//;        # remove last line-continuation
	    $rawInvokeConfigureSERIAL =~ s/\\\n$//;     # remove last line-continuation
	                
	    # create and copy MPI invoke configure
	    if (defined $options{'MPI_DIR'} && defined $options{'MPI_DIR'}[0]) {
	    
            # open INVOKE_CONFIGURE_MPI for writing
            open (INVOKE_CONFIGURE_MPI, ">$options{'TRILINOS_DIR'}[0]/testharness/temp/invoke-configure-mpi")
                or die "$! error trying to open file";
            print INVOKE_CONFIGURE_MPI ".././configure ";
            print INVOKE_CONFIGURE_MPI $rawInvokeConfigureMPI;
            close INVOKE_CONFIGURE_MPI;
            
            # move invoke-configure file into place
            my $command;                     
            $command = "cp $options{'TRILINOS_DIR'}[0]/testharness/temp/invoke-configure-mpi ";
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
            print INVOKE_CONFIGURE_SERIAL $rawInvokeConfigureSERIAL;
            close INVOKE_CONFIGURE_SERIAL;
            
            # move invoke-configure file into place
            my $command;                     
            $command = "cp $options{'TRILINOS_DIR'}[0]/testharness/temp/invoke-configure-serial ";
            $command .= "$options{'TRILINOS_DIR'}[0]/$options{'SERIAL_DIR'}[0]/invoke-configure";
            system $command;
            
            # set invoke-configure permissions
            system "chmod a+rx $options{'TRILINOS_DIR'}[0]/$options{'SERIAL_DIR'}[0]/invoke-configure";
        }
        
    } # prepareInvokeConfigure()

    ############################################################################
    # fixInvokeConfigure()
    #
    # Fix the invoke configure script based on the broken package indicated 
    # in the log file.
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: $log $invokeConfigure
    #   - returns: 

    sub fixInvokeConfigure {   
        my $log = $_[0];    
        my $invokeConfigure = $_[1];
        	
        # open configure/build log for reading
        open (LOG, "<$log")
            or die "can't open $log";
            
        # read configure/build log and extract last instance
        undef $/;           # undefine input record separator
        my $file=<LOG>;     # copy entire file
        $/ = "\n";          # restore it to default newline
	    close LOG;
	    
	    # extract and fix up name of broken package
	    my $brokenPackage;
	    if ($log =~ m/configure/) {
    	    $file =~ m/.*^Running(.*?)Configure Script/ms;
	        if (defined $1) { $brokenPackage = $1; }
    	} elsif ($log =~ m/build/) {    
    	    $file =~ m/.*\/packages\/(\w*)\W/ms;
	        if (defined $1) { $brokenPackage = $1; }
        }
	    
	    if (defined $brokenPackage) {
            $brokenPackage =~ s/^\s*//;             # trim leading spaces
            $brokenPackage =~ s/\s*$//;             # trim trailing spaces
    	    $brokenPackage = lc($brokenPackage);    # convert to lowercase
            $brokenPackage =~ s/ /_/g;              # convert spaces to underscores
        } else {
            printEvent("error fixing invoke-configure--can't detect package\n");
            return ("error", "error");
        }
        
        system "cp $invokeConfigure $invokeConfigure-broken";
        
        # open invoke-configure for reading
        open (INVOKE_CONFIGURE, "<$invokeConfigure")
            or die "can't open $invokeConfigure";
            
        # parse it, extracting lines to keep
        my $newInvokeConfigure = "";
        my $changeMade = 0;
        while (my $line = <INVOKE_CONFIGURE>) {
            $line =~ s/^\s*//;  # trim leading spaces
            $line =~ s/\s*$//;  # trim trailing spaces
            
            # compare lines to package dependencies
            my $dropLine = 0;
            my $lastElementIndex = $#{$dependencies{$brokenPackage}};
            for (my $i=0; $i<=$lastElementIndex; $i++) { 
                if ($line =~ m/$dependencies{$brokenPackage}[$i]\b/i) {
                    $dropLine = 1;   
                }   
            }    
            
            # write line if it isn't a dependency
            if (!$dropLine) {                                        
                $newInvokeConfigure .= "$line\n";
            } else {
                $changeMade = 1;
            }
        } # while ($line)
        close INVOKE_CONFIGURE;
        
        if (!$changeMade) {
            printEvent("error fixing invoke-configure--no changes made ($brokenPackage broke)\n");
            return ("error", "error");
        }
        
        # prepare new invoke-configure
        $newInvokeConfigure =~ s/(.*)\\$/$1/s;  # remove last line-continuation
        chomp ($newInvokeConfigure);            # remove last newline 
        
        # open invoke-configure for writing
        open (INVOKE_CONFIGURE, ">$invokeConfigure")
            or die "can't open $invokeConfigure";
        
        # write new invoke configure
        print INVOKE_CONFIGURE $newInvokeConfigure;
        close INVOKE_CONFIGURE;      
                        
        # create package name corresponding to directory case...grrr...
        my $brokenPackageDir = $brokenPackage;
        if ($brokenPackageDir eq "tsf") { $brokenPackageDir = "TSF" }
        if ($brokenPackageDir eq "tsfcore") { $brokenPackageDir = "TSFCore" }
        if ($brokenPackageDir eq "tsfcoreutils") { $brokenPackageDir = "TSFCoreUtils" }
        if ($brokenPackageDir eq "tsfextended") { $brokenPackageDir = "TSFExtended" }
        
        # return broken package name
        return ($brokenPackage, $brokenPackageDir);
        
    } # fixInvokeConfigure()

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
        	        
        chdir"$options{'TRILINOS_DIR'}[0]/$buildDir";    
        	
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
        	        
        chdir"$options{'TRILINOS_DIR'}[0]/$buildDir";     
    
        my $command = "";
        if (defined $options{'MAKE_FLAGS'} && defined $options{'MAKE_FLAGS'}[0]) {
            $command .= "make $options{'MAKE_FLAGS'}[0] >> $options{'TRILINOS_DIR'}[0]";
            $command .= "/testharness/temp/trilinos_build_log.txt 2>&1";
        } else {
            $command .= "make >> $options{'TRILINOS_DIR'}[0]";
            $command .= "/testharness/temp/trilinos_build_log.txt 2>&1";
        }
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
    # report()
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

    sub report {
        my $code = $_[0];     
        my $message = $_[1];      
        my $comm = $_[2];  
        my $testDir = $_[3];
        my $testName = $_[4];
        $reportCount++;
        
        # preparation ##########################################################
        
        chdir "$options{'TRILINOS_DIR'}[0]/testharness/";
        
        # link to MIME::Lite module
        if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
            use lib "./lib/MIME-Lite-3.01/lib";     # cvs checkin MIME::Lite
            use MIME::Lite;                         # might need Net::SMTP
        }
        
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
        if ($code == $TRILINOS_CONFIGURE_ERROR || $code == $TRILINOS_BUILD_ERROR) {
            push (@{$summary{$code}}, "$comm ($message broke)");
        }        
        if ($code == $TEST_FAILED || $code == $TEST_PASSED) {
            my $package;                
    	    $testDir =~ m/packages\/(.*?)\//;
	        if ($1) { 
	            $package = $1;
                push (@{$summary{$code}}, "$comm - $package - $testName");
            } else {
                push (@{$summary{$code}}, "$comm - $testName");
            }
        }

        # compile list of mail recipients --------------------------------------        
        my $mailTo = "";
        if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
        
            # detect default recipients (if enabled)
            if ($options{'SEND_TO_DEFAULTS'}[0] eq "YES") {
                
                # configure/build-related
                if ($code == $TRILINOS_CONFIGURE_ERROR || $code == $TRILINOS_BUILD_ERROR) {
                    my $package = $message;
                    my $packageRegression = "";   
                    if ($emails{$package}) {
    	                $packageRegression = $emails{$package};
        	            $packageRegression .= "-regression\@software.sandia.gov";
                        $mailTo .= "$packageRegression, ";
                    } else {
            	        printEvent("error - unable to detect package regression list\n");
                	    printEvent("  sending to trilinos-regression instead\n");
                        $mailTo .= "trilinos-regression\@software.sandia.gov, ";
            	    }
                } # if (configure/build-related)
            
                # test-related
                if ($code == $TEST_FAILED || $code == $TEST_PASSED) {
            
                    # extract $scriptOwner email from first line of log$hostOS.txt
                    # (note: the name "log$hostOS.txt" is functional--it is written to 
                    # by the test scripts.
                    my $scriptOwner = "";
                    if (-f "$options{'TRILINOS_DIR'}[0]/log$hostOS.txt") {
                        open 
                        (OWNER, "$options{'TRILINOS_DIR'}[0]/log$hostOS.txt") 
                            or die "$! error trying to open file";
                        $scriptOwner=<OWNER>;   # read first line (email of script owner is on first line)
                        chomp $scriptOwner;     # trim newline
                        close OWNER;
                    }        
                    
                    # script-owner found--append to recipients
                    if ($scriptOwner =~ m/\S+?\@\S+?/i) {    # look for ...@...
                        $mailTo .= "$scriptOwner, ";
                    } 
                    
                    # script-owner not found, use package regression
                    else {    
                        my $error = 0;
                        my $package = "";
                        my $packageRegression = "";                  
                	    $testDir =~ m/packages\/(.*?)\//;
            	        if ($1) { 
            	            $package = $1; 
            	        } else { $error = 1; }
            	        if ($package) { 
            	            if ($emails{$package}) {
            	                $packageRegression = $emails{$package};
                	            $packageRegression .= "-regression\@software.sandia.gov";
                                $mailTo .= "$packageRegression, ";
                            } else { $error = 1; }
                	    } else { $error = 1; }
                	    if ($error) {
                	        printEvent("error - unable to detect package regression list\n");
                	        printEvent("  sending to trilinos-regression instead\n");
                            $mailTo .= "trilinos-regression\@software.sandia.gov, ";
                	    }
                    }            
                } # if (test-related)            
                
                # summary
                if ($code == $SUMMARY) {
                    $mailTo .= "trilinos-regression\@software.sandia.gov, ";
                } 
                
            } # if (SEND_TO_DEFAULTS)
           
            # append ALL_EMAILS       
            if (defined $options{'ALL_EMAILS'} && defined $options{'ALL_EMAILS'}[0]) { 
                my $lastElementIndex = $#{$options{'ALL_EMAILS'}};
                for (my $i=0; $i<=$lastElementIndex; $i++) {
                    if ($mailTo !~ m/$options{'ALL_EMAILS'}[$i]/i) {
                        $mailTo .= "$options{'ALL_EMAILS'}[$i], ";
                    } 
                }
            }
            
            # append SUMMARY_EMAIL
            if ($code == $SUMMARY) { 
                if (defined $options{'SUMMARY_EMAIL'} && defined $options{'SUMMARY_EMAIL'}[0]) {
                    my $lastElementIndex = $#{$options{'SUMMARY_EMAIL'}};
                    for (my $i=0; $i<=$lastElementIndex; $i++) {
                        if ($mailTo !~ m/$options{'SUMMARY_EMAIL'}[$i]/i) {
                            $mailTo .= "$options{'SUMMARY_EMAIL'}[$i], ";
                        }
                    }
                }
            } 
            
            # clean up list of recipients    
            $mailTo =~ s/^\s*//;        # trim leading spaces
            $mailTo =~ s/\s*$//;        # trim trailing spaces
            $mailTo =~ s/,*$//;         # trim trailing comma  
        
        } # compile list of mail recipients (if EMAIL)

        # Create email/report ##################################################

        # subject/filename =====================================================
        
        my $subject = "";
        my $filename;
        
        # subject (email)
        if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
            $subject .= "$hostOS - $hostName - ";
            if (defined $comm) {$subject .= "$comm - ";}
            if (defined $testName) {$subject .= "$testName - ";}
            $subject .= $codes[$code];
        }    
        
        # filename (local report)
        elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") {
            if ($options{'REPORT_NAMES'}[0] eq "ORDER") {
                if (defined $comm) {    
                    $filename = "$reportCount"."-$comm-$codes[$code]";
                } else {
                    $filename = "$reportCount"."-$codes[$code]";
                }
            } elsif ($options{'REPORT_NAMES'}[0] eq "EVENT") {   
                if (defined $comm) {         
                    $filename = "$comm-$codes[$code]";
                } else {
                    $filename = "$codes[$code]";
                }
            }            
        
            if ($code == $TRILINOS_CONFIGURE_ERROR || $code == $TRILINOS_BUILD_ERROR) {
                $filename .= "($message)";
            }
            
            if ($code == $TEST_FAILED || $code == $TEST_PASSED) {
                my $package = "";              
        	    $testDir =~ m/packages\/(.*?)\//;
    	        if ($1) { $package = $1; }
    	        $filename .= "($package-$testName)";
            }
            
            $filename =~ s/ /_/g;       # convert spaces to underscores
            $filename =~ s/\n//g;       # remove any newlines
            $filename = lc($filename);  # convert to lowercase
        }
        
        # construct email/report ===============================================    
        
        # initialize email object
        my $email;
        if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {   
            $email=MIME::Lite->new(
    			From =>     'Trilinos test-harness <trilinos-regression@software.sandia.gov>',
    			To =>       $mailTo,
    			Subject =>  $subject,
    			Type =>     'multipart/mixed'
    		);
		}
            
        # body =================================================================
        
        my $body = "";
                
        # header section
        if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
            $body .= "Host OS:          $hostOS\n";}
        #if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
        #    $body .= "Host OS Release:  $hostOSRelease\n";}    # unnecessary info
        #if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
        #    $body .= "Host OS Version:  $hostOSVersion\n";}    # unnecessary info
        #if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
        #    $body .= "Host Hardware:    $hostHardware\n";}     # unnecessary info
        if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
            $body .= "Host Name:        $hostName\n";}
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
        if ($code == $TRILINOS_CONFIGURE_ERROR || $code == $TRILINOS_BUILD_ERROR) {
            $body .= "Result:           $codes[$code] - $message broke\n";
            $body .= "\n";   
        } elsif ($code != $SUMMARY) {
            $body .= "Result:           $codes[$code]\n";
            $body .= "\n";   
        }             
        
        # fatal error ----------------------------------------------------------
        if ($code == $FILE_SYSTEM_ERROR || $code == $SYSTEM_COMMAND_ERROR ||
            $code == $CONFIG_ERROR || $code == $UPDATE_ERROR) {
        
            $body .= "------------------------------------------------------------\n";
            $body .= "FATAL ERROR\n";
            $body .= "\n";       
            $body .= "This error caused the test-harness to quit prematurely. It is\n";
            $body .= "the sort of error that the test-harness is not designed to\n";
            $body .= "recover from. It is probably either a simple human oversight\n";
            $body .= "or, once fixed, should not occur again on this machine.\n";    
            $body .= "\n";       
            $body .= "Message/exit-status:      $message\n";           
            $body .= "\n";        
        }
        
        # print the summary hash of arrays -------------------------------------        
        if ($code == $SUMMARY) {
            $body .= "------------------------------------------------------------\n";
            $body .= "Summary: \n";
            $body .= "\n";        
            
            for my $id (reverse sort keys %summary) { 
                my $lastElementIndex = $#{$summary{$id}}; 
                $body .= "- $codes[$id] (".($lastElementIndex+1)."): \n";
                $body .= "\n";
                if ($lastElementIndex+1 > 0) {
                    for (my $i=0; $i<=$lastElementIndex; $i++) {
                        $body .= "    $summary{$id}[$i]\n";
                    }
                    $body .= "\n";
                } 
            } 
        } # summary    
        
        # print the test-harness-config hash of arrays -------------------------        
        if ($code == $SUMMARY) {
            $body .= "------------------------------------------------------------\n";
            $body .= "Test-harness Config: \n";
            $body .= "\n";        
            
            for my $key (@optionsOrder) {
                if ($key ne "TRILINOS_DIR" && $key ne "FREQUENCY") {
                    my $lastElementIndex = $#{$options{$key}}; 
                    $body .= "$key (".($lastElementIndex+1)."): ";
                    for (my $j=0; $j<=$lastElementIndex; $j++) {
                        if ($j != 0) { $body .= ", "; }
                        $body .= "$options{$key}[$j]";
                    }
                    $body .= "\n";
                }
            } 
            $body .= "\n";
        } # summary
        
        # print the dependencies hash of arrays --------------------------------        
        if ($code == $SUMMARY) {
            $body .= "------------------------------------------------------------\n";
            $body .= "Package Dependencies: \n";
            $body .= "\n";        
            
            for my $key (sort keys %dependencies) { 
                my $lastElementIndex = $#{$dependencies{$key}}; 
                $body .= "$key (".($lastElementIndex+1)."): \n";
                for (my $i=0; $i<=$lastElementIndex; $i++) {
                    $body .= "    $dependencies{$key}[$i]\n";
                }
                $body .= "\n";
            } 
        } # summary
        
        # attachments ----------------------------------------------------------
        
        my $attachmentText = "";
        my $attachmentsExist = 0;                    
        
        # update failed
        if ($code == $UPDATE_ERROR && -f "update_log.txt") {
            $attachmentsExist = 1;
            my $log = "update_log.txt";     
            my $logPath = "$options{'TRILINOS_DIR'}[0]/testharness/temp/$log";       
            if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                $attachmentText .= "    $log\n";
                $email->attach(Type=>'TEXT', Path=>"$logPath", Disposition=>'attachment');
            } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                $attachmentText .= appendFile($log, $logPath);                
            }
        }
        
        # trilinos configure failed
        if ($code == $TRILINOS_CONFIGURE_ERROR && -f "trilinos_configure_log.txt") {
            $attachmentsExist = 1;
            my $log = "trilinos_configure_log.txt";
            my $logPath = "$options{'TRILINOS_DIR'}[0]/testharness/temp/$log";       
            if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                $attachmentText .= "    $log\n";
                $email->attach(Type=>'TEXT', Path=>"$logPath", Disposition=>'attachment');
            } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                $attachmentText .= appendFile($log, $logPath);                
            }
        }       
        
        # trilinos build failed
        if ($code == $TRILINOS_BUILD_ERROR && -f "trilinos_build_log.txt") {
            $attachmentsExist = 1;
            my $log = "trilinos_build_log.txt";   
            my $logPath = "$options{'TRILINOS_DIR'}[0]/testharness/temp/$log";         
            if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                $attachmentText .= "    $log\n";
                $email->attach(Type=>'TEXT', Path=>"$logPath", Disposition=>'attachment');
            } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") {
                $attachmentText .= appendFile($log, $logPath);                
            }
        }
        
        # test compile log
        if (-f "test_compile_log.txt") {
            $attachmentsExist = 1;
            my $log = "test_compile_log.txt";     
            my $logPath = "$options{'TRILINOS_DIR'}[0]/testharness/temp/$log";       
            if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                $attachmentText .= "    $log\n";
                $email->attach(Type=>'TEXT', Path=>"$logPath", Disposition=>'attachment');
            } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                $attachmentText .= appendFile($log, $logPath);                
            }
        }
        
        # parallel test failed
        if (-f "$options{'TRILINOS_DIR'}[0]/logMpiErrors.txt") {   
            $attachmentsExist = 1;
            my $log = "logMpiErrors.txt";     
            my $logPath = "$options{'TRILINOS_DIR'}[0]/$log";
            if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                $attachmentText .= "    $log\n";
                $email->attach(Type=>'TEXT', Path=>"$logPath", Disposition=>'attachment');
            } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                $attachmentText .= appendFile($log, $logPath);                       
            }
        }
        
        # serial test failed
        if (-f "$options{'TRILINOS_DIR'}[0]/logErrors.txt") {
            $attachmentsExist = 1;
            my $log = "logErrors.txt";     
            my $logPath = "$options{'TRILINOS_DIR'}[0]/$log";       
            if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                $attachmentText .= "    $log\n";
                $email->attach(Type=>'TEXT', Path=>"$logPath", Disposition=>'attachment');
            } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                $attachmentText .= appendFile($log, $logPath);                
            }
        }
        
        # test output
        if (-f "$options{'TRILINOS_DIR'}[0]/log$hostOS.txt") {
            $attachmentsExist = 1;
            my $log = "log$hostOS.txt";     
            my $logPath = "$options{'TRILINOS_DIR'}[0]/$log";       
            if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                $attachmentText .= "    $log\n";
                $email->attach(Type =>'TEXT', Path=>"$logPath", Disposition=>'attachment');
            } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                $attachmentText .= appendFile($log, $logPath);                
            }
        }      
        
        # event log
        if ($code == $SUMMARY && -f "event_log.txt") {
            $attachmentsExist = 1;
            my $log = "event_log.txt";     
            my $logPath = "$options{'TRILINOS_DIR'}[0]/testharness/temp/$log";       
            if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                $attachmentText .= "    $log\n";
                $email->attach(Type=>'TEXT', Path=>"$logPath", Disposition=>'attachment');
            } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                $attachmentText .= appendFile($log, $logPath);                
            }
        }
        
        # invoke-configure attachments -----------------------------------------
        
        # invoke-configure attachments for the summary email
        if ($code == $SUMMARY) {
        
            # mpi invoke-configures
            if (defined $options{'MPI_DIR'}[0]) {
                my $buildPath = "$options{'TRILINOS_DIR'}[0]/$options{'MPI_DIR'}[0]";
            
                # at least one configure/build error
                if (-f "$buildPath/invoke-configure-broken" 
                    && ! -z "$buildPath/invoke-configure-broken") {
                    
                    $attachmentsExist = 1;
                    
                    system "mv invoke-configure-mpi invoke-configure-mpi-original";                    
                    my $log = "invoke-configure-mpi-original";     
                    my $logPath = "$options{'TRILINOS_DIR'}[0]/testharness/temp/$log";       
                    if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                        $attachmentText .= "    $log\n";
                        $email->attach(Type =>'TEXT', Path=>"$logPath", Disposition=>'attachment');
                    } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                        $attachmentText .= appendFile($log, $logPath);                
                    }                    
                    
                    system "cp $buildPath/invoke-configure invoke-configure-mpi-final";
                    $log = "invoke-configure-mpi-final";     
                    $logPath = "$options{'TRILINOS_DIR'}[0]/testharness/temp/$log";       
                    if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                        $attachmentText .= "    $log\n";
                        $email->attach(Type =>'TEXT', Path=>"$logPath", Disposition=>'attachment');
                    } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                        $attachmentText .= appendFile($log, $logPath);                
                    }
                } 
                
                # no configure/build errors
                else {
                    $attachmentsExist = 1;
                    my $log = "invoke-configure-mpi";     
                    my $logPath = "$options{'TRILINOS_DIR'}[0]/testharness/temp/$log";       
                    if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                        $attachmentText .= "    $log\n";
                        $email->attach(Type =>'TEXT', Path=>"$logPath", Disposition=>'attachment');
                    } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                        $attachmentText .= appendFile($log, $logPath);                
                    }               
                }
            }
            
            # serial invoke-configures
            if (defined $options{'SERIAL_DIR'}[0]) {
                my $buildPath = "$options{'TRILINOS_DIR'}[0]/$options{'SERIAL_DIR'}[0]";
            
                # at least one configure/build error
                if (-f "$buildPath/invoke-configure-broken" 
                    && ! -z "$buildPath/invoke-configure-broken") {
                    $attachmentsExist = 1;
                                        
                    system "mv invoke-configure-serial invoke-configure-serial-original";
                    my $log = "invoke-configure-serial-original";     
                    my $logPath = "$options{'TRILINOS_DIR'}[0]/testharness/temp/$log";       
                    if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                        $attachmentText .= "    $log\n";
                        $email->attach(Type =>'TEXT', Path=>"$logPath", Disposition=>'attachment');
                    } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                        $attachmentText .= appendFile($log, $logPath);                
                    }
                    
                    system "cp $buildPath/invoke-configure invoke-configure-serial-final";
                    $log = "invoke-configure-serial-final";     
                    $logPath = "$options{'TRILINOS_DIR'}[0]/testharness/temp/$log";       
                    if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                        $attachmentText .= "    $log\n";
                        $email->attach(Type =>'TEXT', Path=>"$logPath", Disposition=>'attachment');
                    } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                        $attachmentText .= appendFile($log, $logPath);                
                    }
                } 
                
                # no configure/build errors
                else {
                    $attachmentsExist = 1;
                    my $log = "invoke-configure-serial";     
                    my $logPath = "$options{'TRILINOS_DIR'}[0]/testharness/temp/$log";       
                    if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                        $attachmentText .= "    $log\n";
                        $email->attach(Type =>'TEXT', Path=>"$logPath", Disposition=>'attachment');
                    } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                        $attachmentText .= appendFile($log, $logPath);                
                    }                 
                }
            }
        } # invoke-configure attachments for the summary email
        
        # Trilinos configure/build or test compile/pass/fail email
        if ($code == $TRILINOS_CONFIGURE_ERROR || $code == $TRILINOS_BUILD_ERROR 
            || $code == $TEST_FAILED || $code == $TEST_PASSED) {
            
            # build directory
            my $buildDir = "";
            if ($comm eq "mpi") { $buildDir = $options{'MPI_DIR'}[0]; }
            if ($comm eq "serial") { $buildDir = $options{'SERIAL_DIR'}[0]; }
            my $invokeConfigure = "$options{'TRILINOS_DIR'}[0]/$buildDir/invoke-configure";
            
            # configure/build failed--attach broken invoke-configure
            if (($code == $TRILINOS_CONFIGURE_ERROR || $code == $TRILINOS_BUILD_ERROR)
                && -f $invokeConfigure."-broken" && ! -z $invokeConfigure."-broken") {
                $attachmentsExist = 1;
                my $log = "invoke-configure-broken";     
                my $logPath = "$invokeConfigure-broken";       
                if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                    $attachmentText .= "    $log\n";
                    $email->attach(Type =>'TEXT', Path=>"$logPath", Disposition=>'attachment');
                } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                    $attachmentText .= appendFile($log, $logPath);                
                }   
            } 
            
            elsif (($code == $TEST_FAILED || $code == $TEST_PASSED)
                && -f $invokeConfigure && ! -z $invokeConfigure) {
                $attachmentsExist = 1;
                my $log = "invoke-configure";     
                my $logPath = "$invokeConfigure";       
                if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {
                    $attachmentText .= "    $log\n";
                    $email->attach(Type =>'TEXT', Path=>"$logPath", Disposition=>'attachment');
                } elsif ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
                    $attachmentText .= appendFile($log, $logPath);                
                }   
            }
        } # Trilinos configure/build or test compile/pass/fail email
        
        # append attachment section to email body
        if ($attachmentsExist) {
            $body .= "------------------------------------------------------------\n";
            $body .= "Attachments: \n";
            $body .= "\n";            
            $body .= "$attachmentText";
            $body .= "\n";        
        }
        
        # notes ----------------------------------------------------------------          
        
        if ($code == $TEST_FAILED || $code == $TEST_PASSED) {
            $body .= "------------------------------------------------------------\n";
            $body .= "Notes: \n";
            $body .= "\n";        
            
            if ($code == $TEST_FAILED) {
                $body .= "log$hostOS.txt is the output from the test script listed\n";
                $body .= "above (or failed compile attempt). Please note that the -v\n";
                $body .= "option was not selected for this log.\n"; 
            }  
            
            if ($code == $TEST_PASSED) {
                $body .= "log$hostOS.txt is the output from the test script listed\n";
                $body .= "above. Please note that the -v option was not selected for\n";
                $body .= "this log. While no errors occurred during this test, this\n";
                $body .= "log can still be examined to see which tests were run.\n";
            }   
            
            $body .= "\n";        
        }

        # send email / write report ============================================
        
        if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {   
            printEvent("sending email: $mailTo\n");         
    		$email->attach(Type=>'TEXT', Data=>$body);
    		if ($options{'MAIL_METHOD'}[0] eq "sendmail") {
                $email->send();
            } elsif ($options{'MAIL_METHOD'}[0] eq "smtp") {
                $email->send("smtp", $options{'MAIL_METHOD'}[1], Timeout => 30);
            }                
        } else {
            printEvent("generating report: $filename\n"); 
            open (REPORT, ">$options{'TRILINOS_DIR'}[0]/testharness/results/$filename")
                or die "can't open $filename";
            print REPORT $body;
            close REPORT;            
        }
        
        system "rm -f update_log.txt";
        system "rm -f trilinos_configure_log.txt";
        system "rm -f trilinos_build_log.txt";
        system "rm -f test_compile_log.txt";
        system "rm -f $options{'TRILINOS_DIR'}[0]/logErrors.txt";
        system "rm -f $options{'TRILINOS_DIR'}[0]/logMpiErrors.txt";
        system "rm -f $options{'TRILINOS_DIR'}[0]/log$hostOS.txt";
        if ($code == $SUMMARY) {
            system "rm -f event_log.txt";
            system "rm -f invoke-configure-mpi";
            system "rm -f invoke-configure-mpi-original";
            system "rm -f invoke-configure-mpi-final";
            system "rm -f invoke-configure-serial";
            system "rm -f invoke-configure-serial-original";
            system "rm -f invoke-configure-serial-final";
        }
                
    } # report()
    
    ############################################################################
    # appendFile()
    #
    # Prints Test-Harness usage to standart output and exits.
    #   - global variables used: no
    #   - sends mail: no
    #   - args: $log (file to be appended) $logPath (file with absolute path)
    #   - returns: 

    sub appendFile {
        my $log = $_[0];
        my $logPath = $_[1];
    
        my $text = "";
        $text .= "$log:\n";              
        $text .= "\n";
        open (LOG, "<$logPath")
            or die "can't open $log";
        undef $/;           # undefine input record separator
        my $file=<LOG>;     # copy entire file
        $/ = "\n";          # restore it to default newline
        close LOG;
        $text .= "$file\n\n";
        
        return $text;
    
    } # appendFile()
    
    ############################################################################
    # printEvent()
    #
    # Prints an event to standard out and logs it.
    #   - global variables used: no
    #   - sends mail: no
    #   - args: $event to be printed
    #   - returns: 

    sub printEvent {
        my $event = $_[0];
        my $logOnly = $_[1];
    
        my $log = "$options{'TRILINOS_DIR'}[0]/testharness/temp/event_log.txt";
        
        # open log for appending; append; close
        open (LOG, ">>$log")
            or die "can't open $log";            
        print LOG "$event";              
        if (defined $logOnly) { print LOG "$logOnly"; }
        close LOG;
        
        print "$event";
    
    } # printEvent()

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
        print "  -f FILE  : Run test harness normally with given test-harness-config file\n";
        print "\n";
        print "  -nf FILE : Run test harness with given test-harness-config file, but\n";
        print "             don't run tests (for cross-compiling machines).\n";
        print "\n";
        print "  -tf FILE : Run test harness with given test-harness-config file, but\n";
        print "             don't configure and build--only run tests (for cross-\n";
        print "             compiling machines).\n";
        print "\n";
        print "  -p FILE  : Parse given test-harness-config file and exit. This is useful\n";
        print "             for catching errors and inconsistencies without running the\n";
        print "             entire test-harness.\n";
        print "\n";
        print "  -g FILE  : Generate template configuration file (with defaults) named \n";
        print "             FILE and exit\n";
        print "\n";
        print "  -s       : Omit comments from generated configuration file\n";
        print "             (must be of the form: -sg FILE) (do not use -gs)\n";
        print "\n";
        print "  -h       : Print this help page and exit\n";
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
    # parseDependencies()
    #
    # Parses testharness/dependencies and fills global hash of arrays of the 
    # form ({PACKAGE_A, [optionA1, optionA2, ...]}, {PACKAGE_B, [optionB1, ...]})
    #   - global variables used: yes
    #   - sends mail: no
    #   - args: 
    #   - returns: 

    sub parseDependencies {
        my $filename = "dependencies";
        my $line;
        my $name;
        my $value;
        my $continue;
        
        open (IN_FILE, "<$filename")
            or die "can't open $filename";
            
        while ($line = <IN_FILE>) {
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
                
                if (!exists $dependencies{$name}) {  # if there isn't an option with this $name...
                    $dependencies{$name} = ();       # add a new one
                }
                
                while ($value) {
                    $value =~ s/^\s*//;  # trim leading spaces
                    $value =~ s/\s*$//;  # trim trailing spaces
                    # print "\$value: $value\n"; # debugging
                    
                    $value =~ m/^(".*?"|\S+)/;          # grab leftmost value
                    my $v = $1;                         # store temporarily in $v
                    $value =~ s/^$v//;                  # remove $v from remaining list 
                    $v =~ s/"//g;                       # remove any quotes from $v
                    push (@{$dependencies{$name}}, $v); # add $v to %options{$name}     
                }
            } # else (non-comment/blank)
        } # while ($line)
        
        close IN_FILE;
        
        # print the dependencies hash of arrays --------------------------------
        #print "\n\%dependencies:\n\n";                         # debugging
        #for my $name (keys %dependencies) {                    # debugging
        #    my $lastElementIndex = $#{$dependencies{$name}};        # debugging
        #    print "  $name (".($lastElementIndex+1)."): \n";        # debugging
        #    for my $i (0 .. $lastElementIndex) {                    # debugging
        #        print "    $i = $dependencies{$name}[$i]\n";   # debugging
        #    }                                                  # debugging
        #}            
        
    } # parseDependencies()

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
    
        my $configError = 0;
             
        # conversions, fixes, etc. =============================================
             
        # convert <TRILINOS_DIR> psuedo-variable
        for my $name (keys %options) {
            for my $i (0 .. $#{$options{$name}}) {
                $options{$name}[$i] =~ s/<TRILINOS_DIR>/$options{'TRILINOS_DIR'}[0]/;
            }         
        }
        
        # convert <HOST_FILE> psuedo-variable
        for my $name (keys %options) {
            for my $i (0 .. $#{$options{$name}}) {
                if ($options{$name}[$i] =~ m/<HOST_FILE>/) {
                    if (defined $options{'HOST_FILE'} && defined $options{'HOST_FILE'}[0]) {
                        $options{$name}[$i] =~ s/<HOST_FILE>/$options{'HOST_FILE'}[0]/;
                    } else {
                        my $message = "";
                        $message .= "attempting to use <HOST_FILE> value, but HOST_FILE wasn't given\n";
                        if (!$flags{p}) { report($CONFIG_ERROR, $message); }
                        printEvent($message);
                        $configError = 1;
                    }
                } 
            }         
        }      
        
        # if MAKE_FLAGS are specified, but don't begin with a '-', add one
        if (defined $options{'MAKE_FLAGS'} && defined $options{'MAKE_FLAGS'}[0]) {
            if (! $options{'MAKE_FLAGS'}[0] =~ m/^-/) {
                $options{'MAKE_FLAGS'}[0] =~ s/^/-/;
            }
        }
        
        # create results directory if REPORT_METHOD is LOCAL_FILESYSTEM
        if ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM") { 
            system "rm -rf $options{'TRILINOS_DIR'}[0]/testharness/results";
            system "mkdir $options{'TRILINOS_DIR'}[0]/testharness/results";
        }
        
        # delete temp directory and create new one
        chdir "$options{'TRILINOS_DIR'}[0]/testharness";
        system "rm -rf temp";
        system "mkdir temp";
        
        # validations, enforcements, etc. ======================================
        
        # MACHINE_CONFIG_FILE --------------------------------------------------
        
        # enforce MACHINE_CONFIG_FILE requirement
        if (!defined $options{'MACHINE_CONFIG_FILE'} || !defined $options{'MACHINE_CONFIG_FILE'}[0]) {
            my $message = "";
            $message .= "MACHINE_CONFIG_FILE required\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # enforce only one MACHINE_CONFIG_FILE
        elsif (defined $options{'MACHINE_CONFIG_FILE'}[1]) {
            my $message = "";
            $message .= "only one MACHINE_CONFIG_FILE allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # MACHINE_MPI_CONFIG_FILE ----------------------------------------------
        
        # enforce presence of MACHINE_MPI_CONFIG_FILE if MPI_DIR exists
        if (defined $options{'MPI_DIR'} && defined $options{'MPI_DIR'}[0]
            && (!defined $options{'MACHINE_MPI_CONFIG_FILE'} 
            || !defined $options{'MACHINE_MPI_CONFIG_FILE'}[0])) {
            my $message = "";
            $message .= "MACHINE_MPI_CONFIG_FILE must be supplied if MPI_DIR is present\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # enforce only one MACHINE_MPI_CONFIG_FILE
        elsif (defined $options{'MACHINE_MPI_CONFIG_FILE'}[1]) {
            my $message = "";
            $message .= "only one MACHINE_MPI_CONFIG_FILE allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # TRILINOS_CONFIG_FILE -------------------------------------------------
        
        # enforce TRILINOS_CONFIG_FILE requirement
        if (!defined $options{'TRILINOS_CONFIG_FILE'} || !defined $options{'TRILINOS_CONFIG_FILE'}[0]) {
            my $message = "";
            $message .= "TRILINOS_CONFIG_FILE required\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # enforce only one TRILINOS_CONFIG_FILE
        elsif (defined $options{'TRILINOS_CONFIG_FILE'}[1]) {
            my $message = "";
            $message .= "only one TRILINOS_CONFIG_FILE allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # REPORT_METHOD --------------------------------------------------------
        
        # enforce REPORT_METHOD requirement
        if (!defined $options{'REPORT_METHOD'} || !defined $options{'REPORT_METHOD'}[0]) {
            my $message = "";
            $message .= "REPORT_METHOD required\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # enforce only one REPORT_METHOD
        elsif (defined $options{'REPORT_METHOD'}[1]) {
            my $message = "";
            $message .= "only one REPORT_METHOD allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # enforce correct REPORT_METHOD value
        elsif (!($options{'REPORT_METHOD'}[0] eq "EMAIL" 
            || $options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM")) {
            my $message = "";
            $message .= "invalid value for REPORT_METHOD\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # SERIAL_DIR & MPI_DIR -------------------------------------------------
        
        # enforce existence of SERIAL_DIR or MPI_DIR
        if ((!defined $options{'SERIAL_DIR'} || !defined $options{'SERIAL_DIR'}[0])
            && (!defined $options{'MPI_DIR'} || !defined $options{'MPI_DIR'}[0])) {
            my $message = "";
            $message .= "at least one SERIAL_DIR or MPI_DIR required\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
                
        # enforce only one SERIAL_DIR
        elsif (defined $options{'SERIAL_DIR'}[1]) {
            my $message = "";
            $message .= "only one SERIAL_DIR allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
                
        # enforce only one MPI_DIR
        elsif (defined $options{'MPI_DIR'}[1]) {
            my $message = "";
            $message .= "only one MPI_DIR allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # HOST_FILE ------------------------------------------------------------
                
        # enforce only one HOST_FILE
        if (defined $options{'HOST_FILE'}[1]) {
            my $message = "";
            $message .= "only one HOST_FILE allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # MPI_STARTUP_CMD ------------------------------------------------------
                
        # enforce only one MPI_STARTUP_CMD
        if (defined $options{'MPI_STARTUP_CMD'}[1]) {
            my $message = "";
            $message .= "only one MPI_STARTUP_CMD allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # MPI_SHUTDOWN_CMD -----------------------------------------------------
                
        # enforce only one MPI_SHUTDOWN_CMD
        if (defined $options{'MPI_SHUTDOWN_CMD'}[1]) {
            my $message = "";
            $message .= "only one MPI_SHUTDOWN_CMD allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # MAKE_FLAGS -----------------------------------------------------------
                
        # enforce only one MAKE_FLAGS
        if (defined $options{'MAKE_FLAGS'}[1]) {
            my $message = "";
            $message .= "only one MAKE_FLAGS value allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # CVS_UPDATE -----------------------------------------------------------
        
        # enforce CVS_UPDATE requirement
        if (!defined $options{'CVS_UPDATE'} || !defined $options{'CVS_UPDATE'}[0]) {
            my $message = "";
            $message .= "CVS_UPDATE required\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # enforce only one CVS_UPDATE
        elsif (defined $options{'CVS_UPDATE'}[1]) {
            my $message = "";
            $message .= "only one CVS_UPDATE allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # enforce correct CVS_UPDATE value
        elsif (!($options{'CVS_UPDATE'}[0] eq "YES" 
            || $options{'CVS_UPDATE'}[0] eq "NO")) {
            my $message = "";
            $message .= "invalid value for CVS_UPDATE\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # CVS_CMD --------------------------------------------------------------
        
        # enforce presence of CVS_CMD if CVS_UPDATE is set to YES
        if (defined $options{'CVS_UPDATE'} && defined $options{'CVS_UPDATE'}[0]) {
            if ($options{'CVS_UPDATE'}[0] eq "YES"
                && (!defined $options{'CVS_CMD'} 
                || !defined $options{'CVS_CMD'}[0])) {
                my $message = "";
                $message .= "CVS_CMD must be supplied if CVS_UPDATE is set to YES\n";
                if (!$flags{p}) { report($CONFIG_ERROR, $message); }
                printEvent($message);
                $configError = 1;
            }
        }
                
        # enforce only one CVS_CMD
        elsif (defined $options{'CVS_CMD'}[1]) {
            my $message = "";
            $message .= "only one CVS_CMD allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # REPORT_NAMES ---------------------------------------------------------
        
        # enforce presence of REPORT_NAMES if REPORT_METHOD is set to LOCAL_FILESYSTEM
        if (defined $options{'REPORT_METHOD'} && defined $options{'REPORT_METHOD'}[0]) {
            if ($options{'REPORT_METHOD'}[0] eq "LOCAL_FILESYSTEM"
                && (!defined $options{'REPORT_NAMES'} 
                || !defined $options{'REPORT_NAMES'}[0])) {
                my $message = "";
                $message .= "REPORT_NAMES must be supplied if REPORT_METHOD is set to LOCAL_FILESYSTEM\n";
                if (!$flags{p}) { report($CONFIG_ERROR, $message); }
                printEvent($message);
                $configError = 1;
            }
        }
        
        # enforce correct REPORT_NAMES value
        elsif (!($options{'REPORT_NAMES'}[0] eq "ORDER" 
            || $options{'REPORT_NAMES'}[0] eq "EVENT")) {
            my $message = "";
            $message .= "invalid value for REPORT_NAMES\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
                
        # enforce only one REPORT_NAMES
        elsif (defined $options{'REPORT_NAMES'}[1]) {
            my $message = "";
            $message .= "only one REPORT_NAMES value allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # MAIL_METHOD ----------------------------------------------------------
        
        # enforce presence of MAIL_METHOD if REPORT_METHOD is set to EMAIL
        if (defined $options{'REPORT_METHOD'} && defined $options{'REPORT_METHOD'}[0]) {
            if ($options{'REPORT_METHOD'}[0] eq "EMAIL"
                && (!defined $options{'MAIL_METHOD'} 
                || !defined $options{'MAIL_METHOD'}[0])) {
                my $message = "";
                $message .= "MAIL_METHOD must be supplied if REPORT_METHOD is set to EMAIL\n";
                if (!$flags{p}) { report($CONFIG_ERROR, $message); }
                printEvent($message);
                $configError = 1;
            }
        }
                
        # enforce only one MAIL_METHOD if first isn't smtp
        elsif (defined $options{'MAIL_METHOD'}[1] && $options{'MAIL_METHOD'}[0] ne "smtp") {
            my $message = "";
            $message .= "only one MAIL_METHOD value allowed if first value isn't \"smtp\"\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
                
        # enforce only two MAIL_METHOD if first is smtp
        elsif ($options{'MAIL_METHOD'}[0] eq "smtp" && defined $options{'MAIL_METHOD'}[2]) {
            my $message = "";
            $message .= "only two MAIL_METHOD values allowed (if first value is \"smtp\")\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # SUMMARY_EMAIL --------------------------------------------------------
        
            # No enforcable requirements
        
        # ALL_EMAILS -----------------------------------------------------------
        
            # No enforcable requirements
        
        # SEND_TO_DEFAULTS -----------------------------------------------------
        
        # enforce SEND_TO_DEFAULTS requirement
        if (!defined $options{'SEND_TO_DEFAULTS'} || !defined $options{'SEND_TO_DEFAULTS'}[0]) {
            my $message = "";
            $message .= "SEND_TO_DEFAULTS required\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # enforce only one SEND_TO_DEFAULTS
        elsif (defined $options{'SEND_TO_DEFAULTS'}[1]) {
            my $message = "";
            $message .= "only one SEND_TO_DEFAULTS allowed\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # enforce correct SEND_TO_DEFAULTS value
        elsif (!($options{'SEND_TO_DEFAULTS'}[0] eq "YES" 
            || $options{'SEND_TO_DEFAULTS'}[0] eq "NO")) {
            my $message = "";
            $message .= "invalid value for SEND_TO_DEFAULTS\n";
            if (!$flags{p}) { report($CONFIG_ERROR, $message); }
            printEvent($message);
            $configError = 1;
        }
        
        # print the options hash of arrays -------------------------------------
        #print "\n\%options:\n\n";                       # debugging
        #for my $name (keys %options) {                  # debugging
        #    my $lastElementIndex = $#{$options{$name}};      # debugging
        #    print "  $name (".($lastElementIndex+1)."): \n"; # debugging
        #    for my $i (0 .. $lastElementIndex) {             # debugging
        #        print "    $i = $options{$name}[$i]\n"; # debugging
        #    }                                           # debugging
        #}            
        
        # if config error and not just parsing, die
        if ($configError && !$flags{p}) { 
            die " *** test-harness-config error - aborting test-harness ***\n"; 
        }
        
        # report successful parse and validation if parse flag
        elsif (!$configError && $flags{p}) {
            print "test-harness config file \"$flags{p}\" successfully parsed and validated\n"
        }
            
    } # validateOptions()
    
    ############################################################################
    # generateConfig()
    #
    # Generates Test-Harness template config file named in current directory
    # named "test-harness-config" and exits.
    #   - global variables used: no
    #   - sends mail: no
    #   - args: string/boolean filename, boolean short
    #   - returns: 

    sub generateConfig {
        my $filename = $_[0];
        my $short = $_[1];
        my $silent = $_[2];
        my $outFile;
                
        if (!$silent) {
            open (outFile, "> $filename")
                or die "can't open $filename";
        
            print outFile "# test-harness-config\n";  
            print outFile "\n";
        }

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
        
        if (!$silent) {
            print outFile "#===============================================================================\n";
            print outFile "# Frequently changed configuration options\n";
            print outFile "#===============================================================================\n";
            print outFile "\n";
        }
        
        if (!$short) {        
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# The name of the file in \"<TRILINOS_DIR>/testharness/elements-machine\"\n";
            print outFile "# containing the machine-dependent configure options for this machine.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES\n";
            print outFile "\n";
        }
        
        push (@optionsOrder, "MACHINE_CONFIG_FILE");
        if (!$silent) { print outFile "MACHINE_CONFIG_FILE             = \n"; }
        
        if (!$short) {        
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# The name of the file in \"<TRILINOS_DIR>/testharness/elements-machine\"\n";
            print outFile "# containing the machine-independent mpi configure options for this\n";
            print outFile "# configuration of the test-harness.\n";
            print outFile "\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES if MPI_DIR is supplied\n";
            print outFile "\n";
        }
        
        push (@optionsOrder, "MACHINE_MPI_CONFIG_FILE");
        if (!$silent) { print outFile "MACHINE_MPI_CONFIG_FILE         = \n"; }
        
        if (!$short) {        
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# The name of the file in \"<TRILINOS_DIR>/testharness/elements-trilinos\"\n";
            print outFile "# containing the machine-independent configure options for this configuration\n";
            print outFile "# of the test-harness.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES\n";
            print outFile "\n";
        }
        
        push (@optionsOrder, "TRILINOS_CONFIG_FILE");
        if (!$silent) { print outFile "TRILINOS_CONFIG_FILE            = \n"; }
        
        if (!$short) {      
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";         
            print outFile "# Indicate how the test-harness should report results. If the value is EMAIL,\n";              
            print outFile "# then emails will be sent to either the script owner, the package regression\n";
            print outFile "# list, or the Trilinos regression list (see other email options below). If the\n";
            print outFile "# value is LOCAL_FILESYSTEM, then no emails will be sent at all and the results\n";
            print outFile "# will be written to testharness/results/. This directory will be blown away at\n";
            print outFile "# the beginning of each run of the testharness that has REPORT_METHOD set to\n";
            print outFile "# LOCAL_FILESYSTEM.\n";      
            print outFile "#\n";    
            print outFile "# - recognized values: EMAIL LOCAL_FILESYSTEM\n";  
            print outFile "# - multiple values recognized: NO\n"; 
            print outFile "# - value required: YES\n";
            print outFile "\n";
        }
        
        push (@optionsOrder, "REPORT_METHOD");
        if (!$silent) { print outFile "REPORT_METHOD                   = EMAIL\n"; }
        
        if (!$silent) { 
            print outFile "\n";
            print outFile "#===============================================================================\n";
            print outFile "# Less frequently changed configuration options\n";
            print outFile "#===============================================================================\n";
            print outFile "\n";
        }
		
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
        
        push (@optionsOrder, "SERIAL_DIR");
        if (!$silent) { print outFile "SERIAL_DIR                      = SERIAL\n"; }
        
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
        
        push (@optionsOrder, "MPI_DIR");
        if (!$silent) { print outFile "MPI_DIR                         = MPI\n"; }
        
        if (!$short) {    
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";     
            print outFile "# (LAM only) path and filename of the file containing the names of the machines\n"; 
            print outFile "# to be used for parallel jobs. If this file doesn't exist, parallel jobs will\n";
            print outFile "# be run on the local machine only. (If this is the case, make sure\n";  
            print outFile "# <HOSTFILE> is removed from the value of MPI_STARTUP_CMD.)\n";           
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: NO\n";
            print outFile "# - the absolute path of the Trilinos directly above 'testharness' can be\n";
            print outFile "#   referred to with the value <TRILINOS_DIR>\n";          
            print outFile "\n";
        }
        
        push (@optionsOrder, "HOST_FILE");
        if (!$silent) { print outFile "HOST_FILE                       = <TRILINOS_DIR>/hostfile\n"; }
        
        if (!$short) {      
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# Specify the command to start up the MPI implementation on this machine.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: NO\n";
            print outFile "# - the value of the HOST_FILE option can be referred to with the value\n";
            print outFile "#   <HOST_FILE>\n";   
            print outFile "\n";
        }
        
        push (@optionsOrder, "MPI_STARTUP_CMD");
        if (!$silent) { print outFile "MPI_STARTUP_CMD                 = \"lamboot <HOST_FILE> -v\"\n"; }
        
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
        
        push (@optionsOrder, "MPI_SHUTDOWN_CMD");
        if (!$silent) { print outFile "MPI_SHUTDOWN_CMD                = lamhalt\n"; }
        
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
        
        push (@optionsOrder, "MAKE_FLAGS");
        if (!$silent) { print outFile "MAKE_FLAGS                      = \n"; }
        
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
        
        push (@optionsOrder, "CVS_UPDATE");
        if (!$silent) { print outFile "CVS_UPDATE                      = YES \n"; }
        
        if (!$short) {      
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# CVS command on this system. Note that CVS access must not require a password.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES if CVS_UPDATE is set to YES\n";
            print outFile "\n";
        }
        
        push (@optionsOrder, "CVS_CMD");
        if (!$silent) { print outFile "CVS_CMD                         = cvs\n"; }
        
        if (!$short) {    
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";     
            print outFile "# Indicate how report text files should be named. ORDER will result in the files\n";
            print outFile "# being named such that they will sort in the order that they occurred. EVENT\n";
            print outFile "# will result in the files being named such that they will sort according to\n";
            print outFile "# the type of event being reported.\n"; 
            print outFile "#\n";
            print outFile "# - recognized values: ORDER EVENT\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES if REPORT_METHOD is set to LOCAL_FILESYSTEM\n";      
            print outFile "\n";
        }
        
        push (@optionsOrder, "REPORT_NAMES");
        if (!$silent) { print outFile "REPORT_NAMES                    = ORDER\n"; }
        
        if (!$short) {    
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";     
            print outFile "# Indicate how mail should be sent. The unix program sendmail is the default.\n"; 
            print outFile "#\n";
            print outFile "# - recognized values: sendmail\n";
            print outFile "#                      smtp <mail_server>\n";
            print outFile "# - multiple values recognized: NO (except \"smtp\" and \"<mail_server>\")\n";
            print outFile "# - value required: YES if REPORT_METHOD is set to EMAIL\n";      
            print outFile "\n";
        }
        
        push (@optionsOrder, "MAIL_METHOD");
        if (!$silent) { print outFile "MAIL_METHOD                     = sendmail\n"; }
        
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
        
        push (@optionsOrder, "SUMMARY_EMAIL");
        if (!$silent) { print outFile "SUMMARY_EMAIL                   = \n"; }
        
        if (!$short) {        
            print outFile "\n";
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# List the email addresses to which all summaries and error messages will be\n";
            print outFile "# sent regardless of where else they also get sent (script owners, regression\n";
            print outFile "# lists).\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: YES\n";
            print outFile "# - value required: NO\n";
            print outFile "\n";
        }
        
        push (@optionsOrder, "ALL_EMAILS");
        if (!$silent) { print outFile "ALL_EMAILS                      = \n"; }
        
        if (!$short) {        
            print outFile "\n";
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# Use this option to enable/disable sending email to the default recipients.\n";
            print outFile "# YES indicates that emails should be sent to all of the usual places (script\n";
            print outFile "# owners, package lists, and Trilinos lists, in addition to ALL_EMAILS and\n";
            print outFile "# SUMMARY_EMAIL). NO indicates that they should only be sent to ALL_EMAILS\n";
            print outFile "# and SUMMARY_EMAIL\n";
            print outFile "#\n";
            print outFile "# - recognized values: YES NO \n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES\n";
            print outFile "\n";
        }
        
        push (@optionsOrder, "SEND_TO_DEFAULTS");
        if (!$silent) { print outFile "SEND_TO_DEFAULTS                = NO\n"; }
        
        if (!$silent) { 
            print outFile "\n";
            print outFile "# end test-harness-config";
                
            close outFile;
        }
    } # generateConfig()
