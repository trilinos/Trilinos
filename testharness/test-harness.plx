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

# Report variables
my %summary;                # overall test-harness summary
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
my $TEST_COMPILE_ERROR = 6;           # test compile failed
my $TEST_FAILED = 7;                  # test failed
my $TEST_PASSED = 8;                  # test passed
my $SUMMARY = 9;                      # test-harness summary
        
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

# Send summary email
sendMail($SUMMARY);

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
        getopts("f:p:g:snth", \%flags);
        
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
    # prepareVariables()
    #
    # Prepares global varibles.
    #   - global variables used: yes
    #   - sends mail: on error
    #   - args: 
    #   - returns: 

    sub prepareVariables {    
        
        $codes[$FILE_SYSTEM_ERROR] = "file system error";           
        $codes[$SYSTEM_COMMAND_ERROR] = "system command error";        
        $codes[$CONFIG_ERROR] = "test-harness config error";
        $codes[$UPDATE_ERROR] = "cvs udate error";
        $codes[$TRILINOS_CONFIGURE_ERROR] = "Trilinos configure error";
        $codes[$TRILINOS_BUILD_ERROR] = "Trilinos build error";
        $codes[$TEST_COMPILE_ERROR] = "test compile error";
        $codes[$TEST_FAILED] = "test failed";
        $codes[$TEST_PASSED] = "test passed";
        $codes[$SUMMARY] = "test-harness summary";
            
        $summary{$FILE_SYSTEM_ERROR} = ();
        $summary{$SYSTEM_COMMAND_ERROR} = ();
        $summary{$CONFIG_ERROR} = ();
        $summary{$UPDATE_ERROR} = ();
        $summary{$TRILINOS_CONFIGURE_ERROR} = ();
        $summary{$TRILINOS_BUILD_ERROR} = ();
        $summary{$TEST_COMPILE_ERROR} = ();
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
        $emails{"Trilinos"} = "trilinos";
        $emails{"unknown"} = "trilinos";
        
    } # prepareVariables()

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
                sendMail($UPDATE_ERROR);
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
        if (defined $options{'MPI_DIR'} && defined $options{'MPI_DIR'}[0]) { 
            my $command = "$options{'MPI_STARTUP_CMD'}[0]";
            my $commandFailed = system $command; 
            if ($commandFailed) {
                sendMail($SYSTEM_COMMAND_ERROR, $command);
                print $command;          
                die " *** error running system command - aborting test-harness ***\n";
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
        
        if (!$flags{t}) {
            # delete and recreate mpi dir
            system "rm -rf $options{'MPI_DIR'}[0]";
            system "mkdir $options{'MPI_DIR'}[0]";
            
            # delete and recreate serial dir
            system "rm -rf $options{'SERIAL_DIR'}[0]";
            system "mkdir $options{'SERIAL_DIR'}[0]";
        }
        
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
            if ($buildDir[$j] eq $options{'MPI_DIR'}[0]) {
                $comm = "mpi";
            } elsif ($buildDir[$j] eq $options{'SERIAL_DIR'}[0]) {
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
                    sendMail($FILE_SYSTEM_ERROR, $message, $comm);
                    print $message;
                    die " *** file missing - aborting test-harness ***\n";
                } 
                
                # test for read permissions on invoke-configure
                if (!-r "invoke-configure") {
                    my $message = "";
                    $message .= "$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/";
                    $message .= "invoke-configure must be readable\n";
                    sendMail($FILE_SYSTEM_ERROR, $message, $comm);
                    print $message;
                    die " *** file permission wrong - aborting test-harness ***\n";
                } 
                
                # test for executable permission on invoke-configure
                if (!-x "invoke-configure") {
                    my $message = "";
                    $message .= "$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/";
                    $message .= "invoke-configure must be executable\n";
                    sendMail($FILE_SYSTEM_ERROR, $message, $comm);
                    print $message;
                    die " *** file permission wrong - aborting test-harness ***\n";
                }
                
                # should we quit trying? (no invoke-configure left)
                my $quitTrying = 0;
                
                # configure --------------------------------------------------------
                
                # configure passed?
                my $configurePassed = 0;            
                
                my $COUNT = 0;
                while (!$configurePassed && !$quitTrying) {
                    if (++$COUNT >= 3) { $quitTrying = 1; }
                
                    # attempt to configure
                    my $configureOutput = configure($buildDir[$j]);
                    
                    # configure failed
                    if ($configureOutput) {    
                        
                        # fix invoke configure
                        my $log = "$options{'TRILINOS_DIR'}[0]/testharness/temp/trilinos_configure_log.txt";
                        my $invokeConfigure ="$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/invoke-configure";                    
                        my $brokenPackage = fixInvokeConfigure($log, $invokeConfigure);               
                        
                        # remove broken package
                        system "rm -rf $options{'TRILINOS_DIR'}[0]/$buildDir[$j]/packages/$brokenPackage";
                        
                        print "$comm - Trilinos configuration failed--$brokenPackage broke\n";  
                        
                        # there is no invoke-configure left or it's empty
                        if (!-f $invokeConfigure || -z $invokeConfigure) {   
                            $quitTrying = 1;
                        }
                        
                        # send email
                        sendMail($TRILINOS_CONFIGURE_ERROR, $brokenPackage, $comm);
                    } 
                    
                    # configure succeeded
                    else {     
                        $configurePassed = 1;
                        print "$comm - Trilinos configured successfully\n";                      
                    }
                
                    # build --------------------------------------------------------
                    
                    # build passed?            
                    my $buildPassed = 0;            
                    
                    while ($configurePassed && !$buildPassed && !$quitTrying) {
                        if (++$COUNT >= 3) { $quitTrying = 1; }
                    
                        # attempt to build
                        my $buildOutput = build($buildDir[$j]);
                        
                        # build failed
                        if ($buildOutput) {
                                       
                            # fix invoke configure         
                            my $log = "$options{'TRILINOS_DIR'}[0]/testharness/temp/trilinos_build_log.txt";
                            my $invokeConfigure ="$options{'TRILINOS_DIR'}[0]/$buildDir[$j]/invoke-configure";                    
                            my $brokenPackage = fixInvokeConfigure($log, $invokeConfigure);
                            
                            # remove broken package
                            system "rm -rf $options{'TRILINOS_DIR'}[0]/$buildDir[$j]/packages/$brokenPackage";
                            
                            print "$comm - Trilinos build failed--$brokenPackage broke\n";  
                        
                            # there is no invoke-configure left or it's empty
                            if (!-f $invokeConfigure || -z $invokeConfigure) {   
                                $quitTrying = 1;
                            }
                            
                            # force reconfigure
                            $configurePassed = 0;
                            
                            # send email
                            sendMail($TRILINOS_BUILD_ERROR, $brokenPackage, $comm); 
                        } 
                        
                        # build succeeded
                        else {    
                            $buildPassed = 1;
                            print "$comm - Trilinos built successfully\n";                  
                        }
                        
                    } # while (configurePassed && !buildPassed)
                    
                } # while (!configurePassed)                       
            } # if (-t)
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
                                        print "$testNameOnly - Test failed\n";  
                                        sendMail($TEST_FAILED, $testFailed, $comm, $testDir, $potentialScript); 
                                    } else {                                    
                                        print "$testNameOnly - Test passed\n";  
                                        sendMail($TEST_PASSED, $testFailed, $comm, $testDir, $potentialScript);
                                    }
                                    
                                } # if (executable)
                            } # foreach ($potentialScript)  
                        } # if (-f $potentialTestDir)
                                          
                    } # unless (unsupported)
                } # foreach $testDir
            } # if (!-n)
            
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
                sendMail($FILE_SYSTEM_ERROR, $message);
                print $message;
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
            sendMail($FILE_SYSTEM_ERROR, $message);
            print $message;
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
            sendMail($FILE_SYSTEM_ERROR, $message);
            print $message;
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
            print INVOKE_CONFIGURE_SERIAL $rawInvokeConfigureSERIAL;
            close INVOKE_CONFIGURE_SERIAL;
            
            # move invoke-configure file into place
            my $command;                     
            $command = "mv $options{'TRILINOS_DIR'}[0]/testharness/temp/invoke-configure-serial ";
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
    	    $file =~ m/.*\/packages\/(.*?)\//ms;
	        if (defined $1) { $brokenPackage = $1; }
        }
	    
	    if (defined $brokenPackage) {
            $brokenPackage =~ s/^\s*//;             # trim leading spaces
            $brokenPackage =~ s/\s*$//;             # trim trailing spaces
    	    $brokenPackage = lc($brokenPackage);    # convert to lowercase
            $brokenPackage =~ s/ /_/g;              # convert spaces to underscores
        } else {
            return "unknown";
        }
        
        system "cp $invokeConfigure $invokeConfigure-broken";
        
        # open invoke-configure for reading
        open (INVOKE_CONFIGURE, "<$invokeConfigure")
            or die "can't open $invokeConfigure";
            
        # parse it, extracting lines to keep
        my $newInvokeConfigure = "";
        while (my $line = <INVOKE_CONFIGURE>) {
            $line =~ s/^\s*//;  # trim leading spaces
            $line =~ s/\s*$//;  # trim trailing spaces
            
            # compare lines to package dependencies
            my $dropLine = 0;
            my $lastElementIndex = $#{$dependencies{$brokenPackage}};
            for (my $i=0; $i<=$lastElementIndex; $i++) { 
                if ($line =~ m/$dependencies{$brokenPackage}[$i]/) {
                    $dropLine = 1;   
                }   
            }    
            
            # write line if it isn't a dependency
            if (!$dropLine) {                        
                $newInvokeConfigure .= "$line\n";
            }
        } # while ($line)
        close INVOKE_CONFIGURE;
        
        # open invoke-configure for writing
        open (INVOKE_CONFIGURE, ">$invokeConfigure")
            or die "can't open $invokeConfigure";
        
        # write new invoke configure
        chomp ($newInvokeConfigure);        # remove last newline 
        $newInvokeConfigure =~ s/.*\\$//s;  # remove last line-continuation
        print INVOKE_CONFIGURE $newInvokeConfigure;
        close INVOKE_CONFIGURE;    
        
        # return broken package name
        return $brokenPackage;
        
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
        if ($code == $FILE_SYSTEM_ERROR || $code == $SYSTEM_COMMAND_ERROR ||
               $code == $CONFIG_ERROR || $code == $UPDATE_ERROR) {
            push (@{$summary{$code}}, $message);
        }         
        if ($code == $TRILINOS_CONFIGURE_ERROR || $code == $TRILINOS_BUILD_ERROR) {
            push (@{$summary{$code}}, $comm);
        }        
        if ($code == $TEST_COMPILE_ERROR || $code == $TEST_FAILED || $code == $TEST_PASSED) {
            push (@{$summary{$code}}, $testName);
        }

        # compile list of mail recipients --------------------------------------
        
        my $mailTo = "";
    
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
        	        print "error - unable to detect package regression list\n";
            	    print "  sending to trilinos-regression instead\n";
                    $mailTo .= "trilinos-regression\@software.sandia.gov, ";
        	    }
            } # if (configure/build-related)
        
            # test-related
            if ($code == $TEST_COMPILE_ERROR || $code == $TEST_FAILED || $code == $TEST_PASSED) {
        
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
            	        print "error - unable to detect package regression list\n";
            	        print "  sending to trilinos-regression instead\n";
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
        my $lastElementIndex = $#{$options{'ALL_EMAILS'}};
        for (my $i=0; $i<=$lastElementIndex; $i++) {
            if ($mailTo !~ m/$options{'ALL_EMAILS'}[$i]/i) {
                $mailTo .= "$options{'ALL_EMAILS'}[$i], ";
            } 
        }
        
        # append SUMMARY_EMAIL
        if ($code == $SUMMARY) { 
            $lastElementIndex = $#{$options{'SUMMARY_EMAIL'}};
            for (my $i=0; $i<=$lastElementIndex; $i++) {
                if ($mailTo !~ m/$options{'SUMMARY_EMAIL'}[$i]/i) {
                    $mailTo .= "$options{'SUMMARY_EMAIL'}[$i], ";
                }
            }
        } 
        
        # clean up list of recipients    
        $mailTo =~ s/^\s*//;        # trim leading spaces
        $mailTo =~ s/\s*$//;        # trim trailing spaces
        $mailTo =~ s/,*$//;         # trim trailing comma    
        
        print "sending email: $mailTo\n"; 

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
        if ($code == $TRILINOS_CONFIGURE_ERROR || $code == $TRILINOS_BUILD_ERROR) {
            $body .= "Result:           $codes[$code] - $message broke\n";
        } else {
            $body .= "Result:           $codes[$code]\n";
        }
        
        $body .= "\n";        
        
        # fatal error ----------------------------------------------------------
        if ($code == $FILE_SYSTEM_ERROR || $code == $SYSTEM_COMMAND_ERROR ||
            $code == $CONFIG_ERROR || $code == $UPDATE_ERROR) {
        
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
        if ($code == $UPDATE_ERROR && -f "update_log.txt") {
            $attachment = 1;
            $attachmentText .= "    update_log.txt\n";
            $email->attach(Type=>'TEXT', Path=>'update_log.txt', Disposition=>'attachment');
        }
        
        # trilinos configure failed
        if ($code == $TRILINOS_CONFIGURE_ERROR && -f "trilinos_configure_log.txt") {
            $attachment = 1;
            $attachmentText .= "    trilinos_configure_log.txt\n";
            $email->attach(Type=>'TEXT', Path=>'trilinos_configure_log.txt', Disposition=>'attachment');
        }       
        
        # trilinos build failed
        if ($code == $TRILINOS_BUILD_ERROR && -f "trilinos_build_log.txt") {
            $attachment = 1;
            $attachmentText .= "    trilinos_build_log.txt\n";
            $email->attach(Type=>'TEXT', Path=>'trilinos_build_log.txt', Disposition=>'attachment');
        }
        
        # test compile failed
        if ($code == $TEST_COMPILE_ERROR && -f "test_compile_log.txt") {
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
        
        # invoke-configure
        my $buildDir = "";
        if (defined $comm && $comm eq "mpi") {$buildDir = $options{'MPI_DIR'}[0]; }
        if (defined $comm && $comm eq "serial") {$buildDir = $options{'SERIAL_DIR'}[0]; }
        my $invokeConfigure = "$options{'TRILINOS_DIR'}[0]/$buildDir/invoke-configure";
        
        if (($code == $TRILINOS_CONFIGURE_ERROR || $code == $TRILINOS_BUILD_ERROR)
            && -f $invokeConfigure."-broken" && ! -z $invokeConfigure."-broken") {
            $attachment = 1;
            $attachmentText .= "    invoke-configure-broken\n";
            $email->attach(Type =>'TEXT', Path=>"$invokeConfigure-broken", Disposition=>'attachment');
        } 
        
        elsif (!$code == $TRILINOS_CONFIGURE_ERROR && !$code == $TRILINOS_BUILD_ERROR
            && -f $invokeConfigure && ! -z $invokeConfigure) {
            $attachment = 1;
            $attachmentText .= "    invoke-configure\n";
            $email->attach(Type =>'TEXT', Path=>"$invokeConfigure", Disposition=>'attachment');
        }
        
        if ($attachment) {
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
        
        # print the summary hash of arrays -------------------------------------
        
        if ($code == $SUMMARY) {
            $body .= "------------------------------------------------------------\n";
            $body .= "Summary: \n";
            
            for my $id (keys %summary) { 
                my $lastElementIndex = $#{$summary{$id}}; 
                $body .= "\n- $codes[$id] (".($lastElementIndex+1)."): \n\n";
                for my $i (0 .. $lastElementIndex) {
                    $body .= "    $summary{$id}[$i]\n";
                }
            } 
        }           

        # send email ===========================================================
        
        if ($options{'REPORT_METHOD'}[0] eq "EMAIL") {        
    		$email->attach(Type=>'TEXT', Data=>$body);
            $email->send();
        } else {
            print "-- REPRT_METHOD not set to EMAIL - not sending email\n";
            print "-- other report methods not yet implemented\n";
        }
        
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
        
        # if ALL_EMAILS are not specified, assign the empty string
        if (!defined $options{'ALL_EMAILS'} || !defined $options{'ALL_EMAILS'}[0]) {
            $options{'ALL_EMAILS'} = ();
            push (@{$options{'ALL_EMAILS'}}, "");
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
            print outFile "# - value required: YES if MPI_DIR is supplied\n";
            print outFile "\n";
        }
        
        print outFile "MACHINE_MPI_CONFIG_FILE         = \n";
        
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
        
        print outFile "TRILINOS_CONFIG_FILE            = \n";
        
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
        
        print outFile "HOST_FILE                       = <TRILINOS_DIR>/hostfile\n";
        
        if (!$short) {    
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";     
            print outFile "# Indicate how mail should be sent. The unix program sendmail is the default.\n"; 
            print outFile "#\n";
            print outFile "# - recognized values: sendmail smtp::<mail_server> \n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES if REPORT_METHOD is set to EMAIL\n";      
            print outFile "\n";
        }
        
        print outFile "MAIL_METHOD                     = sendmail\n";
        
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
            print outFile "# - value required: YES if CVS_UPDATE is set to YES\n";
            print outFile "\n";
        }
        
        print outFile "CVS_CMD                         = cvs\n";
        
        if (!$short) {      
            print outFile "\n";  
            print outFile "#-------------------------------------------------------------------------------\n";
            print outFile "# Specify the command to start up the MPI implementation on this machine.\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: NO\n";
            print outFile "# - value required: YES MPI_DIR is supplied\n";
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
        
        print outFile "MPI_SHUTDOWN_CMD                = lamhalt\n";
        
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
            print outFile "# List the email addresses to which all summaries and error messages will be\n";
            print outFile "# sent regardless of where else they also get sent (script owners, regression\n";
            print outFile "# lists).\n";
            print outFile "#\n";
            print outFile "# - multiple values recognized: YES\n";
            print outFile "# - value required: NO\n";
            print outFile "\n";
        }
        
        print outFile "ALL_EMAILS                      = \n";
        
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
        
        print outFile "SEND_TO_DEFAULTS                = NO\n";
        
        print outFile "\n";
        print outFile "# end test-harness-config";
            
        close outFile;
    } # genConfigTemp()
