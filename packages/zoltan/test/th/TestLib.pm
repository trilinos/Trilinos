
#use Statistics::TTest;
use POSIX;
package   TestLib;
require   Exporter;
@ISA    = qw(Exporter);
@EXPORT = qw(run_regression_test);

use File::Find;
use Text::ParseWords;

##
## Setup lists of test types, options, frequencies
##
##@tests = qw(regression performance unit memory application);
##@testoptions = qw(param benchmark analysiscode software);
##@test_frequencies = qw(daily weekly monthly smoke all);
#
# Misc
#
$verbose = 0;
#
# Default key info: software/host/configinfo
#
$pwd = `pwd`;
chomp($pwd);
$software = ( split(/\//,$pwd) ) [-1];
$host = `hostname`;
chomp($host);
$configinfo = "unknown";
#
# Configure timing to be more accurate if possible
#
eval("require Time::HiRes");
$hires = 1;
if ($@) {
   require Time::Local;
   $hires=0;
   }

#
# Set the test name
#
sub set_name {
  $software = shift @_;
}


#
# Set the site configuration info name
#
sub init_siteconfig {
  require "siteconfig.pl";
  #$siteconfig = shift @_;
}


sub extract_name {
  my $filename = shift @_;
  open(INPUT,$filename);
  my @lines=<INPUT>;
  close(INPUT);
  my @tmp = grep(/<experimental-study/,@lines);
  my $tmpf;
  foreach $tmpf (@tmp) {
    $tmpf =~ s/(.*)name="(.*)".*/\2/;
    chomp($tmpf);
    #print "HERE :$tmpf:\n";
    return $tmpf;
    }
  return $tmpf;
  }

#
# Setup signal handling mechanism
#
sub my_signal_catcher {
    die "Terminating testdriver due to external signal!\n";
}
$SIG{'INT'} = 'my_signal_catcher';
$SIG{'HUP'} = 'my_signal_catcher';
$SIG{'QUIT'} = 'my_signal_catcher';
$SIG{'TERM'} = 'my_signal_catcher';


##
## Strip out bad xml characters in error messages
## call with a string reference
##
sub scrub {
  my $unsafe='<>&"';
  my $uri="use URI::Escape;";	# someone already wrote this package!
  my $lineref = shift;

  eval $uri;			# but it might not be installed
  if ($@) {
    eval "\$\$lineref =~ tr /$unsafe/+/;";	# not installed
  } else {
    $$lineref = uri_escape($$lineref, $unsafe);
  }
}

##
## Compute the time using 'gettimeofday' if that's available
##
sub mytime {
  if ($hires == 1) {
     return [Time::HiRes::gettimeofday()];
  } else {
     return time;
  }
}

##
## Compute time differences using 'gettimeofday' if that's available
##
sub time_diff {
  my $t0 = shift @_;
  my $t1 = shift @_;
  if ($hires == 1) {
     return Time::HiRes::tv_interval($t0,$t1);
  } else {
     return $t1-$t0;
  }
}

##
## Find extract all '_<name>' options from a list of 
##    <option>=<value> pairs
##
sub process_testoptions {
  my $testoptions = shift @_;
  #
  # Extract test options and put other options in the $expoptions
  # string.
  #
  @auxlist = ();
  @words = &quotewords('\s+', 0, $testoptions);
  while (@words) {
    local $tmp = shift @words;
    push @auxlist, $tmp;
    local @ttmp = split(/=/,$tmp);
    if (substr($ttmp[0],0,1) eq "_") {
       ${@_[0]}{$ttmp[0]} = $ttmp[1];
       }
    }
  #
  # Analyze test options, and add default values where needed
  #
  if ( !defined(${@_[0]}{"_termination"}) ) {
     ${@_[0]}{"_termination"} = "Accuracy";
  }
  if ( defined(${@_[0]}{"_optimum"}) && !defined(${@_[0]}{"_value"}) ) {
     ${@_[0]}{"_value"} = ${@_[0]}{"_optimum"};
  }
  if (!defined(${@_[0]}{"_ctolerance"})) {
     ${@_[0]}{"_ctolerance"} = "0.0";
  }
  #
  # Process the _value parameter
  #
  if (defined(${@_[0]}{"_value"}) &&
      (!( (${@_[0]}{"_value"} eq "Infinity") ||
          (${@_[0]}{"_value"} eq "-Infinity") ||
          (${@_[0]}{"_value"} eq "Infeasible") ) ) ) {
     ${@_[0]}{"_value"} = eval(${@_[0]}{"_value"});
  }
  #
  # Create auxillary options
  #
  $expoptions = "";
  while (@auxlist) {
    local $tmp = shift @auxlist;
    local @ttmp = split(/=/,$tmp);
    if (substr($ttmp[0],0,1) ne "_") {
       if (substr($ttmp[1],0,1) eq "_") {
           $expoptions .= $ttmp[0] . "=" . ${@_[0]}{$ttmp[1]} . " ";
       } else {
           $expoptions .= $tmp . " ";
       }
    }
  }
  return $expoptions;
}


##
## get tokens from a string
##
sub get_tokens {
  my $__tmp = shift @_;
  local @__words = split(/\s+/, $__tmp, length($__tmp));
  local $__lower = 0;
  local $__upper = $#__words;
  if ($__words[0] eq "") {
     $__lower = 1;
     }
  if ($__words[$__upper] eq "") {
     $__upper = $__upper-1;
     }
  return @__words[$__lower .. $__upper];
}

##
## List-of-strings membership test
##
sub is_string_member {
  my $value = shift(@_);
  my $set = shift(@_);
  local $item;
  foreach $item (@$set) {
    if ($verbose) {
       print "$value $item\n";
       }
    if ($value eq $item) { return 1; }
    }
  return 0;
  }

##
## List-of-strings prefix test
##
sub is_string_prefix {
  my $value = shift(@_);
  my $set = shift(@_);
  local $item;
  foreach $item (@$set) {
    if ($verbose) {
       print "$value $item\n";
       }
    if ($value =~ /^$item/) { return 1; }
    }
  return 0;
  }

##
## Update the configinfo variable
##
sub DEPRICATE_update_key {
  if (!(-r "Makefile")) {
     die "ERROR: Missing Makefile!\n";
     }
  local $tmp = `grep "target_canonical =" Makefile`;
  if ($tmp eq "") { $configinfo = "unknown"; }
  else            { $configinfo = ( get_tokens($tmp) ) [-1]; }
  chomp($configinfo);
  local $tmp = `grep "^CCC" Makefile`;
  if ($tmp eq "") { $tmp = "unknown"; }
  else { $tmp = ( get_tokens($tmp) ) [-1]; }
  $configinfo = "$configinfo-$tmp";
  chomp($configinfo);
}

$misc_key_info="";
$scenario_key="";

sub setup_key {
  (open OUTFILE, ">$_[0]") || die "Error opening file $_[0]";
  #
  # Setup the timestamp
  #
  my ($sec,$min,$hr,$mday,$mon,$year,$dts,$foo);
  ($sec,$min,$hr,$mday,$mon,$year,$foo,$foo,$foo) = localtime(time);

  $mon++;			# 1-12 instead of 0-11
  $year += 1900;		# full year

  $_[1] = sprintf("%04d-%02d-%02d %02d:%02d:%02d",
					$year, $mon, $mday, $hr, $min, $sec);
  #
  # Print the key
  #
  print_key_core(\*OUTFILE, $_[1]);
  close OUTFILE;
  #
  # Setup the scenario variable
  #
  $scenario_key = `cat $_[0]`;
  chomp($scenario_key);
  #
  # Setup the prefix variable
  #
  local $prefix = sprintf("%04d%02d%02d#%02d%02d%02d",
					$year, $mon, $mday, $hr, $min, $sec);
  $prefix = $prefix . "#" . `whoami`;
  chop($prefix);
  $prefix = $prefix . "#" . `hostname`;
  chop($prefix);
  $_[2] = $prefix;
}

#
# Print experiment key
#
sub print_key_core {
  my $OUTFILE = shift;
  my $timestamp = shift;

  my $kernal_name=`uname -s`;
  chomp($kernal_name);
  my $nodename=`uname -n`;
  chomp($nodename);
  my $kernal_release=`uname -r`;
  chomp($kernal_release);
  my $kernal_version=`uname -v`;
  chomp($kernal_version);
  my $machine=`uname -m`;
  chomp($machine);
  my $processor=`uname -p`;
  chomp($processor);
  #my $platform=`uname -i`;
  #chomp($platform);
  my $os=`uname -a | awk '{print $1;}'`;
  chomp($os);
  print $OUTFILE "  <Key\n";
  print $OUTFILE "      KernelName=\"$kernal_name\"\n";
  print $OUTFILE "      HostName=\"$nodename\"\n";
  print $OUTFILE "      KernelRelease=\"$kernal_release\"\n";
  print $OUTFILE "      KernelVersion=\"$kernal_version\"\n";
  print $OUTFILE "      Machine=\"$machine\"\n";
  print $OUTFILE "      Processor=\"$processor\"\n";
  #print $OUTFILE "      Platform=\"$platform\"\n";
  print $OUTFILE "      OS=\"$os\"\n";
  if (!($ENV{'RUNEXP_SCENARIO'} eq "")) {
     print $OUTFILE "      Scenario=\"$ENV{'RUNEXP_SCENARIO'}\"\n";
  } else {
     print $OUTFILE "      Scenario=\"$software\"\n";
  }
  if (!($ENV{'RUNEXP_KEYINFO'} eq "")) {
     print $OUTFILE "      $ENV{'RUNEXP_KEYINFO'}\n";
  } elsif (!($misc_key_info eq "")) {
     print $OUTFILE "      $misc_key_info\n";
  }
  print $OUTFILE "      Date=\"$timestamp\"\n";
  print $OUTFILE "  />\n";
}

##
## Print experiment key, using the scenario key file if available
##
sub print_key {
  my $OUTFILE = shift;
  if (!($scenario_key eq "")) {
     print $OUTFILE "$scenario_key\n";
  } else {
     $timestamp = get_timestamp();
     print_key_core(\*$OUTFILE,$timestamp);
  }
}

##
## Returns a timestamp
##
sub get_timestamp {
   my ($sec,$min,$hr,$mday,$mon,$year,$dts,$foo);
   ($sec,$min,$hr,$mday,$mon,$year,$foo,$foo,$foo) = localtime(time);

   $mon++;			# 1-12 instead of 0-11
   $year += 1900;		# full year

   $dts = sprintf("%04d-%02d-%02d %02d:%02d:%02d",
		$year, $mon, $mday, $hr, $min, $sec);
   return $dts;
}

##
## Print xml date string for today
##
sub DEPRICATED_print_date {
  my $OUTFILE = shift;
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  printf $OUTFILE "  <Date>%4d/%02d/%02d</Date>\n",$year+1900,$mon + 1,$mday;
}

##
## Generate canonical description format
##
sub generate_description {
  my $ctr    = shift(@_);
  $i=0;
  while ($i < $ctr) {
    shift(@_);
    $i = $i + 1;
  }
  my $description = "  <Description>\n";
     #$description .= "    <ExperimentName>" . shift(@_) . "</ExperimentName>\n";
     #$description .= "    <Name>" . shift(@_) . "</Name>\n";
  shift(@_);	# The experiment name
  shift(@_);	# The name
  shift(@_);	# The number of replications
  while (@_) {
    local $factor_name = shift(@_);
    local $id = shift(@_);
    local $tmp = shift(@_);
    scrub(\$tmp);
    $description .= "    <Factor name=\"" . $factor_name . "\" level=\"" . $id . "\">" . $tmp . "</Factor>\n";
    }
  $description .= "  </Description>\n";
  return $description;
}


##
## Scan a file for lines that contain warning or error information.
## XML output is sent to $OUTFILE that contains the status of this 
## test.
##
sub scan_file {
  my $OUTFILE = shift(@_);
  my $search = shift(@_);
  my $warning = shift(@_);
  my $ignore = shift(@_);
  my $infile = shift(@_);
  #
  # Read input file until EOF
  #
  (open INPUT_FILE, $infile) || die "ERROR: $!! Cannot open file \"$infile\"\n" ;
  @lines = <INPUT_FILE>;
  close (INPUT_FILE);
  $line = 0;
  while ($line < $#lines) {
    chomp($lines[$line]);
    scrub(\$lines[$line]);
    $line += 1;
    }
  $found = 0;
  #
  # read each line of file and scan for errors
  #
  $line = 0;
  while ($line < $#lines) {
    $line += 1;
    if ($lines[$line-1] =~ m/$search/) {
       #
       # This is an error line
       #
      if ($lines[$line-1] =~ m/$ignore/) {
        ; # Ignore an error
      }
      else { # actual build error
        if ($found == 0) {
           $found = 1;
           print $OUTFILE "  <IntegrityStatus>Fail</IntegrityStatus>\n";
           print $OUTFILE "  <IntegrityFailureInfo>\n";
        }
        print $OUTFILE "    <Explanation line=\"$line\">\n";
	print $OUTFILE "      <Text>";
        print $OUTFILE $lines[$line-1];
	print $OUTFILE "</Text>\n";
        print $OUTFILE "      <PreContext>";
        if (($line-6) >= 0) { print $OUTFILE "$lines[$line-6];\n" }
        if (($line-5) >= 0) { print $OUTFILE "$lines[$line-5];\n" }
        if (($line-4) >= 0) { print $OUTFILE "$lines[$line-4];\n" }
        if (($line-3) >= 0) { print $OUTFILE "$lines[$line-3];\n" }
        if (($line-2) >= 0) { print $OUTFILE "$lines[$line-2];\n" }
        print $OUTFILE "</PreContext>\n";
        print $OUTFILE "      <PostContext>";
        if (($line+0) < $#lines) { print $OUTFILE "$lines[$line+0];\n" }
        if (($line+1) < $#lines) { print $OUTFILE "$lines[$line+1];\n" }
        if (($line+2) < $#lines) { print $OUTFILE "$lines[$line+2];\n" }
        if (($line+3) < $#lines) { print $OUTFILE "$lines[$line+3];\n" }
        if (($line+4) < $#lines) { print $OUTFILE "$lines[$line+4];\n" }
        print $OUTFILE "</PostContext>\n";
        print $OUTFILE "    </Explanation>\n";
      }
    }
  }
  if ($found == 1) {
     print $OUTFILE "  </IntegrityFailureInfo>\n";
  } else {
     print $OUTFILE "  <IntegrityStatus>Pass</IntegrityStatus>\n";
  }
  #
  # read each line of file and scan for warnings
  #
  print $OUTFILE "  <Warnings>\n";
  $line = 0;
  $warnings = 0;
  $warnings_total = 0;
  $package = "packages";
  my %warning_info;
  while ($line < $#lines) {
    $line += 1;
    if ($lines[$line-1] =~ m/=== Running/) {
       #
       # We are starting to build another package
       #
       if ($package ne "packages") {
	  $warning_info{$package} += $warnings;
	  $warnings_total += $warnings;
          }
       $package = (split(/ /,$lines[$line-1]))[4];
       $warnings = 0;
       }
    if ($lines[$line-1] =~ m/$warning/) {
       #
       # This is a warning line
       #
       $warnings += 1;
       print $OUTFILE "    <Explanation line=\"$line\">\n";
       print $OUTFILE "      <Text>";
       print $OUTFILE $lines[$line-1];
       print $OUTFILE "</Text>\n";
       print $OUTFILE "      <PreContext>";
       if (($line-6) >= 0) { print $OUTFILE "$lines[$line-6];\n" }
       if (($line-5) >= 0) { print $OUTFILE "$lines[$line-5];\n" }
       if (($line-4) >= 0) { print $OUTFILE "$lines[$line-4];\n" }
       if (($line-3) >= 0) { print $OUTFILE "$lines[$line-3];\n" }
       if (($line-2) >= 0) { print $OUTFILE "$lines[$line-2];\n" }
       print $OUTFILE "</PreContext>\n";
       print $OUTFILE "      <PostContext>";
       if (($line+0) < $#lines) { print $OUTFILE "$lines[$line+0];\n" }
       if (($line+1) < $#lines) { print $OUTFILE "$lines[$line+1];\n" }
       if (($line+2) < $#lines) { print $OUTFILE "$lines[$line+2];\n" }
       if (($line+3) < $#lines) { print $OUTFILE "$lines[$line+3];\n" }
       if (($line+4) < $#lines) { print $OUTFILE "$lines[$line+4];\n" }
       print $OUTFILE "</PostContext>\n";
       print $OUTFILE "    </Explanation>\n";
       }
  }
  print $OUTFILE "    <Total>$warnings_total</Total>\n";
  $i = 0;
  print $OUTFILE "    <Summary>\n";
  while (($i, $j) = each(%warning_info)) {
    print $OUTFILE "      <Package name=\"$i\">$j</Package>\n";
    }
  print $OUTFILE "    </Summary>\n";
  print $OUTFILE "  </Warnings>\n";
}

##
## Perform configuration
##
sub exec_configure {
  my $outdir = shift(@_);
  $status = 0;
  $logfile = "";
  @config_options = ();
  foreach $val (@_) {
    if (substr($val,0,10) eq "--logfile=") {
      $logfile = substr($val,10,length($val)-10);
    } elsif (substr($val,0,9) eq "--option=") {
      $value = substr($val,9,length($val)-9);
      $val = config_value($value);
      push @config_options, $val;
    } else {
      push @config_options, $val;
    }
  }
  $t0 = mytime();
  if ($logfile eq "") {
     unlink "${outdir}/config.out", "${outdir}/config.xml";
     if (@config_options) {
        $status = 0xffff & system("(which automake; which autoconf; which libtool; autoreconf --install ; ./configure $siteconfig @config_options) > ${outdir}/config.out 2>&1");
     } else {
        $status = 0xffff & system("(which automake; which autoconf; which libtool; autoreconf --install ; ./configure $siteconfig) > ${outdir}/config.out 2>&1");
     }
     if (-e "../test/scenario_key.txt") {
        $scenario_key = `cat ../test/scenario_key.txt`;
        chomp($scenario_key);
     } elsif (-e "bin/software_info.txt") {
        $misc_key_info = `cat bin/software_info.txt`;
        chomp($misc_key_info);
     }
     open LOGFILE, ">${outdir}/config.tmp";
     print LOGFILE "<html><h2>Acro Config Output</h2>\n";
     print LOGFILE "$host ";
     my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
     printf LOGFILE "%4d/%02d/%02d",$year+1900,$mon + 1,$mday;
     print LOGFILE " ($configinfo)<br>\n$software @config_options";
     print LOGFILE "<hr><pre>\n";
     close LOGFILE;
     system("cat ${outdir}/config.out >> ${outdir}/config.tmp");
     open LOGFILE, ">>${outdir}/config.tmp";
     print LOGFILE "</pre></html>\n";
     close LOGFILE;
     system("mv ${outdir}/config.tmp ${outdir}/config.out");
  }
  $elapsed = time_diff($t0, mytime());
  #
  # See if configure executed properly
  #
  if ($status != 0) {
     #
     # Test failed to execute
     #
     $executionstatus   = "Fail";
     $explanation 	= "Configure failed to execute";
  } else {
     $executionstatus   = "Pass";
     $explanation 	= "";
  }
  #
  # Print the XML summary for this configuration
  #
  (open OUTPUT_FILE, ">${outdir}/config.xml")
	|| die "ERROR: cannot open file ${outdir}/config.xml";
  print OUTPUT_FILE "<Configure>\n";
  print_key(\*OUTPUT_FILE);
  #($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  #$hourfrac = sprintf "%05.2f", $hour + ($min + $sec/60.0)/60.0;
  print OUTPUT_FILE "  <Flags>@config_options</Flags>\n";
  print OUTPUT_FILE "  <LogFile>${outdir}/config.out</LogFile>\n";
  print OUTPUT_FILE "  <StartTime>" . get_timestamp() . "</StartTime>\n";
  print OUTPUT_FILE "  <RunTime unit=\"seconds\">$elapsed</RunTime>\n";
  if ($executionstatus eq "Fail") {
     print OUTPUT_FILE "  <ExecutionStatus>$executionstatus</ExecutionStatus>\n";
     print OUTPUT_FILE "  <ExecutionFailureInfo>\n";
     scrub(\$explanation);
     print OUTPUT_FILE "    <Explanation>$explanation</Explanation>\n";
     print OUTPUT_FILE "  </ExecutionFailureInfo>\n";
     }
  else {
     print OUTPUT_FILE "  <ExecutionStatus>$executionstatus</ExecutionStatus>\n";
     #
     # Define strings to search for
     #
     $search = "^\\*\\*\\*| No rule to make | No such file | Error:|configure: error:";
     $warning = "WARNING";
     #
     # Define strings to ignore
     #
     $ignore = "ignored|^diff:";
     #
     # Process input file
     #
     scan_file(\*OUTPUT_FILE, $search, $warning, $ignore, "${outdir}/config.out")
     }
  print OUTPUT_FILE "</Configure>\n";
  close (OUTPUT_FILE);
}

##
## Perform reconfiguration
##
sub exec_reconfigure {
  my $outdir = shift(@_);
  $status = 0;
  $logfile = "";
  foreach $val (@_) {
    if (substr($val,0,10) eq "--logfile=") {
      $logfile = substr($val,10,length($val)-10);
    } else {
      print "ERROR: unknown option for 'reconfigure': $val";
      exit(1);
    }
  }
  $t0 = mytime();
  if ($logfile eq "") {
     unlink "${outdir}/config.out", "${outdir}/config.xml";
     if (@config_options) {
        $status = 0xffff & system("(which automake; which autoconf; which libtool; autoreconf --install ; ./config.status --recheck ) > ${outdir}/config.out 2>&1");
     } else {
        $status = 0xffff & system("(which automake; which autoconf; which libtool; autoreconf --install ; ./config.status --recheck ) > ${outdir}/config.out 2>&1");
     }
     if (-e "../test/scenario_key.txt") {
        $scenario_key = `cat ../test/scenario_key.txt`;
        chomp($scenario_key);
     } elsif (-e "bin/software_info.txt") {
        $misc_key_info = `cat bin/software_info.txt`;
        chomp($misc_key_info);
     }
     open LOGFILE, ">${outdir}/config.tmp";
     print LOGFILE "<html><h2>Acro Config Output</h2>\n";
     print LOGFILE "$host ";
     my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
     printf LOGFILE "%4d/%02d/%02d",$year+1900,$mon + 1,$mday;
     print LOGFILE " ($configinfo)<br>\n$software @config_options";
     print LOGFILE "<hr><pre>\n";
     close LOGFILE;
     system("cat ${outdir}/config.out >> ${outdir}/config.tmp");
     open LOGFILE, ">>${outdir}/config.tmp";
     print LOGFILE "</pre></html>\n";
     close LOGFILE;
     system("mv ${outdir}/config.tmp ${outdir}/config.out");
  }
  $elapsed = time_diff($t0, mytime());
  #
  # See if configure executed properly
  #
  if ($status != 0) {
     #
     # Test failed to execute
     #
     $executionstatus   = "Fail";
     $explanation 	= "Configure failed to execute";
  } else {
     $executionstatus   = "Pass";
     $explanation 	= "";
  }
  #
  # Print the XML summary for this configuration
  #
  (open OUTPUT_FILE, ">${outdir}/config.xml")
	|| die "ERROR: cannot open file ${outdir}/config.xml";
  print OUTPUT_FILE "<Configure>\n";
  print_key(\*OUTPUT_FILE);
  #($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  #$hourfrac = sprintf "%05.2f", $hour + ($min + $sec/60.0)/60.0;
  print OUTPUT_FILE "  <Flags>@config_options</Flags>\n";
  print OUTPUT_FILE "  <LogFile>${outdir}/config.out</LogFile>\n";
  print OUTPUT_FILE "  <StartTime>" . get_timestamp() . "</StartTime>\n";
  print OUTPUT_FILE "  <RunTime unit=\"seconds\">$elapsed</RunTime>\n";
  if ($executionstatus eq "Fail") {
     print OUTPUT_FILE "  <ExecutionStatus>$executionstatus</ExecutionStatus>\n";
     print OUTPUT_FILE "  <ExecutionFailureInfo>\n";
     scrub(\$explanation);
     print OUTPUT_FILE "    <Explanation>$explanation</Explanation>\n";
     print OUTPUT_FILE "  </ExecutionFailureInfo>\n";
     }
  else {
     print OUTPUT_FILE "  <ExecutionStatus>$executionstatus</ExecutionStatus>\n";
     #
     # Define strings to search for
     #
     $search = "^\\*\\*\\*| No rule to make | No such file | Error:|configure: error:";
     $warning = "WARNING";
     #
     # Define strings to ignore
     #
     $ignore = "ignored|^diff:";
     #
     # Process input file
     #
     scan_file(\*OUTPUT_FILE, $search, $warning, $ignore, "${outdir}/config.out")
     }
  print OUTPUT_FILE "</Configure>\n";
  close (OUTPUT_FILE);
}

##
## Perform build
##
sub exec_build {
  my $outdir = shift(@_);
  $t0 = mytime();
  $status = 0;
  $clean = "";
  $logfile = "";
  foreach $val (@_) {
    if ($val eq "--clean") { $clean = "make clean; "; }
    elsif ($val eq "--check") { $checkflag = "check; "; }
    elsif (substr($val,0,10) eq "--logfile=") { $logfile = substr($val,10,length($val)-10); }
    else
       { die "ERROR: Unknown build option: $val\n"; }
    }
  if ($logfile eq "") {
     unlink "${outdir}/build.out", "${outdir}/build.xml";
     open LOGFILE, ">${outdir}/build.out";
     print LOGFILE "<html><h2>Acro Build Output</h2>\n";
     print LOGFILE "$host ";
     my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
     printf LOGFILE "%4d/%02d/%02d",$year+1900,$mon + 1,$mday;
     print LOGFILE " ($configinfo)<br>\n$software @_";
     print LOGFILE "<hr><pre>\n";
     close LOGFILE;
     $status = 0xffff & system("($clean make $checkflag ) >> ${outdir}/build.out 2>&1");
     $logfile = "${outdir}/build.out";
     open LOGFILE, ">>${outdir}/build.out";
     print LOGFILE "</pre></html>\n";
     close LOGFILE;
     }
  $elapsed = time_diff($t0, mytime());
  #
  # See if make executed properly
  # 
  if ($status != 0) {
     #
     # Test failed to execute
     #
     $executionstatus   = "Fail";
     $explanation 	= "Make failed to execute (status = $status)";
  } else {
     $executionstatus   = "Pass";
     $explanation 	= "";
  }
  if (-e "../test/scenario_key.txt") {
     $scenario_key = `cat ../test/scenario_key.txt`;
     chomp($scenario_key);
  } elsif (-e "bin/software_info.txt") {
     $misc_key_info = `cat bin/software_info.txt`;
     chomp($misc_key_info);
  }
  (open OUTPUT_FILE, ">${outdir}/build.xml") || die "ERROR: cannot open file ${outdir}/build.xml";
  print OUTPUT_FILE "<Build>\n";
  print_key(\*OUTPUT_FILE);
  print OUTPUT_FILE "  <Flags>@_</Flags>\n";
  print OUTPUT_FILE "  <LogFile>$logfile</LogFile>\n";
  print OUTPUT_FILE "  <StartTime>" . get_timestamp() . "</StartTime>\n";
  print OUTPUT_FILE "  <RunTime unit=\"seconds\">$elapsed</RunTime>\n";
  print OUTPUT_FILE "  <ExecutionStatus>$executionstatus</ExecutionStatus>\n";
  if ($executionstatus eq "Fail") {
     print OUTPUT_FILE "  <ExecutionFailureInfo>\n";
     scrub(\$explanation);
     print OUTPUT_FILE "    <Explanation>$explanation</Explanation>\n";
     print OUTPUT_FILE "  </ExecutionFailureInfo>\n";
  }
  #
  # Define strings to search for
  #
  $search = "ERROR |\tERROR| Error |\tError| Error:| error:";
  $warning = " warning:| Warning:";
  #
  # Define strings to ignore
  #
  $ignore = "ignored";
  #
  # Process input file
  #
  scan_file(\*OUTPUT_FILE, $search, $warning, $ignore, $logfile);
  print OUTPUT_FILE "</Build>\n";
  close (OUTPUT_FILE);
}

#
# Print summary of test results
#
sub print_experiment_summary {
  my $XMLFILE = shift @_;
  my $category = shift @_;
  my $name = shift @_;
  my $executionstatus = shift @_;
  my $explanation = shift @_;
  my $ntrials = shift @_;
  my $cmdline = shift @_;
  my $measurementSub = shift @_;
  $outname = $category . "-" . $name;
  ##
  ## Print out results
  ##
  print $XMLFILE "<Experiment LogFile=\"$outname\">\n";
  print $XMLFILE "  <Category>$category</Category>\n";
  print $XMLFILE "  <Name>$name</Name>\n";
  print_key(\*$XMLFILE);
  #($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  #$hourfrac = sprintf "%05.2f", $hour + ($min + $sec/60.0)/60.0;
  print $XMLFILE "  <StartTime>" . get_timestamp() . "</StartTime>\n";
  print $XMLFILE "  <RunTime unit=\"seconds\">$elapsed</RunTime>\n";
  print $XMLFILE "  <CommandLine>$cmdline</CommandLine>\n";
  if ($executionstatus eq "Fail") {
     #
     # The testing scripts had some sort of problem
     #
     print $XMLFILE "  <ExecutionStatus>$executionstatus</ExecutionStatus>\n";
     print $XMLFILE "  <ExecutionFailureInfo>\n";
     scrub(\$explanation);
     print $XMLFILE "    <Explanation>$explanation</Explanation>\n";
     print $XMLFILE "  </ExecutionFailureInfo>\n";

  } else {
     print $XMLFILE "  <ExecutionStatus>$executionstatus</ExecutionStatus>\n";
     &$measurementSub($XMLFILE);
     }
  print $XMLFILE "</Experiment>\n";
}

##
## A dummy function
##
sub dummyMeasurementSub {
}


##
## Execute and process a performance test
##
sub run_experiments {
  my $XMLFILE = shift @_;
  my $category = shift @_;
  my $name = shift @_;
  my $cmdline = shift @_;
  my $ntrials = shift @_;
  my $analysisSub = shift @_;
  my $measurementSub = shift @_;
  my $expoptions = shift @_;
  my(%options) = %{(shift)};

  #
  # Get scenario key
  #
  if (-e "../../test/scenario_key.txt") {
     $scenario_key = `cat ../../test/scenario_key.txt`;
     chomp($scenario_key);
  }
  $executionstatus = "Pass";
  $explanation = "";
  #
  # Get seeds
  #
  $seedfile = $category . ".exp.seeds";
  #print "SEEDFILE $seedfile\n";
  @seeds = ();
  if (-r $seedfile) {
     open (INPUT, $seedfile) || die "ERROR: cannot open file \"$seedfile\"!\n";
     @lines = <INPUT>;
     close(INPUT);
     foreach $line (@lines) {
       chomp($line);
       push @seeds, $line;
     }
  }
  #
  # Run experiment
  #
  $i = 1;
  while (($i <= $ntrials) && ($executionstatus eq "Pass")) {
    $logname = $category . "-" . $name . "." . $i . ".log";
    $measurementname = $category . "-" . $name . "." . $i . ".out";
    $tmp_cmdline = sprintf($cmdline, $measurementname, $logname);
    #
    # Setup random number seed
    #
    if ($#seeds >= 0) {
       $ENV{'PSEUDORANDOM_SEED'} = $seeds[$i-1];
       $seed = $seeds[$i-1];
    } else {
       $seed = POSIX::ceil(rand()*1000000);
       $ENV{'PSEUDORANDOM_SEED'} = $seed;
    }
    #
    # Setup header for logfile
    #
    open (OUTPUT, ">${logname}.xml");
    print OUTPUT "<TrialOutput>\n";
    print OUTPUT "<TrialOptions>\n";
    print OUTPUT "  $expoptions";
    print OUTPUT "</TrialOptions>\n";
    print OUTPUT "<CommandLine>\n";
    print OUTPUT "$tmp_cmdline (\$seed = $seed)\n";
    print OUTPUT "</CommandLine>\n";
    print OUTPUT "<CommandLog>\n";
    close(OUTPUT);
    #
    # Launch the system command
    #
    $lcmd = "\$status = 0xffff & system(\"( testlib_exec $tmp_cmdline ) >> ${logname}.xml 2>&1\")";
    $t0 = mytime();
    eval($lcmd);
    $elapsed = time_diff($t0, mytime());
    #
    # See if make executed properly
    #
    if ($status != 0) {
       #
       # Test failed to execute
       #
       $executionstatus   = "Fail";
       $explanation 	= "Failed execution";
    } else {
       #
       # Test executed properly
       #
       if (! (-e "$logname")) {
          #
          # The test log file appears to be missing
          #
          $executionstatus  = "Fail";
          $explanation = "Missing log file $logname";
       } elsif ((get_tokens(`wc $logname`))[1] eq "0") {
          #
          # Check for an empty file
          #
          $executionstatus  = "Fail";
          $explanation = "Empty log file";
       } else {
          open (OUTPUT, ">>$logname");
          print OUTPUT "Seed: $seed\n";
          close (OUTPUT);
          &$analysisSub("$measurementname","$logname",$seed,$executionstatus);
          if ($executionstatus eq "Fail") {
             $explanation = "Bad log format in file $logname";
          }
          open (OUTPUT, ">>${logname}.xml");
          print OUTPUT "</CommandLog>\n";
          print OUTPUT "<Output>\n";
          `cat ${logname} >> ${logname}.xml`;
          print OUTPUT "</Output>\n";
          print OUTPUT "</TrialOutput>\n";
          close (OUTPUT);
       }
    }
    $i = $i + 1;
  }
  #
  # Print results info
  #
  print_experiment_summary($XMLFILE,$category,$name,$executionstatus,$explanation,$i,$cmdline,$measurementSub);
}

##
## Perform tests
##
## NOTE: we might want to re-architect this to distribute jobs on
## other machines.  For now, I'm assuming that all jobs run on the
## current machine.
##
sub exec_tests {
  my $outdir = shift(@_);
  #
  # Process test arguements
  #
  @tags = ();
  $coverage=0;
  $covdir="lcov";
  while (@_) {
    $arg = shift @_;
    if (substr($arg,0,5) eq "--tag") {
       @tags = split(/,/, (split(/=/,$arg,2))[1], $#arg);
    } elsif (substr($arg,0,5) eq "--cov") {
       $coverage=1;
       if ($arg =~ /=(.+)/) {
         $covdir = $1;
       }
# make sure the coverage directory exists and is cleaned up
       system("mkdir test/lcov");
       unlink "test/lcov/acro.info";
       unlink "test/lcov/final.info";
       unlink "test/lcov/lcov.out";
       unlink "test/lcov/genhtml.out";
       system("find packages -name '*.gcda' | xargs --no-run-if-empty /bin/rm");
       system("mkdir ~/public_html/$covdir");   # expected to run on webserver -syc
    } else {
      die "ERROR: Bad test argument: \"$arg\"\n";
    }
  }
  $testtags = "";
  foreach $tag (@tags) {
    $testtags .= " --tag=" . $tag;
  }
  if ($testtags eq "") {
     $testtags = "--tag=smoke";
  }
  #
  # Execute experiments in *.study.xml files
  #
  $result = `(find packages -name '*.study.xml') 2>&1`;
  chomp($result);
  foreach $testfile (get_tokens($result)) {
    local @tmp = split(/\//,$testfile);
    $testname = $tmp[-1];
    $testname = substr($testname,0,$#testname-7);
    #print "      Considering category \"$category\"\n";
    local $dir = join( '/', @tmp[0 .. ($#tmp - 1)]);
    chdir $dir or die "ERROR: Cannot move to directory $dir!\n";
    system("(python ../../../bin/exact --scenario=$software $testtags $tmp[-1]) 2>&1");
    $testname = extract_name($tmp[-1]);
    if (-d $testname) {
       system("(find $testname -name '*.out' -exec /bin/cp {} $outdir \\; 2>&1) > /dev/null");
       system("(find $testname -name '*.log.xml' -exec /bin/cp {} $outdir \\; 2>&1) > /dev/null");
       system("cp $tmp[-1] $outdir");
    }
    chdir $pwd;
  }
  system("(find packages -name '*.results.xml' -exec /bin/cp {} $outdir \\; 2>&1) > /dev/null");
  system("(find packages -name '*.analysis.xml' -exec /bin/cp {} $outdir \\; 2>&1) > /dev/null");
  #
  # Apply coverage method
  #
  if ($coverage) {
    my @dirs = `find packages -name '*.gcda' | xargs -n1 dirname | sort | uniq`;
    print "Gathering coverage data\n";
    for my $d (@dirs) {
      chomp $d;
      print "  $d\n";
      system("lcov -d $d -c -f >> test/lcov/acro.info 2>> test/lcov/lcov.out");
    }
    print "Generating coverage html\n";
    #
    # remove system headerfiles randomly covered
    #
    system("lcov -r test/lcov/acro.info '/usr/include/*' > test/lcov/final.info 2>> test/lcov/lcov.out");
    #
    # genhtml.pl using coverage directory as output (overwrites existing, but keeps
    # any extra files that aren't covered.  Might result in dangling webfiles,
    # but that might be useful to know.
    #
    system("(cd test/lcov; genhtml.pl -o ~/public_html/$covdir -k final.info -t 'Acro' > genhtml.out) 2>&1");
  }
}

##
## Perform code check
##
sub exec_code_checks {
  my $outdir = shift(@_);
  my $status = 0;
  my $check_command="";
  my $name="default";
  foreach $val (@_) {
    if (substr($val,0,7) eq "--exec=") {
	$check_command = substr($val,7,length($val)-7);
    } elsif (substr($val,0,7) eq "--name=") {
	$name = substr($val,7,length($val)-7);
    } else {
	die "ERROR: Unknown check option: $val\n"; }
    }
  if ($name eq "all") {
     foreach $file (glob("qa/checks/*")) {
       $check=(split(/\//,$file))[2];
       if ($check ne "CVS") {
          print "      $check\n";
          exec_code_checks($outdir,"--name=$check");
       }
     }
     return;
  }
  if ($check_command eq "") {
     $check_command = "qa/checks/$name";
     }
  unlink "${outdir}/check-${name}.out", "${outdir}/check-${name}.xml";
  #open LOGFILE, ">${outdir}/${name}.out";
  #print LOGFILE "<html><h2>Acro Code Check Output</h2>\n";
  #print LOGFILE "$host ";
  #my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  #printf LOGFILE "%4d/%02d/%02d",$year+1900,$mon + 1,$mday;
  #print LOGFILE "<hr><pre>\n";
  #close LOGFILE;
  #$status = 0xffff & system("( $check_command ) >> ${outdir}/${name}.out 2>&1");
  #$logfile = "${outdir}/${name}.out";
  #open LOGFILE, ">>${outdir}/${name}.out";
  #print LOGFILE "</pre></html>\n";
  #close LOGFILE;
  #if (-e "../test/scenario_key.txt") {
     #$scenario_key = `cat ../test/scenario_key.txt`;
     #chomp($scenario_key);
  #} elsif (-e "bin/software_info.txt") {
     #$misc_key_info = `cat bin/software_info.txt`;
     #chomp($misc_key_info);
  #}


  #
  # Setup the xml file.
  #
  (open OUTPUT_FILE, ">${outdir}/check-${name}.xml") || die "ERROR: cannot open file ${outdir}/check-${name}.xml";
  print OUTPUT_FILE "<CodeCheck name=\"${name}\">\n";
  #print_key(\*OUTPUT_FILE);
  print OUTPUT_FILE "  <Date>" . get_timestamp() . "</Date>\n";
  #print OUTPUT_FILE "  <Flags>@_</Flags>\n";
  print OUTPUT_FILE "  <LogFile>check-${name}.out</LogFile>\n";
  print OUTPUT_FILE "  <StartTime>" . get_timestamp() . "</StartTime>\n";
  print OUTPUT_FILE "  <CheckInfo>\n";
  close (OUTPUT_FILE);
  #
  # Run the command
  #
  my $t0 = mytime();
  $status = 0xffff & system("( $check_command ${outdir}/check-${name}.out ) >> ${outdir}/check-${name}.xml 2>&1");
  $elapsed = time_diff($t0, mytime());
  #
  # See if make executed properly
  # 
  if ($status != 0) {
     #
     # Test failed to execute
     #
     $executionstatus   = "Fail";
     $explanation 	= "Make failed to execute (status = $status)";
  } else {
     $executionstatus   = "Pass";
     $explanation 	= "";
  }
  (open OUTPUT_FILE, ">>${outdir}/check-${name}.xml") || die "ERROR: cannot open file ${outdir}/check-${name}.xml";
  
  print OUTPUT_FILE "  </CheckInfo>\n";
  print OUTPUT_FILE "  <RunTime unit=\"seconds\">$elapsed</RunTime>\n";
  print OUTPUT_FILE "  <ExecutionStatus>$executionstatus</ExecutionStatus>\n";
  if ($executionstatus eq "Fail") {
     print OUTPUT_FILE "  <ExecutionFailureInfo>\n";
     scrub(\$explanation);
     print OUTPUT_FILE "    <Explanation>$explanation</Explanation>\n";
     print OUTPUT_FILE "  </ExecutionFailureInfo>\n";
     }
  print OUTPUT_FILE "</CodeCheck>\n";
  close (OUTPUT_FILE);
}

##
## Execute a command
##
sub exec_command {
  my $command = shift(@_);
  #
  # Perform configuration operations
  #
  if ("$command" eq "configure") {
    print "   Performing configure\n";
    exec_configure(@_);
  #
  # Perform re-configuration operations
  #
  } elsif ("$command" eq "reconfigure") {
    print "   Performing reconfigure\n";
    exec_reconfigure(@_);
  #
  # Perform build/make operations
  #
  } elsif (("$command" eq "build") || ("$command" eq "make")) {
    print "   Performing build\n";
    exec_build(@_);
  #
  # Perform test operations
  #
  } elsif ("$command" eq "test") {
    print "   Performing tests\n";
    exec_tests(@_);
  #
  # Perform check operations
  #
  } elsif ("$command" eq "check") {
    print "   Performing code checks\n";
    exec_code_checks(@_);
  #
  # ELSE ERROR
  #
  } else {
    die "ERROR: Bad command \"$command\"\n";
    }
}
 

##
## TODO
##
sub DEPRECATED_process_testfile {
  my $testfile = shift @_;
  my $testtype = shift @_;
  my $daily    = shift @_;
  my $weekly   = shift @_;
  my $monthly  = shift @_;
  #
  # Get test category from test filename
  #
  $category = substr($testfile,0,length($testfile)-6);
  print "      Checking $category tests\n";
  #
  # Remove old test output files
  #
  system("rm -f $category.xml $category.*.out $category.*.diff");
  #
  # Run tests
  #
  open (TESTFILE, $testfile) || die "ERROR: cannot open file \"$testfile\"!\n";
  open (XMLFILE, ">$category.tests.xml")
  		|| die "ERROR: cannot open file \"$category.tests.xml\"!\n";
  print XMLFILE  "<ExperimentalResults Category=\"$category\">\n";
  $lineno = 0;
  while (<TESTFILE>) {
    $lineno++;
    #
    # Ignore comment lines
    #
    if (/^[\t ]*#.*/) { next; }
    $line = $_;
    chomp($line);
    @words = get_tokens($line);

    $executeflag = 0;
    $name = shift @words;
    $localtesttype = $testtype;
    $benchmark = "";
    $analysiscode = "";
    #print "$line\n";
    #print "$name\n";

    while (@words) {
      local $tmp = shift @words;
      if (substr($tmp,0,2) ne "--") {
         die "ERROR: bad test option \"$tmp\" at line $lineno in file $testfile. (1)\n";
         }
      #
      # Parse an option and validate it
      #
      $option = (split(/=/,$tmp))[0];
      $option = substr($option,2,length($option));
      $value  = (split(/=/,$tmp))[1];
      #print "$option $value\n";
      if ((is_string_member($option,\@tests) == 0) &&
          (is_string_member($option,\@testoptions) == 0)) {
         die "ERROR: bad test option \"$tmp\" at line $lineno in file $testfile. (2)\n";
         }
      #
      # Process option
      #
      if (($option eq $testtype) || 
          (($executeflag==0) && ($testtype eq "smoke") && (is_string_member($option,\@tests)==1))) {
         #
         # The option is the type of command we're executing
         #                        OR
         # The we're processing smoke tests ... so we consider all
         # command options
         #
         local @tmp = split /,/, $value;
	 foreach $val (@tmp) {
           #print "VAL: $val @tmp $value\n";
	   if (is_string_member($val,\@test_frequencies) == 0)
              { die "ERROR: bad test frequency: $val\n"; }
	   if (($testtype eq "smoke") && ($val eq "smoke")) {
		$localtesttype = $option;
		$executeflag = 1; last;
	   } elsif (($testtype ne "smoke") &&
	  		 (($val eq "all")
			 || ($daily && ($val eq "daily"))
			 || ($weekly && ($val eq "weekly"))
			 || ($monthly && ($val eq "monthly")))) {
		$executeflag = 1; last;
	   }
	 }
      }
      if ($option eq "benchmark") {
	 $benchmark = $value;
      } elsif ($option eq "analysiscode") {
	 $analysiscode = $value;
      } elsif ($option eq "software") {
         $flag = 0;
         local @tmp = split /,/, $value;
	 foreach $val (@tmp) {
           if ($software eq $val) {
	      $flag = 1;
	      last;
	      }
           }
	 if ($flag == 0) {
            $executeflag = 0;
	    goto NEXTLINE;
            }
      }
    }
    #
    # Make sure we have what we need to execute the command
    #
    if (($localtesttype eq "regression") && ($benchmark eq "")) {
       die "ERROR: missing benchmark for test \"$name\"\n";
    } elsif (($localtesttype eq "performance") && ($analysiscode eq "")) {
       die "ERROR: missing analysiscode for test \"$name\"\n";
    }
    NEXTLINE:
    $line = <TESTFILE>;
    $lineno++;
    #
    # Execute the command
    #
    if ($executeflag == 1) {
       run_test(\*XMLFILE,$category, $name, $line, $localtesttype,
						$benchmark, $analysiscode);
    }
  }
  print XMLFILE "</ExperimentalResults>\n";
}

eval 1;
