#!/usr/bin/perl
# This script was originally developed to get Intel compilers to configure
# Trilinos.  The major problem addressed was the format required by Intel 
# compilers to name the object file, if the name was different than the source
# file.  This perl script parses the command line and converts only the 
# necessary arguments, passing the rest through.  In fact, the compiler is
# an argument to this script, so it can be used elsewhere if needed.
#
# Author : Heidi Thornquist
# Date : April 7, 2003
#

$output = "";
$i = 0;
$numinputs = @ARGV;
while ($i < $numinputs) {
	# First things first, fix the -o problem!
	if (index($ARGV[$i],"-o")==0 && index($ARGV[$i+1],"exe")<0) {
		$output = $output . "-Fo" . $ARGV[$i+1] . " ";
		$i++;
	# Whoops, we have a little problem with windows paths :)
	} elsif (index($ARGV[$i],"\\")>=0) {
		$wherebeg=0; $where=0;	
		$newfilename = "";
		while ($where >=0) {
	        	$where = index($ARGV[$i],"\\", $wherebeg );
			if ($where < 0) {
				$part = substr($ARGV[$i], $wherebeg, 100);
				$newfilename = $newfilename . $part;
			} else {
				$part = substr($ARGV[$i], $wherebeg, $where-$wherebeg);
				$newfilename = $newfilename . $part . "\\\\";
			}
			$wherebeg = $where + 1;
		}
		$output = $output . $newfilename . " ";
	# Otherwise, just pass the argument through unchanged.
	} else {
		$output = $output . $ARGV[$i] . " ";
	}
	$i++;
}
# Let the user know what the command line will be, then execute it.
print ("$output","\n");
exec ("$output") || die "Cannot run reformatted compiler line";

