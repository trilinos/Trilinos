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
# Revised : March 31, 2005
# Note:  This script was updated to be compatable with Intel v8 compilers.

$output = "";
$i = 0;
$numinputs = @ARGV;
while ($i < $numinputs) {
  my $currArg = $ARGV[$i];
  if (index($ARGV[$i],"-o")==0 && index($ARGV[$i+1],"exe")<0) {
    # First things first, fix the -o problem!
    $output = $output . "-Fo" . $ARGV[$i+1] . " ";
    $i++;
  } elsif ($currArg eq "-MT") {
    # Fix the depandancy generation option
    $output = $output . "-QM";
    $i++;
  } elsif ($currArg eq "-MD" || $currArg eq "-MP" || $currArg eq "-MF" ) {
    # Ignore these options?
    $i++;
  }  elsif (index($ARGV[$i],"\\")>=0) {
    # Whoops, we have a little problem with windows paths :)
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
  } elsif (index($ARGV[$i],"-l")==0) {
    # Fix the -l problem
    $newlib = "";
    $arg_length = length($ARGV[$i]);
    if ($arg_length > 2) {
      # Check rest of argument, library is tacked on (ex. -lteuchos)
      $templib = substr($ARGV[$i], 2, 100);
      if (index($templib,".lib")>=0) {
        $newlib = $templib;
      } else {
        $newlib = "lib" . $templib . ".a"; # lib<library_name>.a      
        #$newlib = $newlib . " lib" . $templib . ".so";  # lib<library_name>.so      
      }
    } else {
      $newlib = $ARGV[$i+1];
      $i++; 
    }
    $output = $output . " " . $newlib . " ";
  } elsif (index($ARGV[$i],"-L")==0) {
    # Fix the -L problem
    $newlibpath = "";
    $pathpart = substr($ARGV[$i], 2, 100);
    # Check to see if this is a cygwin path (ex. /cygdrive/c/Trilinos/...)
    $wherebeg=index($ARGV[$i],"/cygdrive");
    if ($wherebeg >= 0) {
      # Grab the name of the disk, which is expected to be the next directory
      # after cygdrive (ex. /cygdrive/c/Trilinos/..., "c" is the disk name)
      $part = "";
      $where=index($ARGV[$i],"/",$wherebeg+10);
      $newlibpath=$newlibpath . substr($ARGV[$i], $wherebeg+10, $where-$wherebeg-10);
      $wherebeg=$where + 1;

      # Add colon and backslashes before appending directories
      $newlibpath= $newlibpath . ":\\\\";

      # Find directories and insert in windows' style path
      while ($where >=0) {
        $where = index($ARGV[$i],"/",$wherebeg);
        if ($where < 0) {
          $part = substr($ARGV[$i], $wherebeg, 100);
          $newlibpath = $newlibpath . $part;
        } else {
          $part = substr($ARGV[$i], $wherebeg, $where-$wherebeg);
          $newlibpath = $newlibpath . $part . "\\\\";
        }
        $wherebeg = $where + 1;
      }
    } else {
      # Assume already have windows path
      $newlibpath = $pathpart;
    }
    $output = $output . "/link /libpath:" . $newlibpath . " ";
  } elsif (index($ARGV[$i],"-g")==0) {
    # Fix the -g problem
    #Do nothing for now -g only generates debugging information for GDB
    # Fix the problem when absolute cygwin paths are used (-I/cygdrive/c/...)
  } elsif (index($ARGV[$i],"/cygdrive")>=0) {
    #Grab the part of the argument before "cygdrive" and preserve it.
    $wherebeg=index($ARGV[$i],"/cygdrive");
    $newpathname=substr($ARGV[$i],0,$wherebeg);
    # Grab the name of the disk, which is expected to be the next directory
    # after cygdrive (ex. /cygdrive/c/Trilinos/..., "c" is the disk name)
    $part = "";
    $where=index($ARGV[$i],"/",$wherebeg+10);
    $newpathname=$newpathname . substr($ARGV[$i], $wherebeg+10, $where-$wherebeg-10);
    $wherebeg=$where + 1;
    # Add colon and backslashes before appending directories
    $newpathname= $newpathname . ":\\\\";
    # Find directories and insert in windows' style path
    while ($where >=0) {
      $where = index($ARGV[$i],"/",$wherebeg);
      if ($where < 0) {
        $part = substr($ARGV[$i], $wherebeg, 100);
        $newpathname = $newpathname . $part;
      } else {
        $part = substr($ARGV[$i], $wherebeg, $where-$wherebeg);
        $newpathname = $newpathname . $part . "\\\\";
      }
      $wherebeg = $where + 1;
    }
    $output = $output . $newpathname . " ";
    #print("$ARGV[$i]","\t","$newpathname","\n");
    #Otherwise, just pass the argument through unchanged.
  } else {
    $output = $output . $ARGV[$i] . " ";
  }
  $i++;
}
# Let the user know what the command line will be, then execute it.
print ("$output","\n");
exec ("$output") || die "Cannot run reformatted compiler line";

