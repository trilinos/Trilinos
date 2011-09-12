#!/usr/bin/perl -w

# @HEADER
# ************************************************************************
#
#        Phalanx: A Partial Differential Equation Field Evaluation
#       Kernel for Flexible Management of Complex Dependency Chains
#                    Copyright 2008 Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
# National Laboratories.
#
# ************************************************************************
# @HEADER

#
# This perl script replaces a string with another string.
# Here it is allowd for file_in and file_out to be the
# same file.
#
use strict;
#
my $g_use_msg =
  "Use: string-replace.pl find_string replacement_string file_in file_out\n";
if( scalar(@ARGV) < 3 ) {
  print STDERR $g_use_msg;
  exit(-1);
}
#
my $find_string        = shift;
my $replacement_string = shift;
my $file_in_name       = shift;
my $file_out_name      =$file_in_name;
##my $file_out_name      = shift;
#
#
if($file_in_name=~/CVS/) {
#  print "Do not replace in CVS\n";
  exit;
}
#
open FILE_IN, "<$file_in_name" || die "The file $file_in_name could not be opended for input\n";
my @file_in_array = <FILE_IN>;
close FILE_IN;
#
my @file_out_array;
my $did_replacement = 0;
foreach(@file_in_array) {
  #print $_;
  $did_replacement = 1 if $_=~s/$find_string/$replacement_string/g;
  #print $_;
  push @file_out_array, $_;
}
if($did_replacement || $file_out_name ne $file_in_name) {
  open FILE_OUT, ">$file_out_name" || die "The file $file_out_name could not be opended for output\n";
  print FILE_OUT @file_out_array;
  close FILE_OUT;
}
