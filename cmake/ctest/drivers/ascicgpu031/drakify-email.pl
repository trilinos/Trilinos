#!/usr/bin/perl -w

# $file = "./sample-email.txt";
# 
# open(FILE, "<$file");
# @LINES = <FILE>;
# close(FILE);

my @LINES;
while (<STDIN>)
{
  push @LINES, $_;
}

$numLines = scalar @LINES;

# If you want to reuse this script somewhere else,
# this should be the only line you need to change
@packages=("Tpetra", "Amesos2","Belos","Ifpack2","Xpetra","MueLu","Zoltan2","Thyra");

$num_packages = scalar @packages;

use Class::Struct;
struct Entry => {
  name => '$',  
  BuildWarnings => '@',
  BuildErrors => '@',
  TestPasses => '@',  
  TestFailures => '@',  
  LibBuildWarnings => '@',
  LibBuildErrors => '@',  
};

sub extractLeadingNumber
{
  my $str = shift;
#   print "str is $str.\n";
  # maybe not the best way to do this, but this is what I know how to do without Google's help
  my @fields = split /\s+/, $str;
  return $fields[0];
}

sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };

# my $trimTest = "  string with trailing spaces   ";
# my $trimmed = trim $trimTest;
# print "Trimmed string: \"$trimmed\"\n";

my @entries, my @entriesNightly, my @entriesExperimental, my@entriesSpecialized;
my $gitFailure = "";

my %testType;
$mytestType = "Nightly";
my $pastSummary = 0;


for ($i=0; $i<$numLines; $i++)
{
  my $line = $LINES[$i];
  #print $line;
  if ($pastSummary == 0)
  {
#     print "NOT past summary.\n";
    if ($line =~ /git update FAILED/)
    {
      $gitFailure = "<font color=\"red\">*** git update FAILED ***</font>";
    }
    if ($line =~ /\++ Nightly \++/)
    {
      $mytestType = "Nightly";
    }
    elsif ($line =~ /\++ Experimental \++/)
    {
      $mytestType = "Experimental";
    }
    elsif ($line =~ /\++ Specialized \++/)
    {
      $mytestType = "Specialized";
    }
    
    if ($line =~ /(.*)\s+\.\.\. [pF]/)
    {
      $testName = trim $1;
      $testType{$testName} = $mytestType;
      #print "testType[$testName] = $mytestType\n"
      
    }
    else
    {
      #print "line did not match the test listing regex: $line.\n";
    }
    if ($line =~ /-----------------------------------------------------------------------------/)
    {
      $pastSummary = 1;
      #print "FOUND THE END OF THE SUMMARY";
      $i++; # blank line follows the --- line: we ignore this one.
      $i++;
    }
  }
  
  if ($pastSummary == 1)
  {
#entries are of the following form:
#     OPENMPI_1.10.0_RELEASE_DEV_MueLu_NO_SERIAL             
#     passed  ,  14.3 seconds
#              Xpetra |   0 warnings |   0 errors | 0/0 passed
#          MueLu_libs |   0 warnings |   0 errors
#               MueLu |   0 warnings |   0 errors | 0/0 passed
#         Xpetra_libs |   0 warnings |   0 errors
    my $entry = Entry->new();
    my $name = trim $LINES[$i];
#     print "name: $name.\n";
    $entry->name($name); #trim gets rid of any leading/trailing whitespace
    $i++;
    $line = $LINES[$i]; # timing/overall result line -- for now, we ignore
    my $numEntryLines = 2 * $num_packages;
    for (my $j=0; $j<$numEntryLines; $j++) {
        $i++;
        $line = $LINES[$i];
        my @fields = split /\s*\|\s*/, $line;
#       print "field 0: " . $fields[0] . "\n";
#       print "num fields: " . scalar @fields . "\n";
        my $whichTests = trim $fields[0];
        
        $have_match=0;
        for(my $k=0; $k<$num_packages; $k++) {
            $libsstr=$packages[$k] . '_libs';
            if($whichTests =~ /$libsstr/) {
                ### Check for "libs", e.g. MueLu_libs ###
                # get numWarnings from $fields[1] = "\d+ warnings"
                $numWarnings = extractLeadingNumber $fields[1];
                $entry->LibBuildWarnings($k,$numWarnings);
                # get numErrors from $fields[2] = "\d+ errors"
                $numErrors = extractLeadingNumber $fields[2];
                $entry->LibBuildErrors($k,$numErrors);    
                $have_match=1;
            }
            elsif ( $whichTests =~ /$packages[$k]/) {
                ### Regular libraries, e.g. MueLu ###
                # get numWarnings from $fields[1] = "\d+ warnings"
                $numWarnings = extractLeadingNumber $fields[1];
                $entry->BuildWarnings($k,$numWarnings);
                # get numErrors from $fields[2] = "\d+ errors"
                $numErrors = extractLeadingNumber $fields[2];
                $entry->BuildErrors($k,$numErrors);
                # get pass/total
                my ($pass,$totalString) = split /\//, $fields[3];
                # print "getting total from totalString: $totalString\n";
                my $total = extractLeadingNumber $totalString;
                my $fail = $total - $pass;
                $entry->TestPasses($k,$pass);
                $entry->TestFailures($k,$fail);
                print $packages[$k] . ": test pass/fail :" . $pass . "/" . $fail . "\n";
                $have_match=1;
            }  
            if($have_match == 1) {
                last;
            }
        }#end for num_packages

        if($have_match == 0) {
            print "UNMATCHED whichTests: $whichTests.\n";
        }
    }# end for numEntryLines

    push @entries, $entry;
    $name = $entry->name;
    if ($testType{$name} eq "Nightly")
    {
      push @entriesNightly, $entry;
    }
    elsif ($testType{$name} eq "Experimental")
    {
      push @entriesExperimental, $entry;
    }
    elsif ($testType{$name} eq "Specialized")
    {
      push @entriesSpecialized, $entry;
    }
  }
}

# print "There are " . (scalar @entries) . " entries.\n";

my $theDate = $ARGV[0];
my $cdashDate = $ARGV[1];
my $cdashMachine = $ARGV[2];
my $capMachine = ucfirst($cdashMachine);
my $senderName = $ARGV[3];
$theDate =~ s/_/ /g;

print <<EOF;
From: $senderName\@sandia.gov
Subject: $cdashMachine test summary, $theDate
Mime-Version: 1.0
Content-Type: text/html

<!DOCTYPE html>
<html lang="en-US">

  <head>
    
    <title>Dashboard</title>
    <style>
      .thintable {
        border: 1px solid black;
        border-collapse: collapse;
      }

      .grptr {
        border-bottom: 4px solid black
      }
      .midtr {
        border-bottom: 1px solid black
      }

      .grpth {
        border-left: 4px solid black;
      }
      .midth {
        border-style: none;
        border-left: 1px solid black;
        text-align: center;
      }

      .grptdw {
        border-left: 4px solid black;
        text-align: center;
      }
      .midtdw {
        border-style: none;
        border-left: 1px solid black;
        text-align: center;
      }
      .grptdg {
        border-left: 4px solid black;
        text-align: center;
        background-color: lightgreen;
        color: black;
      }
      .midtdg {
        border-style: none;
        border-left: 1px solid black;
        text-align: center;
        background-color: lightgreen;
        color: black;
      }
      .grptdr {
        border-left: 4px solid black;
        text-align: center;
        background-color: tomato;
        color: black;
      }
      .midtdr {
        border-style: none;
        border-left: 1px solid black;
        text-align: center;
        background-color: tomato;
        color: black;
      }
      .grptdy {
        border-left: 4px solid black;
        text-align: center;
        background-color: yellow;
        color: black;
      }
      .midtdy {
        border-style: none;
        border-left: 1px solid black;
        text-align: center;
        background-color: yellow;
        color: black;
      }
      .grptdc {
        border-left: 4px solid black;
        text-align: center;
        background-color: cyan;
        color: black;
      }
      .midtdc {
        border-style: none;
        border-left: 1px solid black;
        text-align: center;
        background-color: cyan;
        color: black;
      }
      .grptdh {
        border-left: 4px solid black;
        text-align: center;
        background-color: wheat;
        color: black;
      }
      .midtdh {
        border-style: none;
        border-left: 1px solid black;
        text-align: center;
        background-color: wheat;
        color: black;
      }
      .grptdm {
        border-left: 4px solid black;
        text-align: center;
        background-color: magenta;
        color: black;
      }
      .midtdm {
        border-style: none;
        border-left: 1px solid black;
        text-align: center;
        background-color: magenta;
        color: black;
      }
    </style>
  </head>
  <body>

Go to the <a
href="http://testing.sandia.gov/cdash/index.php?project=Trilinos&filtercount=1&showfilters=1&field1=site&compare1=63&value1=$cdashMachine&date=$cdashDate">full report</a>
<br>


<h2>$capMachine Test Summary</h2>
<h2>$gitFailure</h2>
EOF

sub printTableHeader
{
    print '<table class="thintable">'; print "\n";
    print '<tr style="border-style: none;">'; print "\n";
    print '<th></th>'; print "\n";
        
    for(my $k=0; $k<$num_packages; $k++) {
print <<"EOF";
          <th class="grpth" colspan="4" style="text-align: center;">$packages[$k]</th>
          <th class="grpth" colspan="2" style="text-align: center;">$packages[$k] libs</th>
EOF

    }

print <<"EOF";
        </tr>
        <tr class="grptr">
        <th>Test</th>
EOF

   for(my $k=0; $k<$num_packages; $k++) {
print <<"EOF";
            <th class="grpth">pass</th>
            <th class="midth">fail/notrun</th>
            <th class="midth">errors</th>
            <th class="midth">warnings</th>
            <th class="grpth">errors</th>
            <th class="midth">warnings</th>
EOF
     }

print <<"EOF";
</tr>
EOF

}

sub printTableFooter
{
 print <<EOF;
</table>
EOF
}

sub printEntries {
  my @entriesToPrint = @{$_[0]};
  for $entry (@entriesToPrint)
  {

    $name = $entry->name;
    print <<"EOF";
<tr class="midtr">
<td>$name</td>
EOF
    for(my $k=0; $k<$num_packages; $k++) {
        my $TestPassColor = 'g';
        my $TestFailureColor = 'g';
        if ($entry->TestPasses($k) + $entry->TestFailures($k) == 0) {
            $TestPassColor = 'y';
            $TestFailureColor = 'y';
        }
        elsif ($entry->TestFailures($k) > 0) {
            $TestFailureColor = 'r';
        }

        my $BuildErrorColor      = ($entry->BuildErrors($k)      > 0) ? 'r' : 'g';
        my $BuildWarningColor    = ($entry->BuildWarnings($k)    > 0) ? 'y' : 'g';
        my $LibBuildErrorColor   = ($entry->LibBuildErrors($k)   > 0) ? 'r' : 'g';
        my $LibBuildWarningColor = ($entry->LibBuildWarnings($k) > 0) ? 'y' : 'g';

        #print $packages[$k] . ": " . $TestPassColor . " " . $entry->TestPasses($k) . "/" . $entry->TestFailures($k). "\n";

        my $TestPasses       = $entry->TestPasses($k);
        my $TestFailures     = $entry->TestFailures($k);
        my $BuildErrors      = $entry->BuildErrors($k);
        my $BuildWarnings    = $entry->BuildWarnings($k);
        my $LibBuildErrors   = $entry->LibBuildErrors($k);
        my $LibBuildWarnings = $entry->LibBuildWarnings($k);
print <<"EOF";  
<td class="grptd$TestPassColor"><b>$TestPasses</b></td>
<td class="midtd$TestFailureColor"><b>$TestFailures</b></td>
<td class="midtd$BuildErrorColor"><b>$BuildErrors</b></td>
<td class="midtd$BuildWarningColor"><b>$BuildWarnings</b></td>
<td class="grptd$LibBuildErrorColor"><b>$LibBuildErrors</b></td>
<td class="midtd$LibBuildWarningColor"><b>$LibBuildWarnings</b></td>      
EOF
   }
  }
}


print "<h3>Nightly Tests</h3>\n";
if (@entriesNightly > 0)
{
  printTableHeader;
  printEntries(\@entriesNightly);
  printTableFooter;
}
else
{
  print "no results reported\n";
}

print "<h3>Specialized Tests</h3>\n";
if (@entriesSpecialized > 0)
{
  printTableHeader;
  printEntries(\@entriesSpecialized);
  printTableFooter;
}
else
{
  print "no results reported\n";
}

print "<h3>Experimental Tests</h3>\n";
if (@entriesExperimental > 0)
{
  printTableHeader;
  printEntries(\@entriesExperimental);
  printTableFooter;
}
else
{
  print "no results reported\n";
}

print <<EOF;
<br>
<hr>

<!-- Original email text follows -->
<pre>
EOF

foreach $line (@LINES)
{
  print $line;
}

print <<"EOF";
</pre>
</body>
</html>
EOF
