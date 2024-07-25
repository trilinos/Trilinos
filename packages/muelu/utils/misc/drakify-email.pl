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

use Class::Struct;
struct Entry => {
  name => '$',
  XpetraBuildWarnings => '$',
  XpetraBuildErrors => '$',
  XpetraTestPasses => '$',  
  XpetraTestFailures => '$',  
  MueLuBuildWarnings => '$',
  MueLuBuildErrors => '$',
  MueLuTestPasses => '$',  
  MueLuTestFailures => '$',
  XpetraLibBuildWarnings => '$',
  XpetraLibBuildErrors => '$',  
  MueLuLibBuildWarnings => '$',
  MueLuLibBuildErrors => '$'
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
$testType = "Nightly";
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
      $testType = "Nightly";
    }
    elsif ($line =~ /\++ Experimental \++/)
    {
      $testType = "Experimental";
    }
    elsif ($line =~ /\++ Specialized \++/)
    {
      $testType = "Specialized";
    }
    
    if ($line =~ /(.*)\s+\.\.\. [pF]/)
    {
      $testName = trim $1;
      $testType{$testName} = $testType;
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
    my $numEntryLines = 4;
    for (my $j=0; $j<$numEntryLines; $j++)
    {
      $i++;
      $line = $LINES[$i];
      my @fields = split /\s*\|\s*/, $line;
#       print "field 0: " . $fields[0] . "\n";
#       print "num fields: " . scalar @fields . "\n";
      my $whichTests = trim $fields[0];
      if ($whichTests =~ /Xpetra_libs/)
      {
        # get numWarnings from $fields[1] = "\d+ warnings"
        $numWarnings = extractLeadingNumber $fields[1];
        $entry->XpetraLibBuildWarnings($numWarnings);
        # get numErrors from $fields[2] = "\d+ errors"
        $numErrors = extractLeadingNumber $fields[2];
        $entry->XpetraLibBuildErrors($numErrors);    
      }
      elsif ($whichTests =~ /MueLu_libs/)
      {
        # get numWarnings from $fields[1] = "\d+ warnings"
        $numWarnings = extractLeadingNumber $fields[1];
        $entry->MueLuLibBuildWarnings($numWarnings);
        # get numErrors from $fields[2] = "\d+ errors"
        $numErrors = extractLeadingNumber $fields[2];
        $entry->MueLuLibBuildErrors($numErrors); 
      }
      elsif ( $whichTests =~ /Xpetra/)
      {
        # get numWarnings from $fields[1] = "\d+ warnings"
#         print "getting numWarnings from fields 1: $fields[1]\n";
        $numWarnings = extractLeadingNumber $fields[1];
        $entry->XpetraBuildWarnings($numWarnings);
        # get numErrors from $fields[2] = "\d+ errors"
        # print "getting numErrors from fields 2: $fields[2]\n";
        $numErrors = extractLeadingNumber $fields[2];
        $entry->XpetraBuildErrors($numErrors);
        # get pass/total
        my ($pass,$totalString) = split /\//, $fields[3];
        # print "getting total from totalString: $totalString\n";
        my $total = extractLeadingNumber $totalString;
        my $fail = $total - $pass;
        $entry->XpetraTestPasses($pass);
        $entry->XpetraTestFailures($fail);
      }
      elsif ($whichTests =~ /MueLu/)
      {
        # get numWarnings from $fields[1] = "\d+ warnings"
        $numWarnings = extractLeadingNumber $fields[1];
        $entry->MueLuBuildWarnings($numWarnings);
        # get numErrors from $fields[2] = "\d+ errors"
        $numErrors = extractLeadingNumber $fields[2];
        $entry->MueLuBuildErrors($numErrors);
        # get pass/total
        my ($pass,$totalString) = split /\//, $fields[3];
        my $total = extractLeadingNumber $totalString;
        my $fail = $total - $pass;
        $entry->MueLuTestPasses($pass);
        $entry->MueLuTestFailures($fail);
      }
      else
      {
        print "UNMATCHED whichTests: $whichTests.\n";
      }
    }
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
href="http://testing.sandia.gov/cdash/index.php?project=Trilinos&subproject=MueLu&filtercount=1&showfilters=1&field1=site&compare1=63&value1=$cdashMachine&date=$cdashDate">full report</a>
<br>


<h2>$capMachine Test Summary</h2>
<h2>$gitFailure</h2>
EOF

sub printTableHeader
{
  print <<EOF;
<table class="thintable">
<tr style="border-style: none;">
<th></th>
<th class="grpth" colspan="4" style="text-align: center;">MueLu</th>
<th class="grpth" colspan="2" style="text-align: center;">MueLu_libs</th>
<th class="grpth" colspan="4" style="text-align: center;">Xpetra</th>
<th class="grpth" colspan="2" style="text-align: center;">Xpetra_libs</th>
<!-- <th class="grpth" colspan="7">Execution Status</th> -->
</tr>
<tr class="grptr">
<th>Test</th>
<th class="grpth">pass</th>
<th class="midth">fail/notrun</th>
<th class="midth">errors</th>
<th class="midth">warnings</th>
<th class="grpth">errors</th>
<th class="midth">warnings</th>
<th class="grpth">pass</th>
<th class="midth">fail</th>
<th class="midth">errors</th>
<th class="midth">warnings</th>
<th class="grpth">errors</th>
<th class="midth">warnings</th>
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
    $XpetraBuildWarnings = $entry->XpetraBuildWarnings;
    $XpetraBuildErrors = $entry->XpetraBuildErrors;
    $XpetraTestPasses = $entry->XpetraTestPasses;
    $XpetraTestFailures = $entry->XpetraTestFailures;
    $XpetraLibBuildWarnings = $entry->XpetraLibBuildWarnings;
    $XpetraLibBuildErrors = $entry->XpetraLibBuildErrors;
    $MueLuBuildWarnings = $entry->MueLuBuildWarnings;
    $MueLuBuildErrors = $entry->MueLuBuildErrors;
    $MueLuTestPasses = $entry->MueLuTestPasses;
    $MueLuTestFailures = $entry->MueLuTestFailures;
    $MueLuLibBuildWarnings = $entry->MueLuLibBuildWarnings;
    $MueLuLibBuildErrors = $entry->MueLuLibBuildErrors;
    my $MueLuTestPassColor = 'g';
    my $MueLuTestFailureColor = 'g';
    if ($MueLuTestPasses + $MueLuTestFailures == 0)
    {
      $MueLuTestPassColor = 'y';
      $MueLuTestFailureColor = 'y';
    }
    elsif ($MueLuTestFailures > 0)
    {
      $MueLuTestFailureColor = 'r';
    }
    my $MueLuBuildErrorColor = ($MueLuBuildErrors > 0) ? 'r' : 'g';
    my $MueLuBuildWarningColor = ($MueLuBuildWarnings > 0) ? 'y' : 'g';
    my $MueLuLibBuildErrorColor = ($MueLuLibBuildErrors > 0) ? 'r' : 'g';
    my $MueLuLibBuildWarningColor = ($MueLuLibBuildWarnings > 0) ? 'y' : 'g';
    my $XpetraTestPassColor = 'g';
    my $XpetraTestFailureColor = 'g';
    if ($XpetraTestPasses + $XpetraTestFailures == 0)
    {
      $XpetraTestPassColor = 'y';
      $XpetraTestFailureColor = 'y';
    }
    elsif ($XpetraTestFailures > 0)
    {
      $XpetraTestFailureColor = 'r';
    }
    my $XpetraBuildErrorColor = ($XpetraBuildErrors > 0) ? 'r' : 'g';
    my $XpetraBuildWarningColor = ($XpetraBuildWarnings > 0) ? 'y' : 'g';
    my $XpetraLibBuildErrorColor = ($XpetraLibBuildErrors > 0) ? 'r' : 'g';
    my $XpetraLibBuildWarningColor = ($XpetraLibBuildWarnings > 0) ? 'y' : 'g';
    print <<EOE;
<tr class="midtr">
<td>$name</td>
<td class="grptd$MueLuTestPassColor"><b>$MueLuTestPasses</b></td>
<td class="midtd$MueLuTestFailureColor"><b>$MueLuTestFailures</b></td>
<td class="midtd$MueLuBuildErrorColor"><b>$MueLuBuildErrors</b></td>
<td class="midtd$MueLuBuildWarningColor"><b>$MueLuBuildWarnings</b></td>
<td class="grptd$MueLuLibBuildErrorColor"><b>$MueLuLibBuildErrors</b></td>
<td class="midtd$MueLuLibBuildWarningColor"><b>$MueLuLibBuildWarnings</b></td>
<td class="grptd$XpetraTestPassColor"><b>$XpetraTestPasses</b></td>
<td class="midtd$XpetraTestFailureColor"><b>$XpetraTestFailures</b></td>
<td class="midtd$XpetraBuildErrorColor"><b>$XpetraBuildErrors</b></td>
<td class="midtd$XpetraBuildWarningColor"><b>$XpetraBuildWarnings</b></td>
<td class="grptd$XpetraLibBuildErrorColor"><b>$XpetraLibBuildErrors</b></td>
<td class="midtd$XpetraLibBuildWarningColor"><b>$XpetraLibBuildWarnings</b></td>
</tr>
EOE
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
