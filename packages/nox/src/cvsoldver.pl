#!/usr/bin/perl

if ($#ARGV < 0) {
    die "Usage: cvsoldver.pl <filename> [revision number]\n";
}

$filename = $ARGV[0];

if ($#ARGV < 1) {

    $command = "cvs -q log -h $filename | grep head";
    print "Executing: $command \n";
    $_ = `$command`;
    s/^head:\s(.*)\n/$1/;
    $rev = $_;
    print "Revision Number: ~$rev~\n";
}
else {
    $rev = $ARGV[1];
}

# Create command to generate old file
$command = "cvs update -r $rev -p $filename 1> $filename.~$rev~ 2> /dev/null";
print "Executing: $command \n";
`$command`;



