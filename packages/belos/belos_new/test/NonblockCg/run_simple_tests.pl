#!/bin/perl -w

use strict;

my $program = "../../nbl_pcg_tsfcore_simple.exe";

my $dim = 16;
my $opCondNum = 16;
my $precType=0;
my $combinedOpCondNum=1.0;
my $numRhs = 1;
my $blockSize = 1;
my $maxNumIters = 16;
my $baseTol = 1e-3;
my $varTolMag = 3;
my $outputPrec = 2;
my $useNativeStatusTest = 1;

my $success = 1;

#
# Executable statements
#

mkdir("runs") || die "What: $!";
chdir("runs") || die "What: $!";

$numRhs = 32;

$blockSize = 1;
run_case();

$blockSize = 8;
run_case();

$blockSize = 10;
run_case();

$blockSize = 32;
run_case();

$numRhs = 8;
$blockSize = 8;
run_case();

$dim=1000;
$precType=-1; # Left preconditioning!
$opCondNum = 1e+4;
$combinedOpCondNum=5.0;
run_case();

$precType=+1; # Right preconditioning!
$opCondNum = 1e+4;
$combinedOpCondNum=5.0;
run_case();

$precType=+2; # Left and right preconditioning!
$maxNumIters = 80;
$opCondNum = 1e+4;
$combinedOpCondNum=5.0;
$useNativeStatusTest = 0;
run_case();

$precType=0; # No preconditioning!
$maxNumIters = 500;
$useNativeStatusTest = 1;
run_case();

chdir "..";

exit ( $success ? 0 : 1 );

#
# Subroutines
#

sub run_case {
	my $baseDir = "$precType-$numRhs-$blockSize-$useNativeStatusTest-$dim-$opCondNum-$combinedOpCondNum"
		."-$maxNumIters-$baseTol-$varTolMag";
	#my $baseDir = "dim=$dim.op-cond-num=$opCondNum.prec-type=$precType.combined-op-cond-num=$combinedOpCondNum"
	#	.".num-rhs=$numRhs.block-size=$blockSize.max-num-iters=$maxNumIters"
	#	.".base-tol=$baseTol.var-tol-mag=$varTolMag.use-native-status-test=$useNativeStatusTest";
	mkdir($baseDir) || die "What: $!";
	chdir($baseDir) || die "What: $!";
	my $cmnd;
	$cmnd = "$program";
	$cmnd .= " --dim=$dim --op-cond-num=$opCondNum";
	$cmnd .= " --prec-type=\"none\"" if($precType==(0));
	$cmnd .= " --prec-type=\"left\"" if($precType==(-1));
	$cmnd .= " --prec-type=\"right\"" if($precType==(+1));
	$cmnd .= " --prec-type=\"left-right\"" if($precType==(+2));
	$cmnd .= " --combined-op-cond-num=$combinedOpCondNum";
	$cmnd .= " --num-rhs=$numRhs --block-size=$blockSize --max-num-iters=$maxNumIters"
		." --base-tol=$baseTol --var-tol-mag=$varTolMag --output-prec=$outputPrec";
	$cmnd .= " --use-native-status-test" if($useNativeStatusTest);
	$cmnd .= " --use-non-native-status-test" if(!$useNativeStatusTest);
	$cmnd .= " > console.out";
	print "\n$cmnd\n\n";
	my $result = 0;
	$result = system($cmnd);
	if($result!=0) { $success=0; }
	$cmnd = "grep \'The test\' console.out";
	system($cmnd);
	chdir("..") || die "What: $!";
}
