#!/bin/perl -w

use strict;

my $program = "../../nbl_pgmres_tsfcore_simple.exe";

my $dim = 16;
my $opCondNum = 2.0;
my $precType=0;
my $overallOpCondNum=1.0;
my $numRhs = 1;
my $maxNumIters = 16;
my $maxKrylovDim = 16;
my $baseTol = 1e-5;
my $varTolMag = 0.0;
my $outputPrec = 1;
my $useNativeStatusTest = 1;

my $success = 1;

#
# Executable statements
#

mkdir("runs") || die "What: $!";
chdir("runs") || die "What: $!";

$dim = 16;

$numRhs = 1;
run_case();

$maxKrylovDim=4;
$maxNumIters = 32;
run_case();

$useNativeStatusTest = 0;
run_case();

$numRhs = 4;
$maxKrylovDim=$dim;
$maxNumIters = 16;
$useNativeStatusTest = 1;
run_case();

$maxKrylovDim=4;
$maxNumIters = 32;
run_case();

$useNativeStatusTest = 0;
run_case();

$precType=-1; # Left preconditioning!
$numRhs = 1;
$maxKrylovDim=$dim;
$maxNumIters = 16;
$opCondNum=100.0;
$overallOpCondNum=5.0;
$useNativeStatusTest = 0;
run_case();

$precType=-1; # Left preconditioning!
$numRhs = 1;
$maxKrylovDim=$dim;
$maxNumIters = 16;
$opCondNum=100.0;
$overallOpCondNum=5.0;
$useNativeStatusTest = 1;
run_case();

$precType=+1; # Right preconditioning!
$useNativeStatusTest = 1;
run_case();

$precType=+2; # Left and right preconditioning!
$useNativeStatusTest = 0;
run_case();

$precType=0; # No preconditioning
$maxNumIters = 100;
run_case();

chdir "..";

exit ( $success ? 0 : 1 );

#
# Subroutines
#

sub run_case {
	my $baseDir = "dim=$dim.op-cond-num=$opCondNum.prec-type=$precType.overall-op-cond-num=$overallOpCondNum"
		.".num-rhs=$numRhs.max-num-iters=$maxNumIters"
		.".max-krylov-dim=$maxKrylovDim.base-tol=$baseTol.use-native-status-test=$useNativeStatusTest";
	mkdir($baseDir) || die "What: $!";
	chdir($baseDir) || die "What: $!";
	my $cmnd;
	$cmnd = "$program --dim=$dim --op-cond-num=$opCondNum";
	$cmnd .= " --prec-type=\"none\"" if($precType==(0));
	$cmnd .= " --prec-type=\"left\"" if($precType==(-1));
	$cmnd .= " --prec-type=\"right\"" if($precType==(+1));
	$cmnd .= " --prec-type=\"left-right\"" if($precType==(+2));
	$cmnd .= " --num-rhs=$numRhs --max-num-iters=$maxNumIters --max-krylov-dim=$maxKrylovDim"
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
