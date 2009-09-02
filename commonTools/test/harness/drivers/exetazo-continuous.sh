source ~/.keychain/regression.sandia.gov-sh > /dev/null
export CVS_RSH=ssh
export LD_LIBRARY_PATH=/usr/alt/lib:/space/jmwille/contInt/Trilinos/MPI/packages/PyTrilinos/shared

if test -a /home/jmwille/cronOutput/contIntTestsRunning ; then

echo "Tests are already running, contInt tests will not start."
mail  -s "integration - tests running" jmwille@sandia.gov < /home/jmwille/cronOutput/contInt-mpi.txt 
else

#Make sure other tests don't start while these are running
touch /home/jmwille/cronOutput/contIntTestsRunning

cd /space/jmwille/contInt/Trilinos/
cvs -q update -dP | grep -v "? " >& /home/jmwille/cronOutput/contIntCVS.out
diff /home/jmwille/cronOutput/contIntCVS.out /space/jmwille/contInt/Trilinos/commonTools/test/harness/drivers/cvs-update-out.txt > /space/jmwille/contInt/Trilinos/commonTools/test/harness/drivers/cvs-update-out-diff.txt
if test -s /space/jmwille/contInt/Trilinos/commonTools/test/harness/drivers/cvs-update-out-diff.txt ; then

cd /space/jmwille/contInt/Trilinos/commonTools/test/harness

echo "Starting continuous integration build now."
perl runharness --trilinos-dir=/space/jmwille/contInt/Trilinos --build-name=exetazo-mpi-cont >& /home/jmwille/cronOutput/harness-contInt-mpi.txt
#perl runharness --trilinos-dir=/space/jmwille/contInt/Trilinos --build-name=exetazo-mpi-cont --short-circuit >& /home/jmwille/cronOutput/harness-contInt-mpi.txt
more /home/jmwille/cronOutput/harness-contInt-mpi.txt |grep FAILED > /home/jmwille/cronOutput/failedContBuildCheck.txt
if test -s /home/jmwille/cronOutput/failedContBuildCheck.txt ; then
echo "Pausing for failure."
mail  -s "FAILURE in integration tests" trilinos-continuous-integration@software.sandia.gov < /home/jmwille/cronOutput/harness-contInt-mpi.txt
sleep 1800
else
mail  -s "integration tests PASSED" jmwille@sandia.gov < /home/jmwille/cronOutput/harness-contInt-mpi.txt
fi # Failure check
rm -f /home/jmwille/cronOutput/failedContBuildCheck.txt

else

echo "No updates, skipping continuous integration build."
mail  -s "integration - no update" jmwille@sandia.gov < /home/jmwille/cronOutput/contInt-mpi.txt

fi # no update check

rm -f /space/jmwille/contInt/Trilinos/commonTools/test/harness/drivers/cvs-update-out-diff.txt /home/jmwille/cronOutput/contIntCVS.out

# Wait a bit before starting a new test.  Eventually we should only wait if
# there was an error.
#sleep 500
#sleep 20

rm -f /home/jmwille/cronOutput/contIntTestsRunning
 
fi # already running
