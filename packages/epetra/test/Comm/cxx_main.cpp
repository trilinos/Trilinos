// Epetra_Comm Test routine
#ifdef EPETRA_MPI
#include <mpi.h>
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialComm.h"

#include "Epetra_Time.h"
#include "Epetra_Util.h"
#include "../epetra_test_err.h"

int main(int argc, char *argv[]) {
  int ierr = 0, forierr = 0;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

#endif

  bool verbose = false;
  bool verbose1 = false;
  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose1 = true;


#ifdef EPETRA_MPI

  Epetra_MpiComm petracomm( MPI_COMM_WORLD );
  int MyPID =  petracomm.MyPID();
  int NumProc =  petracomm.NumProc();
  EPETRA_TEST_ERR(!(petracomm.MyPID()==rank),ierr);
  EPETRA_TEST_ERR(!(petracomm.NumProc()==size),ierr);
  if (verbose1) verbose = (MyPID==0);
  petracomm.Barrier();

  MPI_Comm MPIComm1 = petracomm.Comm();
  int size1, rank1;
  MPI_Comm_size(MPIComm1, &size1);
  MPI_Comm_rank(MPIComm1, &rank1);
  if (verbose1) cout << petracomm <<  ".  Using MPI_Comm from Petra_Comm: "
                    << "Processor "<< rank1 <<" of " << size1
		          << " (should be the same)."<<endl;

  EPETRA_TEST_ERR(!(rank1==rank),ierr);
  EPETRA_TEST_ERR(!(size1==size),ierr);

  // Do some timing to test barrier function
  
  Epetra_Time before_barrier(petracomm);
  Epetra_Time after_barrier(petracomm);
  // Give each processor rank+1 amount of work
  // Time before barrier should increase roughly linearly
  // Time after barrier should be same for all processors
  double sum = 0.0;
  for (int j=0; j<rank+1; j++)
    for (int i=0; i<1000000; i++) sum += drand48();
  sum /= rank+1;
  if (verbose1) cout << "Processor "<<MyPID
		    <<" Time to reach barrier: "
		    << before_barrier.ElapsedTime() << endl;
  petracomm.Barrier();
  if (verbose1) cout << "Processor "<<MyPID << " Sum result  "
		    << sum <<" Time to get beyond barrier: "
		    << after_barrier.ElapsedTime() << endl;
 
  petracomm.Barrier();

// Some vars needed for the following tests
  int count = 4;
  int *iInputs = new int[count]; // General array for int type tests
  for (int i=0; i<count; i++)
    iInputs[i] = 10*(i + rank - 2) + rank; // if these values are changed, the expected maxs, mins, sums, etc must also change.  NOTE: Broadcst() does not use these values.  The lines that need to be changed are located in the "Values for ****** tests" sections directly below.
  double *dInputs = new double[count]; // General array for double type tests
  for (int i=0; i<count; i++)
    dInputs[i] = pow(2.0,i-rank);// if these values are changed, the expected maxs, mins, sums, etc must also change.  NOTE: Broadcst() does not use these values.  The lines that need to be changed are located in the "Values for ****** tests" sections directly below.

  // Values for Broadcast tests
  int *iVals = new int[count];
  if (rank == 0) {
     for (int i=0; i<count; i++)
       iVals[i] = i; // if these values are changed, the values in iBVals must also be changed
  }
  int *iBVals = new int[count]; // Values to be checked against the values broadcast to the non root processors
  for (int i=0; i<count; i++)
    iBVals[i] = i; // if these values are changed, the values in iVals must also be changed
  double *dVals = new double[count];
  if (rank == 0) {
     for (int i=0; i<count; i++)
       dVals[i] = double(i); // if these values are changed, the values in dBVals must also be changed
  }    
  double *dBVals = new double[count];// Values to be checked against the values broadcast to the non root processors
  for (int i=0; i<count; i++)
    dBVals[i] = i; // if these values are changed, the values in dVals must also be changed

  // Values for MaxAll tests
  int *iMyGlobalMaxs = new int[count];
  for (int i=0; i<count; i++)
    iMyGlobalMaxs[i]=10 * (i + NumProc-1 -2) +  NumProc-1; // if these values are changed, iInput must be changed as well as all other values dependent on iInput
  double *dMyGlobalMaxs = new double[count];
  for (int i=0; i<count; i++)
    dMyGlobalMaxs[i]= pow(2.0,i); //if these values are changed, dInput must be changed as well as all other values dependent on dInput

  // Values for MinAll tests
  int *iMyGlobalMins = new int[count];
  for (int i=0; i<count; i++)
    iMyGlobalMins[i]= 10 * (i - 2); // if these values are changed, iInput must be changed as well as all other values dependent on iInput
  double *dMyGlobalMins = new double[count];
  for (int i=0; i<count; i++)
    dMyGlobalMins[i]= pow(2.0,i-(NumProc-1)); //if these values are changed, dInput must be changed as well as all other values dependent on dInput

  // Values for SumAll tests
  int *iMyGlobalSums = new int[count];
  for (int i=0; i<count; i++){
    iMyGlobalSums[i]=0;
    for (int j=0; j<NumProc; j++)
      iMyGlobalSums[i] += 10*(i+j-2) + j;// if these values are changed, iInput must be changed as well as all other values dependent on iInput
  }
  double *dMyGlobalSums = new double[count];
  for (int i=0; i<count; i++){
    dMyGlobalSums[i]=0;
    for (int j=0; j<NumProc; j++)
      dMyGlobalSums[i] += pow(2.0,i-j);// if these values are changed, dInput must be changed as well as all other values dependent on dInput
  }

  // Values for ScanSum tests
  int *iMyScanSums = new int[count];
  for (int i=0; i<count; i++)
    iMyScanSums[i] = int((rank+1)*(10*(2*i+rank-4)+rank)*.5);// if these values are changed, iInput must be changed as well as all other values dependent on iInput
  double *dMyScanSums = new double[count];
  for (int i=0; i<count; i++) {
    dMyScanSums[i] = 0;
    for (int j=0; j<=rank; j++)
      dMyScanSums[i] += pow(2.0,i-j); //if these values are changed, dInput must be changed as well as all other values dependent on dInput
  }
  // Values for Gather tests
  int totalVals = count*NumProc;
  int *iMyOrderedVals = new int[totalVals];
  double *dMyOrderedVals = new double[totalVals];
  int k=0;
  for (int j=0; j<NumProc; j++) {
    for (int i=0; i<count; i++) {
      iMyOrderedVals[k] = 10*(i + j - 2) + j;; // if these values are changed, iInput must be changed as well as all other values dependent on iInput
      dMyOrderedVals[k] = pow(2.0,i-j); // if these values are changed, dInput must be changed as well as all other values dependent on dInput
      k++;
    }
  }
  petracomm.Barrier();

  // Method testing section
  // Test the Broadcast functions
  EPETRA_TEST_ERR(petracomm.Broadcast(iVals,count,0),ierr);
  if (verbose1) {
    if (rank == 0) 
      cout << "The values on the root processor are: ";
    else
      cout << "The values on processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << iVals[i] << " ";
    cout << endl;
  }
  // ierr = 0; need to track errors the whole way through the file - this line of code seems like a bad idea 
  forierr = 0;
  for (int i=0; i<count; i++)
    forierr += !(iVals[i] == iBVals[i]); // otherwise Broadcast didn't occur properly
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] iVals;
  delete [] iBVals;
  petracomm.Barrier();
  if (verbose) cout << endl << "Broadcast (type int) test passed!" << endl << endl;// If test gets to here the test passed, only output on one node
  petracomm.Barrier();

  EPETRA_TEST_ERR(petracomm.Broadcast(dVals,count,0),ierr);
  if (verbose1) {
    if (rank == 0)
      cout << "The values on the root processor are: ";
    else
      cout << "The values on processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << dVals[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (int i=0; i<count; i++)
    forierr += !(dVals[i] == dBVals[i]); // otherwise Broadcast didn't occur properly
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] dVals;
  delete [] dBVals;
  petracomm.Barrier();
  if (verbose) cout << endl << "Broadcast (type double) test passed!" << endl << endl;// If test gets to here the test passed, only output on one node
  petracomm.Barrier();

 // Test the MaxAll functions
  int *iGlobalMaxs = new int[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << iInputs[i] << " ";
    cout << endl;
  }
  EPETRA_TEST_ERR(petracomm.MaxAll(iInputs,iGlobalMaxs,count),ierr);
  petracomm.Barrier();
  
  if (verbose1) {
    cout << "The max values according to processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << iGlobalMaxs[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (int i=0; i<count; i++) 
    forierr += !(iMyGlobalMaxs[i] == iGlobalMaxs[i]);
  EPETRA_TEST_ERR(forierr,ierr);

  delete [] iGlobalMaxs;
  delete [] iMyGlobalMaxs;
  petracomm.Barrier();
  if (verbose) cout << endl << "MaxAll (type int) test passed!" << endl << endl;// If test gets to here the test passed, only output on one node
  petracomm.Barrier();

  double *dGlobalMaxs = new double[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << dInputs[i] << " ";
    cout << endl;
  }
  EPETRA_TEST_ERR(petracomm.MaxAll(dInputs,dGlobalMaxs,count),ierr);
  petracomm.Barrier();
  
  if (verbose1) {
    cout << "The max values according to processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << dGlobalMaxs[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (int i=0; i<count; i++)
    forierr += !(Epetra_Util::Chop(dMyGlobalMaxs[i] - dGlobalMaxs[i]) == 0);
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] dGlobalMaxs;
  delete [] dMyGlobalMaxs;
  petracomm.Barrier();
  if (verbose) cout << endl << "MaxAll (type double) test passed!" << endl << endl;// If test gets to here the test passed, only output on one node
  petracomm.Barrier();

 // Test the MinAll functions
  int *iGlobalMins = new int[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << iInputs[i] << " ";
    cout << endl;
  }
  EPETRA_TEST_ERR(petracomm.MinAll(iInputs,iGlobalMins,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The min values according to processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << iGlobalMins[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (int i=0; i<count; i++) 
    forierr += !(iMyGlobalMins[i] == iGlobalMins[i]); // otherwise calculated min is wrong
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] iGlobalMins;
  delete [] iMyGlobalMins;
  petracomm.Barrier();
  if (verbose) cout << endl << "MinAll (type int) test passed!" << endl << endl;// If test gets to here the test passed, only output on one node
  petracomm.Barrier();

  double *dGlobalMins = new double[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << dInputs[i] << " ";
    cout << endl;
  }
  EPETRA_TEST_ERR(petracomm.MinAll(dInputs,dGlobalMins,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The min values according to processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << dGlobalMins[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (int i=0; i<count; i++)
    forierr += !(Epetra_Util::Chop(dMyGlobalMins[i] - dGlobalMins[i]) == 0); // otherwise calculated min is wrong
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] dGlobalMins;
  delete [] dMyGlobalMins;
  petracomm.Barrier();
  if (verbose) cout << endl << "MinAll (type double) test passed!" << endl << endl;// If test gets to here the test passed, only output on one node
  petracomm.Barrier();

 // Test the SumAll functions
  int *iGlobalSums = new int[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << iInputs[i] << " ";
    cout << endl;
  }
  EPETRA_TEST_ERR(petracomm.SumAll(iInputs,iGlobalSums,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The sums of all values according to processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << iGlobalSums[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (int i=0; i<count; i++)
    forierr += !(iMyGlobalSums[i] == iGlobalSums[i]); // otherwise calculated sum is wrong
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] iGlobalSums;
  delete [] iMyGlobalSums;
  petracomm.Barrier();
  if (verbose) cout << endl << "SumAll (type int) test passed!" << endl << endl;// If test gets to here the test passed, only output on one node
  petracomm.Barrier();

  double *dGlobalSums = new double[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << dInputs[i] << " ";
    cout << endl;
  }
  EPETRA_TEST_ERR(petracomm.SumAll(dInputs,dGlobalSums,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The sums of all values according to processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << dGlobalSums[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (int i=0; i<count; i++)
    forierr += !(Epetra_Util::Chop(dMyGlobalSums[i] - dGlobalSums[i]) == 0); // otherwise calculated sum is wrong
  EPETRA_TEST_ERR(forierr,ierr);

  delete [] dGlobalSums;
  delete [] dMyGlobalSums;
  petracomm.Barrier();
  if (verbose) cout << endl << "SumAll (type double) test passed!" << endl << endl;// If test gets to here the test passed, only output on one node
  petracomm.Barrier();

 // Test the ScanSum functions
  int *iScanSums = new int[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << iInputs[i] << " ";
    cout << endl;
  }
  
  EPETRA_TEST_ERR(petracomm.ScanSum(iInputs,iScanSums,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The sums of all values on processors 0 - " << rank << " are: ";
    for (int i=0; i<count; i++) {
      cout << iScanSums[i] << " ";
    }
    cout << endl;
  }
  forierr = 0;
  for (int i=0; i<count; i++)
    forierr += !(iMyScanSums[i] == iScanSums[i]);
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] iScanSums;
  delete [] iMyScanSums;
  petracomm.Barrier();
  if (verbose) cout << endl << "ScanSum (type int) test passed!" << endl << endl;// If test gets to here the test passed, only output on one node
  petracomm.Barrier();

  double *dScanSums = new double[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (int i=0; i<count; i++)
      cout << dInputs[i] << " ";
    cout << endl;
  }

  EPETRA_TEST_ERR(petracomm.ScanSum(dInputs,dScanSums,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The sums of all values on processors 0 - " << rank << " are: ";
    for (int i=0; i<count; i++) {
      cout << dScanSums[i] << " ";
    }
    cout << endl;
  }
  forierr = 0;
  for (int i=0; i<count; i++)
    forierr += !(Epetra_Util::Chop(dMyScanSums[i] - dScanSums[i])== 0);
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] dScanSums;
  delete [] dMyScanSums;
  petracomm.Barrier();
  if (verbose) cout << endl << "ScanSum (type double) test passed!" << endl << endl;// If test gets to here the test passed, only output on one node
  petracomm.Barrier();

 // Test the Gather functions
  int *iOrderedVals = new int[totalVals];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << iInputs[i] << " ";
    cout << endl;
  }
  
  EPETRA_TEST_ERR(petracomm.GatherAll(iInputs,iOrderedVals,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The combined list of all values from all processors according to processor " << rank << " is: ";
    for (int i=0; i<totalVals; i++) {
      cout << iOrderedVals[i] << " ";
    }
    cout << endl;
  }
  forierr = 0;
  for (int i=0; i<totalVals; i++)
    forierr += !(iMyOrderedVals[i] == iOrderedVals[i]);
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] iOrderedVals;
  delete [] iMyOrderedVals;
  petracomm.Barrier();
  if (verbose) cout << endl << "GatherAll (type int) test passed!" << endl << endl;// If test gets to here the test passed, only output on one node
  petracomm.Barrier();

  double *dOrderedVals = new double[totalVals];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (int i=0; i<count; i++) 
      cout << dInputs[i] << " ";
    cout << endl;
  }
  
  EPETRA_TEST_ERR(petracomm.GatherAll(dInputs,dOrderedVals,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The combined list of all values from all processors according to processor " << rank << " is: ";
    for (int i=0; i<totalVals; i++) {
      cout << dOrderedVals[i] << " ";
    }
    cout << endl;
  }
  forierr = 0;
  for (int i=0; i<totalVals; i++)
    forierr += !(Epetra_Util::Chop(dMyOrderedVals[i] - dOrderedVals[i]) == 0);
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] dOrderedVals;
  delete [] dMyOrderedVals;
  petracomm.Barrier();
  if (verbose) cout << endl << "GatherAll (type double) test passed!" << endl << endl;// If test gets to here the test passed, only output on one node
  petracomm.Barrier();

  delete [] dInputs;
  delete [] iInputs;
#endif

  // Test serial interface first
  Epetra_SerialComm comm;
  int MyPID1 = comm.MyPID();
  int NumProc1 = comm.NumProc();
  if (verbose1) cout << comm << endl;

  EPETRA_TEST_ERR(!(MyPID1==0),ierr);
  EPETRA_TEST_ERR(!(NumProc1==1),ierr);
  comm.Barrier();
  if (verbose1) cout << comm << " is past serial barrier." << endl << flush;

#ifdef EPETRA_MPI
  petracomm.Barrier();
#endif

  if (verbose1) cout << endl << " Epetra_Comm Check OK." << endl;

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
  return ierr;
}

/*
  end of file main.cc
*/
