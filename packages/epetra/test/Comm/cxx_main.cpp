//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <limits.h>

// Epetra_Comm Test routine
#include "../epetra_test_err.h"
#include "Epetra_Time.h"
#include "Epetra_Util.h"
#include "Epetra_Distributor.h"
#include "Epetra_SerialComm.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Version.h"

#ifdef EPETRA_MPI
#include <mpi.h>
#include "Epetra_MpiComm.h"

int checkMpiDataClass(bool verbose);
#endif

int checkSerialDataClass(bool verbose);
int checkCommMethods(Epetra_Comm& petracomm,
                     bool verbose, bool verbose1,
                     int& NumProc, int& rank);
int checkRankAndSize(Epetra_Comm& petracomm, bool verbose, int rank, int size);
void checkBarrier(Epetra_Comm& petracomm, bool verbose, int rank);

int checkDistributor(Epetra_Distributor* distr,
                     Epetra_Comm& Comm);

int main(int argc, char* argv[]) {
  bool verbose = false;  // used to set verbose false on non-root processors
  bool verbose1 = false; // user's command-line argument
  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose1 = true;

  int ierr = 0;
  int returnierr = 0;
  int size = 1;
  int rank = 0;

  if (verbose1)
    cout << Epetra_Version() << endl << endl;

	// Test Epetra_SerialComm
	if(verbose1) cout << "Testing Epetra_SerialComm..." << endl;
	Epetra_SerialComm serialcomm;
  if (verbose1) cout << serialcomm << endl;
	ierr = checkRankAndSize(serialcomm, verbose1, rank, size);
  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose1 && ierr==0) cout << "Checked OK\n\n" <<endl;
	// method testing
	int numProc = serialcomm.NumProc();
	ierr = checkCommMethods(serialcomm, verbose, verbose1, numProc, rank);
  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose1 && ierr==0) cout << "Checked OK\n\n" <<endl;
	// clone
	if(verbose1) cout << "SerialComm Clone.." << endl;
	Epetra_Comm* cloned_serialcomm = serialcomm.Clone();
	ierr = checkCommMethods(*cloned_serialcomm, verbose, verbose1, numProc, rank);
	delete cloned_serialcomm;
  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose1 && ierr==0) cout << "Checked OK\n\n" <<endl;
	// check inner data class
	ierr = checkSerialDataClass(verbose1);
  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose1 && ierr==0) cout << "Checked OK\n\n" <<endl;

  Epetra_Distributor* serialdistr = serialcomm.CreateDistributor();
  ierr = checkDistributor(serialdistr, serialcomm);
  delete serialdistr;
  EPETRA_TEST_ERR(ierr, returnierr);

	// Test Epetra_MpiComm
#ifdef EPETRA_MPI
  // Initialize MPI
	if(verbose1) cout << "Testing Epetra_MpiComm..." << endl;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	Epetra_MpiComm petracomm( MPI_COMM_WORLD );
	ierr = checkRankAndSize(petracomm, verbose1, rank, size);
  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose1 && ierr==0) cout << "Checked OK\n\n" <<endl;
	MPI_Comm MPIComm1 = petracomm.Comm();
	int size1, rank1;
	MPI_Comm_size(MPIComm1, &size1);
	MPI_Comm_rank(MPIComm1, &rank1);
	if (verbose1) cout << petracomm <<  ".  Using MPI_Comm from Petra_Comm:"
                           << " Processor " << rank1 << " of " << size1
                           << " (should be the same)." << endl;
	EPETRA_TEST_ERR(!(rank1==rank),ierr);
	EPETRA_TEST_ERR(!(size1==size),ierr);
	checkBarrier(petracomm, verbose1, rank);

 	// method testing
	numProc = petracomm.NumProc();
	ierr = checkCommMethods(petracomm, verbose, verbose1, numProc, rank);
  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose1 && ierr==0) cout << "Checked OK\n\n" <<endl;

	// clone
	if(verbose1) cout << "MpiComm Clone.." << endl;
	Epetra_Comm* cloned_mpicomm = petracomm.Clone();
	ierr = checkCommMethods(*cloned_mpicomm, verbose, verbose1, numProc, rank);
	delete cloned_mpicomm;
  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose1 && ierr==0) cout << "Checked OK\n\n" <<endl;

	// check inner data class
	petracomm.Barrier();
	ierr = checkMpiDataClass(verbose1);
  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose1 && ierr==0) cout << "Checked OK\n\n" <<endl;

  Epetra_Distributor* plldistr = petracomm.CreateDistributor();
  ierr = checkDistributor(plldistr, petracomm);
  delete plldistr;
  EPETRA_TEST_ERR(ierr, returnierr);

  petracomm.Barrier();
  MPI_Finalize();
#endif

  return(returnierr);
}

//=============================================================================
void checkBarrier(Epetra_Comm& petracomm, bool verbose, int rank) {
  // Do some timing to test barrier function
  int MyPID = petracomm.MyPID();
  Epetra_Time before_barrier(petracomm);
  Epetra_Time after_barrier(petracomm);
  // Give each processor rank+1 amount of work
  // Time before barrier should increase roughly linearly
  // Time after barrier should be same for all processors
  double sum = 0.0;
  for (int j=0; j<rank+1; j++)
    for (int i=0; i<1000000; i++) 
			sum += ((double )rand())/((double) RAND_MAX);
  sum /= rank+1;
  if (verbose) cout << "Processor " << MyPID
										<< " Time to reach barrier: "
										<< before_barrier.ElapsedTime() << endl;
  petracomm.Barrier();
  if (verbose) cout << "Processor " << MyPID << " Sum result  "
										<< sum << " Time to get beyond barrier: "
										<< after_barrier.ElapsedTime() << endl;
 
  petracomm.Barrier();
}

//=============================================================================
int checkRankAndSize(Epetra_Comm& petracomm, bool verbose, int rank, int size) {
	int ierr = 0;
	//if(verbose) cout << "CRS Breakpoint 1" << endl;
  int MyPID = petracomm.MyPID();
	//if(verbose) cout << "CRS Breakpoint 2" << endl;
  int NumProc = petracomm.NumProc();
  EPETRA_TEST_ERR(!(MyPID==rank),ierr);
  EPETRA_TEST_ERR(!(NumProc==size),ierr);
  petracomm.Barrier();
	return(ierr);
}

//=============================================================================
int checkCommMethods(Epetra_Comm& petracomm, bool verbose, bool verbose1, int& NumProc, int& rank) {
	int i,j;
	int forierr = 0;
	int ierr = 0;

  verbose = (petracomm.MyPID() == 0);  // Turn verbose on; 
                                       // it is always false in main.

	// Some vars needed for the following tests
  int count = 4;
  int* iInputs = new int[count]; // General array for int type tests
  for (i=0; i<count; i++)
    iInputs[i] = 10*(i + rank - 2) + rank; 
	  // if these values are changed, the expected maxs, mins, sums, etc must also change.  
	  //NOTE: Broadcst() does not use these values.  The lines that need to be changed are located
	  //in the "Values for ****** tests" sections directly below.

  double* dInputs = new double[count]; // General array for double type tests
  for (i=0; i<count; i++)
    dInputs[i] = pow(2.0,i-rank); 
	  // if these values are changed, the expected maxs, mins, sums, etc must also change.  
	  //NOTE: Broadcst() does not use these values.  The lines that need to be changed are located
	  //in the "Values for ****** tests" sections directly below.


  // Values for Broadcast tests
  int* iVals = new int[count];
  if (rank == 0) {
		for (i=0; i<count; i++)
       iVals[i] = i; // if these values are changed, the values in iBVals must also be changed
  }
  
  int* iBVals = new int[count]; // Values to be checked against the values broadcast to the non root processors
  for (i=0; i<count; i++)
    iBVals[i] = i; // if these values are changed, the values in iVals must also be changed
  double* dVals = new double[count];
  if (rank == 0) {
     for (i=0; i<count; i++)
       dVals[i] = double(i); // if these values are changed, the values in dBVals must also be changed
  }    
  double* dBVals = new double[count]; // Values to be checked against the values broadcast to the non root processors
  for (i=0; i<count; i++)
    dBVals[i] = i; // if these values are changed, the values in dVals must also be changed

  long long* llVals = new long long[count];
  if (rank == 0) {
     for (i=0; i<count; i++)
       llVals[i] = i+INT_MAX; // if these values are changed, the values in llBVals must also be changed
  }    
  long long* llBVals = new long long[count]; // Values to be checked against the values broadcast to the non root processors
  for (i=0; i<count; i++)
    llBVals[i] = i+INT_MAX; // if these values are changed, the values in dVals must also be changed

  const char *cConst = "Heidi, do you want a cookie?";
  int cCount = strlen(cConst)+1;
  char* cVals = new char[cCount];
  if (rank == 0) {
     strcpy(cVals, cConst);        // if these values are changed, 
     cVals[cCount-1] = '\0';       // the values in cBVals must also be changed
  }
  char* cBVals = new char[cCount]; // Values to be checked against the values 
                                   // broadcast to the non root processors
  strcpy(cBVals, cConst);          // if these values are changed, 
  cBVals[cCount-1] = '\0';         // the values in cVals must also be changed

  // Values for MaxAll tests
  int* iMyGlobalMaxs = new int[count];
  for (i=0; i<count; i++)
    iMyGlobalMaxs[i]=10 * (i + NumProc-1 -2) +  NumProc-1; // if these values are changed, iInput must be changed 
	                                                         //as well as all other values dependent on iInput
  double* dMyGlobalMaxs = new double[count];
  for (i=0; i<count; i++)
    dMyGlobalMaxs[i]= pow(2.0,i); //if these values are changed, dInput must be changed 
	                                //as well as all other values dependent on dInput


  // Values for MinAll tests
  int* iMyGlobalMins = new int[count];
  for (i=0; i<count; i++)
    iMyGlobalMins[i]= 10 * (i - 2); //if these values are changed, iInput must be changed 
	                                  //as well as all other values dependent on iInput
  double* dMyGlobalMins = new double[count];
  for (i=0; i<count; i++)
    dMyGlobalMins[i]= pow(2.0,i-(NumProc-1)); //if these values are changed, dInput must be changed 
	                                            //as well as all other values dependent on dInput


  // Values for SumAll tests
  int* iMyGlobalSums = new int[count];
  for (i=0; i<count; i++){
    iMyGlobalSums[i]=0;
    for (j=0; j<NumProc; j++)
      iMyGlobalSums[i] += 10*(i+j-2) + j;// if these values are changed, iInput must be changed 		                                     
  }                                      //as well as all other values dependent on iInput

  double* dMyGlobalSums = new double[count];
  for (i=0; i<count; i++){
    dMyGlobalSums[i]=0;
    for (j=0; j<NumProc; j++)
      dMyGlobalSums[i] += pow(2.0,i-j);// if these values are changed, dInput must be changed 
	}                                    //as well as all other values dependent on dInput


  // Values for ScanSum tests
  int* iMyScanSums = new int[count];
  for (i=0; i<count; i++)
    iMyScanSums[i] = int((rank+1)*(10*(2*i+rank-4)+rank)*.5);// if these values are changed, 
	                                                           //iInput must be changed as well as 
	                                                           //all other values dependent on iInput
  double* dMyScanSums = new double[count];
  for (i=0; i<count; i++) {
    dMyScanSums[i] = 0;
    for (j=0; j<=rank; j++)
      dMyScanSums[i] += pow(2.0,i-j); //if these values are changed, dInput must be changed 
	}	                                  //as well as all other values dependent on dInput


  // Values for Gather tests
  int totalVals = count*NumProc;
  int* iMyOrderedVals = new int[totalVals];
  double* dMyOrderedVals = new double[totalVals];
  int k=0;
  for (j=0; j<NumProc; j++) {
    for (i=0; i<count; i++) {
      iMyOrderedVals[k] = 10*(i + j - 2) + j;; // if these values are changed, iInput must be changed 
			                                         //as well as all other values dependent on iInput
      dMyOrderedVals[k] = pow(2.0,i-j); // if these values are changed, dInput must be changed 
			                                  //as well as all other values dependent on dInput
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
    for (i=0; i<count; i++) 
      cout << iVals[i] << " ";
    cout << endl;
  }
  // ierr = 0; need to track errors the whole way through the file - this line of code seems like a bad idea 
  forierr = 0;
  for (i=0; i<count; i++)
    forierr += !(iVals[i] == iBVals[i]); // otherwise Broadcast didn't occur properly
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] iVals;
  delete [] iBVals;
  petracomm.Barrier();
  if (verbose) cout << endl << "Broadcast (type int) test passed!" << endl << endl;// If test gets to here the test passed, 
	                                                                                 //only output on one node
  petracomm.Barrier();

  EPETRA_TEST_ERR(petracomm.Broadcast(dVals,count,0),ierr);
  if (verbose1) {
    if (rank == 0)
      cout << "The values on the root processor are: ";
    else
      cout << "The values on processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << dVals[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (i=0; i<count; i++)
    forierr += !(dVals[i] == dBVals[i]); // otherwise Broadcast didn't occur properly
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] dVals;
  delete [] dBVals;
  petracomm.Barrier();
  if (verbose) cout << endl << "Broadcast (type double) test passed!" << endl << endl;// If test gets to here the test passed, 
	                                                                                    //only output on one node
  petracomm.Barrier();

  EPETRA_TEST_ERR(petracomm.Broadcast(llVals,count,0),ierr);
  if (verbose1) {
    if (rank == 0)
      cout << "The values on the root processor are: ";
    else
      cout << "The values on processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << llVals[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (i=0; i<count; i++)
    forierr += !(llVals[i] == llBVals[i]); // otherwise Broadcast didn't occur properly
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] llVals;
  delete [] llBVals;
  petracomm.Barrier();
  if (verbose) cout << endl << "Broadcast (type long long) test passed!" << endl << endl;// If test gets to here the test passed, 
	                                                                                    //only output on one node
  petracomm.Barrier();

  EPETRA_TEST_ERR(petracomm.Broadcast(cVals,cCount,0),ierr);
  if (verbose1) {
    if (rank == 0)
      cout << "The values on the root processor are: " << cVals << endl;
    else
      cout << "The values on processor " << rank << " are: " << cVals << endl;
  }
  forierr = 0;
  for (i=0; i<cCount; i++)
    forierr += !(cVals[i] == cBVals[i]); // otherwise Broadcast didn't work.
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] cVals;
  delete [] cBVals;
  petracomm.Barrier();
  // If test gets to here the test passed, 
  if (verbose) 
    cout << endl << "Broadcast (type char) test passed!" << endl << endl;
	                                                                                    //only output on one node
  petracomm.Barrier();

 // Test the MaxAll functions
  int* iGlobalMaxs = new int[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << iInputs[i] << " ";
    cout << endl;
  }
  EPETRA_TEST_ERR(petracomm.MaxAll(iInputs,iGlobalMaxs,count),ierr);
  petracomm.Barrier();
  
  if (verbose1) {
    cout << "The max values according to processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << iGlobalMaxs[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (i=0; i<count; i++) 
    forierr += !(iMyGlobalMaxs[i] == iGlobalMaxs[i]);
  EPETRA_TEST_ERR(forierr,ierr);

  delete [] iGlobalMaxs;
  delete [] iMyGlobalMaxs;
  petracomm.Barrier();
  if (verbose) cout << endl << "MaxAll (type int) test passed!" << endl << endl;// If test gets to here the test passed, 
	                                                                              //only output on one node
  petracomm.Barrier();

  double* dGlobalMaxs = new double[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << dInputs[i] << " ";
    cout << endl;
  }
  EPETRA_TEST_ERR(petracomm.MaxAll(dInputs,dGlobalMaxs,count),ierr);
  petracomm.Barrier();
  
  if (verbose1) {
    cout << "The max values according to processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << dGlobalMaxs[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (i=0; i<count; i++)
    forierr += !(Epetra_Util::Chop(dMyGlobalMaxs[i] - dGlobalMaxs[i]) == 0);
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] dGlobalMaxs;
  delete [] dMyGlobalMaxs;
  petracomm.Barrier();
  if (verbose) cout << endl << "MaxAll (type double) test passed!" << endl << endl;// If test gets to here the test passed, 
	                                                                                 //only output on one node
  petracomm.Barrier();


 // Test the MinAll functions
  int* iGlobalMins = new int[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << iInputs[i] << " ";
    cout << endl;
  }
  EPETRA_TEST_ERR(petracomm.MinAll(iInputs,iGlobalMins,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The min values according to processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << iGlobalMins[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (i=0; i<count; i++) 
    forierr += !(iMyGlobalMins[i] == iGlobalMins[i]); // otherwise calculated min is wrong
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] iGlobalMins;
  delete [] iMyGlobalMins;
  petracomm.Barrier();
  if (verbose) cout << endl << "MinAll (type int) test passed!" << endl << endl;// If test gets to here the test passed, 
	                                                                              //only output on one node
  petracomm.Barrier();

  double* dGlobalMins = new double[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << dInputs[i] << " ";
    cout << endl;
  }
  EPETRA_TEST_ERR(petracomm.MinAll(dInputs,dGlobalMins,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The min values according to processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << dGlobalMins[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (i=0; i<count; i++)
    forierr += !(Epetra_Util::Chop(dMyGlobalMins[i] - dGlobalMins[i]) == 0); // otherwise calculated min is wrong
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] dGlobalMins;
  delete [] dMyGlobalMins;
  petracomm.Barrier();
  if (verbose) cout << endl << "MinAll (type double) test passed!" << endl << endl;// If test gets to here the test passed, 
	                                                                                 //only output on one node
  petracomm.Barrier();


 // Test the SumAll functions
  int* iGlobalSums = new int[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << iInputs[i] << " ";
    cout << endl;
  }
  EPETRA_TEST_ERR(petracomm.SumAll(iInputs,iGlobalSums,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The sums of all values according to processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << iGlobalSums[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (i=0; i<count; i++)
    forierr += !(iMyGlobalSums[i] == iGlobalSums[i]); // otherwise calculated sum is wrong
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] iGlobalSums;
  delete [] iMyGlobalSums;
  petracomm.Barrier();
  if (verbose) cout << endl << "SumAll (type int) test passed!" << endl << endl;// If test gets to here the test passed, 
	                                                                              //only output on one node
  petracomm.Barrier();

  double* dGlobalSums = new double[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << dInputs[i] << " ";
    cout << endl;
  }
  EPETRA_TEST_ERR(petracomm.SumAll(dInputs,dGlobalSums,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The sums of all values according to processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << dGlobalSums[i] << " ";
    cout << endl;
  }
  forierr = 0;
  for (i=0; i<count; i++)
    forierr += !(Epetra_Util::Chop(dMyGlobalSums[i] - dGlobalSums[i]) == 0); // otherwise calculated sum is wrong
  EPETRA_TEST_ERR(forierr,ierr);

  delete [] dGlobalSums;
  delete [] dMyGlobalSums;
  petracomm.Barrier();
  if (verbose) cout << endl << "SumAll (type double) test passed!" << endl << endl;// If test gets to here the test passed, 
	                                                                                 //only output on one node
  petracomm.Barrier();


 // Test the ScanSum functions
  int* iScanSums = new int[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << iInputs[i] << " ";
    cout << endl;
  }
  
  EPETRA_TEST_ERR(petracomm.ScanSum(iInputs,iScanSums,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The sums of all values on processors 0 - " << rank << " are: ";
    for (i=0; i<count; i++) {
      cout << iScanSums[i] << " ";
    }
    cout << endl;
  }
  forierr = 0;
  for (i=0; i<count; i++)
    forierr += !(iMyScanSums[i] == iScanSums[i]);
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] iScanSums;
  delete [] iMyScanSums;
  petracomm.Barrier();
  if (verbose) cout << endl << "ScanSum (type int) test passed!" << endl << endl;// If test gets to here the test passed, 
	                                                                               //only output on one node
  petracomm.Barrier();

  double* dScanSums = new double[count];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (i=0; i<count; i++)
      cout << dInputs[i] << " ";
    cout << endl;
  }

  EPETRA_TEST_ERR(petracomm.ScanSum(dInputs,dScanSums,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The sums of all values on processors 0 - " << rank << " are: ";
    for (i=0; i<count; i++) {
      cout << dScanSums[i] << " ";
    }
    cout << endl;
  }
  forierr = 0;
  for (i=0; i<count; i++)
    forierr += !(Epetra_Util::Chop(dMyScanSums[i] - dScanSums[i])== 0);
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] dScanSums;
  delete [] dMyScanSums;
  petracomm.Barrier();
  if (verbose) cout << endl << "ScanSum (type double) test passed!" << endl << endl;// If test gets to here the test passed, 
	                                                                                  //only output on one node
  petracomm.Barrier();


 // Test the Gather functions
  int* iOrderedVals = new int[totalVals];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << iInputs[i] << " ";
    cout << endl;
  }
  
  EPETRA_TEST_ERR(petracomm.GatherAll(iInputs,iOrderedVals,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The combined list of all values from all processors according to processor " << rank << " is: ";
    for (i=0; i<totalVals; i++) {
      cout << iOrderedVals[i] << " ";
    }
    cout << endl;
  }
  forierr = 0;
  for (i=0; i<totalVals; i++)
    forierr += !(iMyOrderedVals[i] == iOrderedVals[i]);
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] iOrderedVals;
  delete [] iMyOrderedVals;
  petracomm.Barrier();
  if (verbose) cout << endl << "GatherAll (type int) test passed!" << endl << endl;// If test gets to here the test passed, 
	                                                                                 //only output on one node
  petracomm.Barrier();

  double* dOrderedVals = new double[totalVals];
  if (verbose1) {
    cout << "The values on processor " << rank << " are: ";
    for (i=0; i<count; i++) 
      cout << dInputs[i] << " ";
    cout << endl;
  }
  
  EPETRA_TEST_ERR(petracomm.GatherAll(dInputs,dOrderedVals,count),ierr);
  petracomm.Barrier();

  if (verbose1) {
    cout << "The combined list of all values from all processors according to processor " << rank << " is: ";
    for (i=0; i<totalVals; i++) {
      cout << dOrderedVals[i] << " ";
    }
    cout << endl;
  }
  forierr = 0;
  for (i=0; i<totalVals; i++)
    forierr += !(Epetra_Util::Chop(dMyOrderedVals[i] - dOrderedVals[i]) == 0);
  EPETRA_TEST_ERR(forierr,ierr);
  delete [] dOrderedVals;
  delete [] dMyOrderedVals;
  petracomm.Barrier();
  if (verbose) cout << endl << "GatherAll (type double) test passed!" << endl << endl;// If test gets to here the test passed, 
	                                                                                    //only output on one node
  petracomm.Barrier();

  delete[] dInputs;
  delete[] iInputs;

	return(ierr);
}

//=============================================================================
int checkSerialDataClass(bool verbose) {
	int ierr = 0;
	if(verbose) cout << "Testing Reference Counting... ";								
	Epetra_SerialComm c1;
	int c1count = c1.ReferenceCount();
	const Epetra_SerialCommData* c1addr = c1.DataPtr();
	EPETRA_TEST_ERR(!(c1count==1),ierr); // count should be 1
	if(verbose) cout << "Default constructor. \nc1= " << c1count << "  " << c1addr << endl;

	Epetra_SerialComm* c2 = new Epetra_SerialComm(c1);
	int c2count = c2->ReferenceCount();
	const Epetra_SerialCommData* c2addr = c2->DataPtr();
	int c1countold = c1count;
	c1count = c1.ReferenceCount();
	EPETRA_TEST_ERR(!(c2count==c1count && c1count==(c1countold+1)),ierr); // both counts should be 2
	EPETRA_TEST_ERR(!(c1addr==c2addr),ierr); // addresses should be same
	if(verbose) cout << "Copy constructor(heap). \nc1= " << c1count << "  " << c1addr 
										<< "\nc2= " << c2count << "  " << c2addr << endl;
	delete c2;
	c1countold = c1count;
	c1count = c1.ReferenceCount();
	EPETRA_TEST_ERR(!(c1count==c1countold-1),ierr); // count should have decremented (to 1)
	EPETRA_TEST_ERR(!(c1addr==c1.DataPtr()),ierr); // c1addr should be unchanged
	if(verbose) cout << "c2 Deleted. \nc1= " << c1count << "  " << c1addr << endl;
	{ // inside own set of brackets so that c2a will be automatically at end of brackets
		// so that we can test to make sure objects on the stack deallocate correctly
		Epetra_SerialComm c2a(c1);
		c2count = c2a.ReferenceCount();
		c2addr = c2a.DataPtr();
		c1countold = c1count;
		c1count = c1.ReferenceCount();
		EPETRA_TEST_ERR(!(c2count==c1count && c1count==c1countold+1),ierr); // both counts should be 2
		EPETRA_TEST_ERR(!(c1addr==c2addr),ierr); // addresses should be same
		if(verbose) cout << "Copy constructor(stack). \nc1= " << c1count << "  " << c1addr 
											<< "\nc2a= " << c2count << "  " << c2addr << endl;
	}
	c1countold = c1count;
	c1count = c1.ReferenceCount();
	EPETRA_TEST_ERR(!(c1count==c1countold-1),ierr); // count should have decremented (to 1)
	EPETRA_TEST_ERR(!(c1addr==c1.DataPtr()),ierr); // c1addr should be unchanged
	if(verbose) cout << "c2a Destroyed. \nc1= " << c1count << "  " << c1addr << endl;
	if(verbose) cout << "Assignment operator, post construction" << endl;
	Epetra_SerialComm c3;
	int c3count = c3.ReferenceCount();
	const Epetra_SerialCommData* c3addr = c3.DataPtr();
	EPETRA_TEST_ERR(!(c3count==1),ierr); // c3count should be 1 initially
	EPETRA_TEST_ERR(!(c1addr!=c3addr),ierr); // c1 and c3 should have different ptr addresses
	if(verbose)cout << "Prior to assignment: \nc1=" << c1count << "  " << c1addr 
									 << "\nc3=" << c3count << "  " << c3addr << endl;
	c3 = c1;
	c3count = c3.ReferenceCount();
	c3addr = c3.DataPtr();
	c1countold = c1count;
	c1count = c1.ReferenceCount();
	EPETRA_TEST_ERR(!(c3count==c1count && c1count==c1countold+1),ierr); // both counts should be 2
	EPETRA_TEST_ERR(!(c1addr==c3addr),ierr); // addresses should be same
	if(verbose)cout << "After assignment: \nc1=" << c1count << "  " << c1addr 
									 << "\nc3=" << c3count << "  " << c3addr << endl;
	return(ierr);
}

//=============================================================================
#ifdef EPETRA_MPI
int checkMpiDataClass(bool verbose) {
	int ierr = 0;
	if(verbose) cout << "Testing Reference Counting... ";								
	Epetra_MpiComm c1( MPI_COMM_WORLD );
	int c1count = c1.ReferenceCount();
	const Epetra_MpiCommData* c1addr = c1.DataPtr();
	EPETRA_TEST_ERR(!(c1count==1),ierr); // count should be 1
	if(verbose) cout << "Default constructor. \nc1= " << c1count << "  " << c1addr << endl;

	Epetra_MpiComm* c2 = new Epetra_MpiComm(c1);
	int c2count = c2->ReferenceCount();
	const Epetra_MpiCommData* c2addr = c2->DataPtr();
	int c1countold = c1count;
	c1count = c1.ReferenceCount();
	EPETRA_TEST_ERR(!(c2count==c1count && c1count==(c1countold+1)),ierr); // both counts should be 2
	EPETRA_TEST_ERR(!(c1addr==c2addr),ierr); // addresses should be same
	if(verbose) cout << "Copy constructor(heap). \nc1= " << c1count << "  " << c1addr 
										<< "\nc2= " << c2count << "  " << c2addr << endl;
	delete c2;
	c1countold = c1count;
	c1count = c1.ReferenceCount();
	EPETRA_TEST_ERR(!(c1count==c1countold-1),ierr); // count should have decremented (to 1)
	EPETRA_TEST_ERR(!(c1addr==c1.DataPtr()),ierr); // c1addr should be unchanged
	if(verbose) cout << "c2 Deleted. \nc1= " << c1count << "  " << c1addr << endl;
	{ // inside own set of brackets so that c2a will be automatically at end of brackets
		// so that we can test to make sure objects on the stack deallocate correctly
		Epetra_MpiComm c2a(c1);
		c2count = c2a.ReferenceCount();
		c2addr = c2a.DataPtr();
		c1countold = c1count;
		c1count = c1.ReferenceCount();
		EPETRA_TEST_ERR(!(c2count==c1count && c1count==c1countold+1),ierr); // both counts should be 2
		EPETRA_TEST_ERR(!(c1addr==c2addr),ierr); // addresses should be same
		if(verbose) cout << "Copy constructor(stack). \nc1= " << c1count << "  " << c1addr 
											<< "\nc2a= " << c2count << "  " << c2addr << endl;
	}
	c1countold = c1count;
	c1count = c1.ReferenceCount();
	EPETRA_TEST_ERR(!(c1count==c1countold-1),ierr); // count should have decremented (to 1)
	EPETRA_TEST_ERR(!(c1addr==c1.DataPtr()),ierr); // c1addr should be unchanged
	if(verbose) cout << "c2a Destroyed. \nc1= " << c1count << "  " << c1addr << endl;
	if(verbose) cout << "Assignment operator, post construction" << endl;
	Epetra_MpiComm c3( MPI_COMM_WORLD );
	int c3count = c3.ReferenceCount();
	const Epetra_MpiCommData* c3addr = c3.DataPtr();
	EPETRA_TEST_ERR(!(c3count==1),ierr); // c3count should be 1 initially
	EPETRA_TEST_ERR(!(c1addr!=c3addr),ierr); // c1 and c3 should have different ptr addresses
	if(verbose)cout << "Prior to assignment: \nc1=" << c1count << "  " << c1addr 
									 << "\nc3=" << c3count << "  " << c3addr << endl;
	c3 = c1;
	c3count = c3.ReferenceCount();
	c3addr = c3.DataPtr();
	c1countold = c1count;
	c1count = c1.ReferenceCount();
	EPETRA_TEST_ERR(!(c3count==c1count && c1count==c1countold+1),ierr); // both counts should be 2
	EPETRA_TEST_ERR(!(c1addr==c3addr),ierr); // addresses should be same
	if(verbose)cout << "After assignment: \nc1=" << c1count << "  " << c1addr 
									 << "\nc3=" << c3count << "  " << c3addr << endl;
	return(ierr);
}
#endif

int checkDistributor(Epetra_Distributor* distr,
                     Epetra_Comm& Comm)
{
  int numprocs = Comm.NumProc();

  int numExportIDs = numprocs;
  int* exportPIDs = new int[numExportIDs];
  for(int p=0; p<numExportIDs; ++p) {
    exportPIDs[p] = p;
  }

  bool deterministic = true;
  int numRemoteIDs = 0;

  int err = distr->CreateFromSends(numExportIDs, exportPIDs,
                                   deterministic, numRemoteIDs);

  //numRemoteIDs should equal numExportIDs.

  int returnValue = numRemoteIDs == numExportIDs ? 0 : -99;

  delete [] exportPIDs;

  if (returnValue + err != 0) {
    return(returnValue + err);
  }

  int* exportIDs = new int[numExportIDs];
  for(int i=0; i<numExportIDs; ++i) {
    exportIDs[i] = i+1;
  }

  int len_imports = 0;
  char* imports = NULL;

  err = distr->Do((char*)exportIDs, sizeof(int),
                  len_imports, imports);

  delete [] exportIDs;

  if (len_imports > 0) {
    delete [] imports;
  }

  return(err);
}

/*
  end of file cxx_main.cpp
*/
