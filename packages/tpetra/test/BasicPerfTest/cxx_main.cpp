// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Tpetra_ConfigDefs.hpp" // for <iostream> and <stdlib>
#include <Teuchos_Time.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include "Tpetra_Version.cpp"

//function prototypes
void doSumRuns(int size, int numRuns, bool verbose);
void doSum(int const size, double& timeA, double& timeB, double& timeC, double& timeD, double& timeE);

int main(int argc, char* argv[]) {
	// initialize verbose & debug flags
	bool verbose = false;
	bool override = false;
	if(argc > 1)
		if(argv[1][0] == '-') {
			if(argv[1][1] == 'v' || argv[1][2] == 'v')
				verbose = true;
			if(argv[1][1] == 'o' || argv[1][2] == 'o')
				override = true;
		}

 	if (verbose)
		cout << Tpetra::Tpetra_Version() << endl << endl;

	int base = 2;
	int iterations = 20;
	int runsPerSize = 5;

	if(override) {
		cout << "Base? ";
		cin >> base;
		cout << "How many iterations? ";
		cin >> iterations;
		cout << "How many runs per iteration? ";
		cin >> runsPerSize;
	}

	cout.setf(ios::left, ios::adjustfield);
	cout << setw(20) << "Size" << setw(20) << "Pointer Arith." 
			 << setw(20) << "Pointer(optimized)" << setw(20) << "Array offset" 
			 << setw(20) << "STL accumulate" << setw(20) << "STL iterator" 
       << setw(20) << "Row time" << endl;
	
	int size = base;
	for(int i = 0; i < iterations; i++, size *= base)
		doSumRuns(size, runsPerSize, verbose);

	return(0);
}

void doSumRuns(int size, int numRuns, bool verbose) {
	double totalTimeA = 0.0;
	double totalTimeB = 0.0;
	double totalTimeC = 0.0;
	double totalTimeD = 0.0;
  double totalTimeE = 0.0;
	double timeA, timeB, timeC, timeD, timeE;

	//cout.setf(ios::scientific, ios::floatfield);
	//cout.setf(ios::showpoint);
	//cout.setf(ios::left, ios::adjustfield);

	Teuchos::Time wholeloop("allRunsTogether");
	wholeloop.start();

	for(int i = 0; i < numRuns; i++) {
		doSum(size, timeA, timeB, timeC, timeD, timeE);
		totalTimeA += timeA;
		totalTimeB += timeB;
		totalTimeC += timeC;
		totalTimeD += timeD;
    totalTimeE += timeE;

		if(verbose) {
			cout << setw(20) << size 
					 << setw(20) << timeA 
					 << setw(20) << timeB 
					 << setw(20) << timeC
					 << setw(20) << timeD
           << setw(20) << timeE
					 << endl;
		}
	}

	totalTimeA /= static_cast<double> (numRuns);
	totalTimeB /= static_cast<double> (numRuns);
	totalTimeC /= static_cast<double> (numRuns);
	totalTimeD /= static_cast<double> (numRuns);
  totalTimeE /= static_cast<double> (numRuns);

	wholeloop.stop();

	if(verbose) cout << "--------------------------------------------------------------------------------\n";
	cout << setw(20) << size 
			 << setw(20) << totalTimeA 
			 << setw(20) << totalTimeB 
			 << setw(20) << totalTimeC
			 << setw(20) << totalTimeD
       << setw(20) << totalTimeE
			 << setw(20) << wholeloop.totalElapsedTime()
			 << endl;
	if(verbose) cout << "--------------------------------------------------------------------------------\n";

}

void doSum(int const size, double& timeA, double& timeB, double& timeC, double& timeD, double& timeE) {
	// initialize vector of user-specified size, and set elements to random values.
	std::vector<double> vector1(size, 0.0); 
	for(int i = 0; i < size; i++)
		vector1[i] = Teuchos::ScalarTraits<double>::random();

	// create timer objects
	Teuchos::Time timerA("pointer arithmetic");
	Teuchos::Time timerB("pointer arithmetic (optimized)");
	Teuchos::Time timerC("array offset");
	Teuchos::Time timerD("stl accumulate");
  Teuchos::Time timerE("stl iterator");

	// create double variables to store results in
	double sumA = 0.0;
	double sumB = 0.0;
	double sumC = 0.0;
	double sumD = 0.0;
  double sumE = 0.0;

	// Time using pointer arithmetic
	timerA.start();
	double* ptrA = &vector1[0];
	for(int i = 0; i < size; i++) {
		sumA += *ptrA;
		ptrA++;
	}
	timerA.stop();

	// Time using pointer arithmetic (optimized form)
	timerB.start();
	double* ptrB = &vector1[0];
	for(int i = 0; i < size; i++)
		sumB += *ptrB++;
	timerB.stop();

	// Time using array offsets
	timerC.start();
	for(int i = 0; i < size; i++)
		sumC += vector1[i];
	timerC.stop();

	// Time using STL accumulate
	timerD.start();
	sumD = accumulate(vector1.begin(), vector1.end(), 0.0);
	timerD.stop();

  // Time using STL iterator
  timerE.start();
  std::vector<double>::iterator ii = vector1.begin();
  for(int i = 0; i < size; i++)
    sumE += *ii++;
  timerE.stop();
  
  //cout << "sums: " << sumA << ", " << sumB << ", " << sumC << ", " << sumD << ", " << sumE << endl;

	timeA = timerA.totalElapsedTime();
	timeB = timerB.totalElapsedTime();
	timeC = timerC.totalElapsedTime();
	timeD = timerD.totalElapsedTime();
  timeE = timerE.totalElapsedTime();
}
