//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
// ************************************************************************
//@HEADER

// Epetra_IntSerialDense Test routine

#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "../epetra_test_err.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_DataAccess.h"
#include "Epetra_Version.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

// matrix-testing prototypes
int matrixCoverage(bool verbose, bool debug);
int matrixCtr(bool verbose, bool debug);
int matrixCpyCtr(bool verbose, bool debug);
int matrixAssignment(bool verbose, bool debug);
int matrixExceptions(bool verbose, bool debug);
// vector-testing prototypes
int vectorCoverage(bool verbose, bool debug);
int vectorCtr(bool verbose, bool debug);
int vectorCpyCtr(bool verbose, bool debug);
int vectorAssignment(bool verbose, bool debug);
int vectorExceptions(bool verbose, bool debug);
// helper function prototypes
bool identicalSignatures(Epetra_IntSerialDenseMatrix& a, Epetra_IntSerialDenseMatrix& b, bool testLDA = true);
bool seperateData(Epetra_IntSerialDenseMatrix& a, Epetra_IntSerialDenseMatrix& b);
int* getRandArray(int length);
int randomInt();
void printArray(int* array, int length);
void printMat(const char* name, Epetra_IntSerialDenseMatrix& matrix);
void printHeading(const char* heading);

int main(int argc, char *argv[]) {
//============================
// Check for verbose or debug

bool verbose = false;
bool debug = false;
if(argc > 1) {
    if((argv[1][0] == '-') && (argv[1][1] == 'v'))
        verbose = true;
    if((argv[1][0] == '-') && (argv[1][1] == 'd')) {
        debug = true;
        verbose = true;
    }
}

//============================
// Initialize Comm

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  int size = 1; // Serial case (not using MPI)
  int rank = 0;
	Epetra_SerialComm Comm;
#endif

  if (verbose && Comm.MyPID()==0)
    cout << Epetra_Version() << endl << endl;
	
//============================
// other initial setup

  if (!verbose) Comm.SetTracebackMode(0); // This should shut down any error traceback reporting

	int ierr = 0;
	int returnierr = 0;

//============================
// Test vector

	ierr = vectorExceptions(verbose, debug);
	EPETRA_TEST_ERR(ierr, returnierr);

	ierr = vectorCtr(verbose, debug);
	EPETRA_TEST_ERR(ierr, returnierr);
	
	ierr = vectorCpyCtr(verbose, debug);
	EPETRA_TEST_ERR(ierr, returnierr);
	
	ierr = vectorAssignment(verbose, debug);
	EPETRA_TEST_ERR(ierr, returnierr);

	ierr = vectorCoverage(verbose, debug);
	EPETRA_TEST_ERR(ierr, returnierr);

//============================
// Test matrix

	ierr = matrixExceptions(verbose, debug);
	EPETRA_TEST_ERR(ierr, returnierr);

	ierr = matrixCtr(verbose, debug);
	EPETRA_TEST_ERR(ierr, returnierr);

	ierr = matrixCpyCtr(verbose, debug);
	EPETRA_TEST_ERR(ierr, returnierr);

	ierr = matrixAssignment(verbose, debug);
	EPETRA_TEST_ERR(ierr, returnierr);

	ierr = matrixCoverage(verbose, debug);
	EPETRA_TEST_ERR(ierr, returnierr);

//============================
// end of program cleanup

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
	
  return(returnierr);
}

//=========================================================================
//=========================================================================
// Matrix-testing functions
//=========================================================================
//=========================================================================
// call functions we don't have tests for, for sake of code coverage.
// (the functions are called, but their output isn't checked for correctness).
int matrixCoverage(bool verbose, bool debug) {
	if(verbose) printHeading("Testing other matrix functions");

	int* rand1 = getRandArray(9);
	if(debug) printArray(rand1, 9);
	Epetra_IntSerialDenseMatrix m1(Copy, rand1, 3, 3, 3);
	if(debug) printMat("m1",m1);
	delete[] rand1;

	if(verbose) cout << "calling one norm" << endl;
	int onenorm = m1.OneNorm();
	if(debug) cout << "m1.OneNorm() = " << onenorm << endl;

	if(verbose) cout << "calling infinity norm" << endl;
	int infnorm = m1.InfNorm();
	if(debug) cout << "m1.InfNorm() = " << infnorm << endl;

	if(verbose) cout << "calling random" << endl;
	Epetra_IntSerialDenseMatrix m2(6,6);
	if(debug) printMat("m2 (before)",m2);
	m2.Random();
	if(debug) printMat("m2 (after)",m2);

	if(verbose) cout << "Checked OK." << endl;

	return(0);
}

//=========================================================================
// test matrix default constructor, user constructor (copy & view),
// () [] operators (read & write), shape (larger), reshape (smaller)
int matrixCtr(bool verbose, bool debug) {
	const int m1rows = 5;
	const int m1cols = 4;
	const int m1arows = 4;
	const int m1acols = 6;
	const int m2rows = 2;
	const int m2cols = 7;
	const int m3rows = 8;
	const int m3cols = 3;
	const int m3Rrows = 5; // should be smaller than m3rows
	const int m3Rcols = 2; // should be smaller than m3cols

	int ierr = 0;
	int returnierr = 0;
	if(verbose) printHeading("Testing matrix constructors");

	if(verbose) cout << "default constructor" << endl;
	Epetra_IntSerialDenseMatrix m1;
	EPETRA_TEST_ERR(!(m1.CV() == Copy), ierr);
	if(verbose) cout << "shaping" << endl;
	m1.Shape(m1rows, m1cols);
	EPETRA_TEST_ERR(!(m1.M() == m1rows), ierr);
	for(int i = 0; i < m1rows; i++)
		for(int j = 0; j < m1cols; j++)
			EPETRA_TEST_ERR(!(m1(i,j) == 0), ierr);
	if(debug) printMat("m1",m1);
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	if(verbose) cout << "\nmanually setting values" << endl;
	int* m1rand = getRandArray(m1rows * m1cols);
	for(int i = 0; i < m1rows; i++)
		for(int j = 0; j < m1cols; j++)
			m1(i,j) = m1rand[i*m1rows + j];
	for(int i = 0; i < m1rows; i++)
		for(int j = 0; j < m1cols; j++)
		EPETRA_TEST_ERR(!(m1[j][i] == m1rand[i*m1rows + j]), ierr);
	if(debug) {
		printArray(m1rand, m1rows * m1cols);
		printMat("m1",m1);
	}
	delete[] m1rand;
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	
	if(verbose) cout << "\nshaped constructor" << endl;
	Epetra_IntSerialDenseMatrix m1a(m1arows, m1acols);
	EPETRA_TEST_ERR(!(m1a.M() == m1arows), ierr);
	EPETRA_TEST_ERR(!(m1a.N() == m1acols), ierr);
	EPETRA_TEST_ERR(!(m1a.LDA() == m1arows), ierr);
	EPETRA_TEST_ERR(!(m1a.CV() == Copy),ierr);
	if(debug) printMat("m1a", m1a);
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	
	if(verbose) cout << "\nuser data constructor (view)" << endl;
	int* m2rand = getRandArray(m2rows * m2cols);
	if(debug) printArray(m2rand, m2rows * m2cols);
	Epetra_IntSerialDenseMatrix m2(View, m2rand, m2rows, m2rows, m2cols);
	EPETRA_TEST_ERR(!(m2.CV() == View), ierr);
	EPETRA_TEST_ERR(!(m2.M() == m2rows), ierr);
	EPETRA_TEST_ERR(!(m2.N() == m2cols), ierr);
	EPETRA_TEST_ERR(!(m2.LDA() == m2rows), ierr);
	EPETRA_TEST_ERR(!(m2.A() == m2rand), ierr);
	if(debug) printMat("m2",m2);
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nchanging value, checking for correct behavior" << endl;
	int* m2randcopy = new int[m2rows * m2cols]; // backup of original values
	for(int i = 0; i < m2rows * m2cols; i++)
		m2randcopy[i] = m2rand[i];
	m2(0,4) = m2(0,4) + 1; // magic numbers for which element to change
	m2randcopy[4 * m2rows] = m2randcopy[4 * m2rows] + 1;
	for(int i = 0; i < m2rows * m2cols; i++)
		EPETRA_TEST_ERR(!(m2rand[i] == m2randcopy[i]), ierr); // m2rand should have updated correctly
	if(debug) {
		printArray(m2rand, m2rows * m2cols);
		printMat("m2",m2);
	}
	delete[] m2rand;
	delete[] m2randcopy;
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nuser data constructor (copy)" << endl;
	int* m3rand = getRandArray(m3rows * m3cols);
	if(debug) printArray(m3rand, m3rows * m3cols);
	int* m3randcopy = new int[m3rows * m3cols];
	for(int i = 0; i < m3rows * m3cols; i++)
		m3randcopy[i] = m3rand[i];
	Epetra_IntSerialDenseMatrix m3(Copy, m3rand, m3rows, m3rows, m3cols);
	if(debug) printMat("m3",m3);
	if(verbose) cout << "checking to see if data initialized correctly" << endl;
	EPETRA_TEST_ERR(!(m3.CV() == Copy), ierr);
	EPETRA_TEST_ERR(!(m3.M() == m3rows), ierr);
	EPETRA_TEST_ERR(!(m3.N() == m3cols), ierr);
	EPETRA_TEST_ERR(!(m3.LDA() == m3rows), ierr);
	EPETRA_TEST_ERR(!(m3.A() != m3rand), ierr); // should not be a view
	for(int i = 0; i < m3rows; i++)
		for(int j = 0; j < m3cols; j++)
			EPETRA_TEST_ERR(!(m3[j][i] == m3rand[j * m3rows + i]), ierr); // data should be identical to user array
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nmodifying entry" << endl;
	m3[1][5] = m3[1][5] + 3; // magic numbers for which element to change
	for(int i = 0; i < m3rows * m3cols; i++)
		EPETRA_TEST_ERR(!(m3rand[i] == m3randcopy[i]), ierr); // user array should be unchanged
	m3rand[13] = m3rand[13] + 3; // same magic modification performed to user's array
	for(int i = 0; i < m3rows; i++)
		for(int j = 0; j < m3cols; j++)
			EPETRA_TEST_ERR(!(m3[j][i] == m3rand[j * m3rows + i]), ierr); // should equal user array with same modification performed
	if(debug) {
		printArray(m3rand, m3rows * m3cols);
		printMat("m3",m3);
	}
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	
	if(verbose) cout << "\nreshaping" << endl;
	m3.Reshape(m3Rrows, m3Rcols);
	EPETRA_TEST_ERR(!(m3.M() == m3Rrows), ierr);
	EPETRA_TEST_ERR(!(m3.N() == m3Rcols), ierr);
	EPETRA_TEST_ERR(!(m3.LDA() == m3Rrows), ierr);
	for(int i = 0; i < m3Rrows; i++)
		for(int j = 0; j < m3Rcols; j++)
			EPETRA_TEST_ERR(!(m3[j][i] == m3rand[j * m3rows + i]),ierr);
	if(debug) printMat("m3",m3);

	delete[] m3rand;
	delete[] m3randcopy;
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	
	if(verbose) cout << "\nChecking pointer on zero-sized matrix" << endl;
	int* before = m3.A();
	if(verbose) cout << "Reshaping to (4,0)" << endl;
	if(debug) cout << "Before = " << before << endl;
	EPETRA_TEST_ERR(!(before != 0), ierr);
	m3.Reshape(4,0);
	int* after = m3.A();
	EPETRA_TEST_ERR(!(after == 0), ierr);
	if(debug) cout << "After = " << after << endl;
	m3.Shape(3,3);
	before = m3.A();
	if(verbose) cout << "Shaping to (0,3)" << endl;
	if(debug) cout << "Before = " << before << endl;
	EPETRA_TEST_ERR(!(before != 0), ierr);
	m3.Shape(0,3);
	after = m3.A();
	EPETRA_TEST_ERR(!(after == 0), ierr);
	if(debug) cout << "After = " << after << endl;
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	return(returnierr);
}

//=========================================================================
// test matrix copy constructor (copy & view)
int matrixCpyCtr(bool verbose, bool debug) {
	const int m1rows = 5;
	const int m1cols = 4;
	const int m2rows = 2;
	const int m2cols = 6;

	int ierr = 0;
	int returnierr = 0;
	if(verbose) printHeading("Testing matrix copy constructors");

	if(verbose) cout << "checking copy constructor (view)" << endl;
	int* m1rand = getRandArray(m1rows * m1cols);
	if(debug) printArray(m1rand, m1rows * m1cols);
	Epetra_IntSerialDenseMatrix m1(View, m1rand, m1rows, m1rows, m1cols);
	if(debug) {
		cout << "original matrix:" << endl;
		printMat("m1",m1);
	}
	Epetra_IntSerialDenseMatrix m1clone(m1);
	if(debug) {
		cout << "clone matrix:" << endl;
		printMat("m1clone",m1clone);
	}
	if(verbose) cout << "making sure signatures match" << endl;
	EPETRA_TEST_ERR(!identicalSignatures(m1, m1clone), ierr);
	delete[] m1rand;
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	
	if(verbose) cout << "\nchecking copy constructor (copy)" << endl;
	int* m2rand = getRandArray(m2rows * m2cols);
	if(debug) printArray(m2rand, m2rows * m2cols);
	Epetra_IntSerialDenseMatrix m2(Copy, m2rand, m2rows, m2rows, m2cols);
	if(debug) {
		cout << "original matrix:" << endl;
		printMat("m2",m2);
	}
	Epetra_IntSerialDenseMatrix m2clone(m2);
	if(debug) {
		cout << "clone matrix:" << endl;
		printMat("m2clone",m2clone);
	}
	if(verbose) cout << "checking that signatures match" << endl;
	EPETRA_TEST_ERR(!identicalSignatures(m2, m2clone), ierr);
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nmodifying entry in m2, m2clone should be unchanged" << endl;
	EPETRA_TEST_ERR(!seperateData(m2, m2clone), ierr);
	if(debug) {
		printArray(m2rand, m2rows * m2cols);
		cout << "orig:" << endl;
		printMat("m2",m2);
		cout << "clone:" << endl;
		printMat("m2clone",m2clone);
	}
	delete[] m2rand;
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	
	return(returnierr);
}

//=========================================================================
// test matrix error-reporting
int matrixExceptions(bool verbose, bool debug) {
	int returnierr = 0;
	int ierr = 0;
	bool caught = false;
	Epetra_IntSerialDenseMatrix* matrix;

	if(verbose) printHeading("Testing matrix error-reporting.\nExpect error messages if EPETRA_NO_ERROR_REPORTS is not defined.");
	
	// invalid dimension to sized ctr (2 cases)
	try {
		caught = false;
		if(verbose) cout << "Checking Epetra_IntSerialDenseMatrix(-1, 6) - invalid rows";
		matrix = new Epetra_IntSerialDenseMatrix(-1, 6);
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -1, returnierr);
		if(error == -1)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);
	try {
		caught = false;
		if(verbose) cout << "\nChecking Epetra_IntSerialDenseMatrix(3, -5) - invalid cols";
		matrix = new Epetra_IntSerialDenseMatrix(3, -5);
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -1, returnierr);
		if(error == -1)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);

	// invalid dimension to user-data ctr (3 cases)
	try {
		caught = false;
		if(verbose) cout << "\nChecking Epetra_IntSerialDenseMatrix(Copy, int*, -1, 2, 2) - invalid lda";
		int* rand2 = getRandArray(2);
		matrix = new Epetra_IntSerialDenseMatrix(Copy, rand2, -1, 2, 2);
		delete[] rand2;
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -1, returnierr);
		if(error == -1)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);
	try {
		caught = false;
		if(verbose) cout << "\nChecking Epetra_IntSerialDenseMatrix(Copy, int*, 3, -2, 3) - invalid rows";
		int* rand2 = getRandArray(2);
		matrix = new Epetra_IntSerialDenseMatrix(Copy, rand2, 3, -2, 3);
		delete[] rand2;
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -1, returnierr);
		if(error == -1)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);
	try {
		caught = false;
		if(verbose) cout << "\nChecking Epetra_IntSerialDenseMatrix(Copy, int*, 4, 4, -4) - invalid cols";
		int* rand2 = getRandArray(2);
		matrix = new Epetra_IntSerialDenseMatrix(Copy, rand2, -4, 4, -4);
		delete[] rand2;
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -1, returnierr);
		if(error == -1)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);
	
	// null pointer to user-data ctr
	try {
		caught = false;
		if(verbose) cout << "\nChecking Epetra_IntSerialDenseMatrix(Copy, 0, 5, 5, 5) - null pointer";
		matrix = new Epetra_IntSerialDenseMatrix(Copy, 0, 5, 5, 5);
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -3, returnierr);
		if(error == -3)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);

	// invalid parameter to shape (2 cases)
	Epetra_IntSerialDenseMatrix m1;
	if(verbose) cout << "\nChecking Shape(-2, 2) - invalid rows" << endl;
	ierr = m1.Shape(-2, 2);
	EPETRA_TEST_ERR(!(ierr == -1), returnierr);
	if(ierr == -1)
		if(verbose) cout << "Checked OK." << endl;
	if(verbose) cout << "\nChecking Shape(3, -2) - invalid cols" << endl;
	ierr = m1.Shape(3, -2);
	EPETRA_TEST_ERR(!(ierr == -1), returnierr);
	if(ierr == -1)
		if(verbose) cout << "Checked OK." << endl;
	
	// invalid parameter to reshape (2 cases)
	m1.Shape(5, 5);
	if(verbose) cout << "\nChecking Reshape(-4, 3) - invalid rows" << endl;
	ierr = m1.Reshape(-4, 3);
	EPETRA_TEST_ERR(!(ierr == -1), returnierr);
	if(ierr == -1)
		if(verbose) cout << "Checked OK." << endl;
	if(verbose) cout << "\nChecking Reshape(4, -3) - invalid cols" << endl;
	ierr = m1.Reshape(4, -3);
	EPETRA_TEST_ERR(!(ierr == -1), returnierr);
	if(ierr == -1)
		if(verbose) cout << "Checked OK." << endl;

#ifdef HAVE_EPETRA_ARRAY_BOUNDS_CHECK // only test op() and op[] exceptions if macro is defined.
	// out of range index to op() & op[] (6 cases)
	int* rand16 = getRandArray(16);
	Epetra_IntSerialDenseMatrix m2(View, rand16, 4, 4, 4);

	// op() too high
	try {
		caught = false;
		if(verbose) cout << "\nChecking operator () - row index too high";
		ierr = m2(5, 2);
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -1, returnierr);
		if(error == -1)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);
	try {
		caught = false;
		if(verbose) cout << "\nChecking operator () - col index too high";
		ierr = m2(3, 4);
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -2, returnierr);
		if(error == -2)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);

	// op() too low
	try {
		caught = false;
		if(verbose) cout << "\nChecking operator () - row index too low";
		ierr = m2(-1, 0);
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -1, returnierr);
		if(error == -1)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);
	try {
		caught = false;
		if(verbose) cout << "\nChecking operator () - col index too low";
		ierr = m2(0, -1);
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -2, returnierr);
		if(error == -2)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);

	// op[] too high
	try { 
		caught = false;
		if(verbose) cout << "\nChecking operator [] - col index too high";
		ierr = m2[4][2];
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -2, returnierr);
		if(error == -2)
			if(verbose) cout << "Checked OK." << endl;
	}

	// op[] too low
	try {
		caught = false;
		if(verbose) cout << "\nChecking operator [] - col index too low";
		ierr = m2[-1][0];
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -2, returnierr);
		if(error == -2)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);
#endif // end of HAVE_EPETRA_ARRAY_BOUNDS_CHECK conditional
	
	// ISDM = ISDV
	Epetra_IntSerialDenseMatrix m3;
	Epetra_IntSerialDenseVector v1;
	try {
		caught = false;
		if(verbose) cout << "\nChecking op = - assigning ISDV to ISDM";
		m3 = v1;
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -5, returnierr);
		if(error == -5)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);
	
	return(returnierr);
}

//=========================================================================
// test matrix operator= (copy & view)
int matrixAssignment(bool verbose, bool debug) {
	int ierr = 0;
	int returnierr = 0;
	if(verbose) printHeading("Testing matrix operator=");

	// each section is in its own block so we can reuse variable names
	// lhs = left hand side, rhs = right hand side
	
	{
		// copy->copy (more space needed)
		// orig and dup should have same signature
		// modifying orig or dup should have no effect on the other
		if(verbose) cout << "Checking copy->copy (new alloc)" << endl;
		Epetra_IntSerialDenseMatrix lhs(2,2);
		int* rand1 = getRandArray(25);
		Epetra_IntSerialDenseMatrix rhs(Copy, rand1, 5, 5, 5);
		if(debug) {
			cout << "before assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		lhs = rhs;
		if(debug) {
			cout << "after assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		EPETRA_TEST_ERR(!identicalSignatures(rhs,lhs), ierr);
		EPETRA_TEST_ERR(!seperateData(rhs,lhs), ierr);
		delete[] rand1;
	}
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	{
		// copy->copy (have enough space)
		// orig and dup should have same signature
		// modifying orig or dup should have no effect on the other
		if(verbose) cout << "\nChecking copy->copy (no alloc)" << endl;
		int* rand1 = getRandArray(25);
		int* rand2 = getRandArray(20);
		Epetra_IntSerialDenseMatrix lhs(Copy, rand1, 5, 5, 5);
		Epetra_IntSerialDenseMatrix rhs(Copy, rand2, 4, 4, 5);
		int* origA = lhs.A();
		int origLDA = lhs.LDA();
		if(debug) {
			cout << "before assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		lhs = rhs;
		if(debug) {
			cout << "after assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		// in this case, instead of doing a "normal" LDA test in identSig,
		// we do our own. Since we had enough space already, A and LDA should
		// not have been changed by the assignment. (The extra parameter to
		// identicalSignatures tells it not to test LDA).
		EPETRA_TEST_ERR((lhs.A() != origA) || (lhs.LDA() != origLDA), ierr);
		EPETRA_TEST_ERR(!identicalSignatures(rhs,lhs,false), ierr);
		EPETRA_TEST_ERR(!seperateData(rhs,lhs), ierr);
	}
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	{
		// view->copy
		// orig and dup should have same signature
		// modifying orig or dup should have no effect on the other
		if(verbose) cout << "\nChecking view->copy" << endl;
		int* rand1 = getRandArray(25);
		int* rand2 = getRandArray(64);
		Epetra_IntSerialDenseMatrix lhs(View, rand1, 5, 5, 5);
		Epetra_IntSerialDenseMatrix rhs(Copy, rand2, 8, 8, 8);
		if(debug) {
			cout << "before assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		lhs = rhs;
		if(debug) {
			cout << "after assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		EPETRA_TEST_ERR(!identicalSignatures(rhs,lhs), ierr);
		EPETRA_TEST_ERR(!seperateData(rhs,lhs), ierr);
		delete[] rand1;
		delete[] rand2;
	}
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	{
	  // copy->view
		// orig and dup should have same signature
		// modifying orig or dup should change the other
		if(verbose) cout << "\nChecking copy->view" << endl;
		int* rand1 = getRandArray(10);
		Epetra_IntSerialDenseMatrix lhs(4,4);
		Epetra_IntSerialDenseMatrix rhs(View, rand1, 2, 2, 5);
		if(debug) {
			cout << "before assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		lhs = rhs;
		if(debug) {
			cout << "after assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		EPETRA_TEST_ERR(!identicalSignatures(rhs,lhs), ierr);
		EPETRA_TEST_ERR(seperateData(rhs,lhs), ierr);
		delete[] rand1;
	}
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;	
	{
		// view->view
		// orig and dup should have same signature
		// modifying orig or dup should change the other
		if(verbose) cout << "\nChecking view->view" << endl;
		int* rand1 = getRandArray(9);
		int* rand2 = getRandArray(18);
		Epetra_IntSerialDenseMatrix lhs(View, rand1, 3, 3, 3);
		Epetra_IntSerialDenseMatrix rhs(View, rand2, 3, 3, 6);
		if(debug) {
			cout << "before assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		lhs = rhs;
		if(debug) {
			cout << "after assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		EPETRA_TEST_ERR(!identicalSignatures(rhs,lhs), ierr);
		EPETRA_TEST_ERR(seperateData(rhs,lhs), ierr);
		delete[] rand1;
		delete[] rand2;
	}
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	{
		// test MakeViewOf
		// orig and dup should have same signature except for CV_
		// modifying orig or dup should change the other
		if(verbose) cout << "\nChecking CrsGraph's usage of MakeViewOf" << endl;
		int* rand1 = getRandArray(10);
		int* rand2 = getRandArray(10);
		Epetra_IntSerialDenseMatrix lhs(Copy, rand1, 5, 5, 2);
		Epetra_IntSerialDenseMatrix rhs(Copy, rand2, 5, 5, 2);
		if(debug) {
			cout << "before assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		lhs.MakeViewOf(rhs);
		if(debug) {
			cout << "after assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		EPETRA_TEST_ERR(!(lhs.CV() == View), ierr);
		EPETRA_TEST_ERR(!(lhs.M() == rhs.M()), ierr);
		EPETRA_TEST_ERR(!(lhs.N() == rhs.N()), ierr);
		EPETRA_TEST_ERR(!(lhs.LDA() == rhs.LDA()), ierr);
		EPETRA_TEST_ERR(!(lhs.A() == rhs.A()), ierr);
		EPETRA_TEST_ERR(seperateData(rhs,lhs), ierr);
		delete[] rand1;
		delete[] rand2;
	}
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	return(returnierr);
}

//=========================================================================
//=========================================================================
// Vector-testing functions
//=========================================================================
//=========================================================================
// call functions we don't have tests for, for sake of code coverage.
// (the functions are called, but their output isn't checked for correctness).
int vectorCoverage(bool verbose, bool debug) {
	if(verbose) printHeading("Testing other vector functions");

	int* rand1 = getRandArray(9);
	if(debug) printArray(rand1, 9);
	Epetra_IntSerialDenseVector v1(Copy, rand1, 9);
	if(debug) printMat("v1",v1);
	delete[] rand1;

	if(verbose) cout << "calling one norm" << endl;
	int onenorm = v1.OneNorm();
	if(debug) cout << "v1.OneNorm() = " << onenorm << endl;

	if(verbose) cout << "calling infinity norm" << endl;
	int infnorm = v1.InfNorm();
	if(debug) cout << "v1.InfNorm() = " << infnorm << endl;

	if(verbose) cout << "calling random" << endl;
	Epetra_IntSerialDenseVector v2(12);
	if(debug) printMat("v2 (before)",v2);
	v2.Random();
	if(debug) printMat("v2 (after)",v2);

	if(verbose) cout << "Checked OK (I think)." << endl;

	return(0);
}
//=========================================================================
// test vector default constructor, user constructor (copy & view),
// () [] operators (read & write), size (larger), resize (smaller)
int vectorCtr(bool verbose, bool debug) {
	const int v1size = 5;
	const int v1asize = 15;
	const int v2size = 10;
	const int v3size = 8;
	const int v3resize = 5; // should be smaller than v3size

	int ierr = 0;
	int returnierr = 0;
	if(verbose) printHeading("Testing vector constructors");

	if(verbose) cout << "default constructor" << endl;
	Epetra_IntSerialDenseVector v1;
	EPETRA_TEST_ERR(!(v1.CV() == Copy), ierr);
	if(verbose) cout << "sizing" << endl;
	v1.Size(v1size);
	EPETRA_TEST_ERR(!(v1.Length() == v1size), ierr);
	for(int i = 0; i < v1size; i++) {
		EPETRA_TEST_ERR(!(v1(i) == 0), ierr);
	}
	if(debug) printMat("v1",v1);
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	if(verbose) cout << "\nmanually setting values" << endl;
	int* v1rand = getRandArray(v1size);
	for(int i = 0; i < v1size; i++)
		v1(i) = v1rand[i];
	for(int i = 0; i < v1size; i++)
		EPETRA_TEST_ERR(!(v1[i] == v1rand[i]), ierr);
	if(debug) {
		printArray(v1rand, v1size);
		printMat("v1",v1);
	}
	delete[] v1rand;
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nsized constructor" << endl;
	Epetra_IntSerialDenseVector v1a(v1asize);
	EPETRA_TEST_ERR(!(v1a.Length() == v1asize), ierr);
	EPETRA_TEST_ERR(!(v1a.CV() == Copy),ierr);
	if(debug) printMat("v1a", v1a);
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nuser data constructor (view)" << endl;
	int* v2rand = getRandArray(v2size);
	if(debug) printArray(v2rand, v2size);
	Epetra_IntSerialDenseVector v2(View, v2rand, v2size);
	EPETRA_TEST_ERR(!(v2.CV() == View), ierr);
	EPETRA_TEST_ERR(!(v2.Length() == v2size), ierr);
	EPETRA_TEST_ERR(!(v2.Values() == v2rand), ierr);
	if(debug) printMat("v2",v2);
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nchanging value, checking for correct behavior" << endl;
	int* v2randcopy = new int[v2size]; // backup of original values
	for(int i = 0; i < v2size; i++)
		v2randcopy[i] = v2rand[i];
	v2(4) = v2(4) + 1; // magic number for which element to change
	v2randcopy[4] = v2randcopy[4] +1;
	for(int i = 0; i < v2size; i++)
		EPETRA_TEST_ERR(!(v2rand[i] == v2randcopy[i]), ierr); // v2rand should have updated correctly
	if(debug) {
		printArray(v2rand, v2size);
		printMat("v2",v2);
	}
	delete[] v2rand;
	delete[] v2randcopy;
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nuser data constructor (copy)" << endl;
	int* v3rand = getRandArray(v3size);
	if(debug) printArray(v3rand, v3size);
	int* v3randcopy = new int[v3size];
	for(int i = 0; i < v3size; i++)
		v3randcopy[i] = v3rand[i];
	Epetra_IntSerialDenseVector v3(Copy, v3rand, v3size);
	if(debug) printMat("v3",v3);
	if(verbose) cout << "checking to see if data initialized correctly" << endl;
	EPETRA_TEST_ERR(!(v3.CV() == Copy), ierr);
	EPETRA_TEST_ERR(!(v3.Length() == v3size), ierr);
	EPETRA_TEST_ERR(!(v3.Values() != v3rand), ierr); // should not be a view
	for(int i = 0; i < v3size; i++)
		EPETRA_TEST_ERR(!(v3[i] == v3rand[i]), ierr); // data should be identical to user array
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nmodifying entry" << endl;
	v3[5] = v3[5] + 3; // magic number for which element to change
	for(int i = 0; i < v3size; i++)
		EPETRA_TEST_ERR(!(v3rand[i] == v3randcopy[i]), ierr); // user array should be unchanged
	v3rand[5] = v3rand[5] + 3; // same magic modification performed to user's array
	for(int i = 0; i < v3size; i++)
		EPETRA_TEST_ERR(!(v3[i] == v3rand[i]), ierr); // should equal user array with same modification performed
	if(debug) {
		printArray(v3rand, v3size);
		printMat("v3",v3);
	}
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nresizing" << endl;
	v3.Resize(v3resize);
	EPETRA_TEST_ERR(!(v3.Length() == v3resize), ierr);
	for(int i = 0; i < v3resize; i++)
		EPETRA_TEST_ERR(!(v3[i] == v3rand[i]),ierr);
	if(debug) printMat("v3",v3);
	delete[] v3rand;
	delete[] v3randcopy;
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nChecking pointer on zero-sized vector" << endl;
	int* before = v3.Values();
	if(verbose) cout << "Resizing to 0" << endl;
	if(debug) cout << "Before = " << before << endl;
	EPETRA_TEST_ERR(!(before != 0), ierr);
	v3.Resize(0);
	int* after = v3.Values();
	EPETRA_TEST_ERR(!(after == 0), ierr);
	if(debug) cout << "After = " << after << endl;
	v3.Size(3);
	before = v3.Values();
	if(verbose) cout << "Sizing to 0" << endl;
	if(debug) cout << "Before = " << before << endl;
	EPETRA_TEST_ERR(!(before != 0), ierr);
	v3.Size(0);
	after = v3.Values();
	EPETRA_TEST_ERR(!(after == 0), ierr);
	if(debug) cout << "After = " << after << endl;
	returnierr += ierr;
	if(ierr == 0) 
		 if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	return(returnierr);
}
//=========================================================================
// test vector copy constructor (copy & view)
int vectorCpyCtr(bool verbose, bool debug) {
	const int v1size = 15;
	const int v2size = 12;

	int ierr = 0;
	int returnierr = 0;
	if(verbose) printHeading("Testing vector copy constructors");

	if(verbose) cout << "checking copy constructor (view)" << endl;
	int* v1rand = getRandArray(v1size);
	if(debug) printArray(v1rand, v1size);
	Epetra_IntSerialDenseVector v1(View, v1rand, v1size);
	if(debug) {
		cout << "original vector:" << endl;
		printMat("v1",v1);
	}
	Epetra_IntSerialDenseVector v1clone(v1);
	if(debug) {
		cout << "clone vector:" << endl;
		printMat("v1clone",v1clone);
	}
	if(verbose) cout << "making sure signatures match" << endl;
	EPETRA_TEST_ERR(!identicalSignatures(v1, v1clone), ierr);
	delete[] v1rand;
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nchecking copy constructor (copy)" << endl;
	int* v2rand = getRandArray(v2size);
	if(debug) printArray(v2rand, v2size);
	Epetra_IntSerialDenseVector v2(Copy, v2rand, v2size);
	if(debug) {
		cout << "original vector:" << endl;
		printMat("v2",v2);
	}
	Epetra_IntSerialDenseVector v2clone(v2);
	if(debug) {
		cout << "clone vector:" << endl;
		printMat("v2clone",v2clone);
	}
	if(verbose) cout << "checking that signatures match" << endl;
	EPETRA_TEST_ERR(!identicalSignatures(v2, v2clone), ierr);
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nmodifying entry in v2, v2clone should be unchanged" << endl;
	EPETRA_TEST_ERR(!seperateData(v2, v2clone), ierr);
	if(debug) {
		printArray(v2rand, v2size);
		cout << "orig:" << endl;
		printMat("v2",v2);
		cout << "clone:" << endl;
		printMat("v2clone",v2clone);
	}
	delete[] v2rand;
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	return(returnierr);
}
//=========================================================================
// test vector operator= (copy & view)
int vectorAssignment(bool verbose, bool debug) {
	int ierr = 0;
	int returnierr = 0;
	if(verbose) printHeading("Testing vector operator=");

	// each section is in its own block so we can reuse variable names
	// lhs = left hand side, rhs = right hand side
	
	{
		// copy->copy (more space needed)
		// orig and dup should have same signature
		// modifying orig or dup should have no effect on the other
		if(verbose) cout << "Checking copy->copy (new alloc)" << endl;
		Epetra_IntSerialDenseVector lhs(2);
		int* rand1 = getRandArray(10);
		Epetra_IntSerialDenseVector rhs(Copy, rand1, 10);
		if(debug) {
			cout << "before assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		lhs = rhs;
		if(debug) {
			cout << "after assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		EPETRA_TEST_ERR(!identicalSignatures(rhs,lhs), ierr);
		EPETRA_TEST_ERR(!seperateData(rhs,lhs), ierr);
		delete[] rand1;
	}
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	{
		// copy->copy (have enough space)
		// orig and dup should have same signature
		// modifying orig or dup should have no effect on the other
		if(verbose) cout << "\nChecking copy->copy (no alloc)" << endl;
		int* rand1 = getRandArray(20);
		int* rand2 = getRandArray(15);
		Epetra_IntSerialDenseVector lhs(Copy, rand1, 20);
		Epetra_IntSerialDenseVector rhs(Copy, rand2, 15);
		int* origA = lhs.A();
		int origLDA = lhs.LDA();
		if(debug) {
			cout << "before assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		lhs = rhs;
		if(debug) {
			cout << "after assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		// in this case, instead of doing a "normal" LDA test in identSig,
		// we do our own. Since we had enough space already, A and LDA should
		// not have been changed by the assignment. (The extra parameter to
		// identicalSignatures tells it not to test LDA).
		EPETRA_TEST_ERR((lhs.A() != origA) || (lhs.LDA() != origLDA), ierr);
		EPETRA_TEST_ERR(!identicalSignatures(rhs,lhs,false), ierr);
		EPETRA_TEST_ERR(!seperateData(rhs,lhs), ierr);
	}
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	{
		// view->copy
		// orig and dup should have same signature
		// modifying orig or dup should have no effect on the other
		if(verbose) cout << "\nChecking view->copy" << endl;
		int* rand1 = getRandArray(5);
		int* rand2 = getRandArray(8);
		Epetra_IntSerialDenseVector lhs(View, rand1, 5);
		Epetra_IntSerialDenseVector rhs(Copy, rand2, 8);
		if(debug) {
			cout << "before assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		lhs = rhs;
		if(debug) {
			cout << "after assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		EPETRA_TEST_ERR(!identicalSignatures(rhs,lhs), ierr);
		EPETRA_TEST_ERR(!seperateData(rhs,lhs), ierr);
		delete[] rand1;
		delete[] rand2;
	}
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	{
	  // copy->view
		// orig and dup should have same signature
		// modifying orig or dup should change the other
		if(verbose) cout << "\nChecking copy->view" << endl;
		int* rand1 = getRandArray(10);
		Epetra_IntSerialDenseVector lhs(4);
		Epetra_IntSerialDenseVector rhs(View, rand1, 10);
		if(debug) {
			cout << "before assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		lhs = rhs;
		if(debug) {
			cout << "after assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		EPETRA_TEST_ERR(!identicalSignatures(rhs,lhs), ierr);
		EPETRA_TEST_ERR(seperateData(rhs,lhs), ierr);
		delete[] rand1;
	}
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	{
		// view->view
		// orig and dup should have same signature
		// modifying orig or dup should change the other
		if(verbose) cout << "\nChecking view->view" << endl;
		int* rand1 = getRandArray(9);
		int* rand2 = getRandArray(11);
		Epetra_IntSerialDenseVector lhs(View, rand1, 9);
		Epetra_IntSerialDenseVector rhs(View, rand2, 11);
		if(debug) {
			cout << "before assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		lhs = rhs;
		if(debug) {
			cout << "after assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		EPETRA_TEST_ERR(!identicalSignatures(rhs,lhs), ierr);
		EPETRA_TEST_ERR(seperateData(rhs,lhs), ierr);
		delete[] rand1;
		delete[] rand2;
	}
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	{
		// test MakeViewOf
		// orig and dup should have same signature except for CV_
		// modifying orig or dup should change the other
		if(verbose) cout << "\nChecking CrsGraph's usage of MakeViewOf" << endl;
		int* rand1 = getRandArray(10);
		int* rand2 = getRandArray(10);
		Epetra_IntSerialDenseVector lhs(Copy, rand1, 10);
		Epetra_IntSerialDenseVector rhs(Copy, rand2, 10);
		if(debug) {
			cout << "before assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		lhs.MakeViewOf(rhs);
		if(debug) {
			cout << "after assignment:" << endl;
			printMat("rhs",rhs);
			printMat("lhs",lhs);
		}
		EPETRA_TEST_ERR(!(lhs.CV() == View), ierr);
		EPETRA_TEST_ERR(!(lhs.M() == rhs.M()), ierr);
		EPETRA_TEST_ERR(!(lhs.N() == rhs.N()), ierr);
		EPETRA_TEST_ERR(!(lhs.LDA() == rhs.LDA()), ierr);
		EPETRA_TEST_ERR(!(lhs.A() == rhs.A()), ierr);
		EPETRA_TEST_ERR(seperateData(rhs,lhs), ierr);
		delete[] rand1;
		delete[] rand2;
	}
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	return(returnierr);
}

//=========================================================================
// test vector error-reporting
int vectorExceptions(bool verbose, bool debug) {
	int returnierr = 0;
	int ierr = 0;
	bool caught = false;
	Epetra_IntSerialDenseVector* vector;

	if(verbose) printHeading("Testing vector error-reporting.\nExpect error messages if EPETRA_NO_ERROR_REPORTS is not defined.");
	
	try { // invalid dimension to sized ctr
		caught = false;
		if(verbose) cout << "Checking Epetra_IntSerialDenseVector(-1)";
		vector = new Epetra_IntSerialDenseVector(-1);
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -1, returnierr);
		if(error == -1)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);

	try { // invalid dimension to user-data ctr
		caught = false;
		if(verbose) cout << "\nChecking Epetra_IntSerialDenseVector(Copy, int*, -3)";
		int* rand2 = getRandArray(2);
		vector = new Epetra_IntSerialDenseVector(Copy, rand2, -3);
		delete[] rand2;
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -1, returnierr);
		if(error == -1)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);

	try { // null pointer to user-data ctr
		caught = false;
		if(verbose) cout << "\nChecking Epetra_IntSerialDenseVector(Copy, 0, 5)";
		vector = new Epetra_IntSerialDenseVector(Copy, 0, 5);
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -3, returnierr);
		if(error == -3)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);

	// invalid parameter to size
	if(verbose) cout << "\nChecking Size(-2)" << endl;
	Epetra_IntSerialDenseVector v1;
	ierr = v1.Size(-2);
	EPETRA_TEST_ERR(!(ierr == -1), returnierr);
	if(ierr == -1)
		if(verbose) cout << "Checked OK." << endl;

	// invalid parameter to resize
	if(verbose) cout << "\nChecking Resize(-4)" << endl;
	v1.Size(5);
	ierr = v1.Resize(-4);
	EPETRA_TEST_ERR(!(ierr == -1), returnierr);
	if(ierr == -1)
		if(verbose) cout << "Checked OK." << endl;

#ifdef HAVE_EPETRA_ARRAY_BOUNDS_CHECK // only test op() and op[] exceptions if macro is defined.
	// out of range index to op() & op[]
	int* rand17 = getRandArray(17);
	Epetra_IntSerialDenseVector v2(View, rand17, 17);
	try { // op() too high
		caught = false;
		if(verbose) cout << "\nChecking operator () - index too high";
		ierr = v2(17);
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -1, returnierr);
		if(error == -1)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);
	
	try { // op() too low
		caught = false;
		if(verbose) cout << "\nChecking operator () - index too low";
		ierr = v2(-1);
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -1, returnierr);
		if(error == -1)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);

	try { // op[] too high
		caught = false;
		if(verbose) cout << "\nChecking operator [] - index too high";
		ierr = v2[17];
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -1, returnierr);
		if(error == -1)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);
	
	try { // op[] too low
		caught = false;
		if(verbose) cout << "\nChecking operator [] - index too low";
		ierr = v2[-1];
	}
	catch(int error) {
		caught = true;
		EPETRA_TEST_ERR(error != -1, returnierr);
		if(error == -1)
			if(verbose) cout << "Checked OK." << endl;
	}
	EPETRA_TEST_ERR(!caught, returnierr);
#endif // end of HAVE_EPETRA_ARRAY_BOUNDS_CHECK conditional

	// we don't need to check for ISDV = ISDM, as that is a compile-time error
	
	return(returnierr);
}

//=========================================================================
//=========================================================================
// helper functions
//=========================================================================
//=========================================================================
// checks the signatures of two matrices
bool identicalSignatures(Epetra_IntSerialDenseMatrix& a, Epetra_IntSerialDenseMatrix& b, bool testLDA) {
	/*cout << "M: " << a.M() << "  " << b.M() << endl;
	cout << "N: " << a.N() << "  " << b.N() << endl;
	cout << "LDA: " << a.LDA() << "  " << b.LDA() << endl;
	cout << "CV: " << a.CV() << "  " << b.CV() << endl;
	cout << "A: " << a.A() << "  " << b.A() << endl;*/

	if((a.M()  != b.M()  )|| // check properties first
		 (a.N()  != b.N()  )||
		 (a.CV() != b.CV() ))
		return(false);

	if(testLDA == true)      // if we are coming from op= c->c #2 (have enough space)
		if(a.LDA() != b.LDA()) // then we don't check LDA (but we do check it in the test function)
			return(false);

	if(a.CV() == View) { // if we're still here, we need to check the data
		if(a.A() != b.A()) // for a view, this just means checking the pointers
			return(false);   // for a copy, this means checking each element
	}
	else { // CV == Copy
		const int m = a.M();
		const int n = a.N();
		for(int i = 0; i < m; i++)
			for(int j = 0; j < n; j++) {
				if(a(i,j) != b(i,j))
					return(false);
			}
	}

	return(true); // if we're still here, signatures are identical
}
//=========================================================================
// checks if two matrices are independent or not
bool seperateData(Epetra_IntSerialDenseMatrix& a, Epetra_IntSerialDenseMatrix& b) {
	bool seperate;

	int r = EPETRA_MIN(a.M(),b.M()) / 2; // ensures (r,c) is valid
	int c = EPETRA_MIN(a.N(),b.N()) / 2; // in both matrices

	int orig_a = a(r,c);
	int new_value = a(r,c) + 1;
	if(b(r,c) == new_value) // there's a chance b could be independent, but
		new_value++;          // already have new_value in (r,c).
	
	a(r,c) = new_value;
	if(b(r,c) == new_value)
		seperate = false;
	else
		seperate = true;

	a(r,c) = orig_a; // undo change we made to a

	return(seperate);
}
//=========================================================================
// returns a int* array of a given length, with random values on interval [0,100]
int* getRandArray(int length) {
	int* array = new int[length];
	for(int i = 0; i < length; i++)
		array[i] = randomInt();

	return(array);
}
//=========================================================================
// returns a random integer on the interval [0-maxint).
// this is the same generator used in IntSerialDenseMatrix
int randomInt() {
	const int maxint = 100;

  const double a = 16807.0;
	const double BigInt = 2147483647.0;
	double seed = rand(); // Use POSIX standard random function;

	seed = fmod(a * seed, BigInt); // fmod returns remainder of floating point division
	                               // (a * seed) - (floor(a * seed / BigInt) * BigInt)
	double randdouble = (seed / BigInt); // should be [0,1)
	int randint = int(randdouble * maxint);
	
	return(randint);
}
//=========================================================================
// prints int* array with formatting
void printArray(int* array, int length) {
	cout << "user array (size " << length << "): ";
	for(int i = 0; i < length; i++)
		cout << array[i] << "  ";
	cout << endl;
}
//=========================================================================
// prints IntSerialDenseMatrix/Vector with formatting
void printMat(const char* name, Epetra_IntSerialDenseMatrix& matrix) {
	//cout << "--------------------" << endl;
	cout << "*** " << name << " ***" << endl;
	cout << matrix;
	//cout << "--------------------" << endl;
}

//=========================================================================
// prints section heading with spacers/formatting
void printHeading(const char* heading) {
	cout << "\n==================================================================\n";
	cout << heading << endl;
	cout << "==================================================================\n";
}
