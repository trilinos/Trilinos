//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
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

#include "Kokkos_DenseMultiVector.hpp"
#include "Kokkos_DenseVector.hpp"
#include "Kokkos_HbMatrix.hpp"
#include "Kokkos_BaseSparseMultiply.hpp"
#include "Kokkos_PackedSparseMultiply.hpp"
#include "Kokkos_BaseSparseSolve.hpp"
#include "Kokkos_Time.hpp"
#include "Kokkos_Flops.hpp"
#include "GenerateHbProblem.hpp"

using namespace std;
using namespace Kokkos;

#define OTYPE int
#define STYPE double
#define MULTCLASS BaseSparseMultiply
#define SOLVECLASS BaseSparseSolve

typedef MultiVector<OTYPE, STYPE> DMultiVector;
typedef Vector<OTYPE, STYPE> DVector;
typedef CisMatrix<OTYPE, STYPE> DHbMatrix;

template<typename TYPE>
int PrintTestResults(string, TYPE, TYPE, bool);

template<typename OrdinalType, typename ScalarType>
void GenerateHbProblem(bool generateClassicHbMatrix, bool isRowOriented, bool hasImplicitUnitDiagonal,
		       OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
		       OrdinalType * xoff, OrdinalType * yoff,
		       Kokkos::CisMatrix<OrdinalType, ScalarType> *& A, 
		       Kokkos::Vector<OrdinalType, ScalarType> *& x, 
		       Kokkos::Vector<OrdinalType, ScalarType> *& b,
		       Kokkos::Vector<OrdinalType, ScalarType> *&xexact,
		       OrdinalType & numEntries);

template<typename OrdinalType, typename ScalarType>
void GenerateHbProblem(bool generateClassicHbMatrix, bool isRowOriented, bool hasImplicitUnitDiagonal,
		       OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
		       OrdinalType * xoff, OrdinalType * yoff, OrdinalType nrhs,
		       Kokkos::CisMatrix<OrdinalType, ScalarType> *& A, 
		       Kokkos::MultiVector<OrdinalType, ScalarType> *& x, 
		       Kokkos::MultiVector<OrdinalType, ScalarType> *& b,
		       Kokkos::MultiVector<OrdinalType, ScalarType> *&xexact,
		       OrdinalType & numEntries);

void runChecks (bool verbose, bool generateClassicHbMatrix, bool isRowOriented, 
		int numEquations, int numEntries, int nrhs,
		bool isLowerTriangular, bool isUpperTriangular, bool hasImplicitUnitDiagonal,
		bool canUseStructure, bool canUseValues,
		DHbMatrix * A, DMultiVector * x, DMultiVector * b, DMultiVector * xx, 
		int & numberFailedTests);

int compareMultiVecs(const DMultiVector & v1, const DMultiVector & v2, bool verbose);
int compareVecs(const DVector & v1, const DVector & v2, bool verbose);

int main(int argc, char* argv[]) 
{

  int i, j;
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  int numberFailedTests = 0;
  int returnCode = 0;
  string testName = "";

  DHbMatrix * A;
  DVector * xexact;
  DMultiVector * xm;
  DMultiVector * bm;
  DMultiVector * xexactm;
  OTYPE nrhs = 5;
  OTYPE nx = 250;
  OTYPE ny = nx;
  OTYPE npointsU = 4;
  OTYPE npointsL = 6;
  OTYPE npointsU1 = 5;
  OTYPE npointsL1 = 7;
 
  OTYPE xoffU[] = {1,  0, 1, 2};
  OTYPE yoffU[] = {0,  1, 2, 1}; 

  OTYPE xoffL[] = {-2, -1, -1,  0,  1, -1};
  OTYPE yoffL[] = {-1, -2, -1, -1, -1,  0};

  OTYPE xoffU1[] = {0, 1,  0, 1, 2};
  OTYPE yoffU1[] = {0, 0,  1, 2, 1};

  OTYPE xoffL1[] = {-2, -1, -1,  0,  1, -1,  0};
  OTYPE yoffL1[] = {-1, -2, -1, -1, -1,  0,  0};

  OTYPE numEquations = nx*ny;
  OTYPE numEntries; 
  if (verbose) cout << "Size of Ordinal Type (in bytes) = " << sizeof(OTYPE) << endl
		    << "Size of Scalar  Type (in bytes) = " << sizeof(STYPE) << endl;

  bool canUseStructure = true;
  bool canUseValues = true;

  bool hasImplicitUnitDiagonal = true;
  bool isUpper;

  bool generateClassicHbMatrix;
  bool isRowOriented;
  for (int ii=0; ii<2; ii++) {
    for (int jj=0; jj<2; jj++) {
      generateClassicHbMatrix = (ii==1);
      isRowOriented = (jj==1);
      KokkosTest::GenerateHbProblem<OTYPE, STYPE>
	problem(generateClassicHbMatrix, isRowOriented, hasImplicitUnitDiagonal,
		nx, ny, npointsU, xoffU, yoffU, 
		nrhs, A, xm, bm, xexactm, numEntries);
       
      isUpper = isRowOriented;
      runChecks(verbose, generateClassicHbMatrix, isRowOriented, numEquations, numEntries, nrhs, !isUpper, 
		isUpper, hasImplicitUnitDiagonal, canUseStructure, canUseValues, A, xm, bm, xexactm, numberFailedTests);
       
    }
  }

  for (int ii=0; ii<2; ii++) {
    for (int jj=0; jj<2; jj++) {
      generateClassicHbMatrix = (ii==1);
      isRowOriented = (jj==1);
      KokkosTest::GenerateHbProblem<OTYPE, STYPE>
	problem(generateClassicHbMatrix, isRowOriented, hasImplicitUnitDiagonal,
		nx, ny, npointsL, xoffL, yoffL, 
		nrhs, A, xm, bm, xexactm, numEntries);
       
      isUpper = !isRowOriented;
      runChecks(verbose, generateClassicHbMatrix, isRowOriented, numEquations, numEntries, nrhs, !isUpper, 
		isUpper, hasImplicitUnitDiagonal, canUseStructure, canUseValues, A, xm, bm, xexactm, numberFailedTests);
    
    }
  }
  hasImplicitUnitDiagonal = false;
  for (int ii=0; ii<2; ii++) {
    for (int jj=0; jj<2; jj++) {
      generateClassicHbMatrix = (ii==1);
      isRowOriented = (jj==1);
      KokkosTest::GenerateHbProblem<OTYPE, STYPE>
	problem(generateClassicHbMatrix, isRowOriented, hasImplicitUnitDiagonal,
		nx, ny, npointsL1, xoffL1, yoffL1, 
		nrhs, A, xm, bm, xexactm, numEntries);
       
      isUpper = !isRowOriented;
      runChecks(verbose, generateClassicHbMatrix, isRowOriented, numEquations, numEntries, nrhs, !isUpper, 
		isUpper, hasImplicitUnitDiagonal, canUseStructure, canUseValues, A, xm, bm, xexactm, numberFailedTests);
       
    }
  }
  for (int ii=0; ii<2; ii++) {
    for (int jj=0; jj<2; jj++) {
      generateClassicHbMatrix = (ii==1);
      isRowOriented = (jj==1);
      KokkosTest::GenerateHbProblem<OTYPE, STYPE>
	problem(generateClassicHbMatrix, isRowOriented, hasImplicitUnitDiagonal,
		nx, ny, npointsU1, xoffU1, yoffU1, 
		nrhs, A, xm, bm, xexactm, numEntries);

      isUpper = isRowOriented;
      runChecks(verbose, generateClassicHbMatrix, isRowOriented, numEquations, numEntries, nrhs, !isUpper, 
		isUpper, hasImplicitUnitDiagonal, canUseStructure, canUseValues, A, xm, bm, xexactm, numberFailedTests);
    
    }
  }
  //
  // If a test failed output the number of failed tests.
  //
  if(numberFailedTests > 0) cout << endl<< "Number of failed tests: " << numberFailedTests << endl;
  else cout << endl << "All Tests Passed!"<<endl;

  return 0;
}
void runChecks (bool verbose, bool generateClassicHbMatrix, bool isRowOriented, 
		int numEquations, int numEntries, int nrhs,
		bool isLowerTriangular, bool isUpperTriangular, bool hasImplicitUnitDiagonal,
		bool canUseStructure, bool canUseValues,
		DHbMatrix * A, DMultiVector * x, DMultiVector * b, DMultiVector * xx, 
		int & numberFailedTests) {
  if (verbose) {
    
    cout <<"********** CHECKING KOKKOS HbMatrix **********" 
	 << " Dim = " << numEquations << " Num Nonzeros = " <<numEntries <<endl<<endl;
    cout <<"  generateClassicHbMatrix = "<< generateClassicHbMatrix << endl;
    cout <<"  isRowOriented           = "<< isRowOriented << endl;
    cout <<"  isLowerTriangular       = "<< isLowerTriangular << endl;
    cout <<"  isUpperTriangular       = "<< isUpperTriangular << endl;
    cout <<"  hasImplicitUnitDiagonal = "<< hasImplicitUnitDiagonal << endl;
  }
  

  Kokkos::HbMatrix<OTYPE, STYPE> * HbA = dynamic_cast<Kokkos::HbMatrix<OTYPE, STYPE> *>(A);
  HbA->setHasDiagonalEntries((!hasImplicitUnitDiagonal));
  HbA->setHasImplicitUnitDiagonal(hasImplicitUnitDiagonal);
  HbA->setIsLowerTriangular(isLowerTriangular);
  HbA->setIsUpperTriangular(isUpperTriangular);
  
  // Check output objects
  if (verbose) cout <<"Checking Attribute accessors ......."<<endl<<endl;
  numberFailedTests += PrintTestResults<OTYPE>("getNumRows()", A->getNumRows(), numEquations, verbose);
  numberFailedTests += PrintTestResults<OTYPE>("getNumCols()", A->getNumCols(), numEquations, verbose);
  numberFailedTests += PrintTestResults<OTYPE>("getNumEntries()", A->getNumEntries(), numEntries, verbose);
  numberFailedTests += PrintTestResults<bool>("getIsRowOriented()", A->getIsRowOriented(), isRowOriented, verbose);
  numberFailedTests += PrintTestResults<bool>("getIsUpperTriangular()", A->getIsUpperTriangular(), isUpperTriangular, verbose);
  numberFailedTests += PrintTestResults<bool>("getIsLowerTriangular()", A->getIsLowerTriangular(), isLowerTriangular, verbose);
  numberFailedTests += PrintTestResults<bool>("getHasImplicitUnitDiagonal()", A->getHasImplicitUnitDiagonal(), hasImplicitUnitDiagonal, verbose);

  Kokkos::MULTCLASS<OTYPE, STYPE> multA;
  multA.initializeStructure(*A, true);
  multA.initializeValues(*A, true);
  
  Kokkos::SOLVECLASS<OTYPE, STYPE> solveA;
  solveA.initializeStructure(*A, true);
  solveA.initializeValues(*A, true);
  
  numberFailedTests += PrintTestResults<bool>("getCanUseStructure()", solveA.getCanUseStructure(), canUseStructure, verbose);
  numberFailedTests += PrintTestResults<bool>("getCanUseValues()", solveA.getCanUseValues(), canUseValues, verbose);
  numberFailedTests += PrintTestResults<OTYPE>("getLeftPermutation().getIsIdentity()", solveA.getLeftPermutation().getIsIdentity(), 
					       true, verbose);
  numberFailedTests += PrintTestResults<OTYPE>("getRightPermutation().getIsIdentity()", solveA.getRightPermutation().getIsIdentity(), 
					       true, verbose);
  Kokkos::Flops multCounter;
  Kokkos::Flops solveCounter;
  multA.setFlopCounter(multCounter);
  solveA.setFlopCounter(solveCounter);
  Kokkos::Time timer;
  int ntrials = 20;

  // First solve for Multiple RHS, then multiply
  
  double start = timer.elapsedTime();
  for (int iii=0; iii<ntrials; iii++) solveA.apply(*b, *x);
  double solveAtime = timer.elapsedTime() - start;
  double solveAflops = solveA.getFlops();
  double mflops = solveAflops/solveAtime/1000000.0;
  if (verbose) cout << "Solve MFLOPS = " << mflops << endl; 
  numberFailedTests += compareMultiVecs(*x, *xx, verbose);

  start = timer.elapsedTime();
  for (int iii=0; iii<ntrials; iii++) multA.apply(*xx, *x);
  double multAtime = timer.elapsedTime() - start;
  double multAflops = multA.getFlops();
  mflops = multAflops/multAtime/1000000.0;
  if (verbose) cout << "Multiply MFLOPS = " << mflops << endl; 
  numberFailedTests += compareMultiVecs(*x, *b, verbose);

  // Next solve for single RHS, then multiply

  DenseVector<OTYPE, STYPE> x1;
  DenseVector<OTYPE, STYPE> b1;
  DenseVector<OTYPE, STYPE> xx1;
  x1.initializeValues(x->getNumRows(), x->getValues(0)); 
  b1.initializeValues(b->getNumRows(), b->getValues(0)); 
  xx1.initializeValues(xx->getNumRows(), xx->getValues(0)); 
  
  start = timer.elapsedTime();
  for (int iii=0; iii<ntrials; iii++) solveA.apply(b1, x1);
  solveAtime = timer.elapsedTime() - start;
  solveAflops = solveA.getFlops(); 
  mflops = solveAflops/solveAtime/1000000.0;
  if (verbose) cout << "Solve MFLOPS = " << mflops << endl; 
  numberFailedTests += compareVecs(x1, xx1, verbose);
  
  start = timer.elapsedTime();
  for (int iii=0; iii<ntrials; iii++) multA.apply(xx1, x1);
  multAtime = timer.elapsedTime() - start;
  multAflops = multA.getFlops();
  mflops = multAflops/multAtime/1000000.0;
  if (verbose) cout << "Multiply MFLOPS = " << mflops << endl; 
  numberFailedTests += compareVecs(x1, b1, verbose); 
  return;
}

int compareMultiVecs(const DMultiVector & v1, const DMultiVector & v2, bool verbose) {
  
  STYPE sum = 0.0;
  for (OTYPE i=0; i<v1.getNumCols(); i++) {
    STYPE * v1v = v1.getValues(i);
    STYPE * v2v = v2.getValues(i);
    for (OTYPE i=0; i<v1.getNumRows(); i++) sum += v1v[i] - v2v[i];
  }
  
  if (verbose) cout << "Difference between exact and computed multivectors = " << sum << endl;
  if (!(abs(sum)<1.0E-4)) {
    if (verbose) cout << "********** Difference too large **********" << endl;
    return(1); 
  }
  return(0);
}

int compareVecs(const DVector & v1, const DVector & v2, bool verbose) {

  STYPE sum = 0.0;
  STYPE * v1v = v1.getValues();
  STYPE * v2v = v2.getValues();
  for (OTYPE i=0; i < v1.getLength(); i++) sum += v1v[i] - v2v[i];

  if (verbose) cout << "Difference between exact and computed vectors = " << sum << endl;  
  if (!(abs(sum)<1.0E-4)) {
    if (verbose) cout << "********** Difference too large **********" << endl;
    return(1); 
  }
  return(0);
}

template<typename TYPE>
int PrintTestResults(string testName, TYPE calculatedResult, TYPE expectedResult, bool verbose)
{
  int result;
  if(calculatedResult == expectedResult)
    {
      if(verbose) cout << testName << " successful." << endl;
      result = 0;
    }
  else
    {
      if(verbose) cout << testName << " unsuccessful." << endl;
      result = 1;
    }
  return result;
}

