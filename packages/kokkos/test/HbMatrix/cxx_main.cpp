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
#include "GenerateHbProblem.hpp"

using namespace std;
using namespace Kokkos;

#define OTYPE int
#define STYPE double

template<typename TYPE>
int PrintTestResults(string, TYPE, TYPE, bool);

int ReturnCodeCheck(string, int, int, bool);

template<typename OrdinalType, typename ScalarType>
void GenerateHbProblem(bool generateClassicHbMatrix, bool isRowOriented,
		       OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
		       OrdinalType * xoff, OrdinalType * yoff,
		       Kokkos::CisMatrix<OrdinalType, ScalarType> *& A, 
		       Kokkos::Vector<OrdinalType, ScalarType> *& x, 
		       Kokkos::Vector<OrdinalType, ScalarType> *& b,
		       Kokkos::Vector<OrdinalType, ScalarType> *&xexact,
		       OrdinalType & numEntries);

template<typename OrdinalType, typename ScalarType>
void GenerateHbProblem(bool generateClassicHbMatrix, bool isRowOriented,
		       OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
		       OrdinalType * xoff, OrdinalType * yoff, OrdinalType nrhs,
		       Kokkos::CisMatrix<OrdinalType, ScalarType> *& A, 
		       Kokkos::MultiVector<OrdinalType, ScalarType> *& x, 
		       Kokkos::MultiVector<OrdinalType, ScalarType> *& b,
		       Kokkos::MultiVector<OrdinalType, ScalarType> *&xexact,
		       OrdinalType & numEntries);

typedef MultiVector<OTYPE, STYPE> DMultiVector;
typedef Vector<OTYPE, STYPE> DVector;
typedef CisMatrix<OTYPE, STYPE> DHbMatrix;

int main(int argc, char* argv[]) 
{

  int i, j;
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  int numberFailedTests = 0;
  int returnCode = 0;
  string testName = "";

  DHbMatrix * A;
  DVector * x;
  DVector * b;
  DVector * xexact;
  OTYPE nx = 10;
  OTYPE ny = nx;
  OTYPE npoints = 7;
 
  OTYPE xoff[] = {-1,  0,  1, -1,  0,  1,  0};
  OTYPE yoff[] = {-1, -1, -1,  0,  0,  0,  1};

  OTYPE numEquations = nx*ny;
  OTYPE numEntries; 
  {
    bool generateClassicHbMatrix = true;
    bool isRowOriented = true;
    KokkosTest::GenerateHbProblem<OTYPE, STYPE>
      problem(generateClassicHbMatrix, isRowOriented, nx, ny, npoints, xoff, yoff, A, x, b, xexact, numEntries);
    
    if (verbose) cout<<endl<<"********** CHECKING KOKKOS  Classic HbMatrix **********"<<endl<<endl;
    
    // Check output objects
    if (verbose) cout <<"Checking Attribute accessors ";
    if ( A->getNumRows()!=numEquations || A->getNumCols()!=numEquations || 
	 A->getIsRowOriented()!=isRowOriented ||A->getNumEntries()!=numEntries) {
      if (verbose) cout << "unsuccessful."<<endl;
      numberFailedTests++;
    } else {
      if (verbose) cout << "successful."<<endl;
    }
  }
  {
    bool generateClassicHbMatrix = false;
    bool isRowOriented = false;
    
    KokkosTest::GenerateHbProblem<OTYPE, STYPE>
      problem(generateClassicHbMatrix, isRowOriented, nx, ny, npoints, xoff, yoff, A, x, b, xexact, numEntries);
    
    if (verbose) cout<<endl<<"********** CHECKING KOKKOS  Generalized HbMatrix **********"<<endl<<endl;
    
    // Check output objects
    if (verbose) cout <<"Checking Attribute accessors ";
    if ( A->getNumRows()!=numEquations || A->getNumCols()!=numEquations || 
	 A->getIsRowOriented()!=isRowOriented ||A->getNumEntries()!=numEntries) {
      if (verbose) cout << "unsuccessful."<<endl;
      numberFailedTests++;
    } else {
      if (verbose) cout << "successful."<<endl;
    }
  }
  //
  // If a test failed output the number of failed tests.
  //
  if(numberFailedTests > 0) cout << "Number of failed tests: " << numberFailedTests << endl;

 return 0;
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

int ReturnCodeCheck(string testName, int returnCode, int expectedResult, bool verbose)
{
  int result;
  if(expectedResult == 0)
    {
      if(returnCode == 0)
	{
	  if(verbose) cout << testName << " test successful." << endl;
	  result = 0;
	}
      else
	{
	  if(verbose) cout << testName << " test unsuccessful. Return code was " << returnCode << "." << endl;
	  result = 1;
	}
    }
  else
    {
      if(returnCode != 0)
	{
	  if(verbose) cout << testName << " test successful -- failed as expected." << endl;
	  result = 0;
	}
      else
	{
	  if(verbose) cout << testName << " test unsuccessful -- did not fail as expected. Return code was " << returnCode << "." << endl;
	  result = 1;
	}
    }
  return result;
}
