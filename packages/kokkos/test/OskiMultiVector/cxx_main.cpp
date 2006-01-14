//@HEADER
// ************************************************************************
// 
//                 Kokkos: A Fast Kernel Package
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

#include "Kokkos_OskiMultiVector.hpp"
//#include "Kokkos_DenseVector.hpp"
#include "Kokkos_Version.hpp"

using namespace std;
using namespace Kokkos;

#define OTYPE int
#define STYPE double

typedef OskiMultiVector<OTYPE, STYPE> OMultiVector;
//typedef DenseVector<OTYPE, STYPE> DVector;

int main(int argc, char* argv[]) 
{

  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (verbose)  
     cout << Kokkos::Kokkos_Version() << endl << endl;  
  
  int numberFailedTests = 0;
  string testName = "";



  if (verbose) cout<<endl<<"********** CHECKING KOKKOS  OSKI MULTIVECTOR **********"<<endl<<endl;

  // OSKI must be initialized
  //extern "C" 
  oski_Init();

  // default constructor test
  OMultiVector A;
  if (verbose) cout <<"default constructor -- construct empty multiVector ";
  if ( A.getNumRows()!=0 || A.getNumCols()!=0 || 
       A.getIsStrided()!=false ||A.getRowInc()!=0 ||
       A.getColInc()!=0) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } else {
	if (verbose) cout << "successful."<<endl;
  }
  
  OTYPE numVectors = 10;
  OTYPE length = 100;

  STYPE ** values = new STYPE*[numVectors];
  STYPE * allValues = new STYPE[length*numVectors];
  for (int i=0; i<numVectors; i++) values[i] = allValues+i*length;
  for (int i=0; i<length*numVectors; i++) allValues[i] = (STYPE) i;

//      for (int i=0; i<numVectors; i++) {
//        cout << values[i] << endl;
//      }
//	cout << *(values[2]) << endl;
//  STYPE * Y = *values;
//      for (int i=0; i<numVectors; i++) {
//        cout << Y[i] << endl;
//      }


  A.initializeValues(numVectors, length, values);
  if (verbose) cout <<"default constructor -- initialize using array of pointers ";
cout << endl << "***" << endl << 
	A.getRowInc() << " 0" << endl <<
	A.getColInc() << " 1" << endl << 
	"***" << endl;

  if ( A.getNumRows()!=numVectors || A.getNumCols()!=length ||
       A.getIsStrided()!=false ||A.getRowInc()!=0 ||
       A.getColInc()!=1) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }

  OMultiVector B;
  B.initializeValues(numVectors, length, allValues, length, 1);
  if (verbose) cout <<"default constructor -- initialize using 2D array ";

cout << endl << "***" << endl <<
	B.getRowInc() << " " << length << endl <<
	B.getColInc() << " 1" << endl;

  if ( B.getNumRows()!=numVectors || B.getNumCols()!=length || 
       B.getIsStrided()!=true ||B.getRowInc()!=length ||
       B.getColInc()!=1) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }


  // constructor 2 (copy constructor)

  OMultiVector C(A);
  if(verbose) cout <<"constructor 2 -- copy constructor "; 

cout << endl << "***" << endl <<
        C.getRowInc() << " " << length << endl <<
        C.getColInc() << " 1" << endl;

  if ( C.getNumRows()!=numVectors || C.getNumCols()!=length || 
       C.getIsStrided()!=false || C.getRowInc()!=0 ||
       C.getColInc()!=1) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }

  bool valueTest = true;
  STYPE * Bvals;
  STYPE * Cvals;
  if(verbose) cout <<"constructor 2 -- testing values "; 
  for (int j=0; j<numVectors; j++) {
    Bvals = B.getValues(j);
    Cvals = C.getValues(j);
//    cout << "My Value Test: " << endl;
//    cout << *Bvals << " " << **(B.getValues()) << endl;
    for (int i=0; i<length; i++)//{
      if (Bvals[i]!=Cvals[i]) valueTest = false;
//      cout << Bvals[i] << "  ";}
  }
  if (!valueTest) {
    if (verbose) cout << "unsuccessful."<<endl;
    numberFailedTests++;
  } 
  else {
    if (verbose) cout << "successful."<<endl;
  }
  
  // OSKI can optionally be shut-down
  //extern "C" 
  oski_Close();

  //
  // If a test failed output the number of failed tests.
  //
  if(numberFailedTests > 0) cout << "Number of failed tests: " << numberFailedTests << endl;

 return numberFailedTests;
}  
