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
#include "Kokkos_Permutation.hpp"
#include "Kokkos_Version.hpp"

using namespace std;
using namespace Kokkos;

#define OTYPE int
#define STYPE double

typedef DenseMultiVector<OTYPE, STYPE> DMultiVector;
typedef DenseVector<OTYPE, STYPE> DVector;
typedef Permutation<OTYPE, STYPE> Perm;

int checkPerm(DMultiVector A, Perm P, DMultiVector AP, DVector v, DVector Pv, bool verbose);
int compareMultiVecs(const DMultiVector & v1, const DMultiVector & v2, bool verbose);
int compareVecs(const DVector & v1, const DVector & v2, bool verbose);


int main(int argc, char* argv[]) 
{

  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (verbose)  
     cout << Kokkos::Kokkos_Version() << endl << endl;  
  
  int numberFailedTests = 0;
  string testName = "";



  if (verbose) cout<<endl<<"********** CHECKING KOKKOS  PERMUTATIONS **********"<<endl<<endl;

  // default constructor test
  DMultiVector A;
  if (verbose) cout <<"default constructor -- construct empty matrix ";
  if ( A.getNumRows()!=0 || A.getNumCols()!=0 || 
       A.getIsStrided()!=false ||A.getRowInc()!=0 ||
       A.getColInc()!=0) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } else {
	if (verbose) cout << "successful."<<endl;
  }
  
  DVector v;
  if (verbose) cout <<"default constructor -- construct empty vector ";
  if ( v.getLength()!=0 || v.getInc()!=0) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } else {
	if (verbose) cout << "successful."<<endl;
  }
  
  OTYPE numVectors = 10;
  OTYPE length = 100;

  STYPE ** values = new STYPE*[numVectors];
  STYPE * allValues = new STYPE[length*numVectors];
  STYPE * allValuesReverse = new STYPE[length*numVectors];
  OTYPE * permIndices = new OTYPE[length];
  for (int i=0; i<numVectors; i++) values[i] = allValues+i*length;
  for (int i=0; i<length*numVectors; i++) allValues[i] = (STYPE) i;
  for (int i=0; i<length*numVectors; i++)  allValuesReverse[i] = allValues[length*numVectors-i-1];
  for (int i=0; i<length; i++)  permIndices[i] = length-i-1; // Reverse indices permutation

  A.initializeValues(length, numVectors, values);
  if (verbose) cout <<"default constructor -- initialize using array of pointers ";
  if ( A.getNumRows()!=length || A.getNumCols()!=numVectors || 
       A.getIsStrided()!=false ||A.getRowInc()!=0 ||
       A.getColInc()!=1) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }

  v.initializeValues(length, allValues);
  if (verbose) cout <<"default vector constructor -- initialize using array  ";
  if ( v.getLength()!=length || v.getInc()!=1) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }

  DMultiVector B;
  B.initializeValues(length, numVectors, allValuesReverse, length, 1);
  if (verbose) cout <<"default multivector constructor -- initialize using 2D array ";
  if ( B.getNumRows()!=length || B.getNumCols()!=numVectors || 
       B.getIsStrided()!=true ||B.getRowInc()!=length ||
       B.getColInc()!=1) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }
  
  DVector w;
  w.initializeValues(length, allValuesReverse+length*(numVectors-1));
  if (verbose) cout <<"default vector constructor -- initialize using array ";
  if ( w.getLength()!=length || w.getInc()!=1) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }
  
  Perm P;
  if (verbose) cout <<"default permutation constructor ";
  if ( P.getIsIdentity()!=true || P.getLength()!=0 || 
       P.getIndices()!=0) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }

  int ierr = checkPerm(A, P, A, v, v, verbose); // Check identity permutation
  if (ierr!=0) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }
 
  if (verbose) cout <<"post-construction initialization ";
  P.initialize(length, permIndices);
  if ( P.getIsIdentity()==true || P.getLength()!=length || 
       P.getIndices()!=permIndices) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }
  ierr += checkPerm(A, P, B, v, w, verbose); // Check reverse permutation
  if (ierr!=0) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }
 
   Perm P1(length, permIndices);
  if (verbose) cout <<"Non-identity constructor ";
  if ( P1.getIsIdentity()==true || P1.getLength()!=length || 
       P1.getIndices()!=permIndices) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }
  ierr += checkPerm(A, P1, B, v, w, verbose); // Check reverse permutation
  if (ierr!=0) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }
  

  //
  // If a test failed output the number of failed tests.
  //
  if(numberFailedTests > 0) cout << "Number of failed tests: " << numberFailedTests << endl;

  delete [] values;
  delete [] allValues;
  delete [] allValuesReverse;
  delete [] permIndices;
 return 0;
}

int checkPerm(DMultiVector A, Perm P, DMultiVector AP, DVector v, DVector Pv, bool verbose) {
  DMultiVector B;
  STYPE * BV = new STYPE[A.getNumRows()*A.getNumCols()];
  B.initializeValues(A.getNumRows(), A.getNumCols(), BV, A.getNumRows());
  P.apply(A, B);
  int ierr = compareMultiVecs(B, AP, verbose);

  DVector x;
  STYPE * xv = new STYPE[v.getLength()];
  x.initializeValues(v.getLength(), xv);
  P.apply(v, x);
  ierr += compareVecs(x, Pv, verbose);

  delete [] BV;
  delete [] xv;
  return(ierr);
}
int compareMultiVecs(const DMultiVector & v1, const DMultiVector & v2, bool verbose) {

  STYPE sum = 0.0;
  for (OTYPE i=0; i<v1.getNumCols(); i++) {
    STYPE * v1v = v1.getValues(i);
    STYPE * v2v = v2.getValues(i);
    for (OTYPE i=0; i<v1.getNumRows(); i++) sum += v1v[i] - v2v[i];
  }

  if (verbose) cout << "Difference between exact and computed multivectors = " << sum << endl;
  if (!((sum*sum)<1.0E-4)) {
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
  if (!((sum*sum)<1.0E-4)) {
    if (verbose) cout << "********** Difference too large **********" << endl;
    return(1);
  }
  return(0);
}

