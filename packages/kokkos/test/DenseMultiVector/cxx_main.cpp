#include <iostream>
#include <string>
#include "Kokkos_DenseMultiVector.hpp"
#include "Kokkos_DenseVector.hpp"

using namespace std;
using namespace Kokkos;

#define OTYPE int
#define STYPE double

typedef DenseMultiVector<OTYPE, STYPE> DMultiVector;
typedef DenseVector<OTYPE, STYPE> DVector;

int main(int argc, char* argv[]) 
{

  int i, j;
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  int numberFailedTests = 0;
  int returnCode = 0;
  string testName = "";



  if (verbose) cout<<endl<<"********** CHECKING KOKKOS  DENSE VECTOR/MULTIVECTOR **********"<<endl<<endl;

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
  
  OTYPE numVectors = 10;
  OTYPE length = 100;

  STYPE ** values = new STYPE*[numVectors];
  STYPE * allValues = new STYPE[length*numVectors];
  for (int i=0; i<numVectors; i++) values[i] = allValues+i*length;
  for (int i=0; i<length*numVectors; i++) allValues[i] = (STYPE) i;

  A.initializeValues(numVectors, length, values);
  if (verbose) cout <<"default constructor -- initialize using array of pointers ";
  if ( A.getNumRows()!=numVectors || A.getNumCols()!=length || 
       A.getIsStrided()!=false ||A.getRowInc()!=0 ||
       A.getColInc()!=1) {
	if (verbose) cout << "unsuccessful."<<endl;
	numberFailedTests++;
  } 
  else {
	if (verbose) cout << "successful."<<endl;
  }

  DMultiVector B;
  B.initializeValues(numVectors, length, allValues, length, 1);
  if (verbose) cout <<"default constructor -- initialize using 2D array ";
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

  DMultiVector C(A);
  if(verbose) cout <<"constructor 2 -- copy constructor "; 
  if ( C.getNumRows()!=numVectors || C.getNumCols()!=length || 
       C.getIsStrided()!=false ||C.getRowInc()!=0 ||
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
    for (int i=0; i<length; i++) 
      if (Bvals[i]!=Cvals[i]) valueTest = false;
  }
  if (!valueTest) {
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

 return 0;
}  
