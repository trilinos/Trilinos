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

// Tpetra Vector tester
// Modified: 06-Feb-2003

#include "Tpetra_ConfigDefs.hpp" // for <iostream> and <stdlib>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_Vector.hpp"

// function prototype
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug);
void checkOutputs();

int main(int argc, char* argv[]) {
	// initialize verbose & debug flags
	bool verbose = false;
	bool debug = false;
	if(argc > 1) {
		if(argv[1][0] == '-' && argv[1][1] == 'v')
			verbose = true;
		if(argv[1][0] == '-' && argv[1][1] == 'd') {
			debug = true;
			verbose = true;
		}
	}

	// call test routine
	int ierr = 0;
	if(verbose) cout << "Starting VectorTest..." << endl;
	ierr += unitTests<int, float>(verbose, debug);
	ierr += unitTests<int, double>(verbose, debug);
    
	if(debug) checkOutputs();

	// finish up
	if(verbose) 
		if(ierr == 0)
			cout << "Vector test successful." << endl;
		else
			cout << "Vector test failed." << endl;
	return(ierr);
}

//======================================================================
void checkOutputs() {
    cout << "Doing checkOutput..." << endl;
    int const length = 5;
	int const indexBase = 0;
    const Tpetra::SerialPlatform<int, int> platformE;
	Tpetra::ElementSpace<int> elementspace(length, indexBase, platformE);
    const Tpetra::SerialPlatform<int, double> platformV;
	Tpetra::VectorSpace<int, double> vectorspace(elementspace, platformV);
	Tpetra::Vector<int, double> v1(vectorspace);
    cout << "Created v1, default constructor" << endl;
    v1.printValues(cout);
    cout << "Setting all to 3.0" << endl;
    v1.setAllToScalar(3.0);
    v1.printValues(cout);
    cout << "Setting all to random" << endl;
    v1.setAllToRandom();
    v1.printValues(cout);
    cout << "Min value: " << v1.minValue() << endl;
    cout << "Max value: " << v1.maxValue() << endl;
    cout << "Mean value: " << v1.meanValue() << endl;
    cout << "1-norm: " << v1.norm1() << endl;
    cout << "2-norm: " << v1.norm2() << endl;
    cout << "Inf-norm: " << v1.normInf() << endl;
    cout << "Creating v2, setting all to random" << endl;
    Tpetra::Vector<int, double> v2(vectorspace);
    v2.setAllToRandom();
    v2.printValues(cout);
    cout << "dot product of v1.v2 : " << v1.dotProduct(v2) << endl;
    cout << "Setting v1 by 2.0 (scale):" << endl;
    cout << "v1(before): "; v1.printValues(cout);
    v1.scale(2.0);
    cout << "v1(after): "; v1.printValues(cout);
    cout << "copying v1 to user array" << endl;
    double* dblArray = new double[5];
    v1.extractCopy(dblArray);
    cout << "values of dblArray: ";
    for(int i = 0; i < length; i++)
        cout << dblArray[i] << " ";
    cout << endl;
    delete[] dblArray;
    cout << "Setting v1 to reciprocal of v2" << endl;
    v2[0] = 0.25;
    v2[1] = 5.0;
    v2[2] = 0.125;
    v2[3] = 2.0;
    v2[4] = 0.75;
    cout << "v2: "; v2.printValues(cout);
    v1.reciprocal(v2);
    cout << "v1: "; v1.printValues(cout);
    
    cout << "Elementwise multiply:" << endl;
    v1[0] = 6.25; v1[1] = 18.0; v1[2] = 0.0; v1[3] = 3.0; v1[4] = 1.0;
    cout << "v3 = v1 @ v2" << endl;
    Tpetra::Vector<int, double> v3(vectorspace);
    v3.elementwiseMultiply(1.0, v1, v2, 1.0);
    cout << "v1: "; v1.printValues(cout);
    cout << "v2: "; v2.printValues(cout);
    cout << "v3: "; v3.printValues(cout);
    cout << "v3 = v3 @ v1" << endl;
    v3.elementwiseMultiply(1.0, v1, v3, 0.0);
    cout << "v3: "; v3.printValues(cout);
    
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug) {
	int ierr = 0;
	int returnierr = 0;
  const Tpetra::SerialPlatform<OrdinalType, OrdinalType> platformE;
	const Tpetra::SerialPlatform<OrdinalType, ScalarType> platformV;
	char const * OTName = Teuchos::OrdinalTraits<OrdinalType>::name();
	char const * STName = Teuchos::ScalarTraits<ScalarType>::name();
	if(verbose) cout << "Starting unit tests for Vector<" << OTName << "," << STName << ">." << endl;

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
	if(verbose) cout << "Starting code coverage section..." << endl;

	// constructors
	if(verbose) cout << "Constructors..." << endl;
	// taking a VectorSpace
	OrdinalType const ESlength = 10;
	OrdinalType const ESindexBase = 2;
	Tpetra::ElementSpace<OrdinalType> elementspace(ESlength, ESindexBase, platformE);
	Tpetra::VectorSpace<OrdinalType, ScalarType> vectorspace(elementspace, platformV);
	Tpetra::Vector<OrdinalType, ScalarType> vector(vectorspace);
	// taking a VectorSpace and a user array of entries
	ScalarType* scalarArray = new ScalarType[ESlength];
	for(OrdinalType i = 0; i < ESlength; i++)
		scalarArray[i] = Teuchos::ScalarTraits<ScalarType>::random();
	Tpetra::Vector<OrdinalType, ScalarType> vector1a(scalarArray, ESlength, vectorspace);
	// cpy ctr
	Tpetra::Vector<OrdinalType, ScalarType> v2(vector);

	// print
	if(verbose) {
		cout << "Overloaded << operator..." << endl;
		cout << vector << endl;
	}

	// attribute access
	if(verbose) cout << "Attribute access methods..." << endl;
	OrdinalType temp = 0;
	temp = vector.getNumGlobalEntries();
	temp = vector.getNumMyEntries();

	// element access - [] operator
	if(verbose) cout << "Element access methods..." << endl;
	ScalarType const temp1 = vector[1];
	vector[0] = temp1;

	// set all to scalar
	if(verbose) cout << "setAllToScalar..." << endl;
	ScalarType scalar3 = Teuchos::ScalarTraits<ScalarType>::one();
	scalar3 += (scalar3 + scalar3); // 1 += (1+1) , should equal 3
	vector.setAllToScalar(scalar3);

	// set all to random
	if(verbose) cout << "setAllToRandom..." << endl;
	vector.setAllToRandom();

	// math functions
	
	// dot product
	if(verbose) cout << "dot product..." << endl;
	vector.dotProduct(v2); // throw away return value
  // absolute value
	if(verbose) cout << "absolute value..." << endl;
	vector.absoluteValue(v2);
  // reciprocal
	if(verbose) cout << "reciprocal..." << endl;
	vector.reciprocal(v2);
	// scale
	if(verbose) cout << "scale..." << endl;
	vector.scale(Teuchos::ScalarTraits<ScalarType>::random());
	// scale #2
	if(verbose) cout << "scale #2..." << endl;
	vector.scale(Teuchos::ScalarTraits<ScalarType>::random(), v2);
	// update
	if(verbose) cout << "update..." << endl;
	vector.update(Teuchos::ScalarTraits<ScalarType>::random(), v2, Teuchos::ScalarTraits<ScalarType>::random());
	// update #2
	if(verbose) cout << "update #2..." << endl;
	vector.update(Teuchos::ScalarTraits<ScalarType>::random(), v2, 
			      Teuchos::ScalarTraits<ScalarType>::random(), vector1a, Teuchos::ScalarTraits<ScalarType>::random());
	// 1-norm
	if(verbose) cout << "1-norm..." << endl;
	vector.norm1(); // throw away return value
	// 2-norm
	if(verbose) cout << "2-norm..." << endl;
	vector.norm2(); // throw away return value
	// Infinity-norm
	if(verbose) cout << "Infinity-norm..." << endl;
	vector.normInf(); // throw away return value
	// Weighted 2-norm
	if(verbose) cout << "Weighted 2-norm (RMS norm)..." << endl;
	vector.normWeighted(v2); // throw away return value
	// min value
	if(verbose) cout << "minValue..." << endl;
	vector.minValue();
	// max value
	if(verbose) cout << "maxValue..." << endl;
	vector.maxValue();
	// mean value
	if(verbose) cout << "meanValue..." << endl;
	vector.meanValue();
	// elementwiseMultiply
	if(verbose) cout << "elementwiseMultiply..." << endl;
	vector.elementwiseMultiply(temp1, vector1a, v2, scalar3);
	// elementwiseReciprocalMultiply
	if(verbose) cout << "elementwiseReciprocalMultiply..." << endl;
	vector.elementwiseMultiply(temp1, vector1a, v2, scalar3);

  // assignment operator
  if(verbose) cout << "assignment operator..." << endl;
  v2 = vector;

	if(verbose) cout << "Code coverage section finished." << endl;

	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

	if(verbose) cout << "Starting actual testing section..." << endl;

	// default ctr initializing to zero
	if(verbose) cout << "Checking to see that default constructor initializes to zeros... ";
	Tpetra::Vector<OrdinalType, ScalarType> testVector1(vectorspace);
	OrdinalType length = testVector1.getNumMyEntries();
	for(OrdinalType i = 0; i < length; i++)
		if(testVector1[i] != 0) {
			if(debug) cout << "element " << i << " = " << testVector1[i] << ", should be zero" << endl;
			ierr++;
		}
	if(verbose)
		if(ierr == 0) 
			cout << "Passed" << endl;
		else
			cout << "Failed" << endl;
	returnierr += ierr;
	ierr = 0;

	// user array ctr initializing correctly
	if(verbose) cout << "Checking to see that user array constructor initializes values correctly... ";
	length = vector1a.getNumMyEntries();
	for(OrdinalType i = 0; i < length; i++)
		if(vector1a[i] != scalarArray[i])
	if(verbose)
		if(ierr == 0) 
			cout << "Passed" << endl;
		else
			cout << "Failed" << endl;
	returnierr += ierr;
	ierr = 0;
	
	// changing data
	if(verbose) cout << "Changing data... ";
	ScalarType value2 = Teuchos::ScalarTraits<ScalarType>::random();
	ScalarType value5 = Teuchos::ScalarTraits<ScalarType>::random();
	ScalarType value0 = Teuchos::ScalarTraits<ScalarType>::random();
	testVector1[2] = value2;
	testVector1[5] = value5;
	testVector1[0] = value0;
	if(testVector1[2] != value2) {
		if(debug) cout << "element 2 = " << testVector1[2] << ", should be " << value2 << endl;
		ierr++;
	}
	if(testVector1[5] != value5) {
		if(debug) cout << "element 5 = " << testVector1[5] << ", should be " << value5 << endl;
		ierr++;
	}
	if(testVector1[0] != value0) {
		if(debug) cout << "element 0 = " << testVector1[0] << ", should be " << value0 << endl;
		ierr++;
	}
	if(verbose)
		if(ierr == 0) 
			cout << "Passed" << endl;
		else
			cout << "Failed" << endl;
	returnierr += ierr;
	ierr = 0;

	// check updates
	Tpetra::Vector<OrdinalType, ScalarType> u1(vectorspace);
	Tpetra::Vector<OrdinalType, ScalarType> u2(vectorspace);
	Tpetra::Vector<OrdinalType, ScalarType> u3(vectorspace);

	u1[0] = -1; u1[1] = 2; u1[2] = -1;
	u1[3] = -1; u1[4] = 2; u1[5] = -1;
	u1[6] = -1; u1[7] = 2; u1[8] = -1;
	u1[9] = 2;

	u2[0] = 2; u2[1] = 3; u2[2] = 4; u2[3] = 5;
	u2[4] = 4; u2[5] = 3; u2[6] = 2; u2[7] = 1;
	u2[8] = 8; u2[9] = 9;

	cout << "before update:" << endl;
	cout << "u1:" << endl << u1 << endl;
	cout << "u2:" << endl << u2 << endl;
	cout << "u3:" << endl << u3 << endl << endl;

	u3.update(1.0, u1, 2.0, u2, 0.0);

	cout << "after update:" << endl;
	cout << "u1:" << endl << u1 << endl;
	cout << "u2:" << endl << u2 << endl;
	cout << "u3:" << endl << u3 << endl << endl;

	// finish up
	if(verbose)
		if(returnierr == 0)
			cout << "VectorTest <" << OTName << ", " << STName << "> passed." << endl;
		else
			cout << "VectorTest <" << OTName << ", " << STName << "> failed." << endl;
	return(returnierr);
}
