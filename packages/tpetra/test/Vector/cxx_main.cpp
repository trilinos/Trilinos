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

	// finish up
	if(ierr == 0)
		cout << "Vector test successful." << endl;
	else
		cout << "Vector test failed." << endl;
	return(ierr);
}

template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug) {
	int ierr = 0;
	int returnierr = 0;
  const Tpetra::SerialPlatform<OrdinalType, OrdinalType> platformE;
	const Tpetra::SerialPlatform<OrdinalType, ScalarType> platformV;
	char const * OTName = Teuchos::OrdinalTraits<OrdinalType>::name();
	char const * STName = Teuchos::ScalarTraits<ScalarType>::name();
	if(verbose) cout << "Starting unit tests for Vector<" << OTName << "," << STName << ">." << endl;

	//
	// code coverage section - just call functions, no testing
	//
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

	// element access
	if(verbose) cout << "Element access methods..." << endl;
	ScalarType const temp1 = vector[1];
	vector[0] = temp1;

	// set all to scalar
	if(verbose) cout << "setAllToScalar..." << endl;
	vector.setAllToScalar(3.0);

	// set all to random
	if(verbose) cout << "setAllToRandom..." << endl;
	vector.setAllToRandom();

	if(verbose) cout << "Code coverage section finished." << endl;

	//
	// actual testing section - affects return code
	//

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

	// finish up
	if(verbose)
		if(returnierr == 0)
			cout << "VectorTest <" << OTName << ", " << STName << "> passed." << endl;
		else
			cout << "VectorTest <" << OTName << ", " << STName << "> failed." << endl;
	return(returnierr);
}
