// Tpetra Vector tester
// Modified: 06-Feb-2003

#include <iostream>
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
		cout << "Vector test successfull." << endl;
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
	Tpetra::ElementSpace<OrdinalType> elementspace(10, 2, platformE);
	Tpetra::VectorSpace<OrdinalType, ScalarType> vectorspace(elementspace, platformV);
	Tpetra::Vector<OrdinalType, ScalarType> vector(vectorspace);
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

	if(verbose) cout << "Code coverage section finished." << endl;
	
	// finish up
	if(verbose)
		if(returnierr == 0)
			cout << "VectorTest <" << OTName << ", " << STName << "> passed." << endl;
		else
			cout << "VectorTest <" << OTName << ", " << STName << ">failed." << endl;
	return(returnierr);
}
