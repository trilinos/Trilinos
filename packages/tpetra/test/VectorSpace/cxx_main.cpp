// Tpetra VectorSpace tester
// Modified: 06-Feb-2003

//#define ORDINALTYPE int
//#define SCALARTYPE float

#include <iostream>
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_Vector.hpp"

int main(int argc, char* argv[]) {
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

	int returnierr = 0;
	int ierr = 0;

  const Tpetra::SerialPlatform<int, int> platformE;
	const Tpetra::SerialPlatform<int, float> platformV;

	//
	// code coverage section - just call functions, no testing
	//

	// ctr and cpy ctr
	const Tpetra::ElementSpace<int> elementspace(10, 2, platformE);
	Tpetra::VectorSpace<int, float> vectorspace(elementspace, platformV);
	Tpetra::VectorSpace<int, float> v2(vectorspace);

	// print
	cout << vectorspace << endl;

	// attribute access
	int temp = 0;
	temp = vectorspace.getNumGlobalEntries();
	temp = vectorspace.getNumMyEntries();
	temp = vectorspace.getIndexBase();
	temp = vectorspace.getMinLocalIndex();
	temp = vectorspace.getMaxLocalIndex();
	temp = vectorspace.getMinGlobalIndex();
	temp = vectorspace.getMaxGlobalIndex();
	temp = 0;
	temp = vectorspace.getGlobalIndex(temp);
	temp = vectorspace.getLocalIndex(temp);

	// accessors to other classes
	v2.platform().printInfo(cout);
	temp = v2.comm().getNumImages();

	//
	// actual testing section - affects return code
	//

	// compatibleVector
	Tpetra::Vector<int, float> vector(vectorspace);
	if(!vectorspace.compatibleVector(vector))
		ierr++;
	const Tpetra::ElementSpace<int> differentES(15, 2, platformE);
	Tpetra::VectorSpace<int, float> differentVS(differentES, platformV);
	if(differentVS.compatibleVector(vector))
		ierr++;
	if(ierr == 0) 
		cout << "Compatibility test passed" << endl;
	else
		cout << "Compatibility test failed" << endl;
	returnierr += ierr;
	ierr = 0;

	// vector creation
	Tpetra::Vector<int, float>* vecptr = vectorspace.createVector();
	temp = vecptr->getNumMyEntries();
	delete vecptr;

	// isSameAs
	bool same = vectorspace.isSameAs(v2);
	if(!same)
		ierr++;
	same = vectorspace.isSameAs(differentVS);
	if(same)
		ierr++;
	if(ierr == 0) 
		cout << "IsSameAs test passed" << endl;
	else
		cout << "IsSameAs test failed" << endl;
	returnierr += ierr;
	ierr = 0;

	// finish up
	if(returnierr == 0)
		cout << "VectorSpaceTest passed." << endl;
	else
		cout << "VectorSpaceTest failed." << endl;
	return(returnierr);
}
