// Tpetra VectorSpace tester
// Modified: 06-Feb-2003

//#define ORDINALTYPE int
//#define SCALARTYPE float

#include <iostream>
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_VectorSpace.hpp" 

//Tpetra::ElementSpace<ORDINALTYPE>* esCreator(int constructorNum, bool verbose, bool debug);

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

  const Tpetra::SerialPlatform<int, int> platformE;
	const Tpetra::ElementSpace<int> elementspace(10, 2, platformE);

  const Tpetra::SerialPlatform<int, float> platformV;
	Tpetra::VectorSpace<int, float> vectorspace(elementspace, platformV);

	return(0);
}

/*Tpetra::ElementSpace<ORDINALTYPE>* esCreator(int constructorNum, bool verbose, bool debug) {
  Tpetra::SerialPlatform<ORDINALTYPE, ORDINALTYPE> platform; 
	ORDINALTYPE eList[10] = {1,4,7,8,9,15,22,54,55,58}; 
	ORDINALTYPE eSize = 10; 
	ORDINALTYPE iB = 0; 
	Tpetra::ElementSpace<ORDINALTYPE>* es = 0; 
	
	switch(constructorNum) { 
	case 1: 
		if(verbose) cout << "Creating es1(contiguous, tpetra-defined)..."; 
		es = new Tpetra::ElementSpace<ORDINALTYPE>(eSize, iB, platform); 
		break; 
	case 2: 
		if(verbose) cout << "Creating es2(contiguous, user-defined)..."; 
		es = new Tpetra::ElementSpace<ORDINALTYPE>(-1, eSize, iB, platform); 
		break; 
	case 3: 
		if(verbose) cout << "Creating es3(noncontiguous)..."; 
		es = new Tpetra::ElementSpace<ORDINALTYPE>(-1, eSize, eList, iB, platform); 
		break; 
	}; 

	if(verbose) cout << "Successful." << endl; 
	if(debug) cout << es << endl; 

	return(es); 
	}*/
