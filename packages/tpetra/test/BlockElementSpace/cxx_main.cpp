/*Paul
03-August-2002 BES tester. Initial writeup.
18-Oct-2002 Modified.
22-Oct-2002 Changed to test out BES/BESData friend wrinkles.
*/

#define ORDINALTYPE int

#include <iostream>
#include <iomanip>
#include "Tpetra_SerialPlatform.h" 
#include "Tpetra_ElementSpace.h"
#include "Tpetra_BlockElementSpace.h"

int main(int argc, char* argv[]) {

	const ORDINALTYPE INDEXBASE = 0;
	const ORDINALTYPE NUMELEMENTS = 5;
	const ORDINALTYPE ELEMENTSIZE = 2;
  
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
  
  // Platform
  if(verbose) cout << "Creating platform" << endl;
  Tpetra::SerialPlatform<ORDINALTYPE, ORDINALTYPE> platform;

  // ElementSpace
	// commented out lines are for creating alternate es objects.
  if(verbose) cout << "Creating es, constructor1" << endl;
  Tpetra::ElementSpace<ORDINALTYPE> es(NUMELEMENTS, INDEXBASE, platform);
  //if(verbose) cout << "Creating es, constructor2" << endl;
  //Tpetra::ElementSpace<ORDINALTYPE> es(-1, NUMELEMENTS, INDEXBASE, platform);
  //if(verbose) cout << "Creating es, constructor3" << endl;
  //ORDINALTYPE gidList[NUMELEMENTS] = {1,4,7,8,9};//,15,22,54,55,58};
  //Tpetra::ElementSpace<ORDINALTYPE> es(-1, NUMELEMENTS, gidList, INDEXBASE, platform);
  //if(debug) cout << es;

  //BlockElementSpace
  if(verbose) cout << "Creating bes, constructor1" << endl;
  Tpetra::BlockElementSpace<ORDINALTYPE> bes(es, ELEMENTSIZE);
  if(debug) cout << bes;

	if(verbose) cout << "Creating bes, constructor2" << endl;
	ORDINALTYPE* esizelist = new ORDINALTYPE[NUMELEMENTS];
	esizelist[0] = 1;
	esizelist[1] = 1;
	esizelist[2] = 2;
	esizelist[3] = 3;
	esizelist[4] = 5;
	Tpetra::BlockElementSpace<ORDINALTYPE> bes2(es, esizelist);
	delete[] esizelist;
	esizelist = 0;
	if(debug) cout << bes2;
	
	cout << "BlockElementSpace testing successful." << endl;
  return(0); 
}

