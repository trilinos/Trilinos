// Tpetra Distributor tester
// Modified: 23-Nov-2002

#define PACKETTYPE float
#define ORDINALTYPE int


#include "Tpetra_SerialDistributor.h"

int main(int argc, char* argv[]) {
	bool verbose = false;
	if (argc>1 && argv[1][0]=='-' && argv[1][1]=='v') 
		verbose = true;

  if(verbose) cout << "Creating SerialDistributor object...";
  Tpetra::SerialDistributor<PACKETTYPE, ORDINALTYPE> distributor;
  //if(debug) cout <<distributor.label() << endl;
	if(verbose) cout << "Successful." << endl;

	//  void createFromSends(const OrdinalType& numExportIDs, const OrdinalType* exportImageIDs,
	//											 const bool& deterministic, OrdinalType& numRemoteIDs ) 

	ORDINALTYPE nEIDs = 2;
	ORDINALTYPE* eIIDs = 0;
	bool determ = false;
	ORDINALTYPE nRIDs = 2;
	//distributor.createFromSends(nEIDs, eIIDs, false, nRIDs);
  
	cout << "Distributor test successful." << endl;
  
  return(0);
}
