/*Paul
03-August-2002 BES tester. Initial writeup.
*/

#define ORDINALTYPE int

#include <iostream>
#include <iomanip>
#include "Tpetra_SerialComm.h" 
#include "Tpetra_ElementSpace.h"
#include "Tpetra_BlockElementSpace.h"

int main(int argc, char* argv[]) {

	const ORDINALTYPE INDEXBASE = 0;
	const ORDINALTYPE NUMELEMENTS = 5;
	const ORDINALTYPE ELEMENTSIZE = 2;
  
  bool verbose = false;
  if((argc > 1) && (argv[1][0] == '-') && (argv[1][1] == 'v'))
      verbose = true;
  
  // Comm
  if(verbose) cout << "==Creating comm" << endl;
  Tpetra::SerialComm<ORDINALTYPE> comm;

  // ElementSpace
  if(verbose) cout << "==Creating es, constructor1" << endl;
  Tpetra::ElementSpace<ORDINALTYPE> es(5, 0, comm);
  //if(verbose) cout << "==Creating es, constructor2" << endl;
  //Tpetra::ElementSpace<ORDINALTYPE> es(-1, NUMELEMENTS, INDEXBASE, comm);
  //if(verbose) cout << "==Creating es, constructor3" << endl;
  //ORDINALTYPE gidList[NUMELEMENTS] = {1,4,7,8,9};//,15,22,54,55,58};
  //Tpetra::ElementSpace<ORDINALTYPE> es(-1, NUMELEMENTS, gidList, INDEXBASE, comm);
  
  //cout << es;

  //BlockElementSpace
  if(verbose) cout << "==Creating bes, constructor1" << endl;
  Tpetra::BlockElementSpace<ORDINALTYPE> bes(es, 2);
	//ORDINALTYPE* tmp = bes.elementSpace().getMyGlobalElements();
	//cout << "past gMGE" << endl;
  cout << bes;
	//cout << bes.elementSpace();

  return(0); 
}

