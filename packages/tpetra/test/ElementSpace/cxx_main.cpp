// Tpetra ElementSpace tester
// Modified: 14-Oct-2002
// Paul says: I'll get around to making this a more thorough tester real soon now.

#define ORDINALTYPE int

#include <iostream>
#include <iomanip>
#include "Tpetra_SerialPlatform.h" 
#include "Tpetra_ElementSpace.h"

void platformTester(bool verbose, bool debug);
void esTester(bool verbose, bool debug);
void isLgetG(ORDINALTYPE low, ORDINALTYPE high, Tpetra::ElementSpace<ORDINALTYPE>& es);
void isGgetL(ORDINALTYPE low, ORDINALTYPE high, Tpetra::ElementSpace<ORDINALTYPE>& es);

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

	platformTester(verbose, debug);
	esTester(verbose, debug);

	return(0);
}

void platformTester(bool verbose, bool debug) {
	if(verbose) cout << "Creating platform...";
  Tpetra::SerialPlatform<ORDINALTYPE, ORDINALTYPE> platform;
	if(verbose) cout << "Successful." << endl;
	if(debug) cout << platform << endl;
}

void esTester(bool verbose, bool debug) {
	ORDINALTYPE eList[10] = {1,4,7,8,9,15,22,54,55,58};
	ORDINALTYPE eSize = 10;
  Tpetra::SerialPlatform<ORDINALTYPE, ORDINALTYPE> platform;
	
	if(verbose) cout << "Creating es1(contiguous, tpetra-defined)...";
	Tpetra::ElementSpace<ORDINALTYPE> es1(10, 2, platform);
	if(verbose) cout << "Successful." << endl;
	if(debug) cout << es1 << endl;
	
	if(verbose) cout << "Creating es2(contiguous, user-defined)...";
	Tpetra::ElementSpace<ORDINALTYPE> es2(-1, 10, 2, platform);
	if(verbose) cout << "Successful." << endl;
	if(debug) cout << es2 << endl;
	
	if(verbose) cout << "Testing isSameAs (contig)...";
	assert(es1.isSameAs(es2) == true);
	Tpetra::ElementSpace<ORDINALTYPE> es2a(10, 10, 3, platform);
	assert(es1.isSameAs(es2a) == false);
	if(verbose) cout << "Successful." << endl;
	
	if(verbose) cout << "Creating es3(noncontiguous)..."; 
	Tpetra::ElementSpace<ORDINALTYPE> es3(-1, eSize, eList, 0, platform);
	if(verbose) cout << "Successful." << endl;
	if(debug) cout << es3 << endl;
	
	if(debug) {
		cout << "Testing isMyLID and getGID" << endl;
		isLgetG(0, 13, es3); 
		cout << "Testing isMyGID and getLID" << endl; 
		isGgetL(0, 60, es3);
	}
	
	if(verbose) cout << "Testing isSameAs (noncontig)...";
	Tpetra::ElementSpace<ORDINALTYPE> es3a(eSize, eSize, eList, es3.getIndexBase(), platform);
	assert(es3.isSameAs(es3a) == true);
	eList[(eSize / 2)] += 2;
	Tpetra::ElementSpace<ORDINALTYPE> es3b(eSize, eSize, eList, es3.getIndexBase(), platform);
	assert(es3.isSameAs(es3b) == false);
	if(verbose) cout << "Successful." << endl;

	if(verbose) cout << "Testing copy constructor...";
	Tpetra::ElementSpace<ORDINALTYPE> es4(es3);
	assert(es3.isSameAs(es4) == true);
	assert(es4.isSameAs(es3b) == false);
	if(verbose) cout << "Successful." << endl;

	if(verbose) cout << "Testing getRemoteIDList...";
	const int len = 5;
	ORDINALTYPE gList[len] = {1,4,22,55,58};
	ORDINALTYPE pList[len] = {5,5,5,5,5};
	ORDINALTYPE lList[len] = {0,0,0,0,0};
	es3.getRemoteIDList(5, gList, pList, lList);
	if(debug) cout << "\nGID PID LID getLID" << endl;
	for(int i = 0; i < len; i++) {
		if(debug) cout << setw(3) << gList[i] << setw(4) << pList[i] << setw(4) << lList[i] << setw(4) << es3.getLID(gList[i]) << endl;
		assert(lList[i] == es3.getLID(gList[i]));
	}
	if(verbose) cout << "Successful." << endl;
	
	cout << "ElementSpace test successful." << endl;
}

void isLgetG(ORDINALTYPE low, ORDINALTYPE high, Tpetra::ElementSpace<ORDINALTYPE>& es) {
	  for(ORDINALTYPE i = low; i < high; i++) {
		  if(es.isMyLID(i))
			  cout << "LID" << setw(3) << i << " getGID? " << es.getGID(i) << endl;
	  }
}

void isGgetL(ORDINALTYPE low, ORDINALTYPE high, Tpetra::ElementSpace<ORDINALTYPE>& es) {
	  for(ORDINALTYPE i = low; i < high; i++) {
		  if(es.isMyGID(i))
			  cout << "GID" << setw(3) << i << " getLID? " << es.getLID(i) << endl;
	  }
}
