#define ORDINALTYPE int

#include <iostream>
#include <iomanip>
#include "Tpetra_SerialComm.h" 
#include "Tpetra_ElementSpace.h"

void commTester(bool verbose);
void esTester(bool verbose);
void isLgetG(ORDINALTYPE low, ORDINALTYPE high, Tpetra::ElementSpace<ORDINALTYPE>& es);
void isGgetL(ORDINALTYPE low, ORDINALTYPE high, Tpetra::ElementSpace<ORDINALTYPE>& es);

int main(int argc, char* argv[]) {

	bool verbose = false;
	if(argc > 1)
		if(argv[1][0] == '-' && argv[1][1] == 'v')
			verbose = true;

	commTester(verbose);
	esTester(verbose);

	return(0);
}

void commTester(bool verbose) {
	cout << "==Creating comm" << endl;
  Tpetra::SerialComm<ORDINALTYPE> comm;
	if(verbose) cout << "==comm.print(cout)" << endl;
	if(verbose) comm.print(cout);

}

void esTester(bool verbose) {
	ORDINALTYPE eList[10] = {1,4,7,8,9,15,22,54,55,58};
	ORDINALTYPE eSize = 10;
  Tpetra::SerialComm<ORDINALTYPE> comm;
	
	// Create ES, const 1
	cout << "==Creating es1" << endl;
	Tpetra::ElementSpace<ORDINALTYPE> es1(10, 2, comm);
	if(verbose) cout << "==es1.print(cout)" << endl;
	if(verbose) es1.print(cout);
	
	cout << "==Creating es2" << endl;
	Tpetra::ElementSpace<ORDINALTYPE> es2(-1, 10, 2, comm);
	if(verbose) cout << "==es2.print(cout)" << endl;
	if(verbose) es2.print(cout);
	
	cout << "==Testing isSameAs (contig)" << endl;
	assert(es1.isSameAs(es2) == true);
	Tpetra::ElementSpace<ORDINALTYPE> es2a(10, 10, 3, comm);
	assert(es1.isSameAs(es2a) == false);
	
	cout << "==Creating es3" << endl; 
	Tpetra::ElementSpace<ORDINALTYPE> es3(-1, eSize, eList, 0, comm);
	if(verbose) cout << "==es3.print(cout)" << endl;
	if(verbose) es3.print(cout);
	
	if(verbose) cout << "==Testing isMyLID and getGID" << endl;
	if(verbose) isLgetG(0, 13, es3); 
	
	if(verbose) cout << "==Testing isMyGID and getLID" << endl; 
	if(verbose) isGgetL(0, 60, es3);
	
	cout << "==Testing isSameAs (noncontig)" << endl;
	Tpetra::ElementSpace<ORDINALTYPE> es3a(eSize, eSize, eList, es3.getIndexBase(), comm);
	assert(es3.isSameAs(es3a) == true);
	eList[(eSize / 2)] += 2;
	Tpetra::ElementSpace<ORDINALTYPE> es3b(eSize, eSize, eList, es3.getIndexBase(), comm);
	assert(es3.isSameAs(es3b) == false);

	cout << "==Testing constructor #4" << endl;
	Tpetra::ElementSpace<ORDINALTYPE> es4(es3);
	assert(es3.isSameAs(es4) == true);
	assert(es4.isSameAs(es3b) == false);

	cout << "==Testing getRemoteIDList" << endl;
	ORDINALTYPE gList[5] = {1,4,22,55,58};
	ORDINALTYPE pList[5] = {5,5,5,5,5};
	ORDINALTYPE lList[5] = {0,0,0,0,0};
	es3.getRemoteIDList(5, gList, pList, lList);
	if(verbose)
		for(int i = 0; i < 5; i++)
			cout << setw(3) << gList[i] << setw(3) << pList[i] << setw(3) << lList[i] << endl;
	
	
	cout << "==ESTester Finished." << endl;
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

