/*Paul
16-July-2002 CommTester.
*/

#define ORDINALTYPE int
#define SCALARTYPE float

#include "Tpetra_SerialComm.h"

template<typename T>
void setArray(T* array, T v1, T v2) {
	array[0] = v1;
	array[1] = v2;
}

int main(int argc, char* argv[]) {

	ORDINALTYPE ordinalIn[2];
	ORDINALTYPE ordinalOut[2];
	SCALARTYPE scalarIn[2];
	SCALARTYPE scalarOut[2];

  Tpetra::SerialComm<ORDINALTYPE, SCALARTYPE> comm;
	cout << "===SerialComm object created." << endl;
	cout << comm.label() << endl;

	cout << "===Testing myPID and NumProc." << endl;
	assert(comm.myPID() == 0);
	assert(comm.numProc() == 1);

	cout << "===Testing barrier." << endl;
	comm.barrier();

	cout << "===Testing SumAll(ordinal)." << endl;
	setArray(ordinalIn, 2, 4);
	setArray(ordinalOut, 0, 0);
	comm.sumAll(ordinalIn, ordinalOut, 2);
	assert(ordinalOut[0] == 2);
	assert(ordinalOut[1] == 4);

	/*
	cout << "===Testing maxAll." << endl;
	setArray(scalarIn, 3.0, 9.8);
	setArray(scalarOut, 0.0, 0.0);
	comm.maxAll(scalarIn, scalarOut, 2);
	assert(scalarIn[0] == SCALARTYPE(3.0));
	assert(scalarOut[1] == SCALARTYPE(9.8));
	*/

	cout << "===Finished." << endl;
	
	return 0;
}
