/*Paul
16-July-2002 CommTester.
21-Sept-2002 Updated for Comm/Platform split.
12-Nov-2002 Updated for new templating scheme (no changes).
*/

#include <iostream>
#include "Tpetra_SerialComm.h"

// function prototypes
template<typename PacketType, typename OrdinalType> void setRandom(PacketType& vals, OrdinalType count);
template<typename PacketType, typename OrdinalType> void setToZero(PacketType& vals, OrdinalType count);
template<typename PacketType, typename OrdinalType> void commTest(Tpetra::SerialComm<PacketType, OrdinalType>& comm, bool verbose);

int main(int argc, char* argv[]) {
	bool verbose = false;
	if (argc>1 && argv[1][0]=='-' && argv[1][1]=='v') 
		verbose = true;

	if(verbose) cout << "Creating SerialComm object...";
	Tpetra::SerialComm<int, int> comm;
  if(verbose) cout << "Successful" << endl;
	commTest(comm, verbose);
	cout << "SerialComm testing successfull." << endl;

  return(0);
}

template<typename PacketType, typename OrdinalType> 
void setValues(PacketType& vals, OrdinalType count) {
	for(OrdinalType i = 0; i < count; i++)
		vals[i] = 5;
}

template<typename PacketType, typename OrdinalType> 
void setToZero(PacketType& vals, OrdinalType count) {
	for(OrdinalType i = 0; i < count; i++)
		vals[i] = 0;
}

template<typename PacketType, typename OrdinalType>
void commTest(Tpetra::SerialComm<PacketType, OrdinalType>& comm, bool verbose) {
	OrdinalType count;
	count = 1;	// eventually, count = randnum
	PacketType inVal;
	PacketType outVal;

	inVal = 7; outVal = 7;
	if(verbose) cout << "Testing broadcast...";
	comm.broadcast(&inVal, count, 0);
	assert(inVal = outVal);
	if(verbose) cout << "Successful" << endl;
  
	inVal = 5; outVal = 0;
	if(verbose) cout << "Testing gatherAll...";
	comm.gatherAll(&inVal, &outVal, count);
	assert(inVal = outVal);
	if(verbose) cout << "Successful" << endl;

	inVal = 2; outVal = 0;
	if(verbose) cout << "Testing sumAll...";
	comm.sumAll(&inVal, &outVal, count);
	assert(inVal = outVal);
	if(verbose) cout << "Successful" << endl;

	inVal = 6; outVal = 0;
	if(verbose) cout << "Testing MaxAll...";
	comm.maxAll(&inVal, &outVal, count);
	assert(inVal = outVal);
	if(verbose) cout << "Successful" << endl;

	inVal = 9; outVal = 0;
	if(verbose) cout << "Testing MinAll...";
	comm.minAll(&inVal, &outVal, count);
	assert(inVal = outVal);
	if(verbose) cout << "Successful" << endl;

	inVal = 1; outVal = 0;
	if(verbose) cout << "Testing scanSum...";
	comm.scanSum(&inVal, &outVal, count);
	assert(inVal = outVal);
	if(verbose) cout << "Successful" << endl;

}
