/*Paul
16-July-2002 CommTester.
*/

#define ORDINALTYPE int
#define SCALARTYPE float

#include "Tpetra_SerialComm.h"

template<class T>
void setArray(T* array, T v1, T v2) {
	array[0] = v1;
	array[1] = v2;
}

int main(int argc, char* argv[]) {

  ORDINALTYPE in[2];
  ORDINALTYPE out[2];
  
  Tpetra::SerialComm<ORDINALTYPE> comm;
  cout << "===SerialComm object created." << endl;
  cout << comm.label() << endl;
  
  cout << "===Testing getImageID and getNumImages." << endl;
  assert(comm.getMyImageID() == 0);
  assert(comm.getNumImages() == 1);
  
  cout << "===Testing barrier." << endl;
  comm.barrier();
  
  cout << "===Testing SumAll(ordinal)." << endl;
  setArray(in, 2, 4);
  setArray(out, 0, 0);
  comm.sumAll(in, out, 2);
  assert(out[0] == 2);
  assert(out[1] == 4);
  
  /*
    cout << "===Testing maxAll." << endl;
    setArray(in, 3.0, 9.8);
    setArray(out, 0.0, 0.0);
    comm.maxAll(in, out, 2);
    assert(in[0] == SCALARTYPE(3.0));
    assert(out[1] == SCALARTYPE(9.8));
  */
  
  cout << "===Finished." << endl;
  
  return(0);
}
