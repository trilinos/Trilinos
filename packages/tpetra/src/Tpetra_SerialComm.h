/*Paul
27-July-2002 Status: Templated, OrdinalType only. All Epetra methods except Distributor.
*/

#ifndef _TPETRA_SERIALCOMM_H_
#define _TPETRA_SERIALCOMM_H_

#include "Tpetra_Comm.h"
#include "Tpetra_Object.h"

namespace Tpetra {

// forward declarations
template<typename OrdinalType> class ElementSpace;
template<typename OrdinalType> class Directory;

template<typename OrdinalType>
class SerialComm : public Object, public virtual Comm<OrdinalType> {
public:
  // constructor
	SerialComm();
  
	// copy constructor
  SerialComm(const SerialComm<OrdinalType>& Comm);
  
  // destructor.
	~SerialComm();
	
	// barrier
	void barrier() const {};

	// broadcast
  void broadcast(OrdinalType* MyVals, OrdinalType count, int root) const {};

	// gather all
  void gatherAll(OrdinalType* myVals, OrdinalType* allVals, OrdinalType count) const;

	// sum
	void sumAll(OrdinalType* partialSums, OrdinalType* globalSums, OrdinalType count) const;
	
	// max/min
  void maxAll(OrdinalType* partialMaxs, OrdinalType* globalMaxs, OrdinalType count) const;
  void minAll(OrdinalType* partialMins, OrdinalType* globalMins, OrdinalType count) const;

	// scansum
  void scanSum(OrdinalType* myVals, OrdinalType* scanSums, OrdinalType count) const;
	
  // myPID
  int myPID() const {return(0);};
  	
  // numProcs
  int numProc() const {return(1);};

	// directory
	Directory<OrdinalType>* createDirectory(const ElementSpace<OrdinalType>& ElementSpace) const;

	// print
	void print(ostream& os) const;
	void printInfo(ostream& os) const {print(os);};

}; // class SerialComm

} // namespace Tpetra

#include "Tpetra_SerialComm.cpp"
#endif // _TPETRA_SERIALCOMM_H_
