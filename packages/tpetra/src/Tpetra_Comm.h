/*Paul
27-July-2002 Status: Templated, OrdinalType only. All Epetra methods except Distributor.
*/

#ifndef _TPETRA_COMM_H_
#define _TPETRA_COMM_H_

#include "Tpetra_Common.h"

namespace Tpetra {

// forward declarations
template<typename OrdinalType> class Directory;
template<typename OrdinalType> class ElementSpace;

template<typename OrdinalType>
class Comm {
  public:

  // Destructor.
  virtual ~Comm() {};

  // myPID
  virtual int myPID() const = 0;
  
  // numProc
  virtual int numProc() const = 0;
	
	// Barrier
  virtual void barrier() const = 0;

	// Broadcast
  virtual void broadcast(OrdinalType* MyVals, OrdinalType count, int root) const = 0;

	// Gather All
  virtual void gatherAll(OrdinalType* myVals, OrdinalType* allVals, OrdinalType count) const = 0;

  // Sum
  virtual void sumAll(OrdinalType* partialSums, OrdinalType* globalSums, OrdinalType count) const = 0;
	
  // Max/Min
  virtual void maxAll(OrdinalType* partialMaxs, OrdinalType* globalMaxs, OrdinalType count) const = 0;
  virtual void minAll(OrdinalType* partialMins, OrdinalType* globalMins, OrdinalType count) const = 0;

  // Scan Sum
  virtual void scanSum(OrdinalType* myVals, OrdinalType* scanSums, OrdinalType count) const = 0;

	// Directory
	virtual Directory<OrdinalType>* createDirectory(const ElementSpace<OrdinalType>& ElementSpace) const = 0;

	// Print Info
	virtual void printInfo(ostream& os) const = 0;
	
}; // class Comm

} // namespace Tpetra

#endif // _TPETRA_COMM_H_
