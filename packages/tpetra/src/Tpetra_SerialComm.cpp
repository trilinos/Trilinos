/*Paul
27-July-2002 Status: Templated, OrdinalType only. All Epetra methods except Distributor.
*/

#include "Tpetra_SerialDirectory.h"

namespace Tpetra {

	// default constructor
	template<typename OrdinalType>
	SerialComm<OrdinalType>::SerialComm() 
		: Object("Tpetra::Comm[Serial]") {}

	// copy constructor
	template<typename OrdinalType>
	SerialComm<OrdinalType>::SerialComm(const SerialComm<OrdinalType>& Comm) 
		: Object(Comm.label()) {}

	// destructor
	template<typename OrdinalType>
	SerialComm<OrdinalType>::~SerialComm() {}

	// gather all
	template<typename OrdinalType>
  void SerialComm<OrdinalType>::gatherAll(OrdinalType* myVals, OrdinalType* allVals, OrdinalType count) const {
		for(OrdinalType i=0; i<count; i++)
			allVals[i] = myVals[i];
	}

	// sum
	template<typename OrdinalType>
	void SerialComm<OrdinalType>::sumAll(OrdinalType* partialSums, OrdinalType* globalSums, OrdinalType count) const {
		for(OrdinalType i=0; i<count; i++)
			globalSums[i] = partialSums[i];
	}
	
	// min/max
	template<typename OrdinalType>
	void SerialComm<OrdinalType>::maxAll(OrdinalType* partialMaxs, OrdinalType* globalMaxs, OrdinalType count) const {
		for(OrdinalType i=0; i<count; i++)
			globalMaxs[i] = partialMaxs[i];
	}
	template<typename OrdinalType>
	void SerialComm<OrdinalType>::minAll(OrdinalType* partialMins, OrdinalType* globalMins, OrdinalType count) const {
		for(OrdinalType i=0; i<count; i++)
			globalMins[i] = partialMins[i];
	}
	
	// scansum
	template<typename OrdinalType>
  void SerialComm<OrdinalType>::scanSum(OrdinalType* myVals, OrdinalType* scanSums, OrdinalType count) const {
		for(OrdinalType i=0; i<count; i++)
			scanSums[i] = myVals[i];
  }

	// directory
	template<typename OrdinalType>
	Directory<OrdinalType>* SerialComm<OrdinalType>::createDirectory(const ElementSpace<OrdinalType>& ElementSpace) const {
		Directory<OrdinalType>* dir = dynamic_cast<Directory<OrdinalType>*>(new SerialDirectory<OrdinalType>(ElementSpace));
		return(dir);
	}

	// print
	template<typename OrdinalType>
	void SerialComm<OrdinalType>::print(ostream& os) const {
		os << "::Processor " << myPID() << " of " << numProc() << " total processors";
	}
	
} // namespace Tpetra
