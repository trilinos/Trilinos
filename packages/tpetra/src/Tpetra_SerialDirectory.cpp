/*Paul
27-July-2002 Templated for OrdinalType.
*/

#include "Tpetra_ElementSpace.h"

namespace Tpetra {

// default constructor
template<typename OrdinalType>
SerialDirectory<OrdinalType>::SerialDirectory(const ElementSpace<OrdinalType>& ElementSpace) 
	: Object("Tpetra::Directory[Serial]") 
	, ElementSpace_(&ElementSpace) {}

// copy constructor
template<typename OrdinalType>
SerialDirectory<OrdinalType>::SerialDirectory(const SerialDirectory<OrdinalType>& Directory) 
	: Object(Directory.label()) 
	, ElementSpace_(Directory.ElementSpace_) {}

// destructor
template<typename OrdinalType>
SerialDirectory<OrdinalType>::~SerialDirectory() {}

// query method
template<typename OrdinalType>
void SerialDirectory<OrdinalType>::getDirectoryEntries(OrdinalType numEntries, const OrdinalType* globalEntries, OrdinalType* procs, OrdinalType* localEntries) const {
	int myPID = ElementSpace_->comm().myPID();
	for(OrdinalType i = 0; i < numEntries; i++) {
		if(!ElementSpace_->isMyGID(globalEntries[i]))
			throw reportError("Global ID " + toString(globalEntries[i]) + " was not found in this ElementSpace.", 1);
		else {
			procs[i] = myPID;
			localEntries[i] = ElementSpace_->getLID(globalEntries[i]);
		}
	}
}
	
} // namespace Tpetra
