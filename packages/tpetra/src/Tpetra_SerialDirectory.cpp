/*Paul
27-July-2002 Templated for OrdinalType.
06-August-2002 Changed to images.
*/

#include "Tpetra_ElementSpace.h"

namespace Tpetra {

// default constructor
template<class OrdinalType>
SerialDirectory<OrdinalType>::SerialDirectory(const ElementSpace<OrdinalType>& ElementSpace) 
	: Object("Tpetra::Directory[Serial]") 
	, ElementSpace_(&ElementSpace) {}

// copy constructor
template<class OrdinalType>
SerialDirectory<OrdinalType>::SerialDirectory(const SerialDirectory<OrdinalType>& Directory) 
	: Object(Directory.label()) 
	, ElementSpace_(Directory.ElementSpace_) {}

// destructor
template<class OrdinalType>
SerialDirectory<OrdinalType>::~SerialDirectory() {}

// query method
template<class OrdinalType>
void SerialDirectory<OrdinalType>::getDirectoryEntries(OrdinalType numEntries, const OrdinalType* globalEntries, OrdinalType* images, OrdinalType* localEntries) const {
	int imageID = ElementSpace_->comm().getMyImageID();
	for(OrdinalType i = 0; i < numEntries; i++) {
		if(!ElementSpace_->isMyGID(globalEntries[i]))
			throw reportError("Global ID " + toString(globalEntries[i]) + " was not found in this ElementSpace.", 1);
		else {
			images[i] = imageID;
			localEntries[i] = ElementSpace_->getLID(globalEntries[i]);
		}
	}
}
	
} // namespace Tpetra
