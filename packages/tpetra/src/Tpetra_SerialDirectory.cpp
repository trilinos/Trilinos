/*Paul
27-July-2002 Templated for OrdinalType.
06-August-2002 Changed to images.
21-Sept-2002 Comm/Platform split
23-Nov-2002 Updated for myImageID moved back to Comm
06-Feb-2003 Updated const syntax
*/

#include "Tpetra_ElementSpace.hpp"

namespace Tpetra {

// default constructor
template<typename OrdinalType>
SerialDirectory<OrdinalType>::SerialDirectory(ElementSpace<OrdinalType> const& elementSpace) 
	: Object("Tpetra::Directory[Serial]") 
	, ElementSpace_(&elementSpace) {}

// copy constructor
template<typename OrdinalType>
SerialDirectory<OrdinalType>::SerialDirectory(const SerialDirectory<OrdinalType>& directory) 
	: Object(directory.label()) 
	, ElementSpace_(directory.ElementSpace_) {}

// destructor
template<typename OrdinalType>
SerialDirectory<OrdinalType>::~SerialDirectory() {}

// query method
template<typename OrdinalType>
void SerialDirectory<OrdinalType>::getDirectoryEntries(OrdinalType numEntries, 
																											 OrdinalType const* globalEntries, 
																											 OrdinalType* images, 
																											 OrdinalType* localEntries) const {
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
