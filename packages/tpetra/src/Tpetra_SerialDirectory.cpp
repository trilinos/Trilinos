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
