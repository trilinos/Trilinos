/*Paul
28-Aug-2002 Initial writeup.
*/

//#include "Tpetra_SerialComm.h"
//#include "Tpetra_SerialDirectory.h"

namespace Tpetra {

	// default constructor
	SerialPlatform::SerialPlatform() : Object("Tpetra::Platform[Serial]") {}

	// destructor
	SerialPlatform::~SerialPlatform() {}

	// Comm instance creator
	//template<typename PacketType, typename OrdinalType>
	//void SerialPlatform::createComm(PacketType dummyP, OrdinalType dummyO) const {
	//PacketType p2 = dummyP + 1;
	//OrdinalType o2 = dummyO + 1;
		//Comm<PacketType, OrdinalType>* comm = static_cast<Comm<PacketType, OrdinalType>*>(new SerialComm<PacketType, OrdinalType>());
		// static_cast casts SerialComm* to Comm*
		//return(comm);
	//}

	// Directory instance creator
	//	template<typename OrdinalType>
	//	Directory<OrdinalType>* SerialPlatform::createDirectory(const ElementSpace<OrdinalType>& ElementSpace) const {
	//	Directory<OrdinalType>* dir = static_cast<Directory<OrdinalType>*>(new SerialDirectory<OrdinalType>(ElementSpace)); 
		// static_cast casts SerialDirectory* to Directory*
	//	return(dir);
	//}
 
	// Distributor instance creator (not yet implemented)
	
	// print
	void SerialPlatform::print(ostream& os) const {
		os << "::Memory Image " << getMyImageID() << " of " << getNumImages() << " total images";
	}

} // namespace Tpetra
