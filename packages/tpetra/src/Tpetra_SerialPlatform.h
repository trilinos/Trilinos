#ifndef _TPETRA_SERIALPLATFORM_H_
#define _TPETRA_SERIALPLATFORM_H_

#include "Tpetra_Object.h"
#include "Tpetra_Platform.h"
#include "Tpetra_SerialComm.h"
#include "Tpetra_SerialDirectory.h"

namespace Tpetra {

// forward definition
template<typename OrdinalType> class ElementSpace;

	//! Tpetra::SerialPlatform: Serial Implementation of the Platform class.

 template<typename PacketType, typename OrdinalType>
	class SerialPlatform : public Object, public virtual Platform<PacketType, OrdinalType> {
	public:

		//@{ \name Constructor/Destructor Methods
		//! Constructor
		SerialPlatform() : Object("Tpetra::Platform[Serial]") {};
		//! Copy constructor
		SerialPlatform(const SerialPlatform<PacketType, OrdinalType>& Platform) : Object(Platform.label()) {};
		//! Destructor
		~SerialPlatform() {};
		//@}

		//@{ \name Platform Info Methods
		//! getMyImageID
		int getMyImageID() const {return(0);};
		//! getNumImages
		int getNumImages() const {return(1);};
		//@}

		//@{ \name Class Creation and Accessor Methods

		//! Comm Instance
		Comm<PacketType, OrdinalType>* createComm() const {
			// static_cast casts SerialComm* to Comm*
			Comm<PacketType, OrdinalType>* comm = static_cast<Comm<PacketType, OrdinalType>*>(new SerialComm<PacketType, OrdinalType>());
			return(comm);
		};

		//! Directory Instance
		Directory<OrdinalType>* createDirectory(const ElementSpace<OrdinalType>& ElementSpace) const {
		  // static_cast casts SerialDirectory* to Directory*
		  Directory<OrdinalType>* dir = static_cast<Directory<OrdinalType>*>(new SerialDirectory<OrdinalType>(ElementSpace)); 
		  return(dir);
		};

		//@}

		//@{ \name I/O Methods
		//! print - implements Tpetra::Object virtual print method.
		void print(ostream& os) const {
			os << "::Memory Image " << getMyImageID() << " of " << getNumImages() << " total images" << endl;
		};

		//! printInfo - implements Tpetra::Platform virtual printInfo method.
		void printInfo(ostream& os) const {print(os);};
		//@}

	}; // SerialPlatform class

} // namespace Tpetra

#endif // _TPETRA_SERIALPLATFORM_H_
