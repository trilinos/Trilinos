#ifndef _TPETRA_SERIALPLATFORM_H_
#define _TPETRA_SERIALPLATFORM_H_

#include "Tpetra_Object.h"
#include "Tpetra_Platform.h"
#include "Tpetra_SerialComm.h"
#include "Tpetra_SerialDirectory.h"

namespace Tpetra {

	template<typename OrdinalType> class ElementSpace;

	//! Tpetra::SerialPlatform: Serial Implementation of the Platform class.

	class SerialPlatform : public Object, public virtual Platform {
	public:

		//@{ \name Constructor/Destructor Methods
		//! Constructor
		SerialPlatform() : Object("Tpetra::Platform[Serial]") {};
		//! Destructor
		~SerialPlatform() {};
		//@}

		//@{ \name Class Creation and Accessor Methods

		//! Comm Instance
		template<typename PacketType, typename OrdinalType>
		Comm<PacketType, OrdinalType>* createComm(PacketType dummyP, OrdinalType dummyO) const {
			// static_cast casts SerialComm* to Comm*
			Comm<PacketType, OrdinalType>* comm = static_cast<Comm<PacketType, OrdinalType>*>(new SerialComm<PacketType, OrdinalType>());
			return(comm);
		};

		//! Directory Instance
		template<typename OrdinalType>
		Directory<OrdinalType>* createDirectory(const ElementSpace<OrdinalType>& ElementSpace) const {
		  // static_cast casts SerialDirectory* to Directory*
		  Directory<OrdinalType>* dir = static_cast<Directory<OrdinalType>*>(new SerialDirectory<OrdinalType>(ElementSpace)); 
		  return(dir);
		};

		//@}

		//@{ \name Platform Info Methods
		//! getMyImageID
		int getMyImageID() const {return(0);};
		//! getNumImages
		int getNumImages() const {return(1);};
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
