#ifndef _TPETRA_PLATFORM_H_
#define _TPETRA_PLATFORM_H_

#include "Tpetra_Common.h"
#include "Tpetra_Comm.h"

namespace Tpetra {

	template<typename OrdinalType> class ElementSpace;
	template<typename OrdinalType> class Directory;

	//! Tpetra::Platform: The Tpetra Platform Abstract Base Class
	/*! Platform is an abstract base class. It should never be called directly.
		Rather, an implementation of Platform, such as SerialPlatform, should be used instead.
		Logically, Platform is a pure virtual class, but due to the ways in which templates and 
		virtual functions work together in C++, it is not actually implemented that way.
	*/

	class Platform {
	public:
		//@{ \name Constructor/Destructor Methods
		//! Constructor
		//! Destructor
		virtual ~Platform();
		//@}


		//@{ \name Class Creation and Accessor Methods
		//! Comm Instance
				template<typename PacketType, typename OrdinalType>
				Comm<PacketType, OrdinalType>* createComm(PacketType dummyP, OrdinalType dummyO) const;
		//! Directory Instance
				template<typename OrdinalType>
				Directory<OrdinalType>* createDirectory(const ElementSpace<OrdinalType>& ElementSpace) const;
		//@}


		//@{ \name Platform Info Methods
		//! getMyImageID
		virtual int getMyImageID() const;
		//! getNumImages
		virtual int getNumImages() const;
		//@}

		//@{ \name I/O Methods
		//! printInfo
		virtual void printInfo(ostream& os) const;
		//@}

	}; // Platform class

} // namespace Tpetra

#endif // _TPETRA_PLATFORM_H_
