/*Paul
12-Oct-2002 Updated for Common->Compiler_Directives renaming.
30-Oct-2002 Updated for Compiler_Directives -> ConfigDefs renaming.
*/

#ifndef _TPETRA_PLATFORM_H_
#define _TPETRA_PLATFORM_H_

#include "Tpetra_ConfigDefs.h"
#include "Tpetra_Comm.h"
#include "Tpetra_Directory.h"

namespace Tpetra {

template<typename OrdinalType> class ElementSpace;

//! Tpetra::Platform: The Tpetra Platform Abstract Base Class
/*! Platform is an abstract base class. It should never be called directly.
		Rather, an implementation of Platform, such as SerialPlatform, should be used instead.
		Logically, Platform is a pure virtual class, but due to the ways in which templates and 
		virtual functions work together in C++, it is not actually implemented that way.
*/

template<typename PacketType, typename OrdinalType>
class Platform {
public:

		//@{ \name Constructor/Destructor Methods
		//! Destructor
		virtual ~Platform() {};
		//@}

		//@{ \name Platform Info Methods
		//! getMyImageID
		virtual int getMyImageID() const = 0;
		//! getNumImages
		virtual int getNumImages() const = 0;
		//@}

		//@{ \name Class Creation and Accessor Methods
		//! Comm Instance
		virtual Comm<PacketType, OrdinalType>* createComm() const = 0;
		//! Directory Instance
		virtual Directory<OrdinalType>* createDirectory(const ElementSpace<OrdinalType>& ElementSpace) const = 0;
		//@}

		//@{ \name I/O Methods
		//! printInfo
		void printInfo(ostream& os) const {cout << "ERR: Platform method called.\n";};
		//@}

	}; // Platform class

} // namespace Tpetra

#endif // _TPETRA_PLATFORM_H_
