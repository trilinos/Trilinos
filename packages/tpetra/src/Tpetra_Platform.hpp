/*Paul
12-Oct-2002 Updated for Common->Compiler_Directives renaming.
30-Oct-2002 Updated for Compiler_Directives -> ConfigDefs renaming.
12-Nov-2002 Rewritten with new templating scheme.
19-Nov-2002 myImageID and numImages moved back to Comm.
23-Nov-2002 Distributor methods added.
*/

#ifndef _TPETRA_PLATFORM_HPP_
#define _TPETRA_PLATFORM_HPP_

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_Directory.hpp"
#include "Tpetra_Distributor.hpp"

namespace Tpetra {

template<typename OrdinalType> class ElementSpace;
// Comm, Directory, and Distributor are not forward declared because they are used as return types. 

//! Tpetra::Platform: The Tpetra Platform Abstract Base Class
/*! Platform is an abstract base class. It should never be called directly.
		Rather, an implementation of Platform, such as SerialPlatform, should be used instead.
		Platform is used to generate Comm, Distributor, and Directory instances. It also manages 
		platform-specific information, such as how inter-image communication is implemented.
		An implementation of Platform, such as SerialPlatform, will create corresponding classes,
		such as SerialComm and SerialDistributor. These will then be cast to their base class,
		and passed back to other Tpetra modules. As a result, other Tpetra modules don't need to know
		anything about the platform they're running on, or any implementation-specific details.
*/

template<typename ScalarType, typename OrdinalType>
class Platform {
public:

		//@{ \name Constructor/Destructor Methods
		//! Destructor
		virtual ~Platform() {};
		//@}

		//@{ \name Class Creation and Accessor Methods
		//! Comm Instances
		virtual Comm<ScalarType, OrdinalType>* createScalarComm() const = 0;
		virtual Comm<OrdinalType, OrdinalType>* createOrdinalComm() const = 0;
		//! Distributor Instances
		virtual Distributor<ScalarType, OrdinalType>* createScalarDistributor() const = 0;
		virtual Distributor<OrdinalType, OrdinalType>* createOrdinalDistributor() const = 0;
		//! Directory Instance
		virtual Directory<OrdinalType>* createDirectory(const ElementSpace<OrdinalType>& ElementSpace) const = 0;
		//@}

		//@{ \name I/O Methods
		//! printInfo
		virtual void printInfo(ostream& os) const = 0;
		//@}

	}; // Platform class

} // namespace Tpetra

#endif // _TPETRA_PLATFORM_HPP_
