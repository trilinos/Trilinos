#ifndef _TPETRA_SERIALPLATFORM_HPP_
#define _TPETRA_SERIALPLATFORM_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_Platform.hpp"
#include "Tpetra_SerialComm.hpp"
#include "Tpetra_SerialDirectory.hpp"
#include "Tpetra_SerialDistributor.hpp"

namespace Tpetra {

// forward definition
template<typename OrdinalType> class ElementSpace;

	//! Tpetra::SerialPlatform: Serial Implementation of the Platform class.

 template<typename OrdinalType, typename ScalarType>
	class SerialPlatform : public Object, public virtual Platform<OrdinalType, ScalarType> {
	public:

		//@{ \name Constructor/Destructor Methods
		//! Constructor
		SerialPlatform() : Object("Tpetra::Platform[Serial]") {};
		//! Copy constructor
		SerialPlatform(SerialPlatform<OrdinalType, ScalarType> const& platform) : Object(platform.label()) {};
		//! Destructor
		~SerialPlatform() {};
		//! Clone constructor
		Platform<OrdinalType, ScalarType>* clone() const {
			Platform<OrdinalType, ScalarType>* platform = static_cast<Platform<OrdinalType, ScalarType>*>
				(new SerialPlatform<OrdinalType, ScalarType>(*this));
			return(platform);
		};
		//@}

		//@{ \name Class Creation and Accessor Methods

		//! Comm Instances
		Comm<ScalarType, OrdinalType>* createScalarComm() const {
			// static_cast casts SerialComm* to Comm*
			Comm<ScalarType, OrdinalType>* comm = static_cast<Comm<ScalarType, OrdinalType>*>(new SerialComm<ScalarType, OrdinalType>());
			return(comm);
		};
		Comm<OrdinalType, OrdinalType>* createOrdinalComm() const {
			// static_cast casts SerialComm* to Comm*
			Comm<OrdinalType, OrdinalType>* comm = static_cast<Comm<OrdinalType, OrdinalType>*>(new SerialComm<OrdinalType, OrdinalType>());
			return(comm);
		};

		//! Distributor Instances
		Distributor<ScalarType, OrdinalType>* createScalarDistributor() const {
			// static_cast casts SerialDistributor* to Distributor*
			Distributor<ScalarType, OrdinalType>* distributor = 
				static_cast<Distributor<ScalarType, OrdinalType>*>(new SerialDistributor<ScalarType, OrdinalType>());
			return(distributor);
		};
		Distributor<OrdinalType, OrdinalType>* createOrdinalDistributor() const {
			// static_cast casts SerialDistributor* to Distributor*
			Distributor<OrdinalType, OrdinalType>* distributor = 
				static_cast<Distributor<OrdinalType, OrdinalType>*>(new SerialDistributor<OrdinalType, OrdinalType>());
			return(distributor);
		};

		//! Directory Instance
		Directory<OrdinalType>* createDirectory(ElementSpace<OrdinalType> const& elementSpace) const {
		  // static_cast casts SerialDirectory* to Directory*
		  Directory<OrdinalType>* dir = static_cast<Directory<OrdinalType>*>(new SerialDirectory<OrdinalType>(elementSpace)); 
			return(dir);
		};

		//@}

		//@{ \name I/O Methods
		//! print - implements Tpetra::Object virtual print method.
		void print(ostream& os) const { os << label() << endl;};

		//! printInfo - implements Tpetra::Platform virtual printInfo method.
		void printInfo(ostream& os) const {print(os);};
		//@}

	}; // SerialPlatform class

} // namespace Tpetra

#endif // _TPETRA_SERIALPLATFORM_HPP_
