#ifndef _TPETRA_VECTOR_HPP_
#define _TPETRA_VECTOR_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_VectorSpace.hpp"
#include <Teuchos_CompObject.hpp>
#include <Teuchos_BLAS.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <vector>

namespace Tpetra {

//! Tpetra::Vector: A class for constructing and using vectors.

/*! Vector is templated on ScalarType for the vector entries, and on OrdinalType 
	  for the vector indices. A VectorSpace object is needed for all Vector objects.

		Vector entries can only be accessed through their local index values. (This 
		is due to the fact that a VectorSpace can be created with either an ElementSpace
		or a BlockElementSpace, and there is no way to globally access a BES Point.)
		Global index values can be converted to local indices by using the 
		VectorSpace::getLocalIndex method.

		Note that for most of the mathematical methods that set \e this to the result of
		an operation on vectors passed as parameters, the \e this vector can be used 
		as one of the parameters (unless otherwise specified).
		
    Vector error codes (positive for non-fatal, negative for fatal):
		<ol>
		<li> +1  Specified vector index not found on this image.
		<li> +2  Vector sizes do not match.
		<li> -1  Invalid number of entries passed to constructor.
		<li> -99 Internal Vector error. Contact developer.
		</ol>
*/

template<typename OrdinalType, typename ScalarType>
class Vector : public Teuchos::CompObject, public Object
{

public:
  
  //@{ \name Constructor/Destructor Methods

  //! Sets all vector entries to zero.
  Vector(VectorSpace<OrdinalType, ScalarType> const& VectorSpace) 
		: Object("Tpetra::Vector")
		, BLAS_()
		, VectorSpace_(VectorSpace)
		, A_(VectorSpace.getNumMyEntries())
	{
		ScalarType const zero = Teuchos::ScalarTraits<ScalarType>::zero();
		OrdinalType const length = getNumMyEntries();
		for(OrdinalType i = 0; i < length; i++)
			A_[i] = zero;
	};
  
  //! Set object values from user array. Throws an exception if an incorrect number of entries are specified.
	Vector(ScalarType* vectorEntries, OrdinalType numEntries, VectorSpace<OrdinalType, ScalarType> const& VectorSpace)
		: Object("Tpetra::Vector")
		, BLAS_()
		, VectorSpace_(VectorSpace)
		, A_(VectorSpace.getNumMyEntries())
	{
		OrdinalType const length = getNumMyEntries();
		if(numEntries != length)
			throw reportError("numEntries = " + toString(numEntries) + ".  Should be = " + toString(length) + ".", -1);
		for(OrdinalType i = 0; i < length; i++)
			A_[i] = vectorEntries[i];
	};

	//! Copy constructor.
  Vector(Vector<OrdinalType, ScalarType> const& Source)
		: Object(Source.label())
		, BLAS_(Source.BLAS_)
		, VectorSpace_(Source.VectorSpace_)
		, A_(Source.A_)
	{};

  //! Destructor.  
  ~Vector() {};

  //@}


	//@{ \name Element access methods

	//! [] operator, nonconst version
	ScalarType& operator[](OrdinalType index) {
		return(A_[index]);
	};

	//! [] operator, const version
	ScalarType const& operator[](OrdinalType index) const {
		return(A_[index]);
	};

	//@}


	//@{ \name Attribute access methods

	//! Returns number of vector entries owned by this image.
	OrdinalType getNumMyEntries() const {
		return(vectorSpace().getNumMyEntries());
	};

	//! Returns number of vector entries across all images.
	OrdinalType getNumGlobalEntries() const {
		return(vectorSpace().getNumGlobalEntries());
	};

	//@}


	//@{ \name I/O methods

	//! Print method, used by overloaded << operator.
	void print(ostream& os) const {
		
		OrdinalType myImageID = vectorSpace().comm().getMyImageID();
		OrdinalType numImages = vectorSpace().comm().getNumImages();
		
		for (int imageCtr = 0; imageCtr < numImages; imageCtr++) {
			if (myImageID == imageCtr) {
				if (myImageID == 0) {
					os << endl << "Number of Global Entries  = " << getNumGlobalEntries() << endl;
				}
				os << endl <<   "ImageID = " << myImageID << endl;
				os <<           "Number of Local Entries   = " << getNumMyEntries() << endl;
				os <<           "Contents: ";
				for(OrdinalType i = 0; i < getNumMyEntries(); i++)
					os << A_[i] << " ";
				os << endl << endl;
			}
		}
	};
	//@}

private:

	//! Returns a const reference to the VectorSpace this Vector belongs to.
	VectorSpace<OrdinalType, ScalarType> const& vectorSpace() const {
		return(VectorSpace_);
	};

	Teuchos::BLAS<OrdinalType, ScalarType> BLAS_;
	VectorSpace<OrdinalType, ScalarType> VectorSpace_;
	std::vector<ScalarType> A_;

}; // class Vector

} // namespace Tpetra

#endif /* _TPETRA_VECTOR_HPP_ */
