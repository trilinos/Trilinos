#ifndef _TPETRA_VECTOR_HPP_
#define _TPETRA_VECTOR_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_VectorSpace.hpp"
#include <Teuchos_CompObject.hpp>
#include <Teuchos_BLAS.hpp>

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
		: BLAS_()
		//, VectorSpace_(VectorSpace)
	{
		cout << "vector constructor called" << endl;
	};
  
  //! Set object values from user array. Throws an exception if an incorrect number of entries are specified.
	//Vector(ScalarType* vectorEntries, VectorSpace<OrdinalType, ScalarType> const& vectorSpace);

	//! Copy constructor.
  Vector(Vector<OrdinalType, ScalarType> const& Vector)
		: BLAS_(Vector.BLAS_)
		//, VectorSpace_(Vector.vectorSpace_)
	{
		cout << "vector copy constructor called" << endl;
	};

  //! Destructor.  
  ~Vector() {
		cout << "vector destructor called" << endl;
	};

  //@}


	//@{ \name Element access methods

	//! [] operator, nonconst version
	ScalarType& operator[](OrdinalType index) {
		cout << "operator[] called, nonconst version" << endl;
	};

	//! [] operator, const version
	ScalarType const& operator[](OrdinalType index) const {
		cout << "operator[] called, const version" << endl;
	};

	//@}


	//@{ \name Attribute access methods

	//! Returns number of vector entries owned by this image.
	OrdinalType getNumMyEntries() const {
		//return(VectorSpace().getNumMyEntries());
		cout << "getNumMyEntries called" << endl;
	};

	//! Returns number of vector entries across all images.
	OrdinalType getNumEntries() const {
		//return(VectorSpace().getNumMyEntries());
		cout << "getNumEntries called" << endl;
	};

	//@}


	//@{ \name I/O methods

	//! Print method, used by overloaded << operator.
	void print(ostream& os) const {
		os << "Vector print function called" << endl;
	};
	//@}

private:

	//VectorSpace<OrdinalType, ScalarType> const& VectorSpace() const {
	//	return(&VectorSpace_);
	//};

	Teuchos::BLAS<OrdinalType, ScalarType> BLAS_;
	//VectorSpace<OrdinalType, ScalarType> VectorSpace_;
	//OrdinalType NumMyEntries_;
	//OrdinalType NumGlobalEntries_;

}; // class Vector

} // namespace Tpetra

#endif /* _TPETRA_VECTOR_HPP_ */
