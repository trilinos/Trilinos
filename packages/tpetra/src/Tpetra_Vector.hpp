#ifndef _TPETRA_VECTOR_HPP_
#define _TPETRA_VECTOR_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_CompObject.hpp"
#include "Tpetra_BLAS.hpp"

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
class Vector : public CompObject, public Object, public BLAS<OrdinalType, ScalarType>
{

public:
  
  //@{ \name Constructor/Destructor Methods

  //! Sets all vector entries to zero.
  Vector(VectorSpace<OrdinalType, ScalarType> const& vectorSpace);
  
  //! Set object values from user array. Throws an exception if an incorrect number of entries are specified.
	Vector(ScalarType* vectorEntries, VectorSpace<OrdinalType, ScalarType> const& vectorSpace);

	//! Copy constructor.
  Vector(Vector<OrdinalType, ScalarType> const& vector);

  //! Destructor.  
  ~Vector();

  //@}


	//@{ \name Post-Construction Modification Routines

	//! Combine entries with given local indices. Exact behavior is set using Tpetra::CombineMode.
	void combineEntries(CombineMode cm, OrdinalType numEntries, OrdinalType* indices, ScalarType* values);

	//! Set all entries to scalarValue.
	void setAllToScalar(ScalarType const value);

	//! Set all entries to random values.
	void setAllToRandom();

	//@}


	//@{ \name Extraction Methods

	//! Put vector entries into user array (copy)
	void extractCopy(ScalarType* userArray) const;

	//! Put pointers to vector entries into user array (view)
	void extractView(ScalarType** userPointerArray) const;

	//@}


	//@{ \name Mathematical Methods

	//! Returns result of dot product, \e this = X.Y.
	void dotProduct(Vector<OrdinalType, ScalarType> const& x, Vector<OrdinalType, ScalarType> const& y) const;

	//! Changes this vector to elementwise absolute values of x.
	void absoluteValue(Vector<OrdinalType, ScalarType> const& x);

  //! Changes this vector to element-wise reciprocal values of x.
  void reciprocal(Vector<OrdinalType, ScalarType> const& x);

  //! Scale the current values of a vector, \e this = scalarThis*\e this.
  void scale(ScalarType scalarThis);

  //! Replace vector values with scaled values of x, \e this = scalarX*x.
  void scale(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x);

  //! Update vector values with scaled values of x, \e this = scalarThis*\e this + scalarX*x.
  void update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, ScalarType scalarThis);

  //! Update vector with scaled values of x and y, \e this = scalarThis*\e this + scalarX*x + scalarY*y.
  void update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, ScalarType scalarY, 
							Vector<OrdinalType, ScalarType> const& y, ScalarType scalarThis);

  //! Compute 1-norm of vector.
  ScalarType norm1() const;

  //! Compute 2-norm of vector.
  ScalarType norm2() const;

  //! Compute Inf-norm of vector.
  ScalarType normInf() const;

  //! Compute Weighted 2-norm (RMS Norm) of vector.
  ScalarType normWeighted(Vector<OrdinalType, ScalarType> const& weights) const;

  //! Compute minimum value of vector.
  ScalarType minValue() const;

  //! Compute maximum value of vector.
  ScalarType maxValue() const;

  //! Compute mean (average) value of vector.
  ScalarType meanValue() const;

	//! Vector multiplication (elementwise) 
	/*! \e this = scalarThis*\e this + scalarXY*x@y, where @ represents elementwise multiplication. */
	void elementwiseMultiply(ScalarType scalarXY, Vector<OrdinalType, ScalarType> const& x, 
													 Vector<OrdinalType, ScalarType> const& y, ScalarType scalarThis);

	//! Reciprocal multiply (elementwise)
	/*! \e this = scalarThis*\e this + scalarXY*y@x, where @ represents elementwise division. */
	void elementwiseReciprocalMultiply(ScalarType scalarXY, Vector<OrdinalType, ScalarType> const& x, 
																		 Vector<OrdinalType, ScalarType> const& y, ScalarType scalarThis);

	//@}


	//@{ \name Random number utilities

	//! Get seed
	ScalarType getSeed() const;

	//! Set seed
	void setSeed(ScalarType seed);

	//@}


	//@{ \name Element access methods

	//! [] operator, nonconst version
	ScalarType& operator[](OrdinalType index);

	//! [] operator, const version
	ScalarType const& operator[](OrdinalType index) const;

	//@}


	//@{ \name Attribute access methods

	//! Returns number of vector entries owned by this image.
	OrdinalType getNumMyEntries() const;

	//! Returns number of vector entries across all images.
	OrdinalType getNumEntries() const;

	//@}


	//@{ \name I/O methods

	//! Print method, used by overloaded << operator.
	void print(ostream& os) const;
	//@}

}; // class Vector

} // namespace Tpetra

// begin Tpetra_Vector.cpp
//=======================================================================

namespace Tpetra {

// Constructor 1: all values zeroed
//=======================================================================
template<typename OrdinalType, typename ScalarType>
Vector<OrdinalType, ScalarType>::Vector(VectorSpace<OrdinalType, ScalarType> const& vectorSpace)
  
// Constructor 2: values from user array
//=======================================================================
template<typename OrdinalType, typename ScalarType>
Vector<OrdinalType, ScalarType>::Vector(ScalarType* vectorEntries, 
																				VectorSpace<OrdinalType, ScalarType> const& vectorSpace)

// Constructor 3: Copy constructor
//=======================================================================
template<typename OrdinalType, typename ScalarType>
Vector<OrdinalType, ScalarType>::Vector(Vector<OrdinalType, ScalarType> const& vector)

// Destructor
//=======================================================================
template<typename OrdinalType, typename ScalarType>
Vector<OrdinalType, ScalarType>::~Vector()

// combine
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::combineEntries(CombineMode cm, OrdinalType numEntries, 
																										 OrdinalType* indices, ScalarType* values)

// set all to scalar
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::setAllToScalar(ScalarType const value)

// Set all to random
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::setAllToRandom()

// extract copy
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::extractCopy(ScalarType* userArray) const

// extract view
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::extractView(ScalarType** userPointerArray) const

// dot product
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::dotProduct(Vector<OrdinalType, ScalarType> const& x, 
																								 Vector<OrdinalType, ScalarType> const& y) const

// absolute value
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::absoluteValue(Vector<OrdinalType, ScalarType> const& x)

// reciprocal
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::reciprocal(Vector<OrdinalType, ScalarType> const& x)

// scale: this = scalarThis * this
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::scale(ScalarType scalarThis)

// scale: this = scalarX * x
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::scale(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x)

// update: this = scalarThis * this + scalarX * x
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, 
																						 ScalarType scalarThis)

// update: x and y
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, 
																						 ScalarType scalarY, Vector<OrdinalType, ScalarType> const& y, 
																						 ScalarType scalarThis)

// 1-norm
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::norm1() const

// 2-norm
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::norm2() const

// Inf-norm
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::normInf() const

// Weighted 2-norm (RMS Norm)
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::normWeighted(Vector<OrdinalType, ScalarType> const& weights) const

// minimum
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::minValue() const

// maximum
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::maxValue() const

// mean
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::meanValue() const

// Vector multiplication (elementwise) 
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::elementwiseMultiply(ScalarType scalarXY, 
																													Vector<OrdinalType, ScalarType> const& x, 
																													Vector<OrdinalType, ScalarType> const& y, 
																													ScalarType scalarThis)

// Reciprocal multiply (elementwise)
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::elementwiseReciprocalMultiply(ScalarType scalarXY, 
																																		Vector<OrdinalType, ScalarType> const& x, 
																																		Vector<OrdinalType, ScalarType> const& y, 
																																		ScalarType scalarThis)

// getSeed
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::getSeed() const

// setSeed
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::setSeed(ScalarType seed)

// [] operator, nonconst version
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType& Vector<OrdinalType, ScalarType>::operator[](OrdinalType index)

// [] operator, const version
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType const& Vector<OrdinalType, ScalarType>::operator[](OrdinalType index) const

// GetNumMyEntries
//=======================================================================
template<typename OrdinalType, typename ScalarType>
OrdinalType Vector<OrdinalType, ScalarType>::getNumMyEntries() const

// getNumGlobalEntries
//=======================================================================
template<typename OrdinalType, typename ScalarType>
OrdinalType Vector<OrdinalType, ScalarType>::getNumEntries() const

// print
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::print(ostream& os) const

} //namespace Tpetra

//=======================================================================
// end Tpetra_Vector.cpp

#endif /* _TPETRA_VECTOR_HPP_ */
