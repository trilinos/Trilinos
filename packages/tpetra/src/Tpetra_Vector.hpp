// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef _TPETRA_VECTOR_HPP_
#define _TPETRA_VECTOR_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_VectorSpace.hpp"
#include <Teuchos_CompObject.hpp>
#include <Teuchos_BLAS.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <vector>
#include <numeric>

namespace Tpetra {

//! Tpetra::Vector: A class for constructing and using sparse vectors.

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
		<li> +3  Cannot perform that operation on an empty vector.
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
		, scalarArray_(VectorSpace.getNumMyEntries())
		, seed_(Teuchos::ScalarTraits<ScalarType>::zero())
	{
		ScalarType const scalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		OrdinalType const length = getNumMyEntries();
		for(OrdinalType i = ordinalZero; i < length; i++)
			scalarArray_[i] = scalarZero;
	}
  
  //! Set object values from user array. Throws an exception if an incorrect number of entries are specified.
	Vector(ScalarType* vectorEntries, OrdinalType numEntries, VectorSpace<OrdinalType, ScalarType> const& VectorSpace)
		: Object("Tpetra::Vector")
		, BLAS_()
		, VectorSpace_(VectorSpace)
		, scalarArray_(VectorSpace.getNumMyEntries())
		, seed_(Teuchos::ScalarTraits<ScalarType>::zero())
	{
		OrdinalType const length = getNumMyEntries();
		if(numEntries != length)
			throw reportError("numEntries = " + toString(numEntries) + ".  Should be = " + toString(length) + ".", -1);
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		for(OrdinalType i = ordinalZero; i < length; i++)
			scalarArray_[i] = vectorEntries[i];
	}

	//! Copy constructor.
  Vector(Vector<OrdinalType, ScalarType> const& Source)
		: Object(Source.label())
		, BLAS_(Source.BLAS_)
		, VectorSpace_(Source.VectorSpace_)
		, scalarArray_(Source.scalarArray_)
		, seed_(Source.seed_)
	{};

  //! Destructor.  
  ~Vector() {};

  //@}

	//@{ \name Post-Construction Modification Routines

	//! Submit entries. Values submitted will be summed with existing values.
	void submitEntries(OrdinalType numEntries, OrdinalType* indices, ScalarType* values) {
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		for(OrdinalType i = 0; i < numEntries; i++)
			scalarArray_[indices[i]] = values[i];
	}

	//! Set all entries to scalarValue.
	void setAllToScalar(ScalarType const value) {
		OrdinalType const max = getNumMyEntries();
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		for(OrdinalType i = ordinalZero; i < max; i++)
			scalarArray_[i] = value;
	}

	//! Set all entries to random values.
	void setAllToRandom() {
		OrdinalType const max = getNumMyEntries();
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		for(OrdinalType i = ordinalZero; i < max; i++)
			scalarArray_[i] = Teuchos::ScalarTraits<ScalarType>::random();
	}

	//@}


	//@{ \name Extraction Methods

	//! Put vector entries into user array (copy)
	void extractCopy(ScalarType* userArray) const {
		OrdinalType const max = getNumMyEntries();
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		for(OrdinalType i = ordinalZero; i < max; i++)
			userArray[i] = scalarArray_[i];
	}

	//! Put pointers to vector entries into user array (view)
	void extractView(ScalarType** userPointerArray) const {
		OrdinalType const max = getNumMyEntries();
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		for(OrdinalType i = ordinalZero; i < max; i++)
			userPointerArray[i] = &scalarArray_[i];
	}

	//@}


	//@{ \name Mathematical Methods

	//! Returns result of dot product, \e this = X.Y.
	void dotProduct(Vector<OrdinalType, ScalarType> const& x, Vector<OrdinalType, ScalarType> const& y) const;

	//! Changes this vector to elementwise absolute values of x.
	void absoluteValue(Vector<OrdinalType, ScalarType> const& x) {
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		OrdinalType const length = getNumMyEntries();
		
		for(OrdinalType i = ordinalZero; i < length; i++)
			scalarArray_[i] = Teuchos::ScalarTraits<ScalarType>::magnitude(x.scalarArray_[i]);
	}

  //! Changes this vector to element-wise reciprocal values of x.
  void reciprocal(Vector<OrdinalType, ScalarType> const& x) {
		if(! vectorSpace().isCompatible(x.vectorSpace()))
			throw reportError("Vector sizes do not match.", 2);

		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		ScalarType const scalarOne = Teuchos::ScalarTraits<ScalarType>::one();
		OrdinalType const length = getNumMyEntries();

		for(OrdinalType i = ordinalZero; i < length; i++)
			scalarArray_[i] = scalarOne / x[i];
	}

  //! Scale the current values of a vector, \e this = scalarThis*\e this.
  void scale(ScalarType scalarThis) {
	  OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
	  OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
	  OrdinalType const length = getNumMyEntries();
	  
	  BLAS_.SCAL(length, scalarThis, &scalarArray_[ordinalZero], ordinalOne);
  }

  //! Replace vector values with scaled values of x, \e this = scalarX*x.
  void scale(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x) {
	  OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
	  OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
	  OrdinalType const length = getNumMyEntries();
	  
	  // this = x
	  scalarArray_ = x.scalarArray_;
	  
	  // this = this * scalarX
	  BLAS_.SCAL(length, scalarX, &scalarArray_[ordinalZero], ordinalOne);
  }

  //! Update vector values with scaled values of x, \e this = scalarThis*\e this + scalarX*x.
  void update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, ScalarType scalarThis) {
	  if(! vectorSpace().isCompatible(x.vectorSpace()))
		  throw reportError("Vector sizes do not match.", 2);
	  
	  OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
	  OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
	  OrdinalType const length = getNumMyEntries();
	  
	  // calculate this *= scalarThis
	  BLAS_.SCAL(length, scalarThis, &scalarArray_[ordinalZero], ordinalOne);
	  
	  // calculate this += scalarX * x
	  BLAS_.AXPY(length, scalarX, &x.scalarArray_[ordinalZero], ordinalOne, &scalarArray_[ordinalZero], ordinalOne);
  }

  //! Update vector with scaled values of x and y, \e this = scalarThis*\e this + scalarX*x + scalarY*y.
  void update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, ScalarType scalarY, 
			  Vector<OrdinalType, ScalarType> const& y, ScalarType scalarThis) {
	  if(!vectorSpace().isCompatible(x.vectorSpace()) ||
	     !vectorSpace().isCompatible(y.vectorSpace()))
		  throw reportError("Vector sizes do not match.", 2);
	  
	OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
	OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
	OrdinalType const length = getNumMyEntries();
		  
	// calculate this *= scalarThis
	BLAS_.SCAL(length, scalarThis, &scalarArray_[ordinalZero], ordinalOne);
		  
	// calculate this += scalarX * x
	BLAS_.AXPY(length, scalarX, &x.scalarArray_[ordinalZero], ordinalOne, &scalarArray_[ordinalZero], ordinalOne);
	
	// calculate this += scalarY * y
	BLAS_.AXPY(length, scalarY, &y.scalarArray_[ordinalZero], ordinalOne, &scalarArray_[ordinalZero], ordinalOne);
  }

  //! Compute 1-norm of vector.
  ScalarType norm1() const;

  //! Compute 2-norm of vector.
  ScalarType norm2() const;

  //! Compute Inf-norm of vector.
  ScalarType normInf() const;

  //! Compute Weighted 2-norm (RMS Norm) of vector.
  ScalarType normWeighted(Vector<OrdinalType, ScalarType> const& weights) const;

  //! Compute minimum value of vector.
  ScalarType minValue() const {
		return(*(min_element(scalarArray_.begin(), scalarArray_.end()))); // use STL min_element, takes constant time
	}

  //! Compute maximum value of vector.
  ScalarType maxValue() const {
		return(*(max_element(scalarArray_.begin(), scalarArray_.end()))); // use STL max_element, takes constant time
	}

  //! Compute mean (average) value of vector.
  ScalarType meanValue() const {
		ScalarType const scalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
		ScalarType length = getNumMyEntries(); // implicit cast from OT to ST
		ScalarType total = accumulate(scalarArray_.begin(), scalarArray_.end(), scalarZero); // use STL accumulate, takes linear time
		return(total / length);
	}

	//! Vector multiplication (elementwise) 
	/*! \e this = scalarThis*\e this + scalarXY*x@y, where @ represents elementwise multiplication. */
	void elementwiseMultiply(ScalarType scalarXY, Vector<OrdinalType, ScalarType> const& x, 
							 Vector<OrdinalType, ScalarType> const& y, ScalarType scalarThis) {
		if(!vectorSpace().isCompatible(x.vectorSpace()) ||
	       !vectorSpace().isCompatible(y.vectorSpace()))
			throw reportError("Vector sizes do not match.", 2);

		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
		OrdinalType const length = getNumMyEntries();
		
		// calculate this *= scalarThis
		BLAS_.SCAL(length, scalarThis, &scalarArray_[ordinalZero], ordinalOne);

		// calculate x@y into temp vector
		vector<ScalarType> xytemp(length);
		transform(x.scalarArray_.begin(), x.scalarArray_.end(), y.scalarArray_.begin(), xytemp.begin(), multiplies<ScalarType>());

		// calculate this = scalarXY * temp + this
		BLAS_.AXPY(length, scalarXY, &xytemp[ordinalZero], ordinalOne, &scalarArray_[ordinalZero], ordinalOne);
	}

	//! Reciprocal multiply (elementwise)
	/*! \e this = scalarThis*\e this + scalarXY*y@x, where @ represents elementwise division. */
	void elementwiseReciprocalMultiply(ScalarType scalarXY, Vector<OrdinalType, ScalarType> const& x, 
									Vector<OrdinalType, ScalarType> const& y, ScalarType scalarThis) {
		if(!vectorSpace().isCompatible(x.vectorSpace()) ||
	       !vectorSpace().isCompatible(y.vectorSpace()))
			throw reportError("Vector sizes do not match.", 2);
		
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
		OrdinalType const length = getNumMyEntries();
		
		// calculate this *= scalarThis
		BLAS_.SCAL(length, scalarThis, &scalarArray_[ordinalZero], ordinalOne);
	
		// calculate y@x into temp vector
		vector<ScalarType> xytemp(length);
		transform(y.scalarArray_.begin(), y.scalarArray_.end(), x.scalarArray_.begin(), xytemp.begin(), divides<ScalarType>());
		
		// calculate this += scalarXY * temp
		BLAS_.AXPY(length, scalarXY, &xytemp[ordinalZero], ordinalOne, &scalarArray_[ordinalZero], ordinalOne);
	}

	//@}


	//@{ \name Random number utilities

	//! Get seed
	ScalarType getSeed() const {
	 return(seed_);
	}

	//! Set seed
	void setSeed(ScalarType seed) {
		seed_ = seed;
	}

	//@}


	//@{ \name Element access methods

	//! [] operator, nonconst version
	ScalarType& operator[](OrdinalType index) {
		return(scalarArray_[index]);
	}

	//! [] operator, const version
	ScalarType const& operator[](OrdinalType index) const {
		return(scalarArray_[index]);
	}

	//@}


	//@{ \name Attribute access methods

	//! Returns number of vector entries owned by this image.
	OrdinalType getNumMyEntries() const {
		return(vectorSpace().getNumMyEntries());
	}

	//! Returns number of vector entries across all images.
	OrdinalType getNumGlobalEntries() const {
		return(vectorSpace().getNumGlobalEntries());
	}

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
					os << scalarArray_[i] << " ";
				os << endl << endl;
			}
		}
	}
	//@}

	//@{ \name Misc. 

	//! Returns a const reference to the VectorSpace this Vector belongs to.
	VectorSpace<OrdinalType, ScalarType> const& vectorSpace() const {
		return(VectorSpace_);
	}

	//@}

private:

	Teuchos::BLAS<OrdinalType, ScalarType> BLAS_;
	VectorSpace<OrdinalType, ScalarType> VectorSpace_;
	std::vector<ScalarType> scalarArray_;
	ScalarType seed_;

}; // class Vector

} // namespace Tpetra

#endif /* _TPETRA_VECTOR_HPP_ */
