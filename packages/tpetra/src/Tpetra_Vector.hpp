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
#include "Tpetra_ConfigDefs.hpp" // for STL vector, algorithm, numeric
#include "Tpetra_VectorSpace.hpp"
#include <Teuchos_CompObject.hpp>
#include <Teuchos_BLAS.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>

namespace Tpetra {

// forward declaration of VectorData, needed to prevent circular inclusions
// actual #include statement is at the end of this file
template<typename OrdinalType, typename ScalarType> class VectorData;

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
		, VectorData_()
	{
		ScalarType const scalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		OrdinalType const length = VectorSpace.getNumMyEntries();

		VectorData_ = Teuchos::rcp(new VectorData<OrdinalType, ScalarType>(VectorSpace, length, scalarZero));

		for(OrdinalType i = ordinalZero; i < length; i++)
			VectorData_->scalarArray_[i] = scalarZero;
	}
  
  //! Set object values from user array. Throws an exception if an incorrect number of entries are specified.
	Vector(ScalarType* vectorEntries, OrdinalType numEntries, VectorSpace<OrdinalType, ScalarType> const& VectorSpace)
		: Object("Tpetra::Vector")
		, VectorData_()
	{
		ScalarType const scalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		OrdinalType const length = VectorSpace.getNumMyEntries();

		VectorData_ = Teuchos::rcp(new VectorData<OrdinalType, ScalarType>(VectorSpace, length, scalarZero));

		if(numEntries != length)
			throw reportError("numEntries = " + toString(numEntries) + ".  Should be = " + toString(length) + ".", -1);
		for(OrdinalType i = ordinalZero; i < length; i++)
			VectorData_->scalarArray_[i] = vectorEntries[i];
	}

	//! Copy constructor.
  Vector(Vector<OrdinalType, ScalarType> const& Source)
		: Object(Source.label())
		, VectorData_(Source.VectorData_)
	{};

  //! Destructor.  
  ~Vector() {};

  //@}

	//@{ \name Post-Construction Modification Routines

	//! Submit entries. Values submitted will be summed with existing values.
	void submitEntries(OrdinalType numEntries, OrdinalType* indices, ScalarType* values) {
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		for(OrdinalType i = 0; i < numEntries; i++)
			VectorData_->scalarArray_[indices[i]] = values[i];
	}

	//! Set all entries to scalarValue.
	void setAllToScalar(ScalarType const value) {
		OrdinalType const max = getNumMyEntries();
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		for(OrdinalType i = ordinalZero; i < max; i++)
			VectorData_->scalarArray_[i] = value;
	}

	//! Set all entries to random values.
	void setAllToRandom() {
		OrdinalType const max = getNumMyEntries();
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		for(OrdinalType i = ordinalZero; i < max; i++)
			VectorData_->scalarArray_[i] = Teuchos::ScalarTraits<ScalarType>::random();
	}

	//@}


	//@{ \name Extraction Methods

	//! Put vector entries into user array (copy)
	void extractCopy(ScalarType* userArray) const {
		OrdinalType const max = getNumMyEntries();
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		for(OrdinalType i = ordinalZero; i < max; i++)
			userArray[i] = VectorData_->scalarArray_[i];
	}

	//! Put pointers to vector entries into user array (view)
	void extractView(ScalarType** userPointerArray) const {
		OrdinalType const max = getNumMyEntries();
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		for(OrdinalType i = ordinalZero; i < max; i++)
			userPointerArray[i] = &VectorData_->scalarArray_[i];
	}

	//@}


	//@{ \name Mathematical Methods

	//! Returns result of dot product, \e result = this.x
	ScalarType dotProduct(Vector<OrdinalType, ScalarType> const& x) const {
		if(! vectorSpace().isCompatible(x.vectorSpace()))
			throw reportError("Vector sizes do not match.", 2);
		
		OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
		OrdinalType const length = getNumMyEntries();
		
		// call BLAS routine to calculate local dot product
		ScalarType localDP = BLAS().DOT(length, scalarPointer(), ordinalOne, x.scalarPointer(), ordinalOne);
		
		// use Comm call to sum all local dot products
		ScalarType globalDP;
		vectorSpace().comm().sumAll(&localDP, &globalDP, ordinalOne);
        
		// update flops counter: 2n-1
		updateFlops(length + length - ordinalOne);
		
		return(globalDP);
	}

	//! Changes this vector to elementwise absolute values of x.
	void absoluteValue(Vector<OrdinalType, ScalarType> const& x) {
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		OrdinalType const length = getNumMyEntries();
		
		for(OrdinalType i = ordinalZero; i < length; i++)
			VectorData_->scalarArray_[i] = Teuchos::ScalarTraits<ScalarType>::magnitude(x[i]);
	}

  //! Changes this vector to element-wise reciprocal values of x.
  void reciprocal(Vector<OrdinalType, ScalarType> const& x) {
		if(! vectorSpace().isCompatible(x.vectorSpace()))
			throw reportError("Vector sizes do not match.", 2);

		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		ScalarType const scalarOne = Teuchos::ScalarTraits<ScalarType>::one();
		OrdinalType const length = getNumMyEntries();

		for(OrdinalType i = ordinalZero; i < length; i++)
			VectorData_->scalarArray_[i] = scalarOne / x[i];
        
		// update flops counter: n
		updateFlops(length);
	}

  //! Scale the current values of a vector, \e this = scalarThis*\e this.
  void scale(ScalarType scalarThis) {
	  OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
	  OrdinalType const length = getNumMyEntries();
	  
	  BLAS().SCAL(length, scalarThis, scalarPointer(), ordinalOne);
      
		// update flops counter: n
		updateFlops(length);
  }

  //! Replace vector values with scaled values of x, \e this = scalarX*x.
  void scale(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x) {
	  OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
	  OrdinalType const length = getNumMyEntries();
	  
	  // this = x
	  scalarArray() = x.scalarArray();
	  // this = this * scalarX
	  BLAS().SCAL(length, scalarX, scalarPointer(), ordinalOne);
      
		// update flops counter: n
		updateFlops(length);
  }

  //! Update vector values with scaled values of x, \e this = scalarThis*\e this + scalarX*x.
  void update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, ScalarType scalarThis) {
	  if(! vectorSpace().isCompatible(x.vectorSpace()))
		  throw reportError("Vector sizes do not match.", 2);
	  
	  OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
	  OrdinalType const length = getNumMyEntries();
	  
	  // calculate this *= scalarThis
	  BLAS().SCAL(length, scalarThis, scalarPointer(), ordinalOne);
	  
	  // calculate this += scalarX * x
	  BLAS().AXPY(length, scalarX, x.scalarPointer(), ordinalOne, scalarPointer(), ordinalOne);
      
		// update flops counter: 3n
		updateFlops(length + length + length);
  }

  //! Update vector with scaled values of x and y, \e this = scalarThis*\e this + scalarX*x + scalarY*y.
  void update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, ScalarType scalarY, 
			  Vector<OrdinalType, ScalarType> const& y, ScalarType scalarThis) {
	  if(!vectorSpace().isCompatible(x.vectorSpace()) ||
	     !vectorSpace().isCompatible(y.vectorSpace()))
		  throw reportError("Vector sizes do not match.", 2);
	  
		OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
		OrdinalType const length = getNumMyEntries();
		  
		// calculate this *= scalarThis
		BLAS().SCAL(length, scalarThis, scalarPointer(), ordinalOne);
		// calculate this += scalarX * x
		BLAS().AXPY(length, scalarX, x.scalarPointer(), ordinalOne, scalarPointer(), ordinalOne);
		// calculate this += scalarY * y
		BLAS().AXPY(length, scalarY, y.scalarPointer(), ordinalOne, scalarPointer(), ordinalOne);
    
    // update flops counter: 5n
    updateFlops(length + length + length + length + length);
  }

  //! Compute 1-norm of vector.
	ScalarType norm1() const {
		// 1-norm = sum of abs. values of vector entries
		OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
		OrdinalType const length = getNumMyEntries();
		
		// compute local 1-norm
		ScalarType localNorm = BLAS().ASUM(length, scalarPointer(), ordinalOne);
		// call comm's sumAll method to compute global 1-norm
		ScalarType globalNorm;
		vectorSpace().comm().sumAll(&localNorm, &globalNorm, ordinalOne);
		
		// update flops counter: n-1
		updateFlops(length - ordinalOne);
        
		return(globalNorm);
	}

  //! Compute 2-norm of vector.
	ScalarType norm2() const {
		// 2-norm = square root of the sum of the squares of the abs. values of vector entries
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
		OrdinalType const length = getNumMyEntries();
		
		// add up squares of entries
		ScalarType localSum = Teuchos::ScalarTraits<ScalarType>::zero();
		for(OrdinalType i = ordinalZero; i < length; i++)
			localSum += scalarArray()[i] * scalarArray()[i];
		
		// calculate global sum
		ScalarType globalSum;
		vectorSpace().comm().sumAll(&localSum, &globalSum, ordinalOne);
        
		// update flops counter: 2n
		updateFlops(length + length);
		
		// return square root of global sum
		return(Teuchos::ScalarTraits<ScalarType>::squareroot(globalSum));
	}

  //! Compute Infinity-norm of vector.
  ScalarType normInf() const {
	  // inf-norm = abs. value of the max value
	  return(Teuchos::ScalarTraits<ScalarType>::magnitude(maxValue()));
  }

  //! Compute Weighted 2-norm (RMS Norm) of vector.
	ScalarType normWeighted(Vector<OrdinalType, ScalarType> const& weights) const {
		if(!vectorSpace().isCompatible(weights.vectorSpace()))
			throw reportError("Vector sizes do not match.", 2);
		
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
		OrdinalType const length = getNumMyEntries();
		
		// add up this[i] * weights[i]
		ScalarType localSum = Teuchos::ScalarTraits<ScalarType>::zero();
		for(OrdinalType i = ordinalZero; i < length; i++) {
			ScalarType temp = scalarArray()[i] * weights[i];
			localSum += temp * temp;
		}
		
		// get global sum
		ScalarType globalSum;
		vectorSpace().comm().sumAll(&localSum, &globalSum, ordinalOne);
		
		// divide by global length, and then take square root of that
		globalSum /= static_cast<ScalarType>(getNumGlobalEntries());
        
		// update flops counter: 3n
		updateFlops(length + length + length);
		
		return(Teuchos::ScalarTraits<ScalarType>::squareroot(globalSum));
	}

  //! Compute minimum value of vector.
  ScalarType minValue() const {
		return(*(min_element(scalarArray().begin(), scalarArray().end()))); // use STL min_element, takes constant time
  }

  //! Compute maximum value of vector.
  ScalarType maxValue() const {
		return(*(max_element(scalarArray().begin(), scalarArray().end()))); // use STL max_element, takes constant time
  }

  //! Compute mean (average) value of vector.
  ScalarType meanValue() const {
		ScalarType const scalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
		ScalarType length = getNumMyEntries(); // implicit cast from OT to ST
		ScalarType total = accumulate(scalarArray().begin(), scalarArray().end(), scalarZero); // use STL accumulate, takes linear time
        
		// update flops counter: n
		updateFlops(length);
        
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

		// calculate x@y into temp vector
		vector<ScalarType> xytemp(length);
		transform(x.scalarArray().begin(), x.scalarArray().end(), y.scalarArray().begin(), xytemp.begin(), multiplies<ScalarType>());
        
        // calculate this *= scalarThis
		BLAS().SCAL(length, scalarThis, scalarPointer(), ordinalOne);

		// calculate this = scalarXY * temp + this
		BLAS().AXPY(length, scalarXY, &xytemp[ordinalZero], ordinalOne, scalarPointer(), ordinalOne);
        
		// update flops counter: n
		updateFlops(length);
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
	
		// calculate y@x into temp vector
		vector<ScalarType> xytemp(length);
		transform(y.scalarArray().begin(), y.scalarArray().end(), x.scalarArray().begin(), xytemp.begin(), divides<ScalarType>());
        
        // calculate this *= scalarThis
		BLAS().SCAL(length, scalarThis, scalarPointer(), ordinalOne);
		
		// calculate this += scalarXY * temp
		BLAS().AXPY(length, scalarXY, &xytemp[ordinalZero], ordinalOne, scalarPointer(), ordinalOne);
        
		// update flops counter: 2n
		updateFlops(length + length);
	}

	//@}


	//@{ \name Random number utilities

	//! Get seed
	ScalarType getSeed() const {
	 return(VectorData_->seed_);
	}

	//! Set seed
	void setSeed(ScalarType seed) {
		VectorData_->seed_ = seed;
	}

	//@}


	//@{ \name Element access methods

	//! [] operator, nonconst version
	ScalarType& operator[](OrdinalType index) {
		return(VectorData_->scalarArray_[index]);
	}

	//! [] operator, const version
	ScalarType const& operator[](OrdinalType index) const {
		return(VectorData_->scalarArray_[index]);
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
		OrdinalType const myImageID = vectorSpace().comm().getMyImageID();
		OrdinalType const numImages = vectorSpace().comm().getNumImages();
		OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();
        OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		
		for (OrdinalType imageCtr = ordinalZero; imageCtr < numImages; imageCtr++) {
			if (myImageID == imageCtr) {
				if (myImageID == ordinalZero) {
					os << "Number of Global Entries  = " << getNumGlobalEntries() << endl;
				}
				os <<   "ImageID = " << myImageID << endl;
				os <<           "Number of Local Entries   = " << getNumMyEntries() << endl;
				os <<           "Contents: ";
				for(OrdinalType i = ordinalZero; i < getNumMyEntries(); i++)
					os << VectorData_->scalarArray_[i] << " ";
				os << endl;
			}
		}
	}
    
    void printValues(ostream& os) const {
        for(OrdinalType i = Teuchos::OrdinalTraits<OrdinalType>::zero(); i < getNumMyEntries(); i++)
            os << VectorData_->scalarArray_[i] << " ";
        os << endl;        
    }
	//@}

	//@{ \name Misc. 

	//! Returns a const reference to the VectorSpace this Vector belongs to.
	VectorSpace<OrdinalType, ScalarType> const& vectorSpace() const {
		return(VectorData_->VectorSpace_);
	}

	//@}

private:

	// Accessor for BLAS
	Teuchos::BLAS<OrdinalType, ScalarType> const& BLAS() const {
		return(VectorData_->BLAS_);
	}

	// Accessors for scalarArray
	std::vector<ScalarType>& scalarArray() {
		return(VectorData_->scalarArray_);
	}
	std::vector<ScalarType>const & scalarArray() const{
		return(VectorData_->scalarArray_);
	}

	// Returns pointer to ScalarType array inside of scalarArray
	ScalarType* scalarPointer() {
		return(&VectorData_->scalarArray_[Teuchos::OrdinalTraits<OrdinalType>::zero()]);
	}
	ScalarType const* scalarPointer() const {
		return(&VectorData_->scalarArray_[Teuchos::OrdinalTraits<OrdinalType>::zero()]);
	}

	Teuchos::RefCountPtr< VectorData<OrdinalType, ScalarType> > VectorData_;

}; // class Vector

} // namespace Tpetra

#include "Tpetra_VectorData.hpp"

#endif /* _TPETRA_VECTOR_HPP_ */
