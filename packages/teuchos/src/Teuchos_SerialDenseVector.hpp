// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER


#ifndef _TEUCHOS_SERIALDENSEVECTOR_HPP_
#define _TEUCHOS_SERIALDENSEVECTOR_HPP_

/*! \file Teuchos_SerialDenseVector.hpp
    \brief Templated serial dense vector class
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Object.hpp" 
#include "Teuchos_SerialDenseMatrix.hpp"

/*! \class Teuchos::SerialDenseVector
    \brief This class creates and provides basic support for dense vectors of templated type as a specialization of Teuchos::SerialDenseMatrix.  Additional methods for the SerialDenseVector class, like mathematical methods, can be found documented in SerialDenseMatrix.
*/
namespace Teuchos {

  template<typename OrdinalType, typename ScalarType>
  class SerialDenseVector : public SerialDenseMatrix<OrdinalType,ScalarType> {
    
  public:
    //! @name Constructor/Destructor methods.
  //@{ 

    //! Default Constructor
    /*! Creates an empty vector of no length.  The Sizing methods should be used to size this matrix.  Values of this matrix should be set using the [] or the () operators.
    */
    SerialDenseVector();

    //! Shaped Constructor
    /*!
	\param length - Number of elements in this vector.
	\param zeroOut - Initializes values to 0 if true (default)

    Creates a shaped vector of length \c length.  All values are initialized to 0 when \c zeroOut is true.
    Values of this matrix should be set using the [] or the () operators.
    */
    SerialDenseVector(OrdinalType length, bool zeroOut = true);

    //! Shaped Constructor with Values
    /*!
	\param CV - Enumerated type set to Teuchos::Copy or Teuchos::View.
	\param values - Pointer to an array of ScalarType of the given \c length.
	\param length - Length of vector to be constructed.
    */
    SerialDenseVector(DataAccess CV, ScalarType* values, OrdinalType length);

    //! Copy Constructor
    SerialDenseVector(const SerialDenseVector<OrdinalType,ScalarType>& Source);

    //! Destructor
    virtual ~SerialDenseVector ();
  //@}

  //! @name Sizing methods.
  //@{ 

    //! Size method for changing the size of a SerialDenseVector, initializing entries to zero.
    /*!
	\param length - The length of the new vector.

	This allows the user to define the length of a SerialDenseVector at any point.
	This method can be called at any point after construction.  Any values previously in
	this object will be destroyed and the resized vector starts with all zero values.
    */
    int size(OrdinalType length_in)
      {return(SerialDenseMatrix<OrdinalType, ScalarType>::shape(length_in, 1));}

    //! Same as <tt>size()</tt> except leaves values uninitialized.
    int sizeUninitialized(OrdinalType length_in)
      {return(SerialDenseMatrix<OrdinalType, ScalarType>::shapeUninitialized(length_in, 1));}

    //! Resizing method for changing the size of a SerialDenseVector, keeping the entries.
    /*!
	\param length - The length of the new vector.
	This allows the user to redefine the length of a SerialDenseVector at any point.
	This method can be called at any point after construction.  Any values previously in
	this object will be copied to the resized vector.
    */	
    int resize(OrdinalType length_in)
      {return(SerialDenseMatrix<OrdinalType,ScalarType>::reshape(length_in, 1));}
  //@}

  //! @name Set methods.
  //@{

    //! Set all values in the matrix to a constant value.
    /*!
      \param value - Value to use;
    */
    SerialDenseVector<OrdinalType, ScalarType>& operator= (const ScalarType value) { this->putScalar(value); return(*this); }
  //@}

  //! @name Comparison methods.
  //@{ 
    //! Equality of two matrices.
    /*! \return True if \e this vector and \c Operand are of the same length and have the same entries, else False will be returned.
    */
    bool operator == (const SerialDenseVector<OrdinalType, ScalarType> &Operand) const;

    //! Inequality of two matrices.
    /*! \return True if \e this vector and \c Operand are not of the same length or do not have the same entries, else False will be returned.
    */
    bool operator != (const SerialDenseVector<OrdinalType, ScalarType> &Operand) const;
  //@}

  //! @name Set methods.
  //@{ 

    //! Copies values from one vector to another.
    /*!
	The operator= copies the values from one existing SerialDenseVector to
	another.  If \c Source is a view (i.e. CV = Teuchos::View), then this
	method will return a view.  Otherwise, it will return a copy of \c Source.
	\e this will be resized if it is not large enough to copy \c Source into.
    */
    SerialDenseVector<OrdinalType,ScalarType>& operator = (const SerialDenseVector<OrdinalType,ScalarType>& Source);
  //@}

  //! @name Accessor methods.
  //@{ 
    //! Element access method (non-const).
    /*! Returns the ith element if x(i) is specified, the expression x[i] will return the same element.
	\return (*this)(index)
	\warning The validity of \c index will only be checked if Teuchos is configured with --enable-teuchos-abc.
    */
    ScalarType& operator () (OrdinalType index);
    
    //! Element access method (const).
    /*! Returns the ith element if x(i) is specified, the expression x[i] will return the same element.
	\return (*this)(index)
	\warning The validity of \c index will only be checked if Teuchos is configured with --enable-teuchos-abc.
    */
    const ScalarType& operator () (OrdinalType index) const;

    //! Element access method (non-const).
    /*! Returns the ith element if x[i] is specified, the expression x(i) will return the same element.
	\return (*this)[index]
	\warning The validity of \c index will only be checked if Teuchos is configured with --enable-teuchos-abc.
    */
    ScalarType& operator [] (OrdinalType index);

    //! Element access method (const).
    /*! Returns the ith element if x[i] is specified, the expression x(i) will return the same element.
    	\return (*this)[index]
	\warning The validity of \c index will only be checked if Teuchos is configured with --enable-teuchos-abc.
    */
    const ScalarType& operator [] (OrdinalType index) const;

  //@}

  //! @name Mathematical methods.
  //@{ 
    //! Compute the dot product of \c this vector and \c x.
    ScalarType dot( const SerialDenseVector<OrdinalType,ScalarType> &x) const;
  //@}

  //! @name Attribute methods.
  //@{ 
    //! Returns the length of this vector.
    OrdinalType length() const {return(this->numRows_);}
  //@}

  //! @name I/O methods.
  //@{ 
    //! Print method.  Define the behavior of the std::ostream << operator inherited from the Object class.
    virtual void print(std::ostream& os) const;
  //@}
};

  template<typename OrdinalType, typename ScalarType>
  SerialDenseVector<OrdinalType, ScalarType>::SerialDenseVector() : SerialDenseMatrix<OrdinalType,ScalarType>() {}

  template<typename OrdinalType, typename ScalarType>
  SerialDenseVector<OrdinalType, ScalarType>::SerialDenseVector( OrdinalType length_in, bool zeroOut ) : SerialDenseMatrix<OrdinalType,ScalarType>( length_in, 1, zeroOut ) {}

  template<typename OrdinalType, typename ScalarType>
  SerialDenseVector<OrdinalType, ScalarType>::SerialDenseVector(DataAccess CV, ScalarType* values_in, OrdinalType length_in) : 
    SerialDenseMatrix<OrdinalType,ScalarType>( CV, values_in, length_in, length_in, 1 ) {}

  template<typename OrdinalType, typename ScalarType>
  SerialDenseVector<OrdinalType, ScalarType>::SerialDenseVector(const SerialDenseVector<OrdinalType, ScalarType> &Source) :
    SerialDenseMatrix<OrdinalType,ScalarType>( Source ) {}

  template<typename OrdinalType, typename ScalarType>
  SerialDenseVector<OrdinalType, ScalarType>::~SerialDenseVector() {}
  
  template<typename OrdinalType, typename ScalarType>
  SerialDenseVector<OrdinalType, ScalarType>& SerialDenseVector<OrdinalType,ScalarType>::operator = (const SerialDenseVector<OrdinalType, ScalarType>& Source) 
  {
    SerialDenseMatrix<OrdinalType,ScalarType>::operator=(Source); 
    return(*this);
  }

  template<typename OrdinalType, typename ScalarType>
  bool SerialDenseVector<OrdinalType, ScalarType>::operator == (const SerialDenseVector<OrdinalType, ScalarType> &Operand) const 
  {
    bool result = 1;
    if(this->numRows_ != Operand.numRows_)
      {
	result = 0;
      }
    else
      {
	OrdinalType i;
	for(i = 0; i < this->numRows_; i++) {
	  if((*this)(i) != Operand(i))
	    {
	      return 0;
	    }
	}
      }
    return result;
  }

  template<typename OrdinalType, typename ScalarType>
  bool SerialDenseVector<OrdinalType, ScalarType>::operator != (const SerialDenseVector<OrdinalType, ScalarType> &Operand) const
  {
    return !((*this)==Operand);
  }

  template<typename OrdinalType, typename ScalarType>
  ScalarType SerialDenseVector<OrdinalType, ScalarType>::dot( const SerialDenseVector<OrdinalType, ScalarType> &x) const
  {
    TEST_FOR_EXCEPTION(this->numRows_!= x.numRows_, std::invalid_argument,
    "SerialDenseVector<T>::dot : " << 
    "Number of rows " << this->numRows_ << " not equal to x.numRows_ "<< x.numRows() );
   
    // Compute the dot product and return the result. 
    return BLAS<OrdinalType, ScalarType>::DOT(this->numRows_, this->values(), 1, x.values(), 1); 
  } 

  template<typename OrdinalType, typename ScalarType>
  void SerialDenseVector<OrdinalType, ScalarType>::print(std::ostream& os) const
  {
    os << std::endl;
    if(this->valuesCopied_)
      os << "Values_copied : yes" << std::endl;
    else
      os << "Values_copied : no" << std::endl;
      os << "Length : " << this->numRows_ << std::endl;
    if(this->numRows_ == 0) {
      os << "(std::vector is empty, no values to display)" << std::endl;
    } else {
      for(OrdinalType i = 0; i < this->numRows_; i++) {
	  os << (*this)(i) << " ";
      }
      os << std::endl;
    }
  }

  //----------------------------------------------------------------------------------------------------
  //   Accessor methods 
  //----------------------------------------------------------------------------------------------------

  template<typename OrdinalType, typename ScalarType>
  inline ScalarType& SerialDenseVector<OrdinalType, ScalarType>::operator () (OrdinalType index)
  {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    this->checkIndex( index );
#endif
    return(this->values_[index]);
  }
  
  template<typename OrdinalType, typename ScalarType>
  inline const ScalarType& SerialDenseVector<OrdinalType, ScalarType>::operator () (OrdinalType index) const
  {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    this->checkIndex( index );
#endif
    return(this->values_[index]);
  }
  
  template<typename OrdinalType, typename ScalarType>
  inline const ScalarType& SerialDenseVector<OrdinalType, ScalarType>::operator [] (OrdinalType index) const
  {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    this->checkIndex( index );
#endif
    return(this->values_[index]);
  }
  
  template<typename OrdinalType, typename ScalarType>
  inline ScalarType& SerialDenseVector<OrdinalType, ScalarType>::operator [] (OrdinalType index)
  {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    this->checkIndex( index );
#endif
    return(this->values_[index]);
  }

} // namespace Teuchos

#endif /* _TEUCHOS_SERIALDENSEVECTOR_HPP_ */
