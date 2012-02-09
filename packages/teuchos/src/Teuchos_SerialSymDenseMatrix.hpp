
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

#ifndef _TEUCHOS_SERIALSYMDENSEMATRIX_HPP_
#define _TEUCHOS_SERIALSYMDENSEMATRIX_HPP_
/*! \file Teuchos_SerialSymDenseMatrix.hpp
  \brief Templated serial, dense, symmetric matrix class.
*/

#include "Teuchos_CompObject.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Assert.hpp"

/*! \class Teuchos::SerialSymDenseMatrix
    \brief This class creates and provides basic support for symmetric, positive-definite dense matrices of templated type.

    The Teuchos_SerialSymDenseMatrix class enables the construction and use of 
    symmetric, positive-definite, dense matrices of templated type.

    The Teuchos::SerialSymDenseMatrix class is intended to provide full-featured support for solving 
    linear and eigen system problems for symmetric positive-definite (SPD) matrices.  It is written on 
    top of BLAS and LAPACK and thus has excellent performance and numerical capabilities.  Using this 
    class, one can either perform simple factorizations and solves or apply all the tricks available in 
    LAPACK to get the best possible solution for very ill-conditioned problems.

    <b>Teuchos::SerialSymDenseMatrix vs. Teuchos::LAPACK</b>

    The Teuchos::LAPACK class provides access to most of the same functionality as Teuchos::SerialSymDenseMatrix.
    The primary difference is that Teuchos::LAPACK is a "thin" layer on top of LAPACK and 
    Teuchos::SerialSymDenseMatrix attempts to provide easy access to the more sophisticated aspects of 
    solving dense linear and eigensystems.
<ul>
<li> When you should use Teuchos::LAPACK:  If you are simply looking for a convenient wrapper around 
the Fortran LAPACK routines and you have a well-conditioned problem, you should probably use Teuchos::LAPACK directly.
<li> When you should use Teuchos::SerialSymDenseMatrix: If you want to (or potentially want to) solve ill-conditioned 
     problems or want to work with a more object-oriented interface, you should probably use Teuchos::SerialSymDenseMatrix.    
</ul>

<b>Constructing Teuchos::SerialSymDenseMatrix Objects</b>

There are three Teuchos::SerialSymDenseMatrix constructors.  The first constructs a zero-sized object 
which should be made to appropriate length using the Shape() or Reshape() functions and then filled with 
the [] or () operators. The second is a constructor that accepts user data as a 2D array, the third is a 
copy constructor. The second constructor has two data access modes (specified by the Teuchos::DataAccess argument):
<ol>
  <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
       user data is not needed after construction.
  <li> View mode - Creates a "view" of the user data. In this case, the
       user data is required to remain intact for the life of the object.
</ol>

\warning View mode is \e extremely dangerous from a data hiding perspective. Therefore, we strongly encourage 
users to develop code using Copy mode first and only use the View mode in a secondary optimization phase.

<b>Extracting Data from Teuchos::SerialSymDenseMatrix Objects</b>

Once a Teuchos::SerialSymDenseMatrix is constructed, it is possible to view the data via access functions.

\warning Use of these access functions cam be \e extremely dangerous from a data hiding perspective.

<b>Vector and Utility Functions</b>

Once a Teuchos::SerialSymDenseMatrix is constructed, several mathematical functions can be applied to
the object.  Specifically:
<ul>
  <li> Multiplication.
  <li> Norms.
</ul>

*/

/** \example DenseMatrix/cxx_main_sym.cpp
    This is an example of how to use the Teuchos::SerialSymDenseMatrix class.
*/

namespace Teuchos {
 
template<typename OrdinalType, typename ScalarType>
class SerialSymDenseMatrix : public CompObject, public Object, public BLAS<OrdinalType,ScalarType>
{
 public:

  //! Typedef for ordinal type
  typedef OrdinalType ordinalType;
  //! Typedef for scalar type
  typedef ScalarType scalarType;

  //! @name Constructor/Destructor Methods
  //@{ 
  //! Default constructor; defines a zero size object.
  /*!
    Teuchos::SerialSymDenseMatrix objects defined by the default constructor 
    should be sized with the Shape() or Reshape() functions.  
    Values should be defined by using the [] or ()operators.

    Note: By default the active part of the matrix is assumed to be the lower triangular part.
    To set the upper part as active, call SetUpper(). See Detailed Description section for further discussion.
   */
  SerialSymDenseMatrix();

  //! Basic constructor; defines a matrix of \c numRowsCols size and (optionally) initializes it.
  /*!
    \param numRowsCols - Number of rows and columns in the matrix.
    \param zeroOut - Initializes values to 0 if true (default)

    Creates a shaped matrix with \c numRowsCols rows and cols.  All values are initialized to 0 when \c zeroOut is true.
    Values of this matrix should be set using the [] or the () operators.

    \note By default the active part of the matrix is assumed to be the lower triangular part.
    To set the upper part as active, call SetUpper(). See Detailed Description section for further discussion.
  */
  SerialSymDenseMatrix(OrdinalType numRowsCols, bool zeroOut = true);

  //! Set object values from two-dimensional array.
  /*!
    \param CV - Enumerated type set to Teuchos::Copy or Teuchos::View.

    \param values - Pointer to an array of ScalarType.  The first column starts at \c values,
		the second at \c values+stride, etc.
    \param stride - The stride between the columns of the matrix in memory.
    \param numRowsCols - Number of rows and columns in the matrix.

    \note By default the active part of the matrix is assumed to be the lower triangular part.
    To set the upper part as active, call SetUpper(). See Detailed Description section for further discussion.
  */
  SerialSymDenseMatrix(DataAccess CV, bool upper, ScalarType* values, OrdinalType stride, OrdinalType numRowsCols);
  
  //! Teuchos::SerialSymDenseMatrix copy constructor.
  SerialSymDenseMatrix(const SerialSymDenseMatrix<OrdinalType, ScalarType> &Source);

  //! Submatrix Copy Constructor
  /*! 
    \param CV - Enumerated type set to Teuchos::Copy or Teuchos::View.
    \param Source - Reference to another dense matrix from which values are to be copied.
    \param numRowCols - The number of rows and columns in this matrix.
    \param startRowCol - The row and column of \c Source from which the submatrix copy should start. 

    Creates a shaped matrix with \c numRowCols rows and columns, which is a submatrix of \c Source.
    If \c startRowCol are not given, then the submatrix is the leading submatrix of \c Source.
  */
  SerialSymDenseMatrix(DataAccess CV, const SerialSymDenseMatrix<OrdinalType, ScalarType> &Source, OrdinalType numRowCols, OrdinalType startRowCol=0);

  //! Teuchos::SerialSymDenseMatrix destructor.  
  virtual ~SerialSymDenseMatrix ();
  //@}

  //! @name Shaping Methods
  //@{ 

  //! Set dimensions of a Teuchos::SerialSymDenseMatrix object; init values to zero.
  /*!
    \param numRowsCols - Number of rows and columns in object.

    Allows user to define the dimensions of a Teuchos::SerialSymDenseMatrix at any point. This function can
    be called at any point after construction.  Any values that were previously in this object are
    destroyed and the resized matrix starts off with all zero values.

    \return Integer error code, set to 0 if successful.
  */
  int shape(OrdinalType numRowsCols); 

  //! Set dimensions of a Teuchos::SerialSymDenseMatrix object; don't initialize values.
  /*!
    \param numRowsCols - Number of rows and columns in object.

    Allows user to define the dimensions of a Teuchos::SerialSymDenseMatrix at any point. This function can
    be called at any point after construction.  Any values that were previously in this object are
    destroyed.  The resized matrix has uninitialized values.

    \return Integer error code, set to 0 if successful.
  */
  int shapeUninitialized(OrdinalType numRowsCols); 

  //! Reshape a Teuchos::SerialSymDenseMatrix object.
  /*!
    \param numRowsCols - Number of rows and columns in object.
    
    Allows user to define the dimensions of a Teuchos::SerialSymDenseMatrix at any point. This function can
    be called at any point after construction.  Any values that were previously in this object are
    copied into the new shape.  If the new shape is smaller than the original, the upper left portion
    of the original matrix (the principal submatrix) is copied to the new matrix.
    
    \return Integer error code, set to 0 if successful.
  */
  int reshape(OrdinalType numRowsCols);

  //! Specify that the lower triangle of the \e this matrix should be used.
  /*! \warning This may necessitate the matrix values be copied from the upper to lower portion of the matrix.
  */
  void setLower();

  //! Specify that the upper triangle of the \e this matrix should be used.
  /*! \warning This may necessitate the matrix values be copied from the lower to upper portion of the matrix.
  */
  void setUpper();

  //@}

  //! @name Set methods.
  //@{ 

  //! Copies values from one matrix to another.
  /*!
    The operator= copies the values from one existing SerialSymDenseMatrix to another. 
    If \c Source is a view (i.e. CV = Teuchos::View), then this method will
    return a view.  Otherwise, it will return a copy of \c Source.  \e this object
    will be resized if it is not large enough to copy \c Source into.
  */	
  SerialSymDenseMatrix<OrdinalType, ScalarType>& operator= (const SerialSymDenseMatrix<OrdinalType, ScalarType>& Source);

  //! Copies values from one matrix to another.
  /*!
    Copies the values from one existing SerialSymDenseMatrix to another if the dimension of both matrices are the same.  
    If not, \e this matrix will be returned unchanged.
  */	
  SerialSymDenseMatrix<OrdinalType, ScalarType>& assign (const SerialSymDenseMatrix<OrdinalType, ScalarType>& Source);

  //! Set all values in the matrix to a constant value.
  /*!
    \param value - Value to use;
  */
  SerialSymDenseMatrix<OrdinalType, ScalarType>& operator= (const ScalarType value) { putScalar(value); return(*this); }

  //! Set all values in the matrix to a constant value.
  /*!
    \param value - Value to use; zero if none specified.
    \param fullMatrix - set full matrix entries to \c value, not just active portion of symmetric matrix.
    \return Integer error code, set to 0 if successful.
  */
  int putScalar( const ScalarType value = Teuchos::ScalarTraits<ScalarType>::zero(), bool fullMatrix = false );

  //! Set all values in the active area (upper/lower triangle) of this matrix to be random numbers.
  /*! \note The diagonal will be the sum of the off diagonal elements, plus a bias, so the matrix is SPD.
   */
  int random( const ScalarType bias = 0.1*Teuchos::ScalarTraits<ScalarType>::one() );

  //@}

  //! @name Accessor methods.
  //@{ 

  //! Element access method (non-const).
  /*! Returns the element in the ith row and jth column if A(i,j) is specified. 

    \return Element from the specified \c rowIndex row and \c colIndex column.
    \note If the requested element lies in the inactive part of the matrix, then A(j,i) will be returned.
    \warning The validity of \c rowIndex and \c colIndex will only be checked if Teuchos is
    configured with --enable-teuchos-abc.
  */
  ScalarType& operator () (OrdinalType rowIndex, OrdinalType colIndex);

  //! Element access method (const).
  /*! Returns the element in the ith row and jth column if A(i,j) is specified.

    \return Element from the specified \c rowIndex row and \c colIndex column.
    \note If the requested element lies in the inactive part of the matrix, then A(j,i) will be returned.
    \warning The validity of \c rowIndex and \c colIndex will only be checked if Teuchos is
    configured with --enable-teuchos-abc.
  */
  const ScalarType& operator () (OrdinalType rowIndex, OrdinalType colIndex) const;

  //! Returns the pointer to the ScalarType data array contained in the object. 
  /*! \note The matrix values are only guaranteed to be stored in the active area of the matrix (upper/lower).
  */
  ScalarType* values() const { return(values_); }

  //@}

  //! @name Query methods
  //@{ 

  //! Returns true if upper triangular part of \e this matrix has and will be used.
  bool upper() const {return(upper_);};

  //! Returns character value of UPLO used by LAPACK routines.
  char UPLO() const {return(UPLO_);};
  //@}

  //! @name Mathematical Methods
  //@{ 

  //! Inplace scalar-matrix product A = \c alpha*A.
  /*! Scale a matrix, entry-by-entry using the value \e alpha.  This method is sensitive to
      the UPLO() parameter.

   \param alpha - Scalar to multiply with A.

  */
  SerialSymDenseMatrix<OrdinalType, ScalarType>& operator*= (const ScalarType alpha);

  //! Add another matrix to \e this matrix.
  /*! Add \c Source to \e this if the dimension of both matrices are the same.  If not, \e this matrix
    will be returned unchanged.
  */
  SerialSymDenseMatrix<OrdinalType, ScalarType>& operator+= (const SerialSymDenseMatrix<OrdinalType, ScalarType>& Source);

  //! Subtract another matrix from \e this matrix.
  /*! Subtract \c Source from \e this if the dimension of both matrices are the same.  If not, \e this matrix
    will be returned unchanged.
  */
  SerialSymDenseMatrix<OrdinalType, ScalarType>& operator-= (const SerialSymDenseMatrix<OrdinalType, ScalarType>& Source);

  //@}

  //! @name Comparison methods.
  //@{ 

  //! Equality of two matrices.
  /*! \return True if \e this matrix and \c Operand are of the same shape (rows / columns) and have
    the same entries in the active (upper / lower triangular) area of the matrix, else False will be returned.
  */
  bool operator== (const SerialSymDenseMatrix<OrdinalType, ScalarType> &Operand) const;

  //! Inequality of two matrices.
  /*! \return True if \e this matrix and \c Operand of not of the same shape (rows / columns) or don't
    have the same entries in the active (upper / lower triangular), else False will be returned.
  */
  bool operator!= (const SerialSymDenseMatrix<OrdinalType, ScalarType> &Operand) const;

  //@}

  //! @name Attribute methods.
  //@{ 

  //! Returns the row dimension of this matrix.
  OrdinalType numRows() const { return(numRowCols_); }

  //! Returns the column dimension of this matrix.
  OrdinalType numCols() const { return(numRowCols_); }

  //! Returns the stride between the columns of this matrix in memory.
  OrdinalType stride() const { return(stride_); }

  //! Returns whether this matrix is empty.
  bool empty() const { return(numRowCols_ == 0); }

  //@}

  //! @name Norm methods.
  //@{ 

  //! Returns the 1-norm of the matrix.
  typename ScalarTraits<ScalarType>::magnitudeType normOne() const;

  //! Returns the Infinity-norm of the matrix.
  typename ScalarTraits<ScalarType>::magnitudeType normInf() const;

  //! Returns the Frobenius-norm of the matrix.
  typename ScalarTraits<ScalarType>::magnitudeType normFrobenius() const; 
  //@}

  //! @name I/O methods.
  //@{ 
  //! Print method.  Defines the behavior of the std::ostream << operator inherited from the Object class.
  virtual void print(std::ostream& os) const;

  //@}

 protected:

  // In-place scaling of matrix.
  void scale( const ScalarType alpha );

  // Copy the values from one matrix to the other.  
  void copyMat(bool inputUpper, ScalarType* inputMatrix, OrdinalType inputStride,
	       OrdinalType numRowCols, bool outputUpper, ScalarType* outputMatrix, 
	       OrdinalType outputStride, OrdinalType startRowCol, 
	       ScalarType alpha = ScalarTraits<ScalarType>::zero() );

  // Copy the values from the active triangle of the matrix to the other to make the matrix full symmetric.
  void copyUPLOMat(bool inputUpper, ScalarType* inputMatrix, 
                   OrdinalType inputStride, OrdinalType inputRows);

  void deleteArrays();
  void checkIndex( OrdinalType rowIndex, OrdinalType colIndex = 0 ) const;
  OrdinalType numRowCols_;
  OrdinalType stride_;
  bool valuesCopied_;
  ScalarType* values_;
  bool upper_;
  char UPLO_;


};

//----------------------------------------------------------------------------------------------------
//  Constructors and Destructor
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
SerialSymDenseMatrix<OrdinalType, ScalarType>::SerialSymDenseMatrix()
  : CompObject(), Object("Teuchos::SerialSymDenseMatrix"), numRowCols_(0), stride_(0), valuesCopied_(false), values_(0), upper_(false), UPLO_('L')
{}

template<typename OrdinalType, typename ScalarType>
SerialSymDenseMatrix<OrdinalType, ScalarType>::SerialSymDenseMatrix(OrdinalType numRowCols_in, bool zeroOut)
  : CompObject(), Object("Teuchos::SerialSymDenseMatrix"), numRowCols_(numRowCols_in), stride_(numRowCols_in), valuesCopied_(false), upper_(false), UPLO_('L')
{
  values_ = new ScalarType[stride_*numRowCols_];
  valuesCopied_ = true;
  if (zeroOut == true)
    putScalar( Teuchos::ScalarTraits<ScalarType>::zero(), true );
}
  
template<typename OrdinalType, typename ScalarType>
SerialSymDenseMatrix<OrdinalType, ScalarType>::SerialSymDenseMatrix(
  DataAccess CV, bool upper_in, ScalarType* values_in, OrdinalType stride_in, OrdinalType numRowCols_in
  )
  : CompObject(), Object("Teuchos::SerialSymDenseMatrix"), numRowCols_(numRowCols_in), stride_(stride_in), valuesCopied_(false), 
    values_(values_in), upper_(upper_in)
{
  if (upper_)
    UPLO_ = 'U';
  else
    UPLO_ = 'L';

  if(CV == Copy)
  {
    stride_ = numRowCols_;
    values_ = new ScalarType[stride_*numRowCols_];
    copyMat(upper_in, values_in, stride_in, numRowCols_, upper_, values_, stride_, 0);
    valuesCopied_ = true;
  }
}
  
template<typename OrdinalType, typename ScalarType>
SerialSymDenseMatrix<OrdinalType, ScalarType>::SerialSymDenseMatrix(const SerialSymDenseMatrix<OrdinalType, ScalarType> &Source) : CompObject(), Object("Teuchos::SerialSymDenseMatrix"), numRowCols_(Source.numRowCols_), stride_(Source.numRowCols_), valuesCopied_(true), values_(0), upper_(Source.upper_), UPLO_(Source.UPLO_) 
{
  values_ = new ScalarType[stride_*numRowCols_];
  copyMat(Source.upper_, Source.values_, Source.stride_, numRowCols_, upper_, values_, stride_, 0);
  valuesCopied_ = true;
}
  
template<typename OrdinalType, typename ScalarType>
SerialSymDenseMatrix<OrdinalType, ScalarType>::SerialSymDenseMatrix(
								    DataAccess CV, const SerialSymDenseMatrix<OrdinalType, 
								    ScalarType> &Source, OrdinalType numRowCols_in, OrdinalType startRowCol )
  : CompObject(), Object("Teuchos::SerialSymDenseMatrix"), numRowCols_(numRowCols_in), stride_(Source.stride_), valuesCopied_(false), upper_(Source.upper_), UPLO_(Source.UPLO_)
{
  if(CV == Copy)
  {
    stride_ = numRowCols_in;
    values_ = new ScalarType[stride_ * numRowCols_in];
    copyMat(Source.upper_, Source.values_, Source.stride_, numRowCols_in, upper_, values_, stride_, startRowCol);
    valuesCopied_ = true;
  }
  else // CV == View
  {
    values_ = Source.values_ + (stride_ * startRowCol) + startRowCol;
  }
}
    
template<typename OrdinalType, typename ScalarType>
SerialSymDenseMatrix<OrdinalType, ScalarType>::~SerialSymDenseMatrix()
{
  deleteArrays();
}

//----------------------------------------------------------------------------------------------------
//  Shape methods 
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
int SerialSymDenseMatrix<OrdinalType, ScalarType>::shape( OrdinalType numRowCols_in )
{
  deleteArrays(); // Get rid of anything that might be already allocated
  numRowCols_ = numRowCols_in;
  stride_ = numRowCols_;
  values_ = new ScalarType[stride_*numRowCols_];
  putScalar( Teuchos::ScalarTraits<ScalarType>::zero(), true );
  valuesCopied_ = true;
  return(0);
}

template<typename OrdinalType, typename ScalarType>
int SerialSymDenseMatrix<OrdinalType, ScalarType>::shapeUninitialized( OrdinalType numRowCols_in )
{
  deleteArrays(); // Get rid of anything that might be already allocated
  numRowCols_ = numRowCols_in;
  stride_ = numRowCols_;
  values_ = new ScalarType[stride_*numRowCols_];
  valuesCopied_ = true;
  return(0);
}
  
template<typename OrdinalType, typename ScalarType>
int SerialSymDenseMatrix<OrdinalType, ScalarType>::reshape( OrdinalType numRowCols_in )
{
  // Allocate space for new matrix
  ScalarType* values_tmp = new ScalarType[numRowCols_in * numRowCols_in];
  ScalarType zero = ScalarTraits<ScalarType>::zero();
  for(OrdinalType k = 0; k < numRowCols_in * numRowCols_in; k++)
  {
    values_tmp[k] = zero;
  }
  OrdinalType numRowCols_tmp = TEUCHOS_MIN(numRowCols_, numRowCols_in);
  if(values_ != 0)
  {
    copyMat(upper_, values_, stride_, numRowCols_tmp, upper_, values_tmp, numRowCols_in, 0); // Copy principal submatrix of A to new A
  }
  deleteArrays(); // Get rid of anything that might be already allocated
  numRowCols_ = numRowCols_in;
  stride_ = numRowCols_;
  values_ = values_tmp; // Set pointer to new A
  valuesCopied_ = true;
  return(0);
}
   
//----------------------------------------------------------------------------------------------------
//   Set methods 
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
void SerialSymDenseMatrix<OrdinalType, ScalarType>::setLower() 
{
  // Do nothing if the matrix is already a lower triangular matrix
  if (upper_ != false) {
    copyUPLOMat( true, values_, stride_, numRowCols_ );
    upper_ = false;
    UPLO_ = 'L';
  }
}

template<typename OrdinalType, typename ScalarType>
void SerialSymDenseMatrix<OrdinalType, ScalarType>::setUpper() 
{
  // Do nothing if the matrix is already an upper triangular matrix
  if (upper_ == false) {
    copyUPLOMat( false, values_, stride_, numRowCols_ );
    upper_ = true;
    UPLO_ = 'U';
  }
}

template<typename OrdinalType, typename ScalarType>
int SerialSymDenseMatrix<OrdinalType, ScalarType>::putScalar( const ScalarType value_in, bool fullMatrix )
{
  // Set each value of the dense matrix to "value".
  if (fullMatrix == true) {
    for(OrdinalType j = 0; j < numRowCols_; j++) 
      {
	for(OrdinalType i = 0; i < numRowCols_; i++) 
	  {
	    values_[i + j*stride_] = value_in;
	  }
      }
  }
  // Set the active upper or lower triangular part of the matrix to "value"
  else {
    if (upper_) {
      for(OrdinalType j = 0; j < numRowCols_; j++) {
	for(OrdinalType i = 0; i <= j; i++) {
	  values_[i + j*stride_] = value_in;
	}
      }
    }
    else {
      for(OrdinalType j = 0; j < numRowCols_; j++) {
	for(OrdinalType i = j; i < numRowCols_; i++) {
	  values_[i + j*stride_] = value_in;
	}
      }
    }
  }
  return 0;
}    
    
template<typename OrdinalType, typename ScalarType>
int SerialSymDenseMatrix<OrdinalType, ScalarType>::random( const ScalarType bias )
{
  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MT;

  // Set each value of the dense matrix to a random value.
  std::vector<MT> diagSum( numRowCols_, Teuchos::ScalarTraits<MT>::zero() );
  if (upper_) {
    for(OrdinalType j = 0; j < numRowCols_; j++) {
      for(OrdinalType i = 0; i < j; i++) {
	values_[i + j*stride_] = ScalarTraits<ScalarType>::random();
	diagSum[i] += Teuchos::ScalarTraits<ScalarType>::magnitude( values_[i + j*stride_] );
	diagSum[j] += Teuchos::ScalarTraits<ScalarType>::magnitude( values_[i + j*stride_] );
      }
    }
  }
  else {
    for(OrdinalType j = 0; j < numRowCols_; j++) {
      for(OrdinalType i = j+1; i < numRowCols_; i++) {
	values_[i + j*stride_] = ScalarTraits<ScalarType>::random();
	diagSum[i] += Teuchos::ScalarTraits<ScalarType>::magnitude( values_[i + j*stride_] );
	diagSum[j] += Teuchos::ScalarTraits<ScalarType>::magnitude( values_[i + j*stride_] );
      }
    }
  }
  
  // Set the diagonal.
  for(OrdinalType i = 0; i < numRowCols_; i++) {
    values_[i + i*stride_] = diagSum[i] + bias;
  }
  return 0;
}

template<typename OrdinalType, typename ScalarType>
SerialSymDenseMatrix<OrdinalType,ScalarType>&
SerialSymDenseMatrix<OrdinalType, ScalarType>::operator=( const SerialSymDenseMatrix<OrdinalType,ScalarType>& Source )
{
  if(this == &Source)
    return(*this); // Special case of source same as target
  if((!valuesCopied_) && (!Source.valuesCopied_) && (values_ == Source.values_)) {
    upper_ = Source.upper_;  // Might have to change the active part of the matrix.
    return(*this); // Special case of both are views to same data.
  }
  
  // If the source is a view then we will return a view, else we will return a copy.
  if (!Source.valuesCopied_) {
    if(valuesCopied_)       { 
      // Clean up stored data if this was previously a copy.
      deleteArrays();
    }
    numRowCols_ = Source.numRowCols_; 
    stride_ = Source.stride_;
    values_ = Source.values_;
    upper_ = Source.upper_;
    UPLO_ = Source.UPLO_;
  }
  else {
    // If we were a view, we will now be a copy.
    if(!valuesCopied_) {
      numRowCols_ = Source.numRowCols_;
      stride_ = Source.numRowCols_;
      upper_ = Source.upper_;
      UPLO_ = Source.UPLO_;
      const OrdinalType newsize = stride_ * numRowCols_;
      if(newsize > 0) {
        values_ = new ScalarType[newsize];
        valuesCopied_ = true;
      }
      else {
        values_ = 0;
      }
    }
    // If we were a copy, we will stay a copy.
    else {
      if((Source.numRowCols_ <= stride_) && (Source.numRowCols_ == numRowCols_)) { // we don't need to reallocate
        numRowCols_ = Source.numRowCols_;
	upper_ = Source.upper_;
	UPLO_ = Source.UPLO_;
      }
      else { // we need to allocate more space (or less space)
        deleteArrays();
        numRowCols_ = Source.numRowCols_;
        stride_ = Source.numRowCols_;
	upper_ = Source.upper_;
	UPLO_ = Source.UPLO_;
        const OrdinalType newsize = stride_ * numRowCols_;
        if(newsize > 0) {
          values_ = new ScalarType[newsize];
          valuesCopied_ = true;
        }
      }
    }
    copyMat(Source.upper_, Source.values_, Source.stride_, Source.numRowCols_, upper_, values_, stride_, 0);
  } 
  return(*this);
}

template<typename OrdinalType, typename ScalarType>
SerialSymDenseMatrix<OrdinalType, ScalarType>& SerialSymDenseMatrix<OrdinalType, ScalarType>::operator*= (const ScalarType alpha)
{
  this->scale(alpha);
  return(*this);
}

template<typename OrdinalType, typename ScalarType>
SerialSymDenseMatrix<OrdinalType, ScalarType>& SerialSymDenseMatrix<OrdinalType, ScalarType>::operator+= (const SerialSymDenseMatrix<OrdinalType,ScalarType>& Source )
{
  // Check for compatible dimensions
  if ((numRowCols_ != Source.numRowCols_))
    {
      TEUCHOS_CHK_REF(*this); // Return *this without altering it.
    }
  copyMat(Source.upper_, Source.values_, Source.stride_, numRowCols_, upper_, values_, stride_, 0, ScalarTraits<ScalarType>::one());
  return(*this);
}

template<typename OrdinalType, typename ScalarType>
SerialSymDenseMatrix<OrdinalType, ScalarType>& SerialSymDenseMatrix<OrdinalType, ScalarType>::operator-= (const SerialSymDenseMatrix<OrdinalType,ScalarType>& Source )
{
  // Check for compatible dimensions
  if ((numRowCols_ != Source.numRowCols_))
  {
    TEUCHOS_CHK_REF(*this); // Return *this without altering it.
  }
  copyMat(Source.upper_, Source.values_, Source.stride_, numRowCols_, upper_, values_, stride_, 0, -ScalarTraits<ScalarType>::one());
  return(*this);
}

template<typename OrdinalType, typename ScalarType>
SerialSymDenseMatrix<OrdinalType,ScalarType>& SerialSymDenseMatrix<OrdinalType, ScalarType>::assign (const SerialSymDenseMatrix<OrdinalType,ScalarType>& Source) {
  if(this == &Source)
    return(*this); // Special case of source same as target
  if((!valuesCopied_) && (!Source.valuesCopied_) && (values_ == Source.values_)) {
    upper_ = Source.upper_; // We may have to change the active part of the matrix.
    return(*this); // Special case of both are views to same data.
  }

  // Check for compatible dimensions
  if ((numRowCols_ != Source.numRowCols_))
  {
    TEUCHOS_CHK_REF(*this); // Return *this without altering it.
  }
  copyMat(Source.upper_, Source.values_, Source.stride_, numRowCols_, upper_, values_, stride_, 0 );
  return(*this);
}

//----------------------------------------------------------------------------------------------------
//   Accessor methods 
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
inline ScalarType& SerialSymDenseMatrix<OrdinalType, ScalarType>::operator () (OrdinalType rowIndex, OrdinalType colIndex)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  checkIndex( rowIndex, colIndex );
#endif
  if ( rowIndex <= colIndex ) {
    // Accessing upper triangular part of matrix
    if (upper_)
      return(values_[colIndex * stride_ + rowIndex]);
    else
      return(values_[rowIndex * stride_ + colIndex]);   
  }
  else {
    // Accessing lower triangular part of matrix
    if (upper_)
      return(values_[rowIndex * stride_ + colIndex]);
    else
      return(values_[colIndex * stride_ + rowIndex]);
  }
}
  
template<typename OrdinalType, typename ScalarType>
inline const ScalarType& SerialSymDenseMatrix<OrdinalType, ScalarType>::operator () (OrdinalType rowIndex, OrdinalType colIndex) const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  checkIndex( rowIndex, colIndex );
#endif
  if ( rowIndex <= colIndex ) {
    // Accessing upper triangular part of matrix
    if (upper_)
      return(values_[colIndex * stride_ + rowIndex]);
    else
      return(values_[rowIndex * stride_ + colIndex]);   
  }
  else {
    // Accessing lower triangular part of matrix
    if (upper_)
      return(values_[rowIndex * stride_ + colIndex]);
    else
      return(values_[colIndex * stride_ + rowIndex]);
  }
}  

//----------------------------------------------------------------------------------------------------
//   Norm methods 
//----------------------------------------------------------------------------------------------------
  
template<typename OrdinalType, typename ScalarType>
typename ScalarTraits<ScalarType>::magnitudeType SerialSymDenseMatrix<OrdinalType, ScalarType>::normOne() const
{
  return(normInf());
}
  
template<typename OrdinalType, typename ScalarType>
typename ScalarTraits<ScalarType>::magnitudeType SerialSymDenseMatrix<OrdinalType, ScalarType>::normInf() const
{
  typedef typename ScalarTraits<ScalarType>::magnitudeType MT;

  OrdinalType i, j;
  
  MT sum, anorm = ScalarTraits<MT>::zero();
  ScalarType* ptr;
  
  if (upper_) {
    for (j=0; j<numRowCols_; j++) {
      sum = ScalarTraits<MT>::zero();
      ptr = values_ + j*stride_;
      for (i=0; i<j; i++) {
        sum += ScalarTraits<ScalarType>::magnitude( *ptr++ );
      }
      ptr = values_ + j + j*stride_;
      for (i=j; i<numRowCols_; i++) {
        sum += ScalarTraits<ScalarType>::magnitude( *ptr );
        ptr += stride_;
      }
      anorm = TEUCHOS_MAX( anorm, sum );
    }
  }
  else {
    for (j=0; j<numRowCols_; j++) {
      sum = ScalarTraits<MT>::zero();
      ptr = values_ + j + j*stride_;
      for (i=j; i<numRowCols_; i++) {
        sum += ScalarTraits<ScalarType>::magnitude( *ptr++ );
      }
      ptr = values_ + j;
      for (i=0; i<j; i++) {
        sum += ScalarTraits<ScalarType>::magnitude( *ptr );
        ptr += stride_;
      }
      anorm = TEUCHOS_MAX( anorm, sum );
    }
  }
  return(anorm);
}
  
template<typename OrdinalType, typename ScalarType>
typename ScalarTraits<ScalarType>::magnitudeType SerialSymDenseMatrix<OrdinalType, ScalarType>::normFrobenius() const
{
  typedef typename ScalarTraits<ScalarType>::magnitudeType MT;

  OrdinalType i, j;
  
  MT sum = ScalarTraits<MT>::zero(), anorm = ScalarTraits<MT>::zero();
 
  if (upper_) { 
    for (j = 0; j < numRowCols_; j++) {
      for (i = 0; i < j; i++) {
        sum += ScalarTraits<ScalarType>::magnitude(2.0*values_[i+j*stride_]*values_[i+j*stride_]);
      }
      sum += ScalarTraits<ScalarType>::magnitude(values_[j + j*stride_]*values_[j + j*stride_]);
    }
  }
  else {
    for (j = 0; j < numRowCols_; j++) {
      sum += ScalarTraits<ScalarType>::magnitude(values_[j + j*stride_]*values_[j + j*stride_]);
      for (i = j+1; i < numRowCols_; i++) {
        sum += ScalarTraits<ScalarType>::magnitude(2.0*values_[i+j*stride_]*values_[i+j*stride_]);
      }
    }
  }       
  anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::squareroot(sum));
  return(anorm);
}
  
//----------------------------------------------------------------------------------------------------
//   Comparison methods 
//----------------------------------------------------------------------------------------------------
  
template<typename OrdinalType, typename ScalarType>
bool SerialSymDenseMatrix<OrdinalType, ScalarType>::operator== (const SerialSymDenseMatrix<OrdinalType, ScalarType> &Operand) const
{
  bool result = 1;
  if((numRowCols_ != Operand.numRowCols_)) {
    result = 0;
  }
  else {
    OrdinalType i, j;
    for(i = 0; i < numRowCols_; i++) {
      for(j = 0; j < numRowCols_; j++) {
	if((*this)(i, j) != Operand(i, j)) {
	  return 0;
	}
      }
    }
  }
  return result;
}

template<typename OrdinalType, typename ScalarType>
bool SerialSymDenseMatrix<OrdinalType, ScalarType>::operator!= (const SerialSymDenseMatrix<OrdinalType, ScalarType> &Operand) const
{
  return !((*this) == Operand);
}

//----------------------------------------------------------------------------------------------------
//   Multiplication method 
//----------------------------------------------------------------------------------------------------  

template<typename OrdinalType, typename ScalarType>
void SerialSymDenseMatrix<OrdinalType, ScalarType>::scale( const ScalarType alpha )
{
  OrdinalType i, j;
  ScalarType* ptr;
  
  if (upper_) {
    for (j=0; j<numRowCols_; j++) {
      ptr = values_ + j*stride_;
      for (i=0; i<= j; i++) { *ptr = alpha * (*ptr); ptr++; }
    }
  }
  else {
    for (j=0; j<numRowCols_; j++) {
      ptr = values_ + j*stride_ + j;
      for (i=j; i<numRowCols_; i++) { *ptr = alpha * (*ptr); ptr++; }
    }
  }
}

/*
template<typename OrdinalType, typename ScalarType>
int SerialSymDenseMatrix<OrdinalType, ScalarType>::scale( const SerialSymDenseMatrix<OrdinalType,ScalarType>& A )
{
  OrdinalType i, j;
  ScalarType* ptr;
    
  // Check for compatible dimensions
  if ((numRowCols_ != A.numRowCols_)) {
    TEUCHOS_CHK_ERR(-1); // Return error
  }    

  if (upper_) {
    for (j=0; j<numRowCols_; j++) {
      ptr = values_ + j*stride_;
      for (i=0; i<= j; i++) { *ptr = A(i,j) * (*ptr); ptr++; }
    }
  }
  else {
    for (j=0; j<numRowCols_; j++) {
      ptr = values_ + j*stride_;
      for (i=j; i<numRowCols_; i++) { *ptr = A(i,j) * (*ptr); ptr++; }
    }
  }

  return(0);
}
*/

template<typename OrdinalType, typename ScalarType>
void SerialSymDenseMatrix<OrdinalType, ScalarType>::print(std::ostream& os) const
{
  os << std::endl;
  if(valuesCopied_)
    os << "Values_copied : yes" << std::endl;
  else
    os << "Values_copied : no" << std::endl;
  os << "Rows / Columns : " << numRowCols_ << std::endl;
  os << "LDA : " << stride_ << std::endl;
  if (upper_) 
    os << "Storage: Upper" << std::endl;
  else
    os << "Storage: Lower" << std::endl;
  if(numRowCols_ == 0) {
    os << "(matrix is empty, no values to display)" << std::endl;
  } else {
    for(OrdinalType i = 0; i < numRowCols_; i++) {
      for(OrdinalType j = 0; j < numRowCols_; j++){
        os << (*this)(i,j) << " ";
      }
      os << std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------------------
//   Protected methods 
//----------------------------------------------------------------------------------------------------  

template<typename OrdinalType, typename ScalarType>
inline void SerialSymDenseMatrix<OrdinalType, ScalarType>::checkIndex( OrdinalType rowIndex, OrdinalType colIndex ) const {
  TEUCHOS_TEST_FOR_EXCEPTION(rowIndex < 0 || rowIndex >= numRowCols_, std::out_of_range,
    "SerialSymDenseMatrix<T>::checkIndex: "
    "Row index " << rowIndex << " out of range [0, "<< numRowCols_ << ")");
  TEUCHOS_TEST_FOR_EXCEPTION(colIndex < 0 || colIndex >= numRowCols_, std::out_of_range,
    "SerialSymDenseMatrix<T>::checkIndex: "
    "Col index " << colIndex << " out of range [0, "<< numRowCols_ << ")");
}

template<typename OrdinalType, typename ScalarType>
void SerialSymDenseMatrix<OrdinalType, ScalarType>::deleteArrays(void)
{
  if (valuesCopied_)
  {
    delete [] values_;
    values_ = 0;
    valuesCopied_ = false;
  }
}

template<typename OrdinalType, typename ScalarType>
void SerialSymDenseMatrix<OrdinalType, ScalarType>::copyMat(
							    bool inputUpper, ScalarType* inputMatrix, 
							    OrdinalType inputStride, OrdinalType numRowCols_in, 
							    bool outputUpper, ScalarType* outputMatrix, 
							    OrdinalType outputStride, OrdinalType startRowCol, 
							    ScalarType alpha
							    )
{
  OrdinalType i, j;
  ScalarType* ptr1 = 0;
  ScalarType* ptr2 = 0;

  for (j = 0; j < numRowCols_in; j++) {
    if (inputUpper == true) {
      // The input matrix is upper triangular, start at the beginning of each column.
      ptr2 = inputMatrix + (j + startRowCol) * inputStride + startRowCol;
      if (outputUpper == true) {
	// The output matrix matches the same pattern as the input matrix.
	ptr1 = outputMatrix + j*outputStride;
	if (alpha != Teuchos::ScalarTraits<ScalarType>::zero() ) {
	  for(i = 0; i <= j; i++) {
	    *ptr1++ += alpha*(*ptr2++);
	  }
	} else {
	  for(i = 0; i <= j; i++) {
	    *ptr1++ = *ptr2++;
	  }
	}
      }
      else {
	// The output matrix has the opposite pattern as the input matrix.
	// Copy down across rows of the output matrix, but down columns of the input matrix.
	ptr1 = outputMatrix + j;
	if (alpha != Teuchos::ScalarTraits<ScalarType>::zero() ) {
	  for(i = 0; i <= j; i++) {
	    *ptr1 += alpha*(*ptr2++);
	    ptr1 += outputStride;
	  }
	} else {
	  for(i = 0; i <= j; i++) {
	    *ptr1 = *ptr2++;
	    ptr1 += outputStride;
	  }
	}
      }
    }
    else {
      // The input matrix is lower triangular, start at the diagonal of each row.
      ptr2 = inputMatrix + (startRowCol+j) * inputStride + startRowCol + j;
      if (outputUpper == true) {
	// The output matrix has the opposite pattern as the input matrix.
	// Copy across rows of the output matrix, but down columns of the input matrix.
	ptr1 = outputMatrix + j*outputStride + j;
	if (alpha != Teuchos::ScalarTraits<ScalarType>::zero() ) {
	  for(i = j; i < numRowCols_in; i++) {
	    *ptr1 += alpha*(*ptr2++);
	    ptr1 += outputStride;
	  }
	} else {
	  for(i = j; i < numRowCols_in; i++) {
	    *ptr1 = *ptr2++;
	    ptr1 += outputStride;
	  }
	}
      }
      else {
	// The output matrix matches the same pattern as the input matrix.
	ptr1 = outputMatrix + j*outputStride + j;
	if (alpha != Teuchos::ScalarTraits<ScalarType>::zero() ) {
	  for(i = j; i < numRowCols_in; i++) {
	    *ptr1++ += alpha*(*ptr2++);
	  }
	} else {
	  for(i = j; i < numRowCols_in; i++) {
	    *ptr1++ = *ptr2++;
	  }
	}
      }
    }
  }
}
  
template<typename OrdinalType, typename ScalarType>
void SerialSymDenseMatrix<OrdinalType, ScalarType>::copyUPLOMat(
								bool inputUpper, ScalarType* inputMatrix, 
								OrdinalType inputStride, OrdinalType inputRows 
								)
{
  OrdinalType i, j;
  ScalarType * ptr1 = 0;
  ScalarType * ptr2 = 0;

  if (inputUpper) {
    for (j=1; j<inputRows; j++) {
      ptr1 = inputMatrix + j;
      ptr2 = inputMatrix + j*inputStride;
      for (i=0; i<j; i++) {
	*ptr1 = *ptr2++;
	ptr1+=inputStride;
      }
    }
  }
  else {
    for (i=1; i<inputRows; i++) {
      ptr1 = inputMatrix + i;
      ptr2 = inputMatrix + i*inputStride;
      for (j=0; j<i; j++) {
	*ptr2++ = *ptr1;
	ptr1+=inputStride;
      }
    }
  }
}
  
} // namespace Teuchos

#endif /* _TEUCHOS_SERIALSYMDENSEMATRIX_HPP_ */
