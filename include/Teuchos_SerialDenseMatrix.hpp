// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _TEUCHOS_SERIALDENSEMATRIX_HPP_
#define _TEUCHOS_SERIALDENSEMATRIX_HPP_
/*! \file Teuchos_SerialDenseMatrix.hpp
  \brief Templated serial dense matrix class
*/

#include "Teuchos_CompObject.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"
#include <cstddef>
#include <utility>

/*!     \class Teuchos::SerialDenseMatrix
        \brief This class creates and provides basic support for dense rectangular matrix of templated type.
*/
/** \example DenseMatrix/cxx_main.cpp
    This is an example of how to use the Teuchos::SerialDenseMatrix class.
*/


namespace Teuchos {

template<typename OrdinalType, typename ScalarType>
class SerialDenseMatrix : public CompObject, public BLAS<OrdinalType, ScalarType>
{
public:

  //! Typedef for ordinal type
  typedef OrdinalType ordinalType;
  //! Typedef for scalar type
  typedef ScalarType scalarType;

  //! @name Constructor/Destructor methods.
  //@{

  //! Default Constructor
  /*! Creates a empty matrix of no dimension.  The Shaping methods should be used to size this matrix.
    Values of this matrix should be set using the [], (), or = operators.
  */
  SerialDenseMatrix() = default;

  //! Shaped Constructor
  /*!
    \param numRows - Number of rows in matrix.
    \param numCols - Number of columns in matrix.
    \param zeroOut - Initializes values to 0 if true (default)

    Creates a shaped matrix with \c numRows rows and \c numCols cols.  All values are initialized to 0 when \c zeroOut is true.
    Values of this matrix should be set using the [] or the () operators.
  */
  SerialDenseMatrix(OrdinalType numRows, OrdinalType numCols, bool zeroOut = true);

  //! Shaped Constructor with Values
  /*!
    \param CV - Enumerated type set to Teuchos::Copy or Teuchos::View.
    \param values - Pointer to an array of ScalarType.  The first column starts at \c values,
                the second at \c values+stride, etc.
    \param stride - The stride between the columns of the matrix in memory.
    \param numRows - Number of rows in matrix.
    \param numCols - Number of columns in matrix.
  */
  SerialDenseMatrix(DataAccess CV, ScalarType* values, OrdinalType stride, OrdinalType numRows, OrdinalType numCols);

  //! Copy Constructor
  /*! \note A deep copy of the \c Source transposed can be obtained if \c trans=Teuchos::TRANS, \c else
    a non-transposed copy of \c Source is made.  There is no storage of the transpose state of the matrix
    within the SerialDenseMatrix class, so this information will not propogate to any operation performed
    on a matrix that has been copy constructed in transpose.
  */
  SerialDenseMatrix(const SerialDenseMatrix<OrdinalType, ScalarType> &Source, ETransp trans = Teuchos::NO_TRANS);

  //! Copy Constructor
  /*! \note Only a non-transposed deep copy or view of \c Source is made with this copy constructor.
  */
  SerialDenseMatrix(DataAccess CV, const SerialDenseMatrix<OrdinalType, ScalarType> &Source);

  //! Submatrix Copy Constructor
  /*!
    \param CV - Enumerated type set to Teuchos::Copy or Teuchos::View.
    \param Source - Reference to another dense matrix from which values are to be copied.
    \param numRows - The number of rows in this matrix.
    \param numCols - The number of columns in this matrix.
    \param startRow - The row of \c Source from which the submatrix copy should start.
    \param startCol - The column of \c Source from which the submatrix copy should start.

    Creates a shaped matrix with \c numRows rows and \c numCols columns, which is a submatrix of \c Source.
    If \c startRow and \c startCol are not given, then the submatrix is the leading submatrix of \c Source.
    Otherwise, the (1,1) entry in the copied matrix is the (\c startRow, \c startCol) entry of \c Source.
  */
  SerialDenseMatrix(DataAccess CV, const SerialDenseMatrix<OrdinalType, ScalarType> &Source, OrdinalType numRows, OrdinalType numCols, OrdinalType startRow=0, OrdinalType startCol=0);

  //! Destructor
  virtual ~SerialDenseMatrix();
  //@}

  //! @name Shaping methods.
  //@{
  //! Shape method for changing the size of a SerialDenseMatrix, initializing entries to zero.
  /*!
    \param numRows - The number of rows in this matrix.
    \param numCols - The number of columns in this matrix.

    This method allows the user to define the dimensions of a SerialDenseMatrix at any point.  This method
    can be called at any point after construction.  Any values previously in this object will be destroyed
    and the resized matrix starts of with all zero values.

    \return Integer error code, set to 0 if successful.
  */
  int shape(OrdinalType numRows, OrdinalType numCols);

  //! Same as <tt>shape()</tt> except leaves uninitialized.
  int shapeUninitialized(OrdinalType numRows, OrdinalType numCols);

  //! Reshaping method for changing the size of a SerialDenseMatrix, keeping the entries.
  /*!
    \param numRows - The number of rows in this matrix.
    \param numCols - The number of columns in this matrix.

    This method allows the user to redefine the dimensions of a SerialDenseMatrix at any point.  This method
    can be called at any point after construction.  Any values previously in this object will be copied into
    the reshaped matrix.

    \return Integer error code, set 0 if successful.
  */
  int reshape(OrdinalType numRows, OrdinalType numCols);


  //@}

  //! @name Set methods.
  //@{

  //! Copies values from one matrix to another.
  /*!
    The operator= copies the values from one existing SerialDenseMatrix to another.
    If \c Source is a view (i.e. CV = Teuchos::View), then this method will
    return a view.  Otherwise, it will return a copy of \c Source.  \e this object
    will be resized if it is not large enough to copy \c Source into.
  */
  SerialDenseMatrix<OrdinalType, ScalarType>& operator= (const SerialDenseMatrix<OrdinalType, ScalarType>& Source);

  //! Copies values from one matrix to another.
  /*!
    Copies the values from one existing SerialDenseMatrix to another
    if the dimension of both matrices are the same.  If not, \e this matrix
    will be returned unchanged.
  */
  SerialDenseMatrix<OrdinalType, ScalarType>& assign (const SerialDenseMatrix<OrdinalType, ScalarType>& Source);

  //! Set all values in the matrix to a constant value.
  /*!
    \param value - Value to use;
  */
  SerialDenseMatrix<OrdinalType, ScalarType>& operator= (const ScalarType value) { putScalar(value); return(*this); }

  //! Set all values in the matrix to a constant value.
  /*!
    \param value - Value to use; zero if none specified.
    \return Integer error code, set to 0 if successful.
  */
  int putScalar( const ScalarType value = Teuchos::ScalarTraits<ScalarType>::zero() );

  //! Swap values between this matrix and incoming matrix.
  /*!
    Swaps pointers and associated state without copying the matrix data.
  */
  void swap (SerialDenseMatrix<OrdinalType, ScalarType> &B);

  //! Set all values in the matrix to be random numbers.
  int random();

  //@}

  //! @name Accessor methods.
  //@{

  //! Element access method (non-const).
  /*! Returns the element in the ith row and jth column if A(i,j) is specified, the
    expression A[j][i] will return the same element.

        \return Element from the specified \c rowIndex row and \c colIndex column.
    \warning The validity of \c rowIndex and \c colIndex will only be checked if Teuchos is
    configured with --enable-teuchos-abc.
  */
  ScalarType& operator () (OrdinalType rowIndex, OrdinalType colIndex);

  //! Element access method (const).
  /*! Returns the element in the ith row and jth column if A(i,j) is specified, the expression
    A[j][i] will return the same element.

        \return Element from the specified \c rowIndex row and \c colIndex column.
    \warning The validity of \c rowIndex and \c colIndex will only be checked if Teuchos is
    configured with --enable-teuchos-abc.
  */
  const ScalarType& operator () (OrdinalType rowIndex, OrdinalType colIndex) const;

  //! Column access method (non-const).
  /*! Returns the pointer to the ScalarType array at the jth column if A[j] is specified, the expression
    A[j][i] will return the same element as A(i,j).

    \return Pointer to the ScalarType array at the \c colIndex column ( \c values_+colIndex*stride_ ).
  */
  ScalarType* operator [] (OrdinalType colIndex);

  //! Column access method (const).
  /*! Returns the pointer to the ScalarType array at the jth column if A[j] is specified, the expression
    A[j][i] will return the same element as A(i,j).

    \return Pointer to the ScalarType array at the \c colIndex column ( \c values_+colIndex*stride_ ).
  */
  const ScalarType* operator [] (OrdinalType colIndex) const;

  //! Data array access method.
  /*! \return Pointer to the ScalarType data array contained in the object. */
  ScalarType* values() const { return values_; }

  //@}

  //! @name Mathematical methods.
  //@{

  //! Add another matrix to \e this matrix.
  /*! Add \c Source to \e this if the dimension of both matrices are the same.  If not, \e this matrix
    will be returned unchanged.
  */
  SerialDenseMatrix<OrdinalType, ScalarType>& operator+= (const SerialDenseMatrix<OrdinalType, ScalarType>& Source);

  //! Subtract another matrix from \e this matrix.
  /*! Subtract \c Source from \e this if the dimension of both matrices are the same.  If not, \e this matrix
    will be returned unchanged.
  */
  SerialDenseMatrix<OrdinalType, ScalarType>& operator-= (const SerialDenseMatrix<OrdinalType, ScalarType>& Source);

  //! Scale \c this matrix by \c alpha; \c *this = \c alpha*\c *this.
  /*!
    \param alpha Scalar to multiply \e this by.
  */
  SerialDenseMatrix<OrdinalType, ScalarType>& operator*= (const ScalarType alpha);

  //! Scale \c this matrix by \c alpha; \c *this = \c alpha*\c *this.
  /*!
    \param alpha Scalar to multiply \e this by.
    \return Integer error code, set to 0 if successful.
  */
  int scale ( const ScalarType alpha );

  //! Point-wise scale \c this matrix by \c A; i.e. *this(i,j) *= A(i,j)
  /*! The values of \c *this matrix will be point-wise scaled by the values in A.
    If A and \c this matrix are not the same dimension \c this will be returned unchanged.

    \param B Teuchos::SerialDenseMatrix used to perform element-wise scaling of \e this.
    \return Integer error code, set to 0 if successful.
  */
  int scale ( const SerialDenseMatrix<OrdinalType, ScalarType>& A );

  //! Multiply \c A * \c B and add them to \e this; \e this = \c beta * \e this + \c alpha*A*B.
  /*!
    \param transa - Use the transpose of \c A if transa = Teuchos::TRANS, else don't use the
    transpose if transa = Teuchos::NO_TRANS.
    \param transb - Use the transpose of \c B if transb = Teuchos::TRANS, else don't use the
    transpose if transb = Teuchos::NO_TRANS.
    \param alpha - The scaling factor for \c A * \c B.
    \param A - SerialDenseMatrix
    \param B - SerialDenseMatrix
    \param beta - The scaling factor for \e this.

    If the matrices \c A and \c B are not of the right dimension, consistent with \e this, then \e this
    matrix will not be altered and -1 will be returned.
    \return Integer error code, set to 0 if successful.
  */
  int multiply (ETransp transa, ETransp transb, ScalarType alpha, const SerialDenseMatrix<OrdinalType, ScalarType> &A, const SerialDenseMatrix<OrdinalType, ScalarType> &B, ScalarType beta);

  //! Multiply \c A and \c B and add them to \e this; \e this = \c beta * \e this + \c alpha*A*B or \e this = \c beta * \e this + \c alpha*B*A.
  /*!
    \param sideA - Which side is A on for the multiplication to B, A*B (Teuchos::LEFT_SIDE) or B*A (Teuchos::RIGHT_SIDE).
    \param alpha - The scaling factor for \c A * \c B, or \c B * \c A.
    \param A - SerialSymDenseMatrix (a serial SPD dense matrix)
    \param B - SerialDenseMatrix (a serial dense matrix)
    \param beta - The scaling factor for \e this.

    If the matrices \c A and \c B are not of the right dimension, consistent with \e this, then \e this
    matrix will not be altered and -1 will be returned.
    \return Integer error code, set to 0 if successful.
  */
  int multiply (ESide sideA, ScalarType alpha, const SerialSymDenseMatrix<OrdinalType, ScalarType> &A, const SerialDenseMatrix<OrdinalType, ScalarType> &B, ScalarType beta);

  //@}

  //! @name Comparison methods.
  //@{

  //! Equality of two matrices.
  /*! \return True if \e this matrix and \c Operand are of the same shape (rows and columns) and have
    the same entries, else False will be returned.
  */
  bool operator== (const SerialDenseMatrix<OrdinalType, ScalarType> &Operand) const;

  //! Inequality of two matrices.
  /*! \return True if \e this matrix and \c Operand of not of the same shape (rows and columns) or don't
    have the same entries, else False will be returned.
  */
  bool operator!= (const SerialDenseMatrix<OrdinalType, ScalarType> &Operand) const;

  //@}

  //! @name Attribute methods.
  //@{

  //! Returns the row dimension of this matrix.
  OrdinalType numRows() const { return(numRows_); }

  //! Returns the column dimension of this matrix.
  OrdinalType numCols() const { return(numCols_); }

  //! Returns the stride between the columns of this matrix in memory.
  OrdinalType stride() const { return(stride_); }

  //! Returns whether this matrix is empty.
  bool empty() const { return(numRows_ == 0 || numCols_ == 0); }
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
  //! Print method.  Defines the behavior of the std::ostream << operator
  virtual std::ostream& print(std::ostream& os) const;

  //@}
protected:
  void copyMat(ScalarType* inputMatrix, OrdinalType strideInput,
    OrdinalType numRows, OrdinalType numCols, ScalarType* outputMatrix,
    OrdinalType strideOutput, OrdinalType startRow, OrdinalType startCol,
    ScalarType alpha = ScalarTraits<ScalarType>::zero() );
  void deleteArrays();
  void checkIndex( OrdinalType rowIndex, OrdinalType colIndex = 0 ) const;

  static ScalarType*
  allocateValues(const OrdinalType numRows,
                 const OrdinalType numCols)
  {
    const size_t size = size_t(numRows) * size_t(numCols);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvla"
    return new ScalarType[size];
#pragma GCC diagnostic pop
  }

  OrdinalType numRows_ = 0;
  OrdinalType numCols_ = 0;
  OrdinalType stride_ = 0;
  bool valuesCopied_ = false;
  ScalarType* values_ = nullptr;
}; // class Teuchos_SerialDenseMatrix

//----------------------------------------------------------------------------------------------------
//  Constructors and Destructor
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
SerialDenseMatrix<OrdinalType, ScalarType>::SerialDenseMatrix(
  OrdinalType numRows_in, OrdinalType numCols_in, bool zeroOut
  )
  : numRows_(numRows_in),
    numCols_(numCols_in),
    stride_(numRows_in),
    valuesCopied_(true),
    values_(allocateValues(numRows_in, numCols_in))
{
  if (zeroOut) {
    putScalar();
  }
}

template<typename OrdinalType, typename ScalarType>
SerialDenseMatrix<OrdinalType, ScalarType>::SerialDenseMatrix(
  DataAccess CV, ScalarType* values_in, OrdinalType stride_in, OrdinalType numRows_in,
  OrdinalType numCols_in
  )
  : numRows_(numRows_in),
    numCols_(numCols_in),
    stride_(stride_in),
    valuesCopied_(false),
    values_(values_in)
{
  if(CV == Copy)
  {
    stride_ = numRows_;
    values_ = allocateValues(stride_, numCols_);
    copyMat(values_in, stride_in, numRows_, numCols_, values_, stride_, 0, 0);
    valuesCopied_ = true;
  }
}

template<typename OrdinalType, typename ScalarType>
SerialDenseMatrix<OrdinalType, ScalarType>::SerialDenseMatrix(const SerialDenseMatrix<OrdinalType, ScalarType> &Source, ETransp trans)
  : valuesCopied_(true)
{
  if ( trans == Teuchos::NO_TRANS )
  {
    numRows_ = Source.numRows_;
    numCols_ = Source.numCols_;

    if (!Source.valuesCopied_)
    {
      stride_ = Source.stride_;
      values_ = Source.values_;
      valuesCopied_ = false;
    }
    else
    {
      stride_ = numRows_;
      if(stride_ > 0 && numCols_ > 0) {
        values_ = allocateValues(stride_, numCols_);
        copyMat(Source.values_, Source.stride_, numRows_, numCols_, values_, stride_, 0, 0);
      }
      else {
        numRows_ = 0; numCols_ = 0; stride_ = 0;
        valuesCopied_ = false;
      }
    }
  }
  else if ( trans == Teuchos::CONJ_TRANS && ScalarTraits<ScalarType>::isComplex )
  {
    numRows_ = Source.numCols_;
    numCols_ = Source.numRows_;
    stride_ = numRows_;
    values_ = allocateValues(stride_, numCols_);
    for (OrdinalType j=0; j<numCols_; j++) {
      for (OrdinalType i=0; i<numRows_; i++) {
        values_[j*stride_ + i] = Teuchos::ScalarTraits<ScalarType>::conjugate(Source.values_[i*Source.stride_ + j]);
      }
    }
  }
  else
  {
    numRows_ = Source.numCols_;
    numCols_ = Source.numRows_;
    stride_ = numRows_;
    values_ = allocateValues(stride_, numCols_);
    for (OrdinalType j=0; j<numCols_; j++) {
      for (OrdinalType i=0; i<numRows_; i++) {
        values_[j*stride_ + i] = Source.values_[i*Source.stride_ + j];
      }
    }
  }
}


template<typename OrdinalType, typename ScalarType>
SerialDenseMatrix<OrdinalType, ScalarType>::SerialDenseMatrix(
  DataAccess CV, const SerialDenseMatrix<OrdinalType, ScalarType> &Source
  )
  : numRows_(Source.numRows_), numCols_(Source.numCols_), stride_(Source.stride_),
    valuesCopied_(false), values_(Source.values_)
{
  if(CV == Copy)
  {
    stride_ = numRows_;
    values_ = allocateValues(stride_, numCols_);
    copyMat(Source.values_, Source.stride_, numRows_, numCols_, values_, stride_, 0, 0);
    valuesCopied_ = true;
  }
}


template<typename OrdinalType, typename ScalarType>
SerialDenseMatrix<OrdinalType, ScalarType>::SerialDenseMatrix(
  DataAccess CV, const SerialDenseMatrix<OrdinalType, ScalarType> &Source,
  OrdinalType numRows_in, OrdinalType numCols_in, OrdinalType startRow,
  OrdinalType startCol
  )
  : CompObject(), numRows_(numRows_in), numCols_(numCols_in), stride_(Source.stride_),
    valuesCopied_(false), values_(Source.values_)
{
  if(CV == Copy)
  {
    stride_ = numRows_in;
    values_ = allocateValues(stride_, numCols_in);
    copyMat(Source.values_, Source.stride_, numRows_in, numCols_in, values_, stride_, startRow, startCol);
    valuesCopied_ = true;
  }
  else // CV == View
  {
    values_ = values_ + (stride_ * startCol) + startRow;
  }
}

template<typename OrdinalType, typename ScalarType>
SerialDenseMatrix<OrdinalType, ScalarType>::~SerialDenseMatrix()
{
  deleteArrays();
}

//----------------------------------------------------------------------------------------------------
//  Shape methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
int SerialDenseMatrix<OrdinalType, ScalarType>::shape(
  OrdinalType numRows_in, OrdinalType numCols_in
  )
{
  deleteArrays(); // Get rid of anything that might be already allocated
  numRows_ = numRows_in;
  numCols_ = numCols_in;
  stride_ = numRows_;
  values_ = allocateValues(stride_, numCols_);
  putScalar();
  valuesCopied_ = true;
  return(0);
}

template<typename OrdinalType, typename ScalarType>
int SerialDenseMatrix<OrdinalType, ScalarType>::shapeUninitialized(
  OrdinalType numRows_in, OrdinalType numCols_in
  )
{
  deleteArrays(); // Get rid of anything that might be already allocated
  numRows_ = numRows_in;
  numCols_ = numCols_in;
  stride_ = numRows_;
  values_ = allocateValues(stride_, numCols_);
  valuesCopied_ = true;
  return(0);
}

template<typename OrdinalType, typename ScalarType>
int SerialDenseMatrix<OrdinalType, ScalarType>::reshape(
  OrdinalType numRows_in, OrdinalType numCols_in
  )
{
  // Allocate space for new matrix
  ScalarType* values_tmp = allocateValues(numRows_in, numCols_in);
  ScalarType zero = ScalarTraits<ScalarType>::zero();
  for(OrdinalType k = 0; k < numRows_in * numCols_in; k++)
  {
    values_tmp[k] = zero;
  }
  OrdinalType numRows_tmp = TEUCHOS_MIN(numRows_, numRows_in);
  OrdinalType numCols_tmp = TEUCHOS_MIN(numCols_, numCols_in);
  if(values_ != 0)
  {
    copyMat(values_, stride_, numRows_tmp, numCols_tmp, values_tmp,
      numRows_in, 0, 0); // Copy principal submatrix of A to new A
  }
  deleteArrays(); // Get rid of anything that might be already allocated
  numRows_ = numRows_in;
  numCols_ = numCols_in;
  stride_ = numRows_;
  values_ = values_tmp; // Set pointer to new A
  valuesCopied_ = true;
  return(0);
}

//----------------------------------------------------------------------------------------------------
//   Set methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
int SerialDenseMatrix<OrdinalType, ScalarType>::putScalar( const ScalarType value_in )
{
  // Set each value of the dense matrix to "value".
  for(OrdinalType j = 0; j < numCols_; j++)
  {
    for(OrdinalType i = 0; i < numRows_; i++)
          {
            values_[i + j*stride_] = value_in;
          }
  }
  return 0;
}

template<typename OrdinalType, typename ScalarType> void
SerialDenseMatrix<OrdinalType, ScalarType>::swap(
  SerialDenseMatrix<OrdinalType, ScalarType> &B)
{
  // Notes:
  // > DefaultBLASImpl::SWAP() uses a deep copy. This fn uses a pointer swap.
  // > this fn covers both Vector and Matrix, such that some care must be
  //   employed to not swap across types (creating a Vector with non-unitary
  //   numCols_)
  // > Inherited data that is not currently swapped (since inactive/deprecated):
  //   >> Teuchos::CompObject:
  //        Flops *flopCounter_ [Note: all SerialDenseMatrix ctors initialize a
  //        NULL flop-counter using CompObject(), such that any flop increments
  //        that are computed are not accumulated.]
  //   >> Teuchos::Object: (now removed from inheritance list)
  //        static int tracebackMode (no swap for statics)
  //        std::string label_ (has been reported as a cause of memory overhead)

  std::swap(values_ ,      B.values_);
  std::swap(numRows_,      B.numRows_);
  std::swap(numCols_,      B.numCols_);
  std::swap(stride_,       B.stride_);
  std::swap(valuesCopied_, B.valuesCopied_);
}

template<typename OrdinalType, typename ScalarType>
int SerialDenseMatrix<OrdinalType, ScalarType>::random()
{
  // Set each value of the dense matrix to a random value.
  for(OrdinalType j = 0; j < numCols_; j++)
  {
    for(OrdinalType i = 0; i < numRows_; i++)
          {
            values_[i + j*stride_] = ScalarTraits<ScalarType>::random();
          }
  }
  return 0;
}

template<typename OrdinalType, typename ScalarType>
SerialDenseMatrix<OrdinalType,ScalarType>&
SerialDenseMatrix<OrdinalType, ScalarType>::operator=(
  const SerialDenseMatrix<OrdinalType,ScalarType>& Source
  )
{
  if(this == &Source)
    return(*this); // Special case of source same as target
  if((!valuesCopied_) && (!Source.valuesCopied_) && (values_ == Source.values_))
    return(*this); // Special case of both are views to same data.

  // If the source is a view then we will return a view, else we will return a copy.
  if (!Source.valuesCopied_) {
    if(valuesCopied_)       {
      // Clean up stored data if this was previously a copy.
      deleteArrays();
    }
    numRows_ = Source.numRows_;
    numCols_ = Source.numCols_;
    stride_ = Source.stride_;
    values_ = Source.values_;
  }
  else {
    // If we were a view, we will now be a copy.
    if(!valuesCopied_) {
      numRows_ = Source.numRows_;
      numCols_ = Source.numCols_;
      stride_ = Source.numRows_;
      if(stride_ > 0 && numCols_ > 0) {
        values_ = allocateValues(stride_, numCols_);
        valuesCopied_ = true;
      }
      else {
        values_ = 0;
      }
    }
    // If we were a copy, we will stay a copy.
    else {
      if((Source.numRows_ <= stride_) && (Source.numCols_ == numCols_)) { // we don't need to reallocate
        numRows_ = Source.numRows_;
        numCols_ = Source.numCols_;
      }
      else { // we need to allocate more space (or less space)
        deleteArrays();
        numRows_ = Source.numRows_;
        numCols_ = Source.numCols_;
        stride_ = Source.numRows_;
        if(stride_ > 0 && numCols_ > 0) {
          values_ = allocateValues(stride_, numCols_);
          valuesCopied_ = true;
        }
      }
    }
    copyMat(Source.values_, Source.stride_, numRows_, numCols_, values_, stride_, 0, 0);
  }
  return(*this);
}

template<typename OrdinalType, typename ScalarType>
SerialDenseMatrix<OrdinalType, ScalarType>& SerialDenseMatrix<OrdinalType, ScalarType>::operator+= (const SerialDenseMatrix<OrdinalType,ScalarType>& Source )
{
  // Check for compatible dimensions
  if ((numRows_ != Source.numRows_) || (numCols_ != Source.numCols_))
  {
    TEUCHOS_CHK_REF(*this); // Return *this without altering it.
  }
  copyMat(Source.values_, Source.stride_, numRows_, numCols_, values_, stride_, 0, 0, ScalarTraits<ScalarType>::one());
  return(*this);
}

template<typename OrdinalType, typename ScalarType>
SerialDenseMatrix<OrdinalType, ScalarType>& SerialDenseMatrix<OrdinalType, ScalarType>::operator-= (const SerialDenseMatrix<OrdinalType,ScalarType>& Source )
{
  // Check for compatible dimensions
  if ((numRows_ != Source.numRows_) || (numCols_ != Source.numCols_))
  {
    TEUCHOS_CHK_REF(*this); // Return *this without altering it.
  }
  copyMat(Source.values_, Source.stride_, numRows_, numCols_, values_, stride_, 0, 0, -ScalarTraits<ScalarType>::one());
  return(*this);
}

template<typename OrdinalType, typename ScalarType>
SerialDenseMatrix<OrdinalType,ScalarType>& SerialDenseMatrix<OrdinalType, ScalarType>::assign (const SerialDenseMatrix<OrdinalType,ScalarType>& Source) {
  if(this == &Source)
    return(*this); // Special case of source same as target
  if((!valuesCopied_) && (!Source.valuesCopied_) && (values_ == Source.values_))
    return(*this); // Special case of both are views to same data.

  // Check for compatible dimensions
  if ((numRows_ != Source.numRows_) || (numCols_ != Source.numCols_))
  {
    TEUCHOS_CHK_REF(*this); // Return *this without altering it.
  }
  copyMat(Source.values_, Source.stride_, numRows_, numCols_, values_, stride_, 0, 0);
  return(*this);
}

//----------------------------------------------------------------------------------------------------
//   Accessor methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
inline ScalarType& SerialDenseMatrix<OrdinalType, ScalarType>::operator () (OrdinalType rowIndex, OrdinalType colIndex)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  checkIndex( rowIndex, colIndex );
#endif
  return(values_[colIndex * stride_ + rowIndex]);
}

template<typename OrdinalType, typename ScalarType>
inline const ScalarType& SerialDenseMatrix<OrdinalType, ScalarType>::operator () (OrdinalType rowIndex, OrdinalType colIndex) const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  checkIndex( rowIndex, colIndex );
#endif
  return(values_[colIndex * stride_ + rowIndex]);
}

template<typename OrdinalType, typename ScalarType>
inline const ScalarType* SerialDenseMatrix<OrdinalType, ScalarType>::operator [] (OrdinalType colIndex) const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  checkIndex( 0, colIndex );
#endif
  return(values_ + colIndex * stride_);
}

template<typename OrdinalType, typename ScalarType>
inline ScalarType* SerialDenseMatrix<OrdinalType, ScalarType>::operator [] (OrdinalType colIndex)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  checkIndex( 0, colIndex );
#endif
  return(values_ + colIndex * stride_);
}

//----------------------------------------------------------------------------------------------------
//   Norm methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
typename ScalarTraits<ScalarType>::magnitudeType SerialDenseMatrix<OrdinalType, ScalarType>::normOne() const
{
  OrdinalType i, j;
  typename ScalarTraits<ScalarType>::magnitudeType anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
  typename ScalarTraits<ScalarType>::magnitudeType absSum = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
  ScalarType* ptr;
  for(j = 0; j < numCols_; j++)
  {
    ScalarType sum = 0;
    ptr = values_ + j * stride_;
    for(i = 0; i < numRows_; i++)
          {
            sum += ScalarTraits<ScalarType>::magnitude(*ptr++);
          }
    absSum = ScalarTraits<ScalarType>::magnitude(sum);
    if(absSum > anorm)
          {
            anorm = absSum;
          }
  }
  // don't compute flop increment unless there is an accumulator
  if (flopCounter_!=0) updateFlops(numRows_ * numCols_);
  return(anorm);
}

template<typename OrdinalType, typename ScalarType>
typename ScalarTraits<ScalarType>::magnitudeType SerialDenseMatrix<OrdinalType, ScalarType>::normInf() const
{
  OrdinalType i, j;
  typename ScalarTraits<ScalarType>::magnitudeType sum, anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());

  for (i = 0; i < numRows_; i++) {
    sum = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
    for (j=0; j< numCols_; j++) {
      sum += ScalarTraits<ScalarType>::magnitude(*(values_+i+j*stride_));
    }
    anorm = TEUCHOS_MAX( anorm, sum );
  }
  // don't compute flop increment unless there is an accumulator
  if (flopCounter_!=0) updateFlops(numRows_ * numCols_);
  return(anorm);
}

template<typename OrdinalType, typename ScalarType>
typename ScalarTraits<ScalarType>::magnitudeType SerialDenseMatrix<OrdinalType, ScalarType>::normFrobenius() const
{
  OrdinalType i, j;
  typename ScalarTraits<ScalarType>::magnitudeType anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
  for (j = 0; j < numCols_; j++) {
    for (i = 0; i < numRows_; i++) {
      anorm += ScalarTraits<ScalarType>::magnitude(values_[i+j*stride_]*values_[i+j*stride_]);
    }
  }
  anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::squareroot(anorm));
  // don't compute flop increment unless there is an accumulator
  if (flopCounter_!=0) updateFlops(numRows_ * numCols_);
  return(anorm);
}

//----------------------------------------------------------------------------------------------------
//   Comparison methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
bool SerialDenseMatrix<OrdinalType, ScalarType>::operator== (const SerialDenseMatrix<OrdinalType, ScalarType> &Operand) const
{
  bool result = 1;
  if((numRows_ != Operand.numRows_) || (numCols_ != Operand.numCols_))
  {
    result = 0;
  }
  else
  {
    OrdinalType i, j;
    for(i = 0; i < numRows_; i++)
          {
            for(j = 0; j < numCols_; j++)
      {
        if((*this)(i, j) != Operand(i, j))
        {
          return 0;
        }
      }
          }
  }
  return result;
}

template<typename OrdinalType, typename ScalarType>
bool SerialDenseMatrix<OrdinalType, ScalarType>::operator!= (const SerialDenseMatrix<OrdinalType, ScalarType> &Operand) const
{
  return !((*this) == Operand);
}

//----------------------------------------------------------------------------------------------------
//   Multiplication method
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
SerialDenseMatrix<OrdinalType, ScalarType>& SerialDenseMatrix<OrdinalType, ScalarType>::operator*= (const ScalarType alpha )
{
  this->scale( alpha );
  return(*this);
}

template<typename OrdinalType, typename ScalarType>
int SerialDenseMatrix<OrdinalType, ScalarType>::scale( const ScalarType alpha )
{
  OrdinalType i, j;
  ScalarType* ptr;

  for (j=0; j<numCols_; j++) {
    ptr = values_ + j*stride_;
    for (i=0; i<numRows_; i++) { *ptr = alpha * (*ptr); ptr++; }
  }
  // don't compute flop increment unless there is an accumulator
  if (flopCounter_!=0) updateFlops( numRows_*numCols_ );
  return(0);
}

template<typename OrdinalType, typename ScalarType>
int SerialDenseMatrix<OrdinalType, ScalarType>::scale( const SerialDenseMatrix<OrdinalType,ScalarType>& A )
{
  OrdinalType i, j;
  ScalarType* ptr;

  // Check for compatible dimensions
  if ((numRows_ != A.numRows_) || (numCols_ != A.numCols_))
  {
    TEUCHOS_CHK_ERR(-1); // Return error
  }
  for (j=0; j<numCols_; j++) {
    ptr = values_ + j*stride_;
    for (i=0; i<numRows_; i++) { *ptr = A(i,j) * (*ptr); ptr++; }
  }
  // don't compute flop increment unless there is an accumulator
  if (flopCounter_!=0) updateFlops( numRows_*numCols_ );
  return(0);
}

template<typename OrdinalType, typename ScalarType>
int  SerialDenseMatrix<OrdinalType, ScalarType>::multiply(ETransp transa, ETransp transb, ScalarType alpha, const SerialDenseMatrix<OrdinalType, ScalarType> &A, const SerialDenseMatrix<OrdinalType, ScalarType> &B, ScalarType beta)
{
  // Check for compatible dimensions
  OrdinalType A_nrows = (ETranspChar[transa]!='N') ? A.numCols() : A.numRows();
  OrdinalType A_ncols = (ETranspChar[transa]!='N') ? A.numRows() : A.numCols();
  OrdinalType B_nrows = (ETranspChar[transb]!='N') ? B.numCols() : B.numRows();
  OrdinalType B_ncols = (ETranspChar[transb]!='N') ? B.numRows() : B.numCols();
  if ((numRows_ != A_nrows) || (A_ncols != B_nrows) || (numCols_ != B_ncols))
  {
    TEUCHOS_CHK_ERR(-1); // Return error
  }
  // Call GEMM function
  this->GEMM(transa, transb, numRows_, numCols_, A_ncols, alpha, A.values(), A.stride(), B.values(), B.stride(), beta, values_, stride_);

  // don't compute flop increment unless there is an accumulator
  if (flopCounter_!=0) {
    double nflops = 2 * numRows_;
    nflops *= numCols_;
    nflops *= A_ncols;
    updateFlops(nflops);
  }
  return(0);
}

template<typename OrdinalType, typename ScalarType>
int SerialDenseMatrix<OrdinalType, ScalarType>::multiply (ESide sideA, ScalarType alpha, const SerialSymDenseMatrix<OrdinalType, ScalarType> &A, const SerialDenseMatrix<OrdinalType, ScalarType> &B, ScalarType beta)
{
  // Check for compatible dimensions
  OrdinalType A_nrows = A.numRows(), A_ncols = A.numCols();
  OrdinalType B_nrows = B.numRows(), B_ncols = B.numCols();

  if (ESideChar[sideA]=='L') {
    if ((numRows_ != A_nrows) || (A_ncols != B_nrows) || (numCols_ != B_ncols)) {
      TEUCHOS_CHK_ERR(-1); // Return error
    }
  } else {
    if ((numRows_ != B_nrows) || (B_ncols != A_nrows) || (numCols_ != A_ncols)) {
      TEUCHOS_CHK_ERR(-1); // Return error
    }
  }

  // Call SYMM function
  EUplo uplo = (A.upper() ? Teuchos::UPPER_TRI : Teuchos::LOWER_TRI);
  this->SYMM(sideA, uplo, numRows_, numCols_, alpha, A.values(), A.stride(), B.values(), B.stride(), beta, values_, stride_);

  // don't compute flop increment unless there is an accumulator
  if (flopCounter_!=0) {
    double nflops = 2 * numRows_;
    nflops *= numCols_;
    nflops *= A_ncols;
    updateFlops(nflops);
  }
  return(0);
}

template<typename OrdinalType, typename ScalarType>
std::ostream& SerialDenseMatrix<OrdinalType, ScalarType>::print(std::ostream& os) const
{
  os << std::endl;
  if(valuesCopied_)
    os << "Values_copied : yes" << std::endl;
  else
    os << "Values_copied : no" << std::endl;
  os << "Rows : " << numRows_ << std::endl;
  os << "Columns : " << numCols_ << std::endl;
  os << "LDA : " << stride_ << std::endl;
  if(numRows_ == 0 || numCols_ == 0) {
    os << "(matrix is empty, no values to display)" << std::endl;
  } else {
    for(OrdinalType i = 0; i < numRows_; i++) {
      for(OrdinalType j = 0; j < numCols_; j++){
        os << (*this)(i,j) << " ";
      }
      os << std::endl;
    }
  }
  return os;
}

//----------------------------------------------------------------------------------------------------
//   Protected methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
inline void SerialDenseMatrix<OrdinalType, ScalarType>::checkIndex( OrdinalType rowIndex, OrdinalType colIndex ) const {
  TEUCHOS_TEST_FOR_EXCEPTION(rowIndex < 0 || rowIndex >= numRows_, std::out_of_range,
    "SerialDenseMatrix<T>::checkIndex: "
    "Row index " << rowIndex << " out of range [0, "<< numRows_ << ")");
  TEUCHOS_TEST_FOR_EXCEPTION(colIndex < 0 || colIndex >= numCols_, std::out_of_range,
    "SerialDenseMatrix<T>::checkIndex: "
    "Col index " << colIndex << " out of range [0, "<< numCols_ << ")");
}

template<typename OrdinalType, typename ScalarType>
void SerialDenseMatrix<OrdinalType, ScalarType>::deleteArrays(void)
{
  if (valuesCopied_)
  {
    delete [] values_;
    values_ = 0;
    valuesCopied_ = false;
  }
}

template<typename OrdinalType, typename ScalarType>
void SerialDenseMatrix<OrdinalType, ScalarType>::copyMat(
  ScalarType* inputMatrix, OrdinalType strideInput, OrdinalType numRows_in,
  OrdinalType numCols_in, ScalarType* outputMatrix, OrdinalType strideOutput,
  OrdinalType startRow, OrdinalType startCol, ScalarType alpha
  )
{
  OrdinalType i, j;
  ScalarType* ptr1 = 0;
  ScalarType* ptr2 = 0;
  for(j = 0; j < numCols_in; j++) {
    ptr1 = outputMatrix + (j * strideOutput);
    ptr2 = inputMatrix + (j + startCol) * strideInput + startRow;
    if (alpha != Teuchos::ScalarTraits<ScalarType>::zero() ) {
      for(i = 0; i < numRows_in; i++)
            {
              *ptr1++ += alpha*(*ptr2++);
            }
    } else {
      for(i = 0; i < numRows_in; i++)
            {
              *ptr1++ = *ptr2++;
            }
    }
  }
}

/// \brief Ostream manipulator for SerialDenseMatrix
template<typename OrdinalType, typename ScalarType>
struct SerialDenseMatrixPrinter {
public:
  const SerialDenseMatrix<OrdinalType,ScalarType> &obj;
  SerialDenseMatrixPrinter(
        const SerialDenseMatrix<OrdinalType,ScalarType> &obj_in)
      : obj(obj_in) {}
};

/// \brief Output SerialDenseMatrix object through its stream manipulator.
template<typename OrdinalType, typename ScalarType>
std::ostream&
operator<<(std::ostream &out,
           const SerialDenseMatrixPrinter<OrdinalType,ScalarType> printer)
{
  printer.obj.print(out);
  return out;
}

/// \brief Return SerialDenseMatrix ostream manipulator Use as:
template<typename OrdinalType, typename ScalarType>
SerialDenseMatrixPrinter<OrdinalType,ScalarType>
printMat(const SerialDenseMatrix<OrdinalType,ScalarType> &obj)
{
  return SerialDenseMatrixPrinter<OrdinalType,ScalarType>(obj);
}


} // namespace Teuchos


#endif /* _TEUCHOS_SERIALDENSEMATRIX_HPP_ */
