// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _TEUCHOS_SERIALTRIDIMATRIX_HPP_
#define _TEUCHOS_SERIALTRIDIMATRIX_HPP_
/*! \file Teuchos_SerialTriDiMatrix.hpp
  \brief Templated serial TriDi matrix class
*/

#include "Teuchos_CompObject.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Assert.hpp"

/*!     \class Teuchos::SerialTriDiMatrix
        \brief This class creates and provides basic support for TriDi matrix of templated type.
*/
/** \example TriDiMatrix/cxx_main.cpp
    This is an example of how to use the Teuchos::SerialTriDiMatrix class.
*/


namespace Teuchos {

template<typename OrdinalType, typename ScalarType>
class SerialTriDiMatrix : public CompObject, public Object, public BLAS<OrdinalType, ScalarType > {
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
  SerialTriDiMatrix();

  //! Shaped Constructor
  /*!
    \param numRows - Number of rows in matrix.
    \param numCols - Number of columns in matrix.
    \param zeroOut - Initializes values to 0 if true (default)

    Creates a shaped matrix with \c numRows rows and \c numCols cols.  All values are initialized to 0 when \c zeroOut is true.
    Values of this matrix should be set using the [] or the () operators.
  */
  SerialTriDiMatrix(OrdinalType numRows, OrdinalType numCols, bool zeroOut = true);

  //! Shaped Constructor with Values
  /*!
    \param CV - Enumerated type set to Teuchos::Copy or Teuchos::View.
    \param values - Pointer to an array of ScalarType.
    \param numRowsCols - Number of rows and columns in matrix. 
  */
  SerialTriDiMatrix(DataAccess CV, ScalarType* values, OrdinalType numRowsCols);

  //! Copy Constructor
  /*! \note A deep copy of the \c Source transposed can be obtained if \c trans=Teuchos::TRANS, \c else
    a non-transposed copy of \c Source is made.  There is no storage of the transpose state of the matrix
    within the SerialTriDiMatrix class, so this information will not propogate to any operation performed
    on a matrix that has been copy constructed in transpose.
  */
  SerialTriDiMatrix(const SerialTriDiMatrix<OrdinalType, ScalarType> &Source, ETransp trans = Teuchos::NO_TRANS);

  //! Submatrix Copy Constructor
  /*!
    \param CV - Enumerated type set to Teuchos::Copy or Teuchos::View.
    \param Source - Reference to another TriDi matrix from which values are to be copied.
    \param numRowsCols - The number of rows and columns in this matrix.

    \param startRowCols - The row and col of \c Source from which the submatrix copy should start.

    Creates a shaped matrix with \c numRowsCols rows and columns, which is a submatrix of \c Source.
    If \c startRowCols, then the submatrix is the leading submatrix of \c Source.
    Otherwise, the (1,1) entry in the copied matrix is the (\c startRow, \c startCol) entry of \c Source.
  */
  SerialTriDiMatrix(DataAccess CV, const SerialTriDiMatrix<OrdinalType, ScalarType> &Source, OrdinalType numRowsCols, OrdinalType startRowCols=0);

  //! Destructor
  virtual ~SerialTriDiMatrix();
  //@}

  //! @name Shaping methods.
  //@{
  //! Shape method for changing the size of a SerialTriDiMatrix, initializing entries to zero.
  /*!
    \param numRowsCols - The number of rows in this matrix.

    This method allows the user to define the dimensions of a SerialTriDiMatrix at any point.  This method
    can be called at any point after construction.  Any values previously in this object will be destroyed
    and the resized matrix starts of with all zero values.

    \return Integer error code, set to 0 if successful.
  */
  int shape(OrdinalType numRows);

  //! Same as <tt>shape()</tt> except leaves uninitialized.
  int shapeUninitialized(OrdinalType numRows);

  //! Reshaping method for changing the size of a SerialTriDiMatrix, keeping the entries.
  /*!
    \param numRowsCols - The number of rows in this matrix.

    This method allows the user to redefine the dimensions of a SerialTriDiMatrix at any point.  This method
    can be called at any point after construction.  Any values previously in this object will be copied into
    the reshaped matrix.

    \return Integer error code, set 0 if successful.
  */
  int reshape(OrdinalType numRowsCols);

  //@}

  //! @name Set methods.
  //@{

  //! Copies values from one matrix to another.
  /*!
    The operator= copies the values from one existing SerialTriDiMatrix to another.
    If \c Source is a view (i.e. CV = Teuchos::View), then this method will
    return a view.  Otherwise, it will return a copy of \c Source.  \e this object
    will be resized if it is not large enough to copy \c Source into.
  */
  SerialTriDiMatrix<OrdinalType, ScalarType>& operator= (const SerialTriDiMatrix<OrdinalType, ScalarType>& Source);

  //! Copies values from one matrix to another.
  /*!
    The operator= copies the values from one existing SerialTriDiMatrix to another
    if the dimension of both matrices are the same.  If not, \e this matrix
    will be returned unchanged.
  */
  SerialTriDiMatrix<OrdinalType, ScalarType>& assign (const SerialTriDiMatrix<OrdinalType, ScalarType>& Source);

  //! Set all values in the matrix to a constant value.
  /*!
    \param value - Value to use;
  */
  SerialTriDiMatrix<OrdinalType, ScalarType>& operator= (const ScalarType value) { putScalar(value); return(*this); }

  //! Set all values in the matrix to a constant value.
  /*!
    \param value - Value to use; zero if none specified.
    \return Integer error code, set to 0 if successful.
  */
  int putScalar( const ScalarType value = Teuchos::ScalarTraits<ScalarType>::zero() );

  //! Set all values in the matrix to be random numbers.
  //  int random();

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
    \warning The validity of \c colIndex will only be checked if Teuchos is     configured with
    --enable-teuchos-abc.
  */
  //  ScalarType* operator [] (OrdinalType colIndex);

  //! Column access method (const).
  /*! Returns the pointer to the ScalarType array at the jth column if A[j] is specified, the expression
    A[j][i] will return the same element as A(i,j).

    \return Pointer to the ScalarType array at the \c colIndex column ( \c values_+colIndex*stride_ ).
    \warning The validity of \c colIndex will only be checked if Teuchos is     configured with
    --enable-teuchos-abc.
  */
  //  const ScalarType* operator [] (OrdinalType colIndex) const;

  //! Data array access method.
  /*! \return Pointer to the ScalarType data array contained in the object. */
  ScalarType* values() const { return(values_); }

  ScalarType* D() const { return D_;}
  ScalarType* DL() const { return DL_;}
  ScalarType* DU() const { return DU_;}
  ScalarType* DU2() const { return DU2_;}

  //@}

  //! @name Mathematical methods.
  //@{

  //! Add another matrix to \e this matrix.
  /*! Add \c Source to \e this if the dimension of both matrices are the same.  If not, \e this matrix
    will be returned unchanged.
  */
  SerialTriDiMatrix<OrdinalType, ScalarType>& operator+= (const SerialTriDiMatrix<OrdinalType, ScalarType>& Source);

  //! Subtract another matrix from \e this matrix.
  /*! Subtract \c Source from \e this if the dimension of both matrices are the same.  If not, \e this matrix
    will be returned unchanged.
  */
  SerialTriDiMatrix<OrdinalType, ScalarType>& operator-= (const SerialTriDiMatrix<OrdinalType, ScalarType>& Source);

  //! Scale \c this matrix by \c alpha; \c *this = \c alpha*\c *this.
  /*!
    \param alpha Scalar to multiply \e this by.
  */
  SerialTriDiMatrix<OrdinalType, ScalarType>& operator*= (const ScalarType alpha);

  //! Scale \c this matrix by \c alpha; \c *this = \c alpha*\c *this.
  /*!
    \param alpha Scalar to multiply \e this by.
    \return Integer error code, set to 0 if successful.
  */
  int scale ( const ScalarType alpha );

  //! Point-wise scale \c this matrix by \c A; i.e. *this(i,j) *= A(i,j)
  /*! The values of \c *this matrix will be point-wise scaled by the values in A.
    If A and \c this matrix are not the same dimension \c this will be returned unchanged.

    \param B Teuchos::SerialTriDiMatrix used to perform element-wise scaling of \e this.
    \return Integer error code, set to 0 if successful.
  */
  int scale ( const SerialTriDiMatrix<OrdinalType, ScalarType>& A );

  //! Multiply \c A * \c B and add them to \e this; \e this = \c beta * \e this + \c alpha*A*B.
  /*!
    \param transa - Use the transpose of \c A if transa = Teuchos::TRANS, else don't use the
    transpose if transa = Teuchos::NO_TRANS.
    \param transb - Use the transpose of \c B if transb = Teuchos::TRANS, else don't use the
    transpose if transb = Teuchos::NO_TRANS.
    \param alpha - The scaling factor for \c A * \c B.
    \param A - SerialTriDiMatrix
    \param B - SerialTriDiMatrix
    \param beta - The scaling factor for \e this.

    If the matrices \c A and \c B are not of the right dimension, consistent with \e this, then \e this
    matrix will not be altered and -1 will be returned.
    \return Integer error code, set to 0 if successful.
  */
  //int multiply (ETransp transa, ETransp transb, ScalarType alpha, const SerialTriDiMatrix<OrdinalType, ScalarType> &A, const SerialTriDiMatrix<OrdinalType, ScalarType> &B, ScalarType beta);

  //! Multiply \c A and \c B and add them to \e this; \e this = \c beta * \e this + \c alpha*A*B or \e this = \c beta * \e this + \c alpha*B*A.
  /*!
    \param sideA - Which side is A on for the multiplication to B, A*B (Teuchos::LEFT_SIDE) or B*A (Teuchos::RIGHT_SIDE).
    \param alpha - The scaling factor for \c A * \c B, or \c B * \c A.
    \param A - SerialSymTriDiMatrix (a serial SPD TriDi matrix)
    \param B - SerialTriDiMatrix (a serial TriDi matrix)
    \param beta - The scaling factor for \e this.

    If the matrices \c A and \c B are not of the right dimension, consistent with \e this, then \e this
    matrix will not be altered and -1 will be returned.
    \return Integer error code, set to 0 if successful.
  */
  //int multiply (ESide sideA, ScalarType alpha, const SerialSymTriDiMatrix<OrdinalType, ScalarType> &A, const SerialTriDiMatrix<OrdinalType, ScalarType> &B, ScalarType beta);

  //@}

  //! @name Comparison methods.
  //@{

  //! Equality of two matrices.
  /*! \return True if \e this matrix and \c Operand are of the same shape (rows and columns) and have
    the same entries, else False will be returned.
  */
  bool operator== (const SerialTriDiMatrix<OrdinalType, ScalarType> &Operand) const;

  //! Inequality of two matrices.
  /*! \return True if \e this matrix and \c Operand of not of the same shape (rows and columns) or don't
    have the same entries, else False will be returned.
  */
  bool operator!= (const SerialTriDiMatrix<OrdinalType, ScalarType> &Operand) const;

  //@}

  //! @name Attribute methods.
  //@{

  //! Returns the row dimension of this matrix.
  OrdinalType numRowsCols() const { return(numRowsCols_); }

  //! Returns the column dimension of this matrix.
  //  OrdinalType numCols() const { return(numRowsCols_); }

  //! Returns the stride between the columns of this matrix in memory.
  //  OrdinalType stride() const { return(stride_); }

  //! Returns whether this matrix is empty.
  bool empty() const { return(numRowsCols_ == 0); }
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
  void copyMat(SerialTriDiMatrix<OrdinalType, ScalarType> matrix,
               OrdinalType startCol,
               ScalarType alpha = ScalarTraits<ScalarType>::zero() );
  void deleteArrays();
  void checkIndex( OrdinalType rowIndex, OrdinalType colIndex = 0 ) const;
  OrdinalType numRowsCols_;

  bool valuesCopied_;
  ScalarType* values_;
  ScalarType* DL_;
  ScalarType* D_;
  ScalarType* DU_;
  ScalarType* DU2_;

}; // class Teuchos_SerialTriDiMatrix

//----------------------------------------------------------------------------------------------------
//  Constructors and Destructor
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
SerialTriDiMatrix<OrdinalType, ScalarType>::SerialTriDiMatrix()
  :
    CompObject(),
    numRowsCols_(0),
    valuesCopied_(false),
    values_(0),
    DL_(NULL),
    D_(NULL),
    DU_(NULL),
    DU2_(NULL)
{}

template<typename OrdinalType, typename ScalarType>
SerialTriDiMatrix<OrdinalType, ScalarType>::SerialTriDiMatrix( OrdinalType numRowsCols_in, OrdinalType /* numCols_in */, bool zeroOut)
  : CompObject(), numRowsCols_(numRowsCols_in) {

  OrdinalType numvals = (numRowsCols_ == 1) ? 1 :  4*(numRowsCols_-1);
  values_ = new ScalarType [numvals];
  DL_  = values_;
  D_   = DL_  + (numRowsCols_-1);
  DU_  = D_   +  numRowsCols_;
  DU2_ = DU_  + (numRowsCols_-1);

  valuesCopied_ = true;
  if (zeroOut == true)
    putScalar();
}

template<typename OrdinalType, typename ScalarType>
SerialTriDiMatrix<OrdinalType, ScalarType>::SerialTriDiMatrix(DataAccess CV, ScalarType* values_in, OrdinalType numRowsCols_in  )
  : CompObject(), numRowsCols_(numRowsCols_in),
    valuesCopied_(false), values_(values_in)
{
  const OrdinalType numvals = (numRowsCols_ == 1) ? 1 :  4*(numRowsCols_-1);
  if(CV == Copy) {
    values_ = new ScalarType[numvals];
    valuesCopied_ = true;
  }
  else //CV == View
  {
    values_ = values_in;
    valuesCopied_ = false;
  }
  DL_ = values_;
  D_  = DL_    + (numRowsCols_-1);
  DU_ = D_     + numRowsCols_;
  DU2_ = DU_   + (numRowsCols_-1);
  if(CV == Copy) {
    for(OrdinalType i = 0 ; i < numRowsCols_ ; ++i )
      values_[i] = values_in[i];
  }
}

template<typename OrdinalType, typename ScalarType>
SerialTriDiMatrix<OrdinalType, ScalarType>::SerialTriDiMatrix(const SerialTriDiMatrix<OrdinalType, ScalarType> &Source, ETransp trans) : CompObject(), BLAS<OrdinalType,ScalarType>(), numRowsCols_(0), valuesCopied_(true), values_(0)
{
  if ( trans == Teuchos::NO_TRANS )  {
    numRowsCols_ = Source.numRowsCols_;

    const OrdinalType numvals = (numRowsCols_ == 1) ? 1 : 4*(numRowsCols_-1);
    values_ = new ScalarType[numvals];
    DL_ = values_;
    D_  = DL_+ (numRowsCols_-1);
    DU_ = D_ + numRowsCols_;
    DU2_ = DU_ + (numRowsCols_-1);

    copyMat(Source, 0, 0);
  }
  else if ( trans == Teuchos::CONJ_TRANS && ScalarTraits<ScalarType>::isComplex )
    {
      numRowsCols_ = Source.numRowsCols_;
      const OrdinalType numvals = (numRowsCols_ == 1) ? 1 :  4*(numRowsCols_-1);
      values_ = new ScalarType[numvals];
      DL_ = values_;
      D_  = DL_+(numRowsCols_-1);
      DU_ = D_ + numRowsCols_;
      DU2_ = DU_ + (numRowsCols_-1);

      OrdinalType min = numRowsCols_;
      if(min > Source.numRowsCols_) min = Source.numRowsCols_;

      for(OrdinalType i = 0 ; i< min ; ++i) {
        D_[i] = Teuchos::ScalarTraits<ScalarType>::conjugate(Source.D_[i]);
        if(i < (min-1)) {
          DL_[i] = Teuchos::ScalarTraits<ScalarType>::conjugate(Source.DL_[i]);
          DU_[i] = Teuchos::ScalarTraits<ScalarType>::conjugate(Source.DU_[i]);
        }
        if(i < (min-2)) {
          DU2_[i] = Teuchos::ScalarTraits<ScalarType>::conjugate(Source.DU2_[i]);
        }
      }
    }
  else
    {
      numRowsCols_ = Source.numRowsCols_;
      const OrdinalType numvals = (numRowsCols_  == 1) ? 1 : 4*(numRowsCols_-1);
      values_ = new ScalarType[numvals];
      OrdinalType min = numRowsCols_;
      if(min > Source.numRowsCols_) min = Source.numRowsCols_;
      for(OrdinalType i = 0 ; i< min ; ++i) {
        D_[i] = Source.D_[i];
        if(i < (min-1)) {
          DL_[i] = Source.DL_[i];
          DU_[i] = Source.DU_[i];
        }
        if(i < (min-2)) {
          DU2_[i] = Source.DU2_[i];
        }
    }
  }
}

template<typename OrdinalType, typename ScalarType>
SerialTriDiMatrix<OrdinalType, ScalarType>::SerialTriDiMatrix(
  DataAccess CV, const SerialTriDiMatrix<OrdinalType, ScalarType> &Source,
  OrdinalType numRowsCols_in, OrdinalType startRow )
  : CompObject(), numRowsCols_(numRowsCols_in),
    valuesCopied_(false), values_(Source.values_) {
  if(CV == Copy)
    {
      const OrdinalType numvals = (numRowsCols_ == 1) ? 1 : 4*(numRowsCols_-1);
      values_ = new ScalarType[numvals];
      copyMat(Source, startRow);
      valuesCopied_ = true;
    }
  else // CV == View
    {
      // \todo ???
      //    values_ = values_ + (stride_ * startCol) + startRow;
    }
}

template<typename OrdinalType, typename ScalarType>
SerialTriDiMatrix<OrdinalType, ScalarType>::~SerialTriDiMatrix()
{
  deleteArrays();
}

//----------------------------------------------------------------------------------------------------
//  Shape methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
int SerialTriDiMatrix<OrdinalType, ScalarType>::shape( OrdinalType numRowsCols_in  )
{
  deleteArrays(); // Get rid of anything that might be already allocated
  numRowsCols_ = numRowsCols_in;
  const OrdinalType numvals = ( numRowsCols_ == 1) ? 1 :  4*(numRowsCols_-1);
  values_ = new ScalarType[numvals];

  putScalar();
  valuesCopied_ = true;
  return(0);
}

template<typename OrdinalType, typename ScalarType>
int SerialTriDiMatrix<OrdinalType, ScalarType>::shapeUninitialized(  OrdinalType numRowsCols_in  )
{
  deleteArrays(); // Get rid of anything that might be already allocated
  numRowsCols_ = numRowsCols_in;
  const OrdinalType numvals = ( numRowsCols_ == 1) ? 1 :  4*(numRowsCols_-1);
  values_ = new ScalarType[numvals];
  valuesCopied_ = true;
  return(0);
}

template<typename OrdinalType, typename ScalarType>
int SerialTriDiMatrix<OrdinalType, ScalarType>::reshape(  OrdinalType numRowsCols_in  )
{
  if(numRowsCols_in <1 ) {
    deleteArrays();
    return 0;
  }
  // Allocate space for new matrix
  const OrdinalType numvals = ( numRowsCols_in == 1) ? 1 :  4*(numRowsCols_in - 1);
  ScalarType* values_tmp = new ScalarType[numvals];
  ScalarType zero = ScalarTraits<ScalarType>::zero();
  for(OrdinalType i= 0; i<numvals ; ++i)
    values_tmp[i] = zero;

  OrdinalType min = TEUCHOS_MIN(numRowsCols_, numRowsCols_in);
  ScalarType* dl = values_tmp;
  ScalarType* d = values_tmp + (numRowsCols_in-1);
  ScalarType* du = d+(numRowsCols_in);
  ScalarType* du2 = du+(numRowsCols_in - 1);

  if(values_ != 0) {
   for(OrdinalType i = 0 ; i< min ; ++i) {
      dl[i] = DL_[i];
      d[i]  = D_[i];
      du[i] = DU_[i];
      du2[i] = DU2_[i];
    }
  }

  deleteArrays(); // Get rid of anything that might be already allocated
  numRowsCols_ = numRowsCols_in;

  values_ = values_tmp; // Set pointer to new A
  DL_ = dl;
  D_  = d;
  DU_ = du;
  DU2_ = du2;

  valuesCopied_ = true;
  return(0);
}

//----------------------------------------------------------------------------------------------------
//   Set methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
int SerialTriDiMatrix<OrdinalType, ScalarType>::putScalar( const ScalarType value_in ) {
  // Set each value of the TriDi matrix to "value".
  const OrdinalType numvals = (numRowsCols_ == 1) ? 1 : 4*(numRowsCols_-1);

  for(OrdinalType i = 0; i<numvals ; ++i)
    {
      values_[i] = value_in;
    }
  return 0;
}

// template<typename OrdinalType, typename ScalarType>
// int SerialTriDiMatrix<OrdinalType, ScalarType>::random()
// {
//   // Set each value of the TriDi matrix to a random value.
//   for(OrdinalType j = 0; j < numCols_; j++)
//   {
//     for(OrdinalType i = 0; i < numRowsCols_; i++)
//        {
//          values_[i + j*stride_] = ScalarTraits<ScalarType>::random();
//        }
//   }
//   return 0;
// }

template<typename OrdinalType, typename ScalarType>
SerialTriDiMatrix<OrdinalType,ScalarType>&
SerialTriDiMatrix<OrdinalType, ScalarType>::operator=(const SerialTriDiMatrix<OrdinalType,ScalarType>& Source )
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
    numRowsCols_ = Source.numRowsCols_;
    values_ = Source.values_;
  }
  else {
    // If we were a view, we will now be a copy.
    if(!valuesCopied_) {
      numRowsCols_ = Source.numRowsCols_;
      const OrdinalType numvals = ( Source.numRowsCols_ == 1) ? 1 :  4*(Source.numRowsCols_ - 1);
      if(numvals > 0) {
        values_ = new ScalarType[numvals];
        valuesCopied_ = true;
      }
      else {
        values_ = 0;
      }
    }
    // If we were a copy, we will stay a copy.
    else {
      // we need to allocate more space (or less space)
      deleteArrays();
      numRowsCols_ = Source.numRowsCols_;
      const OrdinalType numvals = ( Source.numRowsCols_ == 1) ? 1 :  4*(Source.numRowsCols_ - 1);
      if(numvals > 0) {
        values_ = new ScalarType[numvals];
        valuesCopied_ = true;
      }
    }
    DL_ = values_;
    if(values_ != 0) {
      D_  = DL_    + (numRowsCols_-1);
      DU_ = D_     +  numRowsCols_;
      DU2_ = DU_   + (numRowsCols_-1);

      OrdinalType min = TEUCHOS_MIN(numRowsCols_, Source.numRowsCols_);
      for(OrdinalType i = 0 ; i < min ; ++i ) {
        D_[i] = Source.D()[i];
        if(i< (min-1 ) )
          {
            DL_[i] = Source.DL()[i];
            DU_[i] = Source.DU()[i];
          }
        if(i< (min-2))
          DU2_[i] = Source.DU2()[i];
      }

    }
    else {
      D_ = DU_ = DU2_ = 0;
    }
  }
  return(*this);
}

template<typename OrdinalType, typename ScalarType>
SerialTriDiMatrix<OrdinalType, ScalarType>& SerialTriDiMatrix<OrdinalType, ScalarType>::operator+= (const SerialTriDiMatrix<OrdinalType,ScalarType>& Source )
{
  // Check for compatible dimensions
  if ((numRowsCols_ != Source.numRowsCols_) )
    {
      TEUCHOS_CHK_REF(*this); // Return *this without altering it.
    }
  copyMat(Source, 0, ScalarTraits<ScalarType>::one());
  return(*this);
}

template<typename OrdinalType, typename ScalarType>
SerialTriDiMatrix<OrdinalType, ScalarType>& SerialTriDiMatrix<OrdinalType, ScalarType>::operator-= (const SerialTriDiMatrix<OrdinalType,ScalarType>& Source )
{
  // Check for compatible dimensions
  if ((numRowsCols_ != Source.numRowsCols_) )
    {
      TEUCHOS_CHK_REF(*this); // Return *this without altering it.
    }
  copyMat(Source, 0, -ScalarTraits<ScalarType>::one());
  return(*this);
}

template<typename OrdinalType,typename ScalarType>
SerialTriDiMatrix<OrdinalType,ScalarType> & SerialTriDiMatrix<OrdinalType,ScalarType>::assign(const SerialTriDiMatrix<OrdinalType,ScalarType> & Source)
{
  if(this == &Source)
    return(*this); // Special case of source same as target
  if((!valuesCopied_) && (!Source.valuesCopied_) && (values_ == Source.values_))
    return(*this); // Special case of both are views to same data.

  // Check for compatible dimensions
  if ((numRowsCols_ != Source.numRowsCols_) )
    {
      TEUCHOS_CHK_REF(*this); // Return *this without altering it.
    }
  copyMat(Source,  0, 0);
  return(*this);
}

//----------------------------------------------------------------------------------------------------
//   Accessor methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType,typename ScalarType>
inline const ScalarType& SerialTriDiMatrix<OrdinalType,ScalarType>::operator () (OrdinalType rowIndex, OrdinalType colIndex) const
{
  OrdinalType diff = colIndex - rowIndex;

  //#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  checkIndex( rowIndex, colIndex );
  //#endif
  switch (diff) {
  case -1:
    return DL_[colIndex];
  case 0:
    return D_[colIndex];
  case 1:
    return DU_[rowIndex];
  case 2:
    return DU2_[rowIndex];
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::out_of_range,
                               "SerialTriDiMatrix<T>::operator (row,col) "
                               "Index (" << rowIndex <<","<<colIndex<<") out of range ");
  }
}

template<typename OrdinalType,typename ScalarType>
inline ScalarType& Teuchos::SerialTriDiMatrix<OrdinalType,ScalarType>::operator () (OrdinalType rowIndex, OrdinalType colIndex)
{
  OrdinalType diff = colIndex - rowIndex;
  //#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  checkIndex( rowIndex, colIndex );
  //#endif
  switch (diff) {
  case -1:
    return DL_[colIndex];
  case 0:
    return D_[colIndex];
  case 1:
    return DU_[rowIndex];
  case 2:
    return DU2_[rowIndex];
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::out_of_range,
                               "SerialTriDiMatrix<T>::operator (row,col) "
                               "Index (" << rowIndex <<","<<colIndex<<") out of range ");
  }
}

//----------------------------------------------------------------------------------------------------
//   Norm methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType,typename ScalarType>
typename ScalarTraits<ScalarType>::magnitudeType SerialTriDiMatrix<OrdinalType,ScalarType>::normOne() const
{
  OrdinalType i, j;
  typename ScalarTraits<ScalarType>::magnitudeType anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
  typename ScalarTraits<ScalarType>::magnitudeType absSum = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());

  // Fix this for Tri DI

  for(j = 0; j < numRowsCols_; j++)
    {
      ScalarType sum = 0;
      if(j-1>=0) sum +=  ScalarTraits<ScalarType>::magnitude((*this)(j-1,j));
      sum+= ScalarTraits<ScalarType>::magnitude((*this)(j,j));
      if(j+1<numRowsCols_) sum+= ScalarTraits<ScalarType>::magnitude((*this)(j+1,j));
      absSum = ScalarTraits<ScalarType>::magnitude(sum);
      if(absSum > anorm)
        {
          anorm = absSum;
        }
    }
  updateFlops(numRowsCols_ * numRowsCols_);
  return(anorm);
}

template<typename OrdinalType, typename ScalarType>
typename ScalarTraits<ScalarType>::magnitudeType SerialTriDiMatrix<OrdinalType, ScalarType>::normInf() const
{
  OrdinalType i,j;
  typename ScalarTraits<ScalarType>::magnitudeType sum, anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());

  for (i = 0; i < numRowsCols_; i++) {
    sum = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
    for (j=i-1; j<= i+1; j++) {
      if(j >= 0 && j < numRowsCols_) sum += ScalarTraits<ScalarType>::magnitude((*this)(i,j));
    }
    anorm = TEUCHOS_MAX( anorm, sum );
  }
  updateFlops(numRowsCols_ * numRowsCols_);
  return(anorm);
}

template<typename OrdinalType, typename ScalarType>
typename ScalarTraits<ScalarType>::magnitudeType SerialTriDiMatrix<OrdinalType, ScalarType>::normFrobenius() const
{
  // \todo fix this
  OrdinalType i, j;
  typename ScalarTraits<ScalarType>::magnitudeType anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
  for (j = 0; j < numRowsCols_; j++) {
    for (i = j-1; i <= j+1; i++) {
      if(i >= 0 && i < numRowsCols_ ) anorm += ScalarTraits<ScalarType>::magnitude((*this)(i,j));
    }
  }
  anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::squareroot(anorm));
  updateFlops( (numRowsCols_ == 1) ? 1 :  4*(numRowsCols_-1) );
  return(anorm);
}

//----------------------------------------------------------------------------------------------------
//   Comparison methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
bool SerialTriDiMatrix<OrdinalType, ScalarType>::operator== (const SerialTriDiMatrix<OrdinalType, ScalarType> &Operand) const
{
  bool allmatch = true;
  // bool result = 1; // unused
  if((numRowsCols_ != Operand.numRowsCols_) )
    {
      // result = 0; // unused
    }
  else
    {
      OrdinalType numvals = (numRowsCols_ == 1)? 1 : 4*(numRowsCols_ -1 );

      for(OrdinalType i = 0; i< numvals; ++i)
        allmatch &= (Operand.values_[i] == values_[i]);
    }
  return allmatch;
}

template<typename OrdinalType, typename ScalarType>
bool SerialTriDiMatrix<OrdinalType, ScalarType>::operator!= (const SerialTriDiMatrix<OrdinalType, ScalarType> &Operand) const {
  return !((*this) == Operand);
}

//----------------------------------------------------------------------------------------------------
//   Multiplication method
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
SerialTriDiMatrix<OrdinalType, ScalarType>& SerialTriDiMatrix<OrdinalType, ScalarType>::operator*= (const ScalarType alpha )
{
  this->scale( alpha );
  return(*this);
}

template<typename OrdinalType, typename ScalarType>
int SerialTriDiMatrix<OrdinalType, ScalarType>::scale( const ScalarType alpha )
{
  OrdinalType i;
  OrdinalType numvals = (numRowsCols_ == 1)? 1 : 4*(numRowsCols_ -1 );
  for (i=0; i < numvals ; ++i ) {
    values_[i] *= alpha;
  }
  return(0);
}

template<typename OrdinalType, typename ScalarType>
int SerialTriDiMatrix<OrdinalType, ScalarType>::scale( const SerialTriDiMatrix<OrdinalType,ScalarType>& A )
{
  OrdinalType j;

  // Check for compatible dimensions
  if ((numRowsCols_ != A.numRowsCols_) )
    {
      TEUCHOS_CHK_ERR(-1); // Return error
    }
  OrdinalType numvals = (numRowsCols_ == 1)? 1 : 4*(numRowsCols_ -1 );
  for (j=0; j<numvals; j++) {
    values_[j] = A.values_ * values_[j];
  }
  updateFlops( numvals );
  return(0);
}

template<typename OrdinalType, typename ScalarType>
void SerialTriDiMatrix<OrdinalType, ScalarType>::print(std::ostream& os) const
{
  os << std::endl;
  if(valuesCopied_)
    os << "A_Copied: yes" << std::endl;
  else
    os << "A_Copied: no" << std::endl;
  os << "Rows and Columns: " << numRowsCols_ << std::endl;
  if(numRowsCols_ == 0) {
    os << "(matrix is empty, no values to display)" << std::endl;
    return;
  }
  else
    {
      os << "DL: "<<std::endl;
      for(int i=0;i<numRowsCols_-1;++i)
        os << DL_[i]<<" ";
      os << std::endl;
      os << "D: "<<std::endl;
      for(int i=0;i<numRowsCols_;++i)
        os << D_[i]<<" ";
      os << std::endl;
      os << "DU: "<<std::endl;
      for(int i=0;i<numRowsCols_-1;++i)
        os << DU_[i]<<" ";
      os << std::endl;
      os << "DU2: "<<std::endl;
      for(int i=0;i<numRowsCols_-2;++i)
        os << DU2_[i]<<" ";
      os << std::endl;
    }

  os <<" square format:"<<std::endl;
  for(int i=0 ; i < numRowsCols_ ; ++i )  {
    for(int j=0;j<numRowsCols_;++j)  {
      if ( j >= i-1  && j <= i+1) {
        os << (*this)(i,j)<<" ";
      }
      else
        os <<". ";
    }
    os << std::endl;
  }

}

//----------------------------------------------------------------------------------------------------
//   Protected methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
inline void SerialTriDiMatrix<OrdinalType, ScalarType>::checkIndex( OrdinalType rowIndex, OrdinalType colIndex ) const
{
  OrdinalType diff = colIndex - rowIndex;
  TEUCHOS_TEST_FOR_EXCEPTION(rowIndex < 0 || rowIndex >= numRowsCols_, std::out_of_range,
                             "SerialTriDiMatrix<T>::checkIndex: "
                             "Row index " << rowIndex << " out of range [0, "<< numRowsCols_ << "]");
  TEUCHOS_TEST_FOR_EXCEPTION(colIndex < 0 || colIndex >= numRowsCols_,
                             std::out_of_range,
                             "SerialTriDiMatrix<T>::checkIndex: "
                             "Col index " << colIndex << " out of range [0, "<< numRowsCols_ << "]");
  TEUCHOS_TEST_FOR_EXCEPTION(diff > 2 || diff < -1 , std::out_of_range,
                             "SerialTriDiMatrix<T>::checkIndex: "
                             "index difference " << diff << " out of range [-1, 2]");
}

template<typename OrdinalType, typename ScalarType>
void SerialTriDiMatrix<OrdinalType, ScalarType>::deleteArrays(void)
{
  if (valuesCopied_)
    {
      delete [] values_;
      values_ = 0;
      valuesCopied_ = false;
    }
}

template<typename OrdinalType, typename ScalarType>
void SerialTriDiMatrix<OrdinalType, ScalarType>::copyMat(SerialTriDiMatrix<OrdinalType, ScalarType> inputMatrix,
                                                         OrdinalType startRowCol,
                                                         ScalarType alpha )
{
  OrdinalType i;
  OrdinalType max = inputMatrix.numRowsCols_;
  if(max > numRowsCols_ ) max = numRowsCols_;
  if(startRowCol > max ) return; //

  for(i = startRowCol ; i < max ; ++i) {

    if (alpha != Teuchos::ScalarTraits<ScalarType>::zero() ) {
      // diagonal
      D()[i] += inputMatrix.D()[i];
      if(i<(max-1) && (i-1) >= startRowCol) {
        DL()[i] += inputMatrix.DL()[i];
        DU()[i] += inputMatrix.DU()[i];
      }
      if(i<(max-2) && (i-2) >= startRowCol) {
        DU2()[i] += inputMatrix.DU2()[i];
      }
    }
    else {
      // diagonal
      D()[i] = inputMatrix.D()[i];
      if(i<(max-1) && (i-1) >= startRowCol) {
        DL()[i] = inputMatrix.DL()[i];
        DU()[i] = inputMatrix.DU()[i];
      }
      if(i<(max-2) && (i-2) >= startRowCol) {
        DU2()[i] = inputMatrix.DU2()[i];
      }
    }
  }
}

}


#endif /* _TEUCHOS_SERIALTRIDIMATRIX_HPP_ */
