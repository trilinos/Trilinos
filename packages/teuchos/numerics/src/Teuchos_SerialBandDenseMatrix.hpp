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

#ifndef _TEUCHOS_SERIALBANDDENSEMATRIX_HPP_
#define _TEUCHOS_SERIALBANDDENSEMATRIX_HPP_
/*! \file Teuchos_SerialBandDenseMatrix.hpp
  \brief Templated serial dense matrix class
*/

#include "Teuchos_CompObject.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

/*! \class Teuchos::SerialBandDenseMatrix
    \brief This class creates and provides basic support for banded dense matrices of templated type.

    The Teuchos_SerialBandDenseMatrix class enables the construction and use of banded dense matrices of templated type.

    The Teuchos::SerialBandDenseMatrix class is intended to provide full-featured support for solving
    linear and eigen system problems for banded matrices.  It is written on
    top of BLAS and LAPACK and thus has excellent performance and numerical capabilities.  Using this
    class, one can either perform simple factorizations and solves or apply all the tricks available in
    LAPACK to get the best possible solution for very ill-conditioned problems.

    <b>Teuchos::SerialBandDenseMatrix vs. Teuchos::LAPACK</b>

    The Teuchos::LAPACK class provides access to most of the same functionality as Teuchos::SerialBandDenseMatrix.
    The primary difference is that Teuchos::LAPACK is a "thin" layer on top of LAPACK and
    Teuchos::SerialBandDenseMatrix attempts to provide easy access to the more sophisticated aspects of
    solving dense linear and eigensystems.
<ul>
<li> When you should use Teuchos::LAPACK:  If you are simply looking for a convenient wrapper around
the Fortran LAPACK routines and you have a well-conditioned problem, you should probably use Teuchos::LAPACK directly.
<li> When you should use Teuchos::SerialBandDenseMatrix: If you want to (or potentially want to) solve ill-conditioned
     problems or want to work with a more object-oriented interface, you should probably use Teuchos::SerialBandDenseMatrix.
</ul>

<b>Constructing Teuchos::SerialBandDenseMatrix Objects</b>

There are three Teuchos::SerialBandDenseMatrix constructors.  The first constructs a zero-sized object
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

<b>Format of the matrix A</b>

The matrix A is stored as a banded matrix AB according to the LAPACK format. Consider using
the non-member conversion routines generalToBanded and bandedToGeneral if the full Teuchos::SerialDenseMatrix is
already in storage. However, it is more efficient to simply construct the Teuchos::SerialBandDenseMatrix with the desired
parameters and use the provided matrix access operators to assign (legal) entries so that the full
rectangular matrix does not need to be stored. A full description of the LAPACK banded format can be found at
http://www.netlib.org/lapack/lug/node124.html. Briefly, the banded storage format of AB used internally is as follows:
<ul>
<li> The dimension of AB is KL+KU+1 by N, where KL and KU are the lower and upper bandwidths, respectively.
<li> The entries are stored in a compressed format: AB(KU+i-j,j) = A(i,j) for max(0,j-KU)<=i<=min(M-1,j+KL).
</ul>

<b>Extracting Data from Teuchos::SerialBandDenseMatrix Objects</b>

Once a Teuchos::SerialBandDenseMatrix is constructed, it is possible to view the data via access functions.

\warning Use of these access functions cam be \e extremely dangerous from a data hiding perspective.

<b>Vector and Utility Functions</b>

Once a Teuchos::SerialBandDenseMatrix is constructed, several mathematical functions can be applied to
the object.  Specifically:
<ul>
  <li> Multiplication.
  <li> Norms.
</ul>

*/

/** \example DenseMatrix/cxx_main_band.cpp
    This is an example of how to use the Teuchos::SerialBandDenseMatrix class.
*/

namespace Teuchos {

template<typename OrdinalType, typename ScalarType>
class SerialBandDenseMatrix : public CompObject, public Object, public BLAS<OrdinalType, ScalarType>
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
  SerialBandDenseMatrix();

  //! Shaped Constructor
  /*!
    \param numRows - Number of rows in matrix.
    \param numCols - Number of columns in matrix.
    \param kl - Lower bandwidth of matrix.
    \param ku - Upper bandwidth of matrix.
    \param zeroOut - Initializes values to 0 if true (default)

    Creates a shaped matrix with \c numRows rows and \c numCols cols.  All values are initialized to 0 when \c zeroOut is true.
    Values of this matrix should be set using the [] or the () operators.
  */
  SerialBandDenseMatrix(OrdinalType numRows, OrdinalType numCols, OrdinalType kl, OrdinalType ku, bool zeroOut = true);

  //! Shaped Constructor with Values
  /*!
    \param CV - Enumerated type set to Teuchos::Copy or Teuchos::View.
    \param values - Pointer to an array of ScalarType.  The first column starts at \c values,
     the second at \c values+stride, etc.
    \param stride - The stride between the columns of the matrix in memory.
    \param numRows - Number of rows in matrix.
    \param numCols - Number of columns in matrix.
    \param kl - Lower bandwidth of matrix.
    \param ku - Upper bandwidth of matrix.
  */
  SerialBandDenseMatrix(DataAccess CV, ScalarType* values, OrdinalType stride, OrdinalType numRows, OrdinalType numCols, OrdinalType kl, OrdinalType ku);

  //! Copy Constructor
  /*! \note A deep copy of the \c Source transposed can be obtained if \c trans=Teuchos::TRANS, \c else
    a non-transposed copy of \c Source is made.  There is no storage of the transpose state of the matrix
    within the SerialBandDenseMatrix class, so this information will not propogate to any operation performed
    on a matrix that has been copy constructed in transpose.
  */
  SerialBandDenseMatrix(const SerialBandDenseMatrix<OrdinalType, ScalarType> &Source, ETransp trans = Teuchos::NO_TRANS);

  //! Submatrix Copy Constructor
  /*!
    \param CV - Enumerated type set to Teuchos::Copy or Teuchos::View.
    \param Source - Reference to another dense matrix from which values are to be copied.
    \param numRows - The number of rows in this matrix.
    \param numCols - The number of columns in this matrix.
    \param kl - Lower bandwidth of matrix.
    \param ku - Upper bandwidth of matrix.
    \param startRow - The row of \c Source from which the submatrix copy should start.
    \param startCol - The column of \c Source from which the submatrix copy should start.

    Creates a shaped matrix with \c numRows rows and \c numCols columns, which is a submatrix of \c Source.
    If \c startRow and \c startCol are not given, then the submatrix is the leading submatrix of \c Source.
    Otherwise, the (1,1) entry in the copied matrix is the (\c startRow, \c startCol) entry of \c Source.
  */
  SerialBandDenseMatrix(DataAccess CV, const SerialBandDenseMatrix<OrdinalType, ScalarType> &Source, OrdinalType numRows, OrdinalType numCols, OrdinalType startCol=0);

  //! Destructor
  virtual ~SerialBandDenseMatrix();
  //@}

  //! @name Shaping methods.
  //@{
  //! Shape method for changing the size of a SerialBandDenseMatrix, initializing entries to zero.
  /*!
    \param numRows - The number of rows in this matrix.
    \param numCols - The number of columns in this matrix.
    \param kl - Lower bandwidth of matrix.
    \param ku - Upper bandwidth of matrix.

    This method allows the user to define the dimensions of a SerialBandDenseMatrix at any point.  This method
    can be called at any point after construction.  Any values previously in this object will be destroyed
    and the resized matrix starts of with all zero values.

    \return Integer error code, set to 0 if successful.
  */
  int shape(OrdinalType numRows, OrdinalType numCols, OrdinalType kl, OrdinalType ku);

  //! Same as <tt>shape()</tt> except leaves uninitialized.
  int shapeUninitialized(OrdinalType numRows, OrdinalType numCols, OrdinalType kl, OrdinalType ku);

  //! Reshaping method for changing the size of a SerialBandDenseMatrix, keeping the entries.
  /*!
    \param numRows - The number of rows in this matrix.
    \param numCols - The number of columns in this matrix.
    \param kl - Lower bandwidth of matrix.
    \param ku - Upper bandwidth of matrix.

    This method allows the user to redefine the dimensions of a SerialBandDenseMatrix at any point.  This method
    can be called at any point after construction.  Any values previously in this object will be copied into
    the reshaped matrix.

    \return Integer error code, set 0 if successful.
  */
  int reshape(OrdinalType numRows, OrdinalType numCols, OrdinalType kl, OrdinalType ku);


  //@}

  //! @name Set methods.
  //@{

  //! Copies values from one matrix to another.
  /*!
    The operator= copies the values from one existing SerialBandDenseMatrix to another.
    If \c Source is a view (i.e. CV = Teuchos::View), then this method will
    return a view.  Otherwise, it will return a copy of \c Source.  \e this object
    will be resized if it is not large enough to copy \c Source into.
  */
  SerialBandDenseMatrix<OrdinalType, ScalarType>& operator= (const SerialBandDenseMatrix<OrdinalType, ScalarType>& Source);

  //! Copies values from one matrix to another.
  /*!
    The operator= copies the values from one existing SerialBandDenseMatrix to another
    if the dimension of both matrices are the same.  If not, \e this matrix
    will be returned unchanged.
  */
  SerialBandDenseMatrix<OrdinalType, ScalarType>& assign (const SerialBandDenseMatrix<OrdinalType, ScalarType>& Source);

  //! Set all values in the matrix to a constant value.
  /*!
    \param value - Value to use;
  */
  SerialBandDenseMatrix<OrdinalType, ScalarType>& operator= (const ScalarType value) { putScalar(value); return(*this); }

  //! Set all values in the matrix to a constant value.
  /*!
    \param value - Value to use; zero if none specified.
    \return Integer error code, set to 0 if successful.
  */
  int putScalar( const ScalarType value = Teuchos::ScalarTraits<ScalarType>::zero() );

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
    \warning The validity of \c colIndex will only be checked if Teuchos is	configured with
    --enable-teuchos-abc.
  */
  ScalarType* operator [] (OrdinalType colIndex);

  //! Column access method (const).
  /*! Returns the pointer to the ScalarType array at the jth column if A[j] is specified, the expression
    A[j][i] will return the same element as A(i,j).

    \return Pointer to the ScalarType array at the \c colIndex column ( \c values_+colIndex*stride_ ).
    \warning The validity of \c colIndex will only be checked if Teuchos is	configured with
    --enable-teuchos-abc.
  */
  const ScalarType* operator [] (OrdinalType colIndex) const;

  //! Data array access method.
  /*! \return Pointer to the ScalarType data array contained in the object. */
  ScalarType* values() const { return(values_); }

  //@}

  //! @name Mathematical methods.
  //@{

  //! Add another matrix to \e this matrix.
  /*! Add \c Source to \e this if the dimension of both matrices are the same.  If not, \e this matrix
    will be returned unchanged.
  */
  SerialBandDenseMatrix<OrdinalType, ScalarType>& operator+= (const SerialBandDenseMatrix<OrdinalType, ScalarType>& Source);

  //! Subtract another matrix from \e this matrix.
  /*! Subtract \c Source from \e this if the dimension of both matrices are the same.  If not, \e this matrix
    will be returned unchanged.
  */
  SerialBandDenseMatrix<OrdinalType, ScalarType>& operator-= (const SerialBandDenseMatrix<OrdinalType, ScalarType>& Source);

  //! Scale \c this matrix by \c alpha; \c *this = \c alpha*\c *this.
  /*!
    \param alpha Scalar to multiply \e this by.
  */
  SerialBandDenseMatrix<OrdinalType, ScalarType>& operator*= (const ScalarType alpha);

  //! Scale \c this matrix by \c alpha; \c *this = \c alpha*\c *this.
  /*!
    \param alpha Scalar to multiply \e this by.
    \return Integer error code, set to 0 if successful.
  */
  int scale ( const ScalarType alpha );

  //! Point-wise scale \c this matrix by \c A; i.e. *this(i,j) *= A(i,j)
  /*! The values of \c *this matrix will be point-wise scaled by the values in A.
    If A and \c this matrix are not the same dimension \c this will be returned unchanged.

    \param B Teuchos::SerialBandDenseMatrix used to perform element-wise scaling of \e this.
    \return Integer error code, set to 0 if successful.
  */
  int scale ( const SerialBandDenseMatrix<OrdinalType, ScalarType>& A );

  //@}

  //! @name Comparison methods.
  //@{

  //! Equality of two matrices.
  /*! \return True if \e this matrix and \c Operand are of the same shape (rows and columns) and have
    the same entries, else False will be returned.
  */
  bool operator== (const SerialBandDenseMatrix<OrdinalType, ScalarType> &Operand) const;

  //! Inequality of two matrices.
  /*! \return True if \e this matrix and \c Operand of not of the same shape (rows and columns) or don't
    have the same entries, else False will be returned.
  */
  bool operator!= (const SerialBandDenseMatrix<OrdinalType, ScalarType> &Operand) const;

  //@}

  //! @name Attribute methods.
  //@{

  //! Returns the row dimension of this matrix.
  OrdinalType numRows() const { return(numRows_); }

  //! Returns the column dimension of this matrix.
  OrdinalType numCols() const { return(numCols_); }

  //! Returns the lower bandwidth of this matrix.
  OrdinalType lowerBandwidth() const { return(kl_); }

  //! Returns the upper bandwidth of this matrix.
  OrdinalType upperBandwidth() const { return(ku_); }

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
  //! Print method.  Defines the behavior of the std::ostream << operator inherited from the Object class.
  virtual void print(std::ostream& os) const;

  //@}
protected:
  void copyMat(ScalarType* inputMatrix, OrdinalType strideInput,
    OrdinalType numRows, OrdinalType numCols, ScalarType* outputMatrix,
    OrdinalType strideOutput, OrdinalType startCol, ScalarType alpha = ScalarTraits<ScalarType>::zero() );
  void deleteArrays();
  void checkIndex( OrdinalType rowIndex, OrdinalType colIndex = 0 ) const;
  OrdinalType numRows_;
  OrdinalType numCols_;
  OrdinalType stride_;
  OrdinalType kl_;
  OrdinalType ku_;
  bool valuesCopied_;
  ScalarType* values_;

}; // class Teuchos_SerialBandDenseMatrix

//----------------------------------------------------------------------------------------------------
//  Constructors and Destructor
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
SerialBandDenseMatrix<OrdinalType, ScalarType>::SerialBandDenseMatrix()
  : CompObject(), Object("Teuchos::SerialBandDenseMatrix"), numRows_(0), numCols_(0), kl_(0), ku_(0), stride_(0), valuesCopied_(false), values_(0)
{}

template<typename OrdinalType, typename ScalarType>
SerialBandDenseMatrix<OrdinalType, ScalarType>::SerialBandDenseMatrix(
  OrdinalType numRows_in, OrdinalType numCols_in, OrdinalType kl_in, OrdinalType ku_in, bool zeroOut
  )
  : CompObject(), Object("Teuchos::SerialBandDenseMatrix"), numRows_(numRows_in), numCols_(numCols_in), kl_(kl_in), ku_(ku_in), stride_(kl_in+ku_in+1)
{
  values_ = new ScalarType[stride_*numCols_];
  if (zeroOut == true)
    putScalar();
  valuesCopied_ = true;
}

template<typename OrdinalType, typename ScalarType>
SerialBandDenseMatrix<OrdinalType, ScalarType>::SerialBandDenseMatrix(
  DataAccess CV, ScalarType* values_in, OrdinalType stride_in, OrdinalType numRows_in,
  OrdinalType numCols_in, OrdinalType kl_in, OrdinalType ku_in
  )
  : CompObject(), Object("Teuchos::SerialBandDenseMatrix"), numRows_(numRows_in), numCols_(numCols_in), kl_(kl_in), ku_(ku_in), stride_(stride_in), valuesCopied_(false), values_(values_in)
{
  if(CV == Copy) {
    stride_ = kl_+ku_+1;
    values_ = new ScalarType[stride_*numCols_];
    copyMat(values_in, stride_in, numRows_, numCols_, values_, stride_, 0 );
    valuesCopied_ = true;
  }
}

template<typename OrdinalType, typename ScalarType>
SerialBandDenseMatrix<OrdinalType, ScalarType>::SerialBandDenseMatrix(const SerialBandDenseMatrix<OrdinalType, ScalarType> &Source, ETransp trans) : CompObject(), Object("Teuchos::SerialBandDenseMatrix"), numRows_(0), numCols_(0), kl_(0), ku_(0), stride_(0), valuesCopied_(true), values_(0)
{

  if ( trans == Teuchos::NO_TRANS ) {

    numRows_ = Source.numRows_;
    numCols_ = Source.numCols_;
    kl_ = Source.kl_;
    ku_ = Source.ku_;
    stride_ = kl_+ku_+1;
    values_ = new ScalarType[stride_*numCols_];
    copyMat(Source.values_, Source.stride_, numRows_, numCols_, values_, stride_, 0);

  } else if ( trans == Teuchos::CONJ_TRANS && ScalarTraits<ScalarType>::isComplex ) {

    numRows_ = Source.numCols_;
    numCols_ = Source.numRows_;
    kl_ = Source.ku_;
    ku_ = Source.kl_;
    stride_ = kl_+ku_+1;
    values_ = new ScalarType[stride_*numCols_];
    for (OrdinalType j=0; j<numCols_; j++) {
      for (OrdinalType i=TEUCHOS_MAX(0,j-ku_); i<=TEUCHOS_MIN(numRows_-1,j+kl_); i++) {
	values_[j*stride_ + (ku_+i-j)] =
	  Teuchos::ScalarTraits<ScalarType>::conjugate(Source.values_[i*Source.stride_ + (Source.ku_+j-i)]);
      }
    }

  } else {

    numRows_ = Source.numCols_;
    numCols_ = Source.numRows_;
    kl_ = Source.ku_;
    ku_ = Source.kl_;
    stride_ = kl_+ku_+1;
    values_ = new ScalarType[stride_*numCols_];
    for (OrdinalType j=0; j<numCols_; j++) {
      for (OrdinalType i=TEUCHOS_MAX(0,j-ku_); i<=TEUCHOS_MIN(numRows_-1,j+kl_); i++) {
	values_[j*stride_ + (ku_+i-j)] = Source.values_[i*Source.stride_ + (Source.ku_+j-i)];
      }
    }

  }

}

template<typename OrdinalType, typename ScalarType>
SerialBandDenseMatrix<OrdinalType, ScalarType>::SerialBandDenseMatrix(
  DataAccess CV, const SerialBandDenseMatrix<OrdinalType, ScalarType> &Source,
  OrdinalType numRows_in, OrdinalType numCols_in, OrdinalType startCol )
  : CompObject(), Object("Teuchos::SerialBandDenseMatrix"), numRows_(numRows_in), numCols_(numCols_in), kl_(Source.kl_), ku_(Source.ku_), stride_(Source.stride_), valuesCopied_(false), values_(Source.values_)
{
  if(CV == Copy) {
      values_ = new ScalarType[stride_ * numCols_in];
      copyMat(Source.values_, Source.stride_, numRows_in, numCols_in, values_, stride_, startCol);
      valuesCopied_ = true;
  } else { // CV = View
    values_ = values_ + (stride_ * startCol);
  }
}

template<typename OrdinalType, typename ScalarType>
SerialBandDenseMatrix<OrdinalType, ScalarType>::~SerialBandDenseMatrix()
{
  deleteArrays();
}

//----------------------------------------------------------------------------------------------------
//  Shape methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseMatrix<OrdinalType, ScalarType>::shape(
  OrdinalType numRows_in, OrdinalType numCols_in, OrdinalType kl_in, OrdinalType ku_in
  )
{
  deleteArrays(); // Get rid of anything that might be already allocated
  numRows_ = numRows_in;
  numCols_ = numCols_in;
  kl_ = kl_in;
  ku_ = ku_in;
  stride_ = kl_+ku_+1;
  values_ = new ScalarType[stride_*numCols_];
  putScalar();
  valuesCopied_ = true;

  return(0);
}

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseMatrix<OrdinalType, ScalarType>::shapeUninitialized(
  OrdinalType numRows_in, OrdinalType numCols_in, OrdinalType kl_in, OrdinalType ku_in
  )
{
  deleteArrays(); // Get rid of anything that might be already allocated
  numRows_ = numRows_in;
  numCols_ = numCols_in;
  kl_ = kl_in;
  ku_ = ku_in;
  stride_ = kl_+ku_+1;
  values_ = new ScalarType[stride_*numCols_];
  valuesCopied_ = true;

  return(0);
}

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseMatrix<OrdinalType, ScalarType>::reshape(
  OrdinalType numRows_in, OrdinalType numCols_in, OrdinalType kl_in, OrdinalType ku_in
  )
{

  // Allocate space for new matrix
  ScalarType* values_tmp = new ScalarType[(kl_in+ku_in+1) * numCols_in];
  ScalarType zero = ScalarTraits<ScalarType>::zero();
  for(OrdinalType k = 0; k < (kl_in+ku_in+1) * numCols_in; k++) {
    values_tmp[k] = zero;
  }
  OrdinalType numRows_tmp = TEUCHOS_MIN(numRows_, numRows_in);
  OrdinalType numCols_tmp = TEUCHOS_MIN(numCols_, numCols_in);
  if(values_ != 0) {
    copyMat(values_, stride_, numRows_tmp, numCols_tmp, values_tmp,
	    kl_in+ku_in+1, 0); // Copy principal submatrix of A to new A
  }
  deleteArrays(); // Get rid of anything that might be already allocated
  numRows_ = numRows_in;
  numCols_ = numCols_in;
  kl_ = kl_in;
  ku_ = ku_in;
  stride_ = kl_+ku_+1;
  values_ = values_tmp; // Set pointer to new A
  valuesCopied_ = true;

  return(0);
}

//----------------------------------------------------------------------------------------------------
//   Set methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseMatrix<OrdinalType, ScalarType>::putScalar( const ScalarType value_in )
{

  // Set each value of the dense matrix to "value".
  for(OrdinalType j = 0; j < numCols_; j++) {
    for (OrdinalType i=TEUCHOS_MAX(0,j-ku_); i<=TEUCHOS_MIN(numRows_-1,j+kl_); i++) {
      values_[(ku_+i-j) + j*stride_] = value_in;
    }
  }

  return 0;
}

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseMatrix<OrdinalType, ScalarType>::random()
{

  // Set each value of the dense matrix to a random value.
  for(OrdinalType j = 0; j < numCols_; j++) {
    for (OrdinalType i=TEUCHOS_MAX(0,j-ku_); i<=TEUCHOS_MIN(numRows_-1,j+kl_); i++) {
      values_[(ku_+i-j) + j*stride_] = ScalarTraits<ScalarType>::random();
    }
  }

  return 0;
}

template<typename OrdinalType, typename ScalarType>
SerialBandDenseMatrix<OrdinalType,ScalarType>&
SerialBandDenseMatrix<OrdinalType, ScalarType>::operator=(
  const SerialBandDenseMatrix<OrdinalType,ScalarType>& Source
  )
{

  if(this == &Source)
    return(*this); // Special case of source same as target
  if((!valuesCopied_) && (!Source.valuesCopied_) && (values_ == Source.values_))
    return(*this); // Special case of both are views to same data.

  // If the source is a view then we will return a view, else we will return a copy.
  if (!Source.valuesCopied_) {
    if(valuesCopied_) {
      // Clean up stored data if this was previously a copy.
      deleteArrays();
    }
    numRows_ = Source.numRows_;
    numCols_ = Source.numCols_;
    kl_ = Source.kl_;
    ku_ = Source.ku_;
    stride_ = Source.stride_;
    values_ = Source.values_;
  } else {
    // If we were a view, we will now be a copy.
    if(!valuesCopied_) {
      numRows_ = Source.numRows_;
      numCols_ = Source.numCols_;
      kl_ = Source.kl_;
      ku_ = Source.ku_;
      stride_ = kl_+ku_+1;
      const OrdinalType newsize = stride_ * numCols_;
      if(newsize > 0) {
	values_ = new ScalarType[newsize];
	valuesCopied_ = true;
      } else {
	values_ = 0;
      }
    } else {
      // If we were a copy, we will stay a copy.
      if((Source.numRows_ <= stride_) && (Source.numCols_ == numCols_)) { // we don't need to reallocate
	numRows_ = Source.numRows_;
	numCols_ = Source.numCols_;
	kl_ = Source.kl_;
	ku_ = Source.ku_;
      } else {
	// we need to allocate more space (or less space)
	deleteArrays();
	numRows_ = Source.numRows_;
	numCols_ = Source.numCols_;
	kl_ = Source.kl_;
	ku_ = Source.ku_;
	stride_ = kl_+ku_+1;
	const OrdinalType newsize = stride_ * numCols_;
	if(newsize > 0) {
	  values_ = new ScalarType[newsize];
	  valuesCopied_ = true;
	}
      }
    }
    copyMat(Source.values_, Source.stride_, numRows_, numCols_, values_, stride_, 0);
  }
  return(*this);

}

template<typename OrdinalType, typename ScalarType>
SerialBandDenseMatrix<OrdinalType, ScalarType>& SerialBandDenseMatrix<OrdinalType, ScalarType>::operator+= (const SerialBandDenseMatrix<OrdinalType,ScalarType>& Source )
{

  // Check for compatible dimensions
  if ((numRows_ != Source.numRows_) || (numCols_ != Source.numCols_) || (kl_ != Source.kl_) || (ku_ != Source.ku_)) {
    TEUCHOS_CHK_REF(*this); // Return *this without altering it.
  }
  copyMat(Source.values_, Source.stride_, numRows_, numCols_, values_, stride_, 0, ScalarTraits<ScalarType>::one());
  return(*this);

}

template<typename OrdinalType, typename ScalarType>
SerialBandDenseMatrix<OrdinalType, ScalarType>& SerialBandDenseMatrix<OrdinalType, ScalarType>::operator-= (const SerialBandDenseMatrix<OrdinalType,ScalarType>& Source )
{

  // Check for compatible dimensions
  if ((numRows_ != Source.numRows_) || (numCols_ != Source.numCols_) || (kl_ != Source.kl_) || (ku_ != Source.ku_)) {
    TEUCHOS_CHK_REF(*this); // Return *this without altering it.
  }
  copyMat(Source.values_, Source.stride_, numRows_, numCols_, values_, stride_, 0, -ScalarTraits<ScalarType>::one());
  return(*this);

}

template<typename OrdinalType, typename ScalarType>
SerialBandDenseMatrix<OrdinalType,ScalarType>& SerialBandDenseMatrix<OrdinalType, ScalarType>::assign (const SerialBandDenseMatrix<OrdinalType,ScalarType>& Source) {

  if(this == &Source)
    return(*this); // Special case of source same as target
  if((!valuesCopied_) && (!Source.valuesCopied_) && (values_ == Source.values_))
    return(*this); // Special case of both are views to same data.

  // Check for compatible dimensions
  if ((numRows_ != Source.numRows_) || (numCols_ != Source.numCols_) || (kl_ != Source.kl_) || (ku_ != Source.ku_)) {
    TEUCHOS_CHK_REF(*this); // Return *this without altering it.
  }
  copyMat(Source.values_, Source.stride_, numRows_, numCols_, values_, stride_, 0);
  return(*this);

}

//----------------------------------------------------------------------------------------------------
//   Accessor methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
inline ScalarType& SerialBandDenseMatrix<OrdinalType, ScalarType>::operator () (OrdinalType rowIndex, OrdinalType colIndex)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  checkIndex( rowIndex, colIndex );
#endif
  return(values_[colIndex * stride_ + ku_+rowIndex-colIndex]);
}

template<typename OrdinalType, typename ScalarType>
inline const ScalarType& SerialBandDenseMatrix<OrdinalType, ScalarType>::operator () (OrdinalType rowIndex, OrdinalType colIndex) const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  checkIndex( rowIndex, colIndex );
#endif
  return(values_[colIndex * stride_ + ku_+rowIndex-colIndex]);
}

template<typename OrdinalType, typename ScalarType>
inline const ScalarType* SerialBandDenseMatrix<OrdinalType, ScalarType>::operator [] (OrdinalType colIndex) const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  checkIndex( 0, colIndex );
#endif
  return(values_ + colIndex * stride_);
}

template<typename OrdinalType, typename ScalarType>
inline ScalarType* SerialBandDenseMatrix<OrdinalType, ScalarType>::operator [] (OrdinalType colIndex)
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
typename ScalarTraits<ScalarType>::magnitudeType SerialBandDenseMatrix<OrdinalType, ScalarType>::normOne() const
{
  OrdinalType i, j;
  typename ScalarTraits<ScalarType>::magnitudeType anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
  typename ScalarTraits<ScalarType>::magnitudeType absSum = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());

  ScalarType* ptr;
  for(j = 0; j < numCols_; j++) {
    ScalarType sum = 0;
    ptr = values_ + j * stride_ + TEUCHOS_MAX(0, ku_-j);
    for (i=TEUCHOS_MAX(0,j-ku_); i<=TEUCHOS_MIN(numRows_-1,j+kl_); i++) {
      sum += ScalarTraits<ScalarType>::magnitude(*ptr++);
    }
    absSum = ScalarTraits<ScalarType>::magnitude(sum);
    if(absSum > anorm) {
      anorm = absSum;
    }
  }
  updateFlops((kl_+ku_+1) * numCols_);

  return(anorm);
}

template<typename OrdinalType, typename ScalarType>
typename ScalarTraits<ScalarType>::magnitudeType SerialBandDenseMatrix<OrdinalType, ScalarType>::normInf() const
{
  OrdinalType i, j;
  typename ScalarTraits<ScalarType>::magnitudeType sum, anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());

  for (i = 0; i < numRows_; i++) {
    sum = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
    for (j=TEUCHOS_MAX(0,i-kl_); j<=TEUCHOS_MIN(numCols_-1,i+ku_); j++) {
      sum += ScalarTraits<ScalarType>::magnitude(*(values_+(ku_+i-j)+j*stride_));
    }
    anorm = TEUCHOS_MAX( anorm, sum );
  }
  updateFlops((kl_+ku_+1) * numCols_);

  return(anorm);
}

template<typename OrdinalType, typename ScalarType>
typename ScalarTraits<ScalarType>::magnitudeType SerialBandDenseMatrix<OrdinalType, ScalarType>::normFrobenius() const
{
  OrdinalType i, j;
  typename ScalarTraits<ScalarType>::magnitudeType anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());

  for (j = 0; j < numCols_; j++) {
    for (i=TEUCHOS_MAX(0,j-ku_); i<=TEUCHOS_MIN(numRows_-1,j+kl_); i++) {
      anorm += ScalarTraits<ScalarType>::magnitude(values_[(ku_+i-j)+j*stride_]*values_[(ku_+i-j)+j*stride_]);
    }
  }
  anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::squareroot(anorm));
  updateFlops((kl_+ku_+1) * numCols_);

  return(anorm);
}

//----------------------------------------------------------------------------------------------------
//   Comparison methods
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
bool SerialBandDenseMatrix<OrdinalType, ScalarType>::operator== (const SerialBandDenseMatrix<OrdinalType, ScalarType> &Operand) const
{
  bool result = 1;

  if((numRows_ != Operand.numRows_) || (numCols_ != Operand.numCols_) || (kl_ != Operand.kl_) || (ku_ != Operand.ku_)) {
    result = 0;
  } else {
    OrdinalType i, j;
    for(j = 0; j < numCols_; j++) {
      for (i=TEUCHOS_MAX(0,j-ku_); i<=TEUCHOS_MIN(numRows_-1,j+kl_); i++) {
	if((*this)(i, j) != Operand(i, j)) {
	  return 0;
	}
      }
    }
  }

  return result;
}

template<typename OrdinalType, typename ScalarType>
bool SerialBandDenseMatrix<OrdinalType, ScalarType>::operator!= (const SerialBandDenseMatrix<OrdinalType, ScalarType> &Operand) const
{
  return !((*this) == Operand);
}

//----------------------------------------------------------------------------------------------------
//   Multiplication method
//----------------------------------------------------------------------------------------------------

template<typename OrdinalType, typename ScalarType>
SerialBandDenseMatrix<OrdinalType, ScalarType>& SerialBandDenseMatrix<OrdinalType, ScalarType>::operator*= (const ScalarType alpha )
{
  this->scale( alpha );
  return(*this);
}

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseMatrix<OrdinalType, ScalarType>::scale( const ScalarType alpha )
{

  OrdinalType i, j;
  ScalarType* ptr;

  for (j=0; j<numCols_; j++) {
    ptr = values_ + j*stride_ + TEUCHOS_MAX(0, ku_-j);
    for (i=TEUCHOS_MAX(0,j-ku_); i<=TEUCHOS_MIN(numRows_-1,j+kl_); i++) {
      *ptr = alpha * (*ptr); ptr++;
    }
  }
  updateFlops( (kl_+ku_+1)*numCols_ );

  return(0);
}

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseMatrix<OrdinalType, ScalarType>::scale( const SerialBandDenseMatrix<OrdinalType,ScalarType>& A )
{

  OrdinalType i, j;
  ScalarType* ptr;

  // Check for compatible dimensions
  if ((numRows_ != A.numRows_) || (numCols_ != A.numCols_) || (kl_ != A.kl_) || (ku_ != A.ku_)) {
    TEUCHOS_CHK_ERR(-1); // Return error
  }
  for (j=0; j<numCols_; j++) {
    ptr = values_ + j*stride_ + TEUCHOS_MAX(0, ku_-j);
    for (i=TEUCHOS_MAX(0,j-ku_); i<=TEUCHOS_MIN(numRows_-1,j+kl_); i++) {
      *ptr = A(i,j) * (*ptr); ptr++;
    }
  }
  updateFlops( (kl_+ku_+1)*numCols_ );

  return(0);
}

template<typename OrdinalType, typename ScalarType>
void SerialBandDenseMatrix<OrdinalType, ScalarType>::print(std::ostream& os) const
{
  os << std::endl;
  if(valuesCopied_)
    os << "Values_copied : yes" << std::endl;
  else
    os << "Values_copied : no" << std::endl;
  os << "Rows : " << numRows_ << std::endl;
  os << "Columns : " << numCols_ << std::endl;
  os << "Lower Bandwidth : " << kl_ << std::endl;
  os << "Upper Bandwidth : " << ku_ << std::endl;
  os << "LDA : " << stride_ << std::endl;
  if(numRows_ == 0 || numCols_ == 0) {
    os << "(matrix is empty, no values to display)" << std::endl;
  } else {

    for(OrdinalType i = 0; i < numRows_; i++) {
      for (OrdinalType j=TEUCHOS_MAX(0,i-kl_); j<=TEUCHOS_MIN(numCols_-1,i+ku_); j++) {
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
inline void SerialBandDenseMatrix<OrdinalType, ScalarType>::checkIndex( OrdinalType rowIndex, OrdinalType colIndex ) const {

  TEUCHOS_TEST_FOR_EXCEPTION(rowIndex < 0 || rowIndex >= numRows_ ||
			     rowIndex < TEUCHOS_MAX(0,colIndex-ku_) || rowIndex > TEUCHOS_MIN(numRows_-1,colIndex+kl_),
			     std::out_of_range,
			     "SerialBandDenseMatrix<T>::checkIndex: "
			     "Row index " << rowIndex << " out of range [0, "<< numRows_ << ")");
  TEUCHOS_TEST_FOR_EXCEPTION(colIndex < 0 || colIndex >= numCols_, std::out_of_range,
			     "SerialBandDenseMatrix<T>::checkIndex: "
			     "Col index " << colIndex << " out of range [0, "<< numCols_ << ")");

}

template<typename OrdinalType, typename ScalarType>
void SerialBandDenseMatrix<OrdinalType, ScalarType>::deleteArrays(void)
{
  if (valuesCopied_) {
    delete [] values_;
    values_ = 0;
    valuesCopied_ = false;
  }
}

template<typename OrdinalType, typename ScalarType>
void SerialBandDenseMatrix<OrdinalType, ScalarType>::copyMat(
  ScalarType* inputMatrix, OrdinalType strideInput, OrdinalType numRows_in,
  OrdinalType numCols_in, ScalarType* outputMatrix, OrdinalType strideOutput, OrdinalType startCol, ScalarType alpha
  )
{
  OrdinalType i, j;
  ScalarType* ptr1 = 0;
  ScalarType* ptr2 = 0;

  for(j = 0; j < numCols_in; j++) {
    ptr1 = outputMatrix + (j * strideOutput) + TEUCHOS_MAX(0, ku_-j);
    ptr2 = inputMatrix + (j + startCol) * strideInput + TEUCHOS_MAX(0, ku_-j);
    if (alpha != Teuchos::ScalarTraits<ScalarType>::zero() ) {
      for (i=TEUCHOS_MAX(0,j-ku_); i<=TEUCHOS_MIN(numRows_in-1,j+kl_); i++) {
	*ptr1++ += alpha*(*ptr2++);
      }
    } else {
      for (i=TEUCHOS_MAX(0,j-ku_); i<=TEUCHOS_MIN(numRows_in-1,j+kl_); i++) {
	*ptr1++ = *ptr2++;
      }
    }
  }
}

} // namespace Teuchos


#endif /* _TEUCHOS_SERIALBANDDENSEMATRIX_HPP_ */
