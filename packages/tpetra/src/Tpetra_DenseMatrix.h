// 27-May-2002 General cleanup. Changed method names to fit namingConvention (already done).

#ifndef _TPETRA_DENSEMATRIX_H_
#define _TPETRA_DENSEMATRIX_H_

#include "Tpetra_CompObject.h"
#include "Tpetra_BLAS.h"

namespace Tpetra
{

//! Tpetra::DenseMatrix: A class for constructing and using template<scalarType> general dense matrices.

/*! The Tpetra::DenseMatrix class enables the construction and use of general template<scalarType>
    dense matrices.  It is built on the BLAS, and derives from the TPetra_BLAS.  Of course, it's also built
    on ScalarTraits, which derives from LAPACK. In any case, 

The Tpetra::DenseMatrix class is intended to provide very basic support for dense rectangular matrices.


<b>Constructing Tpetra::DenseMatrix Objects</b>

There are three Tpetra::DenseMatrix constructors.  The first constructs a zero-sized object which should be made
to appropriate length using the Shape() or Reshape() functions and then filled with the [] or () operators. 
The second is a constructor that accepts user
data as a 2D array, the third is a copy constructor. The second constructor has
two data access modes (specified by the Petra_DataAccess argument):
<ol>
  <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
       user data is not needed after construction.
  <li> View mode - Creates a "view" of the user data. In this case, the
       user data is required to remain intact for the life of the object.
</ol>

\warning View mode is \e extremely dangerous from a data hiding perspective.
Therefore, we strongly encourage users to develop code using Copy mode first and 
only use the View mode in a secondary optimization phase.

<b>Extracting Data from Tpetra::DenseMatrix Objects</b>

Once a Tpetra::DenseMatrix is constructed, it is possible to view the data via access functions.

\warning Use of these access functions cam be \e extremely dangerous from a data hiding perspective.


<b>Vector and Utility Functions</b>

Once a Tpetra::DenseMatrix is constructed, several mathematical functions can be applied to
the object.  Specifically:
<ul>
  <li> Multiplication.
  <li> Norms.
</ul>

The final useful function is flops().  Each Tpetra::DenseMatrix object keep track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The flops() function returns this number as a double precision number.  Using this 
information, in conjunction with the Tpetra_Time class, one can get accurate parallel performance
numbers.


*/
//=========================================================================
template<class scalarType>
class DenseMatrix : public CompObject, public Object, public BLAS<scalarType>
{

  public:
  
  //@{ \name Constructor/Destructor Methods
  //! Default constructor; defines a zero size object.
  /*!
    Tpetra::DenseMatrix objects defined by the default constructor should be sized with the 
    Shape() or Reshape functions.  
    Values should be defined by using the [] or () operators.
   */
  DenseMatrix(void);
  
  //! Set object values from two-dimensional array.
  /*!
    \param In 
           Petra_DataAccess - Enumerated type set to Copy or View.
    \param In
           A - Pointer to an array of template<scalarType> numbers.  The first vector starts at A.
	   The second vector starts at A+stride, the third at A+2*stride, and so on.
    \param In
           stride - The "Leading Dimension", or stride between vectors in memory.
    \param In 
           NumRows - Number of rows in object.
    \param In 
           NumCols - Number of columns in object.

	   See Detailed Description section for further discussion.
  */
  DenseMatrix(Tpetra_DataAccess CV, scalarType *A, int stride, int numRows, int numCols);
  
  //! Tpetra::DenseMatrix copy constructor.
  
  DenseMatrix(const DenseMatrix<scalarType>& Source);

  //! DenseMatrix destructor.  
  virtual ~DenseMatrix ();
  //@}

  //@{ \name Shaping/sizing Methods
  //! Set dimensions of a Tpetra::DenseMatrix object; init values to zero.
  /*!
    \param In 
           numRows - number of rows in object.
    \param In 
           numCols - number of columns in object.

	   Allows user to define the dimensions of a Tpetra::DenseMatrix at any point. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   destroyed and the resized matrix starts off with all zero values.

    \return Integer error code, set to 0 if successful.
  */
  int shape(int numRows, int numCols);
  
  //! Reshape a Tpetra::DenseMatrix object.
  /*!
    \param In 
           numRows - number of rows in object.
    \param In 
           numCols - number of columns in object.

	   Allows user to define the dimensions of a Tpetra::DenseMatrix at any point. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   copied into the new shape.  If the new shape is smaller than the original, the upper left portion
	   of the original matrix (the principal submatrix) is copied to the new matrix.

    \return Integer error code, set to 0 if successful.
  */
  int reshape(int numRows, int numCols);
  //@}
  
  //! Computes the 1-Norm of the \e this matrix.
  /*!
    \return Integer error code, set to 0 if successful.
  */
  ScalarTraits<scalarType>::magnitudeType oneNorm();

  //! Matrix-Matrix multiplication, \e this = Scalar*\e this + ScalarAB*A*B.
  /*! This function performs a variety of matrix-matrix multiply operations.

  \param In
         TransA - Operate with the transpose of A if = 'T', else no transpose if = 'N'.
  \param In
         TransB - Operate with the transpose of B if = 'T', else no transpose if = 'N'.

  \param In
         ScalarAB - Scalar to multiply with A*B.
  \param In
         A - Dense Matrix.
  \param In
         B - Dense Matrix.
  \param In
         Scalar - Scalar to multiply with \e this.

    \return Integer error code, set to 0 if successful.
	 
  */
  int  multiply (char TransA, char TransB, scalarType ScalarAB, 
                 const DenseMatrix<scalarType>& A, 
                 const DenseMatrix<scalarType>& B,
                 scalarType Scalar );

  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Element from the specified row and column.
  */
    scalarType &operator () (int RowIndex, int ColIndex);

  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Element from the specified row and column.
  */
    const scalarType &operator () (int RowIndex, int ColIndex) const;

  //! Column access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Pointer to address of specified column.

    \warning No bounds checking can be done for the index i in the expression A[j][i].
  */
    scalarType *operator [] (int ColIndex);

  //! Column access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Pointer to address of specified column.

    \warning No bounds checking can be done for the index i in the expression A[j][i].
  */
    const scalarType *operator [] (int ColIndex) const;
    
  //! Returns row dimension of system.
  int numRows() const
  {
    return(numRows_);
  };

  //! Returns column dimension of system.
  int numCols() const 
  {
    return(numCols_);
  };

  //! Returns pointer to the \e this matrix.
  scalarType * values() const 
  {
    return(values_);
  };

  //! Returns the leading dimension of the \e this matrix.
  int stride() const 
  {
    return(stride_);
  };

  virtual void print(ostream& os) const;

 protected:

  void copyMat(scalarType * inputMatrix, int strideInput, int numRows, int numCols, scalarType * outputMatrix, int strideOutput);
  void deleteArrays(void);

  int numRows_;
  int numCols_;
  int stride_;
  bool valuesCopied_;
  scalarType * values_;


}; // class Tpetra_DenseMatrix

} // namespace Tpetra

template<class scalarType>
inline ostream& operator<<(ostream& os, const Tpetra::DenseMatrix<scalarType>& obj)
{
  // os << obj.Label();
  obj.print(os);
  return os;
}

#include "Tpetra_DenseMatrix.cpp"

#endif /* _TPETRA_DENSEMATRIX_H_ */
