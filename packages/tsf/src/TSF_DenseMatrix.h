#ifndef _TSF_DENSEMATRIX_H_
#define _TSF_DENSEMATRIX_H_

#include "TSF_Object.h"

namespace TSF {

//! The TSF::DenseMatrix class provides basic support for dense rectangular matrices with elements <scalarType>.

template<class scalarType>
class DenseMatrix  {

  public:
  

  //! TSF::DenseMatrix destructor.  
  virtual ~DenseMatrix (void){};
  //@}

  //@{ \name Shaping/sizing Methods
  //! Set dimensions of a TSF::DenseMatrix object; init values to zero.
  /*!
    \param In 
           numRows - Number of rows in object.
    \param In 
           numCols - Number of columns in object.

	   Allows user to define the dimensions of a TSF::DenseMatrix at any point. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   destroyed and the resized matrix starts off with all zero values.

    \return Integer error code, set to 0 if successful.
  */
  virtual int shape(int numRows, int numCols) = 0;
  
  //! Reshape a TSF::DenseMatrix object.
  /*!
    \param In 
           numRows - Number of rows in object.
    \param In 
           numCols - Number of columns in object.

	   Allows user to define the dimensions of a TSF::DenseMatrix at any point. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   copied into the new shape.  If the new shape is smaller than the original, the upper left portion
	   of the original matrix (the principal submatrix) is copied to the new matrix.

    \return Integer error code, set to 0 if successful.
  */
  virtual int reshape(int numRows, int numCols) = 0;

  //! Returns number of rows in matrix.
  virtual int numRows() = 0;

  //! Returns number of columns in matrix.
  virtual int numCols() = 0;
  //@}
  
  //@{ \name Strided matrix accessor methods.
  //! Strided matrix query; returns true if matrix has provided strided memory access.
  /*! The methods in this section depend on the matrix being strided.  This corresponds
      to a Fortran-style 2D array of values.  If the class that implements this interface
      supports a strided data structure, the isStrided() method should be defined to return true, 
      otherwise it should be defined to return false.
  */
  virtual bool isStrided() const = 0;

  //! Return location of (0,0) entry in this dense matrix.
  /*! Note:  By requiring the following three methods we are necessarily restricting our dense
             matrix data structure to provide a Fortran-style 2D copy of itself.
	     This restriction is purposefully done because so much dense matrix linear algebra
	     is provided for this data structure.
  */
  virtual scalarType * & getValues() const = 0;

  //! Return location of (0,0) entry in this dense matrix (non-const version).
  virtual scalarType * & getValues() = 0;

  //! Returns the so-called "leading dimension" of matrix.
  /*! This method returns the distance between columns of the dense matrix.  Thus, if the
      first element in the first column is pointed to by p = getValues(), the first element
      in the second column is at p[getStride()].
  */
  virtual int getStride() = 0;
  //@}

  //@{ \name Mathematical methods
  //! Computes the 1-Norm of the \e this matrix.
  /*!
    \return Integer error code, set to 0 if successful.
  */
  virtual scalarType norm1() = 0;

  //! Computes the inf-Norm of the \e this matrix.
  /*!
    \return Integer error code, set to 0 if successful.
  */
  virtual scalarType normInf() = 0;

  //! Matrix-Matrix multiplication, C = ScalarC*C + ScalarAB*A*B.
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
         ScalarC - Scalar to multiply with C.
  \param InOut
         C - Dense Matrix.

    \return Integer error code, set to 0 if successful.
	 
  */
  virtual int  multiply (char transA, char transB, scalarType scalarAB, 
                 const TSF::DenseMatrix& A, 
                 const TSF::DenseMatrix& B,
                 scalarType scalar ) = 0;
  //@}

  //@{ \name Element access methods
  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Element from the specified row and column.
  */
  virtual scalarType& operator () (int rowIndex, int colIndex) = 0;

  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Element from the specified row and column.
  */
  virtual const scalarType& operator () (int rowIndex, int colIndex) const = 0;
  //@}
};
} // TSF namespace
#endif /* _TSF_DENSEMATRIX_H_ */
