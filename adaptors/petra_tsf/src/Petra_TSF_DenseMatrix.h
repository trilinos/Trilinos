#ifndef _PETRA_TSF_DENSEMATRIX_H_
#define _PETRA_TSF_DENSEMATRIX_H_

#include "Petra_RDP_MultiVector.h"
#include "Petra_Map.h"
#include "TSF_DenseMatrix.h"

namespace Petra_TSF {

//! The TSF::DenseMatrix class provides basic support for dense rectangular matrices with elements <scalarType>.

template<class scalarType>
class DenseMatrix public virtual TSF::DenseMatrix, public Petra_BLAS_DGE_Matrix  {

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
  virtual int shape(int numRows, int numCols){PETRA_CHK_ERR(Shape(numRows,numCols));};
  
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
  virtual int reshape(int numRows, int numCols){PETRA_CHK_ERR(Reshape(numRows,numCols));};

  //! Returns number of rows in matrix.
  virtual int numRows(){return(M());};

  //! Returns number of columns in matrix.
  virtual int numCols(){return(M());};
  //@}
  
  //@{ \name Strided matrix accessor methods.
  //! Strided matrix query; returns true if matrix has provided strided memory access.
  /*! The Petra_BLAS_DGE_Matrix class supports strided matrix data storage references.  This query will always
      be true.
  */
  virtual bool isStrided(){return(true);};


  //! Return location of the first element in the first column of this dense matrix.
  virtual scalarType * & getValues() {return(A());};

  //! Return location of (0,0) entry in this dense matrix (non-const version).
  virtual scalarType * & getValues() {return(A());};

  //! Returns the so-called "leading dimension" of matrix.
  /*! This method returns the distance between columns of the dense matrix.  Thus, if the
      first element in the first column is pointed to by p = getValues(), the first element
      in the second column is at p[stride()].
  */
  virtual int stride() {return(LDA());};
  //@}

  //@{ \name Mathematical methods

  //! Computes the 1-Norm of the \e this matrix.
  /*!
    \return Integer error code, set to 0 if successful.
  */
  virtual scalarType norm1() {return(OneNorm());};

  //! Computes the inf-Norm of the \e this matrix.
  /*!
    \return Integer error code, set to 0 if successful.
  */
  virtual scalarType NormInf() {return(InfNorm());};

  //! Matrix-Matrix multiplication, C = Scalar*C + ScalarAB*A*B.
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
                 scalarType scalar) {
    const Petra_RDP_MultiVector *A1 = dynamic_cast<const Petra_BLAS_DGE_Matrix *>(&A);
    if (A1==0) PETRA_CHK_ERR(-1);
    const Petra_RDP_MultiVector *B1 = dynamic_cast<const Petra_BLAS_DGE_Matrix *>(&B);
    if (B1==0) PETRA_CHK_ERR(-2);
    PETRA_CHK_ERR(Multiply(transA, transB, scalarAB, *A1, *B1, scalar));
  };
  //@}

  //@{ \name Element access methods
  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Element from the specified row and column.
  */
  virtual scalarType& operator () (int rowIndex, int colIndex) {return((*this)(rowIndex,colIndex));};

  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Element from the specified row and column.
  */
  virtual const scalarType& operator () (int rowIndex, int colIndex) const {return((*this)(rowIndex,colIndex));};
  //@}

};

// Constructor implementations

//=============================================================================
// Petra_BlockMap Constructor

template<class scalarType>
DenseMatrix<scalarType>::DenseMatrix()
  : Petra_BLAS_DGE_Matrix(){}

//==========================================================================
// Copy Constructor


template<class scalarType>
DenseMatrix<scalarType>::DenseMatrix(const Petra_TSF::DenseMatrix<scalarType>& Source)
  :Petra_RDP_DenseMatrix(Source) {}

//==========================================================================
// This constructor copies in or makes view of a standard Fortran array

template<>
DenseMatrix<double>::DenseMatrix(Petra_DataAccess CV, double *A, int LDA, int NumRows, int NumCols)
  : Petra_RDP_DenseMatrix(CV, A, LDA, NumRows, NumCols) {}

} // Petra_TSF namespace
#endif /* _TSF_DENSEMATRIX_H_ */
