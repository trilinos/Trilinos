#ifndef _TPETRA_BLAS_H_
#define _TPETRA_BLAS_H_

//! TPetra_BLAS:  The Petra Templated BLAS Class.
/*! The TPetra_BLAS class provides functionality similar to the BLAS
    (Basic Linear Algebra Subprograms).  The BLAS provide portable, high-
    performance implementations of kernels such as dense vectoer multiplication,
    dot products, dense matrix-vector multiplication and dense matrix-matrix
    multiplication.

    The standard BLAS interface is Fortran-specific.  Unfortunately, the 
    interface between C++ and Fortran is not standard across all computer
    platforms.  The TPetra_BLAS class provides C++ bindings for the BLAS
    kernels in order to insulate the rest of Petra from the details of 
    C++ to Fortran translation.

    In addition to giving access the standard BLAS functionality.
    TPetra_BLAS also provide functionality for any <scalarType> class that
    defines the +, - * and / operators.

    TPetra_BLAS is a single memory image interface only.  This is appropriate 
    since the standard 
    BLAS are only specified for serial execution (or shared memory parallel).
*/


namespace TPetra {

template<class scalarType>
class BLAS {
    
  public:
  //! TPetra::BLAS Constructor.
  /*! Builds an instance of a serial BLAS object.
   */
  BLAS(void);


  //! TPetra::BLAS Copy Constructor.
  /*! Makes an exact copy of an existing TPetra::BLAS instance.
  */
  BLAS(const TPetra::BLAS<scalarType>& BLAS);

  //! TPetra::BLAS Destructor.
  virtual ~BLAS(void);
  
  //! TPetra::BLAS matrix-matrix multiply function (GEMM)
  void GEMM(char TRANSA, char TRANSB, int M, int N, int K,
	    scalarType ALPHA, scalarType * A, int LDA, scalarType * B,
	    int LDB, scalarType BETA, scalarType * C, int LDC) const;
};

} // namespace TPetra
#include "TPetra_BLAS.cpp"
#endif /* _TPETRA_BLAS_H_ */
