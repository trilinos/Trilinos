/*Paul
// 27-May-2002 General cleanup. Checked for newNamingConvention (nothing changed).
06-August-2002 Changed to images (nothing changed).
*/


#ifndef _TPETRA_BLAS_H_
#define _TPETRA_BLAS_H_

namespace Tpetra
{
//! Tpetra_BLAS: The Templated Petra BLAS Class.
/*! The Tpetra_BLAS class provides functionality similar to the BLAS
    (Basic Linear Algebra Subprograms).  The BLAS provide portable, high-
    performance implementations of kernels such as dense vectoer multiplication,
    dot products, dense matrix-vector multiplication and dense matrix-matrix
    multiplication.

    The standard BLAS interface is Fortran-specific.  Unfortunately, the 
    interface between C++ and Fortran is not standard across all computer
    platforms.  The Tpetra_BLAS class provides C++ bindings for the BLAS
    kernels in order to insulate the rest of Petra from the details of 
    C++ to Fortran translation.

    In addition to giving access the standard BLAS functionality.
    Tpetra_BLAS also provide functionality for any <scalarType> class that
    defines the +, - * and / operators.

    Tpetra_BLAS is a single memory image interface only.  This is appropriate 
    since the standard BLAS are only specified for serial execution 
    (or shared memory parallel).
*/

  template<class scalarType>
  class BLAS
  {
    public:
    
    //! Tpetra::BLAS Constructor.
    /*! Builds an instance of a serial BLAS object.
    */
    BLAS(void);
    
    //! Tpetra::BLAS Copy Constructor.
    /*! Makes an exact copy of an existing Tpetra::BLAS instance.
    */
    BLAS(const BLAS<scalarType>& BLAS);
    
    //! Tpetra::BLAS Destructor.
    virtual ~BLAS(void);
    
    //! Tpetra::BLAS matrix-matrix multiply function (GEMM)
    void GEMM(char TRANSA, char TRANSB, int M, int N, int K, scalarType ALPHA, 
        scalarType * A, int LDA, scalarType * B, int LDB, scalarType BETA, 
        scalarType * C, int LDC) const;
  };
} // end of namespace Tpetra

#include "Tpetra_BLAS.cpp"

#endif // end of _TPETRA_BLAS_H_
