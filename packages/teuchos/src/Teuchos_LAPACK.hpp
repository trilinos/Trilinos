// Kris
// 06.11.03 -- Format cleanup
// 06.13.03 -- Initial templatization (Tpetra_LAPACK.cpp is no longer needed)
// 06.17.03 -- Changed LAPACK wrapper calls from XXXXX_F77() to xxxxx_()
//          -- Added warning messages to default function calls
//          -- Added LAPY2 and GEES by request
// 06.18.03 -- Renamed GEES parameters for default implementation
//          -- Changed LAPACK wrapper calls back to XXXXX_F77() from xxxxx_() (Oops!)
//          -- Streamlined character macro stuff
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_LAPACK_HPP_
#define _TEUCHOS_LAPACK_HPP_

/* For INTEL_CXML, the second arg may need to be changed to 'one'.  If so
the appropriate declaration of one will need to be added back into
functions that include the macro: */

#if defined (INTEL_CXML)  
  unsigned int one = 1;
#endif

#ifdef CHAR_MACRO
#undef CHAR_MACRO
#endif
#if defined (INTEL_CXML)
#define CHAR_MACRO(char_var) &char_var, one
#else
#define CHAR_MACRO(char_var) &char_var
#endif

#include "Teuchos_Object.cpp"
#include "Teuchos_LAPACK_wrappers.hpp"

namespace Teuchos
{

  // BEGIN GENERAL TEMPLATE DECLARATION //

  /**
   * Templated LAPACK wrapper class. Somebody needs to document the methods...
   */
  template<typename OrdinalType, typename ScalarType>
  class LAPACK
  {    
  public:
    inline LAPACK(void) {};
    inline LAPACK(const LAPACK& LAPACK) {};
    inline virtual ~LAPACK(void) {};
    void POTRF(char, int, ScalarType*, int, int*) const;
    void POTRS(char, int, int, ScalarType*, int, ScalarType*, int, int*) const;
    void POTRI(char, int, ScalarType*, int, int*) const;
    void POCON(char, int, ScalarType*, int, ScalarType, ScalarType*, ScalarType*, int*, int*) const;
    void POSV(char, int, int, ScalarType*, int, ScalarType*, int, int*) const;
    void POEQU(int, ScalarType*, int, ScalarType*, ScalarType*, ScalarType*, int*) const;
    void PORFS(char, int, int, ScalarType*, int, ScalarType*, int, ScalarType*, int, ScalarType*, int, ScalarType*, ScalarType*, ScalarType*, int*, int*) const;
    void POSVX(char, char, int, int, ScalarType*, int, ScalarType*, int, char, ScalarType*, ScalarType*, int, ScalarType*, int, ScalarType*, ScalarType*, ScalarType*, ScalarType*, int*, int*) const;
    void GELS(char, int, int, int, ScalarType*, int, ScalarType*, int, ScalarType*, int, int) const;
    void GETRF(int, int, ScalarType*, int, int*, int*) const;
    void GETRS(char, int, int, ScalarType*, int, int*, ScalarType*, int, int*) const;
    void GETRI(int, ScalarType*, int, int*, ScalarType*, int*, int*) const;
    void GECON(char, int, ScalarType*, int, ScalarType, ScalarType*, ScalarType*, int*, int*) const;
    void GESV(int, int, ScalarType*, int, int*, ScalarType*, int, int*) const;
    void GEEQU(int, int, ScalarType*, int, ScalarType*, ScalarType*, ScalarType*, ScalarType*, ScalarType*, int*) const;
    void GERFS(char, int, int, ScalarType*, int, ScalarType*, int, int*, ScalarType*, int, ScalarType*, int, ScalarType*, ScalarType*, ScalarType*, int*, int*) const;
    void GESVX(char, char, int, int, ScalarType*, int, ScalarType*, int, int*, char, ScalarType*, ScalarType*, ScalarType*, int, ScalarType*, int, ScalarType*, ScalarType*, ScalarType*, ScalarType*, int*, int*) const;
    void GEHRD(int, int, int, ScalarType*, int, ScalarType*, ScalarType*, int, int*) const;
    void HSEQR(char, char, int, int, int, ScalarType*, int, ScalarType*, ScalarType*, ScalarType*, int, ScalarType*, int, int*) const;
    void ORGHR(int, int, int, ScalarType*, int, ScalarType*, ScalarType*, int, int*) const;
    void ORMHR(char, char, int, int, int, int, ScalarType*, int, ScalarType*, ScalarType*, int, ScalarType*, int, int*) const;
    void TREVC(char, char, int*, int, ScalarType*, int, ScalarType*, int, ScalarType*, int, int, int*, ScalarType*, int*) const;
    void TREXC(char, int, ScalarType*, int, ScalarType*, int, int, int, ScalarType*, int*) const;
    ScalarType LAMCH(char) const;
    void GEES(char, char, int*, int, ScalarType*, int, int, ScalarType*, ScalarType*, ScalarType*, int, ScalarType*, int, int*, int*);
    ScalarType LAPY2(ScalarType, ScalarType) const;
  };

  // END GENERAL TEMPLATE DECLARATION //

  // BEGIN GENERAL TEMPLATE IMPLEMENTATION //

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POTRF(char UPLO, int N, ScalarType* A, int LDA, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::POTRF() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POTRS(char UPLO, int N, int NRHS, ScalarType* A, int LDA, ScalarType* X, int LDX, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::POTRS() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POTRI(char UPLO, int N, ScalarType* A, int LDA, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::POTRI() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POCON(char UPLO, int N, ScalarType* A, int LDA, ScalarType ANORM, ScalarType* RCOND, ScalarType* WORK, int* IWORK, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::POCON() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POSV(char UPLO, int N, int NRHS, ScalarType* A, int LDA, ScalarType* X, int LDX, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::POSV() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POEQU(int N, ScalarType* A, int LDA, ScalarType* S, ScalarType* SCOND, ScalarType* AMAX, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::POEQU() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::PORFS(char UPLO, int N, int NRHS, ScalarType* A, int LDA, ScalarType* AF, int LDAF, ScalarType* B, int LDB, ScalarType* X, int LDX, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, int* IWORK, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::PORFS() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POSVX(char FACT, char UPLO, int N, int NRHS, ScalarType* A, int LDA, ScalarType* AF, int LDAF, char EQUED, ScalarType* S, ScalarType* B, int LDB, ScalarType* X, int LDX, ScalarType* RCOND, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, int* IWORK, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::POSVX() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GELS(char TRANS, int m, int n, int numrhs, ScalarType* a, int lda, ScalarType* b, int ldb, ScalarType* work, int lwork, int info) const
  {
    std::cout << "Warning: default LAPACK::GELS() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GETRF(int M, int N, ScalarType* A, int LDA, int* IPIV, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::GETRF() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GETRS(char TRANS, int N, int NRHS, ScalarType* A, int LDA, int* IPIV, ScalarType* X, int LDX, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::GETRS() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GETRI(int N, ScalarType* A, int LDA, int* IPIV, ScalarType* WORK, int* LWORK, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::GETRI() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GECON(char NORM, int N, ScalarType* A, int LDA, ScalarType ANORM, ScalarType* RCOND, ScalarType* WORK, int* IWORK, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::GECON() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GESV(int N, int NRHS, ScalarType* A, int LDA, int* IPIV, ScalarType* X, int LDX, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::GESV() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GEEQU(int M, int N, ScalarType* A, int LDA, ScalarType* R, ScalarType* C, ScalarType* ROWCND, ScalarType* COLCND, ScalarType* AMAX, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::GEEQU() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GERFS(char TRANS, int N, int NRHS, ScalarType* A, int LDA, ScalarType* AF, int LDAF, int* IPIV, ScalarType* B, int LDB, ScalarType* X, int LDX, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, int* IWORK, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::GERFS() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GESVX(char FACT, char TRANS, int N, int NRHS, ScalarType* A, int LDA, ScalarType* AF, int LDAF, int* IPIV, char EQUED, ScalarType* R, ScalarType* C, ScalarType* B, int LDB, ScalarType* X, int LDX, ScalarType* RCOND, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, int* IWORK, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::GESVX() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GEHRD(int N, int ILO, int IHI, ScalarType* A, int LDA, ScalarType* TAU, ScalarType* WORK, int LWORK, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::GEHRD() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::HSEQR(char JOB, char COMPZ, int N, int ILO, int IHI, ScalarType* H, int LDH, ScalarType* WR, ScalarType* WI, ScalarType* Z, int LDZ, ScalarType* WORK, int LWORK, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::HSEQR() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::ORGHR(int N, int ILO, int IHI, ScalarType* A, int LDA, ScalarType* TAU, ScalarType* WORK, int LWORK, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::ORGHR() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::ORMHR(char SIDE, char TRANS, int M, int N, int ILO, int IHI, ScalarType* A, int LDA, ScalarType* TAU, ScalarType* C, int LDC, ScalarType* WORK, int LWORK, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::ORMHR() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::TREVC(char SIDE, char HOWMNY, int* SELECT, int N, ScalarType* T, int LDT, ScalarType *VL, int LDVL, ScalarType* VR, int LDVR, int MM, int* M, ScalarType* WORK, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::TREVC() not yet implemented" << std::endl;
//     if(HOWMNY=='S')
//       {
// 	*INFO = -3; // We do not support 'S' since it requires a logical array
//       }
//     else
//       {

//       }
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::TREXC(char COMPQ, int N, ScalarType* T, int LDT, ScalarType* Q, int LDQ, int IFST, int ILST, ScalarType* WORK, int* INFO) const
  {
    std::cout << "Warning: default LAPACK::TREXC() not yet implemented" << std::endl;
  }

  template<typename OrdinalType, typename ScalarType>
  ScalarType LAPACK<OrdinalType, ScalarType>::LAMCH(char CMACH) const
  {
    std::cout << "Warning: default LAPACK::LAMCH() not yet implemented" << std::endl;
    ScalarType dummy;
    return dummy;
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GEES(char JOBVS, char SORT, int* SELECT, int N, ScalarType* A, int LDA, int SDIM, ScalarType* WR, ScalarType* WI, ScalarType* VS, int LDVS, ScalarType* WORK, int LWORK, int* BWORK, int* INFO)
  {
    std::cout << "Warning: default LAPACK::GEES() not yet implemented" << std::endl;
  }

  template<typename OrdinalType, typename ScalarType>
  ScalarType LAPACK<OrdinalType, ScalarType>::LAPY2(ScalarType x, ScalarType y) const
  {
    std::cout << "Warning: default LAPACK::LAPY2() not yet implemented" << std::endl;
  }

  // END GENERAL TEMPLATE IMPLEMENTATION //

  // BEGIN FLOAT PARTIAL SPECIALIZATION DECLARATION //
  
  template<typename OrdinalType>
  class LAPACK<OrdinalType, float>
  {    
  public:
    inline LAPACK(void) {};
    inline LAPACK(const LAPACK& LAPACK) {};
    inline virtual ~LAPACK(void) {};
    void POTRF(char, int, float*, int, int*) const;
    void POTRS(char, int, int, float*, int, float*, int, int*) const;
    void POTRI(char, int, float*, int, int*) const;
    void POCON(char, int, float*, int, float, float*, float*, int*, int*) const;
    void POSV(char, int, int, float*, int, float*, int, int*) const;
    void POEQU(int, float*, int, float*, float*, float*, int*) const;
    void PORFS(char, int, int, float*, int, float*, int, float*, int, float*, int, float*, float*, float*, int*, int*) const;
    void POSVX(char, char, int, int, float*, int, float*, int, char, float*, float*, int, float*, int, float*, float*, float*, float*, int*, int*) const;
    void GELS(char, int, int, int, float*, int, float*, int, float*, int, int) const;
    void GETRF(int, int, float*, int, int*, int*) const;
    void GETRS(char, int, int, float*, int, int*, float*, int, int*) const;
    void GETRI(int, float*, int, int*, float*, int*, int*) const;
    void GECON(char, int, float*, int, float, float*, float*, int*, int*) const;
    void GESV(int, int, float*, int, int*, float*, int, int*) const;
    void GEEQU(int, int, float*, int, float*, float*, float*, float*, float*, int*) const;
    void GERFS(char, int, int, float*, int, float*, int, int*, float*, int, float*, int, float*, float*, float*, int*, int*) const;
    void GESVX(char, char, int, int, float*, int, float*, int, int*, char, float*, float*, float*, int, float*, int, float*, float*, float*, float*, int*, int*) const;
    void GEHRD(int, int, int, float*, int, float*, float*, int, int*) const;
    void HSEQR(char, char, int, int, int, float*, int, float*, float*, float*, int, float*, int, int*) const;
    void ORGHR(int, int, int, float*, int, float*, float*, int, int*) const;
    void ORMHR(char, char, int, int, int, int, float*, int, float*, float*, int, float*, int, int*) const;
    void TREVC(char, char, int*, int, float*, int, float*, int, float*, int, int, int*, float*, int*) const;
    void TREXC(char, int, float*, int, float*, int, int, int, float*, int*) const;
    float LAMCH(char) const;
    void GEES(char, char, int*, int, float*, int, int, float*, float*, float*, int, float*, int, int*, int*) const;
    float LAPY2(float, float) const;
  };

  // END FLOAT PARTIAL SPECIALIZATION DECLARATION //

  // BEGIN FLOAT PARTIAL SPECIALIZATION IMPLEMENTATION //

  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POTRF(char UPLO, int N, float* A, int LDA, int* INFO) const
  {
    SPOTRF_F77(CHAR_MACRO(UPLO), &N, A, &LDA, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POTRS(char UPLO, int N, int NRHS, float* A, int LDA, float* X, int LDX, int* INFO) const
  {
    SPOTRS_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POTRI(char UPLO, int N, float* A, int LDA, int* INFO) const
  {
    SPOTRI_F77(CHAR_MACRO(UPLO), &N, A, &LDA, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POCON(char UPLO, int N, float* A, int LDA, float ANORM, float* RCOND, float* WORK, int* IWORK, int* INFO) const
  {
    SPOCON_F77(CHAR_MACRO(UPLO), &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POSV(char UPLO, int N, int NRHS, float* A, int LDA, float* X, int LDX, int* INFO) const
  {
    SPOSV_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POEQU(int N, float* A, int LDA, float* S, float* SCOND, float* AMAX, int* INFO) const
  {
    SPOEQU_F77(&N, A, &LDA, S, SCOND, AMAX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::PORFS(char UPLO, int N, int NRHS, float* A, int LDA, float* AF, int LDAF, float* B, int LDB, float* X, int LDX, float* FERR, float* BERR, float* WORK, int* IWORK, int* INFO) const
  {
    SPORFS_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, AF, &LDAF, B, &LDB, X, &LDX, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POSVX(char FACT, char UPLO, int N, int NRHS, float* A, int LDA, float* AF, int LDAF, char EQUED, float* S, float* B, int LDB, float* X, int LDX, float* RCOND, float* FERR, float* BERR, float* WORK, int* IWORK, int* INFO) const
  {
    SPOSVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, AF, &LDAF, EQUED_CHAR_MACRO, S, B, &LDB, X, &LDX, RCOND, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GELS(char TRANS, int m, int n, int numrhs, float* a, int lda, float* b, int ldb, float* work, int lwork, int info) const
  {
    SGELS_F77(CHAR_MACRO(TRANS), &m, &n, &numrhs, a, &lda, b, &ldb, work, &lwork, &info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GETRF(int M, int N, float* A, int LDA, int* IPIV, int* INFO) const
  {
    SGETRF_F77(&M, &N, A, &LDA, IPIV, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GETRS(char TRANS, int N, int NRHS, float* A, int LDA, int* IPIV, float* X, int LDX, int* INFO) const
  {
    SGETRS_F77(CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GETRI(int N, float* A, int LDA, int* IPIV, float* WORK, int* LWORK, int* INFO) const
  {
    SGETRI_F77(&N, A, &LDA, IPIV, WORK, LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GECON(char NORM, int N, float* A, int LDA, float ANORM, float* RCOND, float* WORK, int* IWORK, int* INFO) const
  {
    SGECON_F77(CHAR_MACRO(NORM), &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GESV(int N, int NRHS, float* A, int LDA, int* IPIV, float* X, int LDX, int* INFO) const
  {
    SGESV_F77(&N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GEEQU(int M, int N, float* A, int LDA, float* R, float* C, float* ROWCND, float* COLCND, float* AMAX, int* INFO) const
  {
    SGEEQU_F77(&M, &N, A, &LDA, R, C, ROWCND, COLCND, AMAX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GERFS(char TRANS, int N, int NRHS, float* A, int LDA, float* AF, int LDAF, int* IPIV, float* B, int LDB, float* X, int LDX, float* FERR, float* BERR, float* WORK, int* IWORK, int* INFO) const
  {
    SGERFS_F77(CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, B, &LDB, X, &LDX, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GESVX(char FACT, char TRANS, int N, int NRHS, float* A, int LDA, float* AF, int LDAF, int* IPIV, char EQUED, float* R, float* C, float* B, int LDB, float* X, int LDX, float* RCOND, float* FERR, float* BERR, float* WORK, int* IWORK, int* INFO) const
  {
    SGESVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, CHAR_MACRO(EQUED), R, C, B, &LDB, X, &LDX, RCOND, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GEHRD(int N, int ILO, int IHI, float* A, int LDA, float* TAU, float* WORK, int LWORK, int* INFO) const
  {
    SGEHRD_F77(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::HSEQR(char JOB, char COMPZ, int N, int ILO, int IHI, float* H, int LDH, float* WR, float* WI, float* Z, int LDZ, float* WORK, int LWORK, int* INFO) const
  {
    SHSEQR_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPZ), &N, &ILO, &IHI, H, &LDH, WR, WI, Z, &LDZ, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::ORGHR(int N, int ILO, int IHI, float* A, int LDA, float* TAU, float* WORK, int LWORK, int* INFO) const
  {
    SORGHR_F77(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::ORMHR(char SIDE, char TRANS, int M, int N, int ILO, int IHI, float* A, int LDA, float* TAU, float* C, int LDC, float* WORK, int LWORK, int* INFO) const
  {
    SORMHR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &M, &N, &ILO, &IHI, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::TREVC(char SIDE, char HOWMNY, int* SELECT, int N, float* T, int LDT, float *VL, int LDVL, float* VR, int LDVR, int MM, int* M, float* WORK, int* INFO) const
  {

    if(HOWMNY=='S')
      {
	*INFO = -3; // We do not support 'S' since it requires a logical array (yuck!)
      }
    else
      {
	STREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), SELECT, &N, T, &LDT, VL, &LDVL, VR, &LDVR, &MM, M, WORK, INFO);
      }
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::TREXC(char COMPQ, int N, float* T, int LDT, float* Q, int LDQ, int IFST, int ILST, float* WORK, int* INFO) const
  {
    STREXC_F77(CHAR_MACRO(COMPQ), &N, T, &LDT, Q, &LDQ, &IFST, &ILST, WORK, INFO);
  }

  template<typename OrdinalType>
  float LAPACK<OrdinalType, float>::LAMCH(char CMACH) const
  {
    return(SLAMCH_F77(CHAR_MACRO(CMACH)));
  }

  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GEES(char JOBVS, char SORT, int* SELECT, int N, float* A, int LDA, int SDIM, float* WR, float* WI, float* VS, int LDVS, float* WORK, int LWORK, int *BWORK, int* INFO ) const
  {
    SGEES_F77(CHAR_CHAR_MACRO(JOBVS), CHAR_CHAR_MACRO(SORT), SELECT, &N, A, &LDA, &SDIM, WR, WI, VS, &LDVS, WORK, &LWORK, BWORK, INFO);
  }

  template<typename OrdinalType>
  float LAPACK<OrdinalType, float>::LAPY2(float x, float y) const
  {
    return SLAPY2_F77(&x, &y);
  }

  // END FLOAT PARTIAL SPECIALIZATION IMPLEMENTATION //

  // BEGIN DOUBLE PARTIAL SPECIALIZATION DECLARATION //

  template<typename OrdinalType>
  class LAPACK<OrdinalType, double>
  {    
  public:
    inline LAPACK(void) {};
    inline LAPACK(const LAPACK& LAPACK) {};
    inline virtual ~LAPACK(void) {};
    void POTRF(char, int, double*, int, int*) const;
    void POTRS(char, int, int, double*, int, double*, int, int*) const;
    void POTRI(char, int, double*, int, int*) const;
    void POCON(char, int, double*, int, double, double*, double*, int*, int*) const;
    void POSV(char, int, int, double*, int, double*, int, int*) const;
    void POEQU(int, double*, int, double*, double*, double*, int*) const;
    void PORFS(char, int, int, double*, int, double*, int, double*, int, double*, int, double*, double*, double*, int*, int*) const;
    void POSVX(char, char, int, int, double*, int, double*, int, char, double*, double*, int, double*, int, double*, double*, double*, double*, int*, int*) const;
    void GELS(char, int, int, int, double*, int, double*, int, double*, int, int) const;
    void GETRF(int, int, double*, int, int*, int*) const;
    void GETRS(char, int, int, double*, int, int*, double*, int, int*) const;
    void GETRI(int, double*, int, int*, double*, int*, int*) const;
    void GECON(char, int, double*, int, double, double*, double*, int*, int*) const;
    void GESV(int, int, double*, int, int*, double*, int, int*) const;
    void GEEQU(int, int, double*, int, double*, double*, double*, double*, double*, int*) const;
    void GERFS(char, int, int, double*, int, double*, int, int*, double*, int, double*, int, double*, double*, double*, int*, int*) const;
   void GESVX(char, char, int, int, double*, int, double*, int, int*, char, double*, double*, double*, int, double*, int, double*, double*, double*, double*, int*, int*) const;
    void GEHRD(int, int, int, double*, int, double*, double*, int, int*) const;
    void HSEQR(char, char, int, int, int, double*, int, double*, double*, double*, int, double*, int, int*) const;
    void ORGHR(int, int, int, double*, int, double*, double*, int, int*) const;
    void ORMHR(char, char, int, int, int, int, double*, int, double*, double*, int, double*, int, int*) const;
    void TREVC(char, char, int*, int, double*, int, double*, int, double*, int, int, int*, double*, int*) const;
    void TREXC(char, int, double*, int, double*, int, int, int, double*, int*) const;
    double LAMCH(char) const;
    void GEES(char, char, int*, int, double*, int, int, double*, double*, double*, int, double*, int, int*, int*) const;
    double LAPY2(double, double) const;
  };

  // END DOUBLE PARTIAL SPECIALIZATION DECLARATION //

  // BEGIN DOUBLE PARTIAL SPECIALIZATION IMPLEMENTATION //

  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POTRF(char UPLO, int N, double* A, int LDA, int* INFO) const
  {
    DPOTRF_F77(CHAR_MACRO(UPLO), &N, A, &LDA, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POTRS(char UPLO, int N, int NRHS, double* A, int LDA, double* X, int LDX, int* INFO) const
  {
    DPOTRS_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POTRI(char UPLO, int N, double* A, int LDA, int* INFO) const
  {
    DPOTRI_F77(CHAR_MACRO(UPLO), &N, A, &LDA, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POCON(char UPLO, int N, double* A, int LDA, double ANORM, double* RCOND, double* WORK, int* IWORK, int* INFO) const
  {
    DPOCON_F77(CHAR_MACRO(UPLO), &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POSV(char UPLO, int N, int NRHS, double* A, int LDA, double* X, int LDX, int* INFO) const
  {
    DPOSV_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POEQU(int N, double* A, int LDA, double* S, double* SCOND, double* AMAX, int* INFO) const
  {
    DPOEQU_F77(&N, A, &LDA, S, SCOND, AMAX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::PORFS(char UPLO, int N, int NRHS, double* A, int LDA, double* AF, int LDAF, double* B, int LDB, double* X, int LDX, double* FERR, double* BERR, double* WORK, int* IWORK, int* INFO) const
  {
    DPORFS_F77(CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, AF, &LDAF, B, &LDB, X, &LDX, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POSVX(char FACT, char UPLO, int N, int NRHS, double* A, int LDA, double* AF, int LDAF, char EQUED, double* S, double* B, int LDB, double* X, int LDX, double* RCOND, double* FERR, double* BERR, double* WORK, int* IWORK, int* INFO) const
  {
    DPOSVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(UPLO), &N, &NRHS, A, &LDA, AF, &LDAF, CHAR_MACRO(EQUED), S, B, &LDB, X, &LDX, RCOND, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GELS(char TRANS, int m, int n, int numrhs, double* a, int lda, double* b, int ldb, double* work, int lwork, int info) const
  {
    DGELS_F77(CHAR_MACRO(TRANS), &m, &n, &numrhs, a, &lda, b, &ldb, work, &lwork, &info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GETRF(int M, int N, double* A, int LDA, int* IPIV, int* INFO) const
  {
    DGETRF_F77(&M, &N, A, &LDA, IPIV, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GETRS(char TRANS, int N, int NRHS, double* A, int LDA, int* IPIV, double* X, int LDX, int* INFO) const
  {
    DGETRS_F77(CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GETRI(int N, double* A, int LDA, int* IPIV, double* WORK, int* LWORK, int* INFO) const
  {
    DGETRI_F77(&N, A, &LDA, IPIV, WORK, LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GECON(char NORM, int N, double* A, int LDA, double ANORM, double* RCOND, double* WORK, int* IWORK, int* INFO) const
  {
    DGECON_F77(CHAR_MACRO(NORM), &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GESV(int N, int NRHS, double* A, int LDA, int* IPIV, double* X, int LDX, int* INFO) const
  {
    DGESV_F77(&N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GEEQU(int M, int N, double* A, int LDA, double* R, double* C, double* ROWCND, double* COLCND, double* AMAX, int* INFO) const
  {
    DGEEQU_F77(&M, &N, A, &LDA, R, C, ROWCND, COLCND, AMAX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GERFS(char TRANS, int N, int NRHS, double* A, int LDA, double* AF, int LDAF, int* IPIV, double* B, int LDB, double* X, int LDX, double* FERR, double* BERR, double* WORK, int* IWORK, int* INFO) const
  {
    DGERFS_F77(CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, B, &LDB, X, &LDX, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GESVX(char FACT, char TRANS, int N, int NRHS, double* A, int LDA, double* AF, int LDAF, int* IPIV, char EQUED, double* R, double* C, double* B, int LDB, double* X, int LDX, double* RCOND, double* FERR, double* BERR, double* WORK, int* IWORK, int* INFO) const
  {
    DGESVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(TRANS), &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, CHAR_MACRO(EQUED), R, C, B, &LDB, X, &LDX, RCOND, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GEHRD(int N, int ILO, int IHI, double* A, int LDA, double* TAU, double* WORK, int LWORK, int* INFO) const
  {
    DGEHRD_F77(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::HSEQR(char JOB, char COMPZ, int N, int ILO, int IHI, double* H, int LDH, double* WR, double* WI, double* Z, int LDZ, double* WORK, int LWORK, int* INFO) const
  {
    DHSEQR_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPZ), &N, &ILO, &IHI, H, &LDH, WR, WI, Z, &LDZ, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::ORGHR(int N, int ILO, int IHI, double* A, int LDA, double* TAU, double* WORK, int LWORK, int* INFO) const
  {
    DORGHR_F77(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::ORMHR(char SIDE, char TRANS, int M, int N, int ILO, int IHI, double* A, int LDA, double* TAU, double* C, int LDC, double* WORK, int LWORK, int* INFO) const
  {
    DORMHR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &M, &N, &ILO, &IHI, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::TREVC(char SIDE, char HOWMNY, int* SELECT, int N, double* T, int LDT, double *VL, int LDVL, double* VR, int LDVR, int MM, int* M, double* WORK, int* INFO) const
  {
    if(HOWMNY=='S')
      {
	*INFO = -3; // We do not support 'S' since it requires a logical array (yuck!)
      }
    else
      {
	DTREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), SELECT, &N, T, &LDT, VL, &LDVL, VR, &LDVR, &MM, M, WORK, INFO);
      }
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::TREXC(char COMPQ, int N, double* T, int LDT, double* Q, int LDQ, int IFST, int ILST, double* WORK, int* INFO) const
  {
    DTREXC_F77(CHAR_MACRO(COMPQ), &N, T, &LDT, Q, &LDQ, &IFST, &ILST, WORK, INFO);
  }

  template<typename OrdinalType>
  double LAPACK<OrdinalType, double>::LAMCH(char CMACH) const
  {
    return(DLAMCH_F77(CHAR_MACRO(CMACH)));
  }

  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GEES(char JOBVS, char SORT, int* SELECT, int N, double* A, int LDA, int SDIM, double* WR, double* WI, double* VS, int LDVS, double* WORK, int LWORK, int* BWORK, int* INFO) const
  {
    DGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(SORT), SELECT, &N, A, &LDA, &SDIM, WR, WI, VS, &LDVS, WORK, &LWORK, BWORK, INFO);
  }

  template<typename OrdinalType>
  double LAPACK<OrdinalType, double>::LAPY2(double x, double y) const
  {
    return DLAPY2_F77(&x, &y);
  }

  // END DOUBLE PARTIAL SPECIALIZATION IMPLEMENTATION //

} // namespace Teuchos

#endif // _TEUCHOS_LAPACK_HPP_
