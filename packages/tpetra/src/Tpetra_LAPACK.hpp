/*Paul
27-May-2002 General cleanup. Changed method names to fit namingConvention (nothing changed).
06-August-2002 Changed to images (nothing changed).
*/

// Kris
// 06.11.03 -- Format cleanup
// 06.13.03 -- Initial templatization (Tpetra_LAPACK.cpp is no longer needed)
// 06.17.03 -- Changed LAPACK wrapper calls from XXXXX_F77() to xxxxx_()
//          -- Added warning messages to default function calls
//          -- Added LAPY2 and GEES by request
//
// To Do:
//     Figure out how we intend to implement Complex numbers
//     Rename GEES parameters
//     Make LAPY2 and GEES const functions

#ifndef _TPETRA_LAPACK_HPP_
#define _TPETRA_LAPACK_HPP_

#if defined (INTEL_CXML)
#define UPLO_MACRO &UPLO, one
#define TRANS_MACRO &TRANS, one
#define FACT_MACRO &FACT, one
#define EQUED_MACRO &EQUED, one
#define NORM_MACRO &NORM, one
#define JOB_MACRO &JOB, one
#define COMPZ_MACRO &COMPZ, one
#define SIDE_MACRO &SIDE, one
#define HOWMNY_MACRO &HOWMNY, one
#define COMPQ_MACRO &COMPQ, one
#define CMACH_MACRO &CMACH, one
#else
#define UPLO_MACRO &UPLO
#define TRANS_MACRO &TRANS
#define FACT_MACRO &FACT
#define EQUED_MACRO &EQUED
#define NORM_MACRO &NORM
#define JOB_MACRO &JOB
#define COMPZ_MACRO &COMPZ
#define SIDE_MACRO &SIDE
#define HOWMNY_MACRO &HOWMNY
#define COMPQ_MACRO &COMPQ
#define CMACH_MACRO &CMACH
#endif

#include "Tpetra_Object.hpp"
#include "Tpetra_LAPACK_wrappers.hpp"

namespace Tpetra
{

  // BEGIN GENERAL TEMPLATE DECLARATION //

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
    void GEES(char, char, int*, int, ScalarType*, int, int, ScalarType*, ScalarType*, ScalarType*, ScalarType*, int, ScalarType*, int, int, int, int);
    ScalarType LAPY2(ScalarType, ScalarType);
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
  void LAPACK<OrdinalType, ScalarType>::GEES(char A, char B, int* C, int D, ScalarType* E, int F, int G, ScalarType* H, ScalarType* I, ScalarType* J, ScalarType* K, int L, ScalarType* M, int N, int O, int P, int Q)
  {
    std::cout << "Warning: default LAPACK::GEES() not yet implemented" << std::endl;
  }

  template<typename OrdinalType, typename ScalarType>
  ScalarType LAPACK<OrdinalType, ScalarType>::LAPY2(ScalarType x, ScalarType y)
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
    void GEES(char, char, int*, int, float*, int, int, float*, float*, float*, float*, int, float*, int, int, int, int);
    float LAPY2(float, float);
  };

  // END FLOAT PARTIAL SPECIALIZATION DECLARATION //

  // BEGIN FLOAT PARTIAL SPECIALIZATION IMPLEMENTATION //

  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POTRF(char UPLO, int N, float* A, int LDA, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    spotrf_(UPLO_MACRO, &N, A, &LDA, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POTRS(char UPLO, int N, int NRHS, float* A, int LDA, float* X, int LDX, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    spotrs_(UPLO_MACRO, &N, &NRHS, A, &LDA, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POTRI(char UPLO, int N, float* A, int LDA, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    spotri_(UPLO_MACRO, &N, A, &LDA, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POCON(char UPLO, int N, float* A, int LDA, float ANORM, float* RCOND, float* WORK, int* IWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    spocon_(UPLO_MACRO, &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POSV(char UPLO, int N, int NRHS, float* A, int LDA, float* X, int LDX, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    sposv_(UPLO_MACRO, &N, &NRHS, A, &LDA, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POEQU(int N, float* A, int LDA, float* S, float* SCOND, float* AMAX, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    spoequ_(&N, A, &LDA, S, SCOND, AMAX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::PORFS(char UPLO, int N, int NRHS, float* A, int LDA, float* AF, int LDAF, float* B, int LDB, float* X, int LDX, float* FERR, float* BERR, float* WORK, int* IWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    sporfs_(UPLO_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, B, &LDB, X, &LDX, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POSVX(char FACT, char UPLO, int N, int NRHS, float* A, int LDA, float* AF, int LDAF, char EQUED, float* S, float* B, int LDB, float* X, int LDX, float* RCOND, float* FERR, float* BERR, float* WORK, int* IWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    sposvx_(FACT_MACRO, UPLO_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, EQUED_MACRO, S, B, &LDB, X, &LDX, RCOND, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GELS(char TRANS, int m, int n, int numrhs, float* a, int lda, float* b, int ldb, float* work, int lwork, int info) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    sgels_(TRANS_MACRO, &m, &n, &numrhs, a, &lda, b, &ldb, work, &lwork, &info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GETRF(int M, int N, float* A, int LDA, int* IPIV, int* INFO) const
  {
    sgetrf_(&M, &N, A, &LDA, IPIV, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GETRS(char TRANS, int N, int NRHS, float* A, int LDA, int* IPIV, float* X, int LDX, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    sgetrs_(TRANS_MACRO, &N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GETRI(int N, float* A, int LDA, int* IPIV, float* WORK, int* LWORK, int* INFO) const
  {
    sgetri_(&N, A, &LDA, IPIV, WORK, LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GECON(char NORM, int N, float* A, int LDA, float ANORM, float* RCOND, float* WORK, int* IWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    sgecon_(NORM_MACRO, &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GESV(int N, int NRHS, float* A, int LDA, int* IPIV, float* X, int LDX, int* INFO) const
  {
    sgesv_(&N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GEEQU(int M, int N, float* A, int LDA, float* R, float* C, float* ROWCND, float* COLCND, float* AMAX, int* INFO) const
  {
    sgeequ_(&M, &N, A, &LDA, R, C, ROWCND, COLCND, AMAX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GERFS(char TRANS, int N, int NRHS, float* A, int LDA, float* AF, int LDAF, int* IPIV, float* B, int LDB, float* X, int LDX, float* FERR, float* BERR, float* WORK, int* IWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    sgerfs_(TRANS_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, B, &LDB, X, &LDX, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GESVX(char FACT, char TRANS, int N, int NRHS, float* A, int LDA, float* AF, int LDAF, int* IPIV, char EQUED, float* R, float* C, float* B, int LDB, float* X, int LDX, float* RCOND, float* FERR, float* BERR, float* WORK, int* IWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    sgesvx_(FACT_MACRO, TRANS_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, EQUED_MACRO, R, C, B, &LDB, X, &LDX, RCOND, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GEHRD(int N, int ILO, int IHI, float* A, int LDA, float* TAU, float* WORK, int LWORK, int* INFO) const
  {
    sgehrd_(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::HSEQR(char JOB, char COMPZ, int N, int ILO, int IHI, float* H, int LDH, float* WR, float* WI, float* Z, int LDZ, float* WORK, int LWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    shseqr_(JOB_MACRO, COMPZ_MACRO, &N, &ILO, &IHI, H, &LDH, WR, WI, Z, &LDZ, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::ORGHR(int N, int ILO, int IHI, float* A, int LDA, float* TAU, float* WORK, int LWORK, int* INFO) const
  {
    sorghr_(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::ORMHR(char SIDE, char TRANS, int M, int N, int ILO, int IHI, float* A, int LDA, float* TAU, float* C, int LDC, float* WORK, int LWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    sormhr_(SIDE_MACRO, TRANS_MACRO, &M, &N, &ILO, &IHI, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::TREVC(char SIDE, char HOWMNY, int* SELECT, int N, float* T, int LDT, float *VL, int LDVL, float* VR, int LDVR, int MM, int* M, float* WORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    if(HOWMNY=='S')
      {
	*INFO = -3; // We do not support 'S' since it requires a logical array (yuck!)
      }
    else
      {
	strevc_(SIDE_MACRO, HOWMNY_MACRO, SELECT, &N, T, &LDT, VL, &LDVL, VR, &LDVR, &MM, M, WORK, INFO);
      }
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::TREXC(char COMPQ, int N, float* T, int LDT, float* Q, int LDQ, int IFST, int ILST, float* WORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    strexc_(COMPQ_MACRO, &N, T, &LDT, Q, &LDQ, &IFST, &ILST, WORK, INFO);
  }

  template<typename OrdinalType>
  float LAPACK<OrdinalType, float>::LAMCH(char CMACH) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    return(slamch_(CMACH_MACRO));
  }

// void AnasaziLAPACK::GEES(char JOBVS, char SORT, int* SELECT, int N, double * A, int LDA,
//                        int SDIM, double* WR, double * WI, double* VS, int LDVS, double* WORK,
//                        int LWORK, int *BWORK, int* INFO ) const {
// DGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(SORT), SELECT, &N, A, &LDA, &SDIM, WR, WI, VS, &LDVS, WORK, &LWORK, BWORK, INFO)$}

  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GEES(char A, char B, int* C, int D, float* E, int F, int G, float* H, float* I, float* J, float* K, int L, float* M, int N, int O, int P, int Q)
  {
    sgees_(&A, &B, C, &D, E, &F, &G, H, I, J, K, &L, M, &N, &O, &P, &Q);
  }

  template<typename OrdinalType>
  float LAPACK<OrdinalType, float>::LAPY2(float x, float y)
  {
    return slapy2_(&x, &y);
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
    void GEES(char, char, int*, int, double*, int, int, double*, double*, double*, double*, int, double*, int, int, int, int);
    double LAPY2(double, double);
  };

  // END DOUBLE PARTIAL SPECIALIZATION DECLARATION //

  // BEGIN DOUBLE PARTIAL SPECIALIZATION IMPLEMENTATION //

  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POTRF(char UPLO, int N, double* A, int LDA, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dpotrf_(UPLO_MACRO, &N, A, &LDA, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POTRS(char UPLO, int N, int NRHS, double* A, int LDA, double* X, int LDX, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dpotrs_(UPLO_MACRO, &N, &NRHS, A, &LDA, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POTRI(char UPLO, int N, double* A, int LDA, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dpotri_(UPLO_MACRO, &N, A, &LDA, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POCON(char UPLO, int N, double* A, int LDA, double ANORM, double* RCOND, double* WORK, int* IWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dpocon_(UPLO_MACRO, &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POSV(char UPLO, int N, int NRHS, double* A, int LDA, double* X, int LDX, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dposv_(UPLO_MACRO, &N, &NRHS, A, &LDA, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POEQU(int N, double* A, int LDA, double* S, double* SCOND, double* AMAX, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dpoequ_(&N, A, &LDA, S, SCOND, AMAX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::PORFS(char UPLO, int N, int NRHS, double* A, int LDA, double* AF, int LDAF, double* B, int LDB, double* X, int LDX, double* FERR, double* BERR, double* WORK, int* IWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dporfs_(UPLO_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, B, &LDB, X, &LDX, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POSVX(char FACT, char UPLO, int N, int NRHS, double* A, int LDA, double* AF, int LDAF, char EQUED, double* S, double* B, int LDB, double* X, int LDX, double* RCOND, double* FERR, double* BERR, double* WORK, int* IWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dposvx_(FACT_MACRO, UPLO_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, EQUED_MACRO, S, B, &LDB, X, &LDX, RCOND, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GELS(char TRANS, int m, int n, int numrhs, double* a, int lda, double* b, int ldb, double* work, int lwork, int info) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dgels_(TRANS_MACRO, &m, &n, &numrhs, a, &lda, b, &ldb, work, &lwork, &info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GETRF(int M, int N, double* A, int LDA, int* IPIV, int* INFO) const
  {
    dgetrf_(&M, &N, A, &LDA, IPIV, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GETRS(char TRANS, int N, int NRHS, double* A, int LDA, int* IPIV, double* X, int LDX, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dgetrs_(TRANS_MACRO, &N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GETRI(int N, double* A, int LDA, int* IPIV, double* WORK, int* LWORK, int* INFO) const
  {
    dgetri_(&N, A, &LDA, IPIV, WORK, LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GECON(char NORM, int N, double* A, int LDA, double ANORM, double* RCOND, double* WORK, int* IWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dgecon_(NORM_MACRO, &N, A, &LDA, &ANORM, RCOND, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GESV(int N, int NRHS, double* A, int LDA, int* IPIV, double* X, int LDX, int* INFO) const
  {
    dgesv_(&N, &NRHS, A, &LDA, IPIV, X, &LDX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GEEQU(int M, int N, double* A, int LDA, double* R, double* C, double* ROWCND, double* COLCND, double* AMAX, int* INFO) const
  {
    dgeequ_(&M, &N, A, &LDA, R, C, ROWCND, COLCND, AMAX, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GERFS(char TRANS, int N, int NRHS, double* A, int LDA, double* AF, int LDAF, int* IPIV, double* B, int LDB, double* X, int LDX, double* FERR, double* BERR, double* WORK, int* IWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dgerfs_(TRANS_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, B, &LDB, X, &LDX, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GESVX(char FACT, char TRANS, int N, int NRHS, double* A, int LDA, double* AF, int LDAF, int* IPIV, char EQUED, double* R, double* C, double* B, int LDB, double* X, int LDX, double* RCOND, double* FERR, double* BERR, double* WORK, int* IWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dgesvx_(FACT_MACRO, TRANS_MACRO, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, EQUED_MACRO, R, C, B, &LDB, X, &LDX, RCOND, FERR, BERR, WORK, IWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GEHRD(int N, int ILO, int IHI, double* A, int LDA, double* TAU, double* WORK, int LWORK, int* INFO) const
  {
    dgehrd_(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::HSEQR(char JOB, char COMPZ, int N, int ILO, int IHI, double* H, int LDH, double* WR, double* WI, double* Z, int LDZ, double* WORK, int LWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dhseqr_(JOB_MACRO, COMPZ_MACRO, &N, &ILO, &IHI, H, &LDH, WR, WI, Z, &LDZ, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::ORGHR(int N, int ILO, int IHI, double* A, int LDA, double* TAU, double* WORK, int LWORK, int* INFO) const
  {
    dorghr_(&N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::ORMHR(char SIDE, char TRANS, int M, int N, int ILO, int IHI, double* A, int LDA, double* TAU, double* C, int LDC, double* WORK, int LWORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dormhr_(SIDE_MACRO, TRANS_MACRO, &M, &N, &ILO, &IHI, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::TREVC(char SIDE, char HOWMNY, int* SELECT, int N, double* T, int LDT, double *VL, int LDVL, double* VR, int LDVR, int MM, int* M, double* WORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    if(HOWMNY=='S')
      {
	*INFO = -3; // We do not support 'S' since it requires a logical array (yuck!)
      }
    else
      {
	dtrevc_(SIDE_MACRO, HOWMNY_MACRO, SELECT, &N, T, &LDT, VL, &LDVL, VR, &LDVR, &MM, M, WORK, INFO);
      }
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::TREXC(char COMPQ, int N, double* T, int LDT, double* Q, int LDQ, int IFST, int ILST, double* WORK, int* INFO) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    dtrexc_(COMPQ_MACRO, &N, T, &LDT, Q, &LDQ, &IFST, &ILST, WORK, INFO);
  }

  template<typename OrdinalType>
  double LAPACK<OrdinalType, double>::LAMCH(char CMACH) const
  {
  #if defined (INTEL_CXML)
    unsigned int one = 1;
  #endif
    return(dlamch_(CMACH_MACRO));
  }

  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GEES(char A, char B, int* C, int D, double* E, int F, int G, double* H, double* I, double* J, double* K, int L, double* M, int N, int O, int P, int Q)
  {
    dgees_(&A, &B, C, &D, E, &F, &G, H, I, J, K, &L, M, &N, &O, &P, &Q);
  }

  template<typename OrdinalType>
  double LAPACK<OrdinalType, double>::LAPY2(double x, double y)
  {
    dlapy2_(&x, &y);
  }

  // END DOUBLE PARTIAL SPECIALIZATION IMPLEMENTATION //

} // namespace Tpetra
#endif /* _TPETRA_LAPACK_HPP_ */
