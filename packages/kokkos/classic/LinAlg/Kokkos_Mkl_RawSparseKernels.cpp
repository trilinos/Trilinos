#include "Kokkos_Mkl_RawSparseKernels.hpp"
#include "Teuchos_ConfigDefs.hpp"

#ifdef HAVE_KOKKOSCLASSIC_MKL
#  include <mkl.h>
#else
#  error "Compiling this file requires building Trilinos with the Intel MKL third-party library."
#  include <Teuchos_StandardCatchMacros.hpp>
#  ifdef MKL_INT
#    undef MKL_INT
#  endif // MKL_INT
#  define MKL_INT int

namespace {
  const char missingMklMsg[] =
    "Calling Intel MKL functions from Trilinos requires having built Trilinos "
    "with the Intel MKL third-party library.  You will need to reconfigure "
    "Trilinos to build with MKL, by setting the following CMake options:"
    "\n"
    "TPL_ENABLE_MKL:BOOL=ON\n"
    "MKL_LIBRARY_DIRS:FILEPATH=<path-to-MKL-libraries>\n"
    "MKL_LIBRARY_NAMES:STRING=\"<list-of-MKL-libraries>\"\n"
    "MKL_INCLUDE_DIRS:FILEPATH=<path-to-MKL-headers>\n"
    "\n"
    "Replace <path-to-MKL-libraries> with the path to MKL libraries, "
    "<list-of-MKL-libraries> with the list of MKL libraries with which "
    "to link (separating library names with a semicolon, and not including "
    "the -l or lib part of the name), and <path-to-MKL-headers> with the "
    "path to MKL header files.";
} // namespace (anonymous)

#endif // HAVE_KOKKOSCLASSIC_MKL


//////////////////////////////////////////////////////////////////////
// Specializations for Scalar=float
//////////////////////////////////////////////////////////////////////

namespace Kokkos {
  namespace Mkl {

    template<>
    void RawSparseKernels<float, MKL_INT>::
    csrmv (const char* const transa,
           const MKL_INT m, // Number of rows in A
           const MKL_INT k, // Number of columns in A
           const float& alpha,
           const char* const matdescra,
           const float* const val,
           const MKL_INT* const ind,
           const MKL_INT* const ptrBegin,
           const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
           const float* const x,
           const float& beta,
           float* const y)
    {
#ifdef HAVE_KOKKOSCLASSIC_MKL
      mkl_scsrmv ((char*) transa, (MKL_INT*) &m, (MKL_INT*) &k,
                  (float*) &alpha, (char*) matdescra,
                  (float*) val, (MKL_INT*) ind, (MKL_INT*) ptrBegin, (MKL_INT*) ptrEnd,
                  (float*) x, (float*) &beta, y);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
    }

    template<>
    void RawSparseKernels<float, MKL_INT>::
    csrmm (const char* const transa,
           const MKL_INT m, // number of rows of A
           const MKL_INT n, // number of columns of C
           const MKL_INT k, // number of columns of A
           const float& alpha,
           const char* const matdescra,
           const float* const val,
           const MKL_INT* const ind,
           const MKL_INT* const ptrBegin,
           const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
           const float* const B,
           const MKL_INT LDB,
           const float& beta,
           float* const C,
           const MKL_INT LDC)
    {
#ifdef HAVE_KOKKOSCLASSIC_MKL
      mkl_scsrmm ((char*) transa, (MKL_INT*) &m, (MKL_INT*) &n, (MKL_INT*) &k,
                  (float*) &alpha, (char*) matdescra,
                  (float*) val, (MKL_INT*) ind, (MKL_INT*) ptrBegin, (MKL_INT*) ptrEnd,
                  (float*) B, (MKL_INT*) &LDB, (float*) &beta, C, (MKL_INT*) &LDC);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
    }

    template<>
    void RawSparseKernels<float, MKL_INT>::
    csrsv (const char* const transa,
           const MKL_INT m,
           const float& alpha,
           const char* const matdescra,
           const float* const val,
           const MKL_INT* const ind,
           const MKL_INT* const ptrBegin,
           const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
           const float* const x,
           float* const y)
    {
#ifdef HAVE_KOKKOSCLASSIC_MKL
      mkl_scsrsv ((char*) transa, (MKL_INT*) &m, (float*) &alpha, (char*) matdescra,
                  (float*) val, (MKL_INT*) ind, (MKL_INT*) ptrBegin, (MKL_INT*) ptrEnd,
                  (float*) x, y);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
    }

    template<>
    void RawSparseKernels<float, MKL_INT>::
    csrsm (const char* const transa,
           const MKL_INT m, // Number of columns in A
           const MKL_INT n, // Number of columns in C
           const float& alpha,
           const char* const matdescra,
           const float* const val,
           const MKL_INT* const ind,
           const MKL_INT* const ptrBegin,
           const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
           const float* const B,
           const MKL_INT LDB,
           float* const C,
           const MKL_INT LDC)
    {
#ifdef HAVE_KOKKOSCLASSIC_MKL
      mkl_scsrsm ((char*) transa, (MKL_INT*) &m, (MKL_INT*) &n,
                  (float*) &alpha, (char*) matdescra,
                  (float*) val, (MKL_INT*) ind, (MKL_INT*) ptrBegin, (MKL_INT*) ptrEnd,
                  (float*) B, (MKL_INT*) &LDB, C, (MKL_INT*) &LDC);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
    }

    //////////////////////////////////////////////////////////////////////
    // Specializations for Scalar=double
    //////////////////////////////////////////////////////////////////////

    template<>
    void RawSparseKernels<double, MKL_INT>::
    csrmv (const char* const transa,
           const MKL_INT m, // Number of rows in A
           const MKL_INT k, // Number of columns in A
           const double& alpha,
           const char* const matdescra,
           const double* const val,
           const MKL_INT* const ind,
           const MKL_INT* const ptrBegin,
           const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
           const double* const x,
           const double& beta,
           double* const y)
    {
#ifdef HAVE_KOKKOSCLASSIC_MKL
      mkl_dcsrmv ((char*) transa, (MKL_INT*) &m, (MKL_INT*) &k,
                  (double*) &alpha, (char*) matdescra,
                  (double*) val, (MKL_INT*) ind, (MKL_INT*) ptrBegin, (MKL_INT*) ptrEnd,
                  (double*) x, (double*) &beta, y);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
    }

    template<>
    void RawSparseKernels<double, MKL_INT>::
    csrmm (const char* const transa,
           const MKL_INT m, // number of rows of A
           const MKL_INT n, // number of columns of C
           const MKL_INT k, // number of columns of A
           const double& alpha,
           const char* const matdescra,
           const double* const val,
           const MKL_INT* const ind,
           const MKL_INT* const ptrBegin,
           const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
           const double* const B,
           const MKL_INT LDB,
           const double& beta,
           double* const C,
           const MKL_INT LDC)
    {
#ifdef HAVE_KOKKOSCLASSIC_MKL
      mkl_dcsrmm ((char*) transa, (MKL_INT*) &m, (MKL_INT*) &n, (MKL_INT*) &k,
                  (double*) &alpha, (char*) matdescra,
                  (double*) val, (MKL_INT*) ind, (MKL_INT*) ptrBegin, (MKL_INT*) ptrEnd,
                  (double*) B, (MKL_INT*) &LDB, (double*) &beta, C, (MKL_INT*) &LDC);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
    }

    template<>
    void RawSparseKernels<double, MKL_INT>::
    csrsv (const char* const transa,
           const MKL_INT m,
           const double& alpha,
           const char* const matdescra,
           const double* const val,
           const MKL_INT* const ind,
           const MKL_INT* const ptrBegin,
           const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
           const double* const x,
           double* const y)
    {
#ifdef HAVE_KOKKOSCLASSIC_MKL
      mkl_dcsrsv ((char*) transa, (MKL_INT*) &m, (double*) &alpha, (char*) matdescra,
                  (double*) val, (MKL_INT*) ind, (MKL_INT*) ptrBegin, (MKL_INT*) ptrEnd,
                  (double*) x, y);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
    }

    template<>
    void RawSparseKernels<double, MKL_INT>::
    csrsm (const char* const transa,
           const MKL_INT m, // Number of columns in A
           const MKL_INT n, // Number of columns in C
           const double& alpha,
           const char* const matdescra,
           const double* const val,
           const MKL_INT* const ind,
           const MKL_INT* const ptrBegin,
           const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
           const double* const B,
           const MKL_INT LDB,
           double* const C,
           const MKL_INT LDC)
    {
#ifdef HAVE_KOKKOSCLASSIC_MKL
      mkl_dcsrsm ((char*) transa, (MKL_INT*) &m, (MKL_INT*) &n,
                  (double*) &alpha, (char*) matdescra,
                  (double*) val, (MKL_INT*) ind, (MKL_INT*) ptrBegin, (MKL_INT*) ptrEnd,
                  (double*) B, (MKL_INT*) &LDB, C, (MKL_INT*) &LDC);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
    }

    // MKL always defines the complex-arithmetic routines, but only build
    // wrappers for them if Teuchos' complex arithmetic support is enabled.
#ifdef HAVE_TEUCHOS_COMPLEX

    //////////////////////////////////////////////////////////////////////
    // Specializations for Scalar=std::complex<float>
    //////////////////////////////////////////////////////////////////////

    template<>
    void RawSparseKernels<std::complex<float>, MKL_INT>::
    csrmv (const char* const transa,
           const MKL_INT m, // Number of rows in A
           const MKL_INT k, // Number of columns in A
           const std::complex<float>& alpha,
           const char* const matdescra,
           const std::complex<float>* const val,
           const MKL_INT* const ind,
           const MKL_INT* const ptrBegin,
           const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
           const std::complex<float>* const x,
           const std::complex<float>& beta,
           std::complex<float>* const y)
    {
#ifdef HAVE_KOKKOSCLASSIC_MKL
      mkl_ccsrmv ((char*) transa,
                  (MKL_INT*) &m,
                  (MKL_INT*) &k,
                  reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (&alpha)),
                  (char*) matdescra,
                  reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (val)),
                  (MKL_INT*) ind,
                  (MKL_INT*) ptrBegin,
                  (MKL_INT*) ptrEnd,
                  reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (x)),
                  reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (&beta)),
                  reinterpret_cast<MKL_Complex8*> (y));
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
    }

    template<>
    void RawSparseKernels<std::complex<float>, MKL_INT>::
    csrmm (const char* const transa,
           const MKL_INT m, // number of rows of A
           const MKL_INT n, // number of columns of C
           const MKL_INT k, // number of columns of A
           const std::complex<float>& alpha,
           const char* const matdescra,
           const std::complex<float>* const val,
           const MKL_INT* const ind,
           const MKL_INT* const ptrBegin,
           const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
           const std::complex<float>* const B,
           const MKL_INT LDB,
           const std::complex<float>& beta,
           std::complex<float>* const C,
           const MKL_INT LDC)
    {
#ifdef HAVE_KOKKOSCLASSIC_MKL
      mkl_ccsrmm ((char*) transa,
                  (MKL_INT*) &m,
                  (MKL_INT*) &n,
                  (MKL_INT*) &k,
                  reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (&alpha)),
                  (char*) matdescra,
                  reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (val)),
                  (MKL_INT*) ind,
                  (MKL_INT*) ptrBegin,
                  (MKL_INT*) ptrEnd,
                  reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (B)), (MKL_INT*) &LDB,
                  reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (&beta)),
                  reinterpret_cast<MKL_Complex8*> (C), (MKL_INT*) &LDC);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
  }

  template<>
  void RawSparseKernels<std::complex<float>, MKL_INT>::
  csrsv (const char* const transa,
         const MKL_INT m,
         const std::complex<float>& alpha,
         const char* const matdescra,
         const std::complex<float>* const val,
         const MKL_INT* const ind,
         const MKL_INT* const ptrBegin,
         const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
         const std::complex<float>* const x,
         std::complex<float>* const y)
  {
#ifdef HAVE_KOKKOSCLASSIC_MKL
    mkl_ccsrsv ((char*) transa,
                (MKL_INT*) &m,
                reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (&alpha)),
                (char*) matdescra,
                reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (val)),
                (MKL_INT*) ind,
                (MKL_INT*) ptrBegin,
                (MKL_INT*) ptrEnd,
                reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (x)),
                reinterpret_cast<MKL_Complex8*> (y));
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
  }

  template<>
  void RawSparseKernels<std::complex<float>, MKL_INT>::
  csrsm (const char* const transa,
         const MKL_INT m, // Number of columns in A
         const MKL_INT n, // Number of columns in C
         const std::complex<float>& alpha,
         const char* const matdescra,
         const std::complex<float>* const val,
         const MKL_INT* const ind,
         const MKL_INT* const ptrBegin,
         const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
         const std::complex<float>* const B,
         const MKL_INT LDB,
         std::complex<float>* const C,
         const MKL_INT LDC)
  {
#ifdef HAVE_KOKKOSCLASSIC_MKL
    mkl_ccsrsm ((char*) transa,
                (MKL_INT*) &m,
                (MKL_INT*) &n,
                reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (&alpha)),
                (char*) matdescra,
                reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (val)),
                (MKL_INT*) ind,
                (MKL_INT*) ptrBegin,
                (MKL_INT*) ptrEnd,
                reinterpret_cast<MKL_Complex8*> (const_cast<std::complex<float>* > (B)), (MKL_INT*) &LDB,
                reinterpret_cast<MKL_Complex8*> (C), (MKL_INT*) &LDC);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
  }

  //////////////////////////////////////////////////////////////////////
  // Specializations for Scalar=std::complex<double>
  //////////////////////////////////////////////////////////////////////

  template<>
  void RawSparseKernels<std::complex<double>, MKL_INT>::
  csrmv (const char* const transa,
         const MKL_INT m, // Number of rows in A
         const MKL_INT k, // Number of columns in A
         const std::complex<double>& alpha,
         const char* const matdescra,
         const std::complex<double>* const val,
         const MKL_INT* const ind,
         const MKL_INT* const ptrBegin,
         const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
         const std::complex<double>* const x,
         const std::complex<double>& beta,
         std::complex<double>* const y)
  {
#ifdef HAVE_KOKKOSCLASSIC_MKL
    mkl_zcsrmv ((char*) transa,
                (MKL_INT*) &m,
                (MKL_INT*) &k,
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (&alpha)),
                (char*) matdescra,
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (val)),
                (MKL_INT*) ind,
                (MKL_INT*) ptrBegin,
                (MKL_INT*) ptrEnd,
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (x)),
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (&beta)),
                reinterpret_cast<MKL_Complex16*> (y));
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
  }

  template<>
  void RawSparseKernels<std::complex<double>, MKL_INT>::
  csrmm (const char* const transa,
         const MKL_INT m, // number of rows of A
         const MKL_INT n, // number of columns of C
         const MKL_INT k, // number of columns of A
         const std::complex<double>& alpha,
         const char* const matdescra,
         const std::complex<double>* const val,
         const MKL_INT* const ind,
         const MKL_INT* const ptrBegin,
         const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
         const std::complex<double>* const B,
         const MKL_INT LDB,
         const std::complex<double>& beta,
         std::complex<double>* const C,
         const MKL_INT LDC)
  {
#ifdef HAVE_KOKKOSCLASSIC_MKL
    mkl_zcsrmm ((char*) transa,
                (MKL_INT*) &m,
                (MKL_INT*) &n,
                (MKL_INT*) &k,
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (&alpha)),
                (char*) matdescra,
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (val)),
                (MKL_INT*) ind,
                (MKL_INT*) ptrBegin,
                (MKL_INT*) ptrEnd,
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (B)), (MKL_INT*) &LDB,
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (&beta)),
                reinterpret_cast<MKL_Complex16*> (C), (MKL_INT*) &LDC);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
  }

  template<>
  void RawSparseKernels<std::complex<double>, MKL_INT>::
  csrsv (const char* const transa,
         const MKL_INT m,
         const std::complex<double>& alpha,
         const char* const matdescra,
         const std::complex<double>* const val,
         const MKL_INT* const ind,
         const MKL_INT* const ptrBegin,
         const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
         const std::complex<double>* const x,
         std::complex<double>* const y)
  {
#ifdef HAVE_KOKKOSCLASSIC_MKL
    mkl_zcsrsv ((char*) transa,
                (MKL_INT*) &m,
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (&alpha)),
                (char*) matdescra,
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (val)),
                (MKL_INT*) ind,
                (MKL_INT*) ptrBegin,
                (MKL_INT*) ptrEnd,
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (x)),
                reinterpret_cast<MKL_Complex16*> (y));
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
  }

  template<>
  void RawSparseKernels<std::complex<double>, MKL_INT>::
  csrsm (const char* const transa,
         const MKL_INT m, // Number of columns in A
         const MKL_INT n, // Number of columns in C
         const std::complex<double>& alpha,
         const char* const matdescra,
         const std::complex<double>* const val,
         const MKL_INT* const ind,
         const MKL_INT* const ptrBegin,
         const MKL_INT* const ptrEnd, // hint: ptrEnd = &ptrBegin[1]
         const std::complex<double>* const B,
         const MKL_INT LDB,
         std::complex<double>* const C,
         const MKL_INT LDC)
  {
#ifdef HAVE_KOKKOSCLASSIC_MKL
    mkl_zcsrsm ((char*) transa,
                (MKL_INT*) &m,
                (MKL_INT*) &n,
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (&alpha)),
                (char*) matdescra,
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (val)),
                (MKL_INT*) ind,
                (MKL_INT*) ptrBegin,
                (MKL_INT*) ptrEnd,
                reinterpret_cast<MKL_Complex16*> (const_cast<std::complex<double>* > (B)), (MKL_INT*) &LDB,
                reinterpret_cast<MKL_Complex16*> (C), (MKL_INT*) &LDC);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, missingMklMsg);
#endif // HAVE_KOKKOSCLASSIC_MKL
  }

} // namespace Mkl
} // namespace Kokkos

#endif // HAVE_TEUCHOS_COMPLEX


