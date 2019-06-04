/// \file  Tacho_Lapack_External.hpp
/// \brief Lapack wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "ShyLU_NodeTacho_config.h"
#include "Tacho_Lapack_External.hpp"

#include "Kokkos_Core.hpp"

extern "C" {
  ///
  /// POTRF
  ///
  void F77_BLAS_MANGLE(spotrf,SPOTRF)( const char*, 
                                       const int*, 
                                       float*, const int*,
                                       int* );
  void F77_BLAS_MANGLE(dpotrf,DPOTRF)( const char*, 
                                       const int*,
                                       double*, const int*,
                                       int*);
  void F77_BLAS_MANGLE(cpotrf,CPOTRF)( const char*, 
                                       const int*,
                                       Kokkos::complex<float>*, const int*,
                                       int*);
  void F77_BLAS_MANGLE(zpotrf,ZPOTRF)( const char*, 
                                       const int*, 
                                       Kokkos::complex<double>*, const int*,
                                       int*);

  void F77_BLAS_MANGLE(ssytrf,SSYTRF)( const char*, 
                                       const int*, 
                                       float*, const int*,
                                       int*,
                                       float*, int*,
                                       int*);
  void F77_BLAS_MANGLE(dsytrf,DSYTRF)( const char*, 
                                       const int*,
                                       double*, const int*,
                                       int*,
                                       double*, int*,
                                       int*);
  void F77_BLAS_MANGLE(csytrf,CSYTRF)( const char*, 
                                       const int*,
                                       Kokkos::complex<float>*, const int*,
                                       int*,
                                       Kokkos::complex<float>*, int *,
                                       int*);
  void F77_BLAS_MANGLE(zsytrf,ZSYTRF)( const char*, 
                                       const int*, 
                                       Kokkos::complex<double>*, const int*,
                                       int*,
                                       Kokkos::complex<double>*, int*,
                                       int*);
}

namespace Tacho {

#define F77_FUNC_SPOTRF F77_BLAS_MANGLE(spotrf,SPOTRF)
#define F77_FUNC_DPOTRF F77_BLAS_MANGLE(dpotrf,DPOTRF)
#define F77_FUNC_CPOTRF F77_BLAS_MANGLE(cpotrf,CPOTRF)
#define F77_FUNC_ZPOTRF F77_BLAS_MANGLE(zpotrf,ZPOTRF)

#define F77_FUNC_SSYTRF F77_BLAS_MANGLE(ssytrf,SSYTRF)
#define F77_FUNC_DSYTRF F77_BLAS_MANGLE(dsytrf,DSYTRF)
#define F77_FUNC_CSYTRF F77_BLAS_MANGLE(csytrf,CSYTRF)
#define F77_FUNC_ZSYTRF F77_BLAS_MANGLE(zsytrf,ZSYTRF)

    template<>
    void 
    Lapack<float>::potrf(const char uplo,
                         const int m,
                         float *a, const int lda,
                         int *info) {
      F77_FUNC_SPOTRF(&uplo,
                      &m,
                      a, &lda,
                      info);
    }
    template<>
    void 
    Lapack<float>::sytrf(const char uplo,
                         const int m,
                         float *a, const int lda,
                         int *ipiv,
                         float *work, int lwork,
                         int *info) {
      F77_FUNC_SSYTRF(&uplo,
                      &m,
                      a, &lda,
                      ipiv,
                      work, &lwork,
                      info);
    }
    
    template<>
    void
    Lapack<double>::potrf(const char uplo,
                          const int m,
                          double *a, const int lda,
                          int *info) {
      F77_FUNC_DPOTRF(&uplo,
                      &m,
                      a, &lda,
                      info);
    }
    template<>
    void 
    Lapack<double>::sytrf(const char uplo,
                          const int m,
                          double *a, const int lda,
                          int *ipiv,
                          double *work, int lwork,
                          int *info) {
      F77_FUNC_DSYTRF(&uplo,
                      &m,
                      a, &lda,
                      ipiv,
                      work, &lwork,
                      info);
    }

    template<>
    void 
    Lapack<Kokkos::complex<float> >::potrf(const char uplo,
                                           const int m,
                                           Kokkos::complex<float> *a, const int lda,
                                           int *info) {
      F77_FUNC_CPOTRF(&uplo,
                      &m,
                      a, &lda,
                      info);
    }
    
    template<>
    void 
    Lapack<Kokkos::complex<float> >::sytrf(const char uplo,
                                           const int m,
                                           Kokkos::complex<float> *a, const int lda,
                                           int *ipiv,
                                           Kokkos::complex<float> *work, int lwork,
                                           int *info) {
      F77_FUNC_CSYTRF(&uplo,
                      &m,
                      a, &lda,
                      ipiv,
                      work, &lwork,
                      info);
    }

    template<>
    void 
    Lapack<Kokkos::complex<double> >::potrf(const char uplo,
                                            const int m,
                                            Kokkos::complex<double> *a, const int lda,
                                            int *info) {
      F77_FUNC_ZPOTRF(&uplo,
                      &m,
                      a, &lda,
                      info);
    }
    template<>
    void 
    Lapack<Kokkos::complex<double> >::sytrf(const char uplo,
                                            const int m,
                                            Kokkos::complex<double> *a, const int lda,
                                            int *ipiv,
                                            Kokkos::complex<double>* work, int lwork,
                                            int *info) {
      F77_FUNC_ZSYTRF(&uplo,
                      &m,
                      a, &lda,
                      ipiv,
                      work, &lwork,
                      info);
    }

    template<>
    void 
    Lapack<std::complex<float> >::potrf(const char uplo,
                                        const int m,
                                        std::complex<float> *a, const int lda,
                                        int *info) {
      F77_FUNC_CPOTRF(&uplo,
                      &m,
                      (Kokkos::complex<float>*)a, &lda,
                      info);
    }
    
    template<>
    void 
    Lapack<std::complex<float> >::sytrf(const char uplo,
                                        const int m,
                                        std::complex<float> *a, const int lda,
                                        int *ipiv,
                                        std::complex<float> *work, int lwork,
                                        int *info) {
      F77_FUNC_CSYTRF(&uplo,
                      &m,
                      (Kokkos::complex<float>*)a, &lda,
                      ipiv,
                      (Kokkos::complex<float>*)work, &lwork,
                      info);
    }

    template<>
    void 
    Lapack<std::complex<double> >::potrf(const char uplo,
                                         const int m,
                                         std::complex<double> *a, const int lda,
                                         int *info) {
      F77_FUNC_ZPOTRF(&uplo,
                      &m,
                      (Kokkos::complex<double>*)a, &lda,
                      info);
    }
    template<>
    void 
    Lapack<std::complex<double> >::sytrf(const char uplo,
                                         const int m,
                                         std::complex<double> *a, const int lda,
                                         int *ipiv,
                                         std::complex<double>* work, int lwork,
                                         int *info) {
      F77_FUNC_ZSYTRF(&uplo,
                      &m,
                      (Kokkos::complex<double>*)a, &lda,
                      ipiv,
                      (Kokkos::complex<double>*)work, &lwork,
                      info);
    }
    

}
