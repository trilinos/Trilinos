// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <vector>
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_LAPACK_wrappers.hpp"
#ifdef HAVE_TEUCHOSCORE_QUADMATH
#  include "Teuchos_Details_Lapack128.hpp" // impl for __float128
#endif // HAVE_TEUCHOSCORE_QUADMATH
#ifdef HAVE_TEUCHOS_LONG_DOUBLE
#  include "Teuchos_Details_LapackLongDouble.hpp" // impl for long double
#endif // HAVE_TEUCHOS_LONG_DOUBLE
#include "Teuchos_TestForException.hpp"

/* for INTEL_CXML, the second arg may need to be changed to 'one'.  If so
the appropriate declaration of one will need to be added back into
functions that include the macro:
*/

#ifdef CHAR_MACRO
#undef CHAR_MACRO
#endif
#if defined (INTEL_CXML)
#define CHAR_MACRO(char_var) &char_var, one
#else
#define CHAR_MACRO(char_var) &char_var
#endif

#ifdef CHARPTR_MACRO
#undef CHARPR_MACRO
#endif
#if defined (INTEL_CXML)
#define CHARPTR_MACRO(charptr_var) charptr_var, one
#else
#define CHARPTR_MACRO(charptr_var) charptr_var
#endif

namespace {

#if defined (INTEL_CXML)
        unsigned int one=1;
#endif

// Use a wrapper function to handle calling ILAENV().  This removes
// duplicaiton and avoid name lookup problems with member functions called
// ILAENV() trying to call nonmember functions called ILAENV() (which does not
// work on Intel compiler on Windows, see Trilinos bug 5762).
inline
int ilaenv_wrapper(
  const int* ispec, const char* name, const unsigned int& name_length,
  const char* opts, const unsigned int& opts_length,
  const int* N1, const int* N2, const int* N3, const int* N4 )
{
#if defined (INTEL_CXML)
    return ILAENV_F77(ispec, name, name_length, opts, opts_length, N1, N2, N3, N4 );
#else
    return ILAENV_F77(ispec, name, opts, N1, N2, N3, N4, name_length, opts_length );
#endif
}

} // namespace



extern "C" {


typedef int (*gees_nullfptr_t)(double*,double*);


} // extern "C"


namespace Teuchos
{
  // BEGIN INT, FLOAT SPECIALIZATION IMPLEMENTATION //

  void LAPACK<int, float>::PTTRF(const int& n, float* d, float* e, int* info) const
  { SPTTRF_F77(&n,d,e,info); }


  void LAPACK<int, float>::PTTRS(const int& n, const int& nrhs, const float* d, const float* e, float* B, const int& ldb, int* info) const
  { SPTTRS_F77(&n,&nrhs,d,e,B,&ldb,info); }


  void LAPACK<int, float>::POTRF(const char& UPLO, const int& n, float* A, const int& lda, int* info) const
  { SPOTRF_F77(CHAR_MACRO(UPLO), &n, A, &lda, info); }


  void LAPACK<int, float>::POTRS(const char& UPLO, const int& n, const int& nrhs, const float* A, const int& lda, float* B, const int& ldb, int* info) const
  { SPOTRS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, B, &ldb, info); }


  void LAPACK<int, float>::POTRI(const char& UPLO, const int& n, float* A, const int& lda, int* info) const
  { SPOTRI_F77(CHAR_MACRO(UPLO), &n, A, &lda, info); }


  void LAPACK<int, float>::POCON(const char& UPLO, const int& n, const float* A, const int& lda, const float& anorm, float* rcond, float* WORK, int* IWORK, int* info) const
  { SPOCON_F77(CHAR_MACRO(UPLO), &n, A, &lda, &anorm, rcond, WORK, IWORK, info); }


  void LAPACK<int, float>::POSV(const char& UPLO, const int& n, const int& nrhs, float* A, const int& lda, float* B, const int& ldb, int* info) const
  { SPOSV_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, B, &ldb, info); }


  void LAPACK<int, float>::POEQU(const int& n, const float* A, const int& lda, float* S, float* scond, float* amax, int* info) const
  { SPOEQU_F77(&n, A, &lda, S, scond, amax, info); }

  void LAPACK<int, float>::PORFS(const char& UPLO, const int& n, const int& nrhs, const float* A, const int& lda, const float* AF, const int& ldaf, const float* B, const int& ldb, float* X, const int& ldx, float* FERR, float* BERR, float* WORK, int* IWORK, int* info) const
  { SPORFS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, AF, &ldaf, B, &ldb, X, &ldx, FERR, BERR, WORK, IWORK, info); }

  void LAPACK<int, float>::POSVX(const char& FACT, const char& UPLO, const int& n, const int& nrhs, float* A, const int& lda, float* AF, const int& ldaf, char* EQUED, float* S, float* B, const int& ldb, float* X, const int& ldx, float* rcond, float* FERR, float* BERR, float* WORK, int* IWORK, int* info) const
  { SPOSVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, AF, &ldaf, CHARPTR_MACRO(EQUED), S, B, &ldb, X, &ldx, rcond, FERR, BERR, WORK, IWORK, info); }


  void LAPACK<int,float>::GELS(const char& TRANS, const int& m, const int& n, const int& nrhs, float* A, const int& lda, float* B, const int& ldb, float* WORK, const int& lwork, int* info) const
  { SGELS_F77(CHAR_MACRO(TRANS), &m, &n, &nrhs, A, &lda, B, &ldb, WORK, &lwork, info); }

  void LAPACK<int,float>::GELSS (const int& m, const int& n, const int& nrhs, float* A, const int& lda, float* B, const int& ldb, float* S, const float& rcond, int* rank, float* WORK, const int& lwork, float* rwork, int* info) const
  {
    (void) rwork;
    SGELSS_F77(&m, &n, &nrhs, A, &lda, B, &ldb, S, &rcond, rank, WORK, &lwork, info);
  }

  void LAPACK<int,float>::GELSS(const int& m, const int& n, const int& nrhs, float* A, const int& lda, float* B, const int& ldb, float* S, const float& rcond, int* rank, float* WORK, const int& lwork, int* info) const
  { SGELSS_F77(&m, &n, &nrhs, A, &lda, B, &ldb, S, &rcond, rank, WORK, &lwork, info); }


  void LAPACK<int,float>::GGLSE(const int& m, const int& n, const int& p, float* A, const int& lda, float* B, const int& ldb, float* C, float* D, float* X, float* WORK, const int& lwork, int* info) const
  { SGGLSE_F77(&m, &n, &p, A, &lda, B, &ldb, C, D, X, WORK, &lwork, info); }


  void LAPACK<int,float>::GEQRF( const int& m, const int& n, float* A, const int& lda, float* TAU, float* WORK, const int& lwork, int* info) const
  { SGEQRF_F77(&m, &n, A, &lda, TAU, WORK, &lwork, info); }

  void LAPACK<int,float>::GEQR2 (const int& m, const int& n, float* A, const int& lda, float* TAU, float* WORK, int* const info) const
  {
    SGEQR2_F77(&m, &n, A, &lda, TAU, WORK, info);
  }

  void LAPACK<int,float>::GETRF(const int& m, const int& n, float* A, const int& lda, int* IPIV, int* info) const
  { SGETRF_F77(&m, &n, A, &lda, IPIV, info); }


  void LAPACK<int,float>::GETRS(const char& TRANS, const int& n, const int& nrhs, const float* A, const int& lda, const int* IPIV, float* B, const int& ldb, int* info) const
  { SGETRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, IPIV, B, &ldb, info); }


  void LAPACK<int,float>::LASCL(const char& TYPE, const int& kl, const int& ku, const float& cfrom, const float& cto, const int& m, const int& n, float* A, const int& lda, int* info) const
  { SLASCL_F77(CHAR_MACRO(TYPE), &kl, &ku, &cfrom, &cto, &m, &n, A, &lda, info); }

  void LAPACK<int,float>::GEQP3 (const int& m, const int& n, float* A, const int& lda, int* jpvt, float* TAU, float* WORK, const int& lwork, float* RWORK, int* info) const
  {
    (void) RWORK;
    SGEQP3_F77(&m, &n, A, &lda, jpvt, TAU, WORK, &lwork, info);
  }

  void LAPACK<int, float>::LASWP (const int& N, float* A, const int& LDA, const int& K1, const int& K2, const int* IPIV, const int& INCX) const
  {
    SLASWP_F77(&N, A, &LDA, &K1, &K2, IPIV, &INCX);
  }

  void LAPACK<int,float>::GBTRF(const int& m, const int& n, const int& kl, const int& ku, float* A, const int& lda, int* IPIV, int* info) const
  { SGBTRF_F77(&m, &n, &kl, &ku, A, &lda, IPIV, info); }


  void LAPACK<int,float>::GBTRS(const char& TRANS, const int& n, const int& kl, const int& ku, const int& nrhs, const float* A, const int& lda, const int* IPIV, float* B, const int& ldb, int* info) const
  { SGBTRS_F77(CHAR_MACRO(TRANS), &n, &kl, &ku, &nrhs, A, &lda, IPIV, B, &ldb, info); }


  void LAPACK<int,float>::GTTRF(const int& n, float* dl, float* d, float* du, float* du2, int* IPIV, int* info) const
  { SGTTRF_F77(&n, dl, d, du, du2, IPIV, info); }


  void LAPACK<int,float>::GTTRS(const char& TRANS, const int& n, const int& nrhs, const float* dl, const float* d, const float* du, const float* du2, const int* IPIV, float* B, const int& ldb, int* info) const
  { SGTTRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, dl, d, du, du2, IPIV, B, &ldb, info); }


  void LAPACK<int,float>::GETRI(const int& n, float* A, const int& lda, const int* IPIV, float* WORK, const int& lwork, int* info) const
  { SGETRI_F77(&n, A, &lda, IPIV, WORK, &lwork, info); }

  void LAPACK<int, float>::LATRS (const char& UPLO, const char& TRANS, const char& DIAG, const char& NORMIN, const int& N, const float* A, const int& LDA, float* X, float* SCALE, float* CNORM, int* INFO) const
  {
    SLATRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), CHAR_MACRO(NORMIN), &N, A, &LDA, X, SCALE, CNORM, INFO);
  }

  void LAPACK<int,float>::GECON(const char& NORM, const int& n, const float* A, const int& lda, const float& anorm, float* rcond, float* WORK, int* IWORK, int* info) const
  { SGECON_F77(CHAR_MACRO(NORM), &n, A, &lda, &anorm, rcond, WORK, IWORK, info); }


  void LAPACK<int,float>::GBCON(const char& NORM, const int& n, const int& kl, const int& ku, const float* A, const int& lda, const int* IPIV, const float& anorm, float* rcond, float* WORK, int* IWORK, int* info) const
  { SGBCON_F77(CHAR_MACRO(NORM), &n, &kl, &ku, A, &lda, IPIV, &anorm, rcond, WORK, IWORK, info); }


  float LAPACK<int,float>::LANGB(const char& NORM, const int& n, const int& kl, const int& ku, const float* A, const int& lda, float* WORK) const
  { return( SLANGB_F77(CHAR_MACRO(NORM), &n, &kl, &ku, A, &lda, WORK) ); }


  void LAPACK<int,float>::GESV(const int& n, const int& nrhs, float* A, const int& lda, int* IPIV, float* B, const int& ldb, int* info) const
  { SGESV_F77(&n, &nrhs, A, &lda, IPIV, B, &ldb, info); }


  void LAPACK<int,float>::GEEQU(const int& m, const int& n, const float* A, const int& lda, float* R, float* C, float* rowcond, float* colcond, float* amax, int* info) const
  { SGEEQU_F77(&m, &n, A, &lda, R, C, rowcond, colcond, amax, info); }


  void LAPACK<int,float>::GERFS(const char& TRANS, const int& n, const int& nrhs, const float* A, const int& lda, const float* AF, const int& ldaf, const int* IPIV, const float* B, const int& ldb, float* X, const int& ldx, float* FERR, float* BERR, float* WORK, int* IWORK, int* info) const
  { SGERFS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, AF, &ldaf, IPIV, B, &ldb, X, &ldx, FERR, BERR, WORK, IWORK, info); }


  void LAPACK<int,float>::GBEQU(const int& m, const int& n, const int& kl, const int& ku, const float* A, const int& lda, float* R, float* C, float* rowcond, float* colcond, float* amax, int* info) const
  { SGBEQU_F77(&m, &n, &kl, &ku, A, &lda, R, C, rowcond, colcond, amax, info); }


  void LAPACK<int,float>::GBRFS(const char& TRANS, const int& n, const int& kl, const int& ku, const int& nrhs, const float* A, const int& lda, const float* AF, const int& ldaf, const int* IPIV, const float* B, const int& ldb, float* X, const int& ldx, float* FERR, float* BERR, float* WORK, int* IWORK, int* info) const
  { SGBRFS_F77(CHAR_MACRO(TRANS), &n, &kl, &ku, &nrhs, A, &lda, AF, &ldaf, IPIV, B, &ldb, X, &ldx, FERR, BERR, WORK, IWORK, info); }

  void LAPACK<int,float>::GESVX(const char& FACT, const char& TRANS, const int& n, const int& nrhs, float* A, const int& lda, float* AF, const int& ldaf, int* IPIV, char* EQUED, float* R, float* C, float* B, const int& ldb, float* X, const int& ldx, float* rcond, float* FERR, float* BERR, float* WORK, int* IWORK, int* info) const
  { SGESVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, AF, &ldaf, IPIV, CHARPTR_MACRO(EQUED), R, C, B, &ldb, X, &ldx, rcond, FERR, BERR, WORK, IWORK, info); }


  void LAPACK<int,float>::SYTRD(const char& UPLO, const int& n, float* A, const int& lda, float* D, float* E, float* TAU, float* WORK, const int& lwork, int* info) const
  { SSYTRD_F77(CHAR_MACRO(UPLO), &n, A, &lda, D, E, TAU, WORK, &lwork, info); }


  void LAPACK<int,float>::GEHRD(const int& n, const int& ilo, const int& ihi, float* A, const int& lda, float* TAU, float* WORK, const int& lwork, int* info) const
  { SGEHRD_F77(&n, &ilo, &ihi, A, &lda, TAU, WORK, &lwork, info); }


  void LAPACK<int,float>::TRTRS(const char& UPLO, const char& TRANS, const char& DIAG, const int& n, const int& nrhs, const float* A, const int& lda, float* B, const int& ldb, int* info) const
  { STRTRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), &n, &nrhs, A, &lda, B, &ldb, info); }


  void LAPACK<int,float>::TRTRI(const char& UPLO, const char& DIAG, const int& n, float* A, const int& lda, int* info) const
  { STRTRI_F77(CHAR_MACRO(UPLO), CHAR_MACRO(DIAG), &n, A, &lda, info); }


  void LAPACK<int,float>::SPEV(const char& JOBZ, const char& UPLO, const int& n, float* AP, float* W, float* Z, const int& ldz, float* WORK, int* info) const
  { SSPEV_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, AP, W, Z, &ldz, WORK, info); }


  void LAPACK<int,float>::SYEV(const char& JOBZ, const char& UPLO, const int& n, float* A, const int& lda, float* W, float* WORK, const int& lwork, int* info) const
  { SSYEV_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, A, &lda, W, WORK, &lwork, info); }


  void LAPACK<int,float>::SYGV(const int& itype, const char& JOBZ, const char& UPLO, const int& n, float* A, const int& lda, float* B, const int& ldb, float* W, float* WORK, const int& lwork, int* info) const
  { SSYGV_F77(&itype, CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, A, &lda, B, &ldb, W, WORK, &lwork, info); }


  void LAPACK<int,float>::HEEV(const char& JOBZ, const char& UPLO, const int& n, float* A, const int& lda, float* W, float* WORK, const int& lwork, float* /* RWORK */, int* info) const
  { SSYEV_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, A, &lda, W, WORK, &lwork, info); }


  void LAPACK<int,float>::HEGV(const int& itype, const char& JOBZ, const char& UPLO, const int& n, float* A, const int& lda, float* B, const int& ldb, float* W, float* WORK, const int& lwork, float* /* RWORK */, int* info) const
  { SSYGV_F77(&itype, CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, A, &lda, B, &ldb, W, WORK, &lwork, info); }


  void LAPACK<int,float>::STEQR(const char& COMPZ, const int& n, float* D, float* E, float* Z, const int& ldz, float* WORK, int* info) const
  { SSTEQR_F77(CHAR_MACRO(COMPZ), &n, D, E, Z, &ldz, WORK, info); }


  void LAPACK<int,float>::PTEQR(const char& COMPZ, const int& n, float* D, float* E, float* Z, const int& ldz, float* WORK, int* info) const
  { SPTEQR_F77(CHAR_MACRO(COMPZ), &n, D, E, Z, &ldz, WORK, info); }


  void LAPACK<int, float>::HSEQR(const char& JOB, const char& COMPZ, const int& n, const int& ilo, const int& ihi, float* H, const int& ldh, float* WR, float* WI, float* Z, const int& ldz, float* WORK, const int& lwork, int* info) const
  { SHSEQR_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPZ), &n, &ilo, &ihi, H, &ldh, WR, WI, Z, &ldz, WORK, &lwork, info); }


  void LAPACK<int, float>::GEES(const char& JOBVS, const char& SORT, int (*ptr2func)(float*, float*), const int& n, float* A, const int& lda, int* sdim, float* WR, float* WI, float* VS, const int& ldvs, float* WORK, const int& lwork, int* BWORK, int* info) const
  { SGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(SORT), ptr2func, &n, A, &lda, sdim, WR, WI, VS, &ldvs, WORK, &lwork, BWORK, info); }


  void LAPACK<int, float>::GEES(const char& JOBVS, const int& n, float* A, const int& lda, int* sdim, float* WR, float* WI, float* VS, const int& ldvs, float* WORK, const int& lwork, float* /* RWORK */, int* BWORK, int* info) const
  {
    int (*nullfptr)(float*,float*) = NULL;
    const char sort = 'N';
    SGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(sort), nullfptr, &n, A, &lda, sdim, WR, WI, VS, &ldvs, WORK, &lwork, BWORK, info);
  }


  void LAPACK<int, float>::GEEV(const char& JOBVL, const char& JOBVR, const int& n, float* A, const int& lda, float* WR, float* WI, float* VL, const int& ldvl, float* VR, const int& ldvr, float* WORK, const int& lwork, int* info) const
  { SGEEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &n, A, &lda, WR, WI, VL, &ldvl, VR, &ldvr, WORK, &lwork, info); }

  void LAPACK<int, float>::GEEV(const char& JOBVL, const char& JOBVR, const int& n, float* A, const int& lda, float* WR, float* WI, float* VL, const int& ldvl, float* VR, const int& ldvr, float* WORK, const int& lwork, float* /* RWORK */, int* info) const
  {
    GEEV (JOBVL, JOBVR, n, A, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, lwork, info);
  }


  void LAPACK<int, float>::GESVD(const char& JOBU, const char& JOBVT, const int& m, const int& n, float* A, const int& lda, float* S, float* U, const int& ldu, float* V, const int& ldv, float* WORK, const int& lwork, float* /* RWORK */, int* info) const
  { SGESVD_F77(CHAR_MACRO(JOBU), CHAR_MACRO(JOBVT), &m, &n, A, &lda, S, U, &ldu, V, &ldv, WORK, &lwork, info); }


  void LAPACK<int,float>::GEEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, float* A, const int& lda, float* WR, float* WI, float* VL, const int& ldvl, float* VR, const int& ldvr, int* ilo, int* ihi, float* SCALE, float* abnrm, float* RCONDE, float* RCONDV, float* WORK, const int& lwork, int* IWORK, int* info) const
  { SGEEVX_F77(CHAR_MACRO(BALANC), CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SENSE), &n, A, &lda, WR, WI, VL, &ldvl, VR, &ldvr, ilo, ihi, SCALE, abnrm, RCONDE, RCONDV, WORK, &lwork, IWORK, info); }


  void LAPACK<int,float>::GGEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, float* A, const int& lda, float* B, const int& ldb, float* ALPHAR, float* ALPHAI, float* BETA, float* VL, const int& ldvl, float* VR, const int& ldvr, int* ilo, int* ihi, float* lscale, float* rscale, float* abnrm, float* bbnrm, float* RCONDE, float* RCONDV, float* WORK, const int& lwork, int* IWORK, int* BWORK, int* info) const
  { SGGEVX_F77(CHAR_MACRO(BALANC), CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SENSE), &n, A, &lda, B, &ldb, ALPHAR, ALPHAI, BETA, VL, &ldvl, VR, &ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, RCONDE, RCONDV, WORK, &lwork, IWORK, BWORK, info); }

  void LAPACK<int,float>::GGEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, float* A, const int& lda, float* B, const int& ldb, float* ALPHAR, float* ALPHAI, float* BETA, float* VL, const int& ldvl, float* VR, const int& ldvr, int* ilo, int* ihi, float* lscale, float* rscale, float* abnrm, float* bbnrm, float* RCONDE, float* RCONDV, float* WORK, const int& lwork, float* /* RWORK */, int* IWORK, int* BWORK, int* info) const
  {
    GGEVX(BALANC, JOBVL, JOBVR, SENSE, n, A, lda, B, ldb, ALPHAR, ALPHAI, BETA, VL, ldvl, VR, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, RCONDE, RCONDV, WORK, lwork, IWORK, BWORK, info);
  }

  void LAPACK<int, float>::GGEV(const char& JOBVL, const char& JOBVR, const int& n, float* A, const int& lda, float* B, const int& ldb, float* ALPHAR, float* ALPHAI, float* BETA, float* VL, const int& ldvl, float* VR, const int& ldvr, float* WORK, const int& lwork, int* info) const
  { SGGEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &n, A, &lda, B, &ldb, ALPHAR, ALPHAI, BETA, VL, &ldvl, VR, &ldvr, WORK, &lwork, info); }


  void LAPACK<int, float>::TRSEN(const char& JOB, const char& COMPQ, const int* SELECT, const int& n, float* T, const int& ldt, float* Q, const int& ldq, float* WR, float* WI, int* M, float* S, float* SEP, float* WORK, const int& lwork, int* IWORK, const int& liwork, int* info ) const
  { STRSEN_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPQ), SELECT, &n, T, &ldt, Q, &ldq, WR, WI, M, S, SEP, WORK, &lwork, IWORK, &liwork, info); }


  void LAPACK<int, float>::TGSEN(const int& ijob, const int& wantq, const int& wantz, const int* SELECT, const int& n, float* A, const int& lda, float* B, const int& ldb, float* ALPHAR, float* ALPHAI, float* BETA, float* Q, const int& ldq, float* Z, const int& ldz, int* M, float* PL, float* PR, float* DIF, float* WORK, const int& lwork, int* IWORK, const int& liwork, int* info ) const
  { STGSEN_F77(&ijob, &wantq, &wantz, SELECT, &n, A, &lda, B, &ldb, ALPHAR, ALPHAI, BETA, Q, &ldq, Z, &ldz, M, PL, PR, DIF, WORK, &lwork, IWORK, &liwork, info); }


  void LAPACK<int, float>::GGES(const char& JOBVL, const char& JOBVR, const char& SORT, int (*ptr2func)(float* , float* , float* ), const int& n, float* A, const int& lda, float* B, const int& ldb, int* sdim, float* ALPHAR, float* ALPHAI, float* BETA, float* VL, const int& ldvl, float* VR, const int& ldvr, float* WORK, const int& lwork, int* BWORK, int* info ) const
  { SGGES_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SORT), ptr2func, &n, A, &lda, B, &ldb, sdim, ALPHAR, ALPHAI, BETA, VL, &ldvl, VR, &ldvr, WORK, &lwork, BWORK, info); }


  void LAPACK<int, float>::ORMQR(const char& SIDE, const char& TRANS, const int& m, const int& n, const int& k, const float* A, const int& lda, const float* TAU, float* C, const int& ldc, float* WORK, const int& lwork, int* info) const
  { SORMQR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &m, &n, &k, A, &lda, TAU, C, &ldc, WORK, &lwork, info); }


  void LAPACK<int, float>::ORM2R(const char& SIDE, const char& TRANS, const int& m, const int& n, const int& k, const float* A, const int& lda, const float* TAU, float* C, const int& ldc, float* WORK, int* const info) const
  { SORM2R_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &m, &n, &k, A, &lda, TAU, C, &ldc, WORK, info); }


  void LAPACK<int, float>::UNMQR(const char& SIDE, const char& TRANS, const int& m, const int& n, const int& k, const float* A, const int& lda, const float* TAU, float* C, const int& ldc, float* WORK, const int& lwork, int* info) const
  {
    // LAPACK only defines UNMQR for Z and C (complex*8 resp. complex*16), but logically, UNMQR means the same thing as ORMQR for real arithmetic.
    ORMQR (SIDE, TRANS, m, n, k, A, lda, TAU, C, ldc, WORK, lwork, info);
  }

  void LAPACK<int, float>::UNM2R (const char& SIDE, const char& TRANS, const int& M, const int& N, const int& K, const float* A, const int& LDA, const float* TAU, float* C, const int& LDC, float* WORK, int* const INFO) const
  {
    // LAPACK only defines UNM2R for Z and C (complex*8
    // resp. complex*16), but logically, UNM2R means the same thing as
    // ORM2R for real arithmetic.
    ORM2R (SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO);
  }


  void LAPACK<int, float>::ORGQR(const int& m, const int& n, const int& k, float* A, const int& lda, const float* TAU, float* WORK, const int& lwork, int* info) const
  { SORGQR_F77( &m, &n, &k, A, &lda, TAU, WORK, &lwork, info); }


  void LAPACK<int, float>::UNGQR(const int& m, const int& n, const int& k, float* A, const int& lda, const float* TAU, float* WORK, const int& lwork, int* info) const
  { SORGQR_F77( &m, &n, &k, A, &lda, TAU, WORK, &lwork, info); }


  void LAPACK<int, float>::ORGHR(const int& n, const int& ilo, const int& ihi, float* A, const int& lda, const float* TAU, float* WORK, const int& lwork, int* info) const
  { SORGHR_F77(&n, &ilo, &ihi, A, &lda, TAU, WORK, &lwork, info); }


  void LAPACK<int, float>::ORMHR(const char& SIDE, const char& TRANS, const int& m, const int& n, const int& ilo, const int& ihi, const float* A, const int& lda, const float* TAU, float* C, const int& ldc, float* WORK, const int& lwork, int* info) const
  { SORMHR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &m, &n, &ilo, &ihi, A, &lda, TAU, C, &ldc, WORK, &lwork, info); }


  void LAPACK<int, float>::TREVC(const char& SIDE, const char& HOWMNY, int* select, const int& n, const float* T, const int& ldt, float* VL, const int& ldvl, float* VR, const int& ldvr, const int& mm, int* m, float* WORK, int* info) const
  { STREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), select, &n, T, &ldt, VL, &ldvl, VR, &ldvr, &mm, m, WORK, info); }


  void LAPACK<int, float>::TREVC(const char& SIDE, const int& n, const float* T, const int& ldt, float* VL, const int& ldvl, float* VR, const int& ldvr, const int& mm, int* m, float* WORK, float* /* RWORK */, int* info) const
  {
    std::vector<int> select(1);
    const char whch = 'A';
    STREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(whch), &select[0], &n, T, &ldt, VL, &ldvl, VR, &ldvr, &mm, m, WORK, info);
  }

  void LAPACK<int, float>::TREXC(const char& COMPQ, const int& n, float* T, const int& ldt, float* Q, const int& ldq, int* ifst, int* ilst, float* WORK, int* info) const
  { STREXC_F77(CHAR_MACRO(COMPQ), &n, T, &ldt, Q, &ldq, ifst, ilst, WORK, info); }


  void LAPACK<int, float>::TGEVC(const char& SIDE, const char& HOWMNY, const int* SELECT, const int& n, const float* S, const int& lds, const float* P, const int& ldp, float* VL, const int& ldvl, float* VR, const int& ldvr, const int& mm, int* M, float* WORK, int* info) const
  { STGEVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), SELECT, &n, S, &lds, P, &ldp, VL, &ldvl, VR, &ldvr, &mm, M, WORK, info); }


  void LAPACK<int, float>::LARTG( const float& f, const float& g, float* c, float* s, float* r ) const
  { SLARTG_F77(&f, &g, c, s, r); }


  void LAPACK<int, float>::LARFG( const int& n, float* alpha, float* x, const int& incx, float* tau ) const
  { SLARFG_F77(&n, alpha, x, &incx, tau); }

  void LAPACK<int, float>::GEBAL(const char& JOBZ, const int& n, float* A, const int& lda, int* ilo, int* ihi, float* scale, int* info) const
  { SGEBAL_F77(CHAR_MACRO(JOBZ),&n, A, &lda, ilo, ihi, scale, info); }


  void LAPACK<int, float>::GEBAK(const char& JOBZ, const char& SIDE, const int& n, const int& ilo, const int& ihi, const float* scale, const int& m, float* V, const int& ldv, int* info) const
  { SGEBAK_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(SIDE), &n, &ilo, &ihi, scale, &m, V, &ldv, info); }

#ifdef HAVE_TEUCHOS_LAPACKLARND
  float LAPACK<int, float>::LARND( const int& idist, int* seed ) const
  { return(SLARND_F77(&idist, seed)); }
#endif

  void LAPACK<int, float>::LARNV( const int& idist, int* seed, const int& n, float* v ) const
  { SLARNV_F77(&idist, seed, &n, v); }


  float LAPACK<int, float>::LAMCH(const char& CMACH) const
  { return(SLAMCH_F77(CHAR_MACRO(CMACH))); }


  int LAPACK<int, float>::ILAENV( const int& ispec, const std::string& NAME, const std::string& OPTS, const int& N1, const int& N2, const int& N3, const int& N4 ) const
  {
    unsigned int opts_length = OPTS.length();
    // if user queries a Hermitian routine, change it to a symmetric routine
    std::string temp_NAME = "s" + NAME;
    if (temp_NAME.substr(1,2) == "he") {
      temp_NAME.replace(1,2,"sy");
    }
    unsigned int name_length = temp_NAME.length();
    return ilaenv_wrapper(&ispec, &temp_NAME[0], name_length, &OPTS[0], opts_length, &N1, &N2, &N3, &N4);
  }


  float LAPACK<int, float>::LAPY2(const float& x, const float& y) const
  {
#if defined(HAVE_TEUCHOS_BLASFLOAT)
    return SLAPY2_F77(&x, &y);
#else
    typedef ScalarTraits<float> ST;
    const float xabs = ST::magnitude(x);
    const float yabs = ST::magnitude(y);
    const float w = TEUCHOS_MAX(xabs, yabs);
    const float z = TEUCHOS_MIN(xabs, yabs);
    if ( z == 0.0 ) {
      return w;
    }
    const float z_over_w = z/w;
    return w*ST::squareroot( 1.0+(z_over_w*z_over_w));
#endif
  }

  // END INT, FLOAT SPECIALIZATION IMPLEMENTATION //

  // BEGIN INT, DOUBLE SPECIALIZATION IMPLEMENTATION //

  void LAPACK<int, double>::PTTRF(const int& n, double* d, double* e, int* info) const
  { DPTTRF_F77(&n,d,e,info); }


  void LAPACK<int, double>::PTTRS(const int& n, const int& nrhs, const double* d, const double* e, double* B, const int& ldb, int* info) const
  { DPTTRS_F77(&n,&nrhs,d,e,B,&ldb,info); }


  void LAPACK<int, double>::POTRF(const char& UPLO, const int& n, double* A, const int& lda, int* info) const
  { DPOTRF_F77(CHAR_MACRO(UPLO), &n, A, &lda, info); }


  void LAPACK<int, double>::POTRS(const char& UPLO, const int& n, const int& nrhs, const double* A, const int& lda, double* B, const int& ldb, int* info) const
  { DPOTRS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, B, &ldb, info); }


  void LAPACK<int, double>::POTRI(const char& UPLO, const int& n, double* A, const int& lda, int* info) const
  { DPOTRI_F77(CHAR_MACRO(UPLO), &n, A, &lda, info); }


  void LAPACK<int, double>::POCON(const char& UPLO, const int& n, const double* A, const int& lda, const double& anorm, double* rcond, double* WORK, int* IWORK, int* info) const
  { DPOCON_F77(CHAR_MACRO(UPLO), &n, A, &lda, &anorm, rcond, WORK, IWORK, info); }


  void LAPACK<int, double>::POSV(const char& UPLO, const int& n, const int& nrhs, double* A, const int& lda, double* B, const int& ldb, int* info) const
  { DPOSV_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, B, &ldb, info); }


  void LAPACK<int, double>::POEQU(const int& n, const double* A, const int& lda, double* S, double* scond, double* amax, int* info) const
  { DPOEQU_F77(&n, A, &lda, S, scond, amax, info); }


  void LAPACK<int, double>::PORFS(const char& UPLO, const int& n, const int& nrhs, const double* A, const int& lda, const double* AF, const int& ldaf, const double* B, const int& ldb, double* X, const int& ldx, double* FERR, double* BERR, double* WORK, int* IWORK, int* info) const
  { DPORFS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, AF, &ldaf, B, &ldb, X, &ldx, FERR, BERR, WORK, IWORK, info); }

  void LAPACK<int, double>::POSVX(const char& FACT, const char& UPLO, const int& n, const int& nrhs, double* A, const int& lda, double* AF, const int& ldaf, char* EQUED, double* S, double* B, const int& ldb, double* X, const int& ldx, double* rcond, double* FERR, double* BERR, double* WORK, int* IWORK, int* info) const
  { DPOSVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, AF, &ldaf, CHARPTR_MACRO(EQUED), S, B, &ldb, X, &ldx, rcond, FERR, BERR, WORK, IWORK, info); }


  void LAPACK<int,double>::GELS(const char& TRANS, const int& m, const int& n, const int& nrhs, double* A, const int& lda, double* B, const int& ldb, double* WORK, const int& lwork, int* info) const
  { DGELS_F77(CHAR_MACRO(TRANS), &m, &n, &nrhs, A, &lda, B, &ldb, WORK, &lwork, info); }


  void LAPACK<int,double>::GELSS(const int& m, const int& n, const int& nrhs, double* A, const int& lda, double* B, const int& ldb, double* S, const double& rcond, int* rank, double* WORK, const int& lwork, double* rwork, int* info) const
  {
    (void) rwork;
    DGELSS_F77(&m, &n, &nrhs, A, &lda, B, &ldb, S, &rcond, rank, WORK, &lwork, info);
  }


  void LAPACK<int,double>::GELSS(const int& m, const int& n, const int& nrhs, double* A, const int& lda, double* B, const int& ldb, double* S, const double& rcond, int* rank, double* WORK, const int& lwork, int* info) const
  { DGELSS_F77(&m, &n, &nrhs, A, &lda, B, &ldb, S, &rcond, rank, WORK, &lwork, info); }


  void LAPACK<int,double>::GGLSE(const int& m, const int& n, const int& p, double* A, const int& lda, double* B, const int& ldb, double* C, double* D, double* X, double* WORK, const int& lwork, int* info) const
  { DGGLSE_F77(&m, &n, &p, A, &lda, B, &ldb, C, D, X, WORK, &lwork, info); }


  void LAPACK<int,double>::GEQRF( const int& m, const int& n, double* A, const int& lda, double* TAU, double* WORK, const int& lwork, int* info) const
  { DGEQRF_F77(&m, &n, A, &lda, TAU, WORK, &lwork, info); }

  void LAPACK<int,double>::GEQR2 (const int& m, const int& n, double* A, const int& lda, double* TAU, double* WORK, int* const info) const
  {
    DGEQR2_F77(&m, &n, A, &lda, TAU, WORK, info);
  }

  void LAPACK<int,double>::GETRF(const int& m, const int& n, double* A, const int& lda, int* IPIV, int* info) const
  { DGETRF_F77(&m, &n, A, &lda, IPIV, info); }


  void LAPACK<int,double>::GETRS(const char& TRANS, const int& n, const int& nrhs, const double* A, const int& lda, const int* IPIV, double* B, const int& ldb, int* info) const
  { DGETRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, IPIV, B, &ldb, info); }


  void LAPACK<int,double>::LASCL(const char& TYPE, const int& kl, const int& ku, const double& cfrom, const double& cto, const int& m, const int& n, double* A, const int& lda, int* info) const
  { DLASCL_F77(CHAR_MACRO(TYPE), &kl, &ku, &cfrom, &cto, &m, &n, A, &lda, info); }

  void LAPACK<int,double>::GEQP3(const int& m, const int& n, double* A, const int& lda, int* jpvt, double* TAU, double* WORK, const int& lwork, double* RWORK, int* info ) const
  {
    (void) RWORK;
    DGEQP3_F77(&m, &n, A, &lda, jpvt, TAU, WORK, &lwork, info);
  }

  void LAPACK<int, double>::LASWP (const int& N, double* A, const int& LDA, const int& K1, const int& K2, const int* IPIV, const int& INCX) const
  { DLASWP_F77(&N, A, &LDA, &K1, &K2, IPIV, &INCX); }

  void LAPACK<int,double>::GBTRF(const int& m, const int& n, const int& kl, const int& ku, double* A, const int& lda, int* IPIV, int* info) const
  { DGBTRF_F77(&m, &n, &kl, &ku, A, &lda, IPIV, info); }


  void LAPACK<int,double>::GBTRS(const char& TRANS, const int& n, const int& kl, const int& ku, const int& nrhs, const double* A, const int& lda, const int* IPIV, double* B, const int& ldb, int* info) const
  { DGBTRS_F77(CHAR_MACRO(TRANS), &n, &kl, &ku, &nrhs, A, &lda, IPIV, B, &ldb, info); }


  void LAPACK<int,double>::GTTRF(const int& n, double* dl, double* d, double* du, double* du2, int* IPIV, int* info) const
  { DGTTRF_F77(&n, dl, d, du, du2, IPIV, info); }


  void LAPACK<int,double>::GTTRS(const char& TRANS, const int& n, const int& nrhs, const double* dl, const double* d, const double* du, const double* du2, const int* IPIV, double* B, const int& ldb, int* info) const
  { DGTTRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, dl, d, du, du2, IPIV, B, &ldb, info); }


  void LAPACK<int,double>::GETRI(const int& n, double* A, const int& lda, const int* IPIV, double* WORK, const int& lwork, int* info) const
  { DGETRI_F77(&n, A, &lda, IPIV, WORK, &lwork, info); }

  void LAPACK<int, double>::LATRS (const char& UPLO, const char& TRANS, const char& DIAG, const char& NORMIN, const int& N, const double* A, const int& LDA, double* X, double* SCALE, double* CNORM, int* INFO) const
  {
    DLATRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), CHAR_MACRO(NORMIN), &N, A, &LDA, X, SCALE, CNORM, INFO);
  }

  void LAPACK<int,double>::GECON(const char& NORM, const int& n, const double* A, const int& lda, const double& anorm, double* rcond, double* WORK, int* IWORK, int* info) const
  { DGECON_F77(CHAR_MACRO(NORM), &n, A, &lda, &anorm, rcond, WORK, IWORK, info); }


  void LAPACK<int,double>::GBCON(const char& NORM, const int& n, const int& kl, const int& ku, const double* A, const int& lda, const int* IPIV, const double& anorm, double* rcond, double* WORK, int* IWORK, int* info) const
  { DGBCON_F77(CHAR_MACRO(NORM), &n, &kl, &ku, A, &lda, IPIV, &anorm, rcond, WORK, IWORK, info); }


  double LAPACK<int,double>::LANGB(const char& NORM, const int& n, const int& kl, const int& ku, const double* A, const int& lda, double* WORK) const
  { return( DLANGB_F77(CHAR_MACRO(NORM), &n, &kl, &ku, A, &lda, WORK) ); }


  void LAPACK<int,double>::GESV(const int& n, const int& nrhs, double* A, const int& lda, int* IPIV, double* B, const int& ldb, int* info) const
  { DGESV_F77(&n, &nrhs, A, &lda, IPIV, B, &ldb, info); }


  void LAPACK<int,double>::GEEQU(const int& m, const int& n, const double* A, const int& lda, double* R, double* C, double* rowcond, double* colcond, double* amax, int* info) const
  { DGEEQU_F77(&m, &n, A, &lda, R, C, rowcond, colcond, amax, info); }


  void LAPACK<int,double>::GERFS(const char& TRANS, const int& n, const int& nrhs, const double* A, const int& lda, const double* AF, const int& ldaf, const int* IPIV, const double* B, const int& ldb, double* X, const int& ldx, double* FERR, double* BERR, double* WORK, int* IWORK, int* info) const
  { DGERFS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, AF, &ldaf, IPIV, B, &ldb, X, &ldx, FERR, BERR, WORK, IWORK, info); }


  void LAPACK<int,double>::GBEQU(const int& m, const int& n, const int& kl, const int& ku, const double* A, const int& lda, double* R, double* C, double* rowcond, double* colcond, double* amax, int* info) const
  { DGBEQU_F77(&m, &n, &kl, &ku, A, &lda, R, C, rowcond, colcond, amax, info); }


  void LAPACK<int,double>::GBRFS(const char& TRANS, const int& n, const int& kl, const int& ku, const int& nrhs, const double* A, const int& lda, const double* AF, const int& ldaf, const int* IPIV, const double* B, const int& ldb, double* X, const int& ldx, double* FERR, double* BERR, double* WORK, int* IWORK, int* info) const
  { DGBRFS_F77(CHAR_MACRO(TRANS), &n, &kl, &ku, &nrhs, A, &lda, AF, &ldaf, IPIV, B, &ldb, X, &ldx, FERR, BERR, WORK, IWORK, info); }

  void LAPACK<int,double>::GESVX(const char& FACT, const char& TRANS, const int& n, const int& nrhs, double* A, const int& lda, double* AF, const int& ldaf, int* IPIV, char* EQUED, double* R, double* C, double* B, const int& ldb, double* X, const int& ldx, double* rcond, double* FERR, double* BERR, double* WORK, int* IWORK, int* info) const
  { DGESVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, AF, &ldaf, IPIV, CHARPTR_MACRO(EQUED), R, C, B, &ldb, X, &ldx, rcond, FERR, BERR, WORK, IWORK, info); }


  void LAPACK<int,double>::SYTRD(const char& UPLO, const int& n, double* A, const int& lda, double* D, double* E, double* TAU, double* WORK, const int& lwork, int* info) const
  { DSYTRD_F77(CHAR_MACRO(UPLO), &n, A, &lda, D, E, TAU, WORK, &lwork, info); }


  void LAPACK<int, double>::GEHRD(const int& n, const int& ilo, const int& ihi, double* A, const int& lda, double* TAU, double* WORK, const int& lwork, int* info) const
  { DGEHRD_F77(&n, &ilo, &ihi, A, &lda, TAU, WORK, &lwork, info); }


  void LAPACK<int,double>::TRTRS(const char& UPLO, const char& TRANS, const char& DIAG, const int& n, const int& nrhs, const double* A, const int& lda, double* B, const int& ldb, int* info) const
  { DTRTRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), &n, &nrhs, A, &lda, B, &ldb, info); }


  void LAPACK<int,double>::TRTRI(const char& UPLO, const char& DIAG, const int& n, double* A, const int& lda, int* info) const
  { DTRTRI_F77(CHAR_MACRO(UPLO), CHAR_MACRO(DIAG), &n, A, &lda, info); }


  void LAPACK<int,double>::SPEV(const char& JOBZ, const char& UPLO, const int& n, double* AP, double* W, double* Z, const int& ldz, double* WORK, int* info) const
  { DSPEV_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, AP, W, Z, &ldz, WORK, info); }


  void LAPACK<int,double>::SYEV(const char& JOBZ, const char& UPLO, const int& n, double* A, const int& lda, double* W, double* WORK, const int& lwork, int* info) const
  {
    DSYEV_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, A, &lda, W, WORK, &lwork, info);
  }


  void LAPACK<int,double>::SYGV(const int& itype, const char& JOBZ, const char& UPLO, const int& n, double* A, const int& lda, double* B, const int& ldb, double* W, double* WORK, const int& lwork, int* info) const
  {
    DSYGV_F77(&itype, CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, A, &lda, B, &ldb, W, WORK, &lwork, info);
  }


  void LAPACK<int,double>::HEEV(const char& JOBZ, const char& UPLO, const int& n, double* A, const int& lda, double* W, double* WORK, const int& lwork, double* /* RWORK */, int* info) const
  {
    DSYEV_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, A, &lda, W, WORK, &lwork, info);
  }


  void LAPACK<int,double>::HEGV(const int& itype, const char& JOBZ, const char& UPLO, const int& n, double* A, const int& lda, double* B, const int& ldb, double* W, double* WORK, const int& lwork, double* /* RWORK */, int* info) const
  {
    DSYGV_F77(&itype, CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, A, &lda, B, &ldb, W, WORK, &lwork, info);
  }


  void LAPACK<int,double>::STEQR(const char& COMPZ, const int& n, double* D, double* E, double* Z, const int& ldz, double* WORK, int* info) const
  { DSTEQR_F77(CHAR_MACRO(COMPZ), &n, D, E, Z, &ldz, WORK, info); }

  
  void LAPACK<int,double>::PTEQR(const char& COMPZ, const int& n, double* D, double* E, double* Z, const int& ldz, double* WORK, int* info) const
  { DPTEQR_F77(CHAR_MACRO(COMPZ), &n, D, E, Z, &ldz, WORK, info); }


  void LAPACK<int, double>::HSEQR(const char& JOB, const char& COMPZ, const int& n, const int& ilo, const int& ihi, double* H, const int& ldh, double* WR, double* WI, double* Z, const int& ldz, double* WORK, const int& lwork, int* info) const
  {
    DHSEQR_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPZ), &n, &ilo, &ihi, H, &ldh, WR, WI, Z, &ldz, WORK, &lwork, info);
  }


  void LAPACK<int, double>::GEES(const char& JOBVS, const char& SORT, int (*ptr2func)(double*, double*), const int& n, double* A, const int& lda, int* sdim, double* WR, double* WI, double* VS, const int& ldvs, double* WORK, const int& lwork, int* BWORK, int* info) const
  {
    DGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(SORT), ptr2func, &n, A, &lda, sdim, WR, WI, VS, &ldvs, WORK, &lwork, BWORK, info);
  }


  void LAPACK<int, double>::GEES(const char& JOBVS, const int& n, double* A, const int& lda, int* sdim, double* WR, double* WI, double* VS, const int& ldvs, double* WORK, const int& lwork, double* /* RWORK */, int* BWORK, int* info) const
  {
    //int (*nullfptr)(double*,double*) = NULL;
    gees_nullfptr_t nullfptr = 0;
    const char sort = 'N';
    DGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(sort), nullfptr, &n, A, &lda, sdim, WR, WI, VS, &ldvs, WORK, &lwork, BWORK, info);
  }


  void LAPACK<int, double>::GEEV(const char& JOBVL, const char& JOBVR, const int& n, double* A, const int& lda, double* WR, double* WI, double* VL, const int& ldvl, double* VR, const int& ldvr, double* WORK, const int& lwork, int* info) const
  {
    DGEEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &n, A, &lda, WR, WI, VL, &ldvl, VR, &ldvr, WORK, &lwork, info);
  }

  void LAPACK<int, double>::GEEV(const char& JOBVL, const char& JOBVR, const int& n, double* A, const int& lda, double* WR, double* WI, double* VL, const int& ldvl, double* VR, const int& ldvr, double* WORK, const int& lwork, double* /* RWORK */, int* info) const
  {
    GEEV (JOBVL, JOBVR, n, A, lda, WR, WI, VL, ldvl, VR, ldvr, WORK, lwork, info);
  }


  void LAPACK<int, double>::GESVD(const char& JOBU, const char& JOBVT, const int& m, const int& n, double* A, const int& lda, double* S, double* U, const int& ldu, double* V, const int& ldv, double* WORK, const int& lwork, double* /* RWORK */, int* info) const {
    DGESVD_F77(CHAR_MACRO(JOBU), CHAR_MACRO(JOBVT), &m, &n, A, &lda, S, U, &ldu, V, &ldv, WORK, &lwork, info);
  }


  void LAPACK<int,double>::GEEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, double* A, const int& lda, double* WR, double* WI, double* VL, const int& ldvl, double* VR, const int& ldvr, int* ilo, int* ihi, double* SCALE, double* abnrm, double* RCONDE, double* RCONDV, double* WORK, const int& lwork, int* IWORK, int* info) const
  {
    DGEEVX_F77(CHAR_MACRO(BALANC), CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SENSE), &n, A, &lda, WR, WI, VL, &ldvl, VR, &ldvr, ilo, ihi, SCALE, abnrm, RCONDE, RCONDV, WORK, &lwork, IWORK, info);
  }


  void LAPACK<int, double>::GGEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, double* A, const int& lda, double* B, const int& ldb, double* ALPHAR, double* ALPHAI, double* BETA, double* VL, const int& ldvl, double* VR, const int& ldvr, int* ilo, int* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, double* RCONDE, double* RCONDV, double* WORK, const int& lwork, int* IWORK, int* BWORK, int* info) const
  {
    DGGEVX_F77(CHAR_MACRO(BALANC), CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SENSE), &n, A, &lda, B, &ldb, ALPHAR, ALPHAI, BETA, VL, &ldvl, VR, &ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, RCONDE, RCONDV, WORK, &lwork, IWORK, BWORK, info);
  }

  void LAPACK<int, double>::GGEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, double* A, const int& lda, double* B, const int& ldb, double* ALPHAR, double* ALPHAI, double* BETA, double* VL, const int& ldvl, double* VR, const int& ldvr, int* ilo, int* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, double* RCONDE, double* RCONDV, double* WORK, const int& lwork, double* /* RWORK */, int* IWORK, int* BWORK, int* info) const
  {
    GGEVX(BALANC, JOBVL, JOBVR, SENSE, n, A, lda, B, ldb, ALPHAR, ALPHAI, BETA, VL, ldvl, VR, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, RCONDE, RCONDV, WORK, lwork, IWORK, BWORK, info);
  }

  void LAPACK<int, double>::GGEV(const char& JOBVL, const char& JOBVR, const int& n, double* A, const int& lda, double* B, const int& ldb, double* ALPHAR, double* ALPHAI, double* BETA, double* VL, const int& ldvl, double* VR, const int& ldvr, double* WORK, const int& lwork, int* info) const
  {
    DGGEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &n, A, &lda, B, &ldb, ALPHAR, ALPHAI, BETA, VL, &ldvl, VR, &ldvr, WORK, &lwork, info);
  }

  void LAPACK<int, double>::TRSEN(const char& JOB, const char& COMPQ, const int* SELECT, const int& n, double* T, const int& ldt, double* Q, const int& ldq, double* WR, double* WI, int* M, double* S, double* SEP, double* WORK, const int& lwork, int* IWORK, const int& liwork, int* info ) const
  { DTRSEN_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPQ), SELECT, &n, T, &ldt, Q, &ldq, WR, WI, M, S, SEP, WORK, &lwork, IWORK, &liwork, info); }


  void LAPACK<int, double>::TGSEN(const int& ijob, const int& wantq, const int& wantz, const int* SELECT, const int& n, double* A, const int& lda, double* B, const int& ldb, double* ALPHAR, double* ALPHAI, double* BETA, double* Q, const int& ldq, double* Z, const int& ldz, int* M, double* PL, double* PR, double* DIF, double* WORK, const int& lwork, int* IWORK, const int& liwork, int* info ) const
  { DTGSEN_F77(&ijob, &wantq, &wantz, SELECT, &n, A, &lda, B, &ldb, ALPHAR, ALPHAI, BETA, Q, &ldq, Z, &ldz, M, PL, PR, DIF, WORK, &lwork, IWORK, &liwork, info); }


  void LAPACK<int, double>::GGES(const char& JOBVL, const char& JOBVR, const char& SORT, int (*ptr2func)(double* , double* , double* ), const int& n, double* A, const int& lda, double* B, const int& ldb, int* sdim, double* ALPHAR, double* ALPHAI, double* BETA, double* VL, const int& ldvl, double* VR, const int& ldvr, double* WORK, const int& lwork, int* BWORK, int* info ) const
  { DGGES_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SORT), ptr2func, &n, A, &lda, B, &ldb, sdim, ALPHAR, ALPHAI, BETA, VL, &ldvl, VR, &ldvr, WORK, &lwork, BWORK, info); }


  void LAPACK<int, double>::ORMQR(const char& SIDE, const char& TRANS, const int& m, const int& n, const int& k, const double* A, const int& lda, const double* TAU, double* C, const int& ldc, double* WORK, const int& lwork, int* info) const
  {
    DORMQR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &m, &n, &k, A, &lda, TAU, C, &ldc, WORK, &lwork, info);
  }

  void LAPACK<int, double>::ORM2R(const char& SIDE, const char& TRANS, const int& m, const int& n, const int& k, const double* A, const int& lda, const double* TAU, double* C, const int& ldc, double* WORK, int* const info) const
  {
    DORM2R_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &m, &n, &k, A, &lda, TAU, C, &ldc, WORK, info);
  }

  void LAPACK<int, double>::UNMQR(const char& SIDE, const char& TRANS, const int& m, const int& n, const int& k, const double* A, const int& lda, const double* TAU, double* C, const int& ldc, double* WORK, const int& lwork, int* info) const
  {
    // LAPACK only defines UNMQR for Z and C (complex*8 resp. complex*16), but logically, UNMQR means the same thing as ORMQR for real arithmetic.
    ORMQR (SIDE, TRANS, m, n, k, A, lda, TAU, C, ldc, WORK, lwork, info);
  }

  void LAPACK<int, double>::UNM2R (const char& SIDE, const char& TRANS, const int& M, const int& N, const int& K, const double* A, const int& LDA, const double* TAU, double* C, const int& LDC, double* WORK, int* const INFO) const
  {
    // LAPACK only defines UNM2R for Z and C (complex*8
    // resp. complex*16), but logically, UNM2R means the same thing as
    // ORM2R for real arithmetic.
    ORM2R (SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO);
  }

  void LAPACK<int, double>::ORGQR(const int& m, const int& n, const int& k, double* A, const int& lda, const double* TAU, double* WORK, const int& lwork, int* info) const
  {
    DORGQR_F77( &m, &n, &k, A, &lda, TAU, WORK, &lwork, info);
  }


  void LAPACK<int, double>::UNGQR(const int& m, const int& n, const int& k, double* A, const int& lda, const double* TAU, double* WORK, const int& lwork, int* info) const
  {
    DORGQR_F77( &m, &n, &k, A, &lda, TAU, WORK, &lwork, info);
  }


  void LAPACK<int, double>::ORGHR(const int& n, const int& ilo, const int& ihi, double* A, const int& lda, const double* TAU, double* WORK, const int& lwork, int* info) const
  {
    DORGHR_F77(&n, &ilo, &ihi, A, &lda, TAU, WORK, &lwork, info);
  }


  void LAPACK<int, double>::ORMHR(const char& SIDE, const char& TRANS, const int& m, const int& n, const int& ilo, const int& ihi, const double* A, const int& lda, const double* TAU, double* C, const int& ldc, double* WORK, const int& lwork, int* info) const
  {
    DORMHR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &m, &n, &ilo, &ihi, A, &lda, TAU, C, &ldc, WORK, &lwork, info);
  }


  void LAPACK<int, double>::TREVC(const char& SIDE, const char& HOWMNY, int* select, const int& n, const double* T, const int& ldt, double* VL, const int& ldvl, double* VR, const int& ldvr, const int& mm, int* m, double* WORK, int* info) const
  {
    DTREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), select, &n, T, &ldt, VL, &ldvl, VR, &ldvr, &mm, m, WORK, info);
  }


  void LAPACK<int, double>::TREVC(const char& SIDE, const int& n, const double* T, const int& ldt, double* VL, const int& ldvl, double* VR, const int& ldvr, const int& mm, int* m, double* WORK, double* /* RWORK */, int* info) const
  {
    std::vector<int> select(1);
    const char whch = 'A';
    DTREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(whch), &select[0], &n, T, &ldt, VL, &ldvl, VR, &ldvr, &mm, m, WORK, info);
  }

  void LAPACK<int, double>::TREXC(const char& COMPQ, const int& n, double* T, const int& ldt, double* Q, const int& ldq, int* ifst, int* ilst, double* WORK, int* info) const
  {
    DTREXC_F77(CHAR_MACRO(COMPQ), &n, T, &ldt, Q, &ldq, ifst, ilst, WORK, info);
  }


  void LAPACK<int, double>::TGEVC(const char& SIDE, const char& HOWMNY, const int* SELECT, const int& n, const double* S, const int& lds, const double* P, const int& ldp, double* VL, const int& ldvl, double* VR, const int& ldvr, const int& mm, int* M, double* WORK, int* info) const
  { DTGEVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), SELECT, &n, S, &lds, P, &ldp, VL, &ldvl, VR, &ldvr, &mm, M, WORK, info); }


  void LAPACK<int, double>::LARTG( const double& f, const double& g, double* c, double* s, double* r ) const
  {
    DLARTG_F77(&f, &g, c, s, r);
  }


  void LAPACK<int, double>::LARFG( const int& n, double* alpha, double* x, const int& incx, double* tau ) const
  {
    DLARFG_F77(&n, alpha, x, &incx, tau);
  }

  void LAPACK<int, double>::GEBAL(const char& JOBZ, const int& n, double* A, const int& lda, int* ilo, int* ihi, double* scale, int* info) const
  {
    DGEBAL_F77(CHAR_MACRO(JOBZ),&n, A, &lda, ilo, ihi, scale, info);
  }


  void LAPACK<int, double>::GEBAK(const char& JOBZ, const char& SIDE, const int& n, const int& ilo, const int& ihi, const double* scale, const int& m, double* V, const int& ldv, int* info) const
  {
    DGEBAK_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(SIDE), &n, &ilo, &ihi, scale, &m, V, &ldv, info);
  }


#ifdef HAVE_TEUCHOS_LAPACKLARND
  double LAPACK<int, double>::LARND( const int& idist, int* seed ) const
  {
    return(DLARND_F77(&idist, seed));
  }
#endif

  void LAPACK<int, double>::LARNV( const int& idist, int* seed, const int& n, double* v ) const
  {
    DLARNV_F77(&idist, seed, &n, v);
  }


  double LAPACK<int, double>::LAMCH(const char& CMACH) const
  {
    return(DLAMCH_F77(CHAR_MACRO(CMACH)));
  }


  int LAPACK<int, double>::ILAENV( const int& ispec, const std::string& NAME, const std::string& OPTS, const int& N1, const int& N2, const int& N3, const int& N4 ) const
  {
    unsigned int opts_length = OPTS.length();
    // if user queries a Hermitian routine, change it to a symmetric routine
    std::string temp_NAME = "d" + NAME;
    if (temp_NAME.substr(1,2) == "he") {
      temp_NAME.replace(1,2,"sy");
    }
    unsigned int name_length = temp_NAME.length();
    return ilaenv_wrapper(&ispec, &temp_NAME[0], name_length, &OPTS[0], opts_length, &N1, &N2, &N3, &N4);
  }


  double LAPACK<int, double>::LAPY2(const double& x, const double& y) const
  {
    return DLAPY2_F77(&x, &y);
  }

  // END INT, DOUBLE SPECIALIZATION IMPLEMENTATION //

#ifdef HAVE_TEUCHOS_COMPLEX

  // BEGIN INT, COMPLEX<FLOAT> SPECIALIZATION IMPLEMENTATION //


  void LAPACK<int, std::complex<float> >::PTTRF(const int& n, float* d, std::complex<float>* e, int* info) const
  {
    CPTTRF_F77(&n,d,e,info);
  }


  void LAPACK<int, std::complex<float> >::PTTRS(const char& UPLO, const int& n, const int& nrhs, const float* d, const std::complex<float>* e, std::complex<float>* B, const int& ldb, int* info) const
  {
    CPTTRS_F77(CHAR_MACRO(UPLO),&n,&nrhs,d,e,B,&ldb,info);
  }


  void LAPACK<int, std::complex<float> >::POTRF(const char& UPLO, const int& n, std::complex<float>* A, const int& lda, int* info) const
  {
    CPOTRF_F77(CHAR_MACRO(UPLO), &n, A, &lda, info);
  }


  void LAPACK<int, std::complex<float> >::POTRS(const char& UPLO, const int& n, const int& nrhs, const std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, int* info) const
  {
    CPOTRS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, B, &ldb, info);
  }


  void LAPACK<int, std::complex<float> >::POTRI(const char& UPLO, const int& n, std::complex<float>* A, const int& lda, int* info) const
  {
    CPOTRI_F77(CHAR_MACRO(UPLO), &n, A, &lda, info);
  }


  void LAPACK<int, std::complex<float> >::POCON(const char& UPLO, const int& n, const std::complex<float>* A, const int& lda, const float& anorm, float* rcond, std::complex<float>* WORK, float* RWORK, int* info) const
  {
    CPOCON_F77(CHAR_MACRO(UPLO), &n, A, &lda, &anorm, rcond, WORK, RWORK, info);
  }


  void LAPACK<int, std::complex<float> >::POSV(const char& UPLO, const int& n, const int& nrhs, std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, int* info) const
  {
    CPOSV_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, B, &ldb, info);
  }


  void LAPACK<int, std::complex<float> >::POEQU(const int& n, const std::complex<float>* A, const int& lda, float* S, float* scond, float* amax, int* info) const
  {
    CPOEQU_F77(&n, A, &lda, S, scond, amax, info);
  }


  void LAPACK<int, std::complex<float> >::PORFS(const char& UPLO, const int& n, const int& nrhs, const std::complex<float>* A, const int& lda, const std::complex<float>* AF, const int& ldaf, const std::complex<float>* B, const int& ldb, std::complex<float>* X, const int& ldx, float* FERR, float* BERR, std::complex<float>* WORK, float* RWORK, int* info) const
  {
    CPORFS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, AF, &ldaf, B, &ldb, X, &ldx, FERR, BERR, WORK, RWORK, info);
  }

  void LAPACK<int, std::complex<float> >::POSVX(const char& FACT, const char& UPLO, const int& n, const int& nrhs, std::complex<float>* A, const int& lda, std::complex<float>* AF, const int& ldaf, char* EQUED, float* S, std::complex<float>* B, const int& ldb, std::complex<float>* X, const int& ldx, float* rcond, float* FERR, float* BERR, std::complex<float>* WORK, float* RWORK, int* info) const
  {
    CPOSVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, AF, &ldaf, CHARPTR_MACRO(EQUED), S, B, &ldb, X, &ldx, rcond, FERR, BERR, WORK, RWORK, info);
  }


  void LAPACK<int,std::complex<float> >::GELS(const char& TRANS, const int& m, const int& n, const int& nrhs, std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, std::complex<float>* WORK, const int& lwork, int* info) const
  {
    CGELS_F77(CHAR_MACRO(TRANS), &m, &n, &nrhs, A, &lda, B, &ldb, WORK, &lwork, info);
  }

  void LAPACK<int, std::complex<float> >::GELSS(const int& m, const int& n, const int& nrhs, std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, float* S, const float& rcond, int* rank, std::complex<float>* WORK, const int& lwork, float* rwork, int* info) const
  {
    CGELSS_F77(&m, &n, &nrhs, A, &lda, B, &ldb, S, &rcond, rank, WORK, &lwork, rwork, info);
  }

  void LAPACK<int,std::complex<float> >::GEQRF( const int& m, const int& n, std::complex<float>* A, const int& lda, std::complex<float>* TAU, std::complex<float>* WORK, const int& lwork, int* info) const
  {
    CGEQRF_F77(&m, &n, A, &lda, TAU, WORK, &lwork, info);
  }

  void LAPACK<int,std::complex<float> >::GEQR2 (const int& m, const int& n, std::complex<float>* A, const int& lda, std::complex<float>* TAU, std::complex<float>* WORK, int* const info) const
  {
    CGEQR2_F77(&m, &n, A, &lda, TAU, WORK, info);
  }

  void LAPACK<int,std::complex<float> >::UNGQR(const int& m, const int& n, const int& k, std::complex<float>* A, const int& lda, const std::complex<float>* TAU, std::complex<float>* WORK, const int& lwork, int* info) const
  {
    CUNGQR_F77( &m, &n, &k, A, &lda, TAU, WORK, &lwork, info);
  }

  void LAPACK<int,std::complex<float> >::UNMQR(const char& SIDE, const char& TRANS, const int& m, const int& n, const int& k, const std::complex<float>* A, const int& lda, const std::complex<float>* TAU, std::complex<float>* C, const int& ldc, std::complex<float>* WORK, const int& lwork, int* info) const
  {
    CUNMQR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &m, &n, &k, A, &lda, TAU, C, &ldc, WORK, &lwork, info);
  }

  void LAPACK<int,std::complex<float> >::UNM2R (const char& SIDE, const char& TRANS, const int& M, const int& N, const int& K, const std::complex<float>* A, const int& LDA, const std::complex<float>* TAU, std::complex<float>* C, const int& LDC, std::complex<float>* WORK, int* const INFO) const
  {
    CUNM2R_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, INFO);
  }

  void LAPACK<int,std::complex<float> >::GETRF(const int& m, const int& n, std::complex<float>* A, const int& lda, int* IPIV, int* info) const
  {
    CGETRF_F77(&m, &n, A, &lda, IPIV, info);
  }

  void LAPACK<int,std::complex<float> >::GETRS(const char& TRANS, const int& n, const int& nrhs, const std::complex<float>* A, const int& lda, const int* IPIV, std::complex<float>* B , const int& ldb, int* info) const
  {
    CGETRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, IPIV, B, &ldb, info);
  }

  void LAPACK<int,std::complex<float> >::LASCL(const char& TYPE, const int& kl, const int& ku, const float& cfrom, const float& cto, const int& m, const int& n, std::complex<float>* A, const int& lda, int* info) const
  { CLASCL_F77(CHAR_MACRO(TYPE), &kl, &ku, &cfrom, &cto, &m, &n, A, &lda, info); }

  void LAPACK<int,std::complex<float> >::GEQP3(const int& m, const int& n, std::complex<float>* A, const int& lda, int* jpvt, std::complex<float>* TAU, std::complex<float>* WORK, const int& lwork, float* RWORK, int* info ) const
  {
    CGEQP3_F77(&m, &n, A, &lda, jpvt, TAU, WORK, &lwork, RWORK, info);
  }

  void LAPACK<int, std::complex<float> >::
  LASWP (const int& N,
         std::complex<float>* A,
         const int& LDA,
         const int& K1,
         const int& K2,
         const int* IPIV,
         const int& INCX) const
  {
    CLASWP_F77(&N, A, &LDA, &K1, &K2, IPIV, &INCX);
  }

  void LAPACK<int,std::complex<float> >::GBTRF(const int& m, const int& n, const int& kl, const int& ku, std::complex<float>* A, const int& lda, int* IPIV, int* info) const
  {
    CGBTRF_F77(&m, &kl, &ku, &n, A, &lda, IPIV, info);
  }


  void LAPACK<int,std::complex<float> >::GBTRS(const char& TRANS, const int& n, const int& kl, const int& ku, const int& nrhs, const std::complex<float>* A, const int& lda, const int* IPIV, std::complex<float>* B , const int& ldb, int* info) const
  {
    CGBTRS_F77(CHAR_MACRO(TRANS), &n, &kl, &ku, &nrhs, A, &lda, IPIV, B, &ldb, info);
  }


  void LAPACK<int,std::complex<float> >::GTTRF(const int& n, std::complex<float>* dl, std::complex<float>* d, std::complex<float>* du, std::complex<float>* du2, int* IPIV, int* info) const
  {
    CGTTRF_F77(&n, dl, d, du, du2, IPIV, info);
  }


  void LAPACK<int,std::complex<float> >::GTTRS(const char& TRANS, const int& n, const int& nrhs, const std::complex<float>* dl, const std::complex<float>* d, const std::complex<float>* du, const std::complex<float>* du2, const int* IPIV, std::complex<float>* B, const int& ldb, int* info) const
  {
    CGTTRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, dl, d, du, du2, IPIV, B, &ldb, info);
  }


  void LAPACK<int,std::complex<float> >::GETRI(const int& n, std::complex<float>* A, const int& lda, const int* IPIV, std::complex<float>* WORK, const int& lwork, int* info) const
  {
    CGETRI_F77(&n, A, &lda, IPIV, WORK, &lwork, info);
  }


  void LAPACK<int, std::complex<float> >::LATRS (const char& UPLO, const char& TRANS, const char& DIAG, const char& NORMIN, const int& N, const std::complex<float>* A, const int& LDA, std::complex<float>* X, float* SCALE, float* CNORM, int* INFO) const
  {
    CLATRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), CHAR_MACRO(NORMIN), &N, A, &LDA, X, SCALE, CNORM, INFO);
  }


  void LAPACK<int,std::complex<float> >::GECON(const char& NORM, const int& n, const std::complex<float>* A, const int& lda, const float& anorm, float* rcond, std::complex<float>* WORK, float* RWORK, int* info) const
  {
    CGECON_F77(CHAR_MACRO(NORM), &n, A, &lda, &anorm, rcond, WORK, RWORK, info);
  }


  void LAPACK<int,std::complex<float> >::GBCON(const char& NORM, const int& n, const int& kl, const int& ku, const std::complex<float>* A, const int& lda, const int* IPIV, const float& anorm, float* rcond, std::complex<float>* WORK, float* RWORK, int* info) const
  {
    CGBCON_F77(CHAR_MACRO(NORM), &n, &kl, &ku, A, &lda, IPIV, &anorm, rcond, WORK, RWORK, info);
  }


  float LAPACK<int,std::complex<float> >::LANGB(const char& NORM, const int& n, const int& kl, const int& ku, const std::complex<float>* A, const int& lda, float* WORK) const
  {
    return( CLANGB_F77(CHAR_MACRO(NORM), &n, &kl, &ku, A, &lda, WORK) );
  }


  void LAPACK<int,std::complex<float> >::GESV(const int& n, const int& nrhs, std::complex<float>* A, const int& lda, int* IPIV, std::complex<float>* B, const int& ldb, int* info) const
  {
    CGESV_F77(&n, &nrhs, A, &lda, IPIV, B, &ldb, info);
  }


  void LAPACK<int,std::complex<float> >::GEEQU(const int& m, const int& n, const std::complex<float>* A, const int& lda, float* R, float* C, float* rowcond, float* colcond, float* amax, int* info) const
  {
    CGEEQU_F77(&m, &n, A, &lda, R, C, rowcond, colcond, amax, info);
  }


  void LAPACK<int,std::complex<float> >::GERFS(const char& TRANS, const int& n, const int& nrhs, const std::complex<float>* A, const int& lda, const std::complex<float>* AF, const int& ldaf, const int* IPIV, const std::complex<float>* B, const int& ldb, std::complex<float>* X, const int& ldx, float* FERR, float* BERR, std::complex<float>* WORK, float* RWORK, int* info) const
  {
    CGERFS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, AF, &ldaf, IPIV, B, &ldb, X, &ldx, FERR, BERR, WORK, RWORK, info);
  }


  void LAPACK<int,std::complex<float> >::GBEQU(const int& m, const int& n, const int& kl, const int& ku, const std::complex<float>* A, const int& lda, float* R, float* C, float* rowcond, float* colcond, float* amax, int* info) const
  {
    CGBEQU_F77(&m, &n, &kl, &ku, A, &lda, R, C, rowcond, colcond, amax, info);
  }


  void LAPACK<int,std::complex<float> >::GBRFS(const char& TRANS, const int& n, const int& kl, const int& ku, const int& nrhs, const std::complex<float>* A, const int& lda, const std::complex<float>* AF, const int& ldaf, const int* IPIV, const std::complex<float>* B, const int& ldb, std::complex<float>* X, const int& ldx, float* FERR, float* BERR, std::complex<float>* WORK, float* RWORK, int* info) const
  {
    CGBRFS_F77(CHAR_MACRO(TRANS), &n, &kl, &ku, &nrhs, A, &lda, AF, &ldaf, IPIV, B, &ldb, X, &ldx, FERR, BERR, WORK, RWORK, info);
  }

  void LAPACK<int,std::complex<float> >::GESVX(const char& FACT, const char& TRANS, const int& n, const int& nrhs, std::complex<float>* A, const int& lda, std::complex<float>* AF, const int& ldaf, int* IPIV, char* EQUED, float* R, float* C, std::complex<float>* B, const int& ldb, std::complex<float>* X, const int& ldx, float* rcond, float* FERR, float* BERR, std::complex<float>* WORK, float* RWORK, int* info) const
  {
    CGESVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, AF, &ldaf, IPIV, CHARPTR_MACRO(EQUED), R, C, B, &ldb, X, &ldx, rcond, FERR, BERR, WORK, RWORK, info);
  }


  void LAPACK<int,std::complex<float> >::GEHRD(const int& n, const int& ilo, const int& ihi, std::complex<float>* A, const int& lda, std::complex<float>* TAU, std::complex<float>* WORK, const int& lwork, int* info) const
  {
    CGEHRD_F77(&n, &ilo, &ihi, A, &lda, TAU, WORK, &lwork, info);
  }


  void LAPACK<int,std::complex<float> >::TRTRS(const char& UPLO, const char& TRANS, const char& DIAG, const int& n, const int& nrhs, const std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, int* info) const
  {
    CTRTRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), &n, &nrhs, A, &lda, B, &ldb, info);
  }


  void LAPACK<int,std::complex<float> >::TRTRI(const char& UPLO, const char& DIAG, const int& n, std::complex<float>* A, const int& lda, int* info) const
  {
    CTRTRI_F77(CHAR_MACRO(UPLO), CHAR_MACRO(DIAG), &n, A, &lda, info);
  }


  void LAPACK<int,std::complex<float> >::STEQR(const char& COMPZ, const int& n, float* D, float* E, std::complex<float>* Z, const int& ldz, float* WORK, int* info) const
  {
    CSTEQR_F77(CHAR_MACRO(COMPZ), &n, D, E, Z, &ldz, WORK, info);
  }


  void LAPACK<int,std::complex<float> >::PTEQR(const char& COMPZ, const int& n, float* D, float* E, std::complex<float>* Z, const int& ldz, float* WORK, int* info) const
  {
    CPTEQR_F77(CHAR_MACRO(COMPZ), &n, D, E, Z, &ldz, WORK, info);
  }


  void LAPACK<int,std::complex<float> >::HEEV(const char& JOBZ, const char& UPLO, const int& n, std::complex<float> * A, const int& lda, float*  W, std::complex<float> * WORK, const int& lwork, float* RWORK, int* info) const
  {
    CHEEV_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, A, &lda, W, WORK, &lwork, RWORK, info);
  }


  void LAPACK<int,std::complex<float> >::HEGV(const int& itype, const char& JOBZ, const char& UPLO, const int& n, std::complex<float> * A, const int& lda, std::complex<float> * B, const int& ldb, float*  W, std::complex<float> * WORK, const int& lwork, float* RWORK, int* info) const
  {
    CHEGV_F77(&itype, CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, A, &lda, B, &ldb, W, WORK, &lwork, RWORK, info);
  }


  void LAPACK<int, std::complex<float> >::HSEQR(const char& JOB, const char& COMPZ, const int& n, const int& ilo, const int& ihi, std::complex<float>* H, const int& ldh, std::complex<float>* W, std::complex<float>* Z, const int& ldz, std::complex<float>* WORK, const int& lwork, int* info) const
  {
    CHSEQR_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPZ), &n, &ilo, &ihi, H, &ldh, W, Z, &ldz, WORK, &lwork, info);
  }


  void LAPACK<int, std::complex<float> >::GEES(const char& JOBVS, const char& SORT, int (*ptr2func)(std::complex<float>*), const int& n, std::complex<float>* A, const int& lda, int* sdim, std::complex<float>* W, std::complex<float>* VS, const int& ldvs, std::complex<float>* WORK, const int& lwork, float* RWORK, int* BWORK, int* info) const
  {
    CGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(SORT), ptr2func, &n, A, &lda, sdim, W, VS, &ldvs, WORK, &lwork, RWORK, BWORK, info);
  }


  void LAPACK<int, std::complex<float> >::GEES(const char& JOBVS, const int& n, std::complex<float>* A, const int& lda, int* sdim, float* WR, float* WI, std::complex<float>* VS, const int& ldvs, std::complex<float>* WORK, const int& lwork, float* RWORK, int* BWORK, int* info) const
  {
    int (*nullfptr)(std::complex<float>*) = NULL;
    std::vector< std::complex<float> > W(n);
    const char sort = 'N';
    CGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(sort), nullfptr, &n, A, &lda, sdim, &W[0], VS, &ldvs, WORK, &lwork, RWORK, BWORK, info);
    for (int i=0; i<n; i++) {
      WR[i] = W[i].real();
      WI[i] = W[i].imag();
    }
  }


  void LAPACK<int, std::complex<float> >::GEEV(const char& JOBVL, const char& JOBVR, const int& n, std::complex<float>* A, const int& lda, std::complex<float>* W, std::complex<float>* VL, const int& ldvl, std::complex<float>* VR, const int& ldvr, std::complex<float>* WORK, const int& lwork, float* RWORK, int* info) const
  {
    CGEEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &n, A, &lda, W, VL, &ldvl, VR, &ldvr, WORK, &lwork, RWORK, info);
  }

  void LAPACK<int, std::complex<float> >::GEEV(const char& JOBVL, const char& JOBVR, const int& n, std::complex<float>* A, const int& lda, float* WR, float* WI, std::complex<float>* VL, const int& ldvl, std::complex<float>* VR, const int& ldvr, std::complex<float>* WORK, const int& lwork, float* RWORK, int* info) const
  {
    std::vector<std::complex<float> > w (n);
    std::complex<float>* w_rawPtr = (n == 0) ? NULL : &w[0];
    GEEV (JOBVL, JOBVR, n, A, lda, w_rawPtr, VL, ldvl, VR, ldvr, WORK, lwork, RWORK, info);
    if (*info == 0) {
      // The eigenvalues are only valid on output if INFO is zero.
      // Otherwise, we shouldn't even write to WR or WI.
      for (int k = 0; k < n; ++k) {
        WR[k] = w[k].real ();
        WI[k] = w[k].imag ();
      }
    }
  }

  void LAPACK<int, std::complex<float> >::GESVD(const char& JOBU, const char& JOBVT, const int& m, const int& n, std::complex<float> * A, const int& lda, float* S, std::complex<float> * U, const int& ldu, std::complex<float> * V, const int& ldv, std::complex<float> * WORK, const int& lwork, float* RWORK, int* info) const {
    CGESVD_F77(CHAR_MACRO(JOBU), CHAR_MACRO(JOBVT), &m, &n, A, &lda, S, U, &ldu, V, &ldv, WORK, &lwork, RWORK, info);
  }


  void LAPACK<int, std::complex<float> >::GEEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, std::complex<float>* A, const int& lda, std::complex<float>* W, std::complex<float>* VL, const int& ldvl, std::complex<float>* VR, const int& ldvr, int* ilo, int* ihi, float* SCALE, float* abnrm, float* RCONDE, float* RCONDV, std::complex<float>* WORK, const int& lwork, float* RWORK, int* info) const
  {
    CGEEVX_F77(CHAR_MACRO(BALANC), CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SENSE), &n, A, &lda, W, VL, &ldvl, VR, &ldvr, ilo, ihi, SCALE, abnrm, RCONDE, RCONDV, WORK, &lwork, RWORK, info);
  }


  void LAPACK<int, std::complex<float> >::GGEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, std::complex<float>* ALPHA, std::complex<float>* BETA, std::complex<float>* VL, const int& ldvl, std::complex<float>* VR, const int& ldvr, int* ilo, int* ihi, float* lscale, float* rscale, float* abnrm, float* bbnrm, float* RCONDE, float* RCONDV, std::complex<float>* WORK, const int& lwork, float* RWORK, int* IWORK, int* BWORK, int* info) const
  {
    CGGEVX_F77(CHAR_MACRO(BALANC), CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SENSE), &n, A, &lda, B, &ldb, ALPHA, BETA, VL, &ldvl, VR, &ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, RCONDE, RCONDV, WORK, &lwork, RWORK, IWORK, BWORK, info);
  }

  void LAPACK<int, std::complex<float> >::GGEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, float* ALPHAR, float* ALPHAI, std::complex<float>* BETA, std::complex<float>* VL, const int& ldvl, std::complex<float>* VR, const int& ldvr, int* ilo, int* ihi, float* lscale, float* rscale, float* abnrm, float* bbnrm, float* RCONDE, float* RCONDV, std::complex<float>* WORK, const int& lwork, float* RWORK, int* IWORK, int* BWORK, int* info) const
  {
    std::vector<std::complex<float> > w (n);
    std::complex<float>* w_rawPtr = (n == 0) ? NULL : &w[0];
    GGEVX(BALANC, JOBVL, JOBVR, SENSE, n, A, lda, B, ldb, w_rawPtr, BETA, VL, ldvl, VR, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, RCONDE, RCONDV, WORK, lwork, RWORK, IWORK, BWORK, info);
    if (*info == 0) {
      // The eigenvalues are only valid on output if INFO is zero.
      // Otherwise, we shouldn't even write to WR or WI.
      for (int k = 0; k < n; ++k) {
        ALPHAR[k] = w[k].real ();
        ALPHAI[k] = w[k].imag ();
      }
    }
  }

  void LAPACK<int, std::complex<float> >::GGEV(const char& JOBVL, const char& JOBVR, const int& n, std::complex<float> *A, const int& lda, std::complex<float> *B, const int& ldb, std::complex<float>* ALPHA, std::complex<float>* BETA, std::complex<float>* VL, const int& ldvl, std::complex<float>* VR, const int& ldvr, std::complex<float> *WORK, const int& lwork, float* RWORK, int* info) const
  {
    CGGEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &n, A, &lda, B, &ldb, ALPHA, BETA, VL, &ldvl, VR, &ldvr, WORK, &lwork, RWORK, info); 
  }

  void LAPACK<int, std::complex<float> >::GGES(const char& JOBVL, const char& JOBVR, const char& SORT, int (*ptr2func)(std::complex<float>*, std::complex<float>*), const int& n, std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, int* sdim, std::complex<float>* ALPHA, std::complex<float>* BETA, std::complex<float>* VL, const int& ldvl, std::complex<float>* VR, const int& ldvr, std::complex<float>* WORK, const int& lwork, float* rwork, int* bwork, int* info ) const
  {
    CGGES_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SORT), ptr2func, &n, A, &lda, B, &ldb, sdim, ALPHA, BETA, VL, &ldvl, VR, &ldvr, WORK, &lwork, rwork, bwork, info); 
  }
  void LAPACK<int, std::complex<float> >::TGSEN(const int& ijob, const int& wantq, const int& wantz, const int* SELECT, const int& n, std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, std::complex<float>* ALPHA, std::complex<float>* BETA, std::complex<float>* Q, const int& ldq, std::complex<float>* Z, const int& ldz, int* M, float* PL, float* PR, float* DIF, std::complex<float>* WORK, const int& lwork, int* IWORK, const int& liwork, int* info ) const
  {
    CTGSEN_F77(&ijob, &wantq, &wantz, SELECT, &n, A, &lda, B, &ldb, ALPHA, BETA, Q, &ldq, Z, &ldz, M, PL, PR, DIF, WORK, &lwork, IWORK, &liwork, info); 
  }

  void LAPACK<int, std::complex<float> >::TREVC(const char& SIDE, const char& HOWMNY, int* select, const int& n, const std::complex<float>* T, const int& ldt, std::complex<float>* VL, const int& ldvl, std::complex<float>* VR, const int& ldvr, const int& mm, int* m, std::complex<float>* WORK, float* RWORK, int* info) const
  {
    CTREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), select, &n, T, &ldt, VL, &ldvl, VR, &ldvr, &mm, m, WORK, RWORK, info);
  }


  void LAPACK<int, std::complex<float> >::TREVC(const char& SIDE, const int& n, const std::complex<float>* T, const int& ldt, std::complex<float>* VL, const int& ldvl, std::complex<float>* VR, const int& ldvr, const int& mm, int* m, std::complex<float>* WORK, float* RWORK, int* info) const
  {
    std::vector<int> select(1);
    const char& whch = 'A';
    CTREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(whch), &select[0], &n, T, &ldt, VL, &ldvl, VR, &ldvr, &mm, m, WORK, RWORK, info);
  }

  void LAPACK<int, std::complex<float> >::TREXC(const char& COMPQ, const int& n, std::complex<float>* T, const int& ldt, std::complex<float>* Q, const int& ldq, int* ifst, int* ilst, std::complex<float>* WORK, int* info) const
  {
    CTREXC_F77(CHAR_MACRO(COMPQ), &n, T, &ldt, Q, &ldq, ifst, ilst, info);
  }


  void LAPACK<int, std::complex<float> >::LARTG( const std::complex<float> f, const std::complex<float> g, float* c, std::complex<float>* s, std::complex<float>* r ) const
  {
    CLARTG_F77(&f, &g, c, s, r);
  }


  void LAPACK<int, std::complex<float> >::LARFG( const int& n, std::complex<float>* alpha, std::complex<float>* x, const int& incx, std::complex<float>* tau ) const
  {
    CLARFG_F77(&n, alpha, x, &incx, tau);
  }

  void LAPACK<int, std::complex<float> >::GEBAL(const char& JOBZ, const int& n, std::complex<float>* A, const int& lda, int* ilo, int* ihi, float* scale, int* info) const
  {
    CGEBAL_F77(CHAR_MACRO(JOBZ),&n, A, &lda, ilo, ihi, scale, info);
  }


  void LAPACK<int, std::complex<float> >::GEBAK(const char& JOBZ, const char& SIDE, const int& n, const int& ilo, const int& ihi, const float* scale, const int& m, std::complex<float>* V, const int& ldv, int* info) const
  {
    CGEBAK_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(SIDE), &n, &ilo, &ihi, scale, &m, V, &ldv, info);
  }


#ifdef HAVE_TEUCHOS_LAPACKLARND
  std::complex<float> LAPACK<int, std::complex<float> >::LARND( const int& idist, int* seed ) const
  {
    float _Complex z = CLARND_F77(&idist, seed);
    return TEUCHOS_LAPACK_CONVERT_COMPLEX_FORTRAN_TO_CXX(float, z);
  }
#endif

  void LAPACK<int, std::complex<float> >::LARNV( const int& idist, int* seed, const int& n, std::complex<float>* v ) const
  {
    CLARNV_F77(&idist, seed, &n, v);
  }


  int LAPACK<int, std::complex<float> >::ILAENV( const int& ispec, const std::string& NAME, const std::string& OPTS, const int& N1, const int& N2, const int& N3, const int& N4 ) const
  {
    unsigned int opts_length = OPTS.length();
    std::string temp_NAME = "c" + NAME;
    unsigned int name_length = temp_NAME.length();
    return ilaenv_wrapper(&ispec, &temp_NAME[0], name_length, &OPTS[0], opts_length, &N1, &N2, &N3, &N4);
  }

  // END INT, COMPLEX<FLOAT> SPECIALIZATION IMPLEMENTATION //

  // BEGIN INT, COMPLEX<DOUBLE> SPECIALIZATION IMPLEMENTATION //


  void LAPACK<int, std::complex<double> >::PTTRF(const int& n, double* d, std::complex<double>* e, int* info) const
  {
    ZPTTRF_F77(&n,d,e,info);
  }


  void LAPACK<int, std::complex<double> >::PTTRS(const char& UPLO, const int& n, const int& nrhs, const double* d, const std::complex<double>* e, std::complex<double>* B, const int& ldb, int* info) const
  {
    ZPTTRS_F77(CHAR_MACRO(UPLO),&n,&nrhs,d,e,B,&ldb,info);
  }


  void LAPACK<int, std::complex<double> >::POTRF(const char& UPLO, const int& n, std::complex<double>* A, const int& lda, int* info) const
  {
    ZPOTRF_F77(CHAR_MACRO(UPLO), &n, A, &lda, info);
  }


  void LAPACK<int, std::complex<double> >::POTRS(const char& UPLO, const int& n, const int& nrhs, const std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, int* info) const
  {
    ZPOTRS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, B, &ldb, info);
  }


  void LAPACK<int, std::complex<double> >::POTRI(const char& UPLO, const int& n, std::complex<double>* A, const int& lda, int* info) const
  {
    ZPOTRI_F77(CHAR_MACRO(UPLO), &n, A, &lda, info);
  }


  void LAPACK<int, std::complex<double> >::POCON(const char& UPLO, const int& n, const std::complex<double>* A, const int& lda, const double& anorm, double* rcond, std::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZPOCON_F77(CHAR_MACRO(UPLO), &n, A, &lda, &anorm, rcond, WORK, RWORK, info);
  }


  void LAPACK<int, std::complex<double> >::POSV(const char& UPLO, const int& n, const int& nrhs, std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, int* info) const
  {
    ZPOSV_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, B, &ldb, info);
  }


  void LAPACK<int, std::complex<double> >::POEQU(const int& n, const std::complex<double>* A, const int& lda, double* S, double* scond, double* amax, int* info) const
  {
    ZPOEQU_F77(&n, A, &lda, S, scond, amax, info);
  }


  void LAPACK<int, std::complex<double> >::PORFS(const char& UPLO, const int& n, const int& nrhs, const std::complex<double>* A, const int& lda, const std::complex<double>* AF, const int& ldaf, const std::complex<double>* B, const int& ldb, std::complex<double>* X, const int& ldx, double* FERR, double* BERR, std::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZPORFS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, AF, &ldaf, B, &ldb, X, &ldx, FERR, BERR, WORK, RWORK, info);
  }

  void LAPACK<int, std::complex<double> >::POSVX(const char& FACT, const char& UPLO, const int& n, const int& nrhs, std::complex<double>* A, const int& lda, std::complex<double>* AF, const int& ldaf, char* EQUED, double* S, std::complex<double>* B, const int& ldb, std::complex<double>* X, const int& ldx, double* rcond, double* FERR, double* BERR, std::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZPOSVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, AF, &ldaf, CHARPTR_MACRO(EQUED), S, B, &ldb, X, &ldx, rcond, FERR, BERR, WORK, RWORK, info);
  }


  void LAPACK<int,std::complex<double> >::GELS(const char& TRANS, const int& m, const int& n, const int& nrhs, std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, std::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZGELS_F77(CHAR_MACRO(TRANS), &m, &n, &nrhs, A, &lda, B, &ldb, WORK, &lwork, info);
  }


  void LAPACK<int, std::complex<double> >::GELSS(const int& m, const int& n, const int& nrhs, std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, double* S, const double& rcond, int* rank, std::complex<double>* WORK, const int& lwork, double* rwork, int* info) const
  {
    ZGELSS_F77(&m, &n, &nrhs, A, &lda, B, &ldb, S, &rcond, rank, WORK, &lwork, rwork, info);
  }


  void LAPACK<int,std::complex<double> >::GEQRF( const int& m, const int& n, std::complex<double>* A, const int& lda, std::complex<double>* TAU, std::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZGEQRF_F77(&m, &n, A, &lda, TAU, WORK, &lwork, info);
  }

  void LAPACK<int,std::complex<double> >::GEQR2 (const int& m, const int& n, std::complex<double>* A, const int& lda, std::complex<double>* TAU, std::complex<double>* WORK, int* const info) const
  {
    ZGEQR2_F77(&m, &n, A, &lda, TAU, WORK, info);
  }

  void LAPACK<int,std::complex<double> >::UNGQR(const int& m, const int& n, const int& k, std::complex<double>* A, const int& lda, const std::complex<double>* TAU, std::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZUNGQR_F77( &m, &n, &k, A, &lda, TAU, WORK, &lwork, info);
  }


  void LAPACK<int,std::complex<double> >::UNMQR(const char& SIDE, const char& TRANS, const int& m, const int& n, const int& k, const std::complex<double>* A, const int& lda, const std::complex<double>* TAU, std::complex<double>* C, const int& ldc, std::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZUNMQR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &m, &n, &k, A, &lda, TAU, C, &ldc, WORK, &lwork, info);
  }

  void LAPACK<int,std::complex<double> >::UNM2R (const char& SIDE, const char& TRANS, const int& M, const int& N, const int& K, const std::complex<double>* A, const int& LDA, const std::complex<double>* TAU, std::complex<double>* C, const int& LDC, std::complex<double>* WORK, int* const INFO) const
  {
    ZUNM2R_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, INFO);
  }

  void LAPACK<int,std::complex<double> >::GETRF(const int& m, const int& n, std::complex<double>* A, const int& lda, int* IPIV, int* info) const
  {
    ZGETRF_F77(&m, &n, A, &lda, IPIV, info);
  }


  void LAPACK<int,std::complex<double> >::GETRS(const char& TRANS, const int& n, const int& nrhs, const std::complex<double>* A, const int& lda, const int* IPIV, std::complex<double>* B, const int& ldb, int* info) const
  {
    ZGETRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, IPIV, B, &ldb, info);
  }


  void LAPACK<int,std::complex<double> >::LASCL(const char& TYPE, const int& kl, const int& ku, const double& cfrom, const double& cto, const int& m, const int& n, std::complex<double>* A, const int& lda, int* info) const
  { ZLASCL_F77(CHAR_MACRO(TYPE), &kl, &ku, &cfrom, &cto, &m, &n, A, &lda, info); }

  void LAPACK<int,std::complex<double> >::GEQP3(const int& m, const int& n, std::complex<double>* A, const int& lda, int* jpvt, std::complex<double>* TAU, std::complex<double>* WORK, const int& lwork, double* RWORK, int* info ) const
  {
    ZGEQP3_F77(&m, &n, A, &lda, jpvt, TAU, WORK, &lwork, RWORK, info);
  }

  void LAPACK<int, std::complex<double> >::LASWP (const int& N, std::complex<double>* A, const int& LDA, const int& K1, const int& K2, const int* IPIV, const int& INCX) const
  {
    ZLASWP_F77(&N, A, &LDA, &K1, &K2, IPIV, &INCX);
  }

  void LAPACK<int,std::complex<double> >::GBTRF(const int& m, const int& n, const int& kl, const int& ku, std::complex<double>* A, const int& lda, int* IPIV, int* info) const
  {
    ZGBTRF_F77(&m, &n, &kl, &ku, A, &lda, IPIV, info);
  }


  void LAPACK<int,std::complex<double> >::GBTRS(const char& TRANS, const int& n, const int& kl, const int& ku, const int& nrhs, const std::complex<double>* A, const int& lda, const int* IPIV, std::complex<double>* B, const int& ldb, int* info) const
  {
    ZGBTRS_F77(CHAR_MACRO(TRANS), &n, &kl, &ku, &nrhs, A, &lda, IPIV, B, &ldb, info);
  }


  void LAPACK<int,std::complex<double> >::GTTRF(const int& n, std::complex<double>* dl, std::complex<double>* d, std::complex<double>* du, std::complex<double>* du2, int* IPIV, int* info) const
  {
    ZGTTRF_F77(&n, dl, d, du, du2, IPIV, info);
  }


  void LAPACK<int,std::complex<double> >::GTTRS(const char& TRANS, const int& n, const int& nrhs, const std::complex<double>* dl, const std::complex<double>* d, const std::complex<double>* du, const std::complex<double>* du2, const int* IPIV, std::complex<double>* B, const int& ldb, int* info) const
  {
    ZGTTRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, dl, d, du, du2, IPIV, B, &ldb, info);
  }


  void LAPACK<int,std::complex<double> >::GETRI(const int& n, std::complex<double>* A, const int& lda, const int* IPIV, std::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZGETRI_F77(&n, A, &lda, IPIV, WORK, &lwork, info);
  }

  void LAPACK<int, std::complex<double> >::LATRS (const char& UPLO, const char& TRANS, const char& DIAG, const char& NORMIN, const int& N, const std::complex<double>* A, const int& LDA, std::complex<double>* X, double* SCALE, double* CNORM, int* INFO) const
  {
    ZLATRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), CHAR_MACRO(NORMIN), &N, A, &LDA, X, SCALE, CNORM, INFO);
  }

  void LAPACK<int,std::complex<double> >::GECON(const char& NORM, const int& n, const std::complex<double>* A, const int& lda, const double& anorm, double* rcond, std::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZGECON_F77(CHAR_MACRO(NORM), &n, A, &lda, &anorm, rcond, WORK, RWORK, info);
  }


  void LAPACK<int,std::complex<double> >::GBCON(const char& NORM, const int& n, const int& kl, const int& ku, const std::complex<double>* A, const int& lda, const int* IPIV, const double& anorm, double* rcond, std::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZGBCON_F77(CHAR_MACRO(NORM), &n, &kl, &ku, A, &lda, IPIV, &anorm, rcond, WORK, RWORK, info);
  }


  double LAPACK<int,std::complex<double> >::LANGB(const char& NORM, const int& n, const int& kl, const int& ku, const std::complex<double>* A, const int& lda, double* WORK) const
  {
    return( ZLANGB_F77(CHAR_MACRO(NORM), &n, &kl, &ku, A, &lda, WORK) );
  }


  void LAPACK<int,std::complex<double> >::GESV(const int& n, const int& nrhs, std::complex<double>* A, const int& lda, int* IPIV, std::complex<double>* B, const int& ldb, int* info) const
  {
    ZGESV_F77(&n, &nrhs, A, &lda, IPIV, B, &ldb, info);
  }


  void LAPACK<int,std::complex<double> >::GEEQU(const int& m, const int& n, const std::complex<double>* A, const int& lda, double* R, double* C, double* rowcond, double* colcond, double* amax, int* info) const
  {
    ZGEEQU_F77(&m, &n, A, &lda, R, C, rowcond, colcond, amax, info);
  }


  void LAPACK<int,std::complex<double> >::GERFS(const char& TRANS, const int& n, const int& nrhs, const std::complex<double>* A, const int& lda, const std::complex<double>* AF, const int& ldaf, const int* IPIV, const std::complex<double>* B, const int& ldb, std::complex<double>* X, const int& ldx, double* FERR, double* BERR, std::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZGERFS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, AF, &ldaf, IPIV, B, &ldb, X, &ldx, FERR, BERR, WORK, RWORK, info);
  }


  void LAPACK<int,std::complex<double> >::GBEQU(const int& m, const int& n, const int& kl, const int& ku, const std::complex<double>* A, const int& lda, double* R, double* C, double* rowcond, double* colcond, double* amax, int* info) const
  {
    ZGBEQU_F77(&m, &n, &kl, &ku, A, &lda, R, C, rowcond, colcond, amax, info);
  }


  void LAPACK<int,std::complex<double> >::GBRFS(const char& TRANS, const int& n, const int& kl, const int& ku, const int& nrhs, const std::complex<double>* A, const int& lda, const std::complex<double>* AF, const int& ldaf, const int* IPIV, const std::complex<double>* B, const int& ldb, std::complex<double>* X, const int& ldx, double* FERR, double* BERR, std::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZGBRFS_F77(CHAR_MACRO(TRANS), &n, &kl, &ku, &nrhs, A, &lda, AF, &ldaf, IPIV, B, &ldb, X, &ldx, FERR, BERR, WORK, RWORK, info);
  }

  void LAPACK<int,std::complex<double> >::GESVX(const char& FACT, const char& TRANS, const int& n, const int& nrhs, std::complex<double>* A, const int& lda, std::complex<double>* AF, const int& ldaf, int* IPIV, char* EQUED, double* R, double* C, std::complex<double>* B, const int& ldb, std::complex<double>* X, const int& ldx, double* rcond, double* FERR, double* BERR, std::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZGESVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, AF, &ldaf, IPIV, CHARPTR_MACRO(EQUED), R, C, B, &ldb, X, &ldx, rcond, FERR, BERR, WORK, RWORK, info);
  }


  void LAPACK<int,std::complex<double> >::GEHRD(const int& n, const int& ilo, const int& ihi, std::complex<double>* A, const int& lda, std::complex<double>* TAU, std::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZGEHRD_F77(&n, &ilo, &ihi, A, &lda, TAU, WORK, &lwork, info);
  }


  void LAPACK<int,std::complex<double> >::TRTRS(const char& UPLO, const char& TRANS, const char& DIAG, const int& n, const int& nrhs, const std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, int* info) const
  {
    ZTRTRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), &n, &nrhs, A, &lda, B, &ldb, info);
  }


  void LAPACK<int,std::complex<double> >::TRTRI(const char& UPLO, const char& DIAG, const int& n, std::complex<double>* A, const int& lda, int* info) const
  {
    ZTRTRI_F77(CHAR_MACRO(UPLO), CHAR_MACRO(DIAG), &n, A, &lda, info);
  }


  void LAPACK<int,std::complex<double> >::STEQR(const char& COMPZ, const int& n, double* D, double* E, std::complex<double>* Z, const int& ldz, double* WORK, int* info) const
  {
    ZSTEQR_F77(CHAR_MACRO(COMPZ), &n, D, E, Z, &ldz, WORK, info);
  }


  void LAPACK<int,std::complex<double> >::PTEQR(const char& COMPZ, const int& n, double* D, double* E, std::complex<double>* Z, const int& ldz, double* WORK, int* info) const
  {
    ZPTEQR_F77(CHAR_MACRO(COMPZ), &n, D, E, Z, &ldz, WORK, info);
  }


  void LAPACK<int,std::complex<double> >::HEEV(const char& JOBZ, const char& UPLO, const int& n, std::complex<double> * A, const int& lda, double*  W, std::complex<double> * WORK, const int& lwork, double* RWORK, int* info) const
  {
    ZHEEV_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, A, &lda, W, WORK, &lwork, RWORK, info);
  }


  void LAPACK<int,std::complex<double> >::HEGV(const int& itype, const char& JOBZ, const char& UPLO, const int& n, std::complex<double> * A, const int& lda, std::complex<double> * B, const int& ldb, double*  W, std::complex<double> * WORK, const int& lwork, double* RWORK, int* info) const
  {
    ZHEGV_F77(&itype, CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, A, &lda, B, &ldb, W, WORK, &lwork, RWORK, info);
  }


  void LAPACK<int, std::complex<double> >::HSEQR(const char& JOB, const char& COMPZ, const int& n, const int& ilo, const int& ihi, std::complex<double>* H, const int& ldh, std::complex<double>* W, std::complex<double>* Z, const int& ldz, std::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZHSEQR_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPZ), &n, &ilo, &ihi, H, &ldh, W, Z, &ldz, WORK, &lwork, info);
  }


  void LAPACK<int, std::complex<double> >::GEES(const char& JOBVS, const char& SORT, int (*ptr2func)(std::complex<double>*), const int& n, std::complex<double>* A, const int& lda, int* sdim, std::complex<double>* W, std::complex<double>* VS, const int& ldvs, std::complex<double>* WORK, const int& lwork, double* RWORK, int* BWORK, int* info) const
  {
    ZGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(SORT), ptr2func, &n, A, &lda, sdim, W, VS, &ldvs, WORK, &lwork, RWORK, BWORK, info);
  }


  void LAPACK<int, std::complex<double> >::GEES(const char& JOBVS, const int& n, std::complex<double>* A, const int& lda, int* sdim, double* WR, double* WI, std::complex<double>* VS, const int& ldvs, std::complex<double>* WORK, const int& lwork, double* RWORK, int* BWORK, int* info) const
  {
    int (*nullfptr)(std::complex<double>*) = NULL;
    std::vector< std::complex<double> > W(n);
    const char sort = 'N';
    ZGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(sort), nullfptr, &n, A, &lda, sdim, &W[0], VS, &ldvs, WORK, &lwork, RWORK, BWORK, info);
    for (int i=0; i<n; i++) {
      WR[i] = W[i].real();
      WI[i] = W[i].imag();
    }
  }


  void LAPACK<int, std::complex<double> >::GEEV(const char& JOBVL, const char& JOBVR, const int& n, std::complex<double>* A, const int& lda, std::complex<double>* W, std::complex<double>* VL, const int& ldvl, std::complex<double>* VR, const int& ldvr, std::complex<double>* WORK, const int& lwork, double* RWORK, int* info) const
  {
    ZGEEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &n, A, &lda, W, VL, &ldvl, VR, &ldvr, WORK, &lwork, RWORK, info);
  }


  void LAPACK<int, std::complex<double> >::GEEV(const char& JOBVL, const char& JOBVR, const int& n, std::complex<double>* A, const int& lda, double* WR, double* WI, std::complex<double>* VL, const int& ldvl, std::complex<double>* VR, const int& ldvr, std::complex<double>* WORK, const int& lwork, double* RWORK, int* info) const
  {
    std::vector<std::complex<double> > w (n);
    std::complex<double>* w_rawPtr = (n == 0) ? NULL : &w[0];
    GEEV (JOBVL, JOBVR, n, A, lda, w_rawPtr, VL, ldvl, VR, ldvr, WORK, lwork, RWORK, info);
    if (*info == 0) {
      // The eigenvalues are only valid on output if INFO is zero.
      // Otherwise, we shouldn't even write to WR or WI.
      for (int k = 0; k < n; ++k) {
        WR[k] = w[k].real ();
        WI[k] = w[k].imag ();
      }
    }
  }


  void LAPACK<int, std::complex<double> >::GESVD(const char& JOBU, const char& JOBVT, const int& m, const int& n, std::complex<double> * A, const int& lda, double* S, std::complex<double> * U, const int& ldu, std::complex<double> * V, const int& ldv, std::complex<double> * WORK, const int& lwork, double* RWORK, int* info) const {
    ZGESVD_F77(CHAR_MACRO(JOBU), CHAR_MACRO(JOBVT), &m, &n, A, &lda, S, U, &ldu, V, &ldv, WORK, &lwork, RWORK, info);
  }

  void LAPACK<int, std::complex<double> >::GEEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, std::complex<double>* A, const int& lda, std::complex<double>* W, std::complex<double>* VL, const int& ldvl, std::complex<double>* VR, const int& ldvr, int* ilo, int* ihi, double* SCALE, double* abnrm, double* RCONDE, double* RCONDV, std::complex<double>* WORK, const int& lwork, double* RWORK, int* info) const
  {
    ZGEEVX_F77(CHAR_MACRO(BALANC), CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SENSE), &n, A, &lda, W, VL, &ldvl, VR, &ldvr, ilo, ihi, SCALE, abnrm, RCONDE, RCONDV, WORK, &lwork, RWORK, info);
  }

  void LAPACK<int, std::complex<double> >::GGEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, std::complex<double>* ALPHA, std::complex<double>* BETA, std::complex<double>* VL, const int& ldvl, std::complex<double>* VR, const int& ldvr, int* ilo, int* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, double* RCONDE, double* RCONDV, std::complex<double>* WORK, const int& lwork, double* RWORK, int* IWORK, int* BWORK, int* info) const
  {
    ZGGEVX_F77(CHAR_MACRO(BALANC), CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SENSE), &n, A, &lda, B, &ldb, ALPHA, BETA, VL, &ldvl, VR, &ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, RCONDE, RCONDV, WORK, &lwork, RWORK, IWORK, BWORK, info);
  }

  void LAPACK<int, std::complex<double> >::GGEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, double* ALPHAR, double* ALPHAI, std::complex<double>* BETA, std::complex<double>* VL, const int& ldvl, std::complex<double>* VR, const int& ldvr, int* ilo, int* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, double* RCONDE, double* RCONDV, std::complex<double>* WORK, const int& lwork, double* RWORK, int* IWORK, int* BWORK, int* info) const
  {
    std::vector<std::complex<double> > w (n);
    std::complex<double>* w_rawPtr = (n == 0) ? NULL : &w[0];
    GGEVX(BALANC, JOBVL, JOBVR, SENSE, n, A, lda, B, ldb, w_rawPtr, BETA, VL, ldvl, VR, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, RCONDE, RCONDV, WORK, lwork, RWORK, IWORK, BWORK, info);
    if (*info == 0) {
      // The eigenvalues are only valid on output if INFO is zero.
      // Otherwise, we shouldn't even write to WR or WI.
      for (int k = 0; k < n; ++k) {
        ALPHAR[k] = w[k].real ();
        ALPHAI[k] = w[k].imag ();
      }
    }
  }

 void LAPACK<int, std::complex<double> >::GGEV(const char& JOBVL, const char& JOBVR, const int& n, std::complex<double> *A, const int& lda, std::complex<double> *B, const int& ldb, std::complex<double>* ALPHA, std::complex<double>* BETA, std::complex<double>* VL, const int& ldvl, std::complex<double>* VR, const int& ldvr, std::complex<double> *WORK, const int& lwork, double* RWORK, int* info) const
  {
    ZGGEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &n, A, &lda, B, &ldb, ALPHA, BETA, VL, &ldvl, VR, &ldvr, WORK, &lwork, RWORK, info); 
  }

  void LAPACK<int, std::complex<double> >::GGES(const char& JOBVL, const char& JOBVR, const char& SORT, int (*ptr2func)(std::complex<double>*, std::complex<double>*), const int& n, std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, int* sdim, std::complex<double>* ALPHA, std::complex<double>* BETA, std::complex<double>* VL, const int& ldvl, std::complex<double>* VR, const int& ldvr, std::complex<double>* WORK, const int& lwork, double* rwork, int* bwork, int* info ) const
  {
    ZGGES_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SORT), ptr2func, &n, A, &lda, B, &ldb, sdim, ALPHA, BETA, VL, &ldvl, VR, &ldvr, WORK, &lwork, rwork, bwork, info); 
  }
  void LAPACK<int, std::complex<double> >::TGSEN(const int& ijob, const int& wantq, const int& wantz, const int* SELECT, const int& n, std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, std::complex<double>* ALPHA, std::complex<double>* BETA, std::complex<double>* Q, const int& ldq, std::complex<double>* Z, const int& ldz, int* M, double* PL, double* PR, double* DIF, std::complex<double>* WORK, const int& lwork, int* IWORK, const int& liwork, int* info ) const
  {
    ZTGSEN_F77(&ijob, &wantq, &wantz, SELECT, &n, A, &lda, B, &ldb, ALPHA, BETA, Q, &ldq, Z, &ldz, M, PL, PR, DIF, WORK, &lwork, IWORK, &liwork, info); 
  }

  void LAPACK<int, std::complex<double> >::TREVC(const char& SIDE, const char& HOWMNY, int* select, const int& n, const std::complex<double>* T, const int& ldt, std::complex<double>* VL, const int& ldvl, std::complex<double>* VR, const int& ldvr, const int& mm, int* m, std::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZTREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), select, &n, T, &ldt, VL, &ldvl, VR, &ldvr, &mm, m, WORK, RWORK, info);
  }


  void LAPACK<int, std::complex<double> >::TREVC(const char& SIDE, const int& n, const std::complex<double>* T, const int& ldt, std::complex<double>* VL, const int& ldvl, std::complex<double>* VR, const int& ldvr, const int& mm, int* m, std::complex<double>* WORK, double* RWORK, int* info) const
  {
    std::vector<int> select(1);
    const char& whch = 'A';
    ZTREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(whch), &select[0], &n, T, &ldt, VL, &ldvl, VR, &ldvr, &mm, m, WORK, RWORK, info);
  }
  
  void LAPACK<int, std::complex<double> >::TREXC(const char& COMPQ, const int& n, std::complex<double>* T, const int& ldt, std::complex<double>* Q, const int& ldq, int* ifst, int* ilst, std::complex<double>* WORK, int* info) const
  {
    ZTREXC_F77(CHAR_MACRO(COMPQ), &n, T, &ldt, Q, &ldq, ifst, ilst, info);
  }

  void LAPACK<int, std::complex<double> >::LARTG( const std::complex<double> f, const std::complex<double> g, double* c, std::complex<double>* s, std::complex<double>* r ) const
  {
    ZLARTG_F77(&f, &g, c, s, r);
  }


  void LAPACK<int, std::complex<double> >::LARFG( const int& n, std::complex<double>* alpha, std::complex<double>* x, const int& incx, std::complex<double>* tau ) const
  {
    ZLARFG_F77(&n, alpha, x, &incx, tau);
  }

  void LAPACK<int, std::complex<double> >::GEBAL(const char& JOBZ, const int& n, std::complex<double>* A, const int& lda, int* ilo, int* ihi, double* scale, int* info) const
  {
    ZGEBAL_F77(CHAR_MACRO(JOBZ),&n, A, &lda, ilo, ihi, scale, info);
  }


  void LAPACK<int, std::complex<double> >::GEBAK(const char& JOBZ, const char& SIDE, const int& n, const int& ilo, const int& ihi, const double* scale, const int& m, std::complex<double>* V, const int& ldv, int* info) const
  {
    ZGEBAK_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(SIDE), &n, &ilo, &ihi, scale, &m, V, &ldv, info);
  }


#ifdef HAVE_TEUCHOS_LAPACKLARND
  std::complex<double> LAPACK<int, std::complex<double> >::LARND( const int& idist, int* seed ) const
  {
    double _Complex z = ZLARND_F77(&idist, seed);
    return TEUCHOS_LAPACK_CONVERT_COMPLEX_FORTRAN_TO_CXX(double, z);
  }
#endif

  void LAPACK<int, std::complex<double> >::LARNV( const int& idist, int* seed, const int& n, std::complex<double>* v ) const
  {
    ZLARNV_F77(&idist, seed, &n, v);
  }


  int LAPACK<int, std::complex<double> >::ILAENV( const int& ispec, const std::string& NAME, const std::string& OPTS, const int& N1, const int& N2, const int& N3, const int& N4 ) const
  {
    unsigned int opts_length = OPTS.length();
    std::string temp_NAME = "z" + NAME;
    unsigned int name_length = temp_NAME.length();
    return ilaenv_wrapper(&ispec, &temp_NAME[0], name_length, &OPTS[0], opts_length, &N1, &N2, &N3, &N4);
  }

  // END INT, COMPLEX<DOUBLE> SPECIALIZATION IMPLEMENTATION //

  // BEGIN INT, KOKKOS::COMPLEX<DOUBLE> SPECIALIZATION IMPLEMENTATION //

  void LAPACK<int, Kokkos::complex<double> >::PTTRF(const int& n, double* d, Kokkos::complex<double>* e, int* info) const
  {
    ZPTTRF_F77(&n, d, reinterpret_cast<std::complex<double>*>(e), info);
  }


  void LAPACK<int, Kokkos::complex<double> >::PTTRS(const char& UPLO, const int& n, const int& nrhs, const double* d, const Kokkos::complex<double>* e, Kokkos::complex<double>* B, const int& ldb, int* info) const
  {
    ZPTTRS_F77(CHAR_MACRO(UPLO), &n, &nrhs, d, reinterpret_cast<const std::complex<double>*>(e), reinterpret_cast<std::complex<double>*>(B), &ldb, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::POTRF(const char& UPLO, const int& n, Kokkos::complex<double>* A, const int& lda, int* info) const
  {
    ZPOTRF_F77(CHAR_MACRO(UPLO), &n, reinterpret_cast<std::complex<double>*>(A), &lda, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::POTRS(const char& UPLO, const int& n, const int& nrhs, const Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* B, const int& ldb, int* info) const
  {
    ZPOTRS_F77(CHAR_MACRO(UPLO), &n, &nrhs, reinterpret_cast<const std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(B), &ldb, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::POTRI(const char& UPLO, const int& n, Kokkos::complex<double>* A, const int& lda, int* info) const
  {
    ZPOTRI_F77(CHAR_MACRO(UPLO), &n, reinterpret_cast<std::complex<double>*>(A), &lda, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::POCON(const char& UPLO, const int& n, const Kokkos::complex<double>* A, const int& lda, const double& anorm, double* rcond, Kokkos::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZPOCON_F77(CHAR_MACRO(UPLO), &n, reinterpret_cast<const std::complex<double>*>(A), &lda, &anorm, rcond, reinterpret_cast<std::complex<double>*>(WORK), RWORK, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::POSV(const char& UPLO, const int& n, const int& nrhs, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* B, const int& ldb, int* info) const
  {
    ZPOSV_F77(CHAR_MACRO(UPLO), &n, &nrhs, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(B), &ldb, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::POEQU(const int& n, const Kokkos::complex<double>* A, const int& lda, double* S, double* scond, double* amax, int* info) const
  {
    ZPOEQU_F77(&n, reinterpret_cast<const std::complex<double>*>(A), &lda, S, scond, amax, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::PORFS(const char& UPLO, const int& n, const int& nrhs, const Kokkos::complex<double>* A, const int& lda, const Kokkos::complex<double>* AF, const int& ldaf, const Kokkos::complex<double>* B, const int& ldb, Kokkos::complex<double>* X, const int& ldx, double* FERR, double* BERR, Kokkos::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZPORFS_F77(CHAR_MACRO(UPLO), &n, &nrhs, reinterpret_cast<const std::complex<double>*>(A), &lda, reinterpret_cast<const std::complex<double>*>(AF), &ldaf, reinterpret_cast<const std::complex<double>*>(B), &ldb, reinterpret_cast<std::complex<double>*>(X), &ldx, FERR, BERR, reinterpret_cast<std::complex<double>*>(WORK), RWORK, info);
  }

  void LAPACK<int, Kokkos::complex<double> >::POSVX(const char& FACT, const char& UPLO, const int& n, const int& nrhs, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* AF, const int& ldaf, char* EQUED, double* S, Kokkos::complex<double>* B, const int& ldb, Kokkos::complex<double>* X, const int& ldx, double* rcond, double* FERR, double* BERR, Kokkos::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZPOSVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(UPLO), &n, &nrhs, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(AF), &ldaf, CHARPTR_MACRO(EQUED), S, reinterpret_cast<std::complex<double>*>(B), &ldb, reinterpret_cast<std::complex<double>*>(X), &ldx, rcond, FERR, BERR, reinterpret_cast<std::complex<double>*>(WORK), RWORK, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::GELS(const char& TRANS, const int& m, const int& n, const int& nrhs, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* B, const int& ldb, Kokkos::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZGELS_F77(CHAR_MACRO(TRANS), &m, &n, &nrhs, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(B), &ldb, reinterpret_cast<std::complex<double>*>(WORK), &lwork, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::GELSS(const int& m, const int& n, const int& nrhs, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* B, const int& ldb, double* S, const double& rcond, int* rank, Kokkos::complex<double>* WORK, const int& lwork, double* rwork, int* info) const
  {
    ZGELSS_F77(&m, &n, &nrhs, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(B), &ldb, S, &rcond, rank, reinterpret_cast<std::complex<double>*>(WORK), &lwork, rwork, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::GEQRF( const int& m, const int& n, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* TAU, Kokkos::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZGEQRF_F77(&m, &n, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(TAU), reinterpret_cast<std::complex<double>*>(WORK), &lwork, info);
  }

  void LAPACK<int,Kokkos::complex<double> >::GEQR2 (const int& m, const int& n, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* TAU, Kokkos::complex<double>* WORK, int* const info) const
  {
    ZGEQR2_F77(&m, &n, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(TAU), reinterpret_cast<std::complex<double>*>(WORK), info);
  }

  void LAPACK<int,Kokkos::complex<double> >::UNGQR(const int& m, const int& n, const int& k, Kokkos::complex<double>* A, const int& lda, const Kokkos::complex<double>* TAU, Kokkos::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZUNGQR_F77( &m, &n, &k, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<const std::complex<double>*>(TAU), reinterpret_cast<std::complex<double>*>(WORK), &lwork, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::UNMQR(const char& SIDE, const char& TRANS, const int& m, const int& n, const int& k, const Kokkos::complex<double>* A, const int& lda, const Kokkos::complex<double>* TAU, Kokkos::complex<double>* C, const int& ldc, Kokkos::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZUNMQR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &m, &n, &k, reinterpret_cast<const std::complex<double>*>(A), &lda, reinterpret_cast<const std::complex<double>*>(TAU), reinterpret_cast<std::complex<double>*>(C), &ldc, reinterpret_cast<std::complex<double>*>(WORK), &lwork, info);
  }

  void LAPACK<int,Kokkos::complex<double> >::UNM2R (const char& SIDE, const char& TRANS, const int& M, const int& N, const int& K, const Kokkos::complex<double>* A, const int& LDA, const Kokkos::complex<double>* TAU, Kokkos::complex<double>* C, const int& LDC, Kokkos::complex<double>* WORK, int* const INFO) const
  {
    ZUNM2R_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &M, &N, &K, reinterpret_cast<const std::complex<double>*>(A), &LDA, reinterpret_cast<const std::complex<double>*>(TAU), reinterpret_cast<std::complex<double>*>(C), &LDC, reinterpret_cast<std::complex<double>*>(WORK), INFO);
  }

  void LAPACK<int,Kokkos::complex<double> >::GETRF(const int& m, const int& n, Kokkos::complex<double>* A, const int& lda, int* IPIV, int* info) const
  {
    ZGETRF_F77(&m, &n, reinterpret_cast<std::complex<double>*>(A), &lda, IPIV, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::GETRS(const char& TRANS, const int& n, const int& nrhs, const Kokkos::complex<double>* A, const int& lda, const int* IPIV, Kokkos::complex<double>* B, const int& ldb, int* info) const
  {
    ZGETRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, reinterpret_cast<const std::complex<double>*>(A), &lda, IPIV, reinterpret_cast<std::complex<double>*>(B), &ldb, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::LASCL(const char& TYPE, const int& kl, const int& ku, const double& cfrom, const double& cto, const int& m, const int& n, Kokkos::complex<double>* A, const int& lda, int* info) const
  { ZLASCL_F77(CHAR_MACRO(TYPE), &kl, &ku, &cfrom, &cto, &m, &n, reinterpret_cast<std::complex<double>*>(A), &lda, info); }

  void LAPACK<int,Kokkos::complex<double> >::GEQP3(const int& m, const int& n, Kokkos::complex<double>* A, const int& lda, int* jpvt, Kokkos::complex<double>* TAU, Kokkos::complex<double>* WORK, const int& lwork, double* RWORK, int* info ) const
  {
    ZGEQP3_F77(&m, &n, reinterpret_cast<std::complex<double>*>(A), &lda, jpvt, reinterpret_cast<std::complex<double>*>(TAU), reinterpret_cast<std::complex<double>*>(WORK), &lwork, RWORK, info);
  }

  void LAPACK<int, Kokkos::complex<double> >::LASWP (const int& N, Kokkos::complex<double>* A, const int& LDA, const int& K1, const int& K2, const int* IPIV, const int& INCX) const
  {
    ZLASWP_F77(&N, reinterpret_cast<std::complex<double>*>(A), &LDA, &K1, &K2, IPIV, &INCX);
  }

  void LAPACK<int,Kokkos::complex<double> >::GBTRF(const int& m, const int& n, const int& kl, const int& ku, Kokkos::complex<double>* A, const int& lda, int* IPIV, int* info) const
  {
    ZGBTRF_F77(&m, &n, &kl, &ku, reinterpret_cast<std::complex<double>*>(A), &lda, IPIV, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::GBTRS(const char& TRANS, const int& n, const int& kl, const int& ku, const int& nrhs, const Kokkos::complex<double>* A, const int& lda, const int* IPIV, Kokkos::complex<double>* B, const int& ldb, int* info) const
  {
    ZGBTRS_F77(CHAR_MACRO(TRANS), &n, &kl, &ku, &nrhs, reinterpret_cast<const std::complex<double>*>(A), &lda, IPIV, reinterpret_cast<std::complex<double>*>(B), &ldb, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::GTTRF(const int& n, Kokkos::complex<double>* dl, Kokkos::complex<double>* d, Kokkos::complex<double>* du, Kokkos::complex<double>* du2, int* IPIV, int* info) const
  {
    ZGTTRF_F77(&n, reinterpret_cast<std::complex<double>*>(dl), reinterpret_cast<std::complex<double>*>(d), reinterpret_cast<std::complex<double>*>(du), reinterpret_cast<std::complex<double>*>(du2), IPIV, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::GTTRS(const char& TRANS, const int& n, const int& nrhs, const Kokkos::complex<double>* dl, const Kokkos::complex<double>* d, const Kokkos::complex<double>* du, const Kokkos::complex<double>* du2, const int* IPIV, Kokkos::complex<double>* B, const int& ldb, int* info) const
  {
    ZGTTRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, reinterpret_cast<const std::complex<double>*>(dl), reinterpret_cast<const std::complex<double>*>(d), reinterpret_cast<const std::complex<double>*>(du), reinterpret_cast<const std::complex<double>*>(du2), IPIV, reinterpret_cast<std::complex<double>*>(B), &ldb, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::GETRI(const int& n, Kokkos::complex<double>* A, const int& lda, const int* IPIV, Kokkos::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZGETRI_F77(&n, reinterpret_cast<std::complex<double>*>(A), &lda, IPIV, reinterpret_cast<std::complex<double>*>(WORK), &lwork, info);
  }

  void LAPACK<int, Kokkos::complex<double> >::LATRS (const char& UPLO, const char& TRANS, const char& DIAG, const char& NORMIN, const int& N, const Kokkos::complex<double>* A, const int& LDA, Kokkos::complex<double>* X, double* SCALE, double* CNORM, int* INFO) const
  {
    ZLATRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), CHAR_MACRO(NORMIN), &N, reinterpret_cast<const std::complex<double>*>(A), &LDA, reinterpret_cast<std::complex<double>*>(X), SCALE, CNORM, INFO);
  }

  void LAPACK<int,Kokkos::complex<double> >::GECON(const char& NORM, const int& n, const Kokkos::complex<double>* A, const int& lda, const double& anorm, double* rcond, Kokkos::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZGECON_F77(CHAR_MACRO(NORM), &n, reinterpret_cast<const std::complex<double>*>(A), &lda, &anorm, rcond, reinterpret_cast<std::complex<double>*>(WORK), RWORK, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::GBCON(const char& NORM, const int& n, const int& kl, const int& ku, const Kokkos::complex<double>* A, const int& lda, const int* IPIV, const double& anorm, double* rcond, Kokkos::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZGBCON_F77(CHAR_MACRO(NORM), &n, &kl, &ku, reinterpret_cast<const std::complex<double>*>(A), &lda, IPIV, &anorm, rcond, reinterpret_cast<std::complex<double>*>(WORK), RWORK, info);
  }


  double LAPACK<int,Kokkos::complex<double> >::LANGB(const char& NORM, const int& n, const int& kl, const int& ku, const Kokkos::complex<double>* A, const int& lda, double* WORK) const
  {
    return( ZLANGB_F77(CHAR_MACRO(NORM), &n, &kl, &ku, reinterpret_cast<const std::complex<double>*>(A), &lda, WORK) );
  }


  void LAPACK<int,Kokkos::complex<double> >::GESV(const int& n, const int& nrhs, Kokkos::complex<double>* A, const int& lda, int* IPIV, Kokkos::complex<double>* B, const int& ldb, int* info) const
  {
    ZGESV_F77(&n, &nrhs, reinterpret_cast<std::complex<double>*>(A), &lda, IPIV, reinterpret_cast<std::complex<double>*>(B), &ldb, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::GEEQU(const int& m, const int& n, const Kokkos::complex<double>* A, const int& lda, double* R, double* C, double* rowcond, double* colcond, double* amax, int* info) const
  {
    ZGEEQU_F77(&m, &n, reinterpret_cast<const std::complex<double>*>(A), &lda, R, C, rowcond, colcond, amax, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::GERFS(const char& TRANS, const int& n, const int& nrhs, const Kokkos::complex<double>* A, const int& lda, const Kokkos::complex<double>* AF, const int& ldaf, const int* IPIV, const Kokkos::complex<double>* B, const int& ldb, Kokkos::complex<double>* X, const int& ldx, double* FERR, double* BERR, Kokkos::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZGERFS_F77(CHAR_MACRO(TRANS), &n, &nrhs, reinterpret_cast<const std::complex<double>*>(A), &lda, reinterpret_cast<const std::complex<double>*>(AF), &ldaf, IPIV, reinterpret_cast<const std::complex<double>*>(B), &ldb, reinterpret_cast<std::complex<double>*>(X), &ldx, FERR, BERR, reinterpret_cast<std::complex<double>*>(WORK), RWORK, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::GBEQU(const int& m, const int& n, const int& kl, const int& ku, const Kokkos::complex<double>* A, const int& lda, double* R, double* C, double* rowcond, double* colcond, double* amax, int* info) const
  {
    ZGBEQU_F77(&m, &n, &kl, &ku, reinterpret_cast<const std::complex<double>*>(A), &lda, R, C, rowcond, colcond, amax, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::GBRFS(const char& TRANS, const int& n, const int& kl, const int& ku, const int& nrhs, const Kokkos::complex<double>* A, const int& lda, const Kokkos::complex<double>* AF, const int& ldaf, const int* IPIV, const Kokkos::complex<double>* B, const int& ldb, Kokkos::complex<double>* X, const int& ldx, double* FERR, double* BERR, Kokkos::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZGBRFS_F77(CHAR_MACRO(TRANS), &n, &kl, &ku, &nrhs, reinterpret_cast<const std::complex<double>*>(A), &lda, reinterpret_cast<const std::complex<double>*>(AF), &ldaf, IPIV, reinterpret_cast<const std::complex<double>*>(B), &ldb, reinterpret_cast<std::complex<double>*>(X), &ldx, FERR, BERR, reinterpret_cast<std::complex<double>*>(WORK), RWORK, info);
  }

  void LAPACK<int,Kokkos::complex<double> >::GESVX(const char& FACT, const char& TRANS, const int& n, const int& nrhs, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* AF, const int& ldaf, int* IPIV, char* EQUED, double* R, double* C, Kokkos::complex<double>* B, const int& ldb, Kokkos::complex<double>* X, const int& ldx, double* rcond, double* FERR, double* BERR, Kokkos::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZGESVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(TRANS), &n, &nrhs, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(AF), &ldaf, IPIV, CHARPTR_MACRO(EQUED), R, C, reinterpret_cast<std::complex<double>*>(B), &ldb, reinterpret_cast<std::complex<double>*>(X), &ldx, rcond, FERR, BERR, reinterpret_cast<std::complex<double>*>(WORK), RWORK, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::GEHRD(const int& n, const int& ilo, const int& ihi, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* TAU, Kokkos::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZGEHRD_F77(&n, &ilo, &ihi, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(TAU), reinterpret_cast<std::complex<double>*>(WORK), &lwork, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::TRTRS(const char& UPLO, const char& TRANS, const char& DIAG, const int& n, const int& nrhs, const Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* B, const int& ldb, int* info) const
  {
    ZTRTRS_F77(CHAR_MACRO(UPLO), CHAR_MACRO(TRANS), CHAR_MACRO(DIAG), &n, &nrhs, reinterpret_cast<const std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(B), &ldb, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::TRTRI(const char& UPLO, const char& DIAG, const int& n, Kokkos::complex<double>* A, const int& lda, int* info) const
  {
    ZTRTRI_F77(CHAR_MACRO(UPLO), CHAR_MACRO(DIAG), &n, reinterpret_cast<std::complex<double>*>(A), &lda, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::STEQR(const char& COMPZ, const int& n, double* D, double* E, Kokkos::complex<double>* Z, const int& ldz, double* WORK, int* info) const
  {
    ZSTEQR_F77(CHAR_MACRO(COMPZ), &n, D, E, reinterpret_cast<std::complex<double>*>(Z), &ldz, WORK, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::PTEQR(const char& COMPZ, const int& n, double* D, double* E, Kokkos::complex<double>* Z, const int& ldz, double* WORK, int* info) const
  {
    ZPTEQR_F77(CHAR_MACRO(COMPZ), &n, D, E, reinterpret_cast<std::complex<double>*>(Z), &ldz, WORK, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::HEEV(const char& JOBZ, const char& UPLO, const int& n, Kokkos::complex<double> * A, const int& lda, double* W, Kokkos::complex<double> * WORK, const int& lwork, double* RWORK, int* info) const
  {
    ZHEEV_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, reinterpret_cast<std::complex<double>*>(A), &lda, W, reinterpret_cast<std::complex<double>*>(WORK), &lwork, RWORK, info);
  }


  void LAPACK<int,Kokkos::complex<double> >::HEGV(const int& itype, const char& JOBZ, const char& UPLO, const int& n, Kokkos::complex<double> * A, const int& lda, Kokkos::complex<double> * B, const int& ldb, double* W, Kokkos::complex<double> * WORK, const int& lwork, double* RWORK, int* info) const
  {
    ZHEGV_F77(&itype, CHAR_MACRO(JOBZ), CHAR_MACRO(UPLO), &n, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(B), &ldb, W, reinterpret_cast<std::complex<double>*>(WORK), &lwork, RWORK, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::HSEQR(const char& JOB, const char& COMPZ, const int& n, const int& ilo, const int& ihi, Kokkos::complex<double>* H, const int& ldh, Kokkos::complex<double>* W, Kokkos::complex<double>* Z, const int& ldz, Kokkos::complex<double>* WORK, const int& lwork, int* info) const
  {
    ZHSEQR_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPZ), &n, &ilo, &ihi, reinterpret_cast<std::complex<double>*>(H), &ldh, reinterpret_cast<std::complex<double>*>(W), reinterpret_cast<std::complex<double>*>(Z), &ldz, reinterpret_cast<std::complex<double>*>(WORK), &lwork, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::GEES(const char& JOBVS, const char& SORT, int (*ptr2func)(std::complex<double>*), const int& n, Kokkos::complex<double>* A, const int& lda, int* sdim, Kokkos::complex<double>* W, Kokkos::complex<double>* VS, const int& ldvs, Kokkos::complex<double>* WORK, const int& lwork, double* RWORK, int* BWORK, int* info) const
  {
    ZGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(SORT), ptr2func, &n, reinterpret_cast<std::complex<double>*>(A), &lda, sdim, reinterpret_cast<std::complex<double>*>(W), reinterpret_cast<std::complex<double>*>(VS), &ldvs, reinterpret_cast<std::complex<double>*>(WORK), &lwork, RWORK, BWORK, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::GEES(const char& JOBVS, const int& n, Kokkos::complex<double>* A, const int& lda, int* sdim, double* WR, double* WI, Kokkos::complex<double>* VS, const int& ldvs, Kokkos::complex<double>* WORK, const int& lwork, double* RWORK, int* BWORK, int* info) const
  {
    int (*nullfptr)(std::complex<double>*) = NULL;
    std::vector< std::complex<double> > W(n);
    const char sort = 'N';
    ZGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(sort), nullfptr, &n, reinterpret_cast<std::complex<double>*>(A), &lda, sdim, &W[0], reinterpret_cast<std::complex<double>*>(VS), &ldvs, reinterpret_cast<std::complex<double>*>(WORK), &lwork, RWORK, BWORK, info);
    for (int i=0; i<n; i++) {
      WR[i] = W[i].real();
      WI[i] = W[i].imag();
    }
  }


  void LAPACK<int, Kokkos::complex<double> >::GEEV(const char& JOBVL, const char& JOBVR, const int& n, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* W, Kokkos::complex<double>* VL, const int& ldvl, Kokkos::complex<double>* VR, const int& ldvr, Kokkos::complex<double>* WORK, const int& lwork, double* RWORK, int* info) const
  {
    ZGEEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &n, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(W), reinterpret_cast<std::complex<double>*>(VL), &ldvl, reinterpret_cast<std::complex<double>*>(VR), &ldvr, reinterpret_cast<std::complex<double>*>(WORK), &lwork, RWORK, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::GEEV(const char& JOBVL, const char& JOBVR, const int& n, Kokkos::complex<double>* A, const int& lda, double* WR, double* WI, Kokkos::complex<double>* VL, const int& ldvl, Kokkos::complex<double>* VR, const int& ldvr, Kokkos::complex<double>* WORK, const int& lwork, double* RWORK, int* info) const
  {
    std::vector<Kokkos::complex<double> > w (n);
    Kokkos::complex<double>* w_rawPtr = (n == 0) ? NULL : &w[0];
    GEEV (JOBVL, JOBVR, n, A, lda, w_rawPtr, VL, ldvl, VR, ldvr, WORK, lwork, RWORK, info);
    if (*info == 0) {
      // The eigenvalues are only valid on output if INFO is zero.
      // Otherwise, we shouldn't even write to WR or WI.
      for (int k = 0; k < n; ++k) {
        WR[k] = w[k].real ();
        WI[k] = w[k].imag ();
      }
    }
  }


  void LAPACK<int, Kokkos::complex<double> >::GESVD(const char& JOBU, const char& JOBVT, const int& m, const int& n, Kokkos::complex<double> * A, const int& lda, double* S, Kokkos::complex<double> * U, const int& ldu, Kokkos::complex<double> * V, const int& ldv, Kokkos::complex<double> * WORK, const int& lwork, double* RWORK, int* info) const {
    ZGESVD_F77(CHAR_MACRO(JOBU), CHAR_MACRO(JOBVT), &m, &n, reinterpret_cast<std::complex<double>*>(A), &lda, S, reinterpret_cast<std::complex<double>*>(U), &ldu, reinterpret_cast<std::complex<double>*>(V), &ldv, reinterpret_cast<std::complex<double>*>(WORK), &lwork, RWORK, info);
  }

  void LAPACK<int, Kokkos::complex<double> >::GEEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* W, Kokkos::complex<double>* VL, const int& ldvl, Kokkos::complex<double>* VR, const int& ldvr, int* ilo, int* ihi, double* SCALE, double* abnrm, double* RCONDE, double* RCONDV, Kokkos::complex<double>* WORK, const int& lwork, double* RWORK, int* info) const
  {
    ZGEEVX_F77(CHAR_MACRO(BALANC), CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SENSE), &n, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(W), reinterpret_cast<std::complex<double>*>(VL), &ldvl, reinterpret_cast<std::complex<double>*>(VR), &ldvr, ilo, ihi, SCALE, abnrm, RCONDE, RCONDV, reinterpret_cast<std::complex<double>*>(WORK), &lwork, RWORK, info);
  }

  void LAPACK<int, Kokkos::complex<double> >::GGEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* B, const int& ldb, Kokkos::complex<double>* ALPHA, Kokkos::complex<double>* BETA, Kokkos::complex<double>* VL, const int& ldvl, Kokkos::complex<double>* VR, const int& ldvr, int* ilo, int* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, double* RCONDE, double* RCONDV, Kokkos::complex<double>* WORK, const int& lwork, double* RWORK, int* IWORK, int* BWORK, int* info) const
  {
    ZGGEVX_F77(CHAR_MACRO(BALANC), CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SENSE), &n, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(B), &ldb, reinterpret_cast<std::complex<double>*>(ALPHA), reinterpret_cast<std::complex<double>*>(BETA), reinterpret_cast<std::complex<double>*>(VL), &ldvl, reinterpret_cast<std::complex<double>*>(VR), &ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, RCONDE, RCONDV, reinterpret_cast<std::complex<double>*>(WORK), &lwork, RWORK, IWORK, BWORK, info);
  }

  void LAPACK<int, Kokkos::complex<double> >::GGEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* B, const int& ldb, double* ALPHAR, double* ALPHAI, Kokkos::complex<double>* BETA, Kokkos::complex<double>* VL, const int& ldvl, Kokkos::complex<double>* VR, const int& ldvr, int* ilo, int* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, double* RCONDE, double* RCONDV, Kokkos::complex<double>* WORK, const int& lwork, double* RWORK, int* IWORK, int* BWORK, int* info) const
  {
    std::vector<Kokkos::complex<double> > w (n);
    Kokkos::complex<double>* w_rawPtr = (n == 0) ? NULL : &w[0];
    GGEVX(BALANC, JOBVL, JOBVR, SENSE, n, A, lda, B, ldb, w_rawPtr, BETA, VL, ldvl, VR, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, RCONDE, RCONDV, WORK, lwork, RWORK, IWORK, BWORK, info);
    if (*info == 0) {
      // The eigenvalues are only valid on output if INFO is zero.
      // Otherwise, we shouldn't even write to WR or WI.
      for (int k = 0; k < n; ++k) {
        ALPHAR[k] = w[k].real ();
        ALPHAI[k] = w[k].imag ();
      }
    }
  }

 void LAPACK<int, Kokkos::complex<double> >::GGEV(const char& JOBVL, const char& JOBVR, const int& n, Kokkos::complex<double> *A, const int& lda, Kokkos::complex<double> *B, const int& ldb, Kokkos::complex<double>* ALPHA, Kokkos::complex<double>* BETA, Kokkos::complex<double>* VL, const int& ldvl, Kokkos::complex<double>* VR, const int& ldvr, Kokkos::complex<double> *WORK, const int& lwork, double* RWORK, int* info) const
  {
    ZGGEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), &n, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(B), &ldb, reinterpret_cast<std::complex<double>*>(ALPHA), reinterpret_cast<std::complex<double>*>(BETA), reinterpret_cast<std::complex<double>*>(VL), &ldvl, reinterpret_cast<std::complex<double>*>(VR), &ldvr, reinterpret_cast<std::complex<double>*>(WORK), &lwork, RWORK, info); 
  }

  void LAPACK<int, Kokkos::complex<double> >::GGES(const char& JOBVL, const char& JOBVR, const char& SORT, int (*ptr2func)(std::complex<double>*, std::complex<double>*), const int& n, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* B, const int& ldb, int* sdim, Kokkos::complex<double>* ALPHA, Kokkos::complex<double>* BETA, Kokkos::complex<double>* VL, const int& ldvl, Kokkos::complex<double>* VR, const int& ldvr, Kokkos::complex<double>* WORK, const int& lwork, double* rwork, int* bwork, int* info ) const
  {
    ZGGES_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(JOBVR), CHAR_MACRO(SORT), ptr2func, &n, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(B), &ldb, sdim, reinterpret_cast<std::complex<double>*>(ALPHA), reinterpret_cast<std::complex<double>*>(BETA), reinterpret_cast<std::complex<double>*>(VL), &ldvl, reinterpret_cast<std::complex<double>*>(VR), &ldvr, reinterpret_cast<std::complex<double>*>(WORK), &lwork, rwork, bwork, info); 
  }
  void LAPACK<int, Kokkos::complex<double> >::TGSEN(const int& ijob, const int& wantq, const int& wantz, const int* SELECT, const int& n, Kokkos::complex<double>* A, const int& lda, Kokkos::complex<double>* B, const int& ldb, Kokkos::complex<double>* ALPHA, Kokkos::complex<double>* BETA, Kokkos::complex<double>* Q, const int& ldq, Kokkos::complex<double>* Z, const int& ldz, int* M, double* PL, double* PR, double* DIF, Kokkos::complex<double>* WORK, const int& lwork, int* IWORK, const int& liwork, int* info ) const
  {
    ZTGSEN_F77(&ijob, &wantq, &wantz, SELECT, &n, reinterpret_cast<std::complex<double>*>(A), &lda, reinterpret_cast<std::complex<double>*>(B), &ldb, reinterpret_cast<std::complex<double>*>(ALPHA), reinterpret_cast<std::complex<double>*>(BETA), reinterpret_cast<std::complex<double>*>(Q), &ldq, reinterpret_cast<std::complex<double>*>(Z), &ldz, M, PL, PR, DIF, reinterpret_cast<std::complex<double>*>(WORK), &lwork, IWORK, &liwork, info); 
  }

  void LAPACK<int, Kokkos::complex<double> >::TREVC(const char& SIDE, const char& HOWMNY, int* select, const int& n, const Kokkos::complex<double>* T, const int& ldt, Kokkos::complex<double>* VL, const int& ldvl, Kokkos::complex<double>* VR, const int& ldvr, const int& mm, int* m, Kokkos::complex<double>* WORK, double* RWORK, int* info) const
  {
    ZTREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), select, &n, reinterpret_cast<const std::complex<double>*>(T), &ldt, reinterpret_cast<std::complex<double>*>(VL), &ldvl, reinterpret_cast<std::complex<double>*>(VR), &ldvr, &mm, m, reinterpret_cast<std::complex<double>*>(WORK), RWORK, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::TREVC(const char& SIDE, const int& n, const Kokkos::complex<double>* T, const int& ldt, Kokkos::complex<double>* VL, const int& ldvl, Kokkos::complex<double>* VR, const int& ldvr, const int& mm, int* m, Kokkos::complex<double>* WORK, double* RWORK, int* info) const
  {
    std::vector<int> select(1);
    const char& whch = 'A';
    ZTREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(whch), &select[0], &n, reinterpret_cast<const std::complex<double>*>(T), &ldt, reinterpret_cast<std::complex<double>*>(VL), &ldvl, reinterpret_cast<std::complex<double>*>(VR), &ldvr, &mm, m, reinterpret_cast<std::complex<double>*>(WORK), RWORK, info);
  }
  
  void LAPACK<int, Kokkos::complex<double> >::TREXC(const char& COMPQ, const int& n, Kokkos::complex<double>* T, const int& ldt, Kokkos::complex<double>* Q, const int& ldq, int* ifst, int* ilst, Kokkos::complex<double>* WORK, int* info) const
  {
    ZTREXC_F77(CHAR_MACRO(COMPQ), &n, reinterpret_cast<std::complex<double>*>(T), &ldt, reinterpret_cast<std::complex<double>*>(Q), &ldq, ifst, ilst, info);
  }

  void LAPACK<int, Kokkos::complex<double> >::LARTG( const Kokkos::complex<double> f, const Kokkos::complex<double> g, double* c, Kokkos::complex<double>* s, Kokkos::complex<double>* r ) const
  {
    ZLARTG_F77(reinterpret_cast<const std::complex<double>*>(&f), reinterpret_cast<const std::complex<double>*>(&g), c, reinterpret_cast<std::complex<double>*>(s), reinterpret_cast<std::complex<double>*>(r));
  }


  void LAPACK<int, Kokkos::complex<double> >::LARFG( const int& n, Kokkos::complex<double>* alpha, Kokkos::complex<double>* x, const int& incx, Kokkos::complex<double>* tau ) const
  {
    ZLARFG_F77(&n, reinterpret_cast<std::complex<double>*>(alpha), reinterpret_cast<std::complex<double>*>(x), &incx, reinterpret_cast<std::complex<double>*>(tau));
  }

  void LAPACK<int, Kokkos::complex<double> >::GEBAL(const char& JOBZ, const int& n, Kokkos::complex<double>* A, const int& lda, int* ilo, int* ihi, double* scale, int* info) const
  {
    ZGEBAL_F77(CHAR_MACRO(JOBZ),&n, reinterpret_cast<std::complex<double>*>(A), &lda, ilo, ihi, scale, info);
  }


  void LAPACK<int, Kokkos::complex<double> >::GEBAK(const char& JOBZ, const char& SIDE, const int& n, const int& ilo, const int& ihi, const double* scale, const int& m, Kokkos::complex<double>* V, const int& ldv, int* info) const
  {
    ZGEBAK_F77(CHAR_MACRO(JOBZ), CHAR_MACRO(SIDE), &n, &ilo, &ihi, scale, &m, reinterpret_cast<std::complex<double>*>(V), &ldv, info);
  }


#ifdef HAVE_TEUCHOS_LAPACKLARND
  Kokkos::complex<double> LAPACK<int, Kokkos::complex<double> >::LARND( const int& idist, int* seed ) const
  {
    double _Complex z = ZLARND_F77(&idist, seed);
    return reinterpret_cast<Kokkos::complex<double>&>(z);
  }
#endif

  void LAPACK<int, Kokkos::complex<double> >::LARNV( const int& idist, int* seed, const int& n, Kokkos::complex<double>* v ) const
  {
    ZLARNV_F77(&idist, seed, &n, reinterpret_cast<std::complex<double>*>(v));
  }


  int LAPACK<int, Kokkos::complex<double> >::ILAENV( const int& ispec, const std::string& NAME, const std::string& OPTS, const int& N1, const int& N2, const int& N3, const int& N4 ) const
  {
    unsigned int opts_length = OPTS.length();
    std::string temp_NAME = "z" + NAME;
    unsigned int name_length = temp_NAME.length();
    return ilaenv_wrapper(&ispec, &temp_NAME[0], name_length, &OPTS[0], opts_length, &N1, &N2, &N3, &N4);
  }

  // END INT, KOKKOS::COMPLEX<DOUBLE> SPECIALIZATION IMPLEMENTATION //

#endif // HAVE_TEUCHOS_COMPLEX


#ifdef HAVE_TEUCHOSCORE_QUADMATH

  // BEGIN int, __float128 SPECIALIZATION IMPLEMENTATION //

  void LAPACK<int, __float128>::
  GEQRF(const int& m, const int& n, __float128* A, const int& lda, __float128* TAU, __float128* WORK, const int& lwork, int* info) const
  {
    Teuchos::Details::Lapack128 lapack;
    lapack.GEQRF (m, n, A, lda, TAU, WORK, lwork, info);
  }

  void LAPACK<int, __float128>::
  GEQR2(const int& m, const int& n, __float128* A, const int& lda, __float128* TAU, __float128* WORK, int* const info) const
  {
    Teuchos::Details::Lapack128 lapack;
    lapack.GEQR2 (m, n, A, lda, TAU, WORK, info);
  }

  void LAPACK<int, __float128>::
  GETRF(const int& m, const int& n, __float128* A, const int& lda, int* IPIV, int* info) const
  {
    Teuchos::Details::Lapack128 lapack;
    lapack.GETRF (m, n, A, lda, IPIV, info);
  }

  void LAPACK<int, __float128>::
  GETRS(const char& TRANS, const int& n, const int& nrhs, const __float128* A, const int& lda, const int* IPIV, __float128* B, const int& ldb, int* info) const
  {
    Teuchos::Details::Lapack128 lapack;
    lapack.GETRS (TRANS, n, nrhs, A, lda, IPIV, B, ldb, info);
  }

  void LAPACK<int, __float128>::
  GETRI (const int& n, __float128* A, const int& lda, const int* IPIV, __float128* WORK, const int& lwork, int* info) const
  {
    Teuchos::Details::Lapack128 lapack;
    lapack.GETRI (n, A, lda, const_cast<int*> (IPIV), WORK, lwork, info);
  }

  void LAPACK<int, __float128>::
  LASWP (const int& N, __float128* A, const int& LDA, const int& K1, const int& K2, const int* IPIV, const int& INCX) const
  {
    Teuchos::Details::Lapack128 lapack;
    lapack.LASWP (N, A, LDA, K1, K2, IPIV, INCX);
  }

  void LAPACK<int, __float128>::
  ORM2R(const char& SIDE, const char& TRANS, const int& m, const int& n, const int& k, const __float128* A, const int& lda, const __float128* TAU, __float128* C, const int& ldc, __float128* WORK, int* const info) const
  {
    Teuchos::Details::Lapack128 lapack;
    lapack.ORM2R (SIDE, TRANS, m, n, k, A, lda, TAU, C, ldc, WORK, info);
  }

  void LAPACK<int, __float128>::
  ORGQR(const int& m, const int& n, const int& k, __float128* A, const int& lda, const __float128* TAU, __float128* WORK, const int& lwork, int* info) const
  {
    Teuchos::Details::Lapack128 lapack;
    lapack.ORGQR (m, n, k, A, lda, TAU, WORK, lwork, info);
  }

  void LAPACK<int, __float128>::
  UNGQR(const int& m, const int& n, const int& k, __float128* A, const int& lda, const __float128* TAU, __float128* WORK, const int& lwork, int* info) const
  {
    Teuchos::Details::Lapack128 lapack;
    lapack.UNGQR (m, n, k, A, lda, TAU, WORK, lwork, info);
  }

  void LAPACK<int, __float128>::
  LARFG( const int& n, __float128* alpha, __float128* x, const int& incx, __float128* tau ) const
  {
    Teuchos::Details::Lapack128 lapack;
    lapack.LARFG (n, alpha, x, incx, tau);
  }

  __float128 LAPACK<int, __float128>::
  LAPY2 (const __float128 x, const __float128 y) const
  {
    Teuchos::Details::Lapack128 lapack;
    return lapack.LAPY2 (x, y);
  }

  void LAPACK<int, __float128>::
  GBTRF (const int& m, const int& n, const int& kl, const int& ku,
         __float128* A, const int& lda, int* IPIV, int* info) const
  {
    Teuchos::Details::Lapack128 lapack;
    return lapack.GBTRF (m, n, kl, ku, A, lda, IPIV, info);
  }

  void LAPACK<int, __float128>::
  GBTRS (const char& TRANS, const int& n, const int& kl, const int& ku,
         const int& nrhs, const __float128* A, const int& lda, const int* IPIV,
         __float128* B, const int& ldb, int* info) const
  {
    Teuchos::Details::Lapack128 lapack;
    return lapack.GBTRS (TRANS, n, kl, ku, nrhs, A, lda, IPIV, B, ldb, info);
  }

  void LAPACK<int, __float128>::
  LASCL (const char& TYPE, const int& kl, const int& ku, const __float128 cfrom,
         const __float128 cto, const int& m, const int& n, __float128* A,
         const int& lda, int* info) const
  {
    Teuchos::Details::Lapack128 lapack;
    return lapack.LASCL (TYPE, kl, ku, cfrom, cto, m, n, A, lda, info);
  }

  // END int, __float128 SPECIALIZATION IMPLEMENTATION //

#endif // HAVE_TEUCHOSCORE_QUADMATH

#ifdef HAVE_TEUCHOS_LONG_DOUBLE

  // BEGIN int, long double SPECIALIZATION IMPLEMENTATION //

  void LAPACK<int, long double>::
  GESV(const int& n, const int& nrhs, long double* A, const int& lda, int* IPIV, long double* B, const int& ldb, int* info) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR in Teuchos::LAPACK: GESV not implemented for long double scalar type!");
  }
  void LAPACK<int, long double>::  
  GTTRS(const char& TRANS, const int& n, const int& nrhs, const long double* dl, const long double* d, const long double* du, const long double* du2, const int* IPIV, long double* B, const int& ldb, int* info) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR in Teuchos::LAPACK: GTTRS not implemented for long double scalar type!");
  }
  void LAPACK<int, long double>::  
  GTTRF(const int& n, long double* dl, long double* d, long double* du, long double* du2, int* IPIV, int* info) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR in Teuchos::LAPACK: GTTRF not implemented for long double scalar type!");
  }
  void LAPACK<int, long double>::  
  SYGV(const int& itype, const char& JOBZ, const char& UPLO, const int& n, long double* A, const int& lda, long double* B, const int& ldb, long double* W, long double* WORK, const int& lwork, int* info) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR in Teuchos::LAPACK: GESV not implemented for long double scalar type!");
  }
  void LAPACK<int, long double>::
  GEEV(const char& JOBVL, const char& JOBVR, const int& n, long double* A, const int& lda, long double* WR, long double* WI, long double* VL, const int& ldvl, long double* VR, const int& ldvr, long double* WORK, const int& lwork, int* info) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR in Teuchos::LAPACK: GEEV not implemented for long double scalar type!");
  }
  void LAPACK<int, long double>::
  GEEV(const char& JOBVL, const char& JOBVR, const int& n, long double* A, const int& lda, long double* WR, long double* WI, long double* VL, const int& ldvl, long double* VR, const int& ldvr, long double* WORK, const int& lwork, long double* /* RWORK */, int* info) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR in Teuchos::LAPACK: GEEV not implemented for long double scalar type!");
  }
  void LAPACK<int,long double>::
  GGEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, long double* A, const int& lda, long double* B, const int& ldb, long double* ALPHAR, long double* ALPHAI, long double* BETA, long double* VL, const int& ldvl, long double* VR, const int& ldvr, int* ilo, int* ihi, long double* lscale, long double* rscale, long double* abnrm, long double* bbnrm, long double* RCONDE, long double* RCONDV, long double* WORK, const int& lwork, int* IWORK, int* BWORK, int* info) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR in Teuchos::LAPACK: GGEVX not implemented for long double scalar type!");
  }
  void LAPACK<int,long double>::
  GGEVX(const char& BALANC, const char& JOBVL, const char& JOBVR, const char& SENSE, const int& n, long double* A, const int& lda, long double* B, const int& ldb, long double* ALPHAR, long double* ALPHAI, long double* BETA, long double* VL, const int& ldvl, long double* VR, const int& ldvr, int* ilo, int* ihi, long double* lscale, long double* rscale, long double* abnrm, long double* bbnrm, long double* RCONDE, long double* RCONDV, long double* WORK, const int& lwork, long double* /* RWORK */, int* IWORK, int* BWORK, int* info) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR in Teuchos::LAPACK: GGEVX not implemented for long double scalar type!");
  }
  void LAPACK<int, long double>::
  PORFS(const char& UPLO, const int& n, const int& nrhs, const long double* A, const int& lda, const long double* AF, const int& ldaf, const long double* B, const int& ldb, long double* X, const int& ldx, long double* FERR, long double* BERR, long double* WORK, int* IWORK, int* info) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR in Teuchos::LAPACK: PORFS not implemented for long double scalar type!");
  }
  void LAPACK<int,long double>::
  PTEQR(const char& COMPZ, const int& n, long double* D, long double* E, long double* Z, const int& ldz, long double* WORK, int* info) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR in Teuchos::LAPACK: PTEQR not implemented for long double scalar type!");
  }
  void LAPACK<int, long double>::
  POTRF(const char& UPLO, const int& n, long double* A, const int& lda, int* info) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR in Teuchos::LAPACK: POTRF not implemented for long double scalar type!");
  }
  void LAPACK<int, long double>::
  POTRS(const char& UPLO, const int& n, const int& nrhs, const long double* A, const int& lda, long double* B, const int& ldb, int* info) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR in Teuchos::LAPACK: POTRS not implemented for long double scalar type!");
  }
  void LAPACK<int, long double>::
  POEQU(const int& n, const long double* A, const int& lda, long double* S, long double* scond, long double* amax, int* info) const 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR in Teuchos::LAPACK: POEQU not implemented for long double scalar type!");
  }
  void LAPACK<int, long double>::
  GEQRF(const int& m, const int& n, long double* A, const int& lda, long double* TAU, long double* WORK, const int& lwork, int* info) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    lapack.GEQRF (m, n, A, lda, TAU, WORK, lwork, info);
  }

  void LAPACK<int, long double>::
  GEQR2(const int& m, const int& n, long double* A, const int& lda, long double* TAU, long double* WORK, int* const info) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    lapack.GEQR2 (m, n, A, lda, TAU, WORK, info);
  }

  void LAPACK<int, long double>::
  GETRF(const int& m, const int& n, long double* A, const int& lda, int* IPIV, int* info) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    lapack.GETRF (m, n, A, lda, IPIV, info);
  }

  void LAPACK<int, long double>::
  GETRS(const char& TRANS, const int& n, const int& nrhs, const long double* A, const int& lda, const int* IPIV, long double* B, const int& ldb, int* info) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    lapack.GETRS (TRANS, n, nrhs, A, lda, IPIV, B, ldb, info);
  }

  void LAPACK<int, long double>::
  GETRI (const int& n, long double* A, const int& lda, const int* IPIV, long double* WORK, const int& lwork, int* info) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    lapack.GETRI (n, A, lda, const_cast<int*> (IPIV), WORK, lwork, info);
  }

  void LAPACK<int, long double>::
  LASWP (const int& N, long double* A, const int& LDA, const int& K1, const int& K2, const int* IPIV, const int& INCX) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    lapack.LASWP (N, A, LDA, K1, K2, IPIV, INCX);
  }

  void LAPACK<int, long double>::
  ORM2R(const char& SIDE, const char& TRANS, const int& m, const int& n, const int& k, const long double* A, const int& lda, const long double* TAU, long double* C, const int& ldc, long double* WORK, int* const info) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    lapack.ORM2R (SIDE, TRANS, m, n, k, A, lda, TAU, C, ldc, WORK, info);
  }

  void LAPACK<int, long double>::
  ORGQR(const int& m, const int& n, const int& k, long double* A, const int& lda, const long double* TAU, long double* WORK, const int& lwork, int* info) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    lapack.ORGQR (m, n, k, A, lda, TAU, WORK, lwork, info);
  }

  void LAPACK<int, long double>::
  UNGQR(const int& m, const int& n, const int& k, long double* A, const int& lda, const long double* TAU, long double* WORK, const int& lwork, int* info) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    lapack.UNGQR (m, n, k, A, lda, TAU, WORK, lwork, info);
  }

  void LAPACK<int, long double>::
  LARFG( const int& n, long double* alpha, long double* x, const int& incx, long double* tau ) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    lapack.LARFG (n, alpha, x, incx, tau);
  }

  long double LAPACK<int, long double>::
  LAPY2 (const long double x, const long double y) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    return lapack.LAPY2 (x, y);
  }

  void LAPACK<int, long double>::
  GBTRF (const int& m, const int& n, const int& kl, const int& ku,
         long double* A, const int& lda, int* IPIV, int* info) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    return lapack.GBTRF (m, n, kl, ku, A, lda, IPIV, info);
  }

  void LAPACK<int, long double>::
  GBTRS (const char& TRANS, const int& n, const int& kl, const int& ku,
         const int& nrhs, const long double* A, const int& lda, const int* IPIV,
         long double* B, const int& ldb, int* info) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    return lapack.GBTRS (TRANS, n, kl, ku, nrhs, A, lda, IPIV, B, ldb, info);
  }

  void LAPACK<int, long double>::
  LASCL (const char& TYPE, const int& kl, const int& ku, const long double cfrom,
         const long double cto, const int& m, const int& n, long double* A,
         const int& lda, int* info) const
  {
    Teuchos::Details::LapackLongDouble lapack;
    return lapack.LASCL (TYPE, kl, ku, cfrom, cto, m, n, A, lda, info);
  }

  // END int, long double SPECIALIZATION IMPLEMENTATION //

#endif // HAVE_TEUCHOS_LONG_DOUBLE


} // namespace Teuchos
