// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// Kris
// 06.11.03 -- Format cleanup
// 06.17.03 -- Added LAPY2 and GEES by request
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_LAPACK_WRAPPERS_HPP_
#define _TEUCHOS_LAPACK_WRAPPERS_HPP_

#include "Teuchos_ConfigDefs.hpp"

/*! \file Teuchos_LAPACK_wrappers.hpp
    \brief The Templated LAPACK wrappers
*/
/* Define fcd (Fortran Teuchos_fcd descriptor) for non-standard situations */

#if defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL)

#if defined(CRAY_T3X)

#include <fortran.h>
#define PREFIX
#define Teuchos_fcd fcd 

#define DGEQRF_F77  F77_FUNC(sgeqrf,SGEQRF)
#define DGETRF_F77  F77_FUNC(sgetrf,SGETRF)
#define DGETRS_F77  F77_FUNC(sgetrs,SGETRS)
#define DGETRI_F77  F77_FUNC(sgetri,SGETRI)
#define DGERFS_F77  F77_FUNC(sgerfs,SGERFS)
#define DGECON_F77  F77_FUNC(sgecon,SGECON)
#define DGESVX_F77  F77_FUNC(sgesvx,SGESVX)
#define DGESV_F77   F77_FUNC(sgesv,SGESV)
#define DGEEQU_F77  F77_FUNC(sgeequ,SGEEQU)
#define DSYTRD_F77  F77_FUNC(ssytrd,SSYTRD)
#define DPOTRF_F77  F77_FUNC(spotrf,SPOTRF)
#define DPOTRS_F77  F77_FUNC(spotrs,SPOTRS)
#define DPOTRI_F77  F77_FUNC(spotri,SPOTRI)
#define DPOCON_F77  F77_FUNC(spocon,SPOCON)
#define DPOSV_F77   F77_FUNC(sposv,SPOSV)
#define DPOEQU_F77  F77_FUNC(spoequ,SPOEQU)
#define DPORFS_F77  F77_FUNC(sporfs,SPORFS)
#define DPOSVX_F77  F77_FUNC(sposvx,SPOSVX)
#define DLAMCH_F77  F77_FUNC(slamch,SLAMCH)
#define DGELS_F77   F77_FUNC(sgels,SGELS)
#define DGEEV_F77   F77_FUNC(sgeev,SGEEV)
#define DGEHRD_F77  F77_FUNC(sgehrd,SGEHRD)
#define DHSEQR_F77  F77_FUNC(shseqr,SHSEQR)
#define DORMQR_F77  F77_FUNC(sormqr,SORMQR)
#define DORGHR_F77  F77_FUNC(sorghr,SORGHR)
#define DORMHR_F77  F77_FUNC(sormhr,SORMHR)
#define DTREVC_F77  F77_FUNC(strevc,STREVC)
#define DTREXC_F77  F77_FUNC(strexc,STREXC)
#define DGEES_F77   F77_FUNC(sgees,SGEES)
#define DSPEV_F77   F77_FUNC(sspev,SSPEV)
#define DSYEV_F77   F77_FUNC(ssyev,SSYEV)
#define DSYGV_F77   F77_FUNC(ssygv,SSYGV)
#define DSTEQR_F77  F77_FUNC(ssteqr,SSTEQR)
#define DLAPY2_F77  F77_FUNC(slapy2,SLAPY2)
#define DLARND_F77  F77_FUNC(slarnd,SLARND)
#define DLARNV_F77  F77_FUNC(slarnv,SLARNV)

#ifdef HAVE_TEUCHOS_COMPLEX

#define ZGEQRF_F77  F77_FUNC(cgeqrf,CGEQRF)
#define ZGETRF_F77  F77_FUNC(cgetrf,CGETRF)
#define ZGETRS_F77  F77_FUNC(cgetrs,CGETRS)
#define ZGETRI_F77  F77_FUNC(cgetri,CGETRI)
#define ZGERFS_F77  F77_FUNC(cgerfs,CGERFS)
#define ZGECON_F77  F77_FUNC(cgecon,CGECON)
#define ZGESVX_F77  F77_FUNC(cgesvx,CGESVX)
#define ZGESV_F77   F77_FUNC(cgesv,CGESV)
#define ZGEEQU_F77  F77_FUNC(cgeequ,CGEEQU)
#define ZPOTRF_F77  F77_FUNC(cpotrf,CPOTRF)
#define ZPOTRS_F77  F77_FUNC(cpotrs,CPOTRS)
#define ZPOTRI_F77  F77_FUNC(cpotri,CPOTRI)
#define ZPOCON_F77  F77_FUNC(cpocon,CPOCON)
#define ZPOSV_F77   F77_FUNC(cposv,CPOSV)
#define ZPOEQU_F77  F77_FUNC(cpoequ,CPOEQU)
#define ZPORFS_F77  F77_FUNC(cporfs,CPORFS)
#define ZPOSVX_F77  F77_FUNC(cposvx,CPOSVX)
#define ZGELS_F77   F77_FUNC(cgels,CGELS)
#define ZGEEV_F77   F77_FUNC(cgeev,CGEEV)
#define ZGEHRD_F77  F77_FUNC(cgehrd,CGEHRD)
#define ZHSEQR_F77  F77_FUNC(chseqr,CHSEQR)
#define ZTREVC_F77  F77_FUNC(ctrevc,CTREVC)
#define ZTREXC_F77  F77_FUNC(ctrexc,CTREXC)
#define ZGEES_F77   F77_FUNC(cgees,CGEES)
#define ZSTEQR_F77  F77_FUNC(csteqr,CSTEQR)
#define ZLARND_F77  F77_FUNC(clarnd,CLARND)
#define ZLARNV_F77  F77_FUNC(clarnv,CLARNV)

#endif

#elif defined(INTEL_CXML)

#define PREFIX __stdcall 
#define Teuchos_fcd const char *, unsigned int 

#define DGEQRF_F77  F77_FUNC(dgeqrf,DGEQRF)
#define DGETRF_F77  F77_FUNC(dgetrf,DGETRF)
#define DGETRS_F77  F77_FUNC(dgetrs,DGETRS)
#define DGETRI_F77  F77_FUNC(dgetri,DGETRI)
#define DGERFS_F77  F77_FUNC(dgerfs,DGERFS)
#define DGECON_F77  F77_FUNC(dgecon,DGECON)
#define DGESVX_F77  F77_FUNC(dgesvx,DGESVX)
#define DGESV_F77   F77_FUNC(dgesv,DGESV)
#define DGEEQU_F77  F77_FUNC(dgeequ,DGEEQU)
#define DSYTRD_F77  F77_FUNC(dsytrd,DSYTRD)
#define DPOTRF_F77  F77_FUNC(dpotrf,DPOTRF)
#define DPOTRS_F77  F77_FUNC(dpotrs,DPOTRS)
#define DPOTRI_F77  F77_FUNC(dpotri,DPOTRI)
#define DPOCON_F77  F77_FUNC(dpocon,DPOCON)
#define DPOSV_F77   F77_FUNC(dposv,DPOSV)
#define DPOEQU_F77  F77_FUNC(dpoequ,DPOEQU)
#define DPORFS_F77  F77_FUNC(dporfs,DPORFS)
#define DPOSVX_F77  F77_FUNC(dposvx,DPOSVX)
#define DLAMCH_F77  F77_FUNC(dlamch,DLAMCH)
#define DGELS_F77   F77_FUNC(dgels,DGELS)
#define DGEEV_F77   F77_FUNC(dgeev,DGEEV)
#define DGEHRD_F77  F77_FUNC(dgehrd,DGEHRD)
#define DHSEQR_F77  F77_FUNC(dhseqr,DHSEQR)
#define DORGHR_F77  F77_FUNC(dorghr,DORGHR)
#define DORMHR_F77  F77_FUNC(dormhr,DORMHR)
#define DORMQR_F77  F77_FUNC(dormqr,DORMQR)
#define DTREVC_F77  F77_FUNC(dtrevc,DTREVC)
#define DTREXC_F77  F77_FUNC(dtrexc,DTREXC)
#define DGEES_F77   F77_FUNC(dgees,DGEES)
#define DSPEV_F77   F77_FUNC(dspev,DSPEV)
#define DSYEV_F77   F77_FUNC(dsyev,DSYEV)
#define DSYGV_F77   F77_FUNC(dsygv,DSYGV)
#define DSTEQR_F77  F77_FUNC(dsteqr,DSTEQR)
#define DLAPY2_F77  F77_FUNC(dlapy2,DLAPY2)
#define DLARND_F77  F77_FUNC(dlarnd,DLARND)
#define DLARNV_F77  F77_FUNC(dlarnv,DLARNV)

#ifdef HAVE_TEUCHOS_COMPLEX

#define ZGEQRF_F77  F77_FUNC(zgeqrf,ZGEQRF)
#define ZGETRF_F77  F77_FUNC(zgetrf,ZGETRF)
#define ZGETRS_F77  F77_FUNC(zgetrs,ZGETRS)
#define ZGETRI_F77  F77_FUNC(zgetri,ZGETRI)
#define ZGERFS_F77  F77_FUNC(zgerfs,ZGERFS)
#define ZGECON_F77  F77_FUNC(zgecon,ZGECON)
#define ZGESVX_F77  F77_FUNC(zgesvx,ZGESVX)
#define ZGESV_F77   F77_FUNC(zgesv,ZGESV)
#define ZGEEQU_F77  F77_FUNC(zgeequ,ZGEEQU)
#define ZPOTRF_F77  F77_FUNC(zpotrf,ZPOTRF)
#define ZPOTRS_F77  F77_FUNC(zpotrs,ZPOTRS)
#define ZPOTRI_F77  F77_FUNC(zpotri,ZPOTRI)
#define ZPOCON_F77  F77_FUNC(zpocon,ZPOCON)
#define ZPOSV_F77   F77_FUNC(zposv,ZPOSV)
#define ZPOEQU_F77  F77_FUNC(zpoequ,ZPOEQU)
#define ZPORFS_F77  F77_FUNC(zporfs,ZPORFS)
#define ZPOSVX_F77  F77_FUNC(zposvx,ZPOSVX)
#define ZGELS_F77   F77_FUNC(zgels,ZGELS)
#define ZGEEV_F77   F77_FUNC(zgeev,ZGEEV)
#define ZGEHRD_F77  F77_FUNC(zgehrd,ZGEHRD)
#define ZHSEQR_F77  F77_FUNC(zhseqr,ZHSEQR)
#define ZTREVC_F77  F77_FUNC(ztrevc,ZTREVC)
#define ZTREXC_F77  F77_FUNC(ztrexc,ZTREXC)
#define ZGEES_F77   F77_FUNC(zgees,ZGEES)
#define ZSTEQR_F77  F77_FUNC(zsteqr,ZSTEQR)
#define ZLARND_F77  F77_FUNC(zlarnd,ZLARND)
#define ZLARNV_F77  F77_FUNC(zlarnv,ZLARNV)

#endif

#elif defined(INTEL_MKL)

#define PREFIX
#define Teuchos_fcd const char *

#define DGEQRF_F77  F77_FUNC(dgeqrf,DGEQRF)
#define DGETRF_F77  F77_FUNC(dgetrf,DGETRF)
#define DGETRS_F77  F77_FUNC(dgetrs,DGETRS)
#define DGETRI_F77  F77_FUNC(dgetri,DGETRI)
#define DGERFS_F77  F77_FUNC(dgerfs,DGERFS)
#define DGECON_F77  F77_FUNC(dgecon,DGECON)
#define DGESVX_F77  F77_FUNC(dgesvx,DGESVX)
#define DGESV_F77   F77_FUNC(dgesv,DGESV)
#define DGEEQU_F77  F77_FUNC(dgeequ,DGEEQU)
#define DSYTRD_F77  F77_FUNC(dsytrd,DSYTRD)
#define DPOTRF_F77  F77_FUNC(dpotrf,DPOTRF)
#define DPOTRS_F77  F77_FUNC(dpotrs,DPOTRS)
#define DPOTRI_F77  F77_FUNC(dpotri,DPOTRI)
#define DPOCON_F77  F77_FUNC(dpocon,DPOCON)
#define DPOSV_F77   F77_FUNC(dposv,DPOSV)
#define DPOEQU_F77  F77_FUNC(dpoequ,DPOEQU)
#define DPORFS_F77  F77_FUNC(dporfs,DPORFS)
#define DPOSVX_F77  F77_FUNC(dposvx,DPOSVX)
#define DLAMCH_F77  F77_FUNC(dlamch,DLAMCH)
#define DGELS_F77   F77_FUNC(dgels,DGELS)
#define DGEEV_F77   F77_FUNC(dgeev,DGEEV)
#define DGEHRD_F77  F77_FUNC(dgehrd,DGEHRD)
#define DHSEQR_F77  F77_FUNC(dhseqr,DHSEQR)
#define DORGHR_F77  F77_FUNC(dorghr,DORGHR)
#define DORMHR_F77  F77_FUNC(dormhr,DORMHR)
#define DORMQR_F77  F77_FUNC(dormqr,DORMQR)
#define DTREVC_F77  F77_FUNC(dtrevc,DTREVC)
#define DTREXC_F77  F77_FUNC(dtrexc,DTREXC)
#define DGEES_F77   F77_FUNC(dgees,DGEES)
#define DSPEV_F77   F77_FUNC(dspev,DSPEV)
#define DSYEV_F77   F77_FUNC(dsyev,DSYEV)
#define DSYGV_F77   F77_FUNC(dsygv,DSYGV)
#define DSTEQR_F77  F77_FUNC(dsteqr,DSTEQR)
#define DLAPY2_F77  F77_FUNC(dlapy2,DLAPY2)
#define DLARND_F77  F77_FUNC(dlarnd,DLARND)
#define DLARNV_F77  F77_FUNC(dlarnv,DLARNV)

#ifdef HAVE_TEUCHOS_COMPLEX

#define ZGEQRF_F77  F77_FUNC(zgeqrf,ZGEQRF)
#define ZGETRF_F77  F77_FUNC(zgetrf,ZGETRF)
#define ZGETRS_F77  F77_FUNC(zgetrs,ZGETRS)
#define ZGETRI_F77  F77_FUNC(zgetri,ZGETRI)
#define ZGERFS_F77  F77_FUNC(zgerfs,ZGERFS)
#define ZGECON_F77  F77_FUNC(zgecon,ZGECON)
#define ZGESVX_F77  F77_FUNC(zgesvx,ZGESVX)
#define ZGESV_F77   F77_FUNC(zgesv,ZGESV)
#define ZGEEQU_F77  F77_FUNC(zgeequ,ZGEEQU)
#define ZPOTRF_F77  F77_FUNC(zpotrf,ZPOTRF)
#define ZPOTRS_F77  F77_FUNC(zpotrs,ZPOTRS)
#define ZPOTRI_F77  F77_FUNC(zpotri,ZPOTRI)
#define ZPOCON_F77  F77_FUNC(zpocon,ZPOCON)
#define ZPOSV_F77   F77_FUNC(zposv,ZPOSV)
#define ZPOEQU_F77  F77_FUNC(zpoequ,ZPOEQU)
#define ZPORFS_F77  F77_FUNC(zporfs,ZPORFS)
#define ZPOSVX_F77  F77_FUNC(zposvx,ZPOSVX)
#define ZGELS_F77   F77_FUNC(zgels,ZGELS)
#define ZGEEV_F77   F77_FUNC(zgeev,ZGEEV)
#define ZGEHRD_F77  F77_FUNC(zgehrd,ZGEHRD)
#define ZHSEQR_F77  F77_FUNC(zhseqr,ZHSEQR)
#define ZTREVC_F77  F77_FUNC(ztrevc,ZTREVC)
#define ZTREXC_F77  F77_FUNC(ztrexc,ZTREXC)
#define ZGEES_F77   F77_FUNC(zgees,ZGEES)
#define ZSTEQR_F77  F77_FUNC(zsteqr,ZSTEQR)
#define ZLARND_F77  F77_FUNC(zlarnd,ZLARND)
#define ZLARNV_F77  F77_FUNC(zlarnv,ZLARNV)

#endif

#endif 

/* All three of these machines use a simple uppercase mangling of Fortran names */

/* if F77_FUNC is defined undefine it because we want to redefine */

#ifdef F77_FUNC
#undef F77_FUNC
#endif

#define F77_FUNC(lcase,UCASE) PREFIX UCASE

#else /* Define Teuchos_fcd for all other machines */

#define PREFIX
#define Teuchos_fcd const char * 

#ifndef HAVE_CONFIG_H

#ifdef F77_FUNC
#undef F77_FUNC
#endif

#ifdef TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE
#define F77_FUNC(lcase,UCASE) lcase
#else /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE not defined*/
#define F77_FUNC(lcase,UCASE) lcase ## _
#endif /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE */

#endif /* HAVE_CONFIG_H */

#define DGEQRF_F77  F77_FUNC(dgeqrf,DGEQRF)
#define DGETRF_F77  F77_FUNC(dgetrf,DGETRF)
#define DGETRS_F77  F77_FUNC(dgetrs,DGETRS)
#define DGETRI_F77  F77_FUNC(dgetri,DGETRI)
#define DGERFS_F77  F77_FUNC(dgerfs,DGERFS)
#define DGECON_F77  F77_FUNC(dgecon,DGECON)
#define DGESVX_F77  F77_FUNC(dgesvx,DGESVX)
#define DGESV_F77   F77_FUNC(dgesv,DGESV)
#define DGEEQU_F77  F77_FUNC(dgeequ,DGEEQU)
#define DSYTRD_F77  F77_FUNC(dsytrd,DSYTRD)
#define DPOTRF_F77  F77_FUNC(dpotrf,DPOTRF)
#define DPOTRS_F77  F77_FUNC(dpotrs,DPOTRS)
#define DPOTRI_F77  F77_FUNC(dpotri,DPOTRI)
#define DPOCON_F77  F77_FUNC(dpocon,DPOCON)
#define DPOSV_F77   F77_FUNC(dposv,DPOSV)
#define DPOEQU_F77  F77_FUNC(dpoequ,DPOEQU)
#define DPORFS_F77  F77_FUNC(dporfs,DPORFS)
#define DPOSVX_F77  F77_FUNC(dposvx,DPOSVX)
#define DLAMCH_F77  F77_FUNC(dlamch,DLAMCH)
#define DGELS_F77   F77_FUNC(dgels,DGELS)
#define DGEEV_F77   F77_FUNC(dgeev,DGEEV)
#define DGEHRD_F77  F77_FUNC(dgehrd,DGEHRD)
#define DHSEQR_F77  F77_FUNC(dhseqr,DHSEQR)
#define DORGHR_F77  F77_FUNC(dorghr,DORGHR)
#define DORMHR_F77  F77_FUNC(dormhr,DORMHR)
#define DORMQR_F77  F77_FUNC(dormqr,DORMQR)
#define DTREVC_F77  F77_FUNC(dtrevc,DTREVC)
#define DTREXC_F77  F77_FUNC(dtrexc,DTREXC)
#define DGEES_F77   F77_FUNC(dgees,DGEES)
#define DSPEV_F77   F77_FUNC(dspev,DSPEV)
#define DSYEV_F77   F77_FUNC(dsyev,DSYEV)
#define DSYGV_F77   F77_FUNC(dsygv,DSYGV)
#define DSTEQR_F77  F77_FUNC(dsteqr,DSTEQR)
#define DLAPY2_F77  F77_FUNC(dlapy2,DLAPY2)
#define DLARND_F77  F77_FUNC(dlarnd,DLARND)
#define DLARNV_F77  F77_FUNC(dlarnv,DLARNV)

#ifdef HAVE_TEUCHOS_COMPLEX

#define ZGEQRF_F77  F77_FUNC(zgeqrf,ZGEQRF)
#define ZGETRF_F77  F77_FUNC(zgetrf,ZGETRF)
#define ZGETRS_F77  F77_FUNC(zgetrs,ZGETRS)
#define ZGETRI_F77  F77_FUNC(zgetri,ZGETRI)
#define ZGERFS_F77  F77_FUNC(zgerfs,ZGERFS)
#define ZGECON_F77  F77_FUNC(zgecon,ZGECON)
#define ZGESVX_F77  F77_FUNC(zgesvx,ZGESVX)
#define ZGESV_F77   F77_FUNC(zgesv,ZGESV)
#define ZGEEQU_F77  F77_FUNC(zgeequ,ZGEEQU)
#define ZPOTRF_F77  F77_FUNC(zpotrf,ZPOTRF)
#define ZPOTRS_F77  F77_FUNC(zpotrs,ZPOTRS)
#define ZPOTRI_F77  F77_FUNC(zpotri,ZPOTRI)
#define ZPOCON_F77  F77_FUNC(zpocon,ZPOCON)
#define ZPOSV_F77   F77_FUNC(zposv,ZPOSV)
#define ZPOEQU_F77  F77_FUNC(zpoequ,ZPOEQU)
#define ZPORFS_F77  F77_FUNC(zporfs,ZPORFS)
#define ZPOSVX_F77  F77_FUNC(zposvx,ZPOSVX)
#define ZGELS_F77   F77_FUNC(zgels,ZGELS)
#define ZGEEV_F77   F77_FUNC(zgeev,ZGEEV)
#define ZGEHRD_F77  F77_FUNC(zgehrd,ZGEHRD)
#define ZHSEQR_F77  F77_FUNC(zhseqr,ZHSEQR)
#define ZTREVC_F77  F77_FUNC(ztrevc,ZTREVC)
#define ZTREXC_F77  F77_FUNC(ztrexc,ZTREXC)
#define ZGEES_F77   F77_FUNC(zgees,ZGEES)
#define ZSTEQR_F77  F77_FUNC(zsteqr,ZSTEQR)
#define ZLARND_F77  F77_FUNC(zlarnd,ZLARND)
#define ZLARNV_F77  F77_FUNC(zlarnv,ZLARNV)

#endif /* HAVE_TEUCHOS_COMPLEX */

#endif

#define SGEQRF_F77  F77_FUNC(sgeqrf,SGEQRF)
#define SGETRF_F77  F77_FUNC(sgetrf,SGETRF)
#define SGETRS_F77  F77_FUNC(sgetrs,SGETRS)
#define SGETRI_F77  F77_FUNC(sgetri,SGETRI)
#define SGERFS_F77  F77_FUNC(sgerfs,SGERFS)
#define SGECON_F77  F77_FUNC(sgecon,SGECON)
#define SGESVX_F77  F77_FUNC(sgesvx,SGESVX)
#define SGESV_F77   F77_FUNC(sgesv,SGESV)
#define SGEEQU_F77  F77_FUNC(sgeequ,SGEEQU)
#define SSYTRD_F77  F77_FUNC(ssytrd,SSYTRD)
#define SPOTRF_F77  F77_FUNC(spotrf,SPOTRF)
#define SPOTRS_F77  F77_FUNC(spotrs,SPOTRS)
#define SPOTRI_F77  F77_FUNC(spotri,SPOTRI)
#define SPOCON_F77  F77_FUNC(spocon,SPOCON)
#define SPOSV_F77   F77_FUNC(sposv,SPOSV)
#define SPOEQU_F77  F77_FUNC(spoequ,SPOEQU)
#define SPORFS_F77  F77_FUNC(sporfs,SPORFS)
#define SPOSVX_F77  F77_FUNC(sposvx,SPOSVX)
#define SGELS_F77   F77_FUNC(sgels,SGELS)
#define SGEEV_F77   F77_FUNC(sgeev,SGEEV)
#define SGEHRD_F77  F77_FUNC(sgehrd,SGEHRD)
#define SHSEQR_F77  F77_FUNC(shseqr,SHSEQR)
#define SORGHR_F77  F77_FUNC(sorghr,SORGHR)
#define SORMHR_F77  F77_FUNC(sormhr,SORMHR)
#define SORMQR_F77  F77_FUNC(sormqr,SORMQR)
#define STREVC_F77  F77_FUNC(strevc,STREVC)
#define STREXC_F77  F77_FUNC(strexc,STREXC)
#define SLAMCH_F77  F77_FUNC(slamch,SLAMCH)
#define SGEES_F77   F77_FUNC(sgees,SGEES)
#define SSPEV_F77   F77_FUNC(sspev,SSPEV)
#define SSYEV_F77   F77_FUNC(ssyev,SSYEV)
#define SSYGV_F77   F77_FUNC(ssygv,SSYGV)
#define SSTEQR_F77  F77_FUNC(ssteqr,SSTEQR)
#define SLAPY2_F77  F77_FUNC(slapy2,SLAPY2)
#define SLARND_F77  F77_FUNC(slarnd,SLARND)
#define SLARNV_F77  F77_FUNC(slarnv,SLARNV)

#ifdef HAVE_TEUCHOS_COMPLEX

#define CGEQRF_F77  F77_FUNC(cgeqrf,CGEQRF)
#define CGETRF_F77  F77_FUNC(cgetrf,CGETRF)
#define CGETRS_F77  F77_FUNC(cgetrs,CGETRS)
#define CGETRI_F77  F77_FUNC(cgetri,CGETRI)
#define CGERFS_F77  F77_FUNC(cgerfs,CGERFS)
#define CGECON_F77  F77_FUNC(cgecon,CGECON)
#define CGESVX_F77  F77_FUNC(cgesvx,CGESVX)
#define CGESV_F77   F77_FUNC(cgesv,CGESV)
#define CGEEQU_F77  F77_FUNC(cgeequ,CGEEQU)
#define CPOTRF_F77  F77_FUNC(cpotrf,CPOTRF)
#define CPOTRS_F77  F77_FUNC(cpotrs,CPOTRS)
#define CPOTRI_F77  F77_FUNC(cpotri,CPOTRI)
#define CPOCON_F77  F77_FUNC(cpocon,CPOCON)
#define CPOSV_F77   F77_FUNC(cposv,CPOSV)
#define CPOEQU_F77  F77_FUNC(cpoequ,CPOEQU)
#define CPORFS_F77  F77_FUNC(cporfs,CPORFS)
#define CPOSVX_F77  F77_FUNC(cposvx,CPOSVX)
#define CGELS_F77   F77_FUNC(cgels,CGELS)
#define CGEEV_F77   F77_FUNC(cgeev,CGEEV)
#define CGEHRD_F77  F77_FUNC(cgehrd,CGEHRD)
#define CHSEQR_F77  F77_FUNC(chseqr,CHSEQR)
#define CTREVC_F77  F77_FUNC(ctrevc,CTREVC)
#define CTREXC_F77  F77_FUNC(ctrexc,CTREXC)
#define CGEES_F77   F77_FUNC(cgees,CGEES)
#define CSTEQR_F77  F77_FUNC(csteqr,CSTEQR)
#define CLARND_F77  F77_FUNC(clarnd,CLARND)
#define CLARNV_F77  F77_FUNC(clarnv,CLARNV)

#endif /* HAVE_TEUCHOS_COMPLEX */

#ifdef __cplusplus
extern "C" {
#endif

// Double precision LAPACK linear solvers
void PREFIX DGELS_F77(Teuchos_fcd ch, const int* m, const int* n, const int* nrhs, double* a, const int* lda, double* b, const int* ldb, double* work, const int* lwork, int* info);
void PREFIX DGEQRF_F77(const int* m, const int* n, double* a, const int* lda, double* tau, double* work, const int* lwork, int* info);
void PREFIX DGETRF_F77(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info); 
void PREFIX DGETRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const double* a, const int* lda,const int* ipiv, double* x , const int* ldx, int* info);
void PREFIX DGETRI_F77(const int* n, double* a, const int* lda, const int* ipiv, double* work , const int* lwork, int* info);
void PREFIX DGECON_F77(Teuchos_fcd norm, const int* n, const double* a, const int* lda, const double* anorm, double* rcond, double* work, int* iwork, int* info); 
void PREFIX DGESV_F77(const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, double* x , const int* ldx, int* info);
void PREFIX DGEEQU_F77(const int* m, const int* n, const double* a, const int* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, int* info); 
void PREFIX DGERFS_F77(Teuchos_fcd, const int* n, const int* nrhs, const double* a, const int* lda, const double* af, const int* ldaf, const int* ipiv, const double* b, const int* ldb, double* x, const int* ldx, double* ferr, double* berr, double* work, int* iwork, int* info);
void PREFIX DGESVX_F77(Teuchos_fcd, Teuchos_fcd, const int* n, const int* nrhs, double* a, const int* lda, double* af, const int* ldaf, int* ipiv, Teuchos_fcd, double* r,
double* c, double* b, const int* ldb, double* x, const int* ldx, double* rcond, double* ferr, double* berr, double* work, int* iwork, int* info);
void PREFIX DSYTRD_F77(Teuchos_fcd, const int* n, double* a, const int* lda, double* D, double* E, double* tau, double* work, const int* lwork, int* info);
void PREFIX DPOTRF_F77(Teuchos_fcd, const int* n, double* a, const int* lda, int* info); 
void PREFIX DPOTRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const double* a, const int* lda, double*x , const int* ldx, int* info);
void PREFIX DPOTRI_F77(Teuchos_fcd, const int* n, double* a, const int* lda, int* info); 
void PREFIX DPOCON_F77(Teuchos_fcd, const int* n, const double* a, const int* lda, const double* anorm, double* rcond, double* work, int* iwork, int* info); 
void PREFIX DPOSV_F77(Teuchos_fcd, const int* n, const int* nrhs, double* a, const int* lda, double*x , const int* ldx, int* info);
void PREFIX DPOEQU_F77(const int* n, const double* a, const int* lda, double* s, double* scond, double* amax, int* info); 
void PREFIX DPORFS_F77(Teuchos_fcd, const int* n, const int* nrhs, double* a, const int* lda, const double* af, const int* ldaf, const double* b, const int* ldb, double* x, const int* ldx, double* ferr, double* berr, double* work, int* iwork, int* info);
void PREFIX DPOSVX_F77(Teuchos_fcd, Teuchos_fcd, const int* n, const int* nrhs, double* a, const int* lda, double* af, const int* ldaf, Teuchos_fcd, double* s, double* b, const int* ldb, double* x, const int* ldx, double* rcond, double* ferr, double* berr, double* work, int* iwork, int* info);

// Single precision LAPACK linear solvers
void PREFIX SGELS_F77(Teuchos_fcd ch, const int* m, const int* n, const int* nrhs, float* a, const int* lda, float* b, const int* ldb, float* work, const int* lwork, int* info);
void PREFIX SGEQRF_F77(const int* m, const int* n, float* a, const int* lda, float* tau, float* work, const int* lwork, int* info);
void PREFIX SGETRF_F77(const int* m, const int* n, float* a, const int* lda, int* ipiv, int* info);
void PREFIX SGETRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const float* a, const int* lda,const int* ipiv, float* x , const int* ldx, int* info);
void PREFIX SGETRI_F77(const int* n, float* a, const int* lda, const int* ipiv, float* work , const int* lwork, int* info);
void PREFIX SGECON_F77(Teuchos_fcd norm, const int* n, const float* a, const int* lda, const float* anorm, float* rcond, float* work, int* iwork, int* info); 
void PREFIX SGESV_F77(const int* n, const int* nrhs, float* a, const int* lda, int* ipiv, float* x , const int* ldx, int* info);
void PREFIX SGEEQU_F77(const int* m, const int* n, const float* a, const int* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, int* info); 
void PREFIX SGERFS_F77(Teuchos_fcd, const int* n, const int* nrhs, const float* a, const int* lda, const float* af, const int* ldaf, const int* ipiv, const float* b, const int* ldb, float* x, const int* ldx, float* ferr, float* berr, float* work, int* iwork, int* info);
void PREFIX SGESVX_F77(Teuchos_fcd, Teuchos_fcd, const int* n, const int* nrhs, float* a, const int* lda, float* af, const int* ldaf, int* ipiv, Teuchos_fcd, float* r,
float* c, float* b, const int* ldb, float* x, const int* ldx, float* rcond, float* ferr, float* berr, float* work, int* iwork, int* info);
void PREFIX SSYTRD_F77(Teuchos_fcd, const int* n, float* a, const int* lda, float* D, float* E, float* tau, float* work, const int* lwork, int* info);
void PREFIX SPOTRF_F77(Teuchos_fcd, const int* n, float* a, const int* lda, int* info); 
void PREFIX SPOTRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const float* a, const int* lda, float*x , const int* ldx, int* info);
void PREFIX SPOTRI_F77(Teuchos_fcd, const int* n, float* a, const int* lda, int* info); 
void PREFIX SPOCON_F77(Teuchos_fcd, const int* n, const float* a, const int* lda, const float* anorm, float* rcond, float* work, int* iwork, int* info); 
void PREFIX SPOSV_F77(Teuchos_fcd, const int* n, const int* nrhs, float* a, const int* lda, float*x , const int* ldx, int* info);
void PREFIX SPOEQU_F77(const int* n, const float* a, const int* lda, float* s, float* scond, float* amax, int* info); 
void PREFIX SPORFS_F77(Teuchos_fcd, const int* n, const int* nrhs, float* a, const int* lda, const float* af, const int* ldaf, const float* b, const int* ldb, float* x, const int* ldx, float* ferr, float* berr, float* work, int* iwork, int* info);
void PREFIX SPOSVX_F77(Teuchos_fcd, Teuchos_fcd, const int* n, const int* nrhs, float* a, const int* lda, float* af, const int* ldaf, Teuchos_fcd, float* s, float* b, const int* ldb, float* x, const int* ldx, float* rcond, float* ferr, float* berr, float* work, int* iwork, int* info);

// Double precision LAPACK eigen solvers
void PREFIX DSPEV_F77(Teuchos_fcd, Teuchos_fcd, const int* n, double* ap, double* w, double* z, const int* ldz, double* work, int* info);
void PREFIX DSYEV_F77(Teuchos_fcd, Teuchos_fcd, const int* n, double* a, const int* lda, double* w, double* work, const int* lwork, int* info);
void PREFIX DSYGV_F77(const int* itype, Teuchos_fcd, Teuchos_fcd, const int* n, double* a, const int* lda, double* B, const int* ldb, double* w, double* work, const int* lwork, int* info);
void PREFIX DSTEQR_F77(Teuchos_fcd, const int* n, double* D, double* E, double* Z, const int* ldz, double* work, int* info);
void PREFIX DGEEV_F77(Teuchos_fcd, Teuchos_fcd, const int* n, double* a, const int* lda, double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr, double* work, const int* lwork, int* info);
void PREFIX DGEHRD_F77(const int* n, const int* ilo, const int* ihi, double* A, const int* lda, double* tau, double* work, const int* lwork, int* info);
void PREFIX DHSEQR_F77(Teuchos_fcd job, Teuchos_fcd, const int* n, const int* ilo, const int* ihi, double* h, const int* ldh, double* wr, double* wi, double* z, const int* ldz, double* work, const int* lwork, int* info);
void PREFIX DGEES_F77(Teuchos_fcd, Teuchos_fcd, int* select, const int* n, double* a, const int* lda, int*sdim, double* wr, double* wi, double* vs, const int* ldvs, double* work, const int* lwork, int* bwork, int* info);
void PREFIX DORGHR_F77(const int* n, const int* ilo, const int* ihi, double* a, const int* lda, double* tau, double* work, int* lwork, int* info);
void PREFIX DORMHR_F77(Teuchos_fcd, Teuchos_fcd, const int* m, const int* n, const int* ilo, const int* ihi, const double* a, const int* lda, const double* tau, double* c, const int* ldc, double* work, int* lwork, int* info);
void PREFIX DORMQR_F77(Teuchos_fcd, Teuchos_fcd, const int* m, const int* n, const int* k, double* a, const int* lda, const double* tau, double* C, const int* ldc, double* work, const int* lwork, int* info);
void PREFIX DTREVC_F77(Teuchos_fcd, Teuchos_fcd, int* select, const int* n, const double* t, const int* ldt, double* vl, const int* ldvl, double* vr, const int* ldvr, const int* mm, int* m, double* work, int* info); 
void PREFIX DTREXC_F77(Teuchos_fcd, const int* n, double* t, const int* ldt, double* q, const int* ldq, int* ifst, int* ilst, double* work, int* info);

// Single precision LAPACK eigen solvers
void PREFIX SSPEV_F77(Teuchos_fcd, Teuchos_fcd, const int* n, float* ap, float* w, float* z, const int* ldz, float* work, int* info);
void PREFIX SSYEV_F77(Teuchos_fcd, Teuchos_fcd, const int* n, float* a, const int* lda, float* w, float* work, const int* lwork, int* info);
void PREFIX SSYGV_F77(const int* itype, Teuchos_fcd, Teuchos_fcd, const int* n, float* a, const int* lda, float* B, const int* ldb, float* w, float* work, const int* lwork, int* info);
void PREFIX SSTEQR_F77(Teuchos_fcd, const int* n, float* D, float* E, float* Z, const int* ldz, float* work, int* info);
void PREFIX SGEEV_F77(Teuchos_fcd, Teuchos_fcd, const int* n, float* a, const int* lda, float* wr, float* wi, float* vl, const int* ldvl, float* vr, const int* ldvr, float* work, const int* lwork, int* info);
void PREFIX SGEHRD_F77(const int* n, const int* ilo, const int* ihi, float* A, const int* lda, float* tau, float* work, const int* lwork, int* info);
void PREFIX SHSEQR_F77(Teuchos_fcd job, Teuchos_fcd, const int* n, const int* ilo, const int* ihi, float* h, const int* ldh, float* wr, float* wi, float* z, const int* ldz, float* work, const int* lwork, int* info);
void PREFIX SGEES_F77(Teuchos_fcd, Teuchos_fcd, int* select, const int* n, float* a, const int* lda, int* sdim, float* wr, float* wi, float* vs, const int* ldvs, float* work, const int* lwork, int* bwork, int* info);
void PREFIX SORGHR_F77(const int* n, const int* ilo, const int* ihi, float* a, const int* lda, float* tau, float* work, int* lwork, int* info);
void PREFIX SORMHR_F77(Teuchos_fcd, Teuchos_fcd, const int* m, const int* n, const int* ilo, const int* ihi, const float* a, const int* lda, const float* tau, float* c, const int* ldc, float* work, int* lwork, int* info);
void PREFIX SORMQR_F77(Teuchos_fcd, Teuchos_fcd, const int* m, const int* n, const int* k, float* a, const int* lda, const float* tau, float* C, const int* ldc, float* work, const int* lwork, int* info);
void PREFIX STREVC_F77(Teuchos_fcd, Teuchos_fcd, int* select, const int* n, const float* t, const int* ldt, float* vl, const int* ldvl, float* vr, const int* ldvr, const int* mm, int* m, float* work, int* info); 
void PREFIX STREXC_F77(Teuchos_fcd, const int* n, float* t, const int* ldt, float* q, const int* ldq, int* ifst, int* ilst, float* work, int* info);

float PREFIX SLARND_F77(const int* idist, int* seed);
double PREFIX DLARND_F77(const int* idist, int* seed);

void PREFIX SLARNV_F77(const int* idist, int* seed, const int* n, float* v);
void PREFIX DLARNV_F77(const int* idist, int* seed, const int* n, double* v);

float PREFIX SLAMCH_F77(Teuchos_fcd);
double PREFIX DLAMCH_F77(Teuchos_fcd);

float PREFIX SLAPY2_F77(const float* x, const float* y);
double PREFIX DLAPY2_F77(const double* x, const double* y);

#ifdef HAVE_TEUCHOS_COMPLEX

// Double precision complex LAPACK linear solvers
void PREFIX ZGELS_F77(Teuchos_fcd ch, const int* m, const int* n, const int* nrhs, complex<double>* a, const int* lda, complex<double>* b, const int* ldb, complex<double>* work, const int* lwork, int* info);
void PREFIX ZGEQRF_F77(const int* m, const int* n, complex<double>* a, const int* lda, complex<double>* tau, complex<double>* work, const int* lwork, int* info);
void PREFIX ZGETRF_F77(const int* m, const int* n, complex<double>* a, const int* lda, int* ipiv, int* info); 
void PREFIX ZGETRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const complex<double>* a, const int* lda,const int* ipiv, complex<double>* x , const int* ldx, int* info);
void PREFIX ZGETRI_F77(const int* n, complex<double>* a, const int* lda, const int* ipiv, complex<double>* work , const int* lwork, int* info);
void PREFIX ZGECON_F77(Teuchos_fcd norm, const int* n, const complex<double>* a, const int* lda, const double* anorm, double* rcond, complex<double>* work, double* rwork, int* info); 
void PREFIX ZGESV_F77(const int* n, const int* nrhs, complex<double>* a, const int* lda, int* ipiv, complex<double>* x , const int* ldx, int* info);
void PREFIX ZGEEQU_F77(const int* m, const int* n, const complex<double>* a, const int* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, int* info); 
void PREFIX ZGERFS_F77(Teuchos_fcd, const int* n, const int* nrhs, const complex<double>* a, const int* lda, const complex<double>* af, const int* ldaf, const int* ipiv, const complex<double>* b, const int* ldb, complex<double>* x, const int* ldx, double* ferr, double* berr, complex<double>* work, double* iwork, int* info);
void PREFIX ZGESVX_F77(Teuchos_fcd, Teuchos_fcd, const int* n, const int* nrhs, complex<double>* a, const int* lda, complex<double>* af, const int* ldaf, int* ipiv, Teuchos_fcd, double* r,
double* c, complex<double>* b, const int* ldb, complex<double>* x, const int* ldx, double* rcond, double* ferr, double* berr, complex<double>* work, double* iwork, int* info);
void PREFIX ZPOTRF_F77(Teuchos_fcd, const int* n, complex<double>* a, const int* lda, int* info); 
void PREFIX ZPOTRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const complex<double>* a, const int* lda, complex<double>*x , const int* ldx, int* info);
void PREFIX ZPOTRI_F77(Teuchos_fcd, const int* n, complex<double>* a, const int* lda, int* info); 
void PREFIX ZPOCON_F77(Teuchos_fcd, const int* n, const complex<double>* a, const int* lda, const double* anorm, double* rcond, complex<double>* work, double* rwork, int* info); 
void PREFIX ZPOSV_F77(Teuchos_fcd, const int* n, const int* nrhs, complex<double>* a, const int* lda, complex<double>*x , const int* ldx, int* info);
void PREFIX ZPOEQU_F77(const int* n, const complex<double>* a, const int* lda, double* s, double* scond, double* amax, int* info); 
void PREFIX ZPORFS_F77(Teuchos_fcd, const int* n, const int* nrhs, complex<double>* a, const int* lda, const complex<double>* af, const int* ldaf, const complex<double>* b, const int* ldb, complex<double>* x, const int* ldx, double* ferr, double* berr, complex<double>* work, double* rwork, int* info);
void PREFIX ZPOSVX_F77(Teuchos_fcd, Teuchos_fcd, const int* n, const int* nrhs, complex<double>* a, const int* lda, complex<double>* af, const int* ldaf, Teuchos_fcd, double* s, complex<double>* b, const int* ldb, complex<double>* x, const int* ldx, double* rcond, double* ferr, double* berr, complex<double>* work, double* rwork, int* info);

// Single precision complex LAPACK linear solvers
void PREFIX CGELS_F77(Teuchos_fcd ch, const int* m, const int* n, const int* nrhs, complex<float>* a, const int* lda, complex<float>* b, const int* ldb, complex<float>* work, const int* lwork, int* info);
void PREFIX CGEQRF_F77(const int* m, const int* n, complex<float>* a, const int* lda, complex<float>* tau, complex<float>* work, const int* lwork, int* info);
void PREFIX CGETRF_F77(const int* m, const int* n, complex<float>* a, const int* lda, int* ipiv, int* info);
void PREFIX CGETRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const complex<float>* a, const int* lda,const int* ipiv, complex<float>* x , const int* ldx, int* info);
void PREFIX CGETRI_F77(const int* n, complex<float>* a, const int* lda, const int* ipiv, complex<float>* work , const int* lwork, int* info);
void PREFIX CGECON_F77(Teuchos_fcd norm, const int* n, const complex<float>* a, const int* lda, const float* anorm, float* rcond, complex<float>* work, float* rwork, int* info); 
void PREFIX CGESV_F77(const int* n, const int* nrhs, complex<float>* a, const int* lda, int* ipiv, complex<float>* x, const int* ldx, int* info);
void PREFIX CGEEQU_F77(const int* m, const int* n, const complex<float>* a, const int* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, int* info); 
void PREFIX CGERFS_F77(Teuchos_fcd, const int* n, const int* nrhs, const complex<float>* a, const int* lda, const complex<float>* af, const int* ldaf, const int* ipiv, const complex<float>* b, const int* ldb, complex<float>* x, const int* ldx, float* ferr, float* berr, complex<float>* work, float* rwork, int* info);
void PREFIX CGESVX_F77(Teuchos_fcd, Teuchos_fcd, const int* n, const int* nrhs, complex<float>* a, const int* lda, complex<float>* af, const int* ldaf, int* ipiv, Teuchos_fcd, float* r,
float* c, complex<float>* b, const int* ldb, complex<float>* x, const int* ldx, float* rcond, float* ferr, float* berr, complex<float>* work, float* rwork, int* info);
void PREFIX CPOTRF_F77(Teuchos_fcd, const int* n, complex<float>* a, const int* lda, int* info); 
void PREFIX CPOTRS_F77(Teuchos_fcd, const int* n, const int* nrhs, const complex<float>* a, const int* lda, complex<float>*x , const int* ldx, int* info);
void PREFIX CPOTRI_F77(Teuchos_fcd, const int* n, complex<float>* a, const int* lda, int* info); 
void PREFIX CPOCON_F77(Teuchos_fcd, const int* n, const complex<float>* a, const int* lda, const float* anorm, float* rcond, complex<float>* work, float* rwork, int* info); 
void PREFIX CPOSV_F77(Teuchos_fcd, const int* n, const int* nrhs, complex<float>* a, const int* lda, complex<float>*x , const int* ldx, int* info);
void PREFIX CPOEQU_F77(const int* n, const complex<float>* a, const int* lda, float* s, float* scond, float* amax, int* info); 
void PREFIX CPORFS_F77(Teuchos_fcd, const int* n, const int* nrhs, complex<float>* a, const int* lda, const complex<float>* af, const int* ldaf, const complex<float>* b, const int* ldb, complex<float>* x, const int* ldx, float* ferr, float* berr, complex<float>* work, float* rwork, int* info);
void PREFIX CPOSVX_F77(Teuchos_fcd, Teuchos_fcd, const int* n, const int* nrhs, complex<float>* a, const int* lda, complex<float>* af, const int* ldaf, Teuchos_fcd, float* s, complex<float>* b, const int* ldb, complex<float>* x, const int* ldx, float* rcond, float* ferr, float* berr, complex<float>* work, float* rwork, int* info);

// Double precision complex LAPACK eigen solvers
void PREFIX ZSTEQR_F77(Teuchos_fcd, const int* n, double* D, double* E, complex<double>* Z, const int* ldz, complex<double>* work, int* info);
void PREFIX ZGEEV_F77(Teuchos_fcd, Teuchos_fcd, const int* n, complex<double>* a, const int* lda, complex<double>* w, complex<double>* vl, const int* ldvl, complex<double>* vr, const int* ldvr, complex<double>* work, const int* lwork, double* rwork, int* info);
void PREFIX ZGEHRD_F77(const int* n, const int* ilo, const int* ihi, complex<double>* A, const int* lda, complex<double>* tau, complex<double>* work, const int* lwork, int* info);
void PREFIX ZHSEQR_F77(Teuchos_fcd job, Teuchos_fcd, const int* n, const int* ilo, const int* ihi, complex<double>* h, const int* ldh, complex<double>* w, complex<double>* z, const int* ldz, complex<double>* work, const int* lwork, int* info);
void PREFIX ZGEES_F77(Teuchos_fcd, Teuchos_fcd, int* select, const int* n, complex<double>* a, const int* lda, int* sdim, complex<double>* w, complex<double>* vs, const int* ldvs, complex<double>* work, const int* lwork, double* rwork, int* bwork, int* info);
void PREFIX ZTREVC_F77(Teuchos_fcd, Teuchos_fcd, int* select, const int* n, const complex<double>* t, const int* ldt, complex<double>* vl, const int* ldvl, complex<double>* vr, const int* ldvr, const int* mm, int* m, complex<double>* work, double* rwork, int* info); 
void PREFIX ZTREXC_F77(Teuchos_fcd, const int* n, complex<double>* t, const int* ldt, complex<double>* q, const int* ldq, int* ifst, int* ilst, int* info);

// Single precision complex LAPACK eigen solvers
void PREFIX CSTEQR_F77(Teuchos_fcd, const int* n, complex<float>* D, complex<float>* E, complex<float>* Z, const int* ldz, complex<float>* work, int* info);
void PREFIX CGEEV_F77(Teuchos_fcd, Teuchos_fcd, const int* n, complex<float>* a, const int* lda, complex<float>* wr, complex<float>* vl, const int* ldvl, complex<float>* vr, const int* ldvr, complex<float>* work, const int* lwork, float* rwork, int* info);
void PREFIX CGEHRD_F77(const int* n, const int* ilo, const int* ihi, complex<float>* A, const int* lda, complex<float>* tau, complex<float>* work, const int* lwork, int* info);
void PREFIX CHSEQR_F77(Teuchos_fcd job, Teuchos_fcd, const int* n, const int* ilo, const int* ihi, complex<float>* h, const int* ldh, complex<float>* w, complex<float>* z, const int* ldz, complex<float>* work, const int* lwork, int* info);
void PREFIX CGEES_F77(Teuchos_fcd, Teuchos_fcd, int* select, const int* n, complex<float>* a, const int* lda, int* sdim, complex<float>* w, complex<float>* vs, const int* ldvs, complex<float>* work, const int* lwork, float* rwork, int* bwork, int* info);
void PREFIX CTREVC_F77(Teuchos_fcd, Teuchos_fcd, int* select, const int* n, const complex<float>* t, const int* ldt, complex<float>* vl, const int* ldvl, complex<float>* vr, const int* ldvr, const int* mm, int* m, complex<float>* work, float* rwork, int* info); 
void PREFIX CTREXC_F77(Teuchos_fcd, const int* n, complex<float>* t, const int* ldt, complex<float>* q, const int* ldq, int* ifst, int* ilst, int* info);

complex<float> PREFIX CLARND_F77(const int* idist, int* seed);
complex<double> PREFIX ZLARND_F77(const int* idist, int* seed);

void PREFIX CLARNV_F77(const int* idist, int* seed, const int* n, complex<float>* v);
void PREFIX ZLARNV_F77(const int* idist, int* seed, const int* n, complex<double>* v);

#endif /* HAVE_TEUCHOS_COMPLEX */

#ifdef __cplusplus
}

#endif

#endif // end of TEUCHOS_LAPACK_WRAPPERS_HPP_
