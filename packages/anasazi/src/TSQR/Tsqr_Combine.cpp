// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
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

#include <Tsqr_Combine.hpp>
#include <Tsqr_FortranCInterface.hpp>
#include <Tsqr_Config.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

extern "C" void FortranCInterface_MODULE_(tsqrcombine,d_apply_inner, TSQRCOMBINE,D_APPLY_INNER)
  (const char* const trans, 
   const int m, 
   const int ncols_C, 
   const int ncols_Q, 
   const double A[], 
   const int lda, 
   const double tau[],
   double C_top[], 
   const int ldc_top, 
   double C_bot[], 
   const int ldc_bot, 
   double work[]);

extern "C" void FortranCInterface_MODULE_(tsqrcombine,s_apply_inner, TSQRCOMBINE,S_APPLY_INNER)
  (const char* const trans, 
   const int m, 
   const int ncols_C, 
   const int ncols_Q, 
   const float A[], 
   const int lda, 
   const float tau[],
   float C_top[], 
   const int ldc_top, 
   float C_bot[], 
   const int ldc_bot, 
   float work[]);

extern "C" void FortranCInterface_MODULE_(tsqrcombine,d_factor_inner, TSQRCOMBINE,D_FACTOR_INNER)
  (const int m,
   const int n,
   double R[],
   const int ldr,
   double A[],
   const int lda, 
   double tau[],
   double work[]);

extern "C" void FortranCInterface_MODULE_(tsqrcombine,s_factor_inner, TSQRCOMBINE,S_FACTOR_INNER)
  (const int m,
   const int n,
   float R[],
   const int ldr,
   float A[],
   const int lda, 
   float tau[],
   float work[]);

extern "C" void FortranCInterface_MODULE_(tsqrcombine,d_factor_pair, TSQRCOMBINE,D_FACTOR_PAIR)
  (const int n, 
   double R_top[], 
   const int ldr_top, 
   double R_bot[], 
   const int ldr_bot, 
   double tau[], 
   double work[]);

extern "C" void FortranCInterface_MODULE_(tsqrcombine,s_factor_pair, TSQRCOMBINE,S_FACTOR_PAIR)
  (const int n, 
   float R_top[], 
   const int ldr_top, 
   float R_bot[], 
   const int ldr_bot, 
   float tau[], 
   float work[]);

extern "C" void FortranCInterface_MODULE_(tsqrcombine,d_apply_pair, TSQRCOMBINE,D_APPLY_PAIR)
(const char* const trans, 
 const int ncols_C, 
 const int ncols_Q, 
 const double R_bot[], 
 const int ldr_bot,
 const double tau[], 
 double C_top[], 
 const int ldc_top, 
 double C_bot[], 
 const int ldc_bot, 
 double work[]);

extern "C" void FortranCInterface_MODULE_(tsqrcombine,s_apply_pair, TSQRCOMBINE,S_APPLY_PAIR)
(const char* const trans, 
 const int ncols_C, 
 const int ncols_Q, 
 const float R_bot[], 
 const int ldr_bot,
 const float tau[], 
 float C_top[], 
 const int ldc_top, 
 float C_bot[], 
 const int ldc_bot, 
 float work[]);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  template<>
  void
  Combine<int, double, false >::
  factor_inner (const int m,
		const int n,
		double R[],
		const int ldr,
		double A[],
		const int lda, 
		double tau[],
		double work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,d_factor_inner, TSQRCOMBINE,D_FACTOR_INNER) 
      (m, n, R, ldr, A, lda, tau, work);
  }

  template<>
  void
  Combine<int, double, false >::
  apply_inner (const char* const trans, 
	       const int m, 
	       const int ncols_C, 
	       const int ncols_Q, 
	       const double A[], 
	       const int lda, 
	       const double tau[],
	       double C_top[], 
	       const int ldc_top, 
	       double C_bot[], 
	       const int ldc_bot, 
	       double work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,d_apply_inner, TSQRCOMBINE,D_APPLY_INNER) 
      (trans, m, ncols_C, ncols_Q, A, lda, tau, C_top, ldc_top, C_bot, ldc_bot, work);
  }

  template<>
  void
  Combine<int, double, false >::
  factor_pair (const int n, 
	       double R_top[], 
	       const int ldr_top, 
	       double R_bot[], 
	       const int ldr_bot, 
	       double tau[], 
	       double work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,d_factor_pair, TSQRCOMBINE,D_FACTOR_PAIR) 
      (n, R_top, ldr_top, R_bot, ldr_bot, tau, work);
  }

  template<>
  void
  Combine<int, double, false >::
  apply_pair (const char* const trans, 
	      const int ncols_C, 
	      const int ncols_Q, 
	      const double R_bot[], 
	      const int ldr_bot,
	      const double tau[], 
	      double C_top[], 
	      const int ldc_top, 
	      double C_bot[], 
	      const int ldc_bot, 
	      double work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,d_apply_pair, TSQRCOMBINE,D_APPLY_PAIR)
      (trans, ncols_C, ncols_Q, R_bot, ldr_bot, tau, C_top, ldc_top, C_bot, ldc_bot, work);
  }

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  template<>
  void
  Combine<int, float, false >::
  factor_inner (const int m,
		const int n,
		float R[],
		const int ldr,
		float A[],
		const int lda, 
		float tau[],
		float work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,s_factor_inner, TSQRCOMBINE,S_FACTOR_INNER) 
      (m, n, R, ldr, A, lda, tau, work);
  }

  template<>
  void
  Combine<int, float, false >::
  apply_inner (const char* const trans, 
	       const int m, 
	       const int ncols_C, 
	       const int ncols_Q, 
	       const float A[], 
	       const int lda, 
	       const float tau[],
	       float C_top[], 
	       const int ldc_top, 
	       float C_bot[], 
	       const int ldc_bot, 
	       float work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,s_apply_inner, TSQRCOMBINE,S_APPLY_INNER) 
      (trans, m, ncols_C, ncols_Q, A, lda, tau, C_top, ldc_top, C_bot, ldc_bot, work);
  }

  template<>
  void
  Combine<int, float, false >::
  factor_pair (const int n, 
	       float R_top[], 
	       const int ldr_top, 
	       float R_bot[], 
	       const int ldr_bot, 
	       float tau[], 
	       float work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,s_factor_pair, TSQRCOMBINE,S_FACTOR_PAIR) 
      (n, R_top, ldr_top, R_bot, ldr_bot, tau, work);
  }

  template<>
  void
  Combine<int, float, false >::
  apply_pair (const char* const trans, 
	      const int ncols_C, 
	      const int ncols_Q, 
	      const float R_bot[], 
	      const int ldr_bot,
	      const float tau[], 
	      float C_top[], 
	      const int ldc_top, 
	      float C_bot[], 
	      const int ldc_bot, 
	      float work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,s_apply_pair, TSQRCOMBINE,S_APPLY_PAIR)
      (trans, ncols_C, ncols_Q, R_bot, ldr_bot, tau, C_top, ldc_top, C_bot, ldc_bot, work);
  }

} // namespace TSQR
