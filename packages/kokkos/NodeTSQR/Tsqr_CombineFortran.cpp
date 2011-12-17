//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <Tsqr_ConfigDefs.hpp>
#include <Tsqr_CombineFortran.hpp>
#include <Tsqr_FortranCInterface.hpp>

#include <complex>

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

extern "C" void FortranCInterface_MODULE_(tsqrcombine,z_apply_inner, TSQRCOMBINE,Z_APPLY_INNER)
  (const char* const trans, 
   const int m, 
   const int ncols_C, 
   const int ncols_Q, 
   const std::complex<double> A[], 
   const int lda, 
   const std::complex<double> tau[],
   std::complex<double> C_top[], 
   const int ldc_top, 
   std::complex<double> C_bot[], 
   const int ldc_bot, 
   std::complex<double> work[]);

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

extern "C" void FortranCInterface_MODULE_(tsqrcombine,c_apply_inner, TSQRCOMBINE,C_APPLY_INNER)
  (const char* const trans, 
   const int m, 
   const int ncols_C, 
   const int ncols_Q, 
   const std::complex<float> A[], 
   const int lda, 
   const std::complex<float> tau[],
   std::complex<float> C_top[], 
   const int ldc_top, 
   std::complex<float> C_bot[], 
   const int ldc_bot, 
   std::complex<float> work[]);

////////////////////////////////////////////////////////////////////////////////

extern "C" void FortranCInterface_MODULE_(tsqrcombine,d_factor_inner, TSQRCOMBINE,D_FACTOR_INNER)
  (const int m,
   const int n,
   double R[],
   const int ldr,
   double A[],
   const int lda, 
   double tau[],
   double work[]);

extern "C" void FortranCInterface_MODULE_(tsqrcombine,z_factor_inner, TSQRCOMBINE,Z_FACTOR_INNER)
  (const int m,
   const int n,
   std::complex<double> R[],
   const int ldr,
   std::complex<double> A[],
   const int lda, 
   std::complex<double> tau[],
   std::complex<double> work[]);

extern "C" void FortranCInterface_MODULE_(tsqrcombine,s_factor_inner, TSQRCOMBINE,S_FACTOR_INNER)
  (const int m,
   const int n,
   float R[],
   const int ldr,
   float A[],
   const int lda, 
   float tau[],
   float work[]);

extern "C" void FortranCInterface_MODULE_(tsqrcombine,c_factor_inner, TSQRCOMBINE,C_FACTOR_INNER)
  (const int m,
   const int n,
   std::complex<float> R[],
   const int ldr,
   std::complex<float> A[],
   const int lda, 
   std::complex<float> tau[],
   std::complex<float> work[]);

////////////////////////////////////////////////////////////////////////////////

extern "C" void FortranCInterface_MODULE_(tsqrcombine,d_factor_pair, TSQRCOMBINE,D_FACTOR_PAIR)
  (const int n, 
   double R_top[], 
   const int ldr_top, 
   double R_bot[], 
   const int ldr_bot, 
   double tau[], 
   double work[]);

extern "C" void FortranCInterface_MODULE_(tsqrcombine,z_factor_pair, TSQRCOMBINE,Z_FACTOR_PAIR)
  (const int n, 
   std::complex<double> R_top[], 
   const int ldr_top, 
   std::complex<double> R_bot[], 
   const int ldr_bot, 
   std::complex<double> tau[], 
   std::complex<double> work[]);

extern "C" void FortranCInterface_MODULE_(tsqrcombine,s_factor_pair, TSQRCOMBINE,S_FACTOR_PAIR)
  (const int n, 
   float R_top[], 
   const int ldr_top, 
   float R_bot[], 
   const int ldr_bot, 
   float tau[], 
   float work[]);

extern "C" void FortranCInterface_MODULE_(tsqrcombine,c_factor_pair, TSQRCOMBINE,C_FACTOR_PAIR)
  (const int n, 
   std::complex<float> R_top[], 
   const int ldr_top, 
   std::complex<float> R_bot[], 
   const int ldr_bot, 
   std::complex<float> tau[], 
   std::complex<float> work[]);

////////////////////////////////////////////////////////////////////////////////

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

extern "C" void FortranCInterface_MODULE_(tsqrcombine,z_apply_pair, TSQRCOMBINE,Z_APPLY_PAIR)
(const char* const trans, 
 const int ncols_C, 
 const int ncols_Q, 
 const std::complex<double> R_bot[], 
 const int ldr_bot,
 const std::complex<double> tau[], 
 std::complex<double> C_top[], 
 const int ldc_top, 
 std::complex<double> C_bot[], 
 const int ldc_bot, 
 std::complex<double> work[]);

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

extern "C" void FortranCInterface_MODULE_(tsqrcombine,c_apply_pair, TSQRCOMBINE,C_APPLY_PAIR)
(const char* const trans, 
 const int ncols_C, 
 const int ncols_Q, 
 const std::complex<float> R_bot[], 
 const int ldr_bot,
 const std::complex<float> tau[], 
 std::complex<float> C_top[], 
 const int ldc_top, 
 std::complex<float> C_bot[], 
 const int ldc_bot, 
 std::complex<float> work[]);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  template<>
  void
  CombineFortran< double, false >::
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
  CombineFortran< double, false >::
  apply_inner (const ApplyType& apply_type,
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
      (apply_type.toString().c_str(), m, ncols_C, ncols_Q, 
       A, lda, tau, C_top, ldc_top, C_bot, ldc_bot, work);
  }

  template<>
  void
  CombineFortran< double, false >::
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
  CombineFortran< double, false >::
  apply_pair (const ApplyType& apply_type,
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
      (apply_type.toString().c_str(), ncols_C, ncols_Q, 
       R_bot, ldr_bot, tau, C_top, ldc_top, C_bot, ldc_bot, work);
  }

  //////////////////////////////////////////////////////////////////////////////  
  //////////////////////////////////////////////////////////////////////////////  

  template<>
  void
  CombineFortran< float, false >::
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
  CombineFortran< float, false >::
  apply_inner (const ApplyType& apply_type,
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
      (apply_type.toString().c_str(), m, ncols_C, ncols_Q, 
       A, lda, tau, C_top, ldc_top, C_bot, ldc_bot, work);
  }

  template<>
  void
  CombineFortran< float, false >::
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
  CombineFortran< float, false >::
  apply_pair (const ApplyType& apply_type,
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
      (apply_type.toString().c_str(), ncols_C, ncols_Q, 
       R_bot, ldr_bot, tau, C_top, ldc_top, C_bot, ldc_bot, work);
  }
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  template<>
  void
  CombineFortran< std::complex<double>, true >::
  factor_inner (const int m,
  		const int n,
  		std::complex<double> R[],
  		const int ldr,
  		std::complex<double> A[],
  		const int lda, 
  		std::complex<double> tau[],
  		std::complex<double> work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,z_factor_inner, TSQRCOMBINE,Z_FACTOR_INNER) 
      (m, n, R, ldr, A, lda, tau, work);
  }

  template<>
  void
  CombineFortran< std::complex<double>, true >::
  apply_inner (const ApplyType& apply_type,
  	       const int m, 
  	       const int ncols_C, 
  	       const int ncols_Q, 
  	       const std::complex<double> A[], 
  	       const int lda, 
  	       const std::complex<double> tau[],
  	       std::complex<double> C_top[], 
  	       const int ldc_top, 
  	       std::complex<double> C_bot[], 
  	       const int ldc_bot, 
  	       std::complex<double> work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,z_apply_inner, TSQRCOMBINE,Z_APPLY_INNER) 
      (apply_type.toString().c_str(), m, ncols_C, ncols_Q, 
       A, lda, tau, C_top, ldc_top, C_bot, ldc_bot, work);
  }

  template<>
  void
  CombineFortran< std::complex<double>, true >::
  factor_pair (const int n, 
  	       std::complex<double> R_top[], 
  	       const int ldr_top, 
  	       std::complex<double> R_bot[], 
  	       const int ldr_bot, 
  	       std::complex<double> tau[], 
  	       std::complex<double> work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,z_factor_pair, TSQRCOMBINE,Z_FACTOR_PAIR) 
      (n, R_top, ldr_top, R_bot, ldr_bot, tau, work);
  }

  template<>
  void
  CombineFortran< std::complex<double>, true >::
  apply_pair (const ApplyType& apply_type,
  	      const int ncols_C, 
  	      const int ncols_Q, 
  	      const std::complex<double> R_bot[], 
  	      const int ldr_bot,
  	      const std::complex<double> tau[], 
  	      std::complex<double> C_top[], 
  	      const int ldc_top, 
  	      std::complex<double> C_bot[], 
  	      const int ldc_bot, 
  	      std::complex<double> work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,z_apply_pair, TSQRCOMBINE,Z_APPLY_PAIR)
      (apply_type.toString().c_str(), ncols_C, ncols_Q, 
       R_bot, ldr_bot, tau, C_top, ldc_top, C_bot, ldc_bot, work);
  }

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  template<>
  void
  CombineFortran< std::complex<float>, true >::
  factor_inner (const int m,
  		const int n,
  		std::complex<float> R[],
  		const int ldr,
  		std::complex<float> A[],
  		const int lda, 
  		std::complex<float> tau[],
  		std::complex<float> work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,c_factor_inner, TSQRCOMBINE,C_FACTOR_INNER) 
      (m, n, R, ldr, A, lda, tau, work);
  }

  template<>
  void
  CombineFortran< std::complex<float>, true >::
  apply_inner (const ApplyType& apply_type,
  	       const int m, 
  	       const int ncols_C, 
  	       const int ncols_Q, 
  	       const std::complex<float> A[], 
  	       const int lda, 
  	       const std::complex<float> tau[],
  	       std::complex<float> C_top[], 
  	       const int ldc_top, 
  	       std::complex<float> C_bot[], 
  	       const int ldc_bot, 
  	       std::complex<float> work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,c_apply_inner, TSQRCOMBINE,C_APPLY_INNER) 
      (apply_type.toString().c_str(), m, ncols_C, ncols_Q, 
       A, lda, tau, C_top, ldc_top, C_bot, ldc_bot, work);
  }

  template<>
  void
  CombineFortran< std::complex<float>, true >::
  factor_pair (const int n, 
  	       std::complex<float> R_top[], 
  	       const int ldr_top, 
  	       std::complex<float> R_bot[], 
  	       const int ldr_bot, 
  	       std::complex<float> tau[], 
  	       std::complex<float> work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,c_factor_pair, TSQRCOMBINE,C_FACTOR_PAIR) 
      (n, R_top, ldr_top, R_bot, ldr_bot, tau, work);
  }

  template<>
  void
  CombineFortran< std::complex<float>, true >::
  apply_pair (const ApplyType& apply_type,
  	      const int ncols_C, 
  	      const int ncols_Q, 
  	      const std::complex<float> R_bot[], 
  	      const int ldr_bot,
  	      const std::complex<float> tau[], 
  	      std::complex<float> C_top[], 
  	      const int ldc_top, 
  	      std::complex<float> C_bot[], 
  	      const int ldc_bot, 
  	      std::complex<float> work[]) const
  {
    FortranCInterface_MODULE_(tsqrcombine,c_apply_pair, TSQRCOMBINE,C_APPLY_PAIR)
      (apply_type.toString().c_str(), ncols_C, ncols_Q, 
       R_bot, ldr_bot, tau, C_top, ldc_top, C_bot, ldc_bot, work);
  }

} // namespace TSQR
