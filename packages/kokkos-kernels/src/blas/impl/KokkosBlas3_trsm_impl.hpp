/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSBLAS3_TRSM_IMPL_HPP_
#define KOKKOSBLAS3_TRSM_IMPL_HPP_

/// \file KokkosBlas3_trsm_impl.hpp
/// \brief Implementation(s) of triangular linear system solve (with multiple RHSs)
/// \brief Sequential fall-back implementation calls the exisiting serial batched TRSM.
/// \brief Two sequential fall-back implementations for conjugate transpose case are
/// \brief also based on the exisiting serial batched TRSM.

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "KokkosBatched_Trsm_Decl.hpp"
#include "KokkosBatched_Trsm_Serial_Impl.hpp"

using namespace KokkosBatched;

namespace KokkosBlas {
namespace Impl {

template<typename ScalarType,
         typename ValueType>
int
SerialTrsmInternalLeftLowerConj(const bool use_unit_diag,
                                const int m, const int n,
                                const ScalarType alpha,
                                const ValueType* KOKKOS_RESTRICT A, const int as0, const int as1,
                                /**/  ValueType* KOKKOS_RESTRICT B, const int bs0, const int bs1) {

  typedef Kokkos::Details::ArithTraits<ValueType> AT;
  
  const ScalarType one(1.0), zero(0.0);
      
  if (alpha == zero)   SerialSetInternal  ::invoke(m, n, zero,  B, bs0, bs1);
  else {
    if (alpha != one)  SerialScaleInternal::invoke(m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;

    for (int p=0;p<m;++p) {
      const int iend = m-p-1, jend = n;
        
      const ValueType
        *KOKKOS_RESTRICT a21 = iend ? A+(p+1)*as0+p*as1 : NULL;
          
      ValueType
        *KOKKOS_RESTRICT b1t =        B+p*bs0,
        *KOKKOS_RESTRICT B2  = iend ? B+(p+1)*bs0 : NULL;
        
      if (!use_unit_diag) {
        const ValueType alpha11 = AT::conj(A[p*as0+p*as1]);
        for (int j=0;j<jend;++j)
          b1t[j*bs1] = b1t[j*bs1] / alpha11;
      }
        
      for (int i=0;i<iend;++i)
        for (int j=0;j<jend;++j)
          B2[i*bs0+j*bs1] -= AT::conj(a21[i*as0]) * b1t[j*bs1];
    }
  }      
  return 0;
}

template<typename ScalarType,
         typename ValueType>
int
SerialTrsmInternalLeftUpperConj(const bool use_unit_diag,
                                const int m, const int n,
                                const ScalarType alpha,
                                const ValueType* KOKKOS_RESTRICT A, const int as0, const int as1,
                                /**/  ValueType* KOKKOS_RESTRICT B, const int bs0, const int bs1) {

  typedef Kokkos::Details::ArithTraits<ValueType> AT;

  const ScalarType one(1.0), zero(0.0);

  if (alpha == zero)  SerialSetInternal  ::invoke(m, n, zero,  B, bs0, bs1);
  else {
    if (alpha != one) SerialScaleInternal::invoke(m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;
      
    ValueType *KOKKOS_RESTRICT B0 = B;
    for (int p=(m-1);p>=0;--p) {
      const int iend = p, jend = n;

      const ValueType* KOKKOS_RESTRICT a01 = A+p*as1;
      ValueType* KOKKOS_RESTRICT b1t = B+p*bs0;

      if (!use_unit_diag) {
        const ValueType alpha11 = AT::conj(A[p*as0+p*as1]);
        for (int j=0;j<n;++j)
          b1t[j*bs1] = b1t[j*bs1] / alpha11;
      }
        
      if (p>0){//Note: A workaround to produce correct results for complex<double> with Intel-18.2.199
        for (int i=0;i<iend;++i)
          for (int j=0;j<jend;++j)
            B0[i*bs0+j*bs1] -= AT::conj(a01[i*as0]) * b1t[j*bs1];
      }
    }
  }
  return 0;
}


template<class AViewType, class BViewType>
void SerialTrsm_Invoke (const char side[],
                  const char uplo[],
                  const char trans[],
                  const char diag[],
                  typename BViewType::const_value_type& alpha,
                  const AViewType& A,
                  const BViewType& B)
{
  //Side::Left, Uplo::Lower, Trans::NoTranspose
  if (((side[0]=='L')||(side[0]=='l'))&&((uplo[0]=='L')||(uplo[0]=='l'))&&((trans[0]=='N')||(trans[0]=='n'))&&((diag[0]=='U')||(diag[0]=='u')))
    SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(Diag::Unit::use_unit_diag,
                                                               B.extent(0), B.extent(1),
                                                               alpha,
                                                               A.data(), A.stride(0), A.stride(1),
                                                               B.data(), B.stride(0), B.stride(1));
  if (((side[0]=='L')||(side[0]=='l'))&&((uplo[0]=='L')||(uplo[0]=='l'))&&((trans[0]=='N')||(trans[0]=='n'))&&((diag[0]=='N')||(diag[0]=='n')))
    SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(Diag::NonUnit::use_unit_diag,
                                                               B.extent(0), B.extent(1),
                                                               alpha,
                                                               A.data(), A.stride(0), A.stride(1),
                                                               B.data(), B.stride(0), B.stride(1));

  //Side::Left, Uplo::Lower, Trans::Transpose
  if (((side[0]=='L')||(side[0]=='l'))&&((uplo[0]=='L')||(uplo[0]=='l'))&&((trans[0]=='T')||(trans[0]=='t'))&&((diag[0]=='U')||(diag[0]=='u')))
    SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(Diag::Unit::use_unit_diag,
                                                               B.extent(0), B.extent(1),
                                                               alpha,
                                                               A.data(), A.stride(1), A.stride(0),
                                                               B.data(), B.stride(0), B.stride(1));
  if (((side[0]=='L')||(side[0]=='l'))&&((uplo[0]=='L')||(uplo[0]=='l'))&&((trans[0]=='T')||(trans[0]=='t'))&&((diag[0]=='N')||(diag[0]=='n')))
    SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(Diag::NonUnit::use_unit_diag,
                                                               B.extent(0), B.extent(1),
                                                               alpha,
                                                               A.data(), A.stride(1), A.stride(0),
                                                               B.data(), B.stride(0), B.stride(1));

  //Side::Left, Uplo::Lower, Trans::ConjTranspose
  if (((side[0]=='L')||(side[0]=='l'))&&((uplo[0]=='L')||(uplo[0]=='l'))&&((trans[0]=='C')||(trans[0]=='c'))&&((diag[0]=='U')||(diag[0]=='u')))
    SerialTrsmInternalLeftUpperConj(Diag::Unit::use_unit_diag,
                                    B.extent(0), B.extent(1),
                                    alpha,
                                    A.data(), A.stride(1), A.stride(0),
                                    B.data(), B.stride(0), B.stride(1));
  if (((side[0]=='L')||(side[0]=='l'))&&((uplo[0]=='L')||(uplo[0]=='l'))&&((trans[0]=='C')||(trans[0]=='c'))&&((diag[0]=='N')||(diag[0]=='n')))
    SerialTrsmInternalLeftUpperConj(Diag::NonUnit::use_unit_diag,
                                    B.extent(0), B.extent(1),
                                    alpha,
                                    A.data(), A.stride(1), A.stride(0),
                                    B.data(), B.stride(0), B.stride(1));

  //Side::Left, Uplo::Upper, Trans::NoTranspose
  if (((side[0]=='L')||(side[0]=='l'))&&((uplo[0]=='U')||(uplo[0]=='u'))&&((trans[0]=='N')||(trans[0]=='n'))&&((diag[0]=='U')||(diag[0]=='u')))
    SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(Diag::Unit::use_unit_diag,
                                                               B.extent(0), B.extent(1),
                                                               alpha,
                                                               A.data(), A.stride(0), A.stride(1),
                                                               B.data(), B.stride(0), B.stride(1));
  if (((side[0]=='L')||(side[0]=='l'))&&((uplo[0]=='U')||(uplo[0]=='u'))&&((trans[0]=='N')||(trans[0]=='n'))&&((diag[0]=='N')||(diag[0]=='n')))
    SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(Diag::NonUnit::use_unit_diag,
                                                               B.extent(0), B.extent(1),
                                                               alpha,
                                                               A.data(), A.stride(0), A.stride(1),
                                                               B.data(), B.stride(0), B.stride(1));

  //Side::Left, Uplo::Upper, Trans::Transpose
  if (((side[0]=='L')||(side[0]=='l'))&&((uplo[0]=='U')||(uplo[0]=='u'))&&((trans[0]=='T')||(trans[0]=='t'))&&((diag[0]=='U')||(diag[0]=='u')))
    SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(Diag::Unit::use_unit_diag,
                                                               B.extent(0), B.extent(1),
                                                               alpha,
                                                               A.data(), A.stride(1), A.stride(0),
                                                               B.data(), B.stride(0), B.stride(1));
  if (((side[0]=='L')||(side[0]=='l'))&&((uplo[0]=='U')||(uplo[0]=='u'))&&((trans[0]=='T')||(trans[0]=='t'))&&((diag[0]=='N')||(diag[0]=='n')))
    SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(Diag::NonUnit::use_unit_diag,
                                                               B.extent(0), B.extent(1),
                                                               alpha,
                                                               A.data(), A.stride(1), A.stride(0),
                                                               B.data(), B.stride(0), B.stride(1));

  //Side::Left, Uplo::Upper, Trans::ConjTranspose
  if (((side[0]=='L')||(side[0]=='l'))&&((uplo[0]=='U')||(uplo[0]=='u'))&&((trans[0]=='C')||(trans[0]=='c'))&&((diag[0]=='U')||(diag[0]=='u')))
    SerialTrsmInternalLeftLowerConj(Diag::Unit::use_unit_diag,
                                    B.extent(0), B.extent(1),
                                    alpha,
                                    A.data(), A.stride(1), A.stride(0),
                                    B.data(), B.stride(0), B.stride(1));
  if (((side[0]=='L')||(side[0]=='l'))&&((uplo[0]=='U')||(uplo[0]=='u'))&&((trans[0]=='C')||(trans[0]=='c'))&&((diag[0]=='N')||(diag[0]=='n')))
    SerialTrsmInternalLeftLowerConj(Diag::NonUnit::use_unit_diag,
                                    B.extent(0), B.extent(1),
                                    alpha,
                                    A.data(), A.stride(1), A.stride(0),
                                    B.data(), B.stride(0), B.stride(1));
  ////
  //Side::Right, Uplo::Lower, Trans::NoTranspose
  if (((side[0]=='R')||(side[0]=='r'))&&((uplo[0]=='L')||(uplo[0]=='l'))&&((trans[0]=='N')||(trans[0]=='n'))&&((diag[0]=='U')||(diag[0]=='u')))
    SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(Diag::Unit::use_unit_diag,
                                                               B.extent(1), B.extent(0),
                                                               alpha,
                                                               A.data(), A.stride(1), A.stride(0),
                                                               B.data(), B.stride(1), B.stride(0));
  if (((side[0]=='R')||(side[0]=='r'))&&((uplo[0]=='L')||(uplo[0]=='l'))&&((trans[0]=='N')||(trans[0]=='n'))&&((diag[0]=='N')||(diag[0]=='n')))
    SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(Diag::NonUnit::use_unit_diag,
                                                               B.extent(1), B.extent(0),
                                                               alpha,
                                                               A.data(), A.stride(1), A.stride(0),
                                                               B.data(), B.stride(1), B.stride(0));

  //Side::Right, Uplo::Lower, Trans::Transpose
  if (((side[0]=='R')||(side[0]=='r'))&&((uplo[0]=='L')||(uplo[0]=='l'))&&((trans[0]=='T')||(trans[0]=='t'))&&((diag[0]=='U')||(diag[0]=='u')))
    SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(Diag::Unit::use_unit_diag,
                                                               B.extent(1), B.extent(0),
                                                               alpha,
                                                               A.data(), A.stride(0), A.stride(1),
                                                               B.data(), B.stride(1), B.stride(0));
  if (((side[0]=='R')||(side[0]=='r'))&&((uplo[0]=='L')||(uplo[0]=='l'))&&((trans[0]=='T')||(trans[0]=='t'))&&((diag[0]=='N')||(diag[0]=='n')))
    SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(Diag::NonUnit::use_unit_diag,
                                                               B.extent(1), B.extent(0),
                                                               alpha,
                                                               A.data(), A.stride(0), A.stride(1),
                                                               B.data(), B.stride(1), B.stride(0));

  //Side::Right, Uplo::Lower, Trans::ConjTranspose
  if (((side[0]=='R')||(side[0]=='r'))&&((uplo[0]=='L')||(uplo[0]=='l'))&&((trans[0]=='C')||(trans[0]=='c'))&&((diag[0]=='U')||(diag[0]=='u')))
    SerialTrsmInternalLeftLowerConj(Diag::Unit::use_unit_diag,
                                    B.extent(1), B.extent(0),
                                    alpha,
                                    A.data(), A.stride(0), A.stride(1),
                                    B.data(), B.stride(1), B.stride(0));
  if (((side[0]=='R')||(side[0]=='r'))&&((uplo[0]=='L')||(uplo[0]=='l'))&&((trans[0]=='C')||(trans[0]=='c'))&&((diag[0]=='N')||(diag[0]=='n')))
    SerialTrsmInternalLeftLowerConj(Diag::NonUnit::use_unit_diag,
                                    B.extent(1), B.extent(0),
                                    alpha,
                                    A.data(), A.stride(0), A.stride(1),
                                    B.data(), B.stride(1), B.stride(0));

  //Side::Right, Uplo::Upper, Trans::NoTranspose
  if (((side[0]=='R')||(side[0]=='r'))&&((uplo[0]=='U')||(uplo[0]=='u'))&&((trans[0]=='N')||(trans[0]=='n'))&&((diag[0]=='U')||(diag[0]=='u')))
    SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(Diag::Unit::use_unit_diag,
                                                               B.extent(1), B.extent(0),
                                                               alpha,
                                                               A.data(), A.stride(1), A.stride(0),
                                                               B.data(), B.stride(1), B.stride(0));
  if (((side[0]=='R')||(side[0]=='r'))&&((uplo[0]=='U')||(uplo[0]=='u'))&&((trans[0]=='N')||(trans[0]=='n'))&&((diag[0]=='N')||(diag[0]=='n')))
    SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(Diag::NonUnit::use_unit_diag,
                                                               B.extent(1), B.extent(0),
                                                               alpha,
                                                               A.data(), A.stride(1), A.stride(0),
                                                               B.data(), B.stride(1), B.stride(0));

  //Side::Right, Uplo::Upper, Trans::Transpose
  if (((side[0]=='R')||(side[0]=='r'))&&((uplo[0]=='U')||(uplo[0]=='u'))&&((trans[0]=='T')||(trans[0]=='t'))&&((diag[0]=='U')||(diag[0]=='u')))
    SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(Diag::Unit::use_unit_diag,
                                                               B.extent(1), B.extent(0),
                                                               alpha,
                                                               A.data(), A.stride(0), A.stride(1),
                                                               B.data(), B.stride(1), B.stride(0));
  if (((side[0]=='R')||(side[0]=='r'))&&((uplo[0]=='U')||(uplo[0]=='u'))&&((trans[0]=='T')||(trans[0]=='t'))&&((diag[0]=='N')||(diag[0]=='n')))
    SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(Diag::NonUnit::use_unit_diag,
                                                               B.extent(1), B.extent(0),
                                                               alpha,
                                                               A.data(), A.stride(0), A.stride(1),
                                                               B.data(), B.stride(1), B.stride(0));
															   
  //Side::Right, Uplo::Upper, Trans::ConjTranspose
  if (((side[0]=='R')||(side[0]=='r'))&&((uplo[0]=='U')||(uplo[0]=='u'))&&((trans[0]=='C')||(trans[0]=='c'))&&((diag[0]=='U')||(diag[0]=='u')))
    SerialTrsmInternalLeftUpperConj(Diag::Unit::use_unit_diag,
                                    B.extent(1), B.extent(0),
                                    alpha,
                                    A.data(), A.stride(0), A.stride(1),
                                    B.data(), B.stride(1), B.stride(0));
  if (((side[0]=='R')||(side[0]=='r'))&&((uplo[0]=='U')||(uplo[0]=='u'))&&((trans[0]=='C')||(trans[0]=='c'))&&((diag[0]=='N')||(diag[0]=='n')))
    SerialTrsmInternalLeftUpperConj(Diag::NonUnit::use_unit_diag,
                                    B.extent(1), B.extent(0),
                                    alpha,
                                    A.data(), A.stride(0), A.stride(1),
                                    B.data(), B.stride(1), B.stride(0));
}

}// namespace Impl
}// namespace KokkosBlas
#endif // KOKKOSBLAS3_TRSM_IMPL_HPP_
