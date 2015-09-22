/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER
*/

#ifndef TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_LOCAL_DEEP_COPY_UQ_PCE_HPP
#define TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_LOCAL_DEEP_COPY_UQ_PCE_HPP

#include "Tpetra_KokkosRefactor_Details_MultiVectorLocalDeepCopy.hpp"
#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"

namespace Tpetra {
namespace Details {

template<class DT, class DL, class DD, class DM,
           class ST, class SL, class SD, class SM,
           class DstWhichVecsType,
           class SrcWhichVecsType>
  void
  localDeepCopy (
    const Kokkos::View<DT,DL,DD,DM,Kokkos::Impl::ViewPCEContiguous>& dst,
    const Kokkos::View<ST,SL,SD,SM,Kokkos::Impl::ViewPCEContiguous>& src,
    const bool dstConstStride, const bool srcConstStride,
    const DstWhichVecsType& dstWhichVecs,
    const SrcWhichVecsType& srcWhichVecs)
  {
    typedef Kokkos::View<DT,DL,DD,DM,Kokkos::Impl::ViewPCEContiguous> DstViewType;
    typedef Kokkos::View<ST,SL,SD,SM,Kokkos::Impl::ViewPCEContiguous> SrcViewType;
    typename DstViewType::flat_array_type dst_flat = dst;
    typename SrcViewType::flat_array_type src_flat = src;
    localDeepCopy( dst_flat, src_flat, dstConstStride, srcConstStride,
                   dstWhichVecs, srcWhichVecs );
  }

  template<class DT, class DL, class DD, class DM,
           class ST, class SL, class SD, class SM>
  void
  localDeepCopyConstStride (
    const Kokkos::View<DT,DL,DD,DM,Kokkos::Impl::ViewPCEContiguous>& dst,
    const Kokkos::View<ST,SL,SD,SM,Kokkos::Impl::ViewPCEContiguous>& src)
  {
    typedef Kokkos::View<DT,DL,DD,DM,Kokkos::Impl::ViewPCEContiguous> DstViewType;
    typedef Kokkos::View<ST,SL,SD,SM,Kokkos::Impl::ViewPCEContiguous> SrcViewType;
    typename DstViewType::flat_array_type dst_flat = dst;
    typename SrcViewType::flat_array_type src_flat = src;
    localDeepCopyConstStride( dst_flat, src_flat );
  }

} // Details namespace
} // Tpetra namespace

#endif // TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_LOCAL_DEEP_COPY_UQ_PCE_HPP
