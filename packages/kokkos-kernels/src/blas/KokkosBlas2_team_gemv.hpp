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

#ifndef KOKKOSBLAS2_TEAM_GEMV_HPP_
#define KOKKOSBLAS2_TEAM_GEMV_HPP_

#include<KokkosBlas2_team_gemv_spec.hpp>

namespace KokkosBlas {
namespace Experimental {

template<class TeamType, class MatrixType, class XVector,class YVector>
void KOKKOS_INLINE_FUNCTION gemv (const TeamType& team,
                 const char trans,
                 const typename MatrixType::non_const_value_type& alpha,
                 const MatrixType& A,
                 const XVector& x,
                 const typename YVector::non_const_value_type& beta,
                 const YVector& y)
{
  if (trans == 'N' || trans == 'n')
    return Impl::TeamGEMV<TeamType,MatrixType,XVector,YVector,0>::team_gemv(team,alpha,A,x,beta,y);
  if (trans == 'T' || trans == 't')
    return Impl::TeamGEMV<TeamType,MatrixType,XVector,YVector,1>::team_gemv(team,alpha,A,x,beta,y);
  if (trans == 'C' || trans == 'c')
    return Impl::TeamGEMV<TeamType,MatrixType,XVector,YVector,2>::team_gemv(team,alpha,A,x,beta,y);
}

}
}

#endif
