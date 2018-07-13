// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include "PanzerDiscFE_config.hpp"
#include "Panzer_PointValues2.hpp"
#include "Panzer_PointValues2_impl.hpp"

#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_Traits.hpp"


namespace panzer {

// do some explicit instantiation so things build faster.

#define POINT_VALUES_INSTANTIATION(SCALAR) \
template class PointValues2<SCALAR>;

#define POINT_VALUES_INSTANTIATION2(SCALAR,SCALAR2)\
template void PointValues2<SCALAR>::copyNodeCoords<PHX::MDField<SCALAR2> >(const PHX::MDField<SCALAR2> & in_node_coords); \
template void PointValues2<SCALAR>::copyNodeCoords<PHX::MDField<SCALAR2,Cell,NODE,Dim> >(const PHX::MDField<SCALAR2,Cell,NODE,Dim> & in_node_coords); \
template void PointValues2<SCALAR>::copyNodeCoords<Kokkos::DynRankView<SCALAR2,PHX::Device> >(const Kokkos::DynRankView<SCALAR2,PHX::Device> & in_node_coords); \
\
template void PointValues2<SCALAR>::copyPointCoords<PHX::MDField<SCALAR2> >(const PHX::MDField<SCALAR2> & in_node_coords); \
template void PointValues2<SCALAR>::copyPointCoords<PHX::MDField<SCALAR2,BASIS,Dim> >(const PHX::MDField<SCALAR2,BASIS,Dim> & in_node_coords); \
template void PointValues2<SCALAR>::copyPointCoords<Kokkos::DynRankView<SCALAR2,PHX::Device> >(const Kokkos::DynRankView<SCALAR2,PHX::Device> & in_node_coords);

// special case for PointGenerator....yikes!
template void PointValues2<double>::copyPointCoords<Kokkos::DynRankView<double> >(const Kokkos::DynRankView<double> & in_node_coords);

POINT_VALUES_INSTANTIATION(panzer::Traits::RealType)
POINT_VALUES_INSTANTIATION(panzer::Traits::FadType)
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
POINT_VALUES_INSTANTIATION(panzer::Traits::HessianType)
#endif

// This is very complicated for reasons I don't fully understand...
POINT_VALUES_INSTANTIATION2(panzer::Traits::RealType,panzer::Traits::RealType)
POINT_VALUES_INSTANTIATION2(panzer::Traits::FadType,panzer::Traits::RealType)
POINT_VALUES_INSTANTIATION2(panzer::Traits::FadType,panzer::Traits::FadType)
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
POINT_VALUES_INSTANTIATION2(panzer::Traits::HessianType,panzer::Traits::RealType)
#endif

} // namespace panzer
