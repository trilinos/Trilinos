// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
// Disabled FAD support due to long build times on cuda (in debug mode
// it takes multiple hours on some platforms). If we need
// sensitivities wrt coordinates, we can reenable.
// POINT_VALUES_INSTANTIATION(panzer::Traits::FadType)
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
POINT_VALUES_INSTANTIATION(panzer::Traits::HessianType)
#endif

// This is very complicated for reasons I don't fully understand...
POINT_VALUES_INSTANTIATION2(panzer::Traits::RealType,panzer::Traits::RealType)
// Disabled FAD support due to long build times on cuda (in debug mode
// it takes multiple hours on some platforms). If we need
// sensitivities wrt coordinates, we can reenable.
// POINT_VALUES_INSTANTIATION2(panzer::Traits::FadType,panzer::Traits::RealType)
// POINT_VALUES_INSTANTIATION2(panzer::Traits::FadType,panzer::Traits::FadType)
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
POINT_VALUES_INSTANTIATION2(panzer::Traits::HessianType,panzer::Traits::RealType)
#endif

} // namespace panzer
