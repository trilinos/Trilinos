// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ProjectionToolsDefHVOL.hpp
    \brief  Header file for the Intrepid2::ProjectionTools
            containing definitions for HVOL projections.
    \author Created by Mauro Perego
 */

#ifndef __INTREPID2_PROJECTIONTOOLSDEFHVOL_HPP__
#define __INTREPID2_PROJECTIONTOOLSDEFHVOL_HPP__

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


namespace Intrepid2 {

template<typename DeviceType>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<DeviceType>::getHVolBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtTargetEPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<DeviceType, typename BasisType::scalarType> * projStruct){

  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,DeviceType> ScalarViewType;
  ordinal_type dim = cellBasis->getBaseCellTopology().getDimension();

  ordinal_type basisCardinality = cellBasis->getCardinality();

  ordinal_type numCells = targetAtTargetEPoints.extent(0);

  auto refTargetEWeights = projStruct->getTargetEvalWeights(dim,0);
  auto targetEPointsRange = projStruct->getTargetPointsRange();

  auto refBasisEWeights = projStruct->getBasisEvalWeights(dim,0);
  auto basisEPointsRange = projStruct->getBasisPointsRange();

  ordinal_type numTargetEPoints = range_size(targetEPointsRange(dim,0));
  ordinal_type numBasisEPoints = range_size(basisEPointsRange(dim,0));

  ScalarViewType basisAtBasisEPoints("basisAtBasisEPoints", 1, basisCardinality, numBasisEPoints);
  ScalarViewType basisAtTargetEPoints("basisAtTargetEPoints", basisCardinality, numTargetEPoints);

  auto basisEPoints = projStruct->getAllEvalPoints(EvalPointsType::BASIS);
  auto targetEPoints = projStruct->getAllEvalPoints(EvalPointsType::TARGET);

  cellBasis->getValues(Kokkos::subview(basisAtBasisEPoints, 0, Kokkos::ALL(), Kokkos::ALL()), basisEPoints);
  cellBasis->getValues(basisAtTargetEPoints, targetEPoints);

  ScalarViewType weightedBasisAtTargetEPoints("weightedBasisAtTargetEPoints_",numCells, basisCardinality, numTargetEPoints);
  ScalarViewType weightedBasisAtBasisEPoints("weightedBasisAtBasisEPoints", 1, basisCardinality, numBasisEPoints);

  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(MemSpaceType(), cellBasis->getAllDofOrdinal());
  auto cellDofs = Kokkos::subview(tagToOrdinal, dim, 0, Kokkos::ALL());

  ScalarViewType
  massMat0("massMat0", 1, basisCardinality, basisCardinality),
  massMat("massMat", numCells, basisCardinality, basisCardinality),
  rhsMat("rhsMat", numCells, basisCardinality  );

  ordinal_type offsetBasis = basisEPointsRange(dim,0).first;
  ordinal_type offsetTarget = targetEPointsRange(dim,0).first;

  using HostSpaceType = Kokkos::DefaultHostExecutionSpace;
  auto hWeightedBasisAtBasisEPoints = Kokkos::create_mirror_view(weightedBasisAtBasisEPoints);
  auto hWeightedBasisAtTargetEPoints = Kokkos::create_mirror_view(weightedBasisAtTargetEPoints);
  auto hBasisAtBasisEPoints = Kokkos::create_mirror_view_and_copy(HostSpaceType(), basisAtBasisEPoints);
  auto hBasisAtTargetEPoints = Kokkos::create_mirror_view_and_copy(HostSpaceType(), basisAtTargetEPoints);

  for(ordinal_type j=0; j <basisCardinality; ++j) {
    ordinal_type idof = cellBasis->getDofOrdinal(dim, 0, j);
      for(ordinal_type iq=0; iq <ordinal_type(refBasisEWeights.extent(0)); ++iq)
        hWeightedBasisAtBasisEPoints(0,j,iq) = hBasisAtBasisEPoints(0,idof,offsetBasis+iq) * refBasisEWeights(iq);
      for(ordinal_type iq=0; iq <ordinal_type(refTargetEWeights.extent(0)); ++iq)
        hWeightedBasisAtTargetEPoints(0,j,iq) = hBasisAtTargetEPoints(idof,offsetTarget+iq)* refTargetEWeights(iq);
  }
  Kokkos::deep_copy(weightedBasisAtBasisEPoints,hWeightedBasisAtBasisEPoints);
  Kokkos::deep_copy(weightedBasisAtTargetEPoints,hWeightedBasisAtTargetEPoints);
  FunctionSpaceTools<DeviceType >::integrate(massMat0, basisAtBasisEPoints, weightedBasisAtBasisEPoints);
  RealSpaceTools<DeviceType>::clone(massMat, Kokkos::subview(massMat0,0,Kokkos::ALL(), Kokkos::ALL()));
  RealSpaceTools<DeviceType>::clone(weightedBasisAtTargetEPoints,  Kokkos::subview(weightedBasisAtTargetEPoints,0,Kokkos::ALL(), Kokkos::ALL()));
  FunctionSpaceTools<DeviceType >::integrate(rhsMat, targetAtTargetEPoints, weightedBasisAtTargetEPoints);

  typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, DeviceType> WorkArrayViewType;
  ScalarViewType t_("t",numCells, basisCardinality);
  WorkArrayViewType w_("w",numCells,basisCardinality);

  ElemSystem cellSystem("cellSystem", true);
  cellSystem.solve(basisCoeffs, massMat, rhsMat, t_, w_, cellDofs, basisCardinality);
}

} // Intrepid2 namespace

#endif

