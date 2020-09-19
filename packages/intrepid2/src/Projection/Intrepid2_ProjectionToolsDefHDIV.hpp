// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_ProjectionToolsDefHDIV.hpp
    \brief  Header file for the Intrepid2::Experimental::ProjectionTools
            containing definitions for HDIV projections.
    \author Created by Mauro Perego
 */

#ifndef __INTREPID2_PROJECTIONTOOLSDEFHDIV_HPP__
#define __INTREPID2_PROJECTIONTOOLSDEFHDIV_HPP__

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


namespace Intrepid2 {
namespace Experimental {


template<typename ViewType1, typename ViewType2, typename ViewType3>
struct ComputeBasisCoeffsOnSides_HDiv {
  const ViewType1 sideBasisNormalAtBasisEPoints_;
  const ViewType1 basisAtBasisEPoints_;
  const ViewType2 basisEWeights_;
  const ViewType1 wBasisDofAtBasisEPoints_;
  const ViewType2 targetEWeights_;
  const ViewType1 basisAtTargetEPoints_;
  const ViewType1 wBasisDofAtTargetEPoints_;
  const ViewType3 tagToOrdinal_;
  const ViewType1 targetAtEPoints_;
  const ViewType1 targetAtTargetEPoints_;
  const ViewType1 refSidesNormal_;
  ordinal_type sideCardinality_;
  ordinal_type offsetBasis_;
  ordinal_type offsetTarget_;
  ordinal_type sideDim_;
  ordinal_type dim_;
  ordinal_type iside_;

  ComputeBasisCoeffsOnSides_HDiv(const ViewType1 sideBasisNormalAtBasisEPoints,
      const ViewType1 basisAtBasisEPoints, const ViewType2 basisEWeights,  const ViewType1 wBasisDofAtBasisEPoints,   const ViewType2 targetEWeights,
      const ViewType1 basisAtTargetEPoints, const ViewType1 wBasisDofAtTargetEvalPoint, const ViewType3 tagToOrdinal,
      const ViewType1 targetAtEPoints, const ViewType1 targetAtTargetEPoints,
      const ViewType1 refSidesNormal, ordinal_type sideCardinality, ordinal_type offsetBasis,
      ordinal_type offsetTarget, ordinal_type sideDim,
      ordinal_type dim, ordinal_type iside) :
        sideBasisNormalAtBasisEPoints_(sideBasisNormalAtBasisEPoints),
        basisAtBasisEPoints_(basisAtBasisEPoints), basisEWeights_(basisEWeights), wBasisDofAtBasisEPoints_(wBasisDofAtBasisEPoints), targetEWeights_(targetEWeights),
        basisAtTargetEPoints_(basisAtTargetEPoints), wBasisDofAtTargetEPoints_(wBasisDofAtTargetEvalPoint),
        tagToOrdinal_(tagToOrdinal), targetAtEPoints_(targetAtEPoints),
        targetAtTargetEPoints_(targetAtTargetEPoints),
        refSidesNormal_(refSidesNormal), sideCardinality_(sideCardinality), offsetBasis_(offsetBasis),
        offsetTarget_(offsetTarget), sideDim_(sideDim), dim_(dim), iside_(iside)
  {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    //Note: we are not considering the jacobian of the orientation map since it is simply a scalar term for the integrals and it does not affect the projection
    for(ordinal_type j=0; j <sideCardinality_; ++j) {
      ordinal_type jdof = tagToOrdinal_(sideDim_, iside_, j);

    for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq) {
      for(ordinal_type d=0; d <dim_; ++d)
        sideBasisNormalAtBasisEPoints_(ic,j,iq) += refSidesNormal_(iside_,d)*basisAtBasisEPoints_(ic,jdof,offsetBasis_+iq,d);
      wBasisDofAtBasisEPoints_(ic,j,iq) = sideBasisNormalAtBasisEPoints_(ic,j,iq) * basisEWeights_(iq);
    }
    for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {
      typename ViewType2::value_type sum=0;
      for(ordinal_type d=0; d <dim_; ++d)
        sum += refSidesNormal_(iside_,d)*basisAtTargetEPoints_(ic,jdof,offsetTarget_+iq,d);
      wBasisDofAtTargetEPoints_(ic,j,iq) = sum * targetEWeights_(iq);
    }
  }

  for(ordinal_type d=0; d <dim_; ++d)
      for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq)
        targetAtTargetEPoints_(ic,iq) += refSidesNormal_(iside_,d)*targetAtEPoints_(ic,offsetTarget_+iq,d);
  }
};


template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5>
struct ComputeBasisCoeffsOnCells_HDiv {
  const ViewType1 basisCoeffs_;
  const ViewType2 negPartialProjAtBasisEPoints_;
  const ViewType2 nonWeightedBasisAtBasisEPoints_;
  const ViewType2 basisAtBasisEPoints_;
  const ViewType3 basisEWeights_;
  const ViewType2 wBasisAtBasisEPoints_;
  const ViewType3 targetEWeights_;
  const ViewType2 basisAtTargetEPoints_;
  const ViewType2 wBasisAtTargetEPoints_;
  const ViewType4 computedDofs_;
  const ViewType5 cellDof_;
  ordinal_type numCellDofs_;
  ordinal_type offsetBasis_;
  ordinal_type offsetTarget_;
  ordinal_type numSideDofs_;

  ComputeBasisCoeffsOnCells_HDiv(const ViewType1 basisCoeffs, ViewType2 negPartialProjAtBasisEPoints,  const ViewType2 nonWeightedBasisAtBasisEPoints,
      const ViewType2 basisAtBasisEPoints, const ViewType3 basisEWeights,  const ViewType2 wBasisAtBasisEPoints,   const ViewType3 targetEWeights,
      const ViewType2 basisAtTargetEPoints, const ViewType2 wBasisAtTargetEPoints, const ViewType4 computedDofs, const ViewType5 cellDof,
      ordinal_type numCellDofs, ordinal_type offsetBasis, ordinal_type offsetTarget, ordinal_type numSideDofs) :
        basisCoeffs_(basisCoeffs), negPartialProjAtBasisEPoints_(negPartialProjAtBasisEPoints), nonWeightedBasisAtBasisEPoints_(nonWeightedBasisAtBasisEPoints),
        basisAtBasisEPoints_(basisAtBasisEPoints), basisEWeights_(basisEWeights), wBasisAtBasisEPoints_(wBasisAtBasisEPoints), targetEWeights_(targetEWeights),
        basisAtTargetEPoints_(basisAtTargetEPoints), wBasisAtTargetEPoints_(wBasisAtTargetEPoints),
        computedDofs_(computedDofs), cellDof_(cellDof),numCellDofs_(numCellDofs), offsetBasis_(offsetBasis),
        offsetTarget_(offsetTarget), numSideDofs_(numSideDofs) {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    for(ordinal_type j=0; j <numCellDofs_; ++j) {
      ordinal_type idof = cellDof_(j);
      for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq) {
        nonWeightedBasisAtBasisEPoints_(ic,j,iq) = basisAtBasisEPoints_(ic,idof,offsetBasis_+iq);
        wBasisAtBasisEPoints_(ic,j,iq) = nonWeightedBasisAtBasisEPoints_(ic,j,iq) * basisEWeights_(iq);
      }
      for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {
        wBasisAtTargetEPoints_(ic,j,iq) = basisAtTargetEPoints_(ic,idof,offsetTarget_+iq)* targetEWeights_(iq);
      }
    }
    for(ordinal_type j=0; j < numSideDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq)
        negPartialProjAtBasisEPoints_(ic,iq) -=  basisCoeffs_(ic,jdof)*basisAtBasisEPoints_(ic,jdof,offsetBasis_+iq);
    }
  }
};


template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5, typename ViewType6>
struct ComputeHCurlBasisCoeffsOnCells_HDiv {
  const ViewType1 basisCoeffs_;
  const ViewType2 negPartialProjAtBasisEPoints_;
  const ViewType2 nonWeightedBasisAtBasisEPoints_;
  const ViewType2 basisAtBasisEPoints_;
  const ViewType2 hcurlBasisCurlAtBasisEPoints_;
  const ViewType3 basisEWeights_;
  const ViewType2 wHCurlBasisAtBasisEPoints_;
  const ViewType3 targetEWeights_;
  const ViewType2 hcurlBasisCurlAtTargetEPoints_;
  const ViewType2 wHCurlBasisAtTargetEPoints_;
  const ViewType4 tagToOrdinal_;
  const ViewType5 computedDofs_;
  const ViewType6 hCurlDof_;
  ordinal_type numCellDofs_;
  ordinal_type offsetBasis_;
  ordinal_type numSideDofs_;
  ordinal_type dim_;

  ComputeHCurlBasisCoeffsOnCells_HDiv(const ViewType1 basisCoeffs, ViewType2 negPartialProjAtBasisEPoints,  const ViewType2 nonWeightedBasisAtBasisEPoints,
      const ViewType2 basisAtBasisEPoints, const ViewType2 hcurlBasisCurlAtBasisEPoints, const ViewType3 basisEWeights,  const ViewType2 wHCurlBasisAtBasisEPoints,   const ViewType3 targetEWeights,
      const ViewType2 hcurlBasisCurlAtTargetEPoints, const ViewType2 wHCurlBasisAtTargetEPoints, const ViewType4 tagToOrdinal, const ViewType5 computedDofs, const ViewType6 hCurlDof,
      ordinal_type numCellDofs, ordinal_type offsetBasis, ordinal_type numSideDofs, ordinal_type dim) :
        basisCoeffs_(basisCoeffs), negPartialProjAtBasisEPoints_(negPartialProjAtBasisEPoints), nonWeightedBasisAtBasisEPoints_(nonWeightedBasisAtBasisEPoints),
        basisAtBasisEPoints_(basisAtBasisEPoints), hcurlBasisCurlAtBasisEPoints_(hcurlBasisCurlAtBasisEPoints), basisEWeights_(basisEWeights), wHCurlBasisAtBasisEPoints_(wHCurlBasisAtBasisEPoints), targetEWeights_(targetEWeights),
        hcurlBasisCurlAtTargetEPoints_(hcurlBasisCurlAtTargetEPoints), wHCurlBasisAtTargetEPoints_(wHCurlBasisAtTargetEPoints),
        tagToOrdinal_(tagToOrdinal), computedDofs_(computedDofs), hCurlDof_(hCurlDof),numCellDofs_(numCellDofs), offsetBasis_(offsetBasis),
        numSideDofs_(numSideDofs), dim_(dim) {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    ordinal_type numBasisEPoints = basisEWeights_.extent(0);

    for(ordinal_type i=0; i<numSideDofs_; ++i) {
      ordinal_type idof = computedDofs_(i);
      for(ordinal_type iq=0; iq<numBasisEPoints; ++iq){
        for(ordinal_type d=0; d <dim_; ++d)
          negPartialProjAtBasisEPoints_(ic,iq,d) -= basisCoeffs_(ic,idof)*basisAtBasisEPoints_(ic,idof,offsetBasis_+iq,d);
      }
    }

    for(ordinal_type i=0; i<numCellDofs_; ++i) {
      ordinal_type idof = tagToOrdinal_(dim_, 0, i);
      for(ordinal_type iq=0; iq<numBasisEPoints; ++iq) {
        for(ordinal_type d=0; d<dim_; ++d)
          nonWeightedBasisAtBasisEPoints_(ic,i,iq,d) = basisAtBasisEPoints_(ic,idof,offsetBasis_+iq,d);
      }
    }

    for(ordinal_type i=0; i<static_cast<ordinal_type>(hCurlDof_.extent(0)); ++i) {
      ordinal_type idof = hCurlDof_(i);
      for(ordinal_type d=0; d<dim_; ++d) {
        for(ordinal_type iq=0; iq<numBasisEPoints; ++iq) {
          wHCurlBasisAtBasisEPoints_(ic,i,iq,d) = hcurlBasisCurlAtBasisEPoints_(idof,iq,d)*basisEWeights_(iq);
        }
        for(ordinal_type iq=0; iq<static_cast<ordinal_type>(targetEWeights_.extent(0)); ++iq) {
          wHCurlBasisAtTargetEPoints_(ic,i,iq,d) = hcurlBasisCurlAtTargetEPoints_(idof,iq,d)*targetEWeights_(iq);
        }
      }
    }
  }
};


template<typename SpT>
template<typename BasisType,
typename ortValueType,       class ...ortProperties>
void
ProjectionTools<SpT>::getHDivEvaluationPoints(typename BasisType::ScalarViewType targetEPoints,
    typename BasisType::ScalarViewType targetDivEPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct,
    const EvalPointsType evalPointType){
  auto refTopologyKey =  Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTopologyKey());
  //const auto cellTopoKey = cellBasis->getBaseCellTopology().getKey();
  ordinal_type dim = cellBasis->getBaseCellTopology().getDimension();
  ordinal_type sideDim = dim-1;

  ordinal_type numSides = cellBasis->getBaseCellTopology().getSideCount();

  ordinal_type numCells = orts.extent(0);

  CellTools<SpT>::setSubcellParametrization();
  typename CellTools<SpT>::subcellParamViewType subcellParamSide;
  if(numSides>0)
    CellTools<SpT>::getSubcellParametrization(subcellParamSide,  sideDim, cellBasis->getBaseCellTopology());

  auto evalPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getPointsRange(evalPointType));

  for(ordinal_type is=0; is<numSides; ++is)  {

    const auto topoKey = refTopologyKey(sideDim,is);
    auto sideBasisEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getEvalPoints(sideDim,is,evalPointType));

    Kokkos::parallel_for
      ("Evaluate Points Sides ",
       Kokkos::RangePolicy<SpT, int> (0, numCells),
       KOKKOS_LAMBDA (const size_t ic) {
      auto sidePointsRange = evalPointsRange(sideDim, is);

      ordinal_type sOrt[6];
      if(dim == 3)
        orts(ic).getFaceOrientation(sOrt, numSides);
      else
        orts(ic).getEdgeOrientation(sOrt, numSides);
      ordinal_type ort = sOrt[is];

      Impl::OrientationTools::mapSubcellCoordsToRefCell(Kokkos::subview(targetEPoints,ic,sidePointsRange,Kokkos::ALL()),
          sideBasisEPoints, subcellParamSide, topoKey, is, ort);
    });
  }

  if(cellBasis->getDofCount(dim,0) <= 0)
    return;

  auto evalDivPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getDerivPointsRange(evalPointType));
  auto cellDivPointsRange = evalDivPointsRange(dim, 0);
  auto cellBasisDivEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getDerivEvalPoints(dim,0,evalPointType));

  RealSpaceTools<SpT>::clone(Kokkos::subview(targetDivEPoints, Kokkos::ALL(), cellDivPointsRange, Kokkos::ALL()), cellBasisDivEPoints);


  if(projStruct->getTargetEvalPoints(dim, 0).data() != NULL)
  {
    auto cellPointsRange = evalPointsRange(dim, 0);
    auto cellBasisEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getEvalPoints(dim,0,evalPointType));
    RealSpaceTools<SpT>::clone(Kokkos::subview(targetEPoints, Kokkos::ALL(), cellPointsRange, Kokkos::ALL()), cellBasisEPoints);
  }
}


template<typename SpT>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<SpT>::getHDivBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtEPoints,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetDivAtDivEPoints,
    const typename BasisType::ScalarViewType targetEPoints,
    const typename BasisType::ScalarViewType targetDivEPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct){
  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> ScalarViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type numTotalEvaluationPoints(targetAtEPoints.extent(1)),
      numTotalDivEvaluationPoints(targetDivAtDivEPoints.extent(1));
  ordinal_type basisCardinality = cellBasis->getCardinality();
  ordinal_type numCells = targetAtEPoints.extent(0);
  const ordinal_type sideDim = dim-1;

  const std::string& name = cellBasis->getName();

  ordinal_type numSides = cellBasis->getBaseCellTopology().getSideCount();
  ScalarViewType refSideNormal("refSideNormal", dim);

  ordinal_type numSideDofs(0);
  for(ordinal_type is=0; is<numSides; ++is)
    numSideDofs += cellBasis->getDofCount(sideDim,is);

  Kokkos::View<ordinal_type*, SpT> computedDofs("computedDofs",numSideDofs);

  const Kokkos::RangePolicy<SpT> policy(0, numCells);

  auto targetEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetPointsRange());
  auto basisEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisPointsRange());

  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(), cellBasis->getAllDofOrdinal());

  ordinal_type numTotalBasisEPoints = projStruct->getNumBasisEvalPoints(), numTotalBasisDivEPoints = projStruct->getNumBasisDerivEvalPoints();
  ScalarViewType basisEPoints_("basisEPoints",numCells,numTotalBasisEPoints, dim);
  ScalarViewType basisDivEPoints("basisDivEPoints",numCells,numTotalBasisDivEPoints, dim);
  getHDivEvaluationPoints(basisEPoints_, basisDivEPoints, orts, cellBasis, projStruct, EvalPointsType::BASIS);

  ScalarViewType basisAtBasisEPoints("basisAtBasisEPoints",numCells,basisCardinality, numTotalBasisEPoints, dim);
  ScalarViewType basisAtTargetEPoints("basisAtTargetEPoints",numCells,basisCardinality, numTotalEvaluationPoints, dim);
  {
    ScalarViewType nonOrientedBasisAtBasisEPoints("nonOrientedBasisAtBasisEPoints",numCells,basisCardinality, numTotalBasisEPoints, dim);
    ScalarViewType nonOrientedBasisAtTargetEPoints("nonOrientedBasisAtTargetEPoints",numCells,basisCardinality, numTotalEvaluationPoints, dim);
    for(ordinal_type ic=0; ic<numCells; ++ic) {
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtTargetEPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(targetEPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtBasisEPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(basisEPoints_, ic, Kokkos::ALL(), Kokkos::ALL()));
    }

    OrientationTools<SpT>::modifyBasisByOrientation(basisAtBasisEPoints, nonOrientedBasisAtBasisEPoints, orts, cellBasis);
    OrientationTools<SpT>::modifyBasisByOrientation(basisAtTargetEPoints, nonOrientedBasisAtTargetEPoints, orts, cellBasis);
  }

  ScalarViewType basisDivAtBasisDivEPoints;
  ScalarViewType basisDivAtTargetDivEPoints;
  if(numTotalDivEvaluationPoints>0) {
    ScalarViewType nonOrientedBasisDivAtTargetDivEPoints, nonOrientedBasisDivAtBasisDivEPoints;
    basisDivAtBasisDivEPoints = ScalarViewType ("basisDivAtBasisDivEPoints",numCells,basisCardinality, numTotalBasisDivEPoints);
    nonOrientedBasisDivAtBasisDivEPoints = ScalarViewType ("nonOrientedBasisDivAtBasisDivEPoints",numCells,basisCardinality, numTotalBasisDivEPoints);
    basisDivAtTargetDivEPoints = ScalarViewType("basisDivAtTargetDivEPoints",numCells,basisCardinality, numTotalDivEvaluationPoints);
    nonOrientedBasisDivAtTargetDivEPoints = ScalarViewType("nonOrientedBasisDivAtTargetDivEPoints",numCells,basisCardinality, numTotalDivEvaluationPoints);

    for(ordinal_type ic=0; ic<numCells; ++ic) {
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisDivAtBasisDivEPoints,ic,Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(basisDivEPoints, ic, Kokkos::ALL(), Kokkos::ALL()),OPERATOR_DIV);
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisDivAtTargetDivEPoints,ic,Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(targetDivEPoints, ic, Kokkos::ALL(), Kokkos::ALL()),OPERATOR_DIV);
    }
    OrientationTools<SpT>::modifyBasisByOrientation(basisDivAtBasisDivEPoints, nonOrientedBasisDivAtBasisDivEPoints, orts, cellBasis);
    OrientationTools<SpT>::modifyBasisByOrientation(basisDivAtTargetDivEPoints, nonOrientedBasisDivAtTargetDivEPoints, orts, cellBasis);
  }

  ScalarViewType refSidesNormal("refSidesNormal", numSides, dim);
  ordinal_type computedDofsCount = 0;
  for(ordinal_type is=0; is<numSides; ++is)  {

    ordinal_type sideCardinality = cellBasis->getDofCount(sideDim,is);

    ordinal_type numTargetEPoints = range_size(targetEPointsRange(sideDim,is));
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(sideDim,is));

    for(ordinal_type i=0; i<sideCardinality; ++i)
      computedDofs(computedDofsCount++) = cellBasis->getDofOrdinal(sideDim, is, i);

    auto sideNormal = Kokkos::subview(refSidesNormal,is,Kokkos::ALL());
    auto sideNormalHost = Kokkos::create_mirror_view(sideNormal);
    CellTools<SpT>::getReferenceSideNormal(sideNormalHost, is, cellTopo);
    Kokkos::deep_copy(sideNormal, sideNormalHost);

    ScalarViewType basisNormalAtBasisEPoints("normalBasisAtBasisEPoints",numCells,sideCardinality, numBasisEPoints);
    ScalarViewType wBasisNormalAtBasisEPoints("wBasisNormalAtBasisEPoints",numCells,sideCardinality, numBasisEPoints);
    ScalarViewType wBasisNormalAtTargetEPoints("wBasisNormalAtTargetEPoints",numCells,sideCardinality, numTargetEPoints);
    ScalarViewType targetNormalAtTargetEPoints("targetNormalAtTargetEPoints",numCells, numTargetEPoints);

    ordinal_type offsetBasis = basisEPointsRange(sideDim,is).first;
    ordinal_type offsetTarget = targetEPointsRange(sideDim,is).first;
    auto targetEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetEvalWeights(sideDim,is));
    auto basisEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisEvalWeights(sideDim,is));


    typedef ComputeBasisCoeffsOnSides_HDiv<ScalarViewType,  decltype(basisEWeights), decltype(tagToOrdinal)> functorTypeSide;
    Kokkos::parallel_for(policy, functorTypeSide(basisNormalAtBasisEPoints, basisAtBasisEPoints,
        basisEWeights,  wBasisNormalAtBasisEPoints, targetEWeights,
        basisAtTargetEPoints, wBasisNormalAtTargetEPoints, tagToOrdinal,
        targetAtEPoints, targetNormalAtTargetEPoints,
        refSidesNormal, sideCardinality, offsetBasis,
        offsetTarget, sideDim,
        dim, is));

    ScalarViewType sideMassMat_("sideMassMat_", numCells, sideCardinality+1, sideCardinality+1),
        sideRhsMat_("rhsMat_", numCells, sideCardinality+1);

    ScalarViewType targetEWeights_("targetEWeights", numCells, 1, targetEWeights.extent(0));
    RealSpaceTools<SpT>::clone(targetEWeights_, targetEWeights);

    range_type range_H(0, sideCardinality);
    range_type range_B(sideCardinality, sideCardinality+1);
    ScalarViewType ones("ones",numCells,1,numBasisEPoints);
    Kokkos::deep_copy(ones,1);

    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(sideMassMat_, Kokkos::ALL(), range_H, range_H), basisNormalAtBasisEPoints, wBasisNormalAtBasisEPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(sideMassMat_, Kokkos::ALL(), range_H, range_B), wBasisNormalAtBasisEPoints, ones);

    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(sideRhsMat_, Kokkos::ALL(), range_H), targetNormalAtTargetEPoints, wBasisNormalAtTargetEPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(sideRhsMat_, Kokkos::ALL(), range_B), targetNormalAtTargetEPoints, targetEWeights_);

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
    ScalarViewType t_("t",numCells, sideCardinality+1);
    WorkArrayViewType w_("w",numCells, sideCardinality+1);

    auto sideDof = Kokkos::subview(tagToOrdinal, sideDim, is, Kokkos::ALL());

    ElemSystem sideSystem("sideSystem", false);
    sideSystem.solve(basisCoeffs, sideMassMat_, sideRhsMat_, t_, w_, sideDof, sideCardinality, 1);
  }


  //Cell
  ordinal_type numCellDofs = cellBasis->getDofCount(dim,0);
  if(numCellDofs==0)
    return;

  Basis<SpT,scalarType,scalarType> *hcurlBasis = NULL;
  if(cellTopo.getKey() == shards::getCellTopologyData<shards::Hexahedron<8> >()->key)
    hcurlBasis = new Basis_HCURL_HEX_In_FEM<SpT,scalarType,scalarType>(cellBasis->getDegree());
  else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key)
    hcurlBasis = new Basis_HCURL_TET_In_FEM<SpT,scalarType,scalarType>(cellBasis->getDegree());
  else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Quadrilateral<4> >()->key)
    hcurlBasis = new Basis_HGRAD_QUAD_Cn_FEM<SpT,scalarType,scalarType>(cellBasis->getDegree());
  else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Triangle<3> >()->key)
    hcurlBasis = new Basis_HGRAD_TRI_Cn_FEM<SpT,scalarType,scalarType>(cellBasis->getDegree());
  else  {
    std::stringstream ss;
    ss << ">>> ERROR (Intrepid2::ProjectionTools::getHDivEvaluationPoints): "
        << "Method not implemented for basis " << name;
    INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
  }
  if(hcurlBasis == NULL) return;


  auto targetDivEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetDerivPointsRange());
  auto basisDivEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisDerivPointsRange());

  ordinal_type numTargetDivEPoints = range_size(targetDivEPointsRange(dim,0));
  ordinal_type numBasisDivEPoints = range_size(basisDivEPointsRange(dim,0));

  auto targetDivEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetDerivEvalWeights(dim,0));
  auto divEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisDerivEvalWeights(dim,0));

  ordinal_type offsetBasisDiv = basisDivEPointsRange(dim,0).first;
  ordinal_type offsetTargetDiv = targetDivEPointsRange(dim,0).first;

  ScalarViewType weightedBasisDivAtBasisEPoints("weightedBasisDivAtBasisEPoints",numCells,numCellDofs, numBasisDivEPoints);
  ScalarViewType weightedBasisDivAtTargetEPoints("weightedBasisDivAtTargetEPoints",numCells, numCellDofs, numTargetDivEPoints);
  ScalarViewType basisDivAtBasisEPoints("basisDivAtBasisEPoints",numCells,numCellDofs, numBasisDivEPoints);
  ScalarViewType targetSideDivAtBasisEPoints("targetSideDivAtBasisEPoints",numCells, numBasisDivEPoints);

  auto cellDofs = Kokkos::subview(tagToOrdinal, dim, 0, Kokkos::ALL());
  typedef ComputeBasisCoeffsOnCells_HDiv<decltype(basisCoeffs), ScalarViewType,  decltype(divEWeights), decltype(computedDofs), decltype(cellDofs)> functorType;
  Kokkos::parallel_for(policy, functorType( basisCoeffs, targetSideDivAtBasisEPoints,  basisDivAtBasisEPoints,
      basisDivAtBasisDivEPoints, divEWeights,  weightedBasisDivAtBasisEPoints, targetDivEWeights, basisDivAtTargetDivEPoints, weightedBasisDivAtTargetEPoints,
      computedDofs, cellDofs, numCellDofs, offsetBasisDiv, offsetTargetDiv, numSideDofs));


  ordinal_type hcurlBasisCardinality = hcurlBasis->getCardinality();
  ordinal_type numCurlInteriorDOFs = hcurlBasis->getDofCount(dim,0);

  range_type range_H(0, numCellDofs);
  range_type range_B(numCellDofs, numCellDofs+numCurlInteriorDOFs);


  ScalarViewType massMat_("massMat_",numCells,numCellDofs+numCurlInteriorDOFs,numCellDofs+numCurlInteriorDOFs),
                 rhsMatTrans("rhsMatTrans",numCells,numCellDofs+numCurlInteriorDOFs);

  FunctionSpaceTools<SpT >::integrate(Kokkos::subview(massMat_, Kokkos::ALL(), range_H,range_H), basisDivAtBasisEPoints, weightedBasisDivAtBasisEPoints);
  FunctionSpaceTools<SpT >::integrate(Kokkos::subview(rhsMatTrans, Kokkos::ALL(), range_H), targetDivAtDivEPoints, weightedBasisDivAtTargetEPoints);
  FunctionSpaceTools<SpT >::integrate(Kokkos::subview(rhsMatTrans, Kokkos::ALL(), range_H), targetSideDivAtBasisEPoints, weightedBasisDivAtBasisEPoints,true);

  if(numCurlInteriorDOFs>0){
    ordinal_type numTargetEPoints = range_size(targetEPointsRange(dim,0));
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(dim,0));

    auto targetEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetEvalWeights(dim,0));
    auto basisEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisEvalWeights(dim,0));

    ordinal_type offsetBasis = basisEPointsRange(dim,0).first;

    auto basisEPoints = Kokkos::subview(basisEPoints_, 0, basisEPointsRange(dim,0), Kokkos::ALL());

    ScalarViewType negPartialProjAtBasisEPoints("targetSideAtBasisEPoints",numCells, numBasisEPoints, dim);
    ScalarViewType nonWeightedBasisAtBasisEPoints("basisAtBasisEPoints",numCells,numCellDofs, numBasisEPoints, dim);
    ScalarViewType hcurlBasisCurlAtBasisEPoints("hcurlBasisCurlAtBasisEPoints",hcurlBasisCardinality, numBasisEPoints,dim);
    ScalarViewType hcurlBasisCurlAtTargetEPoints("hcurlBasisCurlAtTargetEPoints", hcurlBasisCardinality,numTargetEPoints, dim);
    ScalarViewType wHcurlBasisCurlAtBasisEPoints("wHcurlBasisHcurlAtBasisEPoints", numCells, numCurlInteriorDOFs, numBasisEPoints,dim);
    ScalarViewType wHcurlBasisCurlAtTargetEPoints("wHcurlBasisHcurlAtTargetEPoints",numCells, numCurlInteriorDOFs, numTargetEPoints,dim);

    hcurlBasis->getValues(hcurlBasisCurlAtBasisEPoints, basisEPoints, OPERATOR_CURL);
    hcurlBasis->getValues(hcurlBasisCurlAtTargetEPoints, Kokkos::subview(targetEPoints,0,targetEPointsRange(dim,0),Kokkos::ALL()), OPERATOR_CURL);

    auto hCurlTagToOrdinal = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(), hcurlBasis->getAllDofOrdinal());
    auto cellHCurlDof = Kokkos::subview(hCurlTagToOrdinal, dim, 0, range_type(0, numCurlInteriorDOFs));

    typedef ComputeHCurlBasisCoeffsOnCells_HDiv<decltype(basisCoeffs), ScalarViewType,  decltype(divEWeights),
        decltype(tagToOrdinal), decltype(computedDofs), decltype(cellDofs)> functorTypeHCurlCells;
    Kokkos::parallel_for(policy, functorTypeHCurlCells(basisCoeffs, negPartialProjAtBasisEPoints,  nonWeightedBasisAtBasisEPoints,
        basisAtBasisEPoints, hcurlBasisCurlAtBasisEPoints, basisEWeights, wHcurlBasisCurlAtBasisEPoints, targetEWeights,
        hcurlBasisCurlAtTargetEPoints, wHcurlBasisCurlAtTargetEPoints, tagToOrdinal, computedDofs, cellHCurlDof,
        numCellDofs, offsetBasis, numSideDofs, dim));

    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(massMat_, Kokkos::ALL(), range_H,range_B), nonWeightedBasisAtBasisEPoints, wHcurlBasisCurlAtBasisEPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(rhsMatTrans, Kokkos::ALL(), range_B), Kokkos::subview(targetAtEPoints, Kokkos::ALL(), targetEPointsRange(dim,0), Kokkos::ALL()), wHcurlBasisCurlAtTargetEPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(rhsMatTrans, Kokkos::ALL(), range_B), negPartialProjAtBasisEPoints, wHcurlBasisCurlAtBasisEPoints,true);
  }
  delete hcurlBasis;

  typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
  ScalarViewType t_("t",numCells, numCellDofs+numCurlInteriorDOFs);
  WorkArrayViewType w_("w",numCells, numCellDofs+numCurlInteriorDOFs);

  ElemSystem cellSystem("cellSystem", true);
  cellSystem.solve(basisCoeffs, massMat_, rhsMatTrans, t_, w_, cellDofs, numCellDofs, numCurlInteriorDOFs);
}



}
}

#endif

