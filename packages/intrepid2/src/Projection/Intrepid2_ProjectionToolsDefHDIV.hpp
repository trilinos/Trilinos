// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ProjectionToolsDefHDIV.hpp
    \brief  Header file for the Intrepid2::ProjectionTools
            containing definitions for HDIV projections.
    \author Created by Mauro Perego
 */

#ifndef __INTREPID2_PROJECTIONTOOLSDEFHDIV_HPP__
#define __INTREPID2_PROJECTIONTOOLSDEFHDIV_HPP__

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


namespace Intrepid2 {

namespace FunctorsProjectionTools {

template<typename ViewType1, typename ViewType2, typename ViewType3>
struct ComputeWBasisSide_HDiv {
  ViewType1 basisNormalAtBasisEPoints_;
  ViewType1 wBasisNormalAtBasisEPoints_;
  const ViewType1 basisAtBasisEPoints_;
  const ViewType2 basisEWeights_;
  const ViewType1 refSideNormal_;
  const ViewType3 tagToOrdinal_;
  ordinal_type sideDim_;
  ordinal_type iside_;
  ordinal_type offsetBasis_;

  ComputeWBasisSide_HDiv(ViewType1 basisNormalAtBasisEPoints, ViewType1 wBasisNormalAtBasisEPoints,
      ViewType1 basisAtBasisEPoints, ViewType2 basisEWeights,  ViewType1 refSideNormal, ViewType3 tagToOrdinal,  
      ordinal_type sideDim, ordinal_type iside, ordinal_type offsetBasis) :
        basisNormalAtBasisEPoints_(basisNormalAtBasisEPoints), wBasisNormalAtBasisEPoints_(wBasisNormalAtBasisEPoints), basisAtBasisEPoints_(basisAtBasisEPoints),
        basisEWeights_(basisEWeights), refSideNormal_(refSideNormal), tagToOrdinal_(tagToOrdinal),
        sideDim_(sideDim), iside_(iside), offsetBasis_(offsetBasis) {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type j, const ordinal_type iq) const {
    ordinal_type jdof = tagToOrdinal_(sideDim_, iside_, j);
    for(ordinal_type d=0; d <ordinal_type(refSideNormal_.extent(0)); ++d)
      basisNormalAtBasisEPoints_(0,j,iq) += refSideNormal_(d)*basisAtBasisEPoints_(jdof,offsetBasis_+iq,d);
    wBasisNormalAtBasisEPoints_(0,j,iq) = basisNormalAtBasisEPoints_(0,j,iq)*basisEWeights_(iq);
  }
};

template<typename ViewType1, typename ViewType2, typename ViewType3, typename ViewType4>
struct ComputeBasisCoeffsOnSides_HDiv {
  const ViewType2 targetEWeights_;
  const ViewType1 basisAtTargetEPoints_;
  const ViewType1 wBasisDofAtTargetEPoints_;
  const ViewType3 tagToOrdinal_;
  const ViewType4 targetAtEPoints_;
  const ViewType1 targetAtTargetEPoints_;
  const ViewType1 refSideNormal_;
  ordinal_type sideCardinality_;
  ordinal_type offsetTarget_;
  ordinal_type sideDim_;
  ordinal_type dim_;
  ordinal_type iside_;

  ComputeBasisCoeffsOnSides_HDiv(const ViewType2 targetEWeights,
      const ViewType1 basisAtTargetEPoints, const ViewType1 wBasisDofAtTargetEvalPoint, const ViewType3 tagToOrdinal,
      const ViewType4 targetAtEPoints, const ViewType1 targetAtTargetEPoints,
      const ViewType1 refSideNormal, ordinal_type sideCardinality,
      ordinal_type offsetTarget, ordinal_type sideDim,
      ordinal_type dim, ordinal_type iside) :
        targetEWeights_(targetEWeights),
        basisAtTargetEPoints_(basisAtTargetEPoints), wBasisDofAtTargetEPoints_(wBasisDofAtTargetEvalPoint),
        tagToOrdinal_(tagToOrdinal), targetAtEPoints_(targetAtEPoints),
        targetAtTargetEPoints_(targetAtTargetEPoints),
        refSideNormal_(refSideNormal), sideCardinality_(sideCardinality),
        offsetTarget_(offsetTarget), sideDim_(sideDim), dim_(dim), iside_(iside)
  {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    typename ViewType1::value_type tmp =0;
    for(ordinal_type j=0; j <sideCardinality_; ++j) {
      ordinal_type jdof = tagToOrdinal_(sideDim_, iside_, j);
      for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {
        tmp=0;
        for(ordinal_type d=0; d <dim_; ++d)
          tmp += refSideNormal_(d)*basisAtTargetEPoints_(jdof,offsetTarget_+iq,d);
        wBasisDofAtTargetEPoints_(ic,j,iq) = tmp * targetEWeights_(iq);
      }
    }

    for(ordinal_type d=0; d <dim_; ++d)
      for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq)
        targetAtTargetEPoints_(ic,iq) += refSideNormal_(d)*targetAtEPoints_(ic,offsetTarget_+iq,d);
  }
};


template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4>
struct ComputeBasisCoeffsOnCells_HDiv {
  const ViewType1 basisCoeffs_;
  const ViewType1 negPartialProjAtBasisEPoints_;
  const ViewType1 nonWeightedBasisAtBasisEPoints_;
  const ViewType1 basisAtBasisEPoints_;
  const ViewType2 basisEWeights_;
  const ViewType1 wBasisAtBasisEPoints_;
  const ViewType2 targetEWeights_;
  const ViewType1 basisAtTargetEPoints_;
  const ViewType1 wBasisAtTargetEPoints_;
  const ViewType3 computedDofs_;
  const ViewType4 cellDof_;
  ordinal_type numCellDofs_;
  ordinal_type offsetBasis_;
  ordinal_type offsetTarget_;
  ordinal_type numSideDofs_;

  ComputeBasisCoeffsOnCells_HDiv(const ViewType1 basisCoeffs, ViewType1 negPartialProjAtBasisEPoints,  const ViewType1 nonWeightedBasisAtBasisEPoints,
      const ViewType1 basisAtBasisEPoints, const ViewType2 basisEWeights,  const ViewType1 wBasisAtBasisEPoints,   const ViewType2 targetEWeights,
      const ViewType1 basisAtTargetEPoints, const ViewType1 wBasisAtTargetEPoints, const ViewType3 computedDofs, const ViewType4 cellDof,
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
        nonWeightedBasisAtBasisEPoints_(0,j,iq) = basisAtBasisEPoints_(idof,offsetBasis_+iq);
        wBasisAtBasisEPoints_(ic,j,iq) = nonWeightedBasisAtBasisEPoints_(0,j,iq) * basisEWeights_(iq);
      }
      for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {
        wBasisAtTargetEPoints_(ic,j,iq) = basisAtTargetEPoints_(idof,offsetTarget_+iq)* targetEWeights_(iq);
      }
    }
    for(ordinal_type j=0; j < numSideDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq)
        negPartialProjAtBasisEPoints_(ic,iq) -=  basisCoeffs_(ic,jdof)*basisAtBasisEPoints_(jdof,offsetBasis_+iq);
    }
  }
};


template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5>
struct ComputeHCurlBasisCoeffsOnCells_HDiv {
  const ViewType1 basisCoeffs_;
  const ViewType1 negPartialProjAtBasisEPoints_;
  const ViewType1 nonWeightedBasisAtBasisEPoints_;
  const ViewType1 basisAtBasisEPoints_;
  const ViewType1 hcurlBasisCurlAtBasisEPoints_;
  const ViewType2 basisEWeights_;
  const ViewType1 wHCurlBasisAtBasisEPoints_;
  const ViewType2 targetEWeights_;
  const ViewType1 hcurlBasisCurlAtTargetEPoints_;
  const ViewType1 wHCurlBasisAtTargetEPoints_;
  const ViewType3 tagToOrdinal_;
  const ViewType4 computedDofs_;
  const ViewType5 hCurlDof_;
  ordinal_type numCellDofs_;
  ordinal_type offsetBasis_;
  ordinal_type numSideDofs_;
  ordinal_type dim_;

  ComputeHCurlBasisCoeffsOnCells_HDiv(const ViewType1 basisCoeffs, ViewType1 negPartialProjAtBasisEPoints, const ViewType1 nonWeightedBasisAtBasisEPoints,
      const ViewType1 basisAtBasisEPoints, const ViewType1 hcurlBasisCurlAtBasisEPoints, const ViewType2 basisEWeights,  const ViewType1 wHCurlBasisAtBasisEPoints,   const ViewType2 targetEWeights,
      const ViewType1 hcurlBasisCurlAtTargetEPoints, const ViewType1 wHCurlBasisAtTargetEPoints, const ViewType3 tagToOrdinal, const ViewType4 computedDofs, const ViewType5 hCurlDof,
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
          negPartialProjAtBasisEPoints_(ic,iq,d) -= basisCoeffs_(ic,idof)*basisAtBasisEPoints_(idof,offsetBasis_+iq,d);
      }
    }

    for(ordinal_type i=0; i<numCellDofs_; ++i) {
      ordinal_type idof = tagToOrdinal_(dim_, 0, i);
      for(ordinal_type iq=0; iq<numBasisEPoints; ++iq) {
        for(ordinal_type d=0; d<dim_; ++d)
          nonWeightedBasisAtBasisEPoints_(0,i,iq,d) = basisAtBasisEPoints_(idof,offsetBasis_+iq,d);
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
}  // FunctorsProjectionTools namespace


template<typename DeviceType>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<DeviceType>::getHDivBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtEPoints,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetDivAtDivEPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<DeviceType, typename BasisType::scalarType> * projStruct){
  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,DeviceType> ScalarViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type basisCardinality = cellBasis->getCardinality();
  ordinal_type numCells = targetAtEPoints.extent(0);
  const ordinal_type sideDim = dim-1;

  const std::string& name = cellBasis->getName();

  ordinal_type numSides = cellBasis->getBaseCellTopology().getSideCount();

  ordinal_type numSideDofs(0);
  for(ordinal_type is=0; is<numSides; ++is)
    numSideDofs += cellBasis->getDofCount(sideDim,is);

  Kokkos::View<ordinal_type*, DeviceType> computedDofs("computedDofs",numSideDofs);

  const Kokkos::RangePolicy<ExecSpaceType> policy(0, numCells);

  auto targetEPointsRange = projStruct->getTargetPointsRange();
  auto basisEPointsRange = projStruct->getBasisPointsRange();

  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(MemSpaceType(), cellBasis->getAllDofOrdinal());

  ordinal_type numTotalBasisEPoints = projStruct->getNumBasisEvalPoints(), numTotalBasisDivEPoints = projStruct->getNumBasisDerivEvalPoints();
  auto basisEPoints = projStruct->getAllEvalPoints(EvalPointsType::BASIS);
  auto basisDivEPoints = projStruct->getAllDerivEvalPoints(EvalPointsType::BASIS);

  ordinal_type numTotalTargetEPoints = projStruct->getNumTargetEvalPoints(), numTotalTargetDivEPoints = projStruct->getNumTargetDerivEvalPoints();
  auto targetEPoints = projStruct->getAllEvalPoints(EvalPointsType::TARGET);
  auto targetDivEPoints = projStruct->getAllDerivEvalPoints(EvalPointsType::TARGET);

  ScalarViewType basisAtBasisEPoints("basisAtBasisEPoints",basisCardinality, numTotalBasisEPoints, dim);
  ScalarViewType basisAtTargetEPoints("basisAtTargetEPoints",basisCardinality, numTotalTargetEPoints, dim);
  cellBasis->getValues(basisAtTargetEPoints, targetEPoints);
  cellBasis->getValues(basisAtBasisEPoints, basisEPoints);

  ScalarViewType basisDivAtBasisDivEPoints;
  ScalarViewType basisDivAtTargetDivEPoints;
  if(numTotalTargetDivEPoints>0) {
    basisDivAtBasisDivEPoints = ScalarViewType ("basisDivAtBasisDivEPoints",basisCardinality, numTotalBasisDivEPoints);
    basisDivAtTargetDivEPoints = ScalarViewType("basisDivAtTargetDivEPoints",basisCardinality, numTotalTargetDivEPoints);
    cellBasis->getValues(basisDivAtBasisDivEPoints, basisDivEPoints, OPERATOR_DIV);
    cellBasis->getValues(basisDivAtTargetDivEPoints, targetDivEPoints, OPERATOR_DIV);
  }

  ScalarViewType refBasisCoeffs("refBasisCoeffs", basisCoeffs.extent(0), basisCoeffs.extent(1));
  
  ordinal_type computedDofsCount = 0;
  for(ordinal_type is=0; is<numSides; ++is)  {

    ordinal_type sideCardinality = cellBasis->getDofCount(sideDim,is);

    ordinal_type numTargetEPoints = range_size(targetEPointsRange(sideDim,is));
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(sideDim,is));

    auto sideDofs = Kokkos::subview(tagToOrdinal, sideDim, is, range_type(0,sideCardinality));
    auto computedSideDofs = Kokkos::subview(computedDofs, range_type(computedDofsCount,computedDofsCount+sideCardinality));
    deep_copy(computedSideDofs, sideDofs);
    computedDofsCount += sideCardinality;

    ScalarViewType refSideNormal("refSideNormal", dim);
    CellTools<DeviceType>::getReferenceSideNormal(refSideNormal, is, cellTopo);

    ScalarViewType basisNormalAtBasisEPoints("normalBasisAtBasisEPoints",1,sideCardinality, numBasisEPoints);
    ScalarViewType wBasisNormalAtBasisEPoints("wBasisNormalAtBasisEPoints",1,sideCardinality, numBasisEPoints);
    ScalarViewType wBasisNormalAtTargetEPoints("wBasisNormalAtTargetEPoints",numCells,sideCardinality, numTargetEPoints);
    ScalarViewType targetNormalAtTargetEPoints("targetNormalAtTargetEPoints",numCells, numTargetEPoints);

    ordinal_type offsetBasis = basisEPointsRange(sideDim,is).first;
    ordinal_type offsetTarget = targetEPointsRange(sideDim,is).first;
    auto targetEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetEvalWeights(sideDim,is));
    auto basisEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisEvalWeights(sideDim,is));

    using functorTypeWBasisEdge = FunctorsProjectionTools::ComputeWBasisSide_HDiv<ScalarViewType,  decltype(basisEWeights), decltype(tagToOrdinal)>;
    Kokkos::parallel_for(Kokkos::MDRangePolicy<ExecSpaceType, Kokkos::Rank<2> >({0,0}, {sideCardinality,numBasisEPoints}), 
      functorTypeWBasisEdge(basisNormalAtBasisEPoints,wBasisNormalAtBasisEPoints,basisAtBasisEPoints,basisEWeights,refSideNormal,tagToOrdinal,sideDim,is,offsetBasis));

    using functorTypeSide = FunctorsProjectionTools::ComputeBasisCoeffsOnSides_HDiv<ScalarViewType,  decltype(targetEWeights), decltype(tagToOrdinal), decltype(targetAtEPoints)>;
    Kokkos::parallel_for(policy, functorTypeSide(targetEWeights,
        basisAtTargetEPoints, wBasisNormalAtTargetEPoints, tagToOrdinal,
        targetAtEPoints, targetNormalAtTargetEPoints,
        refSideNormal, sideCardinality,
        offsetTarget, sideDim,
        dim, is));

    ScalarViewType sideMassMat_("sideMassMat_", 1, sideCardinality+1, sideCardinality+1),
        sideRhsMat_("rhsMat_", numCells, sideCardinality+1);

    ScalarViewType targetEWeights_("targetEWeights", numCells, 1, targetEWeights.extent(0));
    RealSpaceTools<DeviceType>::clone(targetEWeights_, targetEWeights);

    range_type range_H(0, sideCardinality);
    range_type range_B(sideCardinality, sideCardinality+1);
    ScalarViewType ones("ones",1,1,numBasisEPoints);
    Kokkos::deep_copy(ones,1);

    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(sideMassMat_, Kokkos::ALL(), range_H, range_H), basisNormalAtBasisEPoints, wBasisNormalAtBasisEPoints);
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(sideMassMat_, Kokkos::ALL(), range_H, range_B), wBasisNormalAtBasisEPoints, ones);

    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(sideRhsMat_, Kokkos::ALL(), range_H), targetNormalAtTargetEPoints, wBasisNormalAtTargetEPoints);
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(sideRhsMat_, Kokkos::ALL(), range_B), targetNormalAtTargetEPoints, targetEWeights_);

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, DeviceType> WorkArrayViewType;
    ScalarViewType t_("t",numCells, sideCardinality+1);
    WorkArrayViewType w_("w",numCells, sideCardinality+1);

    auto sideDof = Kokkos::subview(tagToOrdinal, sideDim, is, Kokkos::ALL());

    ElemSystem sideSystem("sideSystem", true);
    sideSystem.solve(refBasisCoeffs, sideMassMat_, sideRhsMat_, t_, w_, sideDof, sideCardinality, 1);
  }


  //Cell
  ordinal_type numCellDofs = cellBasis->getDofCount(dim,0);
  if(numCellDofs!=0) {
    Basis<DeviceType,scalarType,scalarType> *hcurlBasis = NULL;
    if(cellTopo.getKey() == shards::getCellTopologyData<shards::Hexahedron<8> >()->key)
      hcurlBasis = new Basis_HCURL_HEX_In_FEM<DeviceType,scalarType,scalarType>(cellBasis->getDegree());
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key)
      hcurlBasis = new Basis_HCURL_TET_In_FEM<DeviceType,scalarType,scalarType>(cellBasis->getDegree());
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Wedge<6> >()->key)
      hcurlBasis = new typename DerivedNodalBasisFamily<DeviceType,scalarType,scalarType>::HCURL_WEDGE(cellBasis->getDegree());
   // TODO: uncomment the next two lines once H(curl) pyramid implemented
   //  else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Pyramid<5> >()->key)
   //    hcurlBasis = new typename DerivedNodalBasisFamily<DeviceType,scalarType,scalarType>::HCURL_PYR(cellBasis->getDegree());
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Quadrilateral<4> >()->key)
      hcurlBasis = new Basis_HGRAD_QUAD_Cn_FEM<DeviceType,scalarType,scalarType>(cellBasis->getDegree());
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Triangle<3> >()->key)
      hcurlBasis = new Basis_HGRAD_TRI_Cn_FEM<DeviceType,scalarType,scalarType>(cellBasis->getDegree());
    else  {
      std::stringstream ss;
      ss << ">>> ERROR (Intrepid2::ProjectionTools::getHDivBasisCoeffs): "
          << "Method not implemented for basis " << name;
      INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
    }
    if(hcurlBasis == NULL) return;


    auto targetDivEPointsRange = projStruct->getTargetDerivPointsRange();
    auto basisDivEPointsRange = projStruct->getBasisDerivPointsRange();

    ordinal_type numTargetDivEPoints = range_size(targetDivEPointsRange(dim,0));
    ordinal_type numBasisDivEPoints = range_size(basisDivEPointsRange(dim,0));

    auto targetDivEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetDerivEvalWeights(dim,0));
    auto divEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisDerivEvalWeights(dim,0));

    ordinal_type offsetBasisDiv = basisDivEPointsRange(dim,0).first;
    ordinal_type offsetTargetDiv = targetDivEPointsRange(dim,0).first;

    ScalarViewType weightedBasisDivAtBasisEPoints("weightedBasisDivAtBasisEPoints",numCells,numCellDofs, numBasisDivEPoints);
    ScalarViewType weightedBasisDivAtTargetEPoints("weightedBasisDivAtTargetEPoints",numCells, numCellDofs, numTargetDivEPoints);
    ScalarViewType basisDivAtBasisEPoints("basisDivAtBasisEPoints",1,numCellDofs, numBasisDivEPoints);
    ScalarViewType targetSideDivAtBasisEPoints("targetSideDivAtBasisEPoints",numCells, numBasisDivEPoints);

    auto cellDofs = Kokkos::subview(tagToOrdinal, dim, 0, Kokkos::ALL());
    using functorType = FunctorsProjectionTools::ComputeBasisCoeffsOnCells_HDiv<ScalarViewType,  decltype(divEWeights), decltype(computedDofs), decltype(cellDofs)>;
    Kokkos::parallel_for(policy, functorType( refBasisCoeffs, targetSideDivAtBasisEPoints,  basisDivAtBasisEPoints,
        basisDivAtBasisDivEPoints, divEWeights,  weightedBasisDivAtBasisEPoints, targetDivEWeights, basisDivAtTargetDivEPoints, weightedBasisDivAtTargetEPoints,
        computedDofs, cellDofs, numCellDofs, offsetBasisDiv, offsetTargetDiv, numSideDofs));


    ordinal_type hcurlBasisCardinality = hcurlBasis->getCardinality();
    ordinal_type numCurlInteriorDOFs = hcurlBasis->getDofCount(dim,0);

    range_type range_H(0, numCellDofs);
    range_type range_B(numCellDofs, numCellDofs+numCurlInteriorDOFs);


    ScalarViewType massMat_("massMat_",1,numCellDofs+numCurlInteriorDOFs,numCellDofs+numCurlInteriorDOFs),
                  rhsMatTrans("rhsMatTrans",numCells,numCellDofs+numCurlInteriorDOFs);

    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(massMat_, Kokkos::ALL(), range_H,range_H), basisDivAtBasisEPoints, Kokkos::subview(weightedBasisDivAtBasisEPoints, std::make_pair(0,1), Kokkos::ALL(), Kokkos::ALL()) );
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(rhsMatTrans, Kokkos::ALL(), range_H), targetDivAtDivEPoints, weightedBasisDivAtTargetEPoints);
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(rhsMatTrans, Kokkos::ALL(), range_H), targetSideDivAtBasisEPoints, weightedBasisDivAtBasisEPoints,true);

    if(numCurlInteriorDOFs>0){
      ordinal_type numTargetEPoints = range_size(targetEPointsRange(dim,0));
      ordinal_type numBasisEPoints = range_size(basisEPointsRange(dim,0));

      auto targetEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetEvalWeights(dim,0));
      auto basisEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisEvalWeights(dim,0));

      ordinal_type offsetBasis = basisEPointsRange(dim,0).first;

      ScalarViewType negPartialProjAtBasisEPoints("targetSideAtBasisEPoints",numCells, numBasisEPoints, dim);
      ScalarViewType nonWeightedBasisAtBasisEPoints("basisAtBasisEPoints",1,numCellDofs, numBasisEPoints, dim);
      ScalarViewType hcurlBasisCurlAtBasisEPoints("hcurlBasisCurlAtBasisEPoints",hcurlBasisCardinality, numBasisEPoints,dim);
      ScalarViewType hcurlBasisCurlAtTargetEPoints("hcurlBasisCurlAtTargetEPoints", hcurlBasisCardinality,numTargetEPoints, dim);
      ScalarViewType wHcurlBasisCurlAtBasisEPoints("wHcurlBasisHcurlAtBasisEPoints", numCells, numCurlInteriorDOFs, numBasisEPoints,dim);
      ScalarViewType wHcurlBasisCurlAtTargetEPoints("wHcurlBasisHcurlAtTargetEPoints",numCells, numCurlInteriorDOFs, numTargetEPoints,dim);

      hcurlBasis->getValues(hcurlBasisCurlAtBasisEPoints, Kokkos::subview(basisEPoints, basisEPointsRange(dim,0), Kokkos::ALL()), OPERATOR_CURL);
      hcurlBasis->getValues(hcurlBasisCurlAtTargetEPoints, Kokkos::subview(targetEPoints,targetEPointsRange(dim,0),Kokkos::ALL()), OPERATOR_CURL);

      auto hCurlTagToOrdinal = Kokkos::create_mirror_view_and_copy(MemSpaceType(), hcurlBasis->getAllDofOrdinal());
      auto cellHCurlDof = Kokkos::subview(hCurlTagToOrdinal, dim, 0, range_type(0, numCurlInteriorDOFs));

      using functorTypeHCurlCells = FunctorsProjectionTools::ComputeHCurlBasisCoeffsOnCells_HDiv<ScalarViewType,  decltype(divEWeights),
          decltype(tagToOrdinal), decltype(computedDofs), decltype(cellDofs)>;
      Kokkos::parallel_for(policy, functorTypeHCurlCells(refBasisCoeffs, negPartialProjAtBasisEPoints,  nonWeightedBasisAtBasisEPoints,
          basisAtBasisEPoints, hcurlBasisCurlAtBasisEPoints, basisEWeights, wHcurlBasisCurlAtBasisEPoints, targetEWeights,
          hcurlBasisCurlAtTargetEPoints, wHcurlBasisCurlAtTargetEPoints, tagToOrdinal, computedDofs, cellHCurlDof,
          numCellDofs, offsetBasis, numSideDofs, dim));

      FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(massMat_, Kokkos::ALL(), range_H,range_B), nonWeightedBasisAtBasisEPoints, Kokkos::subview(wHcurlBasisCurlAtBasisEPoints, std::make_pair(0,1), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) );
      FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(rhsMatTrans, Kokkos::ALL(), range_B), Kokkos::subview(targetAtEPoints, Kokkos::ALL(), targetEPointsRange(dim,0), Kokkos::ALL()), wHcurlBasisCurlAtTargetEPoints);
      FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(rhsMatTrans, Kokkos::ALL(), range_B), negPartialProjAtBasisEPoints, wHcurlBasisCurlAtBasisEPoints,true);
    }
    delete hcurlBasis;

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, DeviceType> WorkArrayViewType;
    ScalarViewType t_("t",numCells, numCellDofs+numCurlInteriorDOFs);
    WorkArrayViewType w_("w",numCells, numCellDofs+numCurlInteriorDOFs);

    ElemSystem cellSystem("cellSystem", true);
    cellSystem.solve(refBasisCoeffs, massMat_, rhsMatTrans, t_, w_, cellDofs, numCellDofs, numCurlInteriorDOFs);
  }

  OrientationTools<DeviceType>::modifyBasisByOrientationInverse(basisCoeffs, refBasisCoeffs, orts, cellBasis, true);

}

}  // Intrepid2 namespace

#endif

