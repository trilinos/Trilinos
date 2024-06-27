// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ProjectionToolsDefHGRAD.hpp
    \brief  Header file for the Intrepid2::ProjectionTools
            containing definitions for HGRAD projections.
    \author Created by Mauro Perego
 */

#ifndef __INTREPID2_PROJECTIONTOOLSDEFHGRAD_HPP__
#define __INTREPID2_PROJECTIONTOOLSDEFHGRAD_HPP__

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Intrepid2_TestUtils.hpp"

namespace Intrepid2 {

namespace FunctorsProjectionTools {

template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5>
struct ComputeBasisCoeffsOnVertices_HGRAD {
  ViewType1 basisCoeffs_;
  const ViewType2 tagToOrdinal_;
  const ViewType3 targetEPointsRange_;
  const ViewType4 targetAtTargetEPoints_;
  const ViewType5 basisAtTargetEPoints_;
  ordinal_type numVertices_;


  ComputeBasisCoeffsOnVertices_HGRAD(ViewType1 basisCoeffs, ViewType2 tagToOrdinal, ViewType3 targetEPointsRange,
      ViewType4  targetAtTargetEPoints, ViewType5 basisAtTargetEPoints, ordinal_type numVertices) :
        basisCoeffs_(basisCoeffs), tagToOrdinal_(tagToOrdinal), targetEPointsRange_(targetEPointsRange),
        targetAtTargetEPoints_(targetAtTargetEPoints), basisAtTargetEPoints_(basisAtTargetEPoints), numVertices_(numVertices) {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {


    for(ordinal_type iv=0; iv<numVertices_; ++iv) {
      ordinal_type idof = tagToOrdinal_(0, iv, 0);
      ordinal_type pt = targetEPointsRange_(0,iv).first;
      //the value of the basis at the vertex might be different than 1; HGrad basis, so the function is scalar
      basisCoeffs_(ic,idof) = targetAtTargetEPoints_(ic,pt)/basisAtTargetEPoints_(idof,pt,0);
    }
  }
};

template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5>
struct ComputeBasisCoeffsOnEdges_HGRAD {
  const ViewType1 basisCoeffs_;
  const ViewType1 negPartialProjGrad_;
  const ViewType1 basisTanAtEPoints_;
  const ViewType1 basisGradAtBasisGradEPoints_;
  const ViewType2 basisGradEWeights_;
  const ViewType1 wBasisAtBasisGradEPoints_;
  const ViewType2 targetGradEWeights_;
  const ViewType1 basisGradAtTargetGradEPoints_;
  const ViewType1 wBasisAtTargetGradEPoints_;
  const ViewType3 computedDofs_;
  const ViewType4 tagToOrdinal_;
  const ViewType1 targetGradTanAtTargetGradEPoints_;
  const ViewType5 targetGradAtTargetGradEPoints_;
  const ViewType1 refEdgeTan_;
  ordinal_type edgeCardinality_;
  ordinal_type offsetBasis_;
  ordinal_type offsetTarget_;
  ordinal_type numVertexDofs_;
  ordinal_type edgeDim_;
  ordinal_type dim_;
  ordinal_type iedge_;

  ComputeBasisCoeffsOnEdges_HGRAD(const ViewType1 basisCoeffs, ViewType1 negPartialProjGrad,  const ViewType1 basisTanAtEPoints,
      const ViewType1 basisGradAtBasisGradEPoints, const ViewType2 basisGradEWeights,  const ViewType1 wBasisAtBasisGradEPoints,   const ViewType2 targetGradEWeights,
      const ViewType1 basisGradAtTargetGradEPoints, const ViewType1 wBasisAtTargetGradEPoints, const ViewType3 computedDofs, const ViewType4 tagToOrdinal,
      const ViewType1 targetGradTanAtTargetGradEPoints, const ViewType5 targetGradAtTargetGradEPoints, const ViewType1 refEdgeTan,
      ordinal_type edgeCardinality, ordinal_type offsetBasis,
      ordinal_type offsetTarget, ordinal_type numVertexDofs, ordinal_type edgeDim, ordinal_type dim, ordinal_type iedge) :
        basisCoeffs_(basisCoeffs), negPartialProjGrad_(negPartialProjGrad), basisTanAtEPoints_(basisTanAtEPoints),
        basisGradAtBasisGradEPoints_(basisGradAtBasisGradEPoints), basisGradEWeights_(basisGradEWeights), wBasisAtBasisGradEPoints_(wBasisAtBasisGradEPoints), targetGradEWeights_(targetGradEWeights),
        basisGradAtTargetGradEPoints_(basisGradAtTargetGradEPoints), wBasisAtTargetGradEPoints_(wBasisAtTargetGradEPoints),
        computedDofs_(computedDofs), tagToOrdinal_(tagToOrdinal), targetGradTanAtTargetGradEPoints_(targetGradTanAtTargetGradEPoints),
        targetGradAtTargetGradEPoints_(targetGradAtTargetGradEPoints), refEdgeTan_(refEdgeTan),
        edgeCardinality_(edgeCardinality), offsetBasis_(offsetBasis),
        offsetTarget_(offsetTarget), numVertexDofs_(numVertexDofs), edgeDim_(edgeDim), dim_(dim), iedge_(iedge)
  {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    ordinal_type numBasisGradEPoints = basisGradEWeights_.extent(0);
    ordinal_type numTargetGradEPoints = targetGradEWeights_.extent(0);

    //this loop could be moved outside of cell loop
    for(ordinal_type j=0; j <edgeCardinality_; ++j) {
      ordinal_type jdof = tagToOrdinal_(edgeDim_, iedge_, j);

      for(ordinal_type iq=0; iq <numBasisGradEPoints; ++iq) {
        typename ViewType1::value_type tmp =0;
        for(ordinal_type d=0; d <dim_; ++d)
          tmp += basisGradAtBasisGradEPoints_(jdof,offsetBasis_+iq,d)*refEdgeTan_(d);
        basisTanAtEPoints_(0,j,iq) = tmp;
        wBasisAtBasisGradEPoints_(ic,j,iq) = tmp*basisGradEWeights_(iq);
      }
      for(ordinal_type iq=0; iq <numTargetGradEPoints; ++iq) {
        for(ordinal_type d=0; d <dim_; ++d)
          wBasisAtTargetGradEPoints_(ic,j,iq) += basisGradAtTargetGradEPoints_(jdof,offsetTarget_+iq,d)*refEdgeTan_(d)*targetGradEWeights_(iq);
      }
    }

    for(ordinal_type iq=0; iq <numTargetGradEPoints; ++iq)
      for(ordinal_type d=0; d <dim_; ++d)
        targetGradTanAtTargetGradEPoints_(ic,iq) += targetGradAtTargetGradEPoints_(ic,offsetTarget_+iq,d)*refEdgeTan_(d);

    for(ordinal_type j=0; j <numVertexDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <numBasisGradEPoints; ++iq)
        for(ordinal_type d=0; d <dim_; ++d)
          negPartialProjGrad_(ic,iq) -=  basisCoeffs_(ic,jdof)*basisGradAtBasisGradEPoints_(jdof,offsetBasis_+iq,d)*refEdgeTan_(d);
    }
  }
};

template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5>
struct ComputeBasisCoeffsOnFaces_HGRAD {
  const ViewType1 basisCoeffs_;
  const ViewType1 negPartialProjGrad_;
  const ViewType1 faceBasisGradAtGradEPoints_;
  const ViewType1 basisGradAtBasisGradEPoints_;
  const ViewType2 basisGradEWeights_;
  const ViewType1 wBasisGradAtGradEPoints_;
  const ViewType2 targetGradEWeights_;
  const ViewType1 basisGradAtTargetGradEPoints_;
  const ViewType1 wBasisGradBasisAtTargetGradEPoints_;
  const ViewType3 computedDofs_;
  const ViewType4 tagToOrdinal_;
  const ViewType1 targetGradTanAtTargetGradEPoints_;
  const ViewType5 targetGradAtTargetGradEPoints_;
  const ViewType1 refFaceNormal_;
  ordinal_type faceCardinality_;
  ordinal_type offsetBasisGrad_;
  ordinal_type offsetTargetGrad_;
  ordinal_type numVertexEdgeDofs_;
  ordinal_type numFaces_;
  ordinal_type faceDim_;
  ordinal_type dim_;
  ordinal_type iface_;

  ComputeBasisCoeffsOnFaces_HGRAD(const ViewType1 basisCoeffs, ViewType1 negPartialProjGrad,  const ViewType1 faceBasisGradAtGradEPoints,
      const ViewType1 basisGradAtBasisGradEPoints, const ViewType2 basisGradEWeights,  const ViewType1 wBasisGradAtGradEPoints,   const ViewType2 targetGradEWeights,
      const ViewType1 basisGradAtTargetGradEPoints, const ViewType1 wBasisGradBasisAtTargetGradEPoints, const ViewType3 computedDofs, const ViewType4 tagToOrdinal,
      const ViewType1 targetGradTanAtTargetGradEPoints, const ViewType5 targetGradAtTargetGradEPoints,
      const ViewType1 refFaceNormal, ordinal_type faceCardinality, ordinal_type offsetBasisGrad,
      ordinal_type offsetTargetGrad, ordinal_type numVertexEdgeDofs, ordinal_type numFaces, ordinal_type faceDim,
      ordinal_type dim, ordinal_type iface) :
        basisCoeffs_(basisCoeffs), negPartialProjGrad_(negPartialProjGrad), faceBasisGradAtGradEPoints_(faceBasisGradAtGradEPoints),
        basisGradAtBasisGradEPoints_(basisGradAtBasisGradEPoints), basisGradEWeights_(basisGradEWeights), wBasisGradAtGradEPoints_(wBasisGradAtGradEPoints),
        targetGradEWeights_(targetGradEWeights),
        basisGradAtTargetGradEPoints_(basisGradAtTargetGradEPoints), wBasisGradBasisAtTargetGradEPoints_(wBasisGradBasisAtTargetGradEPoints),
        computedDofs_(computedDofs), tagToOrdinal_(tagToOrdinal), targetGradTanAtTargetGradEPoints_(targetGradTanAtTargetGradEPoints),
        targetGradAtTargetGradEPoints_(targetGradAtTargetGradEPoints), refFaceNormal_(refFaceNormal),
        faceCardinality_(faceCardinality), offsetBasisGrad_(offsetBasisGrad),
        offsetTargetGrad_(offsetTargetGrad), numVertexEdgeDofs_(numVertexEdgeDofs), numFaces_(numFaces),
        faceDim_(faceDim), dim_(dim), iface_(iface)
  {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    ordinal_type numBasisGradEPoints = basisGradEWeights_.extent(0);
    ordinal_type numTargetGradEPoints = targetGradEWeights_.extent(0);

    //normal
    typename ViewType2::value_type n[3] = {refFaceNormal_(0), refFaceNormal_(1), refFaceNormal_(2)};

    for(ordinal_type j=0; j <faceCardinality_; ++j) {
      ordinal_type jdof = tagToOrdinal_(faceDim_, iface_, j);
      for(ordinal_type d=0; d <dim_; ++d) {
        ordinal_type dp1 = (d+1) % dim_;
        ordinal_type dp2 = (d+2) % dim_;
        for(ordinal_type iq=0; iq <numBasisGradEPoints; ++iq) {
          faceBasisGradAtGradEPoints_(0,j,iq,d) = basisGradAtBasisGradEPoints_(jdof,offsetBasisGrad_+iq,dp1)*n[dp2] - basisGradAtBasisGradEPoints_(jdof,offsetBasisGrad_+iq,dp2)*n[dp1];        
          wBasisGradAtGradEPoints_(ic,j,iq,d) = faceBasisGradAtGradEPoints_(0,j,iq,d) * basisGradEWeights_(iq);
        }
        for(ordinal_type iq=0; iq <numTargetGradEPoints; ++iq) {
          wBasisGradBasisAtTargetGradEPoints_(ic,j,iq,d) = (basisGradAtTargetGradEPoints_(jdof,offsetTargetGrad_+iq,dp1)*n[dp2] - basisGradAtTargetGradEPoints_(jdof,offsetTargetGrad_+iq,dp2)*n[dp1]) * targetGradEWeights_(iq);
        }
      }
    }

    for(ordinal_type d=0; d <dim_; ++d) {
      ordinal_type dp1 = (d+1) % dim_;
      ordinal_type dp2 = (d+2) % dim_;
      for(ordinal_type iq=0; iq <numTargetGradEPoints; ++iq)
        targetGradTanAtTargetGradEPoints_(ic,iq,d) = targetGradAtTargetGradEPoints_(ic,offsetTargetGrad_+iq,dp1)*n[dp2] - targetGradAtTargetGradEPoints_(ic,offsetTargetGrad_+iq,dp2)*n[dp1];
    }

    for(ordinal_type j=0; j <numVertexEdgeDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type d=0; d <dim_; ++d) {
        ordinal_type dp1 = (d+1) % dim_;
        ordinal_type dp2 = (d+2) % dim_;
        for(ordinal_type iq=0; iq <numBasisGradEPoints; ++iq)
          negPartialProjGrad_(ic,iq,d) -= (basisGradAtBasisGradEPoints_(jdof,offsetBasisGrad_+iq,dp1)*n[dp2] - basisGradAtBasisGradEPoints_(jdof,offsetBasisGrad_+iq,dp2)*n[dp1]) * basisCoeffs_(ic,jdof);
      }
    }
  }
};


template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4>
struct ComputeBasisCoeffsOnCells_HGRAD {
  const ViewType1 basisCoeffs_;
  const ViewType1 negPartialProjGrad_;
  const ViewType1 cellBasisGradAtGradEPoints_;
  const ViewType1 basisGradAtBasisGradEPoints_;
  const ViewType2 basisGradEWeights_;
  const ViewType1 wBasisGradAtGradEPoints_;
  const ViewType2 targetGradEWeights_;
  const ViewType1 basisGradAtTargetGradEPoints_;
  const ViewType1 wBasisGradBasisAtTargetGradEPoints_;
  const ViewType3 computedDofs_;
  const ViewType4 elemDof_;
  ordinal_type dim_;
  ordinal_type numElemDofs_;
  ordinal_type offsetBasisGrad_;
  ordinal_type offsetTargetGrad_;
  ordinal_type numVertexEdgeFaceDofs_;

  ComputeBasisCoeffsOnCells_HGRAD(const ViewType1 basisCoeffs, ViewType1 negPartialProjGrad,  const ViewType1 cellBasisGradAtGradEPoints,
      const ViewType1 basisGradAtBasisGradEPoints, const ViewType2 basisGradEWeights,  const ViewType1 wBasisGradAtGradEPoints,   const ViewType2 targetGradEWeights,
      const ViewType1 basisGradAtTargetGradEPoints, const ViewType1 wBasisGradBasisAtTargetGradEPoints, const ViewType3 computedDofs, const ViewType4 elemDof,
      ordinal_type dim, ordinal_type numElemDofs, ordinal_type offsetBasisGrad, ordinal_type offsetTargetGrad, ordinal_type numVertexEdgeFaceDofs) :
        basisCoeffs_(basisCoeffs), negPartialProjGrad_(negPartialProjGrad), cellBasisGradAtGradEPoints_(cellBasisGradAtGradEPoints),
        basisGradAtBasisGradEPoints_(basisGradAtBasisGradEPoints), basisGradEWeights_(basisGradEWeights), wBasisGradAtGradEPoints_(wBasisGradAtGradEPoints), targetGradEWeights_(targetGradEWeights),
        basisGradAtTargetGradEPoints_(basisGradAtTargetGradEPoints), wBasisGradBasisAtTargetGradEPoints_(wBasisGradBasisAtTargetGradEPoints),
        computedDofs_(computedDofs), elemDof_(elemDof), dim_(dim), numElemDofs_(numElemDofs), offsetBasisGrad_(offsetBasisGrad),
        offsetTargetGrad_(offsetTargetGrad), numVertexEdgeFaceDofs_(numVertexEdgeFaceDofs) {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {
    ordinal_type numBasisGradEPoints = basisGradEWeights_.extent(0);
    ordinal_type numTargetGradEPoints = targetGradEWeights_.extent(0);
    for(ordinal_type j=0; j <numElemDofs_; ++j) {
      ordinal_type idof = elemDof_(j);
      for(ordinal_type d=0; d <dim_; ++d) {
        for(ordinal_type iq=0; iq <numBasisGradEPoints; ++iq) {
          cellBasisGradAtGradEPoints_(0,j,iq,d) = basisGradAtBasisGradEPoints_(idof,offsetBasisGrad_+iq,d);
          wBasisGradAtGradEPoints_(ic,j,iq,d) = cellBasisGradAtGradEPoints_(0,j,iq,d) * basisGradEWeights_(iq);
        }
        for(ordinal_type iq=0; iq <numTargetGradEPoints; ++iq) {
          wBasisGradBasisAtTargetGradEPoints_(ic,j,iq,d )= basisGradAtTargetGradEPoints_(idof,offsetTargetGrad_+iq,d) * targetGradEWeights_(iq);
        }
      }
    }
    for(ordinal_type j=0; j <numVertexEdgeFaceDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <numBasisGradEPoints; ++iq)
        for(ordinal_type d=0; d <dim_; ++d) {
          negPartialProjGrad_(ic,iq,d) -=  basisCoeffs_(ic,jdof)*basisGradAtBasisGradEPoints_(jdof,offsetBasisGrad_+iq,d);
        }
    }
  }
};

} // FunctorsProjectionTools namespace


template<typename DeviceType>
template<class BasisCoeffsViewType, class TargetValueViewType, class TargetGradViewType, class BasisType, class OrientationViewType>
void
ProjectionTools<DeviceType>::getHGradBasisCoeffs(BasisCoeffsViewType basisCoeffs,
                                          const TargetValueViewType targetAtTargetEPoints,
                                          const TargetGradViewType targetGradAtTargetGradEPoints,
                                          const OrientationViewType orts,
                                          const BasisType* cellBasis,
                                          ProjectionStruct<DeviceType, typename BasisType::scalarType> * projStruct)
{
  using funValsValueType = typename TargetValueViewType::value_type;
  static_assert(std::is_same<funValsValueType,typename TargetGradViewType::value_type>::value,
                "targetGradAtTargetGradEPoints and targetAtTargetEPoints must agree on their value type" );
  
  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,DeviceType> ScalarViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type basisCardinality = cellBasis->getCardinality();
  ordinal_type numCells = targetAtTargetEPoints.extent(0);
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;

  ordinal_type numVertices = (cellBasis->getDofCount(0, 0) > 0) ? cellTopo.getVertexCount() : 0;
  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;

  ScalarViewType refEdgeTan("refTan", dim);
  ScalarViewType refFaceNormal("refNormal", dim);

  ordinal_type numVertexDofs = numVertices;

  ordinal_type numEdgeDofs(0);
  for(ordinal_type ie=0; ie<numEdges; ++ie)
    numEdgeDofs += cellBasis->getDofCount(edgeDim,ie);

  ordinal_type numFaceDofs(0);
  for(ordinal_type iface=0; iface<numFaces; ++iface)
    numFaceDofs += cellBasis->getDofCount(faceDim,iface);

  Kokkos::View<ordinal_type*, DeviceType> computedDofs("computedDofs",numVertexDofs+numEdgeDofs+numFaceDofs);

  ordinal_type numTotalBasisGradEPoints = projStruct->getNumBasisDerivEvalPoints();
  auto basisGradEPoints = projStruct->getAllDerivEvalPoints(EvalPointsType::BASIS);

  ordinal_type numTotalTargetEPoints = projStruct->getNumTargetEvalPoints(), numTotalTargetGradEPoints = projStruct->getNumTargetDerivEvalPoints();
  auto targetEPoints = projStruct->getAllEvalPoints(EvalPointsType::TARGET);
  auto targetGradEPoints = projStruct->getAllDerivEvalPoints(EvalPointsType::TARGET);


  ScalarViewType basisAtTargetEPoints("basisAtTargetEPoints",basisCardinality, numTotalTargetEPoints);
  cellBasis->getValues(basisAtTargetEPoints, targetEPoints);

  ScalarViewType basisGradAtBasisGradEPoints;
  ScalarViewType basisGradAtTargetGradEPoints;
  if(numTotalBasisGradEPoints>0) {
    basisGradAtBasisGradEPoints = ScalarViewType ("basisGradAtBasisGradEPoints",basisCardinality, numTotalBasisGradEPoints, dim);
    basisGradAtTargetGradEPoints = ScalarViewType("basisGradAtTargetGradEPoints",basisCardinality, numTotalTargetGradEPoints, dim);

    cellBasis->getValues(basisGradAtBasisGradEPoints, basisGradEPoints, OPERATOR_GRAD);
    cellBasis->getValues(basisGradAtTargetGradEPoints, targetGradEPoints, OPERATOR_GRAD);
  }

  auto targetEPointsRange = Kokkos::create_mirror_view_and_copy(MemSpaceType(), projStruct->getTargetPointsRange());
  auto targetGradEPointsRange  = projStruct->getTargetDerivPointsRange();
  auto basisGradEPointsRange  = projStruct->getBasisDerivPointsRange();

  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(MemSpaceType(), cellBasis->getAllDofOrdinal());

  typename RefSubcellParametrization<DeviceType>::ConstViewType  subcellParamFace;
  if(numFaces>0)
    subcellParamFace = RefSubcellParametrization<DeviceType>::get(faceDim, cellBasis->getBaseCellTopology().getKey());


  Kokkos::parallel_for("Compute Dofs ", Kokkos::RangePolicy<ExecSpaceType, int> (0, numVertices),
      KOKKOS_LAMBDA (const size_t iv) {
    computedDofs(iv) = tagToOrdinal(0, iv, 0);
  });
  ordinal_type computedDofsCount = numVertices;

  ScalarViewType refBasisCoeffs("refBasisCoeffs", basisCoeffs.extent(0), basisCoeffs.extent(1));

  const Kokkos::RangePolicy<ExecSpaceType> policy(0, numCells);
  using functorType = FunctorsProjectionTools::ComputeBasisCoeffsOnVertices_HGRAD<ScalarViewType, decltype(tagToOrdinal), decltype(targetEPointsRange),
      decltype(targetAtTargetEPoints), decltype(basisAtTargetEPoints)>;
  Kokkos::parallel_for(policy, functorType(refBasisCoeffs, tagToOrdinal, targetEPointsRange,
      targetAtTargetEPoints, basisAtTargetEPoints, numVertices));
  
  for(ordinal_type ie=0; ie<numEdges; ++ie)  {

    ordinal_type edgeCardinality = cellBasis->getDofCount(edgeDim,ie);
    ordinal_type offsetBasis = basisGradEPointsRange(edgeDim, ie).first;
    ordinal_type offsetTarget = targetGradEPointsRange(edgeDim, ie).first;
    ordinal_type numBasisGradEPoints = range_size(basisGradEPointsRange(edgeDim, ie));
    ordinal_type numTargetGradEPoints = range_size(targetGradEPointsRange(edgeDim, ie));
    auto basisGradEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisDerivEvalWeights(edgeDim,ie));
    auto targetGradEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetDerivEvalWeights(edgeDim,ie));

    CellTools<DeviceType>::getReferenceEdgeTangent(refEdgeTan, ie, cellTopo);

    ScalarViewType basisTanAtEPoints("basisTanAtEPoints",1,edgeCardinality, numBasisGradEPoints);
    ScalarViewType targetGradTanAtTargetGradEPoints("tanBasisAtTargetGradEPoints",numCells, numTargetGradEPoints);
    ScalarViewType wBasisAtBasisGradEPoints("wTanBasisAtBasisGradEPoints",numCells,edgeCardinality, numBasisGradEPoints);
    ScalarViewType wBasisAtTargetGradEPoints("wTanBasisAtTargetGradEPoints",numCells,edgeCardinality, numTargetGradEPoints);
    ScalarViewType negPartialProjGrad("negPartialProjGrad", numCells, numBasisGradEPoints);

    using functorTypeEdge = FunctorsProjectionTools::ComputeBasisCoeffsOnEdges_HGRAD<ScalarViewType,  decltype(basisGradEWeights),
        decltype(computedDofs), decltype(tagToOrdinal), decltype(targetGradAtTargetGradEPoints)>;

    Kokkos::parallel_for(policy, functorTypeEdge(refBasisCoeffs, negPartialProjGrad,  basisTanAtEPoints,
        basisGradAtBasisGradEPoints, basisGradEWeights,  wBasisAtBasisGradEPoints,   targetGradEWeights,
        basisGradAtTargetGradEPoints, wBasisAtTargetGradEPoints, computedDofs, tagToOrdinal,
        targetGradTanAtTargetGradEPoints, targetGradAtTargetGradEPoints, refEdgeTan,
        edgeCardinality, offsetBasis,
        offsetTarget, numVertexDofs, edgeDim, dim, ie));

    ScalarViewType edgeMassMat_("edgeMassMat_", 1, edgeCardinality, edgeCardinality),
        edgeRhsMat_("rhsMat_", numCells, edgeCardinality);

    FunctionSpaceTools<DeviceType >::integrate(edgeMassMat_, basisTanAtEPoints, Kokkos::subview(wBasisAtBasisGradEPoints, std::make_pair(0,1), Kokkos::ALL(), Kokkos::ALL()) );
    FunctionSpaceTools<DeviceType >::integrate(edgeRhsMat_, targetGradTanAtTargetGradEPoints, wBasisAtTargetGradEPoints);
    FunctionSpaceTools<DeviceType >::integrate(edgeRhsMat_, negPartialProjGrad, wBasisAtBasisGradEPoints,true);

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, DeviceType> WorkArrayViewType;
    ScalarViewType t_("t",numCells, edgeCardinality);
    WorkArrayViewType w_("w",numCells, edgeCardinality);

    auto edgeDofs = Kokkos::subview(tagToOrdinal, edgeDim, ie, Kokkos::ALL());

    ElemSystem edgeSystem("edgeSystem", true);
    edgeSystem.solve(refBasisCoeffs, edgeMassMat_, edgeRhsMat_, t_, w_, edgeDofs, edgeCardinality);

    Kokkos::parallel_for("Compute Dofs ", Kokkos::RangePolicy<ExecSpaceType, int> (0, edgeCardinality),
        KOKKOS_LAMBDA (const size_t i) {
      computedDofs(computedDofsCount+i) = tagToOrdinal(edgeDim, ie, i);
    });
    computedDofsCount += edgeCardinality;

  }

  for(ordinal_type iface=0; iface<numFaces; ++iface) {

    ordinal_type faceCardinality = cellBasis->getDofCount(faceDim,iface);

    ordinal_type numGradEPoints = range_size(basisGradEPointsRange(faceDim, iface));
    ordinal_type numTargetGradEPoints = range_size(targetGradEPointsRange(faceDim, iface));
    ordinal_type offsetBasisGrad = basisGradEPointsRange(faceDim, iface).first;
    ordinal_type offsetTargetGrad = targetGradEPointsRange(faceDim, iface).first;
    auto basisGradEWeights =  Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisDerivEvalWeights(faceDim,iface));
    auto targetGradEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetDerivEvalWeights(faceDim,iface));

    CellTools<DeviceType>::getReferenceFaceNormal(refFaceNormal, iface, cellTopo);

    ScalarViewType faceBasisGradAtGradEPoints("faceBasisGradAtGradEPoints",1,faceCardinality, numGradEPoints,dim);
    ScalarViewType wBasisGradAtGradEPoints("wNormalBasisGradAtGradEPoints",numCells,faceCardinality, numGradEPoints,dim);
    ScalarViewType wBasisGradBasisAtTargetGradEPoints("wNormalBasisGradAtTargetGradEPoints",numCells,faceCardinality, numTargetGradEPoints,dim);
    ScalarViewType targetGradTanAtTargetGradEPoints("targetGradTanAtTargetGradEPoints",numCells, numTargetGradEPoints,dim);
    ScalarViewType negPartialProjGrad("mNormalComputedProjection", numCells,numGradEPoints,dim);

    using functorTypeFace_HGRAD = FunctorsProjectionTools::ComputeBasisCoeffsOnFaces_HGRAD<ScalarViewType,  decltype(basisGradEWeights),
        decltype(computedDofs), decltype(tagToOrdinal), decltype(targetGradAtTargetGradEPoints)>;

    Kokkos::parallel_for(policy, functorTypeFace_HGRAD(refBasisCoeffs, negPartialProjGrad, faceBasisGradAtGradEPoints,
        basisGradAtBasisGradEPoints, basisGradEWeights, wBasisGradAtGradEPoints, targetGradEWeights,
        basisGradAtTargetGradEPoints,wBasisGradBasisAtTargetGradEPoints, computedDofs, tagToOrdinal,
        targetGradTanAtTargetGradEPoints,targetGradAtTargetGradEPoints,
        refFaceNormal, faceCardinality, offsetBasisGrad,
        offsetTargetGrad, numVertexDofs+numEdgeDofs, numFaces, faceDim,
        dim, iface));

    ScalarViewType faceMassMat_("faceMassMat_", 1, faceCardinality, faceCardinality),
        faceRhsMat_("rhsMat_", numCells, faceCardinality);

    FunctionSpaceTools<DeviceType >::integrate(faceMassMat_, faceBasisGradAtGradEPoints, Kokkos::subview(wBasisGradAtGradEPoints, std::make_pair(0,1), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) );
    FunctionSpaceTools<DeviceType >::integrate(faceRhsMat_, targetGradTanAtTargetGradEPoints, wBasisGradBasisAtTargetGradEPoints);
    FunctionSpaceTools<DeviceType >::integrate(faceRhsMat_, negPartialProjGrad, wBasisGradAtGradEPoints,true);

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, DeviceType> WorkArrayViewType;
    ScalarViewType t_("t",numCells, faceCardinality);
    WorkArrayViewType w_("w",numCells, faceCardinality);

    auto faceDofs = Kokkos::subview(tagToOrdinal, faceDim, iface, Kokkos::ALL());

    ElemSystem faceSystem("faceSystem", true);
    faceSystem.solve(refBasisCoeffs, faceMassMat_, faceRhsMat_, t_, w_, faceDofs, faceCardinality);

    Kokkos::parallel_for("Compute Face Dofs ", Kokkos::RangePolicy<ExecSpaceType, int> (0, faceCardinality),
        KOKKOS_LAMBDA (const size_t i) {
      computedDofs(computedDofsCount+i) = tagToOrdinal(faceDim, iface, i);
    });
    computedDofsCount += faceCardinality;
  }

  ordinal_type numElemDofs = cellBasis->getDofCount(dim,0);
  if(numElemDofs>0) {

    range_type cellTargetGradEPointsRange = targetGradEPointsRange(dim, 0);
    ordinal_type numTargetGradEPoints = range_size(cellTargetGradEPointsRange);
    ordinal_type numGradEPoints = range_size(basisGradEPointsRange(dim, 0));
    ordinal_type offsetBasisGrad = basisGradEPointsRange(dim, 0).first;
    ordinal_type offsetTargetGrad = cellTargetGradEPointsRange.first;
    auto targetGradEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetDerivEvalWeights(dim,0));
    auto basisGradEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisDerivEvalWeights(dim,0));

    ScalarViewType cellBasisGradAtGradEPoints("internalBasisGradAtEPoints",1,numElemDofs, numGradEPoints, dim);
    ScalarViewType negPartialProjGrad("negPartialProjGrad", numCells, numGradEPoints, dim);
    ScalarViewType wBasisGradAtGradEPoints("wBasisGradAtGradEPoints",numCells,numElemDofs, numGradEPoints,dim);
    ScalarViewType wBasisGradBasisAtTargetGradEPoints("wBasisGradAtTargetGradEPoints",numCells,numElemDofs, numTargetGradEPoints,dim);

    auto elemDof = Kokkos::subview(tagToOrdinal, dim, 0, Kokkos::ALL());
    using functorTypeCell_HGRAD = FunctorsProjectionTools::ComputeBasisCoeffsOnCells_HGRAD<ScalarViewType,  decltype(basisGradEWeights), decltype(computedDofs), decltype(elemDof)>;
    Kokkos::parallel_for(policy, functorTypeCell_HGRAD(refBasisCoeffs, negPartialProjGrad,  cellBasisGradAtGradEPoints,
        basisGradAtBasisGradEPoints, basisGradEWeights,  wBasisGradAtGradEPoints,   targetGradEWeights,
        basisGradAtTargetGradEPoints, wBasisGradBasisAtTargetGradEPoints, computedDofs, elemDof,
        dim, numElemDofs, offsetBasisGrad, offsetTargetGrad, numVertexDofs+numEdgeDofs+numFaceDofs));
    ScalarViewType cellMassMat_("cellMassMat_", 1, numElemDofs, numElemDofs),
        cellRhsMat_("rhsMat_", numCells, numElemDofs);

    FunctionSpaceTools<DeviceType >::integrate(cellMassMat_, cellBasisGradAtGradEPoints, Kokkos::subview(wBasisGradAtGradEPoints, std::make_pair(0,1), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) );
    FunctionSpaceTools<DeviceType >::integrate(cellRhsMat_, Kokkos::subview(targetGradAtTargetGradEPoints,Kokkos::ALL(),cellTargetGradEPointsRange,Kokkos::ALL()), wBasisGradBasisAtTargetGradEPoints);
    FunctionSpaceTools<DeviceType >::integrate(cellRhsMat_, negPartialProjGrad, wBasisGradAtGradEPoints, true);

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, DeviceType> WorkArrayViewType;
    ScalarViewType t_("t",numCells, numElemDofs);
    WorkArrayViewType w_("w",numCells, numElemDofs);

    auto cellDofs = Kokkos::subview(tagToOrdinal, dim, 0, Kokkos::ALL());
    ElemSystem cellSystem("cellSystem", true);
    cellSystem.solve(refBasisCoeffs, cellMassMat_, cellRhsMat_, t_, w_, cellDofs, numElemDofs);
  }

  OrientationTools<DeviceType>::modifyBasisByOrientationInverse(basisCoeffs, refBasisCoeffs, orts, cellBasis, true);

}

}  // Intrepid2 namespace

#endif

