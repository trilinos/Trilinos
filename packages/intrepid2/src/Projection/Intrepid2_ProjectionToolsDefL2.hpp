// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ProjectionToolsDefL2.hpp
    \brief  Header file for the Intrepid2::ProjectionTools
            containing definitions for L2 projections.
    \author Created by Mauro Perego
 */

#ifndef __INTREPID2_PROJECTIONTOOLSDEFL2_HPP__
#define __INTREPID2_PROJECTIONTOOLSDEFL2_HPP__

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


namespace Intrepid2 {

namespace FunctorsProjectionTools {

template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4>
struct ComputeBasisCoeffsOnVertices_L2 {
  ViewType1 basisCoeffs_;
  const ViewType2 tagToOrdinal_;
  const ViewType3 targetEPointsRange_;
  const ViewType4 targetAtTargetEPoints_;
  const ViewType1 basisAtTargetEPoints_;
  ordinal_type numVertices_;


  ComputeBasisCoeffsOnVertices_L2(ViewType1 basisCoeffs, ViewType2 tagToOrdinal, ViewType3 targetEPointsRange,
      ViewType4  targetAtTargetEPoints, ViewType1 basisAtTargetEPoints, ordinal_type numVertices) :
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
struct ComputeBasisCoeffsOnEdges_L2 {
  const ViewType1 basisCoeffs_;
  const ViewType1 negPartialProj_;
  const ViewType1 basisDofDofAtBasisEPoints_;
  const ViewType1 basisAtBasisEPoints_;
  const ViewType2 basisEWeights_;
  const ViewType1 wBasisDofAtBasisEPoints_;
  const ViewType2 targetEWeights_;
  const ViewType1 basisAtTargetEPoints_;
  const ViewType1 wBasisDofAtTargetEPoints_;
  const ViewType3 computedDofs_;
  const ViewType4 tagToOrdinal_;
  const ViewType5 targetAtTargetEPoints_;
  const ViewType1 targetTanAtTargetEPoints_;
  const ViewType1 refEdgesVec_;
  ordinal_type fieldDim_;
  ordinal_type edgeCardinality_;
  ordinal_type offsetBasis_;
  ordinal_type offsetTarget_;
  ordinal_type numVertexDofs_;
  ordinal_type edgeDim_;
  ordinal_type iedge_;

  ComputeBasisCoeffsOnEdges_L2(const ViewType1 basisCoeffs, ViewType1 negPartialProj,  const ViewType1 basisDofDofAtBasisEPoints,
      const ViewType1 basisAtBasisEPoints, const ViewType2 basisEWeights,  const ViewType1 wBasisDofAtBasisEPoints,   const ViewType2 targetEWeights,
      const ViewType1 basisAtTargetEPoints, const ViewType1 wBasisDofAtTargetEPoints, const ViewType3 computedDofs, const ViewType4 tagToOrdinal,
      const ViewType5 targetAtTargetEPoints, const ViewType1 targetTanAtTargetEPoints, const ViewType1 refEdgesVec,
      ordinal_type fieldDim, ordinal_type edgeCardinality, ordinal_type offsetBasis,
      ordinal_type offsetTarget, ordinal_type numVertexDofs, ordinal_type edgeDim, ordinal_type iedge) :
        basisCoeffs_(basisCoeffs), negPartialProj_(negPartialProj), basisDofDofAtBasisEPoints_(basisDofDofAtBasisEPoints),
        basisAtBasisEPoints_(basisAtBasisEPoints), basisEWeights_(basisEWeights), wBasisDofAtBasisEPoints_(wBasisDofAtBasisEPoints), targetEWeights_(targetEWeights),
        basisAtTargetEPoints_(basisAtTargetEPoints), wBasisDofAtTargetEPoints_(wBasisDofAtTargetEPoints),
        computedDofs_(computedDofs), tagToOrdinal_(tagToOrdinal), targetAtTargetEPoints_(targetAtTargetEPoints),
        targetTanAtTargetEPoints_(targetTanAtTargetEPoints), refEdgesVec_(refEdgesVec),
        fieldDim_(fieldDim), edgeCardinality_(edgeCardinality), offsetBasis_(offsetBasis),
        offsetTarget_(offsetTarget), numVertexDofs_(numVertexDofs), edgeDim_(edgeDim), iedge_(iedge)
  {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {
    for(ordinal_type j=0; j <edgeCardinality_; ++j) {
      ordinal_type jdof =tagToOrdinal_(edgeDim_, iedge_, j);
      for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq) {
        typename ViewType1::value_type tmp =0;
        for(ordinal_type d=0; d <fieldDim_; ++d)
          tmp += basisAtBasisEPoints_(jdof,offsetBasis_+iq,d)*refEdgesVec_(iedge_,d);
        basisDofDofAtBasisEPoints_(0,j,iq) = tmp;
        wBasisDofAtBasisEPoints_(ic,j,iq) = tmp*basisEWeights_(iq);
      }
      for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {
        for(ordinal_type d=0; d <fieldDim_; ++d)
          wBasisDofAtTargetEPoints_(ic,j,iq) += basisAtTargetEPoints_(jdof,offsetTarget_+iq,d)*refEdgesVec_(iedge_,d)*targetEWeights_(iq);
      }
    }

    for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq)
      for(ordinal_type d=0; d <fieldDim_; ++d)
        targetTanAtTargetEPoints_(ic,iq) += targetAtTargetEPoints_(ic,offsetTarget_+iq,d)*refEdgesVec_(iedge_,d);

    for(ordinal_type j=0; j <numVertexDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq)
        for(ordinal_type d=0; d <fieldDim_; ++d)
          negPartialProj_(ic,iq) -=  basisCoeffs_(ic,jdof)*basisAtBasisEPoints_(jdof,offsetBasis_+iq,d)*refEdgesVec_(iedge_,d);
    }
  }
};

template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5>
struct ComputeBasisCoeffsOnFaces_L2 {
  const ViewType1 basisCoeffs_;
  const ViewType1 negPartialProj_;
  const ViewType1 faceBasisDofAtBasisEPoints_;
  const ViewType1 basisAtBasisEPoints_;
  const ViewType2 basisEWeights_;
  const ViewType1 wBasisDofAtBasisEPoints_;
  const ViewType2 targetEWeights_;
  const ViewType1 basisAtTargetEPoints_;
  const ViewType1 wBasisDofAtTargetEPoints_;
  const ViewType3 computedDofs_;
  const ViewType4 tagToOrdinal_;
  const ViewType5 targetAtTargetEPoints_;
  const ViewType1 targetDofAtTargetEPoints_;
  const ViewType1 refSideNormal_;
  ordinal_type fieldDim_;
  ordinal_type faceCardinality_;
  ordinal_type offsetBasis_;
  ordinal_type offsetTarget_;
  ordinal_type numVertexEdgeDofs_;
  ordinal_type numFaces_;
  ordinal_type faceDim_;
  ordinal_type dim_;
  ordinal_type iface_;
  bool isHCurlBasis_, isHDivBasis_;

  ComputeBasisCoeffsOnFaces_L2(const ViewType1 basisCoeffs, ViewType1 negPartialProj,  const ViewType1 faceBasisDofAtBasisEPoints,
      const ViewType1 basisAtBasisEPoints, const ViewType2 basisEWeights,  const ViewType1 wBasisDofAtBasisEPoints,   const ViewType2 targetEWeights,
      const ViewType1 basisAtTargetEPoints, const ViewType1 wBasisDofAtTargetEPoints, const ViewType3 computedDofs, const ViewType4 tagToOrdinal,
      const ViewType5 targetAtTargetEPoints, const ViewType1 targetDofAtTargetEPoints,
      const ViewType1 refSideNormal, ordinal_type fieldDim, ordinal_type faceCardinality, ordinal_type offsetBasis,
      ordinal_type offsetTarget, ordinal_type numVertexEdgeDofs, ordinal_type numFaces, ordinal_type faceDim,
      ordinal_type dim, ordinal_type iface, bool isHCurlBasis, bool isHDivBasis) :
        basisCoeffs_(basisCoeffs), negPartialProj_(negPartialProj), faceBasisDofAtBasisEPoints_(faceBasisDofAtBasisEPoints),
        basisAtBasisEPoints_(basisAtBasisEPoints), basisEWeights_(basisEWeights), wBasisDofAtBasisEPoints_(wBasisDofAtBasisEPoints), targetEWeights_(targetEWeights),
        basisAtTargetEPoints_(basisAtTargetEPoints), wBasisDofAtTargetEPoints_(wBasisDofAtTargetEPoints),
        computedDofs_(computedDofs), tagToOrdinal_(tagToOrdinal), targetAtTargetEPoints_(targetAtTargetEPoints),
        targetDofAtTargetEPoints_(targetDofAtTargetEPoints), refSideNormal_(refSideNormal),
        fieldDim_(fieldDim), faceCardinality_(faceCardinality), offsetBasis_(offsetBasis),
        offsetTarget_(offsetTarget), numVertexEdgeDofs_(numVertexEdgeDofs), numFaces_(numFaces),
        faceDim_(faceDim), dim_(dim), iface_(iface),
        isHCurlBasis_(isHCurlBasis), isHDivBasis_(isHDivBasis)
  {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    if(isHCurlBasis_) {
      //normal
      typename ViewType1::value_type n[3] = {refSideNormal_(0), refSideNormal_(1), refSideNormal_(2)};

      for(ordinal_type j=0; j <faceCardinality_; ++j) {
        ordinal_type jdof = tagToOrdinal_(faceDim_, iface_, j);
        for(ordinal_type d=0; d <fieldDim_; ++d) {
          ordinal_type dp1 = (d+1) % dim_;
          ordinal_type dp2 = (d+2) % dim_;
          for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq) {

            faceBasisDofAtBasisEPoints_(0,j,iq,d) = basisAtBasisEPoints_(jdof,offsetBasis_+iq,dp1)*n[dp2] - basisAtBasisEPoints_(jdof,offsetBasis_+iq,dp2)*n[dp1];
            // basis \times n
            wBasisDofAtBasisEPoints_(ic,j,iq,d) = faceBasisDofAtBasisEPoints_(0,j,iq,d) * basisEWeights_(iq);
          }
          for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {
            // basis \times n
            wBasisDofAtTargetEPoints_(ic,j,iq,d) = (basisAtTargetEPoints_(jdof,offsetTarget_+iq,dp1)*n[dp2] - basisAtTargetEPoints_(jdof,offsetTarget_+iq,dp2)*n[dp1]) * targetEWeights_(iq);
          }
        }
      }

      for(ordinal_type d=0; d <fieldDim_; ++d) {
        ordinal_type dp1 = (d+1) % dim_;
        ordinal_type dp2 = (d+2) % dim_;
        for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {     
          // target \times n
          targetDofAtTargetEPoints_(ic,iq,d) = targetAtTargetEPoints_(ic,offsetTarget_+iq,dp1)*n[dp2] - targetAtTargetEPoints_(ic,offsetTarget_+iq,dp2)*n[dp1];
        }
      }

      for(ordinal_type j=0; j <numVertexEdgeDofs_; ++j) {
        ordinal_type jdof = computedDofs_(j);
        for(ordinal_type d=0; d <fieldDim_; ++d) {
          ordinal_type dp1 = (d+1) % dim_;
          ordinal_type dp2 = (d+2) % dim_;
          for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq) {
            // basis \times n
            negPartialProj_(ic,iq,d) -= (basisAtBasisEPoints_(jdof,offsetBasis_+iq,dp1)*n[dp2] - basisAtBasisEPoints_(jdof,offsetBasis_+iq,dp2)*n[dp1])*basisCoeffs_(ic,jdof);
          } 
        }
      }
    } else { // isHDivBasis_ || isHGradBasis_
      typename ViewType1::value_type coeff[3] = {1,0,0}; //only need first component for HGrad basis
      if (isHDivBasis_) {
        for (ordinal_type d=0; d<3; d++)
          coeff[d] =  refSideNormal_(d);
      }

      for(ordinal_type j=0; j <faceCardinality_; ++j) {
        ordinal_type jdof = tagToOrdinal_(faceDim_, iface_, j);
        for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq) {
          typename ViewType1::value_type tmp =0;
          for(ordinal_type d=0; d <fieldDim_; ++d)
            tmp += coeff[d]*basisAtBasisEPoints_(jdof,offsetBasis_+iq,d);
          faceBasisDofAtBasisEPoints_(0,j,iq,0) = tmp;
          wBasisDofAtBasisEPoints_(ic,j,iq,0) = tmp * basisEWeights_(iq);
        }
        for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {
          typename ViewType2::value_type sum=0;
          for(ordinal_type d=0; d <fieldDim_; ++d)
            sum += coeff[d]*basisAtTargetEPoints_(jdof,offsetTarget_+iq,d);
          wBasisDofAtTargetEPoints_(ic,j,iq,0) = sum * targetEWeights_(iq);
        }
      }

      for(ordinal_type d=0; d <fieldDim_; ++d)
        for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq)
          targetDofAtTargetEPoints_(ic,iq,0) += coeff[d]*targetAtTargetEPoints_(ic,offsetTarget_+iq,d);

      for(ordinal_type j=0; j <numVertexEdgeDofs_; ++j) {
        ordinal_type jdof = computedDofs_(j);
        for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq)
          for(ordinal_type d=0; d <fieldDim_; ++d)
            negPartialProj_(ic,iq,0) -=  basisCoeffs_(ic,jdof)*coeff[d]*basisAtBasisEPoints_(jdof,offsetBasis_+iq,d);
      }
    }
  }
};


template<typename ViewType1, typename ViewType2, typename ViewType3, typename ViewType4>
struct ComputeBasisCoeffsOnCells_L2 {
  const ViewType1 basisCoeffs_;
  const ViewType1 negPartialProj_;
  const ViewType1 internalBasisAtBasisEPoints_;
  const ViewType1 basisAtBasisEPoints_;
  const ViewType2 basisEWeights_;
  const ViewType1 wBasisAtBasisEPoints_;
  const ViewType2 targetEWeights_;
  const ViewType1 basisAtTargetEPoints_;
  const ViewType1 wBasisDofAtTargetEPoints_;
  const ViewType3 computedDofs_;
  const ViewType4 elemDof_;
  ordinal_type fieldDim_;
  ordinal_type numElemDofs_;
  ordinal_type offsetBasis_;
  ordinal_type offsetTarget_;
  ordinal_type numVertexEdgeFaceDofs_;

  ComputeBasisCoeffsOnCells_L2(const ViewType1 basisCoeffs, ViewType1 negPartialProj,  const ViewType1 internalBasisAtBasisEPoints,
      const ViewType1 basisAtBasisEPoints, const ViewType2 basisEWeights,  const ViewType1 wBasisAtBasisEPoints,   const ViewType2 targetEWeights,
      const ViewType1 basisAtTargetEPoints, const ViewType1 wBasisDofAtTargetEPoints, const ViewType3 computedDofs, const ViewType4 elemDof,
      ordinal_type fieldDim, ordinal_type numElemDofs, ordinal_type offsetBasis, ordinal_type offsetTarget, ordinal_type numVertexEdgeFaceDofs) :
        basisCoeffs_(basisCoeffs), negPartialProj_(negPartialProj), internalBasisAtBasisEPoints_(internalBasisAtBasisEPoints),
        basisAtBasisEPoints_(basisAtBasisEPoints), basisEWeights_(basisEWeights), wBasisAtBasisEPoints_(wBasisAtBasisEPoints), targetEWeights_(targetEWeights),
        basisAtTargetEPoints_(basisAtTargetEPoints), wBasisDofAtTargetEPoints_(wBasisDofAtTargetEPoints),
        computedDofs_(computedDofs), elemDof_(elemDof), fieldDim_(fieldDim), numElemDofs_(numElemDofs), offsetBasis_(offsetBasis),
        offsetTarget_(offsetTarget), numVertexEdgeFaceDofs_(numVertexEdgeFaceDofs) {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    for(ordinal_type j=0; j <numElemDofs_; ++j) {
      ordinal_type idof = elemDof_(j);
      for(ordinal_type d=0; d <fieldDim_; ++d) {
        for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq) {
          internalBasisAtBasisEPoints_(0,j,iq,d) = basisAtBasisEPoints_(idof,offsetBasis_+iq,d);
          wBasisAtBasisEPoints_(ic,j,iq,d) = internalBasisAtBasisEPoints_(0,j,iq,d) * basisEWeights_(iq);
        }
        for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {
          wBasisDofAtTargetEPoints_(ic,j,iq,d) = basisAtTargetEPoints_(idof,offsetTarget_+iq,d)* targetEWeights_(iq);
        }
      }
    }
    for(ordinal_type j=0; j < numVertexEdgeFaceDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq)
        for(ordinal_type d=0; d <fieldDim_; ++d) {
          negPartialProj_(ic,iq,d) -=  basisCoeffs_(ic,jdof)*basisAtBasisEPoints_(jdof,offsetBasis_+iq,d);
        }
    }
  }
};

template<typename ViewType1, typename ViewType2>
struct MultiplyBasisByWeights {
  const ViewType1 basisAtBasisEPoints_;
  const ViewType2 basisEWeights_;
  const ViewType1 wBasisAtBasisEPoints_;
  const ViewType2 targetEWeights_;
  const ViewType1 basisAtTargetEPoints_;
  const ViewType1 wBasisDofAtTargetEPoints_;
  ordinal_type fieldDim_;
  ordinal_type numElemDofs_;

  MultiplyBasisByWeights(const ViewType1 basisAtBasisEPoints, const ViewType2 basisEWeights,  const ViewType1 wBasisAtBasisEPoints,   const ViewType2 targetEWeights,
      const ViewType1 basisAtTargetEPoints, const ViewType1 wBasisDofAtTargetEPoints,
      ordinal_type fieldDim, ordinal_type numElemDofs) :
        basisAtBasisEPoints_(basisAtBasisEPoints), basisEWeights_(basisEWeights), wBasisAtBasisEPoints_(wBasisAtBasisEPoints), targetEWeights_(targetEWeights),
        basisAtTargetEPoints_(basisAtTargetEPoints), wBasisDofAtTargetEPoints_(wBasisDofAtTargetEPoints),
        fieldDim_(fieldDim), numElemDofs_(numElemDofs) {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    for(ordinal_type j=0; j <numElemDofs_; ++j) {
      for(ordinal_type d=0; d <fieldDim_; ++d) {
        for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq) {
          wBasisAtBasisEPoints_(ic,j,iq,d) = basisAtBasisEPoints_(j,iq,d) * basisEWeights_(iq);
        }
        for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {
          wBasisDofAtTargetEPoints_(ic,j,iq,d) = basisAtTargetEPoints_(j,iq,d)* targetEWeights_(iq);
        }
      }
    }
  }
};

} // FunctorsProjectionTools namespace


template<typename DeviceType>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<DeviceType>::getL2BasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtTargetEPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<DeviceType, typename BasisType::scalarType> * projStruct){

  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,DeviceType> ScalarViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type basisCardinality = cellBasis->getCardinality();
  ordinal_type numCells = targetAtTargetEPoints.extent(0);
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;
  const ordinal_type fieldDim = (targetAtTargetEPoints.rank()==2) ? 1 : targetAtTargetEPoints.extent(2);

  ordinal_type numVertices = (cellBasis->getDofCount(0, 0) > 0) ? cellTopo.getVertexCount() : 0;
  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;

  ScalarViewType refEdgesVec("refEdgesVec", numEdges, dim);
  ScalarViewType refFacesTangents("refFaceTangents", numFaces, dim, 2);
  ScalarViewType refFacesNormal("refFaceNormal", numFaces, dim);

  ordinal_type numVertexDofs = numVertices;

  ordinal_type numEdgeDofs(0);
  for(ordinal_type ie=0; ie<numEdges; ++ie)
    numEdgeDofs += cellBasis->getDofCount(edgeDim,ie);

  ordinal_type numFaceDofs(0);
  for(ordinal_type iface=0; iface<numFaces; ++iface)
    numFaceDofs += cellBasis->getDofCount(faceDim,iface);

  Kokkos::View<ordinal_type*, DeviceType> computedDofs("computedDofs", numVertexDofs+numEdgeDofs+numFaceDofs);

  auto basisEPointsRange = projStruct->getBasisPointsRange();

  ordinal_type numTotalBasisEPoints = projStruct->getNumBasisEvalPoints();
  auto basisEPoints = projStruct->getAllEvalPoints(EvalPointsType::BASIS);

  ordinal_type numTotalTargetEPoints = projStruct->getNumTargetEvalPoints();
  auto targetEPoints = projStruct->getAllEvalPoints(EvalPointsType::TARGET);

  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(MemSpaceType(), cellBasis->getAllDofOrdinal());

  ScalarViewType basisAtBasisEPoints("basisAtBasisEPoints",basisCardinality, numTotalBasisEPoints, fieldDim);
  ScalarViewType basisAtTargetEPoints("basisAtTargetEPoints",basisCardinality, numTotalTargetEPoints, fieldDim);
  {
    if(fieldDim == 1) {
      cellBasis->getValues(Kokkos::subview(basisAtTargetEPoints, Kokkos::ALL(), Kokkos::ALL(), 0), targetEPoints);
      cellBasis->getValues(Kokkos::subview(basisAtBasisEPoints, Kokkos::ALL(), Kokkos::ALL(), 0), basisEPoints);
    }
    else {
      cellBasis->getValues(basisAtTargetEPoints, targetEPoints);
      cellBasis->getValues(basisAtBasisEPoints, basisEPoints);
    }
  }

  {
    auto hostComputedDofs = Kokkos::create_mirror_view(computedDofs);
    ordinal_type computedDofsCount = 0;
    for(ordinal_type iv=0; iv<numVertices; ++iv)
      hostComputedDofs(computedDofsCount++) = cellBasis->getDofOrdinal(0, iv, 0);

    for(ordinal_type ie=0; ie<numEdges; ++ie)  {
      ordinal_type edgeCardinality = cellBasis->getDofCount(edgeDim,ie);
      for(ordinal_type i=0; i<edgeCardinality; ++i)
        hostComputedDofs(computedDofsCount++) = cellBasis->getDofOrdinal(edgeDim, ie, i);
    }

    for(ordinal_type iface=0; iface<numFaces; ++iface)  {
      ordinal_type faceCardinality = cellBasis->getDofCount(faceDim,iface);
      for(ordinal_type i=0; i<faceCardinality; ++i)
        hostComputedDofs(computedDofsCount++) = cellBasis->getDofOrdinal(faceDim, iface, i);
    }
    Kokkos::deep_copy(computedDofs,hostComputedDofs);
  }

  bool isHGradBasis = (cellBasis->getFunctionSpace() == FUNCTION_SPACE_HGRAD);
  bool isHCurlBasis = (cellBasis->getFunctionSpace() == FUNCTION_SPACE_HCURL);
  bool isHDivBasis = (cellBasis->getFunctionSpace() == FUNCTION_SPACE_HDIV);
  ordinal_type faceDofDim = isHCurlBasis ? 3 : 1;
  ScalarViewType edgeCoeff("edgeCoeff", fieldDim);


  const Kokkos::RangePolicy<ExecSpaceType> policy(0, numCells);
  ScalarViewType refBasisCoeffs("refBasisCoeffs", basisCoeffs.extent(0), basisCoeffs.extent(1));

  if(isHGradBasis) {

    auto targetEPointsRange = Kokkos::create_mirror_view_and_copy(MemSpaceType(), projStruct->getTargetPointsRange());
    using functorType = FunctorsProjectionTools::ComputeBasisCoeffsOnVertices_L2<ScalarViewType, decltype(tagToOrdinal), decltype(targetEPointsRange),
        decltype(targetAtTargetEPoints)>;
    Kokkos::parallel_for(policy, functorType(refBasisCoeffs, tagToOrdinal, targetEPointsRange,
        targetAtTargetEPoints, basisAtTargetEPoints, numVertices));
  }

  auto targetEPointsRange = projStruct->getTargetPointsRange();
  for(ordinal_type ie=0; ie<numEdges; ++ie)  {
    auto edgeVec = Kokkos::subview(refEdgesVec, ie, Kokkos::ALL());
    //auto edgeVecHost = Kokkos::create_mirror_view(edgeVec);

    if(isHCurlBasis) {
      CellTools<DeviceType>::getReferenceEdgeTangent(edgeVec, ie, cellTopo);
    } else if(isHDivBasis) {
      CellTools<DeviceType>::getReferenceSideNormal(edgeVec, ie, cellTopo);
    } else {
       deep_copy(edgeVec, 1.0);
    }

    ordinal_type edgeCardinality = cellBasis->getDofCount(edgeDim,ie);
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(edgeDim, ie));
    ordinal_type numTargetEPoints = range_size(targetEPointsRange(edgeDim, ie));

    ScalarViewType basisDofAtBasisEPoints("BasisDofAtBasisEPoints",1,edgeCardinality, numBasisEPoints);
    ScalarViewType tragetDofAtTargetEPoints("TargetDofAtTargetEPoints",numCells, numTargetEPoints);
    ScalarViewType weightedBasisAtBasisEPoints("weightedTanBasisAtBasisEPoints",numCells,edgeCardinality, numBasisEPoints);
    ScalarViewType weightedBasisAtTargetEPoints("weightedTanBasisAtTargetEPoints",numCells,edgeCardinality, numTargetEPoints);
    ScalarViewType negPartialProj("negPartialProj", numCells, numBasisEPoints);

    auto targetEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetEvalWeights(edgeDim,ie));
    auto basisEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisEvalWeights(edgeDim,ie));

    //Note: we are not considering the jacobian of the orientation map since it is simply a scalar term for the integrals and it does not affect the projection
    ordinal_type offsetBasis = basisEPointsRange(edgeDim, ie).first;
    ordinal_type offsetTarget = targetEPointsRange(edgeDim, ie).first;


    using functorTypeEdge = FunctorsProjectionTools::ComputeBasisCoeffsOnEdges_L2<ScalarViewType,  decltype(basisEWeights),
        decltype(computedDofs), decltype(tagToOrdinal), decltype(targetAtTargetEPoints)>;

    Kokkos::parallel_for(policy, functorTypeEdge(refBasisCoeffs, negPartialProj, basisDofAtBasisEPoints,
        basisAtBasisEPoints, basisEWeights, weightedBasisAtBasisEPoints, targetEWeights,
        basisAtTargetEPoints, weightedBasisAtTargetEPoints, computedDofs, tagToOrdinal,
        targetAtTargetEPoints,tragetDofAtTargetEPoints, refEdgesVec, fieldDim,
        edgeCardinality, offsetBasis, offsetTarget, numVertexDofs, edgeDim, ie));


    ScalarViewType edgeMassMat_("edgeMassMat_", 1, edgeCardinality, edgeCardinality),
        edgeRhsMat_("rhsMat_", numCells, edgeCardinality);

    FunctionSpaceTools<DeviceType >::integrate(edgeMassMat_, basisDofAtBasisEPoints, Kokkos::subview(weightedBasisAtBasisEPoints, std::make_pair(0,1), Kokkos::ALL(),Kokkos::ALL()) );
    FunctionSpaceTools<DeviceType >::integrate(edgeRhsMat_, tragetDofAtTargetEPoints, weightedBasisAtTargetEPoints);
    FunctionSpaceTools<DeviceType >::integrate(edgeRhsMat_, negPartialProj, weightedBasisAtBasisEPoints,true);


    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, DeviceType> WorkArrayViewType;
    ScalarViewType t_("t",numCells, edgeCardinality);
    WorkArrayViewType w_("w",numCells, edgeCardinality);

    auto edgeDof = Kokkos::subview(tagToOrdinal, edgeDim, ie, Kokkos::ALL());

    ElemSystem edgeSystem("edgeSystem", true);
    edgeSystem.solve(refBasisCoeffs, edgeMassMat_, edgeRhsMat_, t_, w_, edgeDof, edgeCardinality);
  }

  for(ordinal_type iface=0; iface<numFaces; ++iface) {
    ordinal_type faceCardinality = cellBasis->getDofCount(faceDim,iface);

    ordinal_type numTargetEPoints = range_size(targetEPointsRange(faceDim, iface));
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(faceDim, iface));

    ScalarViewType faceBasisDofAtBasisEPoints("faceBasisDofAtBasisEPoints",1,faceCardinality, numBasisEPoints,faceDofDim);
    ScalarViewType wBasisDofAtBasisEPoints("weightedBasisDofAtBasisEPoints",numCells,faceCardinality, numBasisEPoints,faceDofDim);

    ScalarViewType faceBasisAtTargetEPoints("faceBasisDofAtTargetEPoints",numCells,faceCardinality, numTargetEPoints,faceDofDim);
    ScalarViewType wBasisDofAtTargetEPoints("weightedBasisDofAtTargetEPoints",numCells,faceCardinality, numTargetEPoints,faceDofDim);

    ScalarViewType targetDofAtTargetEPoints("targetDofAtTargetEPoints",numCells, numTargetEPoints,faceDofDim);
    ScalarViewType negPartialProj("negComputedProjection", numCells,numBasisEPoints,faceDofDim);

    ordinal_type offsetBasis = basisEPointsRange(faceDim, iface).first;
    ordinal_type offsetTarget = targetEPointsRange(faceDim, iface).first;
    auto targetEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetEvalWeights(faceDim,iface));
    auto basisEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisEvalWeights(faceDim,iface));

    ScalarViewType refSideNormal("refSideNormal", dim);
    CellTools<DeviceType>::getReferenceSideNormal(refSideNormal, iface, cellTopo);

    using functorTypeFace = FunctorsProjectionTools::ComputeBasisCoeffsOnFaces_L2<ScalarViewType,  decltype(basisEWeights),
        decltype(computedDofs), decltype(tagToOrdinal), decltype(targetAtTargetEPoints)>;

    Kokkos::parallel_for(policy, functorTypeFace(refBasisCoeffs, negPartialProj,faceBasisDofAtBasisEPoints,
        basisAtBasisEPoints, basisEWeights, wBasisDofAtBasisEPoints, targetEWeights,
        basisAtTargetEPoints, wBasisDofAtTargetEPoints, computedDofs, tagToOrdinal,
        targetAtTargetEPoints,targetDofAtTargetEPoints,
        refSideNormal, fieldDim, faceCardinality, offsetBasis,
        offsetTarget, numVertexDofs+numEdgeDofs, numFaces, faceDim,
        dim, iface, isHCurlBasis, isHDivBasis));

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, DeviceType> WorkArrayViewType;
    ScalarViewType faceMassMat_("faceMassMat_", 1, faceCardinality, faceCardinality),
        faceRhsMat_("rhsMat_", numCells, faceCardinality);

    FunctionSpaceTools<DeviceType >::integrate(faceMassMat_, faceBasisDofAtBasisEPoints, Kokkos::subview(wBasisDofAtBasisEPoints, std::make_pair(0,1), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) );
    FunctionSpaceTools<DeviceType >::integrate(faceRhsMat_, targetDofAtTargetEPoints, wBasisDofAtTargetEPoints);
    FunctionSpaceTools<DeviceType >::integrate(faceRhsMat_, negPartialProj, wBasisDofAtBasisEPoints,true);

    ScalarViewType t_("t",numCells, faceCardinality);
    WorkArrayViewType w_("w",numCells,faceCardinality);

    auto faceDof = Kokkos::subview(tagToOrdinal, faceDim, iface, Kokkos::ALL());

    ElemSystem faceSystem("faceSystem", true);
    faceSystem.solve(refBasisCoeffs, faceMassMat_, faceRhsMat_, t_, w_, faceDof, faceCardinality);
  }

  ordinal_type numElemDofs = cellBasis->getDofCount(dim,0);


  if(numElemDofs>0) {

    auto cellDofs = Kokkos::subview(tagToOrdinal, dim, 0, Kokkos::ALL());

    range_type cellPointsRange = targetEPointsRange(dim, 0);

    ordinal_type numTargetEPoints = range_size(targetEPointsRange(dim,0));
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(dim,0));

    ScalarViewType internalBasisAtBasisEPoints("internalBasisAtBasisEPoints",1,numElemDofs, numBasisEPoints, fieldDim);
    ScalarViewType negPartialProj("negPartialProj", numCells, numBasisEPoints, fieldDim);

    auto targetEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetEvalWeights(dim,0));
    auto basisEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisEvalWeights(dim,0));
    ordinal_type offsetBasis = basisEPointsRange(dim, 0).first;
    ordinal_type offsetTarget = targetEPointsRange(dim, 0).first;

    ScalarViewType wBasisAtBasisEPoints("weightedBasisAtBasisEPoints",numCells,numElemDofs, numBasisEPoints,fieldDim);
    ScalarViewType wBasisDofAtTargetEPoints("weightedBasisAtTargetEPoints",numCells,numElemDofs, numTargetEPoints,fieldDim);

    using functorType = FunctorsProjectionTools::ComputeBasisCoeffsOnCells_L2<ScalarViewType,  decltype(basisEWeights), decltype(computedDofs), decltype(cellDofs)>;
    Kokkos::parallel_for(policy, functorType( refBasisCoeffs, negPartialProj,  internalBasisAtBasisEPoints,
        basisAtBasisEPoints, basisEWeights,  wBasisAtBasisEPoints, targetEWeights, basisAtTargetEPoints, wBasisDofAtTargetEPoints,
        computedDofs, cellDofs, fieldDim, numElemDofs, offsetBasis, offsetTarget, numVertexDofs+numEdgeDofs+numFaceDofs));

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, DeviceType> WorkArrayViewType;
    ScalarViewType cellMassMat_("cellMassMat_", 1, numElemDofs, numElemDofs),
        cellRhsMat_("rhsMat_", numCells, numElemDofs);

    FunctionSpaceTools<DeviceType >::integrate(cellMassMat_, internalBasisAtBasisEPoints, Kokkos::subview(wBasisAtBasisEPoints, std::make_pair(0,1), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) );
    if(fieldDim==1)
      FunctionSpaceTools<DeviceType >::integrate(cellRhsMat_, Kokkos::subview(targetAtTargetEPoints,Kokkos::ALL(),cellPointsRange,Kokkos::ALL()),
          Kokkos::subview(wBasisDofAtTargetEPoints,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),0));
    else
      FunctionSpaceTools<DeviceType >::integrate(cellRhsMat_, Kokkos::subview(targetAtTargetEPoints,Kokkos::ALL(),cellPointsRange,Kokkos::ALL()), wBasisDofAtTargetEPoints);
    FunctionSpaceTools<DeviceType >::integrate(cellRhsMat_, negPartialProj, wBasisAtBasisEPoints, true);

    ScalarViewType t_("t",numCells, numElemDofs);
    WorkArrayViewType w_("w",numCells,numElemDofs);
    ElemSystem cellSystem("cellSystem", true);
    cellSystem.solve(refBasisCoeffs, cellMassMat_, cellRhsMat_, t_, w_, cellDofs, numElemDofs);
  }

  OrientationTools<DeviceType>::modifyBasisByOrientationInverse(basisCoeffs, refBasisCoeffs, orts, cellBasis, true);
}


template<typename DeviceType>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<DeviceType>::getL2DGBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtTargetEPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<DeviceType, typename BasisType::scalarType> * projStruct){
  
  Kokkos::DynRankView<typename BasisType::scalarType,DeviceType>  refBasisCoeffs("refBasisCoeffs", basisCoeffs.extent(0), basisCoeffs.extent(1));
  getL2DGBasisCoeffs(refBasisCoeffs, targetAtTargetEPoints, cellBasis, projStruct);

  OrientationTools<DeviceType>::modifyBasisByOrientationInverse(basisCoeffs, refBasisCoeffs, orts, cellBasis, true);
}


template<typename DeviceType>
template<typename basisViewType, typename targetViewType, typename BasisType>
void
ProjectionTools<DeviceType>::getL2DGBasisCoeffs(basisViewType basisCoeffs,
    const targetViewType targetAtTargetEPoints,
    const BasisType* cellBasis,
    ProjectionStruct<DeviceType, typename BasisType::scalarType> * projStruct){

  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,DeviceType> ScalarViewType;
  const auto cellTopo = cellBasis->getBaseCellTopology();

  ordinal_type dim = cellTopo.getDimension();

  auto basisEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),
      projStruct->getEvalPoints(dim,0,EvalPointsType::BASIS));
  auto targetEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),
      projStruct->getEvalPoints(dim,0,EvalPointsType::TARGET));

  ordinal_type numTotalTargetEPoints(targetAtTargetEPoints.extent(1));
  ordinal_type basisCardinality = cellBasis->getCardinality();
  ordinal_type numCells = targetAtTargetEPoints.extent(0);
  const ordinal_type fieldDim = (targetAtTargetEPoints.rank()==2) ? 1 : targetAtTargetEPoints.extent(2);

  ordinal_type numTotalBasisEPoints = projStruct->getNumBasisEvalPoints();
  ScalarViewType basisAtBasisEPoints("basisAtBasisEPoints",1,basisCardinality, numTotalBasisEPoints, fieldDim);
  ScalarViewType basisAtTargetEPoints("basisAtTargetEPoints",basisCardinality, numTotalTargetEPoints, fieldDim);
  {
    if(fieldDim == 1) {
      cellBasis->getValues(Kokkos::subview(basisAtTargetEPoints,Kokkos::ALL(),Kokkos::ALL(),0), targetEPoints);
      cellBasis->getValues(Kokkos::subview(basisAtBasisEPoints,0,Kokkos::ALL(),Kokkos::ALL(),0), basisEPoints);
    }
    else {
      cellBasis->getValues(basisAtTargetEPoints, targetEPoints);
      cellBasis->getValues(Kokkos::subview(basisAtBasisEPoints,0,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), basisEPoints);
    }
  }

  const Kokkos::RangePolicy<ExecSpaceType> policy(0, numCells);
  ordinal_type numElemDofs = cellBasis->getCardinality();

  auto targetEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetEvalWeights(dim,0));
  auto basisEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisEvalWeights(dim,0));

  ScalarViewType wBasisAtBasisEPoints("weightedBasisAtBasisEPoints",1,numElemDofs, numTotalBasisEPoints,fieldDim);
  ScalarViewType wBasisDofAtTargetEPoints("weightedBasisAtTargetEPoints",numCells,numElemDofs, numTotalTargetEPoints,fieldDim);
  Kokkos::DynRankView<ordinal_type, DeviceType> cellDofs("cellDoFs", numElemDofs);

  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElemDofs),
  KOKKOS_LAMBDA (const int &j) {
    for(ordinal_type d=0; d <fieldDim; ++d) {
      for(ordinal_type iq=0; iq < numTotalBasisEPoints; ++iq)
        wBasisAtBasisEPoints(0,j,iq,d) = basisAtBasisEPoints(0,j,iq,d) * basisEWeights(iq);
      for(ordinal_type iq=0; iq <numTotalTargetEPoints; ++iq) {
        wBasisDofAtTargetEPoints(0,j,iq,d) = basisAtTargetEPoints(j,iq,d)* targetEWeights(iq);
      }
    }
    cellDofs(j) = j;
  });
  RealSpaceTools<DeviceType>::clone(wBasisDofAtTargetEPoints, Kokkos::subview(wBasisDofAtTargetEPoints,0,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()));

  typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, DeviceType> WorkArrayViewType;
  ScalarViewType cellMassMat_("cellMassMat_", 1, numElemDofs, numElemDofs),
      cellRhsMat_("rhsMat_", numCells, numElemDofs);

  FunctionSpaceTools<DeviceType >::integrate(cellMassMat_, basisAtBasisEPoints, wBasisAtBasisEPoints);
  if(fieldDim==1)
    FunctionSpaceTools<DeviceType >::integrate(cellRhsMat_, targetAtTargetEPoints,
        Kokkos::subview(wBasisDofAtTargetEPoints,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),0));
  else
    FunctionSpaceTools<DeviceType >::integrate(cellRhsMat_, targetAtTargetEPoints, wBasisDofAtTargetEPoints);

  ScalarViewType t_("t",1, numElemDofs);
  WorkArrayViewType w_("w",numCells,numElemDofs);
  ElemSystem cellSystem("cellSystem", true);
  cellSystem.solve(basisCoeffs, cellMassMat_, cellRhsMat_, t_, w_, cellDofs, numElemDofs);
}

} //Intrepid2 namespace

#endif

