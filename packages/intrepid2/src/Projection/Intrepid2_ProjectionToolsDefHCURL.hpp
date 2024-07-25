// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ProjectionToolsDefHCURL.hpp
    \brief  Header file for the Intrepid2::ProjectionTools
            containing definitions for HCURL projections.
    \author Created by Mauro Perego
 */

#ifndef __INTREPID2_PROJECTIONTOOLSDEFHCURL_HPP__
#define __INTREPID2_PROJECTIONTOOLSDEFHCURL_HPP__

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"


namespace Intrepid2 {

namespace FunctorsProjectionTools {


template<typename ViewType1, typename ViewType2, typename ViewType3>
struct ComputeWBasisEdge_HCurl {
  ViewType1 basisTanAtBasisEPoints_;
  ViewType1 weightedTanBasisAtBasisEPoints_;
  const ViewType1 basisAtBasisEPoints_;
  const ViewType2 basisEWeights_;
  const ViewType1 refEdgeTangent_;
  const ViewType3 tagToOrdinal_;
  ordinal_type edgeDim_;
  ordinal_type iedge_;
  ordinal_type offsetBasis_;

  ComputeWBasisEdge_HCurl(ViewType1 basisTanAtBasisEPoints, ViewType1 weightedTanBasisAtBasisEPoints,
      ViewType1 basisAtBasisEPoints, ViewType2 basisEWeights,  ViewType1 refEdgeTangent, ViewType3 tagToOrdinal,  
      ordinal_type edgeDim, ordinal_type iedge, ordinal_type offsetBasis) :
        basisTanAtBasisEPoints_(basisTanAtBasisEPoints), weightedTanBasisAtBasisEPoints_(weightedTanBasisAtBasisEPoints), basisAtBasisEPoints_(basisAtBasisEPoints),
        basisEWeights_(basisEWeights), refEdgeTangent_(refEdgeTangent), tagToOrdinal_(tagToOrdinal),
        edgeDim_(edgeDim), iedge_(iedge), offsetBasis_(offsetBasis) {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type j, const ordinal_type iq) const {
    ordinal_type jdof = tagToOrdinal_(edgeDim_, iedge_, j);
    for(ordinal_type d=0; d <ordinal_type(refEdgeTangent_.extent(0)); ++d)
      basisTanAtBasisEPoints_(0,j,iq) += refEdgeTangent_(d)*basisAtBasisEPoints_(jdof,offsetBasis_+iq,d);
    weightedTanBasisAtBasisEPoints_(0,j,iq) = basisTanAtBasisEPoints_(0,j,iq)*basisEWeights_(iq);
  }
};

template<typename ViewType1, typename ViewType2, typename ViewType3, typename ViewType4>
struct ComputeBasisCoeffsOnEdges_HCurl {
  const ViewType2 targetEWeights_;
  const ViewType1 basisAtTargetEPoints_;
  const ViewType1 wTanBasisAtTargetEPoints_;
  const ViewType3 tagToOrdinal_;
  const ViewType4 targetAtTargetEPoints_;
  const ViewType1 targetTanAtTargetEPoints_;
  const ViewType1 refEdgeTangent_;
  ordinal_type edgeCardinality_;
  ordinal_type offsetTarget_;
  ordinal_type edgeDim_;
  ordinal_type dim_;
  ordinal_type iedge_;

  ComputeBasisCoeffsOnEdges_HCurl(
      const ViewType2 targetEWeights,
      const ViewType1 basisAtTargetEPoints, const ViewType1 wTanBasisAtTargetEPoints, const ViewType3 tagToOrdinal,
      const ViewType4 targetAtTargetEPoints, const ViewType1 targetTanAtTargetEPoints,
      const ViewType1 refEdgeTangent, ordinal_type edgeCardinality,
      ordinal_type offsetTarget, ordinal_type edgeDim,
      ordinal_type dim, ordinal_type iedge) :
        targetEWeights_(targetEWeights),
        basisAtTargetEPoints_(basisAtTargetEPoints), wTanBasisAtTargetEPoints_(wTanBasisAtTargetEPoints),
        tagToOrdinal_(tagToOrdinal), targetAtTargetEPoints_(targetAtTargetEPoints),
        targetTanAtTargetEPoints_(targetTanAtTargetEPoints),
        refEdgeTangent_(refEdgeTangent), edgeCardinality_(edgeCardinality),
        offsetTarget_(offsetTarget), edgeDim_(edgeDim), dim_(dim), iedge_(iedge)
  {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    ordinal_type numTargetEPoints = targetEWeights_.extent(0);
    typename ViewType1::value_type  tmp = 0; 
    for(ordinal_type j=0; j <edgeCardinality_; ++j) {
      ordinal_type jdof = tagToOrdinal_(edgeDim_, iedge_, j);
      for(ordinal_type iq=0; iq <numTargetEPoints; ++iq) {
        tmp = 0;
        for(ordinal_type d=0; d <dim_; ++d)
          tmp += refEdgeTangent_(d)*basisAtTargetEPoints_(jdof,offsetTarget_+iq,d);
        wTanBasisAtTargetEPoints_(ic,j,iq) = tmp*targetEWeights_(iq);
      }
    }
    for(ordinal_type iq=0; iq <numTargetEPoints; ++iq)
      for(ordinal_type d=0; d <dim_; ++d)
        targetTanAtTargetEPoints_(ic,iq) += refEdgeTangent_(d)*targetAtTargetEPoints_(ic,offsetTarget_+iq,d);
  }
};


template<typename ViewType1, typename ViewType2, typename ViewType3, typename ViewType4,
typename ViewType5>
struct ComputeBasisCoeffsOnFaces_HCurl {
  const ViewType1 basisCoeffs_;
  const ViewType1 negPartialProjTan_;
  const ViewType1 negPartialProjCurlNormal_;
  const ViewType1 hgradBasisGradAtBasisEPoints_;
  const ViewType1 wHgradBasisGradAtBasisEPoints_;
  const ViewType1 basisCurlAtBasisCurlEPoints_;
  const ViewType1 basisCurlNormalAtBasisCurlEPoints_;
  const ViewType1 basisAtBasisEPoints_;
  const ViewType1 normalTargetCurlAtTargetEPoints_;
  const ViewType1 basisTanAtBasisEPoints_;
  const ViewType1 hgradBasisGradAtTargetEPoints_;
  const ViewType1 wHgradBasisGradAtTargetEPoints_;
  const ViewType1 wNormalBasisCurlAtBasisCurlEPoints_;
  const ViewType1 basisCurlAtTargetCurlEPoints_;
  const ViewType1 wNormalBasisCurlBasisAtTargetCurlEPoints_;
  const ViewType2 targetAtTargetEPoints_;
  const ViewType1 targetTanAtTargetEPoints_;
  const ViewType2 targetCurlAtTargetCurlEPoints_;
  const ViewType3 basisEWeights_;
  const ViewType3 targetEWeights_;
  const ViewType3 basisCurlEWeights_;
  const ViewType3 targetCurlEWeights_;
  const ViewType4 tagToOrdinal_;
  const ViewType4 hGradTagToOrdinal_;
  const ViewType5 computedDofs_;
  const ViewType1 refFaceNormal_;
  ordinal_type offsetBasis_;
  ordinal_type offsetBasisCurl_;
  ordinal_type offsetTarget_;
  ordinal_type offsetTargetCurl_;
  ordinal_type iface_;
  ordinal_type hgradCardinality_;
  ordinal_type numFaces_;
  ordinal_type numFaceDofs_;
  ordinal_type numEdgeDofs_;
  ordinal_type faceDim_;
  ordinal_type dim_;




  ComputeBasisCoeffsOnFaces_HCurl(const ViewType1 basisCoeffs,
      const ViewType1 negPartialProjTan, const ViewType1 negPartialProjCurlNormal,
      const ViewType1 hgradBasisGradAtBasisEPoints, const ViewType1 wHgradBasisGradAtBasisEPoints,
      const ViewType1 basisCurlAtBasisCurlEPoints, const ViewType1 basisCurlNormalAtBasisCurlEPoints,
      const ViewType1 basisAtBasisEPoints,
      const ViewType1 normalTargetCurlAtTargetEPoints,
      const ViewType1 basisTanAtBasisEPoints,
      const ViewType1 hgradBasisGradAtTargetEPoints, const ViewType1 wHgradBasisGradAtTargetEPoints,
      const ViewType1 wNormalBasisCurlAtBasisCurlEPoints, const ViewType1 basisCurlAtTargetCurlEPoints,
      const ViewType1 wNormalBasisCurlBasisAtTargetCurlEPoints, const ViewType2 targetAtTargetEPoints,
      const ViewType1 targetTanAtTargetEPoints, const ViewType2 targetCurlAtTargetCurlEPoints,
      const ViewType3 basisEWeights, const ViewType3 targetEWeights,
      const ViewType3 basisCurlEWeights, const ViewType3 targetCurlEWeights, const ViewType4 tagToOrdinal,
      const ViewType4 hGradTagToOrdinal, const ViewType5 computedDofs, 
      const ViewType1 refFaceNormal , ordinal_type offsetBasis,
      ordinal_type offsetBasisCurl, ordinal_type offsetTarget,
      ordinal_type offsetTargetCurl, ordinal_type iface,
      ordinal_type hgradCardinality, ordinal_type numFaces,
      ordinal_type numFaceDofs, ordinal_type numEdgeDofs,
      ordinal_type faceDim, ordinal_type dim):
        basisCoeffs_(basisCoeffs),
        negPartialProjTan_(negPartialProjTan),  negPartialProjCurlNormal_(negPartialProjCurlNormal),
        hgradBasisGradAtBasisEPoints_(hgradBasisGradAtBasisEPoints),  wHgradBasisGradAtBasisEPoints_(wHgradBasisGradAtBasisEPoints),
        basisCurlAtBasisCurlEPoints_(basisCurlAtBasisCurlEPoints),  basisCurlNormalAtBasisCurlEPoints_(basisCurlNormalAtBasisCurlEPoints),
        basisAtBasisEPoints_(basisAtBasisEPoints),
        normalTargetCurlAtTargetEPoints_(normalTargetCurlAtTargetEPoints),  basisTanAtBasisEPoints_(basisTanAtBasisEPoints),
        hgradBasisGradAtTargetEPoints_(hgradBasisGradAtTargetEPoints),  wHgradBasisGradAtTargetEPoints_(wHgradBasisGradAtTargetEPoints),
        wNormalBasisCurlAtBasisCurlEPoints_(wNormalBasisCurlAtBasisCurlEPoints),  basisCurlAtTargetCurlEPoints_(basisCurlAtTargetCurlEPoints),
        wNormalBasisCurlBasisAtTargetCurlEPoints_(wNormalBasisCurlBasisAtTargetCurlEPoints),  targetAtTargetEPoints_(targetAtTargetEPoints),
        targetTanAtTargetEPoints_(targetTanAtTargetEPoints),  targetCurlAtTargetCurlEPoints_(targetCurlAtTargetCurlEPoints),
        basisEWeights_(basisEWeights),  targetEWeights_(targetEWeights),
        basisCurlEWeights_(basisCurlEWeights), targetCurlEWeights_(targetCurlEWeights),  tagToOrdinal_(tagToOrdinal),
        hGradTagToOrdinal_(hGradTagToOrdinal), computedDofs_(computedDofs),
        refFaceNormal_(refFaceNormal), offsetBasis_(offsetBasis),
        offsetBasisCurl_(offsetBasisCurl), offsetTarget_(offsetTarget),
        offsetTargetCurl_(offsetTargetCurl), iface_(iface),
        hgradCardinality_(hgradCardinality), numFaces_(numFaces),
        numFaceDofs_(numFaceDofs), numEdgeDofs_(numEdgeDofs),
        faceDim_(faceDim), dim_(dim){}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    typename ViewType1::value_type n[3] = {refFaceNormal_(0), refFaceNormal_(1), refFaceNormal_(2)};
    typename ViewType1::value_type tmp=0;

    ordinal_type numBasisEPoints = basisEWeights_.extent(0);
    ordinal_type numTargetEPoints = targetEWeights_.extent(0);
    for(ordinal_type j=0; j <hgradCardinality_; ++j) {
      ordinal_type face_dof = hGradTagToOrdinal_(faceDim_, iface_, j);
      for(ordinal_type iq=0; iq <numBasisEPoints; ++iq) {
        for(ordinal_type d=0; d <dim_; ++d) {
          ordinal_type dp1 = (d+1) % dim_;
          ordinal_type dp2 = (d+2) % dim_;
          // basis \times n
          wHgradBasisGradAtBasisEPoints_(ic,j,iq,d) = (hgradBasisGradAtBasisEPoints_(face_dof,iq,dp1)*n[dp2] - hgradBasisGradAtBasisEPoints_(face_dof,iq,dp2)*n[dp1]) * basisEWeights_(iq);
        }
      }
      
      for(ordinal_type iq=0; iq <numTargetEPoints; ++iq) {
        for(ordinal_type d=0; d <dim_; ++d) {
          ordinal_type dp1 = (d+1) % dim_;
          ordinal_type dp2 = (d+2) % dim_;
          wHgradBasisGradAtTargetEPoints_(ic,j,iq,d) = (hgradBasisGradAtTargetEPoints_(face_dof,iq,dp1)*n[dp2] - hgradBasisGradAtTargetEPoints_(face_dof,iq,dp2)*n[dp1]) * targetEWeights_(iq);
        }
      }
    }

    ordinal_type numBasisCurlEPoints = basisCurlEWeights_.extent(0);
    for(ordinal_type j=0; j <numFaceDofs_; ++j) {
      ordinal_type jdof = tagToOrdinal_(faceDim_, iface_, j);
      for(ordinal_type iq=0; iq <numBasisEPoints; ++iq) {
        for(ordinal_type d=0; d <dim_; ++d) {
          ordinal_type dp1 = (d+1) % dim_;
          ordinal_type dp2 = (d+2) % dim_;
          // basis \times n
          basisTanAtBasisEPoints_(0,j,iq,d) = basisAtBasisEPoints_(jdof,offsetBasis_+iq,dp1)*n[dp2] - basisAtBasisEPoints_(jdof,offsetBasis_+iq,dp2)*n[dp1];
        }
      }
      //Note: we are not considering the jacobian of the orientation map for normals since it is simply a scalar term for the integrals and it does not affect the projection
      for(ordinal_type iq=0; iq <numBasisCurlEPoints; ++iq) {
        tmp=0;
        for(ordinal_type d=0; d <dim_; ++d)
          tmp += n[d]*basisCurlAtBasisCurlEPoints_(jdof,offsetBasisCurl_+iq,d);
        basisCurlNormalAtBasisCurlEPoints_(0,j,iq) = tmp;
        wNormalBasisCurlAtBasisCurlEPoints_(ic,j,iq) = tmp * basisCurlEWeights_(iq);
      }

      ordinal_type numTargetCurlEPoints = targetCurlEWeights_.extent(0);
      for(ordinal_type iq=0; iq <numTargetCurlEPoints; ++iq) {
        tmp=0;
        // target \cdot n  
        for(ordinal_type d=0; d <dim_; ++d)
          tmp += n[d]*basisCurlAtTargetCurlEPoints_(jdof,offsetTargetCurl_+iq,d);
        wNormalBasisCurlBasisAtTargetCurlEPoints_(ic,j,iq) = tmp*targetCurlEWeights_(iq);
      }
    }

    for(ordinal_type j=0; j <numEdgeDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <numBasisEPoints; ++iq)
        for(ordinal_type d=0; d <dim_; ++d) {
          ordinal_type dp1 = (d+1) % dim_;
          ordinal_type dp2 = (d+2) % dim_;
          negPartialProjCurlNormal_(ic,iq) -=  n[d]*basisCoeffs_(ic,jdof)*basisCurlAtBasisCurlEPoints_(jdof,offsetBasisCurl_+iq,d);
          // basis \times n
          negPartialProjTan_(ic,iq,d) -=  (basisAtBasisEPoints_(jdof,offsetBasis_+iq,dp1)*n[dp2] - basisAtBasisEPoints_(jdof,offsetBasis_+iq,dp2)*n[dp1])*basisCoeffs_(ic,jdof);
        }
    }

    ordinal_type numTargetCurlEPoints = targetCurlEWeights_.extent(0);
    for(ordinal_type iq=0; iq <numTargetEPoints; ++iq)
      for(ordinal_type d=0; d <dim_; ++d) {
        ordinal_type dp1 = (d+1) % dim_;
        ordinal_type dp2 = (d+2) % dim_;
        // target \times n
        targetTanAtTargetEPoints_(ic,iq,d) = (targetAtTargetEPoints_(ic,offsetTarget_+iq,dp1)*n[dp2] - targetAtTargetEPoints_(ic,offsetTarget_+iq,dp2)*n[dp1]);
      }

    for(ordinal_type iq=0; iq <numTargetCurlEPoints; ++iq)
      for(ordinal_type d=0; d <dim_; ++d) {
        // target \cdot n
        normalTargetCurlAtTargetEPoints_(ic,iq) += n[d]*targetCurlAtTargetCurlEPoints_(ic,offsetTargetCurl_+iq,d);
      }
  }
};


template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4>
struct ComputeBasisCoeffsOnCell_HCurl {
  const ViewType1 basisCoeffs_;
  const ViewType1 negPartialProj_;
  const ViewType1 negPartialProjCurl_;
  const ViewType1 cellBasisAtBasisEPoints_;
  const ViewType1 cellBasisCurlAtBasisCurlEPoints_;
  const ViewType1 basisAtBasisEPoints_;
  const ViewType1 hgradBasisGradAtBasisEPoints_;
  const ViewType1 basisCurlAtBasisCurlEPoints_;
  const ViewType1 hgradBasisGradAtTargetEPoints_;
  const ViewType1 basisCurlAtTargetCurlEPoints_;
  const ViewType2 basisEWeights_;
  const ViewType2 basisCurlEWeights_;
  const ViewType1 wHgradBasisGradAtBasisEPoints_;
  const ViewType1 wBasisCurlAtBasisCurlEPoints_;
  const ViewType2 targetEWeights_;
  const ViewType2 targetCurlEWeights_;
  const ViewType1 wHgradBasisGradAtTargetEPoints_;
  const ViewType1 wBasisCurlAtTargetCurlEPoints_;
  const ViewType3 computedDofs_;
  const ViewType4 tagToOrdinal_;
  const ViewType4 hGradTagToOrdinal_;
  ordinal_type numCellDofs_;
  ordinal_type hgradCardinality_;
  ordinal_type offsetBasis_;
  ordinal_type offsetBasisCurl_;
  ordinal_type offsetTargetCurl_;
  ordinal_type numEdgeFaceDofs_;
  ordinal_type dim_;
  ordinal_type derDim_;

  ComputeBasisCoeffsOnCell_HCurl(const ViewType1 basisCoeffs, ViewType1 negPartialProj,  ViewType1 negPartialProjCurl,
      const ViewType1 cellBasisAtBasisEPoints, const ViewType1 cellBasisCurlAtBasisCurlEPoints,
      const ViewType1 basisAtBasisEPoints, const ViewType1 hgradBasisGradAtBasisEPoints, const ViewType1 basisCurlAtBasisCurlEPoints,
      const ViewType1 hgradBasisGradAtTargetEPoints,   const ViewType1 basisCurlAtTargetCurlEPoints,
      const ViewType2 basisEWeights,  const ViewType2 basisCurlEWeights,
      const ViewType1 wHgradBasisGradAtBasisEPoints, const ViewType1 wBasisCurlAtBasisCurlEPoints,
      const ViewType2 targetEWeights, const ViewType2 targetCurlEWeights,
      const ViewType1 wHgradBasisGradAtTargetEPoints,
      const ViewType1 wBasisCurlAtTargetCurlEPoints, const ViewType3 computedDofs,
      const ViewType4 tagToOrdinal, const ViewType4 hGradTagToOrdinal,
      ordinal_type numCellDofs, ordinal_type hgradCardinality,
      ordinal_type offsetBasis, ordinal_type offsetBasisCurl, ordinal_type offsetTargetCurl,
      ordinal_type numEdgeFaceDofs, ordinal_type dim, ordinal_type derDim) :
        basisCoeffs_(basisCoeffs), negPartialProj_(negPartialProj), negPartialProjCurl_(negPartialProjCurl),
        cellBasisAtBasisEPoints_(cellBasisAtBasisEPoints), cellBasisCurlAtBasisCurlEPoints_(cellBasisCurlAtBasisCurlEPoints),
        basisAtBasisEPoints_(basisAtBasisEPoints), hgradBasisGradAtBasisEPoints_(hgradBasisGradAtBasisEPoints),
        basisCurlAtBasisCurlEPoints_(basisCurlAtBasisCurlEPoints),
        hgradBasisGradAtTargetEPoints_(hgradBasisGradAtTargetEPoints),
        basisCurlAtTargetCurlEPoints_(basisCurlAtTargetCurlEPoints),
        basisEWeights_(basisEWeights), basisCurlEWeights_(basisCurlEWeights),
        wHgradBasisGradAtBasisEPoints_(wHgradBasisGradAtBasisEPoints),
        wBasisCurlAtBasisCurlEPoints_(wBasisCurlAtBasisCurlEPoints),
        targetEWeights_(targetEWeights), targetCurlEWeights_(targetCurlEWeights),
        wHgradBasisGradAtTargetEPoints_(wHgradBasisGradAtTargetEPoints),
        wBasisCurlAtTargetCurlEPoints_(wBasisCurlAtTargetCurlEPoints),
        computedDofs_(computedDofs), tagToOrdinal_(tagToOrdinal), hGradTagToOrdinal_(hGradTagToOrdinal),
        numCellDofs_(numCellDofs), hgradCardinality_(hgradCardinality),
        offsetBasis_(offsetBasis), offsetBasisCurl_(offsetBasisCurl), offsetTargetCurl_(offsetTargetCurl),
        numEdgeFaceDofs_(numEdgeFaceDofs), dim_(dim), derDim_(derDim) {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    ordinal_type numBasisPoints = basisEWeights_.extent(0);
    ordinal_type numBasisCurlPoints = basisCurlEWeights_.extent(0);
    ordinal_type numTargetPoints = targetEWeights_.extent(0);
    ordinal_type numTargetCurlPoints = targetCurlEWeights_.extent(0);
    for(ordinal_type j=0; j <hgradCardinality_; ++j) {
      ordinal_type idof = hGradTagToOrdinal_(dim_, 0, j);
      for(ordinal_type d=0; d <dim_; ++d) {
        for(ordinal_type iq=0; iq <numBasisPoints; ++iq)
          wHgradBasisGradAtBasisEPoints_(ic,j,iq,d) = hgradBasisGradAtBasisEPoints_(idof,iq,d)*basisEWeights_(iq);
        for(ordinal_type iq=0; iq <numTargetPoints; ++iq)
          wHgradBasisGradAtTargetEPoints_(ic,j,iq,d) = hgradBasisGradAtTargetEPoints_(idof,iq,d)*targetEWeights_(iq);
      }
    }
    for(ordinal_type j=0; j <numCellDofs_; ++j) {
      ordinal_type idof = tagToOrdinal_(dim_, 0, j);
      for(ordinal_type d=0; d <dim_; ++d)
        for(ordinal_type iq=0; iq <numBasisPoints; ++iq)
          cellBasisAtBasisEPoints_(0,j,iq,d)=basisAtBasisEPoints_(idof,offsetBasis_+iq,d);

      for(ordinal_type d=0; d <derDim_; ++d) {
        for(ordinal_type iq=0; iq <numBasisCurlPoints; ++iq) {
          cellBasisCurlAtBasisCurlEPoints_(0,j,iq,d)=basisCurlAtBasisCurlEPoints_(idof,offsetBasisCurl_+iq,d);
          wBasisCurlAtBasisCurlEPoints_(ic,j,iq,d)=cellBasisCurlAtBasisCurlEPoints_(0,j,iq,d)*basisCurlEWeights_(iq);
        }
        for(ordinal_type iq=0; iq <numTargetCurlPoints; ++iq)
          wBasisCurlAtTargetCurlEPoints_(ic,j,iq,d) = basisCurlAtTargetCurlEPoints_(idof,offsetTargetCurl_+iq,d)*targetCurlEWeights_(iq);
      }
    }
    for(ordinal_type j=0; j < numEdgeFaceDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type d=0; d <derDim_; ++d)
        for(ordinal_type iq=0; iq <numBasisCurlPoints; ++iq)
          negPartialProjCurl_(ic,iq,d) -=  basisCoeffs_(ic,jdof)*basisCurlAtBasisCurlEPoints_(jdof,offsetBasisCurl_+iq,d);
      for(ordinal_type d=0; d <dim_; ++d)
        for(ordinal_type iq=0; iq <numBasisPoints; ++iq)
          negPartialProj_(ic,iq,d) -=  basisCoeffs_(ic,jdof)*basisAtBasisEPoints_(jdof,offsetBasis_+iq,d);
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
ProjectionTools<DeviceType>::getHCurlBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtTargetEPoints,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetCurlAtTargetCurlEPoints,
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
  const ordinal_type derDim = dim == 3 ? dim : 1;

  const Kokkos::RangePolicy<ExecSpaceType> policy(0, numCells);

  const std::string& name = cellBasis->getName();

  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;

  ScalarViewType refEdgeTangent("refEdgeTangent", dim);
  ScalarViewType refFaceNormal("refFaceNormal", dim);

  ordinal_type numEdgeDofs(0);
  for(ordinal_type ie=0; ie<numEdges; ++ie)
    numEdgeDofs += cellBasis->getDofCount(edgeDim,ie);

  ordinal_type numTotalFaceDofs(0);
  for(ordinal_type iface=0; iface<numFaces; ++iface)
    numTotalFaceDofs += cellBasis->getDofCount(faceDim,iface);

  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(MemSpaceType(), cellBasis->getAllDofOrdinal());

  Kokkos::View<ordinal_type*, DeviceType> computedDofs("computedDofs",numEdgeDofs+numTotalFaceDofs);

  auto targetEPointsRange  = projStruct->getTargetPointsRange();
  auto targetCurlEPointsRange  = projStruct->getTargetDerivPointsRange();

  auto basisEPointsRange = projStruct->getBasisPointsRange();
  auto basisCurlEPointsRange = projStruct->getBasisDerivPointsRange();

  ordinal_type numTotalBasisEPoints = projStruct->getNumBasisEvalPoints(), numTotalBasisCurlEPoints = projStruct->getNumBasisDerivEvalPoints();
  auto basisEPoints = projStruct->getAllEvalPoints(EvalPointsType::BASIS);
  auto basisCurlEPoints = projStruct->getAllDerivEvalPoints(EvalPointsType::BASIS);

  ordinal_type numTotalTargetEPoints = projStruct->getNumTargetEvalPoints(), numTotalTargetCurlEPoints = projStruct->getNumTargetDerivEvalPoints();
  auto targetEPoints = projStruct->getAllEvalPoints(EvalPointsType::TARGET);
  auto targetCurlEPoints = projStruct->getAllDerivEvalPoints(EvalPointsType::TARGET);

  ScalarViewType basisAtBasisEPoints("basisAtBasisEPoints",basisCardinality, numTotalBasisEPoints, dim);
  ScalarViewType basisAtTargetEPoints("basisAtTargetEPoints",basisCardinality, numTotalTargetEPoints, dim);
  cellBasis->getValues(basisAtTargetEPoints, targetEPoints);
  cellBasis->getValues(basisAtBasisEPoints, basisEPoints);

  ScalarViewType basisCurlAtBasisCurlEPoints;
  ScalarViewType basisCurlAtTargetCurlEPoints;
  if(numTotalBasisCurlEPoints>0) {
    ScalarViewType nonOrientedBasisCurlAtTargetCurlEPoints, nonOrientedBasisCurlAtBasisCurlEPoints;
    if (dim == 3) {
      basisCurlAtBasisCurlEPoints = ScalarViewType ("basisCurlAtBasisCurlEPoints",basisCardinality, numTotalBasisCurlEPoints, dim);
      basisCurlAtTargetCurlEPoints = ScalarViewType("basisCurlAtTargetCurlEPoints",basisCardinality, numTotalTargetCurlEPoints, dim);
    } else {
      basisCurlAtBasisCurlEPoints = ScalarViewType ("basisCurlAtBasisCurlEPoints",basisCardinality, numTotalBasisCurlEPoints);
      basisCurlAtTargetCurlEPoints = ScalarViewType("basisCurlAtTargetCurlEPoints",basisCardinality, numTotalTargetCurlEPoints);
    }

    cellBasis->getValues(basisCurlAtBasisCurlEPoints, basisCurlEPoints,OPERATOR_CURL);
    cellBasis->getValues(basisCurlAtTargetCurlEPoints, targetCurlEPoints,OPERATOR_CURL);
  }

  ScalarViewType refBasisCoeffs("refBasisCoeffs", basisCoeffs.extent(0), basisCoeffs.extent(1));

  ordinal_type computedDofsCount = 0;
  for(ordinal_type ie=0; ie<numEdges; ++ie)  {

    ordinal_type edgeCardinality = cellBasis->getDofCount(edgeDim,ie);
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(edgeDim, ie));
    ordinal_type numTargetEPoints = range_size(targetEPointsRange(edgeDim, ie));

    CellTools<DeviceType>::getReferenceEdgeTangent(refEdgeTangent, ie, cellTopo);

    ScalarViewType basisTanAtBasisEPoints("basisTanAtBasisEPoints",1,edgeCardinality, numBasisEPoints);
    ScalarViewType weightedTanBasisAtBasisEPoints("weightedTanBasisAtBasisEPoints",1,edgeCardinality, numBasisEPoints);
    ScalarViewType weightedTanBasisAtTargetEPoints("weightedTanBasisAtTargetEPoints",numCells,edgeCardinality, numTargetEPoints);
    ScalarViewType targetTanAtTargetEPoints("normalTargetAtTargetEPoints",numCells, numTargetEPoints);

    auto targetEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetEvalWeights(edgeDim,ie));
    auto basisEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisEvalWeights(edgeDim,ie));

    ordinal_type offsetBasis = basisEPointsRange(edgeDim, ie).first;
    ordinal_type offsetTarget = targetEPointsRange(edgeDim, ie).first;

    using functorTypeWBasisEdge = FunctorsProjectionTools::ComputeWBasisEdge_HCurl<ScalarViewType,  decltype(basisEWeights), decltype(tagToOrdinal)>;
    Kokkos::parallel_for(Kokkos::MDRangePolicy<ExecSpaceType, Kokkos::Rank<2> >({0,0}, {edgeCardinality,numBasisEPoints}), 
      functorTypeWBasisEdge(basisTanAtBasisEPoints,weightedTanBasisAtBasisEPoints,basisAtBasisEPoints,basisEWeights,refEdgeTangent,tagToOrdinal,edgeDim,ie,offsetBasis));

    using functorTypeEdge = FunctorsProjectionTools::ComputeBasisCoeffsOnEdges_HCurl<ScalarViewType,  decltype(targetEWeights), decltype(tagToOrdinal), decltype(targetAtTargetEPoints)>;
    Kokkos::parallel_for(policy, functorTypeEdge(targetEWeights,
        basisAtTargetEPoints, weightedTanBasisAtTargetEPoints, tagToOrdinal,
        targetAtTargetEPoints, targetTanAtTargetEPoints,
        refEdgeTangent, edgeCardinality,
        offsetTarget, edgeDim,
        dim, ie));

    ScalarViewType edgeMassMat_("edgeMassMat_", 1, edgeCardinality+1, edgeCardinality+1),
        edgeRhsMat_("rhsMat_", numCells, edgeCardinality+1);

    ScalarViewType eWeights_("eWeights_", 1, 1, basisEWeights.extent(0)), targetEWeights_("targetEWeights", numCells, 1, targetEWeights.extent(0));
    RealSpaceTools<DeviceType>::clone(eWeights_, basisEWeights);
    RealSpaceTools<DeviceType>::clone(targetEWeights_, targetEWeights);

    range_type range_H(0, edgeCardinality);
    range_type range_B(edgeCardinality, edgeCardinality+1);
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(edgeMassMat_,Kokkos::ALL(),range_H,range_H), basisTanAtBasisEPoints, weightedTanBasisAtBasisEPoints);
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(edgeMassMat_,Kokkos::ALL(),range_H,range_B), basisTanAtBasisEPoints, eWeights_);
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(edgeRhsMat_,Kokkos::ALL(),range_H), targetTanAtTargetEPoints, weightedTanBasisAtTargetEPoints);
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(edgeRhsMat_,Kokkos::ALL(),range_B), targetTanAtTargetEPoints, targetEWeights_);

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, DeviceType> WorkArrayViewType;
    ScalarViewType t_("t",numCells, edgeCardinality+1);
    WorkArrayViewType w_("w",numCells, edgeCardinality+1);

    auto edgeDofs = Kokkos::subview(tagToOrdinal, edgeDim, ie, range_type(0,edgeCardinality));
    ElemSystem edgeSystem("edgeSystem", true);
    edgeSystem.solve(refBasisCoeffs, edgeMassMat_, edgeRhsMat_, t_, w_, edgeDofs, edgeCardinality, 1);

    auto computedEdgeDofs = Kokkos::subview(computedDofs, range_type(computedDofsCount,computedDofsCount+edgeCardinality));
    deep_copy(computedEdgeDofs, edgeDofs);
    computedDofsCount += edgeCardinality;

  }

  typename RefSubcellParametrization<DeviceType>::ConstViewType  subcellParamFace;
  Basis<DeviceType,scalarType,scalarType> *hgradBasis = NULL;
  if(numFaces>0) {
    subcellParamFace = RefSubcellParametrization<DeviceType>::get(faceDim, cellBasis->getBaseCellTopology().getKey());
    if(cellTopo.getKey() == shards::getCellTopologyData<shards::Hexahedron<8> >()->key)
      hgradBasis = new Basis_HGRAD_HEX_Cn_FEM<DeviceType,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key)
      hgradBasis = new Basis_HGRAD_TET_Cn_FEM<DeviceType,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Wedge<6> >()->key) {
      hgradBasis = new typename DerivedNodalBasisFamily<DeviceType,scalarType,scalarType>::HGRAD_WEDGE(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    }
    else {
      std::stringstream ss;
      ss << ">>> ERROR (Intrepid2::ProjectionTools::getHCurlBasisCoeffs): "
          << "Method not implemented for basis " << name;
      INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
    }
  }
  for(ordinal_type iface=0; iface<numFaces; ++iface) {

    ordinal_type numTargetEPoints = range_size(targetEPointsRange(faceDim, iface));
    ordinal_type numTargetCurlEPoints = range_size(targetCurlEPointsRange(faceDim, iface));
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(faceDim, iface));
    ordinal_type numBasisCurlEPoints = range_size(basisCurlEPointsRange(faceDim, iface));

    ordinal_type numFaceDofs = cellBasis->getDofCount(faceDim,iface);

    ScalarViewType hgradBasisGradAtBasisEPoints("hgradBasisGradAtBasisEPoints", hgradBasis->getCardinality(), numBasisEPoints, dim);
    ScalarViewType hgradBasisGradAtTargetEPoints("hgradBasisGradAtTargetEPoints", hgradBasis->getCardinality(), numTargetEPoints, dim);

    ordinal_type hgradCardinality = hgradBasis->getDofCount(faceDim,iface);

    hgradBasis->getValues(hgradBasisGradAtBasisEPoints, Kokkos::subview(basisEPoints, basisEPointsRange(faceDim, iface), Kokkos::ALL()), OPERATOR_GRAD);
    hgradBasis->getValues(hgradBasisGradAtTargetEPoints, Kokkos::subview(targetEPoints, targetEPointsRange(faceDim, iface), Kokkos::ALL()), OPERATOR_GRAD);
    
    //no need to orient these basis as they act locally as test functions.

    auto hGradTagToOrdinal = Kokkos::create_mirror_view_and_copy(MemSpaceType(), hgradBasis->getAllDofOrdinal());

    ScalarViewType basisTanAtBasisEPoints("basisTanAtBasisEPoints",1,numFaceDofs, numBasisEPoints,dim);
    ScalarViewType basisCurlNormalAtBasisCurlEPoints("normaBasisCurlAtBasisEPoints",1,numFaceDofs, numBasisCurlEPoints);
    ScalarViewType wNormalBasisCurlAtBasisCurlEPoints("weightedNormalBasisCurlAtBasisEPoints",numCells,numFaceDofs, numBasisCurlEPoints);

    ScalarViewType targetTanAtTargetEPoints("targetTanAtTargetEPoints",numCells, numTargetEPoints, dim);
    ScalarViewType normalTargetCurlAtTargetEPoints("normalTargetCurlAtTargetEPoints",numCells, numTargetCurlEPoints);
    ScalarViewType wNormalBasisCurlBasisAtTargetCurlEPoints("weightedNormalBasisCurlAtTargetCurlEPoints",numCells,numFaceDofs, numTargetCurlEPoints);

    ScalarViewType wHgradBasisGradAtBasisEPoints("wHgradBasisGradAtBasisEPoints",numCells, hgradCardinality, numBasisEPoints, dim);
    ScalarViewType wHgradBasisGradAtTargetEPoints("wHgradBasisGradAtTargetEPoints",numCells, hgradCardinality, numTargetEPoints, dim);
  
    ScalarViewType negPartialProjCurlNormal("mNormalComputedProjection", numCells,numBasisEPoints);
    ScalarViewType negPartialProjTan("negPartialProjTan", numCells,numBasisEPoints,dim);

    ordinal_type offsetBasis = basisEPointsRange(faceDim, iface).first;
    ordinal_type offsetBasisCurl = basisCurlEPointsRange(faceDim, iface).first;
    ordinal_type offsetTarget = targetEPointsRange(faceDim, iface).first;
    ordinal_type offsetTargetCurl = targetCurlEPointsRange(faceDim, iface).first;

    auto basisEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisEvalWeights(faceDim,iface));
    auto targetEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetEvalWeights(faceDim,iface));
    auto targetCurlEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetDerivEvalWeights(faceDim,iface));
    auto basisCurlEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisDerivEvalWeights(faceDim,iface));

    CellTools<DeviceType>::getReferenceFaceNormal(refFaceNormal, iface, cellTopo);

    using functorTypeFaces = FunctorsProjectionTools::ComputeBasisCoeffsOnFaces_HCurl<ScalarViewType,  decltype(targetAtTargetEPoints), decltype(basisEWeights),
        decltype(tagToOrdinal), decltype(computedDofs)>;
    Kokkos::parallel_for(policy, functorTypeFaces(refBasisCoeffs,
        negPartialProjTan, negPartialProjCurlNormal,
        hgradBasisGradAtBasisEPoints, wHgradBasisGradAtBasisEPoints,
        basisCurlAtBasisCurlEPoints, basisCurlNormalAtBasisCurlEPoints,
        basisAtBasisEPoints,
        normalTargetCurlAtTargetEPoints, basisTanAtBasisEPoints,
        hgradBasisGradAtTargetEPoints, wHgradBasisGradAtTargetEPoints,
        wNormalBasisCurlAtBasisCurlEPoints, basisCurlAtTargetCurlEPoints,
        wNormalBasisCurlBasisAtTargetCurlEPoints, targetAtTargetEPoints,
        targetTanAtTargetEPoints, targetCurlAtTargetCurlEPoints,
        basisEWeights, targetEWeights,
        basisCurlEWeights, targetCurlEWeights, tagToOrdinal,
        hGradTagToOrdinal,
        computedDofs, refFaceNormal, offsetBasis,
        offsetBasisCurl, offsetTarget,
        offsetTargetCurl, iface,
        hgradCardinality, numFaces,
        numFaceDofs, numEdgeDofs,
        faceDim, dim));


    ScalarViewType faceMassMat_("faceMassMat_", 1, numFaceDofs+hgradCardinality, numFaceDofs+hgradCardinality),
        faceRhsMat_("rhsMat_", numCells, numFaceDofs+hgradCardinality);
    range_type range_H(0, numFaceDofs);
    range_type range_B(numFaceDofs, numFaceDofs+hgradCardinality);
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(faceMassMat_,Kokkos::ALL(),range_H,range_H), basisCurlNormalAtBasisCurlEPoints, Kokkos::subview(wNormalBasisCurlAtBasisCurlEPoints, std::make_pair(0,1), Kokkos::ALL(), Kokkos::ALL()) );
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(faceMassMat_,Kokkos::ALL(),range_H,range_B), basisTanAtBasisEPoints, Kokkos::subview(wHgradBasisGradAtBasisEPoints, std::make_pair(0,1), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) );

    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(faceRhsMat_,Kokkos::ALL(),range_H), normalTargetCurlAtTargetEPoints, wNormalBasisCurlBasisAtTargetCurlEPoints);
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(faceRhsMat_,Kokkos::ALL(),range_H), negPartialProjCurlNormal, wNormalBasisCurlAtBasisCurlEPoints,true);

    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(faceRhsMat_,Kokkos::ALL(),range_B), targetTanAtTargetEPoints, wHgradBasisGradAtTargetEPoints);
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(faceRhsMat_,Kokkos::ALL(),range_B), negPartialProjTan, wHgradBasisGradAtBasisEPoints,true);


    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, DeviceType> WorkArrayViewType;
    ScalarViewType t_("t",numCells, numFaceDofs+hgradCardinality);
    WorkArrayViewType w_("w",numCells, numFaceDofs+hgradCardinality);

    auto faceDofs = Kokkos::subview(tagToOrdinal, faceDim, iface, range_type(0,numFaceDofs));
    ElemSystem faceSystem( "faceSystem", true);
    faceSystem.solve(refBasisCoeffs, faceMassMat_, faceRhsMat_, t_, w_, faceDofs, numFaceDofs, hgradCardinality);

    auto computedFaceDofs = Kokkos::subview(computedDofs, range_type(computedDofsCount,computedDofsCount+numFaceDofs));
    deep_copy(computedFaceDofs, faceDofs);
    computedDofsCount += numFaceDofs;

  }
  delete hgradBasis;

  ordinal_type numCellDofs = cellBasis->getDofCount(dim,0);
  if(numCellDofs>0) {
    if(cellTopo.getKey() == shards::getCellTopologyData<shards::Hexahedron<8> >()->key)
      hgradBasis = new Basis_HGRAD_HEX_Cn_FEM<DeviceType,scalarType,scalarType>(cellBasis->getDegree());
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key)
      hgradBasis = new Basis_HGRAD_TET_Cn_FEM<DeviceType,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Wedge<6> >()->key)
      hgradBasis = new typename DerivedNodalBasisFamily<DeviceType,scalarType,scalarType>::HGRAD_WEDGE(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Triangle<3> >()->key)
      hgradBasis = new Basis_HGRAD_TRI_Cn_FEM<DeviceType,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Quadrilateral<4> >()->key)
      hgradBasis = new Basis_HGRAD_QUAD_Cn_FEM<DeviceType,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else  {
      std::stringstream ss;
      ss << ">>> ERROR (Intrepid2::ProjectionTools::getHCurlBasisCoeffs): "
          << "Method not implemented for basis " << name;
      INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
    }

    range_type cellPointsRange = targetEPointsRange(dim, 0);
    range_type cellCurlPointsRange = targetCurlEPointsRange(dim, 0);

    ordinal_type numTargetCurlEPoints = range_size(targetCurlEPointsRange(dim,0));
    ordinal_type numBasisCurlEPoints = range_size(basisCurlEPointsRange(dim,0));
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(dim,0));
    ordinal_type numTargetEPoints = range_size(targetEPointsRange(dim,0));

    ScalarViewType hgradBasisGradAtBasisEPoints("hgradBasisGradAtBasisEPoints",hgradBasis->getCardinality(), numBasisEPoints, dim);
    ScalarViewType hgradBasisGradAtTargetEPoints("hgradBasisGradAtTargetEPoints",hgradBasis->getCardinality(), numTargetEPoints, dim);

    ordinal_type hgradCardinality = hgradBasis->getDofCount(dim,0);
    ScalarViewType wHgradBasisGradAtTargetEPoints("wHgradBasisGradAtTargetEPoints",numCells, hgradCardinality, numTargetEPoints, dim);
    ScalarViewType wHgradBasisGradAtBasisEPoints("wHgradBasisGradAtBasisEPoints",numCells, hgradCardinality, numBasisEPoints, dim);

    hgradBasis->getValues(hgradBasisGradAtBasisEPoints,Kokkos::subview(basisEPoints, basisEPointsRange(dim, 0), Kokkos::ALL()), OPERATOR_GRAD);
    hgradBasis->getValues(hgradBasisGradAtTargetEPoints,Kokkos::subview(targetEPoints, targetEPointsRange(dim, 0), Kokkos::ALL()),OPERATOR_GRAD);

    ScalarViewType cellBasisAtBasisEPoints("basisCellAtEPoints",1,numCellDofs, numBasisEPoints, dim);
    ScalarViewType cellBasisCurlAtCurlEPoints("cellBasisCurlAtCurlEPoints",1,numCellDofs, numBasisCurlEPoints, derDim);
    ScalarViewType negPartialProjCurl("negPartialProjCurl", numCells, numBasisEPoints, derDim);
    ScalarViewType negPartialProj("negPartialProj", numCells, numBasisEPoints, dim);
    ScalarViewType wBasisCurlAtCurlEPoints("weightedBasisCurlAtBasisEPoints",numCells,numCellDofs, numBasisCurlEPoints,derDim);
    ScalarViewType wBasisCurlBasisAtTargetCurlEPoints("weightedBasisCurlAtTargetCurlEPoints",numCells,numCellDofs, numTargetCurlEPoints,derDim);

    auto targetEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetEvalWeights(dim,0));
    auto basisEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisEvalWeights(dim,0));
    auto targetCurlEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getTargetDerivEvalWeights(dim,0));
    auto basisCurlEWeights = Kokkos::create_mirror_view_and_copy(MemSpaceType(),projStruct->getBasisDerivEvalWeights(dim,0));
    ordinal_type offsetBasis = basisEPointsRange(dim, 0).first;
    ordinal_type offsetBasisCurl = basisCurlEPointsRange(dim, 0).first;
    ordinal_type offsetTargetCurl = targetCurlEPointsRange(dim, 0).first;


    auto hGradTagToOrdinal = Kokkos::create_mirror_view_and_copy(MemSpaceType(), hgradBasis->getAllDofOrdinal());

    using functorTypeCell = FunctorsProjectionTools::ComputeBasisCoeffsOnCell_HCurl<ScalarViewType,  decltype(basisEWeights),
        decltype(computedDofs), decltype(tagToOrdinal)>;
    Kokkos::parallel_for(policy, functorTypeCell(refBasisCoeffs, negPartialProj, negPartialProjCurl,
        cellBasisAtBasisEPoints, cellBasisCurlAtCurlEPoints,
        basisAtBasisEPoints, hgradBasisGradAtBasisEPoints, basisCurlAtBasisCurlEPoints,
        hgradBasisGradAtTargetEPoints,  basisCurlAtTargetCurlEPoints,
        basisEWeights, basisCurlEWeights, wHgradBasisGradAtBasisEPoints,
        wBasisCurlAtCurlEPoints, targetEWeights, targetCurlEWeights,
        wHgradBasisGradAtTargetEPoints,
        wBasisCurlBasisAtTargetCurlEPoints, computedDofs,
        tagToOrdinal, hGradTagToOrdinal,
        numCellDofs, hgradCardinality,
        offsetBasis, offsetBasisCurl,  offsetTargetCurl,
        numEdgeDofs+numTotalFaceDofs, dim, derDim));

    ScalarViewType cellMassMat_("cellMassMat_", 1, numCellDofs+hgradCardinality, numCellDofs+hgradCardinality),
        cellRhsMat_("rhsMat_", numCells, numCellDofs+hgradCardinality);

    range_type range_H(0, numCellDofs);
    range_type range_B(numCellDofs, numCellDofs+hgradCardinality);
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(cellMassMat_,Kokkos::ALL(),range_H,range_H), cellBasisCurlAtCurlEPoints, Kokkos::subview(wBasisCurlAtCurlEPoints, std::make_pair(0,1), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) );
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(cellMassMat_,Kokkos::ALL(),range_H,range_B), cellBasisAtBasisEPoints, Kokkos::subview(wHgradBasisGradAtBasisEPoints, std::make_pair(0,1), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) );
    if(dim==3)
      FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_H), Kokkos::subview(targetCurlAtTargetCurlEPoints,Kokkos::ALL(),cellCurlPointsRange,Kokkos::ALL()), wBasisCurlBasisAtTargetCurlEPoints);
    else
      FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_H), Kokkos::subview(targetCurlAtTargetCurlEPoints,Kokkos::ALL(),cellCurlPointsRange), Kokkos::subview(wBasisCurlBasisAtTargetCurlEPoints,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),0));
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_H), negPartialProjCurl, wBasisCurlAtCurlEPoints, true);

    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_B), Kokkos::subview(targetAtTargetEPoints,Kokkos::ALL(),cellPointsRange,Kokkos::ALL()), wHgradBasisGradAtTargetEPoints);
    FunctionSpaceTools<DeviceType >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_B), negPartialProj, wHgradBasisGradAtBasisEPoints, true);

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, DeviceType> WorkArrayViewType;
    ScalarViewType t_("t",numCells, numCellDofs+hgradCardinality);
    WorkArrayViewType w_("w",numCells, numCellDofs+hgradCardinality);

    auto cellDofs = Kokkos::subview(tagToOrdinal, dim, 0, Kokkos::ALL());
    ElemSystem cellSystem( "cellSystem", true);
    cellSystem.solve(refBasisCoeffs, cellMassMat_, cellRhsMat_, t_, w_, cellDofs, numCellDofs, hgradCardinality);

    delete hgradBasis;
  }

  OrientationTools<DeviceType>::modifyBasisByOrientationInverse(basisCoeffs, refBasisCoeffs, orts, cellBasis, true);
}

}   // Intrepid2 namespace

#endif

