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

/** \file   Intrepid2_ProjectionToolsDefHCURL.hpp
    \brief  Header file for the Intrepid2::Experimental::ProjectionTools
            containing definitions for HCURL projections.
    \author Created by Mauro Perego
 */

#ifndef __INTREPID2_PROJECTIONTOOLSDEFHCURL_HPP__
#define __INTREPID2_PROJECTIONTOOLSDEFHCURL_HPP__

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


namespace Intrepid2 {
namespace Experimental {

template<typename ViewType1, typename ViewType2, typename ViewType3, typename ViewType4>
struct ComputeBasisCoeffsOnEdges_HCurl {
  const ViewType1 basisTanAtBasisEPoints_;
  const ViewType1 basisAtBasisEPoints_;
  const ViewType2 basisEWeights_;
  const ViewType1 wTanBasisAtBasisEPoints_;
  const ViewType2 targetEWeights_;
  const ViewType1 basisAtTargetEPoints_;
  const ViewType1 wTanBasisAtTargetEPoints_;
  const ViewType3 tagToOrdinal_;
  const ViewType4 targetAtTargetEPoints_;
  const ViewType1 targetTanAtTargetEPoints_;
  const ViewType1 refEdgesTangent_;
  ordinal_type edgeCardinality_;
  ordinal_type offsetBasis_;
  ordinal_type offsetTarget_;
  ordinal_type edgeDim_;
  ordinal_type dim_;
  ordinal_type iedge_;

  ComputeBasisCoeffsOnEdges_HCurl(const ViewType1 basisTanAtBasisEPoints,
      const ViewType1 basisAtBasisEPoints, const ViewType2 basisEWeights,  const ViewType1 wTanBasisAtBasisEPoints,   const ViewType2 targetEWeights,
      const ViewType1 basisAtTargetEPoints, const ViewType1 wTanBasisAtTargetEPoints, const ViewType3 tagToOrdinal,
      const ViewType4 targetAtTargetEPoints, const ViewType1 targetTanAtTargetEPoints,
      const ViewType1 refEdgesTangent, ordinal_type edgeCardinality, ordinal_type offsetBasis,
      ordinal_type offsetTarget, ordinal_type edgeDim,
      ordinal_type dim, ordinal_type iedge) :
        basisTanAtBasisEPoints_(basisTanAtBasisEPoints),
        basisAtBasisEPoints_(basisAtBasisEPoints), basisEWeights_(basisEWeights), wTanBasisAtBasisEPoints_(wTanBasisAtBasisEPoints), targetEWeights_(targetEWeights),
        basisAtTargetEPoints_(basisAtTargetEPoints), wTanBasisAtTargetEPoints_(wTanBasisAtTargetEPoints),
        tagToOrdinal_(tagToOrdinal), targetAtTargetEPoints_(targetAtTargetEPoints),
        targetTanAtTargetEPoints_(targetTanAtTargetEPoints),
        refEdgesTangent_(refEdgesTangent), edgeCardinality_(edgeCardinality), offsetBasis_(offsetBasis),
        offsetTarget_(offsetTarget), edgeDim_(edgeDim), dim_(dim), iedge_(iedge)
  {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    ordinal_type numBasisEPoints = basisEWeights_.extent(0);
    ordinal_type numTargetEPoints = targetEWeights_.extent(0);
    for(ordinal_type j=0; j <edgeCardinality_; ++j) {
      ordinal_type jdof = tagToOrdinal_(edgeDim_, iedge_, j);
      for(ordinal_type iq=0; iq <numBasisEPoints; ++iq) {
        for(ordinal_type d=0; d <dim_; ++d)
          basisTanAtBasisEPoints_(ic,j,iq) += refEdgesTangent_(iedge_,d)*basisAtBasisEPoints_(ic,jdof,offsetBasis_+iq,d);
        wTanBasisAtBasisEPoints_(ic,j,iq) = basisTanAtBasisEPoints_(ic,j,iq)*basisEWeights_(iq);
      }

      for(ordinal_type iq=0; iq <numTargetEPoints; ++iq) {
        typename ViewType2::value_type  tmp = 0;
        for(ordinal_type d=0; d <dim_; ++d)
          tmp += refEdgesTangent_(iedge_,d)*basisAtTargetEPoints_(ic,jdof,offsetTarget_+iq,d);
        wTanBasisAtTargetEPoints_(ic,j,iq) = tmp*targetEWeights_(iq);
      }
    }
    for(ordinal_type iq=0; iq <numTargetEPoints; ++iq)
      for(ordinal_type d=0; d <dim_; ++d)
        targetTanAtTargetEPoints_(ic,iq) += refEdgesTangent_(iedge_,d)*targetAtTargetEPoints_(ic,offsetTarget_+iq,d);
  }
};


template<typename ViewType1, typename ViewType2, typename ViewType3, typename ViewType4,
typename ViewType5, typename ViewType6, typename ViewType7, typename ViewType8, typename ViewType9>
struct ComputeBasisCoeffsOnFaces_HCurl {
  const ViewType1 basisCoeffs_;
  const ViewType2 orts_;
  const ViewType3 negPartialProjTan_;
  const ViewType3 negPartialProjCurlNormal_;
  const ViewType3 hgradBasisGradAtBasisEPoints_;
  const ViewType3 wHgradBasisGradAtBasisEPoints_;
  const ViewType3 basisCurlAtBasisCurlEPoints_;
  const ViewType3 basisCurlNormalAtBasisCurlEPoints_;
  const ViewType3 basisAtBasisEPoints_;
  const ViewType3 normalTargetCurlAtTargetEPoints_;
  const ViewType3 basisTanAtBasisEPoints_;
  const ViewType3 hgradBasisGradAtTargetEPoints_;
  const ViewType3 wHgradBasisGradAtTargetEPoints_;
  const ViewType3 wNormalBasisCurlAtBasisCurlEPoints_;
  const ViewType3 basisCurlAtTargetCurlEPoints_;
  const ViewType3 wNormalBasisCurlBasisAtTargetCurlEPoints_;
  const ViewType4 targetAtTargetEPoints_;
  const ViewType3 targetTanAtTargetEPoints_;
  const ViewType4 targetCurlAtTargetCurlEPoints_;
  const ViewType5 basisEWeights_;
  const ViewType5 targetEWeights_;
  const ViewType5 basisCurlEWeights_;
  const ViewType5 targetCurlEWeights_;
  const ViewType6 tagToOrdinal_;
  const ViewType6 hGradTagToOrdinal_;
  const ViewType7 refTopologyKey_;
  const ViewType8 faceParametrization_;
  const ViewType9 computedDofs_;
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
      const ViewType2 orts, const ViewType3 negPartialProjTan, const ViewType3 negPartialProjCurlNormal,
      const ViewType3 hgradBasisGradAtBasisEPoints, const ViewType3 wHgradBasisGradAtBasisEPoints,
      const ViewType3 basisCurlAtBasisCurlEPoints, const ViewType3 basisCurlNormalAtBasisCurlEPoints,
      const ViewType3 basisAtBasisEPoints,
      const ViewType3 normalTargetCurlAtTargetEPoints,
      const ViewType3 basisTanAtBasisEPoints,
      const ViewType3 hgradBasisGradAtTargetEPoints, const ViewType3 wHgradBasisGradAtTargetEPoints,
      const ViewType3 wNormalBasisCurlAtBasisCurlEPoints, const ViewType3 basisCurlAtTargetCurlEPoints,
      const ViewType3 wNormalBasisCurlBasisAtTargetCurlEPoints, const ViewType4 targetAtTargetEPoints,
      const ViewType3 targetTanAtTargetEPoints, const ViewType4 targetCurlAtTargetCurlEPoints,
      const ViewType5 basisEWeights, const ViewType5 targetEWeights,
      const ViewType5 basisCurlEWeights, const ViewType5 targetCurlEWeights, const ViewType6 tagToOrdinal,
      const ViewType6 hGradTagToOrdinal, const ViewType7 refTopologyKey,
      const ViewType8 faceParametrization,
      const ViewType9 computedDofs, ordinal_type offsetBasis,
      ordinal_type offsetBasisCurl, ordinal_type offsetTarget,
      ordinal_type offsetTargetCurl, ordinal_type iface,
      ordinal_type hgradCardinality, ordinal_type numFaces,
      ordinal_type numFaceDofs, ordinal_type numEdgeDofs,
      ordinal_type faceDim, ordinal_type dim):
        basisCoeffs_(basisCoeffs),
        orts_(orts),  negPartialProjTan_(negPartialProjTan),  negPartialProjCurlNormal_(negPartialProjCurlNormal),
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
        hGradTagToOrdinal_(hGradTagToOrdinal),  refTopologyKey_(refTopologyKey),
        faceParametrization_(faceParametrization),
        computedDofs_(computedDofs), offsetBasis_(offsetBasis),
        offsetBasisCurl_(offsetBasisCurl), offsetTarget_(offsetTarget),
        offsetTargetCurl_(offsetTargetCurl), iface_(iface),
        hgradCardinality_(hgradCardinality), numFaces_(numFaces),
        numFaceDofs_(numFaceDofs), numEdgeDofs_(numEdgeDofs),
        faceDim_(faceDim), dim_(dim){}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    ordinal_type fOrt[6];
    orts_(ic).getFaceOrientation(fOrt, numFaces_);

    ordinal_type ort = fOrt[iface_];

    typename ViewType3::value_type data[3*3];
    auto tangentsAndNormal = ViewType3(data, dim_, dim_);

    Impl::OrientationTools::getRefSideTangentsAndNormal(tangentsAndNormal, faceParametrization_,refTopologyKey_(faceDim_,iface_), iface_, ort);

    ordinal_type numBasisEPoints = basisEWeights_.extent(0);
    ordinal_type numTargetEPoints = targetEWeights_.extent(0);
    for(ordinal_type j=0; j <hgradCardinality_; ++j) {
      ordinal_type face_dof = hGradTagToOrdinal_(faceDim_, 0, j);
      for(ordinal_type d=0; d <faceDim_; ++d) {
        for(ordinal_type iq=0; iq <numBasisEPoints; ++iq)
          wHgradBasisGradAtBasisEPoints_(ic, j, iq, d) = hgradBasisGradAtBasisEPoints_(face_dof,iq,d) * basisEWeights_(iq);
        for(ordinal_type iq=0; iq <numTargetEPoints; ++iq)
          wHgradBasisGradAtTargetEPoints_(ic,j,iq,d)= hgradBasisGradAtTargetEPoints_(face_dof,iq,d) * targetEWeights_(iq);
      }
    }

    //Note: we are not considering the jacobian of the orientation map for normals since it is simply a scalar term for the integrals and it does not affect the projection
    ordinal_type numBasisCurlEPoints = basisCurlEWeights_.extent(0);
    for(ordinal_type j=0; j <numFaceDofs_; ++j) {
      ordinal_type jdof = tagToOrdinal_(faceDim_, iface_, j);
      for(ordinal_type iq=0; iq <numBasisEPoints; ++iq) {
        for(ordinal_type d=0; d <dim_; ++d) {
          for(ordinal_type itan=0; itan <dim_-1; ++itan)
            basisTanAtBasisEPoints_(ic,j,iq,itan) += tangentsAndNormal(itan,d)*basisAtBasisEPoints_(ic,jdof,offsetBasis_+iq,d);
        }
      }
      for(ordinal_type iq=0; iq <numBasisCurlEPoints; ++iq) {
        for(ordinal_type d=0; d <dim_; ++d)
          basisCurlNormalAtBasisCurlEPoints_(ic,j,iq) += tangentsAndNormal(dim_-1,d)*basisCurlAtBasisCurlEPoints_(ic,jdof,offsetBasisCurl_+iq,d);
        wNormalBasisCurlAtBasisCurlEPoints_(ic,j,iq) = basisCurlNormalAtBasisCurlEPoints_(ic,j,iq) * basisCurlEWeights_(iq);
      }

      ordinal_type numTargetCurlEPoints = targetCurlEWeights_.extent(0);
      for(ordinal_type iq=0; iq <numTargetCurlEPoints; ++iq) {
        typename ViewType3::value_type tmp=0;
        for(ordinal_type d=0; d <dim_; ++d)
          tmp += tangentsAndNormal(dim_-1,d)*basisCurlAtTargetCurlEPoints_(ic,jdof,offsetTargetCurl_+iq,d);
        wNormalBasisCurlBasisAtTargetCurlEPoints_(ic,j,iq) = tmp*targetCurlEWeights_(iq);
      }
    }

    for(ordinal_type j=0; j <numEdgeDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <numBasisEPoints; ++iq)
        for(ordinal_type d=0; d <dim_; ++d) {
          negPartialProjCurlNormal_(ic,iq) -=  tangentsAndNormal(dim_-1,d)*basisCoeffs_(ic,jdof)*basisCurlAtBasisCurlEPoints_(ic,jdof,offsetBasisCurl_+iq,d);
          for(ordinal_type itan=0; itan <dim_-1; ++itan)
              negPartialProjTan_(ic,iq,itan) -=  tangentsAndNormal(itan,d)*basisCoeffs_(ic,jdof)*basisAtBasisEPoints_(ic,jdof,offsetBasis_+iq,d);
        }
    }

    ordinal_type numTargetCurlEPoints = targetCurlEWeights_.extent(0);
    for(ordinal_type iq=0; iq <numTargetEPoints; ++iq)
      for(ordinal_type d=0; d <dim_; ++d)
        for(ordinal_type itan=0; itan <dim_-1; ++itan)
            targetTanAtTargetEPoints_(ic,iq,itan) += tangentsAndNormal(itan,d)*targetAtTargetEPoints_(ic,offsetTarget_+iq,d);

    for(ordinal_type iq=0; iq <numTargetCurlEPoints; ++iq)
      for(ordinal_type d=0; d <dim_; ++d)
        normalTargetCurlAtTargetEPoints_(ic,iq) += tangentsAndNormal(dim_-1,d)*targetCurlAtTargetCurlEPoints_(ic,offsetTargetCurl_+iq,d);
  }
};


template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5>
struct ComputeBasisCoeffsOnCell_HCurl {
  const ViewType1 basisCoeffs_;
  const ViewType2 negPartialProj_;
  const ViewType2 negPartialProjCurl_;
  const ViewType2 cellBasisAtBasisEPoints_;
  const ViewType2 cellBasisCurlAtBasisCurlEPoints_;
  const ViewType2 basisAtBasisEPoints_;
  const ViewType2 hgradBasisGradAtBasisEPoints_;
  const ViewType2 basisCurlAtBasisCurlEPoints_;
  const ViewType2 hgradBasisGradAtTargetEPoints_;
  const ViewType2 basisCurlAtTargetCurlEPoints_;
  const ViewType3 basisEWeights_;
  const ViewType3 basisCurlEWeights_;
  const ViewType2 wHgradBasisGradAtBasisEPoints_;
  const ViewType2 wBasisCurlAtBasisCurlEPoints_;
  const ViewType3 targetEWeights_;
  const ViewType3 targetCurlEWeights_;
  const ViewType2 wHgradBasisGradAtTargetEPoints_;
  const ViewType2 wBasisCurlAtTargetCurlEPoints_;
  const ViewType4 computedDofs_;
  const ViewType5 tagToOrdinal_;
  const ViewType5 hGradTagToOrdinal_;
  ordinal_type numCellDofs_;
  ordinal_type hgradCardinality_;
  ordinal_type offsetBasis_;
  ordinal_type offsetBasisCurl_;
  ordinal_type offsetTargetCurl_;
  ordinal_type numEdgeFaceDofs_;
  ordinal_type dim_;
  ordinal_type derDim_;

  ComputeBasisCoeffsOnCell_HCurl(const ViewType1 basisCoeffs, ViewType2 negPartialProj,  ViewType2 negPartialProjCurl,
      const ViewType2 cellBasisAtBasisEPoints, const ViewType2 cellBasisCurlAtBasisCurlEPoints,
      const ViewType2 basisAtBasisEPoints, const ViewType2 hgradBasisGradAtBasisEPoints, const ViewType2 basisCurlAtBasisCurlEPoints,
      const ViewType2 hgradBasisGradAtTargetEPoints,   const ViewType2 basisCurlAtTargetCurlEPoints,
      const ViewType3 basisEWeights,  const ViewType3 basisCurlEWeights,
      const ViewType2 wHgradBasisGradAtBasisEPoints, const ViewType2 wBasisCurlAtBasisCurlEPoints,
      const ViewType3 targetEWeights, const ViewType3 targetCurlEWeights,
      const ViewType2 wHgradBasisGradAtTargetEPoints,
      const ViewType2 wBasisCurlAtTargetCurlEPoints, const ViewType4 computedDofs,
      const ViewType5 tagToOrdinal, const ViewType5 hGradTagToOrdinal,
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
          cellBasisAtBasisEPoints_(ic,j,iq,d)=basisAtBasisEPoints_(ic,idof,offsetBasis_+iq,d);

      for(ordinal_type d=0; d <derDim_; ++d) {
        for(ordinal_type iq=0; iq <numBasisCurlPoints; ++iq) {
          cellBasisCurlAtBasisCurlEPoints_(ic,j,iq,d)=basisCurlAtBasisCurlEPoints_(ic,idof,offsetBasisCurl_+iq,d);
          wBasisCurlAtBasisCurlEPoints_(ic,j,iq,d)=cellBasisCurlAtBasisCurlEPoints_(ic,j,iq,d)*basisCurlEWeights_(iq);
        }
        for(ordinal_type iq=0; iq <numTargetCurlPoints; ++iq)
          wBasisCurlAtTargetCurlEPoints_(ic,j,iq,d) = basisCurlAtTargetCurlEPoints_(ic,idof,offsetTargetCurl_+iq,d)*targetCurlEWeights_(iq);
      }
    }
    for(ordinal_type j=0; j < numEdgeFaceDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type d=0; d <derDim_; ++d)
        for(ordinal_type iq=0; iq <numBasisCurlPoints; ++iq)
          negPartialProjCurl_(ic,iq,d) -=  basisCoeffs_(ic,jdof)*basisCurlAtBasisCurlEPoints_(ic,jdof,offsetBasisCurl_+iq,d);
      for(ordinal_type d=0; d <dim_; ++d)
        for(ordinal_type iq=0; iq <numBasisPoints; ++iq)
          negPartialProj_(ic,iq,d) -=  basisCoeffs_(ic,jdof)*basisAtBasisEPoints_(ic,jdof,offsetBasis_+iq,d);
    }
  }
};


template<typename SpT>
template<typename BasisType,
typename ortValueType,       class ...ortProperties>
void
ProjectionTools<SpT>::getHCurlEvaluationPoints(typename BasisType::ScalarViewType targetEPoints,
    typename BasisType::ScalarViewType targetCurlEPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct,
    const EvalPointsType evalPointType) {
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type numCells = targetEPoints.extent(0);
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;

  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;

  CellTools<SpT>::setSubcellParametrization();
  typename CellTools<SpT>::subcellParamViewType subcellParamEdge,  subcellParamFace;
  if(numEdges>0)
    CellTools<SpT>::getSubcellParametrization(subcellParamEdge,  edgeDim, cellTopo);
  if(numFaces>0)
    CellTools<SpT>::getSubcellParametrization(subcellParamFace,  faceDim, cellTopo);

  auto refTopologyKey =  Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTopologyKey());

  auto evalPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getPointsRange(evalPointType));
  auto curlEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getDerivPointsRange(evalPointType));

  for(ordinal_type ie=0; ie<numEdges; ++ie) {

    auto edgePointsRange = evalPointsRange(edgeDim, ie);
    auto edgeEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getEvalPoints(edgeDim,ie,evalPointType));

    Kokkos::parallel_for
    ("Evaluate Points Edges ",
        Kokkos::RangePolicy<SpT, int> (0, numCells),
        KOKKOS_LAMBDA (const size_t ic) {

      ordinal_type eOrt[12];
      orts(ic).getEdgeOrientation(eOrt, numEdges);
      ordinal_type ort = eOrt[ie];

      Impl::OrientationTools::mapSubcellCoordsToRefCell(Kokkos::subview(targetEPoints,ic,edgePointsRange,Kokkos::ALL()),
          edgeEPoints, subcellParamEdge, refTopologyKey(edgeDim, ie), ie, ort);
    });
  }


  for(ordinal_type iface=0; iface<numFaces; ++iface) {
    auto faceCurlPointsRange = curlEPointsRange(faceDim, iface);
    auto faceCurlEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getDerivEvalPoints(faceDim,iface,evalPointType));

    auto facePointsRange = evalPointsRange(faceDim, iface);
    auto faceEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getEvalPoints(faceDim,iface,evalPointType));

    Kokkos::parallel_for
    ("Evaluate Points Faces ",
        Kokkos::RangePolicy<SpT, int> (0, numCells),
        KOKKOS_LAMBDA (const size_t ic) {

      ordinal_type fOrt[6];
      orts(ic).getFaceOrientation(fOrt, numFaces);
      ordinal_type ort = fOrt[iface];

      Impl::OrientationTools::mapSubcellCoordsToRefCell(Kokkos::subview(targetEPoints, ic, facePointsRange, Kokkos::ALL()),
          faceEPoints, subcellParamFace, refTopologyKey(faceDim, iface), iface, ort);

      Impl::OrientationTools::mapSubcellCoordsToRefCell(Kokkos::subview(targetCurlEPoints,  ic, faceCurlPointsRange, Kokkos::ALL()),
          faceCurlEPoints, subcellParamFace, refTopologyKey(faceDim, iface), iface, ort);
    });
  }


  if(cellBasis->getDofCount(dim,0)>0) {
    auto cellEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getEvalPoints(dim,0,evalPointType));
    RealSpaceTools<SpT>::clone(Kokkos::subview(targetEPoints, Kokkos::ALL(), evalPointsRange(dim, 0), Kokkos::ALL()), cellEPoints);

    auto cellCurlEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getDerivEvalPoints(dim,0,evalPointType));
    RealSpaceTools<SpT>::clone(Kokkos::subview(targetCurlEPoints, Kokkos::ALL(), curlEPointsRange(dim, 0), Kokkos::ALL()), cellCurlEPoints);
  }
}


template<typename SpT>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<SpT>::getHCurlBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtTargetEPoints,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetCurlAtTargetCurlEPoints,
    const typename BasisType::ScalarViewType targetEPoints,
    const typename BasisType::ScalarViewType targetCurlEPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct){

  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> ScalarViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type numTotalTargetEPoints(targetAtTargetEPoints.extent(1)),
      numTotalTargetCurlEPoints(targetCurlAtTargetCurlEPoints.extent(1));
  ordinal_type basisCardinality = cellBasis->getCardinality();
  ordinal_type numCells = targetAtTargetEPoints.extent(0);
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;
  const ordinal_type derDim = dim == 3 ? dim : 1;

  const Kokkos::RangePolicy<SpT> policy(0, numCells);

  const std::string& name = cellBasis->getName();

  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;

  ScalarViewType refEdgesTangent("refEdgesTangent", numEdges, dim);
  ScalarViewType refFacesTangents("refFaceTangents", numFaces, dim, 2);
  ScalarViewType refFacesNormal("refFaceNormal", numFaces, dim);

  ordinal_type numEdgeDofs(0);
  for(ordinal_type ie=0; ie<numEdges; ++ie)
    numEdgeDofs += cellBasis->getDofCount(edgeDim,ie);

  ordinal_type numFaceDofs(0);
  for(ordinal_type iface=0; iface<numFaces; ++iface)
    numFaceDofs += cellBasis->getDofCount(faceDim,iface);

  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(), cellBasis->getAllDofOrdinal());

  Kokkos::View<ordinal_type*, SpT> computedDofs("computedDofs",numEdgeDofs+numFaceDofs);

  auto targetEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetPointsRange());
  auto targetCurlEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetDerivPointsRange());

  auto basisEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisPointsRange());
  auto basisCurlEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisDerivPointsRange());

  auto refTopologyKey =  Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTopologyKey());

  ordinal_type numTotalBasisEPoints = projStruct->getNumBasisEvalPoints(), numTotalBasisCurlEPoints = projStruct->getNumBasisDerivEvalPoints();

  ScalarViewType basisEPoints("basisEPoints",numCells,numTotalBasisEPoints, dim);
  ScalarViewType basisCurlEPoints("basisCurlEPoints",numCells,numTotalBasisCurlEPoints, dim);
  getHCurlEvaluationPoints(basisEPoints, basisCurlEPoints, orts, cellBasis, projStruct, EvalPointsType::BASIS);

  ScalarViewType basisAtBasisEPoints("basisAtBasisEPoints",numCells,basisCardinality, numTotalBasisEPoints, dim);
  ScalarViewType basisAtTargetEPoints("basisAtTargetEPoints",numCells,basisCardinality, numTotalTargetEPoints, dim);
  {
    ScalarViewType nonOrientedBasisAtBasisEPoints("nonOrientedBasisAtEPoints",numCells,basisCardinality, numTotalBasisEPoints, dim);
    ScalarViewType nonOrientedBasisAtTargetEPoints("nonOrientedBasisAtTargetEPoints",numCells,basisCardinality, numTotalTargetEPoints, dim);
    for(ordinal_type ic=0; ic<numCells; ++ic) {
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtTargetEPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(targetEPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtBasisEPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(basisEPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
    }

    OrientationTools<SpT>::modifyBasisByOrientation(basisAtBasisEPoints, nonOrientedBasisAtBasisEPoints, orts, cellBasis);
    OrientationTools<SpT>::modifyBasisByOrientation(basisAtTargetEPoints, nonOrientedBasisAtTargetEPoints, orts, cellBasis);
  }

  ScalarViewType basisCurlAtBasisCurlEPoints;
  ScalarViewType basisCurlAtTargetCurlEPoints;
  if(numTotalBasisCurlEPoints>0) {
    ScalarViewType nonOrientedBasisCurlAtTargetCurlEPoints, nonOrientedBasisCurlAtBasisCurlEPoints;
    if (dim == 3) {
      basisCurlAtBasisCurlEPoints = ScalarViewType ("basisCurlAtBasisCurlEPoints",numCells,basisCardinality, numTotalBasisCurlEPoints, dim);
      nonOrientedBasisCurlAtBasisCurlEPoints = ScalarViewType ("nonOrientedBasisCurlAtBasisCurlEPoints",numCells,basisCardinality, numTotalBasisCurlEPoints, dim);
      basisCurlAtTargetCurlEPoints = ScalarViewType("basisCurlAtTargetCurlEPoints",numCells,basisCardinality, numTotalTargetCurlEPoints, dim);
      nonOrientedBasisCurlAtTargetCurlEPoints = ScalarViewType("nonOrientedBasisCurlAtTargetCurlEPoints",numCells,basisCardinality, numTotalTargetCurlEPoints, dim);
    } else {
      basisCurlAtBasisCurlEPoints = ScalarViewType ("basisCurlAtBasisCurlEPoints",numCells,basisCardinality, numTotalBasisCurlEPoints);
      nonOrientedBasisCurlAtBasisCurlEPoints = ScalarViewType ("nonOrientedBasisCurlAtBasisCurlEPoints",numCells,basisCardinality, numTotalBasisCurlEPoints);
      basisCurlAtTargetCurlEPoints = ScalarViewType("basisCurlAtTargetCurlEPoints",numCells,basisCardinality, numTotalTargetCurlEPoints);
      nonOrientedBasisCurlAtTargetCurlEPoints = ScalarViewType("nonOrientedBasisCurlAtTargetCurlEPoints",numCells,basisCardinality, numTotalTargetCurlEPoints);
    }
    for(ordinal_type ic=0; ic<numCells; ++ic) {
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisCurlAtBasisCurlEPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(basisCurlEPoints, ic, Kokkos::ALL(), Kokkos::ALL()),OPERATOR_CURL);
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisCurlAtTargetCurlEPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(targetCurlEPoints, ic, Kokkos::ALL(), Kokkos::ALL()),OPERATOR_CURL);
    }
    OrientationTools<SpT>::modifyBasisByOrientation(basisCurlAtBasisCurlEPoints, nonOrientedBasisCurlAtBasisCurlEPoints, orts, cellBasis);
    OrientationTools<SpT>::modifyBasisByOrientation(basisCurlAtTargetCurlEPoints, nonOrientedBasisCurlAtTargetCurlEPoints, orts, cellBasis);
  }

  ordinal_type computedDofsCount = 0;
  for(ordinal_type ie=0; ie<numEdges; ++ie)  {

    ordinal_type edgeCardinality = cellBasis->getDofCount(edgeDim,ie);
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(edgeDim, ie));
    ordinal_type numTargetEPoints = range_size(targetEPointsRange(edgeDim, ie));

    {
      auto refEdgeTan = Kokkos::subview(refEdgesTangent, ie, Kokkos::ALL());
      auto refEdgeTanHost = Kokkos::create_mirror_view(refEdgeTan);
      CellTools<SpT>::getReferenceEdgeTangent(refEdgeTanHost, ie, cellTopo);
      Kokkos::deep_copy(refEdgeTan,refEdgeTanHost);
    }

    ScalarViewType basisTanAtBasisEPoints("basisTanAtBasisEPoints",numCells,edgeCardinality, numBasisEPoints);
    ScalarViewType basisTanAtTargetEPoints("basisTanAtTargetEPoints",numCells,edgeCardinality, numTargetEPoints);
    ScalarViewType weightedTanBasisAtBasisEPoints("weightedTanBasisAtBasisEPoints",numCells,edgeCardinality, numBasisEPoints);
    ScalarViewType weightedTanBasisAtTargetEPoints("weightedTanBasisAtTargetEPoints",numCells,edgeCardinality, numTargetEPoints);
    ScalarViewType targetTanAtTargetEPoints("normalTargetAtTargetEPoints",numCells, numTargetEPoints);

    auto targetEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetEvalWeights(edgeDim,ie));
    auto basisEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisEvalWeights(edgeDim,ie));

    //Note: we are not considering the jacobian of the orientation map since it is simply a scalar term for the integrals and it does not affect the projection
    ordinal_type offsetBasis = basisEPointsRange(edgeDim, ie).first;
    ordinal_type offsetTarget = targetEPointsRange(edgeDim, ie).first;

    typedef ComputeBasisCoeffsOnEdges_HCurl<ScalarViewType,  decltype(basisEWeights), decltype(tagToOrdinal), decltype(targetAtTargetEPoints)> functorTypeEdge;
    Kokkos::parallel_for(policy, functorTypeEdge(basisTanAtBasisEPoints,basisAtBasisEPoints,basisEWeights,
        weightedTanBasisAtBasisEPoints, targetEWeights,
        basisAtTargetEPoints, weightedTanBasisAtTargetEPoints, tagToOrdinal,
        targetAtTargetEPoints, targetTanAtTargetEPoints,
        refEdgesTangent, edgeCardinality, offsetBasis,
        offsetTarget, edgeDim,
        dim, ie));

    ScalarViewType edgeMassMat_("edgeMassMat_", numCells, edgeCardinality+1, edgeCardinality+1),
        edgeRhsMat_("rhsMat_", numCells, edgeCardinality+1);

    ScalarViewType eWeights_("eWeights_", numCells, 1, basisEWeights.extent(0)), targetEWeights_("targetEWeights", numCells, 1, targetEWeights.extent(0));
    RealSpaceTools<SpT>::clone(eWeights_, basisEWeights);
    RealSpaceTools<SpT>::clone(targetEWeights_, targetEWeights);

    range_type range_H(0, edgeCardinality);
    range_type range_B(edgeCardinality, edgeCardinality+1);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(edgeMassMat_,Kokkos::ALL(),range_H,range_H), basisTanAtBasisEPoints, weightedTanBasisAtBasisEPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(edgeMassMat_,Kokkos::ALL(),range_H,range_B), basisTanAtBasisEPoints, eWeights_);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(edgeRhsMat_,Kokkos::ALL(),range_H), targetTanAtTargetEPoints, weightedTanBasisAtTargetEPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(edgeRhsMat_,Kokkos::ALL(),range_B), targetTanAtTargetEPoints, targetEWeights_);

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
    ScalarViewType t_("t",numCells, edgeCardinality+1);
    WorkArrayViewType w_("w",numCells, edgeCardinality+1);

    auto edgeDofs = Kokkos::subview(tagToOrdinal, edgeDim, ie, Kokkos::ALL());
    ElemSystem edgeSystem("edgeSystem", false);
    edgeSystem.solve(basisCoeffs, edgeMassMat_, edgeRhsMat_, t_, w_, edgeDofs, edgeCardinality, 1);

    for(ordinal_type i=0; i<edgeCardinality; ++i)
      computedDofs(computedDofsCount++) = cellBasis->getDofOrdinal(edgeDim, ie, i);

  }

  CellTools<SpT>::setSubcellParametrization();
  typename CellTools<SpT>::subcellParamViewType  subcellParamFace;
  if(numFaces>0)
    CellTools<SpT>::getSubcellParametrization(subcellParamFace,  faceDim, cellBasis->getBaseCellTopology());

  Basis<SpT,scalarType,scalarType> *hgradBasis = NULL;
  for(ordinal_type iface=0; iface<numFaces; ++iface) {

    if(cellTopo.getKey() == shards::getCellTopologyData<shards::Hexahedron<8> >()->key)
      hgradBasis = new Basis_HGRAD_QUAD_Cn_FEM<SpT,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key)
      hgradBasis = new Basis_HGRAD_TRI_Cn_FEM<SpT,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else  {
      std::stringstream ss;
      ss << ">>> ERROR (Intrepid2::ProjectionTools::getHCurlBasisCoeffs): "
          << "Method not implemented for basis " << name;
      INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
    }

    ordinal_type numTargetEPoints = range_size(targetEPointsRange(faceDim, iface));
    ordinal_type numTargetCurlEPoints = range_size(targetCurlEPointsRange(faceDim, iface));
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(faceDim, iface));
    ordinal_type numBasisCurlEPoints = range_size(basisCurlEPointsRange(faceDim, iface));

    ordinal_type numFaceDofs = cellBasis->getDofCount(faceDim,iface);

    ScalarViewType hgradBasisGradAtBasisEPoints("hgradBasisGradAtBasisEPoints",hgradBasis->getCardinality(), numBasisEPoints, faceDim);
    ScalarViewType hgradBasisGradAtTargetEPoints("hgradBasisGradAtTargetEPoints",hgradBasis->getCardinality(), numTargetEPoints, faceDim);

    ordinal_type hgradCardinality = hgradBasis->getDofCount(faceDim,0);

    auto refBasisEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisEvalPoints(faceDim, iface));
    auto refTargetEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetEvalPoints(faceDim, iface));
    hgradBasis->getValues(hgradBasisGradAtBasisEPoints, refBasisEPoints, OPERATOR_GRAD);
    hgradBasis->getValues(hgradBasisGradAtTargetEPoints, refTargetEPoints, OPERATOR_GRAD);

    ScalarViewType basisTanAtBasisEPoints("basisTanAtBasisEPoints",numCells,numFaceDofs, numBasisEPoints,dim-1);
    ScalarViewType basisTanAtTargetEPoints("basisTanAtTargetEPoints",numCells,numFaceDofs, numTargetEPoints,dim-1);
    ScalarViewType basisCurlNormalAtBasisCurlEPoints("normaBasisCurlAtBasisEPoints",numCells,numFaceDofs, numBasisCurlEPoints);
    ScalarViewType wNormalBasisCurlAtBasisCurlEPoints("weightedNormalBasisCurlAtBasisEPoints",numCells,numFaceDofs, numBasisCurlEPoints);

    ScalarViewType targetTanAtTargetEPoints("targetTanAtTargetEPoints",numCells, numTargetEPoints, dim-1);
    ScalarViewType normalTargetCurlAtTargetEPoints("normalTargetCurlAtTargetEPoints",numCells, numTargetCurlEPoints);
    ScalarViewType wNormalBasisCurlBasisAtTargetCurlEPoints("weightedNormalBasisCurlAtTargetCurlEPoints",numCells,numFaceDofs, numTargetCurlEPoints);

    ScalarViewType wHgradBasisGradAtBasisEPoints("wHgradBasisGradAtBasisEPoints",numCells, hgradCardinality, numBasisEPoints, faceDim);
    ScalarViewType wHgradBasisGradAtTargetEPoints("wHgradBasisGradAtTargetEPoints",numCells, hgradCardinality, numTargetEPoints, faceDim);

    ScalarViewType negPartialProjCurlNormal("mNormalComputedProjection", numCells,numBasisEPoints);
    ScalarViewType negPartialProjTan("negPartialProjTan", numCells,numBasisEPoints,dim-1);


    ordinal_type offsetBasis = basisEPointsRange(faceDim, iface).first;
    ordinal_type offsetBasisCurl = basisCurlEPointsRange(faceDim, iface).first;
    ordinal_type offsetTarget = targetEPointsRange(faceDim, iface).first;
    ordinal_type offsetTargetCurl = targetCurlEPointsRange(faceDim, iface).first;

    auto hGradTagToOrdinal = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(), hgradBasis->getAllDofOrdinal());

    auto basisEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisEvalWeights(faceDim,iface));
    auto targetEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetEvalWeights(faceDim,iface));
    auto targetCurlEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetDerivEvalWeights(faceDim,iface));
    auto basisCurlEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisDerivEvalWeights(faceDim,iface));
    typedef ComputeBasisCoeffsOnFaces_HCurl<decltype(basisCoeffs), decltype(orts), ScalarViewType,  decltype(targetAtTargetEPoints), decltype(basisEWeights),
        decltype(tagToOrdinal), decltype(refTopologyKey),  decltype(subcellParamFace), decltype(computedDofs)> functorTypeFaces;
    Kokkos::parallel_for(policy, functorTypeFaces(basisCoeffs,
        orts, negPartialProjTan, negPartialProjCurlNormal,
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
        hGradTagToOrdinal, refTopologyKey,
        subcellParamFace,
        computedDofs, offsetBasis,
        offsetBasisCurl, offsetTarget,
        offsetTargetCurl, iface,
        hgradCardinality, numFaces,
        numFaceDofs, numEdgeDofs,
        faceDim, dim));


    ScalarViewType faceMassMat_("faceMassMat_", numCells, numFaceDofs+hgradCardinality, numFaceDofs+hgradCardinality),
        faceRhsMat_("rhsMat_", numCells, numFaceDofs+hgradCardinality);
    range_type range_H(0, numFaceDofs);
    range_type range_B(numFaceDofs, numFaceDofs+hgradCardinality);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(faceMassMat_,Kokkos::ALL(),range_H,range_H), basisCurlNormalAtBasisCurlEPoints, wNormalBasisCurlAtBasisCurlEPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(faceMassMat_,Kokkos::ALL(),range_H,range_B), basisTanAtBasisEPoints, wHgradBasisGradAtBasisEPoints);

    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(faceRhsMat_,Kokkos::ALL(),range_H), normalTargetCurlAtTargetEPoints, wNormalBasisCurlBasisAtTargetCurlEPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(faceRhsMat_,Kokkos::ALL(),range_H), negPartialProjCurlNormal, wNormalBasisCurlAtBasisCurlEPoints,true);

    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(faceRhsMat_,Kokkos::ALL(),range_B), targetTanAtTargetEPoints, wHgradBasisGradAtTargetEPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(faceRhsMat_,Kokkos::ALL(),range_B), negPartialProjTan, wHgradBasisGradAtBasisEPoints,true);


    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
    ScalarViewType t_("t",numCells, numFaceDofs+hgradCardinality);
    WorkArrayViewType w_("w",numCells, numFaceDofs+hgradCardinality);

    auto faceDofs = Kokkos::subview(tagToOrdinal, faceDim, iface, Kokkos::ALL());
    ElemSystem faceSystem( "faceSystem", false);
    faceSystem.solve(basisCoeffs, faceMassMat_, faceRhsMat_, t_, w_, faceDofs, numFaceDofs, hgradCardinality);

    for(ordinal_type i=0; i<numFaceDofs; ++i)
      computedDofs(computedDofsCount++) = cellBasis->getDofOrdinal(faceDim, iface, i);

    delete hgradBasis;
  }

  ordinal_type numCellDofs = cellBasis->getDofCount(dim,0);
  if(numCellDofs>0) {
    if(cellTopo.getKey() == shards::getCellTopologyData<shards::Hexahedron<8> >()->key)
      hgradBasis = new Basis_HGRAD_HEX_Cn_FEM<SpT,scalarType,scalarType>(cellBasis->getDegree());
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key)
      hgradBasis = new Basis_HGRAD_TET_Cn_FEM<SpT,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Triangle<3> >()->key)
      hgradBasis = new Basis_HGRAD_TRI_Cn_FEM<SpT,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Quadrilateral<4> >()->key)
      hgradBasis = new Basis_HGRAD_QUAD_Cn_FEM<SpT,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
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

    hgradBasis->getValues(hgradBasisGradAtBasisEPoints,Kokkos::subview(basisEPoints, 0, basisEPointsRange(dim, 0), Kokkos::ALL()), OPERATOR_GRAD);
    hgradBasis->getValues(hgradBasisGradAtTargetEPoints,Kokkos::subview(targetEPoints, 0, targetEPointsRange(dim, 0), Kokkos::ALL()),OPERATOR_GRAD);

    ScalarViewType cellBasisAtBasisEPoints("basisCellAtEPoints",numCells,numCellDofs, numBasisEPoints, dim);
    ScalarViewType cellBasisCurlAtCurlEPoints("cellBasisCurlAtCurlEPoints",numCells,numCellDofs, numBasisCurlEPoints, derDim);
    ScalarViewType negPartialProjCurl("negPartialProjCurl", numCells, numBasisEPoints, derDim);
    ScalarViewType negPartialProj("negPartialProj", numCells, numBasisEPoints, dim);
    ScalarViewType wBasisCurlAtCurlEPoints("weightedBasisCurlAtBasisEPoints",numCells,numCellDofs, numBasisCurlEPoints,derDim);
    ScalarViewType wBasisCurlBasisAtTargetCurlEPoints("weightedBasisCurlAtTargetCurlEPoints",numCells,numCellDofs, numTargetCurlEPoints,derDim);

    auto targetEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetEvalWeights(dim,0));
    auto basisEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisEvalWeights(dim,0));
    auto targetCurlEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetDerivEvalWeights(dim,0));
    auto basisCurlEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisDerivEvalWeights(dim,0));
    ordinal_type offsetBasis = basisEPointsRange(dim, 0).first;
    ordinal_type offsetBasisCurl = basisCurlEPointsRange(dim, 0).first;
    ordinal_type offsetTargetCurl = targetCurlEPointsRange(dim, 0).first;


    auto hGradTagToOrdinal = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(), hgradBasis->getAllDofOrdinal());

    typedef ComputeBasisCoeffsOnCell_HCurl<decltype(basisCoeffs), ScalarViewType,  decltype(basisEWeights),
        decltype(computedDofs), decltype(tagToOrdinal)> functorTypeCell;
    Kokkos::parallel_for(policy, functorTypeCell(basisCoeffs, negPartialProj, negPartialProjCurl,
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
        numEdgeDofs+numFaceDofs, dim, derDim));

    ScalarViewType cellMassMat_("cellMassMat_", numCells, numCellDofs+hgradCardinality, numCellDofs+hgradCardinality),
        cellRhsMat_("rhsMat_", numCells, numCellDofs+hgradCardinality);

    range_type range_H(0, numCellDofs);
    range_type range_B(numCellDofs, numCellDofs+hgradCardinality);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellMassMat_,Kokkos::ALL(),range_H,range_H), cellBasisCurlAtCurlEPoints, wBasisCurlAtCurlEPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellMassMat_,Kokkos::ALL(),range_H,range_B), cellBasisAtBasisEPoints, wHgradBasisGradAtBasisEPoints);
    if(dim==3)
      FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_H), Kokkos::subview(targetCurlAtTargetCurlEPoints,Kokkos::ALL(),cellCurlPointsRange,Kokkos::ALL()), wBasisCurlBasisAtTargetCurlEPoints);
    else
      FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_H), Kokkos::subview(targetCurlAtTargetCurlEPoints,Kokkos::ALL(),cellCurlPointsRange), Kokkos::subview(wBasisCurlBasisAtTargetCurlEPoints,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),0));
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_H), negPartialProjCurl, wBasisCurlAtCurlEPoints, true);

    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_B), Kokkos::subview(targetAtTargetEPoints,Kokkos::ALL(),cellPointsRange,Kokkos::ALL()), wHgradBasisGradAtTargetEPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_B), negPartialProj, wHgradBasisGradAtBasisEPoints, true);

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
    ScalarViewType t_("t",numCells, numCellDofs+hgradCardinality);
    WorkArrayViewType w_("w",numCells, numCellDofs+hgradCardinality);

    auto cellDofs = Kokkos::subview(tagToOrdinal, dim, 0, Kokkos::ALL());
    ElemSystem cellSystem( "cellSystem", true);
    cellSystem.solve(basisCoeffs, cellMassMat_, cellRhsMat_, t_, w_, cellDofs, numCellDofs, hgradCardinality);

    delete hgradBasis;
  }
}
}
}

#endif

