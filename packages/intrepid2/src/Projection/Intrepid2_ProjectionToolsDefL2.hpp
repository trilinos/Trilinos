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

/** \file   Intrepid2_ProjectionToolsDefL2.hpp
    \brief  Header file for the Intrepid2::Experimental::ProjectionTools
            containing definitions for L2 projections.
    \author Created by Mauro Perego
 */

#ifndef __INTREPID2_PROJECTIONTOOLSDEFL2_HPP__
#define __INTREPID2_PROJECTIONTOOLSDEFL2_HPP__

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


namespace Intrepid2 {
namespace Experimental {


template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5>
struct ComputeBasisCoeffsOnVertices_L2 {
  ViewType1 basisCoeffs_;
  const ViewType2 tagToOrdinal_;
  const ViewType3 targetEPointsRange_;
  const ViewType4 targetAtTargetEPoints_;
  const ViewType5 basisAtTargetEPoints_;
  ordinal_type numVertices_;


  ComputeBasisCoeffsOnVertices_L2(ViewType1 basisCoeffs, ViewType2 tagToOrdinal, ViewType3 targetEPointsRange,
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
      basisCoeffs_(ic,idof) = targetAtTargetEPoints_(ic,pt)/basisAtTargetEPoints_(ic,idof,pt,0);
    }
  }
};


template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5, typename ViewType6>
struct ComputeBasisCoeffsOnEdges_L2 {
  const ViewType1 basisCoeffs_;
  const ViewType2 negPartialProj_;
  const ViewType2 basisDofDofAtBasisEPoints_;
  const ViewType2 basisAtBasisEPoints_;
  const ViewType3 basisEWeights_;
  const ViewType2 wBasisDofAtBasisEPoints_;
  const ViewType3 targetEWeights_;
  const ViewType2 basisAtTargetEPoints_;
  const ViewType2 wBasisDofAtTargetEPoints_;
  const ViewType4 computedDofs_;
  const ViewType5 tagToOrdinal_;
  const ViewType6 targetAtTargetEPoints_;
  const ViewType2 targetTanAtTargetEPoints_;
  const ViewType2 refEdgesVec_;
  ordinal_type fieldDim_;
  ordinal_type edgeCardinality_;
  ordinal_type offsetBasis_;
  ordinal_type offsetTarget_;
  ordinal_type numVertexDofs_;
  ordinal_type edgeDim_;
  ordinal_type iedge_;

  ComputeBasisCoeffsOnEdges_L2(const ViewType1 basisCoeffs, ViewType2 negPartialProj,  const ViewType2 basisDofDofAtBasisEPoints,
      const ViewType2 basisAtBasisEPoints, const ViewType3 basisEWeights,  const ViewType2 wBasisDofAtBasisEPoints,   const ViewType3 targetEWeights,
      const ViewType2 basisAtTargetEPoints, const ViewType2 wBasisDofAtTargetEPoints, const ViewType4 computedDofs, const ViewType5 tagToOrdinal,
      const ViewType6 targetAtTargetEPoints, const ViewType2 targetTanAtTargetEPoints, const ViewType2 refEdgesVec,
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
        for(ordinal_type d=0; d <fieldDim_; ++d)
          basisDofDofAtBasisEPoints_(ic,j,iq) += basisAtBasisEPoints_(ic,jdof,offsetBasis_+iq,d)*refEdgesVec_(iedge_,d);
        wBasisDofAtBasisEPoints_(ic,j,iq) = basisDofDofAtBasisEPoints_(ic,j,iq)*basisEWeights_(iq);
      }
      for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {
        for(ordinal_type d=0; d <fieldDim_; ++d)
          wBasisDofAtTargetEPoints_(ic,j,iq) += basisAtTargetEPoints_(ic,jdof,offsetTarget_+iq,d)*refEdgesVec_(iedge_,d)*targetEWeights_(iq);
      }
    }

    for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq)
      for(ordinal_type d=0; d <fieldDim_; ++d)
        targetTanAtTargetEPoints_(ic,iq) += targetAtTargetEPoints_(ic,offsetTarget_+iq,d)*refEdgesVec_(iedge_,d);

    for(ordinal_type j=0; j <numVertexDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq)
        for(ordinal_type d=0; d <fieldDim_; ++d)
          negPartialProj_(ic,iq) -=  basisCoeffs_(ic,jdof)*basisAtBasisEPoints_(ic,jdof,offsetBasis_+iq,d)*refEdgesVec_(iedge_,d);
    }
  }
};

template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5, typename ViewType6, typename ViewType7, typename ViewType8>
struct ComputeBasisCoeffsOnFaces_L2 {
  const ViewType1 basisCoeffs_;
  const ViewType2 negPartialProj_;
  const ViewType2 faceBasisDofAtBasisEPoints_;
  const ViewType2 basisAtBasisEPoints_;
  const ViewType3 basisEWeights_;
  const ViewType2 wBasisDofAtBasisEPoints_;
  const ViewType3 targetEWeights_;
  const ViewType2 basisAtTargetEPoints_;
  const ViewType2 wBasisDofAtTargetEPoints_;
  const ViewType4 computedDofs_;
  const ViewType5 tagToOrdinal_;
  const ViewType6 orts_;
  const ViewType7 targetAtTargetEPoints_;
  const ViewType2 targetDofAtTargetEPoints_;
  const ViewType2 faceCoeff_;
  const ViewType8 faceParametrization_;
  ordinal_type fieldDim_;
  ordinal_type faceCardinality_;
  ordinal_type offsetBasis_;
  ordinal_type offsetTarget_;
  ordinal_type numVertexEdgeDofs_;
  ordinal_type numFaces_;
  ordinal_type faceDim_;
  ordinal_type faceDofDim_;
  ordinal_type dim_;
  ordinal_type iface_;
  unsigned topoKey_;
  bool isHCurlBasis_, isHDivBasis_;

  ComputeBasisCoeffsOnFaces_L2(const ViewType1 basisCoeffs, ViewType2 negPartialProj,  const ViewType2 faceBasisDofAtBasisEPoints,
      const ViewType2 basisAtBasisEPoints, const ViewType3 basisEWeights,  const ViewType2 wBasisDofAtBasisEPoints,   const ViewType3 targetEWeights,
      const ViewType2 basisAtTargetEPoints, const ViewType2 wBasisDofAtTargetEPoints, const ViewType4 computedDofs, const ViewType5 tagToOrdinal,
      const ViewType6 orts, const ViewType7 targetAtTargetEPoints, const ViewType2 targetDofAtTargetEPoints, const ViewType2 faceCoeff,
      const ViewType2 faceParametrization, ordinal_type fieldDim, ordinal_type faceCardinality, ordinal_type offsetBasis,
      ordinal_type offsetTarget, ordinal_type numVertexEdgeDofs, ordinal_type numFaces, ordinal_type faceDim, ordinal_type faceDofDim,
      ordinal_type dim, ordinal_type iface, unsigned topoKey, bool isHCurlBasis, bool isHDivBasis) :
        basisCoeffs_(basisCoeffs), negPartialProj_(negPartialProj), faceBasisDofAtBasisEPoints_(faceBasisDofAtBasisEPoints),
        basisAtBasisEPoints_(basisAtBasisEPoints), basisEWeights_(basisEWeights), wBasisDofAtBasisEPoints_(wBasisDofAtBasisEPoints), targetEWeights_(targetEWeights),
        basisAtTargetEPoints_(basisAtTargetEPoints), wBasisDofAtTargetEPoints_(wBasisDofAtTargetEPoints),
        computedDofs_(computedDofs), tagToOrdinal_(tagToOrdinal), orts_(orts), targetAtTargetEPoints_(targetAtTargetEPoints),
        targetDofAtTargetEPoints_(targetDofAtTargetEPoints), faceCoeff_(faceCoeff),
        faceParametrization_(faceParametrization),
        fieldDim_(fieldDim), faceCardinality_(faceCardinality), offsetBasis_(offsetBasis),
        offsetTarget_(offsetTarget), numVertexEdgeDofs_(numVertexEdgeDofs), numFaces_(numFaces),
        faceDim_(faceDim), faceDofDim_(faceDofDim), dim_(dim), iface_(iface), topoKey_(topoKey),
        isHCurlBasis_(isHCurlBasis), isHDivBasis_(isHDivBasis)
  {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    ordinal_type fOrt[6];
    orts_(ic).getFaceOrientation(fOrt, numFaces_);
    ordinal_type ort = fOrt[iface_];
    //Note: we are not considering the jacobian of the orientation map since it is simply a scalar term for the integrals and it does not affect the projection

    typename ViewType3::value_type data[3*3];
    auto tangentsAndNormal = ViewType3(data, dim_, dim_);

    if(isHCurlBasis_ || isHDivBasis_)
    Impl::OrientationTools::getRefSideTangentsAndNormal(tangentsAndNormal, faceParametrization_,topoKey_, iface_, ort);

    if(isHCurlBasis_) {
      for(ordinal_type d=0; d <dim_; ++d)
        for(ordinal_type itan=0; itan <faceDim_; ++itan) {
          faceCoeff_(ic,d,itan) = tangentsAndNormal(itan,d);
        }
    } else if (isHDivBasis_) {
      for(ordinal_type d=0; d <dim_; ++d)
        faceCoeff_(ic,d,0) = tangentsAndNormal(dim_-1,d);
    } else
      faceCoeff_(ic,0,0) = 1;
    for(ordinal_type j=0; j <faceCardinality_; ++j) {
      ordinal_type jdof = tagToOrdinal_(faceDim_, iface_, j);
      for(ordinal_type itan=0; itan <faceDofDim_; ++itan) {
        for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq) {
          for(ordinal_type d=0; d <fieldDim_; ++d)
            faceBasisDofAtBasisEPoints_(ic,j,iq,itan) += faceCoeff_(ic,d, itan)*basisAtBasisEPoints_(ic,jdof,offsetBasis_+iq,d);
          wBasisDofAtBasisEPoints_(ic,j,iq,itan) = faceBasisDofAtBasisEPoints_(ic,j,iq,itan) * basisEWeights_(iq);
        }
        for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {
          typename ViewType2::value_type sum=0;
          for(ordinal_type d=0; d <fieldDim_; ++d)
            sum += faceCoeff_(ic, d, itan)*basisAtTargetEPoints_(ic,jdof,offsetTarget_+iq,d);
          wBasisDofAtTargetEPoints_(ic,j,iq,itan) = sum * targetEWeights_(iq);
        }
      }
    }

    for(ordinal_type d=0; d <fieldDim_; ++d)
      for(ordinal_type itan=0; itan <faceDofDim_; ++itan) {
        for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq)
          targetDofAtTargetEPoints_(ic,iq,itan) += faceCoeff_(ic, d, itan)*targetAtTargetEPoints_(ic,offsetTarget_+iq,d);
      }

    for(ordinal_type j=0; j <numVertexEdgeDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq)
        for(ordinal_type d=0; d <fieldDim_; ++d)
          for(ordinal_type itan=0; itan <faceDofDim_; ++itan)
            negPartialProj_(ic,iq,itan) -=  basisCoeffs_(ic,jdof)*faceCoeff_(ic, d, itan)*basisAtBasisEPoints_(ic,jdof,offsetBasis_+iq,d);
    }
  }
};


template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5>
struct ComputeBasisCoeffsOnCells_L2 {
  const ViewType1 basisCoeffs_;
  const ViewType2 negPartialProj_;
  const ViewType2 internalBasisAtBasisEPoints_;
  const ViewType2 basisAtBasisEPoints_;
  const ViewType3 basisEWeights_;
  const ViewType2 wBasisAtBasisEPoints_;
  const ViewType3 targetEWeights_;
  const ViewType2 basisAtTargetEPoints_;
  const ViewType2 wBasisDofAtTargetEPoints_;
  const ViewType4 computedDofs_;
  const ViewType5 elemDof_;
  ordinal_type fieldDim_;
  ordinal_type numElemDofs_;
  ordinal_type offsetBasis_;
  ordinal_type offsetTarget_;
  ordinal_type numVertexEdgeFaceDofs_;

  ComputeBasisCoeffsOnCells_L2(const ViewType1 basisCoeffs, ViewType2 negPartialProj,  const ViewType2 internalBasisAtBasisEPoints,
      const ViewType2 basisAtBasisEPoints, const ViewType3 basisEWeights,  const ViewType2 wBasisAtBasisEPoints,   const ViewType3 targetEWeights,
      const ViewType2 basisAtTargetEPoints, const ViewType2 wBasisDofAtTargetEPoints, const ViewType4 computedDofs, const ViewType5 elemDof,
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
          internalBasisAtBasisEPoints_(ic,j,iq,d) = basisAtBasisEPoints_(ic,idof,offsetBasis_+iq,d);
          wBasisAtBasisEPoints_(ic,j,iq,d) = internalBasisAtBasisEPoints_(ic,j,iq,d) * basisEWeights_(iq);
        }
        for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {
          wBasisDofAtTargetEPoints_(ic,j,iq,d) = basisAtTargetEPoints_(ic,idof,offsetTarget_+iq,d)* targetEWeights_(iq);
        }
      }
    }
    for(ordinal_type j=0; j < numVertexEdgeFaceDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <ordinal_type(basisEWeights_.extent(0)); ++iq)
        for(ordinal_type d=0; d <fieldDim_; ++d) {
          negPartialProj_(ic,iq,d) -=  basisCoeffs_(ic,jdof)*basisAtBasisEPoints_(ic,jdof,offsetBasis_+iq,d);
        }
    }
  }
};


template<typename SpT>
template<typename BasisType,
typename ortValueType,       class ...ortProperties>
void
ProjectionTools<SpT>::getL2EvaluationPoints(typename BasisType::ScalarViewType ePoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct,
    const EvalPointsType ePointType) {
  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,ortProperties...> ScalarViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  const auto cellTopo = cellBasis->getBaseCellTopology();
  //const auto cellTopoKey = cellBasis->getBaseCellTopology().getKey();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type numCells = ePoints.extent(0);
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;

  ordinal_type numVertices = (cellBasis->getDofCount(0, 0) > 0) ? cellTopo.getVertexCount() : 0;
  ordinal_type numEdges = (cellBasis->getDofCount(edgeDim, 0) > 0) ? cellTopo.getEdgeCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(faceDim, 0) > 0) ? cellTopo.getFaceCount() : 0;
  ordinal_type numVols = (cellBasis->getDofCount(dim, 0) > 0);

  CellTools<SpT>::setSubcellParametrization();

  auto ePointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getPointsRange(ePointType));

  typename CellTools<SpT>::subcellParamViewType subcellParamEdge,  subcellParamFace;
  if(numEdges>0)
    CellTools<SpT>::getSubcellParametrization(subcellParamEdge,  edgeDim, cellTopo);
  if(numFaces>0)
    CellTools<SpT>::getSubcellParametrization(subcellParamFace,  faceDim, cellTopo);

  auto refTopologyKey =  Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTopologyKey());

  ScalarViewType workView("workView", numCells, projStruct->getMaxNumEvalPoints(ePointType), dim-1);

  if(numVertices>0) {
    for(ordinal_type iv=0; iv<numVertices; ++iv) {
      auto vertexEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getEvalPoints(0,iv,ePointType));
      RealSpaceTools<SpT>::clone(Kokkos::subview(ePoints, Kokkos::ALL(),
          ePointsRange(0, iv), Kokkos::ALL()), vertexEPoints);
    }
  }

  for(ordinal_type ie=0; ie<numEdges; ++ie) {
    auto edgePointsRange = ePointsRange(edgeDim, ie);
    auto edgeRefPointsRange = range_type(0, range_size(edgePointsRange));
    auto edgeEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getEvalPoints(edgeDim,ie,ePointType));

    Kokkos::parallel_for
    ("Evaluate Points",
        Kokkos::RangePolicy<SpT, int> (0, numCells),
        KOKKOS_LAMBDA (const size_t ic) {

      ordinal_type eOrt[12];
      orts(ic).getEdgeOrientation(eOrt, numEdges);
      ordinal_type ort = eOrt[ie];

      auto orientedEdgeEPoints = Kokkos::subview(workView, ic, edgeRefPointsRange, range_type(0,edgeDim));

      Impl::OrientationTools::mapToModifiedReference(orientedEdgeEPoints,edgeEPoints,refTopologyKey(edgeDim,ie),ort);
      CellTools<SpT>::mapToReferenceSubcell(Kokkos::subview(ePoints,ic,edgePointsRange,Kokkos::ALL()), orientedEdgeEPoints, subcellParamEdge, edgeDim, ie, dim);
    });
  }

  for(ordinal_type iface=0; iface<numFaces; ++iface) {
    auto facePointsRange = ePointsRange(faceDim, iface);
    auto faceRefPointsRange = range_type(0, range_size(facePointsRange));
    auto faceEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getEvalPoints(faceDim,iface,ePointType));

    Kokkos::parallel_for
    ("Evaluate Points",
        Kokkos::RangePolicy<SpT, int> (0, numCells),
        KOKKOS_LAMBDA (const size_t ic) {
      ordinal_type fOrt[6];
      orts(ic).getFaceOrientation(fOrt, numFaces);
      ordinal_type ort = fOrt[iface];

      auto orientedFaceEPoints = Kokkos::subview(workView, ic, faceRefPointsRange, Kokkos::ALL());

      Impl::OrientationTools::mapToModifiedReference(orientedFaceEPoints,faceEPoints,refTopologyKey(faceDim,iface),ort);
      CellTools<SpT>::mapToReferenceSubcell(Kokkos::subview(ePoints,  ic, facePointsRange, Kokkos::ALL()), orientedFaceEPoints, subcellParamFace, faceDim, iface, dim);
    });
  }


  if(numVols > 0) {
    auto pointsRange = ePointsRange(dim, 0);
    auto cellEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getEvalPoints(dim,0,ePointType));
    RealSpaceTools<SpT>::clone(Kokkos::subview(ePoints, Kokkos::ALL(), pointsRange, Kokkos::ALL()), cellEPoints);
  }
}

template<typename SpT>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<SpT>::getL2BasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtTargetEPoints,
    const typename BasisType::ScalarViewType targetEPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct){

  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> ScalarViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type numTotalTargetEPoints(targetAtTargetEPoints.extent(1));
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

  Kokkos::View<ordinal_type*, SpT> computedDofs("computedDofs", numVertexDofs+numEdgeDofs+numFaceDofs);

  auto targetEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetPointsRange());
  auto basisEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisPointsRange());

  auto refTopologyKey =  Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTopologyKey());

  ordinal_type numTotalBasisEPoints = projStruct->getNumBasisEvalPoints();
  ScalarViewType basisEPoints("basisEPoints",numCells,numTotalBasisEPoints, dim);
  getL2EvaluationPoints(basisEPoints, orts, cellBasis, projStruct, EvalPointsType::BASIS);

  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(), cellBasis->getAllDofOrdinal());

  ScalarViewType basisAtBasisEPoints("basisAtBasisEPoints",numCells,basisCardinality, numTotalBasisEPoints, fieldDim);
  ScalarViewType basisAtTargetEPoints("basisAtTargetEPoints",numCells,basisCardinality, numTotalTargetEPoints, fieldDim);
  {
    if(fieldDim == 1) {
      ScalarViewType nonOrientedBasisAtBasisEPoints("nonOrientedBasisAtBasisEPoints",numCells,basisCardinality, numTotalBasisEPoints);
      ScalarViewType nonOrientedBasisAtTargetEPoints("nonOrientedBasisAtTargetEPoints",numCells,basisCardinality, numTotalTargetEPoints);
      for(ordinal_type ic=0; ic<numCells; ++ic) {
        cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtTargetEPoints,ic,Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(targetEPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
        cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtBasisEPoints,ic,Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(basisEPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
      }
      OrientationTools<SpT>::modifyBasisByOrientation(Kokkos::subview(basisAtBasisEPoints, Kokkos::ALL(), Kokkos::ALL(),
          Kokkos::ALL(),0), nonOrientedBasisAtBasisEPoints, orts, cellBasis);
      OrientationTools<SpT>::modifyBasisByOrientation(Kokkos::subview(basisAtTargetEPoints, Kokkos::ALL(),
          Kokkos::ALL(), Kokkos::ALL(),0), nonOrientedBasisAtTargetEPoints, orts, cellBasis);
    }
    else {
      ScalarViewType nonOrientedBasisAtBasisEPoints("nonOrientedBasisAtBasisEPoints",numCells,basisCardinality, numTotalBasisEPoints,fieldDim);
      ScalarViewType nonOrientedBasisAtTargetEPoints("nonOrientedBasisAtTargetEPoints",numCells,basisCardinality, numTotalTargetEPoints,fieldDim);
      for(ordinal_type ic=0; ic<numCells; ++ic) {
        cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtTargetEPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(targetEPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
        cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtBasisEPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(basisEPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
      }
      OrientationTools<SpT>::modifyBasisByOrientation(basisAtBasisEPoints, nonOrientedBasisAtBasisEPoints, orts, cellBasis);
      OrientationTools<SpT>::modifyBasisByOrientation(basisAtTargetEPoints, nonOrientedBasisAtTargetEPoints, orts, cellBasis);
    }
  }

  {
    ordinal_type computedDofsCount = 0;
    for(ordinal_type iv=0; iv<numVertices; ++iv)
      computedDofs(computedDofsCount++) = tagToOrdinal(0, iv, 0);

    for(ordinal_type ie=0; ie<numEdges; ++ie)  {
      ordinal_type edgeCardinality = cellBasis->getDofCount(edgeDim,ie);
      for(ordinal_type i=0; i<edgeCardinality; ++i)
        computedDofs(computedDofsCount++) = tagToOrdinal(edgeDim, ie, i);
    }

    for(ordinal_type iface=0; iface<numFaces; ++iface)  {
      ordinal_type faceCardinality = cellBasis->getDofCount(faceDim,iface);
      for(ordinal_type i=0; i<faceCardinality; ++i)
        computedDofs(computedDofsCount++) = tagToOrdinal(faceDim, iface, i);
    }
  }

  bool isHGradBasis = (cellBasis->getFunctionSpace() == FUNCTION_SPACE_HGRAD);
  bool isHCurlBasis = (cellBasis->getFunctionSpace() == FUNCTION_SPACE_HCURL);
  bool isHDivBasis = (cellBasis->getFunctionSpace() == FUNCTION_SPACE_HDIV);
  ordinal_type faceDofDim = isHCurlBasis ? 2 : 1;
  ScalarViewType edgeCoeff("edgeCoeff", fieldDim);


  const Kokkos::RangePolicy<SpT> policy(0, numCells);

  if(isHGradBasis) {

    typedef ComputeBasisCoeffsOnVertices_L2<decltype(basisCoeffs), decltype(tagToOrdinal), decltype(targetEPointsRange),
        decltype(targetAtTargetEPoints), decltype(basisAtTargetEPoints)> functorType;
    Kokkos::parallel_for(policy, functorType(basisCoeffs, tagToOrdinal, targetEPointsRange,
        targetAtTargetEPoints, basisAtTargetEPoints, numVertices));
  }

  for(ordinal_type ie=0; ie<numEdges; ++ie)  {
    auto edgeVec = Kokkos::subview(refEdgesVec, ie, Kokkos::ALL());
    auto edgeVecHost = Kokkos::create_mirror_view(edgeVec);

    if(isHCurlBasis) {
      CellTools<SpT>::getReferenceEdgeTangent(edgeVecHost,ie, cellTopo);
    } else if(isHDivBasis) {
      CellTools<SpT>::getReferenceSideNormal(edgeVecHost, ie, cellTopo);
    } else {
      edgeVecHost(0) = 1;
    }
    Kokkos::deep_copy(edgeVec,edgeVecHost);

    ordinal_type edgeCardinality = cellBasis->getDofCount(edgeDim,ie);
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(edgeDim, ie));
    ordinal_type numTargetEPoints = range_size(targetEPointsRange(edgeDim, ie));

    ScalarViewType basisDofAtBasisEPoints("BasisDofAtBasisEPoints",numCells,edgeCardinality, numBasisEPoints);
    ScalarViewType tragetDofAtTargetEPoints("TargetDofAtTargetEPoints",numCells, numTargetEPoints);
    ScalarViewType weightedBasisAtBasisEPoints("weightedTanBasisAtBasisEPoints",numCells,edgeCardinality, numBasisEPoints);
    ScalarViewType weightedBasisAtTargetEPoints("weightedTanBasisAtTargetEPoints",numCells,edgeCardinality, numTargetEPoints);
    ScalarViewType negPartialProj("negPartialProj", numCells, numBasisEPoints);

    auto targetEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetEvalWeights(edgeDim,ie));
    auto basisEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisEvalWeights(edgeDim,ie));

    //Note: we are not considering the jacobian of the orientation map since it is simply a scalar term for the integrals and it does not affect the projection
    ordinal_type offsetBasis = basisEPointsRange(edgeDim, ie).first;
    ordinal_type offsetTarget = targetEPointsRange(edgeDim, ie).first;


    typedef ComputeBasisCoeffsOnEdges_L2<decltype(basisCoeffs), ScalarViewType,  decltype(basisEWeights),
        decltype(computedDofs), decltype(tagToOrdinal), decltype(targetAtTargetEPoints)> functorTypeEdge;

    Kokkos::parallel_for(policy, functorTypeEdge(basisCoeffs, negPartialProj, basisDofAtBasisEPoints,
        basisAtBasisEPoints, basisEWeights, weightedBasisAtBasisEPoints, targetEWeights,
        basisAtTargetEPoints, weightedBasisAtTargetEPoints, computedDofs, tagToOrdinal,
        targetAtTargetEPoints,tragetDofAtTargetEPoints, refEdgesVec, fieldDim,
        edgeCardinality, offsetBasis, offsetTarget, numVertexDofs, edgeDim, ie));


    ScalarViewType edgeMassMat_("edgeMassMat_", numCells, edgeCardinality, edgeCardinality),
        edgeRhsMat_("rhsMat_", numCells, edgeCardinality);

    FunctionSpaceTools<SpT >::integrate(edgeMassMat_, basisDofAtBasisEPoints, weightedBasisAtBasisEPoints);
    FunctionSpaceTools<SpT >::integrate(edgeRhsMat_, tragetDofAtTargetEPoints, weightedBasisAtTargetEPoints);
    FunctionSpaceTools<SpT >::integrate(edgeRhsMat_, negPartialProj, weightedBasisAtBasisEPoints,true);


    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
    ScalarViewType t_("t",numCells, edgeCardinality);
    WorkArrayViewType w_("w",numCells, edgeCardinality);

    auto edgeDof = Kokkos::subview(tagToOrdinal, edgeDim, ie, Kokkos::ALL());

    ElemSystem edgeSystem("edgeSystem", false);
    edgeSystem.solve(basisCoeffs, edgeMassMat_, edgeRhsMat_, t_, w_, edgeDof, edgeCardinality);
  }

  CellTools<SpT>::setSubcellParametrization();
  typename CellTools<SpT>::subcellParamViewType  subcellParamFace;
  if(numFaces>0)
    CellTools<SpT>::getSubcellParametrization(subcellParamFace,  faceDim, cellBasis->getBaseCellTopology());

  ScalarViewType faceCoeff("faceCoeff", numCells, fieldDim, faceDofDim);
  for(ordinal_type iface=0; iface<numFaces; ++iface) {
    const auto topoKey = refTopologyKey(faceDim,iface);
    ordinal_type faceCardinality = cellBasis->getDofCount(faceDim,iface);

    ordinal_type numTargetEPoints = range_size(targetEPointsRange(faceDim, iface));
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(faceDim, iface));

    ScalarViewType faceBasisDofAtBasisEPoints("normaBasisAtBasisEPoints",numCells,faceCardinality, numBasisEPoints,faceDofDim);
    ScalarViewType wBasisDofAtBasisEPoints("weightedNormalBasisAtBasisEPoints",numCells,faceCardinality, numBasisEPoints,faceDofDim);

    ScalarViewType faceBasisAtTargetEPoints("normalBasisAtTargetEPoints",numCells,faceCardinality, numTargetEPoints,faceDofDim);
    ScalarViewType wBasisDofAtTargetEPoints("weightedNormalBasisAtTargetEPoints",numCells,faceCardinality, numTargetEPoints,faceDofDim);

    ScalarViewType targetDofAtTargetEPoints("targetDofAtTargetEPoints",numCells, numTargetEPoints,faceDofDim);
    ScalarViewType negPartialProj("mNormalComputedProjection", numCells,numBasisEPoints,faceDofDim);

    ordinal_type offsetBasis = basisEPointsRange(faceDim, iface).first;
    ordinal_type offsetTarget = targetEPointsRange(faceDim, iface).first;
    auto targetEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetEvalWeights(faceDim,iface));
    auto basisEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisEvalWeights(faceDim,iface));


    typedef ComputeBasisCoeffsOnFaces_L2<decltype(basisCoeffs), ScalarViewType,  decltype(basisEWeights),
        decltype(computedDofs), decltype(tagToOrdinal), decltype(orts), decltype(targetAtTargetEPoints), decltype(subcellParamFace)> functorTypeFace;

    Kokkos::parallel_for(policy, functorTypeFace(basisCoeffs, negPartialProj,faceBasisDofAtBasisEPoints,
        basisAtBasisEPoints, basisEWeights, wBasisDofAtBasisEPoints, targetEWeights,
        basisAtTargetEPoints, wBasisDofAtTargetEPoints, computedDofs, tagToOrdinal,
        orts, targetAtTargetEPoints,targetDofAtTargetEPoints, faceCoeff,
        subcellParamFace, fieldDim, faceCardinality, offsetBasis,
        offsetTarget, numVertexDofs+numEdgeDofs, numFaces, faceDim,faceDofDim,
        dim, iface, topoKey, isHCurlBasis, isHDivBasis));

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
    ScalarViewType faceMassMat_("faceMassMat_", numCells, faceCardinality, faceCardinality),
        faceRhsMat_("rhsMat_", numCells, faceCardinality);

    FunctionSpaceTools<SpT >::integrate(faceMassMat_, faceBasisDofAtBasisEPoints, wBasisDofAtBasisEPoints);
    FunctionSpaceTools<SpT >::integrate(faceRhsMat_, targetDofAtTargetEPoints, wBasisDofAtTargetEPoints);
    FunctionSpaceTools<SpT >::integrate(faceRhsMat_, negPartialProj, wBasisDofAtBasisEPoints,true);

    ScalarViewType t_("t",numCells, faceCardinality);
    WorkArrayViewType w_("w",numCells,faceCardinality);

    auto faceDof = Kokkos::subview(tagToOrdinal, faceDim, iface, Kokkos::ALL());

    ElemSystem faceSystem("faceSystem", false);
    faceSystem.solve(basisCoeffs, faceMassMat_, faceRhsMat_, t_, w_, faceDof, faceCardinality);
  }

  ordinal_type numElemDofs = cellBasis->getDofCount(dim,0);


  if(numElemDofs>0) {

    auto cellDofs = Kokkos::subview(tagToOrdinal, dim, 0, Kokkos::ALL());

    range_type cellPointsRange = targetEPointsRange(dim, 0);

    ordinal_type numTargetEPoints = range_size(targetEPointsRange(dim,0));
    ordinal_type numBasisEPoints = range_size(basisEPointsRange(dim,0));

    ScalarViewType internalBasisAtBasisEPoints("internalBasisAtBasisEPoints",numCells,numElemDofs, numBasisEPoints, fieldDim);
    ScalarViewType negPartialProj("negPartialProj", numCells, numBasisEPoints, fieldDim);

    auto targetEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetEvalWeights(dim,0));
    auto basisEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisEvalWeights(dim,0));
    ordinal_type offsetBasis = basisEPointsRange(dim, 0).first;
    ordinal_type offsetTarget = targetEPointsRange(dim, 0).first;


    ScalarViewType wBasisAtBasisEPoints("weightedBasisAtBasisEPoints",numCells,numElemDofs, numBasisEPoints,fieldDim);
    ScalarViewType wBasisDofAtTargetEPoints("weightedBasisAtTargetEPoints",numCells,numElemDofs, numTargetEPoints,fieldDim);

    typedef ComputeBasisCoeffsOnCells_L2<decltype(basisCoeffs), ScalarViewType,  decltype(basisEWeights), decltype(computedDofs), decltype(cellDofs)> functorType;
    Kokkos::parallel_for(policy, functorType( basisCoeffs, negPartialProj,  internalBasisAtBasisEPoints,
        basisAtBasisEPoints, basisEWeights,  wBasisAtBasisEPoints, targetEWeights, basisAtTargetEPoints, wBasisDofAtTargetEPoints,
        computedDofs, cellDofs, fieldDim, numElemDofs, offsetBasis, offsetTarget, numVertexDofs+numEdgeDofs+numFaceDofs));

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
    ScalarViewType cellMassMat_("cellMassMat_", numCells, numElemDofs, numElemDofs),
        cellRhsMat_("rhsMat_", numCells, numElemDofs);

    FunctionSpaceTools<SpT >::integrate(cellMassMat_, internalBasisAtBasisEPoints, wBasisAtBasisEPoints);
    if(fieldDim==1)
      FunctionSpaceTools<SpT >::integrate(cellRhsMat_, Kokkos::subview(targetAtTargetEPoints,Kokkos::ALL(),cellPointsRange,Kokkos::ALL()),
          Kokkos::subview(wBasisDofAtTargetEPoints,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),0));
    else
      FunctionSpaceTools<SpT >::integrate(cellRhsMat_, Kokkos::subview(targetAtTargetEPoints,Kokkos::ALL(),cellPointsRange,Kokkos::ALL()), wBasisDofAtTargetEPoints);
    FunctionSpaceTools<SpT >::integrate(cellRhsMat_, negPartialProj, wBasisAtBasisEPoints, true);

    ScalarViewType t_("t",numCells, numElemDofs);
    WorkArrayViewType w_("w",numCells,numElemDofs);
    ElemSystem cellSystem("cellSystem", false);
    cellSystem.solve(basisCoeffs, cellMassMat_, cellRhsMat_, t_, w_, cellDofs, numElemDofs);
  }
}

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
          wBasisAtBasisEPoints_(ic,j,iq,d) = basisAtBasisEPoints_(ic,j,iq,d) * basisEWeights_(iq);
        }
        for(ordinal_type iq=0; iq <ordinal_type(targetEWeights_.extent(0)); ++iq) {
          wBasisDofAtTargetEPoints_(ic,j,iq,d) = basisAtTargetEPoints_(ic,j,iq,d)* targetEWeights_(iq);
        }
      }
    }
  }
};

template<typename SpT>
template<typename BasisType>
void
ProjectionTools<SpT>::getL2DGEvaluationPoints(typename BasisType::ScalarViewType ePoints,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct,
    const EvalPointsType ePointType) {

  ordinal_type dim = cellBasis->getBaseCellTopology().getDimension();
  auto cellEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getEvalPoints(dim,0,ePointType));
  RealSpaceTools<SpT>::clone(ePoints, cellEPoints);
}

template<typename SpT>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<SpT>::getL2DGBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtTargetEPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct){

  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> ScalarViewType;
  const auto cellTopo = cellBasis->getBaseCellTopology();

  ordinal_type dim = cellTopo.getDimension();

  auto basisEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),
      projStruct->getEvalPoints(dim,0,EvalPointsType::BASIS));
  auto targetEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),
      projStruct->getEvalPoints(dim,0,EvalPointsType::TARGET));


  ordinal_type numTotalTargetEPoints(targetAtTargetEPoints.extent(1));
  ordinal_type basisCardinality = cellBasis->getCardinality();
  ordinal_type numCells = targetAtTargetEPoints.extent(0);
  const ordinal_type fieldDim = (targetAtTargetEPoints.rank()==2) ? 1 : targetAtTargetEPoints.extent(2);

  ordinal_type numTotalBasisEPoints = projStruct->getNumBasisEvalPoints();
  ScalarViewType basisAtBasisEPoints("basisAtBasisEPoints",numCells,basisCardinality, numTotalBasisEPoints, fieldDim);
  ScalarViewType basisAtTargetEPoints("basisAtTargetEPoints",numCells,basisCardinality, numTotalTargetEPoints, fieldDim);
  {
    if(fieldDim == 1) {
      ScalarViewType nonOrientedBasisAtBasisEPoints("nonOrientedBasisAtBasisEPoints",numCells,basisCardinality, numTotalBasisEPoints);
      ScalarViewType nonOrientedBasisAtTargetEPoints("nonOrientedBasisAtTargetEPoints",numCells,basisCardinality, numTotalTargetEPoints);
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtTargetEPoints,0,Kokkos::ALL(),Kokkos::ALL()), targetEPoints);
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtBasisEPoints,0,Kokkos::ALL(),Kokkos::ALL()), basisEPoints);

      RealSpaceTools<SpT>::clone(nonOrientedBasisAtTargetEPoints, Kokkos::subview(nonOrientedBasisAtTargetEPoints,0,Kokkos::ALL(),Kokkos::ALL()));
      RealSpaceTools<SpT>::clone(nonOrientedBasisAtBasisEPoints, Kokkos::subview(nonOrientedBasisAtBasisEPoints,0,Kokkos::ALL(),Kokkos::ALL()));
      OrientationTools<SpT>::modifyBasisByOrientation(Kokkos::subview(basisAtBasisEPoints, Kokkos::ALL(), Kokkos::ALL(),
          Kokkos::ALL(),0), nonOrientedBasisAtBasisEPoints, orts, cellBasis);
      OrientationTools<SpT>::modifyBasisByOrientation(Kokkos::subview(basisAtTargetEPoints, Kokkos::ALL(),
          Kokkos::ALL(), Kokkos::ALL(),0), nonOrientedBasisAtTargetEPoints, orts, cellBasis);
    }
    else {
      ScalarViewType nonOrientedBasisAtBasisEPoints("nonOrientedBasisAtBasisEPoints",numCells,basisCardinality, numTotalBasisEPoints,fieldDim);
      ScalarViewType nonOrientedBasisAtTargetEPoints("nonOrientedBasisAtTargetEPoints",numCells,basisCardinality, numTotalTargetEPoints,fieldDim);
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtTargetEPoints,0,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), targetEPoints);
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtBasisEPoints,0,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), basisEPoints);

      RealSpaceTools<SpT>::clone(nonOrientedBasisAtTargetEPoints, Kokkos::subview(nonOrientedBasisAtTargetEPoints,0,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()));
      RealSpaceTools<SpT>::clone(nonOrientedBasisAtBasisEPoints, Kokkos::subview(nonOrientedBasisAtBasisEPoints,0,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()));
      OrientationTools<SpT>::modifyBasisByOrientation(basisAtBasisEPoints, nonOrientedBasisAtBasisEPoints, orts, cellBasis);
      OrientationTools<SpT>::modifyBasisByOrientation(basisAtTargetEPoints, nonOrientedBasisAtTargetEPoints, orts, cellBasis);
    }
  }

  const Kokkos::RangePolicy<SpT> policy(0, numCells);
  ordinal_type numElemDofs = cellBasis->getCardinality();

  auto targetEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetEvalWeights(dim,0));
  auto basisEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisEvalWeights(dim,0));

  ScalarViewType wBasisAtBasisEPoints("weightedBasisAtBasisEPoints",numCells,numElemDofs, numTotalBasisEPoints,fieldDim);
  ScalarViewType wBasisDofAtTargetEPoints("weightedBasisAtTargetEPoints",numCells,numElemDofs, numTotalTargetEPoints,fieldDim);

  typedef MultiplyBasisByWeights<decltype(basisAtBasisEPoints), decltype(basisEWeights)> functorType;
  Kokkos::parallel_for( "Multiply basis by weights", policy, functorType(basisAtBasisEPoints, basisEWeights,
      wBasisAtBasisEPoints, targetEWeights,  basisAtTargetEPoints, wBasisDofAtTargetEPoints, fieldDim, numElemDofs));// )){

  typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
  ScalarViewType cellMassMat_("cellMassMat_", numCells, numElemDofs, numElemDofs),
      cellRhsMat_("rhsMat_", numCells, numElemDofs);

  FunctionSpaceTools<SpT >::integrate(cellMassMat_, basisAtBasisEPoints, wBasisAtBasisEPoints);
  if(fieldDim==1)
    FunctionSpaceTools<SpT >::integrate(cellRhsMat_, targetAtTargetEPoints,
        Kokkos::subview(wBasisDofAtTargetEPoints,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),0));
  else
    FunctionSpaceTools<SpT >::integrate(cellRhsMat_, targetAtTargetEPoints, wBasisDofAtTargetEPoints);

  Kokkos::DynRankView<ordinal_type,Kokkos::HostSpace> hCellDofs("cellDoFs", numElemDofs);
  for(ordinal_type i=0; i<numElemDofs; ++i) hCellDofs(i) = i;
  auto cellDofs = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),hCellDofs);

  ScalarViewType t_("t",numCells, numElemDofs);
  WorkArrayViewType w_("w",numCells,numElemDofs);
  ElemSystem cellSystem("cellSystem", false);
  cellSystem.solve(basisCoeffs, cellMassMat_, cellRhsMat_, t_, w_, cellDofs, numElemDofs);
}

template<typename SpT>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType>
void
ProjectionTools<SpT>::getL2DGBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtTargetEPoints,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct){

  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> ScalarViewType;
  const auto cellTopo = cellBasis->getBaseCellTopology();

  ordinal_type dim = cellTopo.getDimension();

  auto basisEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),
      projStruct->getEvalPoints(dim,0,EvalPointsType::BASIS));
  auto targetEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),
      projStruct->getEvalPoints(dim,0,EvalPointsType::TARGET));

  ordinal_type numTotalTargetEPoints(targetAtTargetEPoints.extent(1));
  ordinal_type basisCardinality = cellBasis->getCardinality();
  ordinal_type numCells = targetAtTargetEPoints.extent(0);
  const ordinal_type fieldDim = (targetAtTargetEPoints.rank()==2) ? 1 : targetAtTargetEPoints.extent(2);

  ordinal_type numTotalBasisEPoints = projStruct->getNumBasisEvalPoints();
  ScalarViewType basisAtBasisEPoints("basisAtBasisEPoints",1,basisCardinality, numTotalBasisEPoints, fieldDim);
  ScalarViewType basisAtTargetEPoints("basisAtTargetEPoints",numCells,basisCardinality, numTotalTargetEPoints, fieldDim);
  {
    if(fieldDim == 1) {
      cellBasis->getValues(Kokkos::subview(basisAtTargetEPoints,0,Kokkos::ALL(),Kokkos::ALL(),0), targetEPoints);
      cellBasis->getValues(Kokkos::subview(basisAtBasisEPoints,0,Kokkos::ALL(),Kokkos::ALL(),0), basisEPoints);

      RealSpaceTools<SpT>::clone(basisAtTargetEPoints, Kokkos::subview(basisAtTargetEPoints,0,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()));
    }
    else {
      cellBasis->getValues(Kokkos::subview(basisAtTargetEPoints,0,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), targetEPoints);
      cellBasis->getValues(Kokkos::subview(basisAtBasisEPoints,0,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), basisEPoints);

      RealSpaceTools<SpT>::clone(basisAtTargetEPoints, Kokkos::subview(basisAtTargetEPoints,0,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()));
    }
  }

  const Kokkos::RangePolicy<SpT> policy(0, numCells);
  ordinal_type numElemDofs = cellBasis->getCardinality();

  auto targetEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetEvalWeights(dim,0));
  auto basisEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisEvalWeights(dim,0));

  ScalarViewType wBasisAtBasisEPoints("weightedBasisAtBasisEPoints",1,numElemDofs, numTotalBasisEPoints,fieldDim);
  ScalarViewType wBasisDofAtTargetEPoints("weightedBasisAtTargetEPoints",numCells,numElemDofs, numTotalTargetEPoints,fieldDim);

  for(ordinal_type j=0; j <numElemDofs; ++j) {
    for(ordinal_type d=0; d <fieldDim; ++d) {
      for(ordinal_type iq=0; iq < numTotalBasisEPoints; ++iq)
        wBasisAtBasisEPoints(0,j,iq,d) = basisAtBasisEPoints(0,j,iq,d) * basisEWeights(iq);
      for(ordinal_type iq=0; iq <numTotalTargetEPoints; ++iq) {
        wBasisDofAtTargetEPoints(0,j,iq,d) = basisAtTargetEPoints(0,j,iq,d)* targetEWeights(iq);
      }
    }
  }
  RealSpaceTools<SpT>::clone(wBasisDofAtTargetEPoints, Kokkos::subview(wBasisDofAtTargetEPoints,0,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()));

  typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
  ScalarViewType cellMassMat_("cellMassMat_", 1, numElemDofs, numElemDofs),
      cellRhsMat_("rhsMat_", numCells, numElemDofs);

  FunctionSpaceTools<SpT >::integrate(cellMassMat_, basisAtBasisEPoints, wBasisAtBasisEPoints);
  if(fieldDim==1)
    FunctionSpaceTools<SpT >::integrate(cellRhsMat_, targetAtTargetEPoints,
        Kokkos::subview(wBasisDofAtTargetEPoints,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),0));
  else
    FunctionSpaceTools<SpT >::integrate(cellRhsMat_, targetAtTargetEPoints, wBasisDofAtTargetEPoints);

  Kokkos::DynRankView<ordinal_type,Kokkos::HostSpace> hCellDofs("cellDoFs", numElemDofs);
  for(ordinal_type i=0; i<numElemDofs; ++i) hCellDofs(i) = i;
  auto cellDofs = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),hCellDofs);

  ScalarViewType t_("t",1, numElemDofs);
  WorkArrayViewType w_("w",numCells,numElemDofs);
  ElemSystem cellSystem("cellSystem", true);
  cellSystem.solve(basisCoeffs, cellMassMat_, cellRhsMat_, t_, w_, cellDofs, numElemDofs);
}

}
}

#endif

