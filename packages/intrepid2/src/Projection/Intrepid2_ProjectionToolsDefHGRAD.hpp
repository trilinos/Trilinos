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

/** \file   Intrepid2_ProjectionToolsDefHGRAD.hpp
    \brief  Header file for the Intrepid2::Experimental::ProjectionTools
            containing definitions for HGRAD projections.
    \author Created by Mauro Perego
 */

#ifndef __INTREPID2_PROJECTIONTOOLSDEFHGRAD_HPP__
#define __INTREPID2_PROJECTIONTOOLSDEFHGRAD_HPP__

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


namespace Intrepid2 {
namespace Experimental {

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
      basisCoeffs_(ic,idof) = targetAtTargetEPoints_(ic,pt)/basisAtTargetEPoints_(ic,idof,pt,0);
    }
  }
};

template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5, typename ViewType6>
struct ComputeBasisCoeffsOnEdges_HGRAD {
  const ViewType1 basisCoeffs_;
  const ViewType2 negPartialProjGrad_;
  const ViewType2 basisTanAtEPoints_;
  const ViewType2 basisGradAtBasisGradEPoints_;
  const ViewType3 basisGradEWeights_;
  const ViewType2 wBasisAtBasisGradEPoints_;
  const ViewType3 targetGradEWeights_;
  const ViewType2 basisGradAtTargetGradEPoints_;
  const ViewType2 wBasisAtTargetGradEPoints_;
  const ViewType4 computedDofs_;
  const ViewType5 tagToOrdinal_;
  const ViewType2 targetGradTanAtTargetGradEPoints_;
  const ViewType6 targetGradAtTargetGradEPoints_;
  const ViewType2 refEdgesTan_;
  ordinal_type edgeCardinality_;
  ordinal_type offsetBasis_;
  ordinal_type offsetTarget_;
  ordinal_type numVertexDofs_;
  ordinal_type edgeDim_;
  ordinal_type dim_;
  ordinal_type iedge_;

  ComputeBasisCoeffsOnEdges_HGRAD(const ViewType1 basisCoeffs, ViewType2 negPartialProjGrad,  const ViewType2 basisTanAtEPoints,
      const ViewType2 basisGradAtBasisGradEPoints, const ViewType3 basisGradEWeights,  const ViewType2 wBasisAtBasisGradEPoints,   const ViewType3 targetGradEWeights,
      const ViewType2 basisGradAtTargetGradEPoints, const ViewType2 wBasisAtTargetGradEPoints, const ViewType4 computedDofs, const ViewType5 tagToOrdinal,
      const ViewType2 targetGradTanAtTargetGradEPoints, const ViewType6 targetGradAtTargetGradEPoints, const ViewType2 refEdgesTan,
      ordinal_type edgeCardinality, ordinal_type offsetBasis,
      ordinal_type offsetTarget, ordinal_type numVertexDofs, ordinal_type edgeDim, ordinal_type dim, ordinal_type iedge) :
        basisCoeffs_(basisCoeffs), negPartialProjGrad_(negPartialProjGrad), basisTanAtEPoints_(basisTanAtEPoints),
        basisGradAtBasisGradEPoints_(basisGradAtBasisGradEPoints), basisGradEWeights_(basisGradEWeights), wBasisAtBasisGradEPoints_(wBasisAtBasisGradEPoints), targetGradEWeights_(targetGradEWeights),
        basisGradAtTargetGradEPoints_(basisGradAtTargetGradEPoints), wBasisAtTargetGradEPoints_(wBasisAtTargetGradEPoints),
        computedDofs_(computedDofs), tagToOrdinal_(tagToOrdinal), targetGradTanAtTargetGradEPoints_(targetGradTanAtTargetGradEPoints),
        targetGradAtTargetGradEPoints_(targetGradAtTargetGradEPoints), refEdgesTan_(refEdgesTan),
        edgeCardinality_(edgeCardinality), offsetBasis_(offsetBasis),
        offsetTarget_(offsetTarget), numVertexDofs_(numVertexDofs), edgeDim_(edgeDim), dim_(dim), iedge_(iedge)
  {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    ordinal_type numBasisGradEPoints = basisGradEWeights_.extent(0);
    ordinal_type numTargetGradEPoints = targetGradEWeights_.extent(0);
    for(ordinal_type j=0; j <edgeCardinality_; ++j) {
      ordinal_type jdof = tagToOrdinal_(edgeDim_, iedge_, j);

      //Note: we are not considering the orientation map since it is simply results in a scalar term for the integrals and it does not affect the projection
      for(ordinal_type iq=0; iq <numBasisGradEPoints; ++iq) {
        for(ordinal_type d=0; d <dim_; ++d)
          basisTanAtEPoints_(ic,j,iq) += basisGradAtBasisGradEPoints_(ic,jdof,offsetBasis_+iq,d)*refEdgesTan_(iedge_, d);
        wBasisAtBasisGradEPoints_(ic,j,iq) = basisTanAtEPoints_(ic,j,iq)*basisGradEWeights_(iq);
      }
      for(ordinal_type iq=0; iq <numTargetGradEPoints; ++iq) {
        for(ordinal_type d=0; d <dim_; ++d)
          wBasisAtTargetGradEPoints_(ic,j,iq) += basisGradAtTargetGradEPoints_(ic,jdof,offsetTarget_+iq,d)*refEdgesTan_(iedge_, d)*targetGradEWeights_(iq);
      }
    }

    for(ordinal_type iq=0; iq <numTargetGradEPoints; ++iq)
      for(ordinal_type d=0; d <dim_; ++d)
        targetGradTanAtTargetGradEPoints_(ic,iq) += targetGradAtTargetGradEPoints_(ic,offsetTarget_+iq,d)*refEdgesTan_(iedge_, d);

    for(ordinal_type j=0; j <numVertexDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <numBasisGradEPoints; ++iq)
        for(ordinal_type d=0; d <dim_; ++d)
          negPartialProjGrad_(ic,iq) -=  basisCoeffs_(ic,jdof)*basisGradAtBasisGradEPoints_(ic,jdof,offsetBasis_+iq,d)*refEdgesTan_(iedge_, d);
    }
  }
};

template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5, typename ViewType6, typename ViewType7, typename ViewType8>
struct ComputeBasisCoeffsOnFaces_HGRAD {
  const ViewType1 basisCoeffs_;
  const ViewType2 negPartialProjGrad_;
  const ViewType2 faceBasisGradAtGradEPoints_;
  const ViewType2 basisGradAtBasisGradEPoints_;
  const ViewType3 basisGradEWeights_;
  const ViewType2 wBasisGradAtGradEPoints_;
  const ViewType3 targetGradEWeights_;
  const ViewType2 basisGradAtTargetGradEPoints_;
  const ViewType2 wBasisGradBasisAtTargetGradEPoints_;
  const ViewType4 computedDofs_;
  const ViewType5 tagToOrdinal_;
  const ViewType6 orts_;
  const ViewType2 targetGradTanAtTargetGradEPoints_;
  const ViewType7 targetGradAtTargetGradEPoints_;
  const ViewType2 faceCoeff_;
  const ViewType8 faceParametrization_;
  ordinal_type faceCardinality_;
  ordinal_type offsetBasisGrad_;
  ordinal_type offsetTargetGrad_;
  ordinal_type numVertexEdgeDofs_;
  ordinal_type numFaces_;
  ordinal_type faceDim_;
  ordinal_type dim_;
  ordinal_type iface_;
  unsigned topoKey_;

  ComputeBasisCoeffsOnFaces_HGRAD(const ViewType1 basisCoeffs, ViewType2 negPartialProjGrad,  const ViewType2 faceBasisGradAtGradEPoints,
      const ViewType2 basisGradAtBasisGradEPoints, const ViewType3 basisGradEWeights,  const ViewType2 wBasisGradAtGradEPoints,   const ViewType3 targetGradEWeights,
      const ViewType2 basisGradAtTargetGradEPoints, const ViewType2 wBasisGradBasisAtTargetGradEPoints, const ViewType4 computedDofs, const ViewType5 tagToOrdinal,
      const ViewType6 orts, const ViewType2 targetGradTanAtTargetGradEPoints, const ViewType7 targetGradAtTargetGradEPoints,
      const ViewType8 faceParametrization, ordinal_type faceCardinality, ordinal_type offsetBasisGrad,
      ordinal_type offsetTargetGrad, ordinal_type numVertexEdgeDofs, ordinal_type numFaces, ordinal_type faceDim,
      ordinal_type dim, ordinal_type iface, unsigned topoKey) :
        basisCoeffs_(basisCoeffs), negPartialProjGrad_(negPartialProjGrad), faceBasisGradAtGradEPoints_(faceBasisGradAtGradEPoints),
        basisGradAtBasisGradEPoints_(basisGradAtBasisGradEPoints), basisGradEWeights_(basisGradEWeights), wBasisGradAtGradEPoints_(wBasisGradAtGradEPoints),
        targetGradEWeights_(targetGradEWeights),
        basisGradAtTargetGradEPoints_(basisGradAtTargetGradEPoints), wBasisGradBasisAtTargetGradEPoints_(wBasisGradBasisAtTargetGradEPoints),
        computedDofs_(computedDofs), tagToOrdinal_(tagToOrdinal), orts_(orts), targetGradTanAtTargetGradEPoints_(targetGradTanAtTargetGradEPoints),
        targetGradAtTargetGradEPoints_(targetGradAtTargetGradEPoints), faceParametrization_(faceParametrization),
        faceCardinality_(faceCardinality), offsetBasisGrad_(offsetBasisGrad),
        offsetTargetGrad_(offsetTargetGrad), numVertexEdgeDofs_(numVertexEdgeDofs), numFaces_(numFaces),
        faceDim_(faceDim), dim_(dim), iface_(iface), topoKey_(topoKey)
  {}

  void
  KOKKOS_INLINE_FUNCTION
  operator()(const ordinal_type ic) const {

    typename ViewType2::value_type tanData[2*3];
    auto tangents = ViewType2(tanData, 2, 3);

    ordinal_type fOrt[6];
    orts_(ic).getFaceOrientation(fOrt, numFaces_);
    ordinal_type ort = fOrt[iface_];

    ordinal_type numBasisGradEPoints = basisGradEWeights_.extent(0);
    ordinal_type numTargetGradEPoints = targetGradEWeights_.extent(0);
    Impl::OrientationTools::getRefSubcellTangents(tangents, faceParametrization_,topoKey_, iface_, ort);

    for(ordinal_type j=0; j <faceCardinality_; ++j) {
      ordinal_type jdof = tagToOrdinal_(faceDim_, iface_, j);
      for(ordinal_type itan=0; itan <faceDim_; ++itan) {
        for(ordinal_type iq=0; iq <numBasisGradEPoints; ++iq) {
          for(ordinal_type d=0; d <dim_; ++d)
            faceBasisGradAtGradEPoints_(ic,j,iq,itan) += tangents(itan,d)*basisGradAtBasisGradEPoints_(ic,jdof,offsetBasisGrad_+iq,d);
          wBasisGradAtGradEPoints_(ic,j,iq,itan) = faceBasisGradAtGradEPoints_(ic,j,iq,itan) * basisGradEWeights_(iq);
        }
        for(ordinal_type iq=0; iq <numTargetGradEPoints; ++iq) {
          typename ViewType2::value_type tmp=0;
          for(ordinal_type d=0; d <dim_; ++d) {
            tmp += tangents(itan, d)*basisGradAtTargetGradEPoints_(ic,jdof,offsetTargetGrad_+iq,d);
          }
          wBasisGradBasisAtTargetGradEPoints_(ic,j,iq,itan) = tmp * targetGradEWeights_(iq);
        }
      }
    }

    for(ordinal_type d=0; d <dim_; ++d)
      for(ordinal_type itan=0; itan <faceDim_; ++itan)
        for(ordinal_type iq=0; iq <numTargetGradEPoints; ++iq)
          targetGradTanAtTargetGradEPoints_(ic,iq,itan) += tangents(itan, d)*targetGradAtTargetGradEPoints_(ic,offsetTargetGrad_+iq,d);

    for(ordinal_type j=0; j <numVertexEdgeDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <numBasisGradEPoints; ++iq)
        for(ordinal_type d=0; d <dim_; ++d)
          for(ordinal_type itan=0; itan <faceDim_; ++itan)
            negPartialProjGrad_(ic,iq,itan) -=  basisCoeffs_(ic,jdof)*tangents(itan, d)*basisGradAtBasisGradEPoints_(ic,jdof,offsetBasisGrad_+iq,d);
    }
  }
};


template<typename ViewType1, typename ViewType2, typename ViewType3,
typename ViewType4, typename ViewType5>
struct ComputeBasisCoeffsOnCells_HGRAD {
  const ViewType1 basisCoeffs_;
  const ViewType2 negPartialProjGrad_;
  const ViewType2 cellBasisGradAtGradEPoints_;
  const ViewType2 basisGradAtBasisGradEPoints_;
  const ViewType3 basisGradEWeights_;
  const ViewType2 wBasisGradAtGradEPoints_;
  const ViewType3 targetGradEWeights_;
  const ViewType2 basisGradAtTargetGradEPoints_;
  const ViewType2 wBasisGradBasisAtTargetGradEPoints_;
  const ViewType4 computedDofs_;
  const ViewType5 elemDof_;
  ordinal_type dim_;
  ordinal_type numElemDofs_;
  ordinal_type offsetBasisGrad_;
  ordinal_type offsetTargetGrad_;
  ordinal_type numVertexEdgeFaceDofs_;

  ComputeBasisCoeffsOnCells_HGRAD(const ViewType1 basisCoeffs, ViewType2 negPartialProjGrad,  const ViewType2 cellBasisGradAtGradEPoints,
      const ViewType2 basisGradAtBasisGradEPoints, const ViewType3 basisGradEWeights,  const ViewType2 wBasisGradAtGradEPoints,   const ViewType3 targetGradEWeights,
      const ViewType2 basisGradAtTargetGradEPoints, const ViewType2 wBasisGradBasisAtTargetGradEPoints, const ViewType4 computedDofs, const ViewType5 elemDof,
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
          cellBasisGradAtGradEPoints_(ic,j,iq,d) = basisGradAtBasisGradEPoints_(ic,idof,offsetBasisGrad_+iq,d);
          wBasisGradAtGradEPoints_(ic,j,iq,d) = cellBasisGradAtGradEPoints_(ic,j,iq,d) * basisGradEWeights_(iq);
        }
        for(ordinal_type iq=0; iq <numTargetGradEPoints; ++iq) {
          wBasisGradBasisAtTargetGradEPoints_(ic,j,iq,d )= basisGradAtTargetGradEPoints_(ic,idof,offsetTargetGrad_+iq,d) * targetGradEWeights_(iq);
        }
      }
    }
    for(ordinal_type j=0; j <numVertexEdgeFaceDofs_; ++j) {
      ordinal_type jdof = computedDofs_(j);
      for(ordinal_type iq=0; iq <numBasisGradEPoints; ++iq)
        for(ordinal_type d=0; d <dim_; ++d) {
          negPartialProjGrad_(ic,iq,d) -=  basisCoeffs_(ic,jdof)*basisGradAtBasisGradEPoints_(ic,jdof,offsetBasisGrad_+iq,d);
        }
    }
  }
};


template<typename SpT>
template<typename BasisType,
typename ortValueType,       class ...ortProperties>
void
ProjectionTools<SpT>::getHGradEvaluationPoints(typename BasisType::ScalarViewType ePoints,
    typename BasisType::ScalarViewType gradEPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct,
    const EvalPointsType evalPointType) {
  const auto cellTopo = cellBasis->getBaseCellTopology();
  //const auto cellTopoKey = cellBasis->getBaseCellTopology().getKey();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type numCells = ePoints.extent(0);
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;

  ordinal_type numVertices = (cellBasis->getDofCount(0, 0) > 0) ? cellTopo.getVertexCount() : 0;
  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;

  CellTools<SpT>::setSubcellParametrization();
  typename CellTools<SpT>::subcellParamViewType subcellParamEdge,  subcellParamFace;
  if(numEdges>0)
    CellTools<SpT>::getSubcellParametrization(subcellParamEdge,  edgeDim, cellBasis->getBaseCellTopology());
  if(numFaces>0)
    CellTools<SpT>::getSubcellParametrization(subcellParamFace,  faceDim, cellBasis->getBaseCellTopology());

  auto refTopologyKey =  Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTopologyKey());

  auto ePointsRange = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getPointsRange(evalPointType));
  auto gradEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getDerivPointsRange(evalPointType));

  if(numVertices>0) {
    for(ordinal_type iv=0; iv<numVertices; ++iv) {
      auto vertexEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getEvalPoints(0,iv,evalPointType));
      RealSpaceTools<SpT>::clone(Kokkos::subview(ePoints, Kokkos::ALL(),
          ePointsRange(0, iv), Kokkos::ALL()), vertexEPoints);
    }
  }

  for(ordinal_type ie=0; ie<numEdges; ++ie) {

    auto edgeGradEPointsRange = gradEPointsRange(edgeDim, ie);
    auto edgeGradEPoints =  Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getDerivEvalPoints(edgeDim,ie,evalPointType));

    const auto topoKey = refTopologyKey(edgeDim,ie);

    Kokkos::parallel_for
    ("Evaluate Points Edges ",
        Kokkos::RangePolicy<SpT, int> (0, numCells),
        KOKKOS_LAMBDA (const size_t ic) {


      ordinal_type eOrt[12];
      orts(ic).getEdgeOrientation(eOrt, numEdges);
      ordinal_type ort = eOrt[ie];

      Impl::OrientationTools::mapSubcellCoordsToRefCell(Kokkos::subview(gradEPoints,ic,edgeGradEPointsRange,Kokkos::ALL()),
          edgeGradEPoints, subcellParamEdge, topoKey, ie, ort);
    });
  }


  for(ordinal_type iface=0; iface<numFaces; ++iface) {

    auto faceGradEPointsRange = gradEPointsRange(faceDim, iface);
    auto refGradEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getDerivEvalPoints(faceDim,iface,evalPointType));


    const auto topoKey = refTopologyKey(faceDim,iface);
    Kokkos::parallel_for
    ("Evaluate Points Faces ",
        Kokkos::RangePolicy<SpT, int> (0, numCells),
        KOKKOS_LAMBDA (const size_t ic) {

      ordinal_type fOrt[6];
      orts(ic).getFaceOrientation(fOrt, numFaces);
      ordinal_type ort = fOrt[iface];

      Impl::OrientationTools::mapSubcellCoordsToRefCell(Kokkos::subview(gradEPoints,  ic, faceGradEPointsRange, Kokkos::ALL()),
          refGradEPoints, subcellParamFace, topoKey, iface, ort);
    });
  }

  if(cellBasis->getDofCount(dim,0)>0) {
    auto gradPointsRange = gradEPointsRange(dim, 0);
    //auto refGradEPointsRange = range_type(0, range_size(gradPointsRange));
    auto refCellGradEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getDerivEvalPoints(dim,0,evalPointType));
    RealSpaceTools<SpT>::clone(Kokkos::subview(gradEPoints, Kokkos::ALL(), gradPointsRange, Kokkos::ALL()), refCellGradEPoints);
  }
}


template<typename SpT>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<SpT>::getHGradBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtTargetEPoints,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetGradAtTargetGradEPoints,
    const typename BasisType::ScalarViewType targetEPoints,
    const typename BasisType::ScalarViewType targetGradEPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct){

  typedef typename Kokkos::Impl::is_space<SpT>::host_mirror_space::execution_space host_space_type;
  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> ScalarViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type numTotalTargetEPoints(targetAtTargetEPoints.extent(1)),
      numTotalTargetGradEPoints(targetGradAtTargetGradEPoints.extent(1));
  ordinal_type basisCardinality = cellBasis->getCardinality();
  ordinal_type numCells = targetAtTargetEPoints.extent(0);
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;

  ordinal_type numVertices = (cellBasis->getDofCount(0, 0) > 0) ? cellTopo.getVertexCount() : 0;
  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;

  ScalarViewType refEdgesTan("refEdgesTan",  numEdges, dim);
  ScalarViewType refFacesTangents("refFacesTangents", numFaces, dim, faceDim);

  ordinal_type numVertexDofs = numVertices;

  ordinal_type numEdgeDofs(0);
  for(ordinal_type ie=0; ie<numEdges; ++ie)
    numEdgeDofs += cellBasis->getDofCount(edgeDim,ie);

  ordinal_type numFaceDofs(0);
  for(ordinal_type iface=0; iface<numFaces; ++iface)
    numFaceDofs += cellBasis->getDofCount(faceDim,iface);

  Kokkos::View<ordinal_type*> computedDofs("computedDofs",numVertexDofs+numEdgeDofs+numFaceDofs);

  ordinal_type numTotalBasisEPoints = projStruct->getNumBasisEvalPoints(),
      numTotalBasisGradEPoints = projStruct->getNumBasisDerivEvalPoints();
  ScalarViewType basisEPoints("basisEPoints",numCells,numTotalBasisEPoints, dim);
  ScalarViewType basisGradEPoints("basisGradEPoints",numCells,numTotalBasisGradEPoints, dim);
  getHGradEvaluationPoints(basisEPoints, basisGradEPoints, orts, cellBasis, projStruct, EvalPointsType::BASIS);

  ScalarViewType basisAtTargetEPoints("basisAtTargetEPoints",numCells,basisCardinality, numTotalTargetEPoints);
  {
    ScalarViewType nonOrientedBasisAtTargetEPoints("nonOrientedBasisAtTargetEPoints",numCells,basisCardinality, numTotalTargetEPoints);
    for(ordinal_type ic=0; ic<numCells; ++ic) {
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtTargetEPoints,ic,Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(targetEPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
    }
    OrientationTools<SpT>::modifyBasisByOrientation(basisAtTargetEPoints, nonOrientedBasisAtTargetEPoints, orts, cellBasis);
  }

  ScalarViewType basisGradAtBasisGradEPoints;
  ScalarViewType basisGradAtTargetGradEPoints;
  if(numTotalBasisGradEPoints>0) {
    ScalarViewType nonOrientedBasisGradAtTargetGradEPoints, nonOrientedBasisGradAtBasisGradEPoints;

    basisGradAtBasisGradEPoints = ScalarViewType ("basisGradAtBasisGradEPoints",numCells,basisCardinality, numTotalBasisGradEPoints, dim);
    nonOrientedBasisGradAtBasisGradEPoints = ScalarViewType ("nonOrientedBasisGradAtBasisGradEPoints",numCells,basisCardinality, numTotalBasisGradEPoints, dim);
    basisGradAtTargetGradEPoints = ScalarViewType("basisGradAtTargetGradEPoints",numCells,basisCardinality, numTotalTargetGradEPoints, dim);
    nonOrientedBasisGradAtTargetGradEPoints = ScalarViewType("nonOrientedBasisGradAtTargetGradEPoints",numCells,basisCardinality, numTotalTargetGradEPoints, dim);

    for(ordinal_type ic=0; ic<numCells; ++ic) {
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisGradAtBasisGradEPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(basisGradEPoints, ic, Kokkos::ALL(), Kokkos::ALL()),OPERATOR_GRAD);
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisGradAtTargetGradEPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(targetGradEPoints, ic, Kokkos::ALL(), Kokkos::ALL()),OPERATOR_GRAD);
    }

    OrientationTools<SpT>::modifyBasisByOrientation(basisGradAtBasisGradEPoints, nonOrientedBasisGradAtBasisGradEPoints, orts, cellBasis);
    OrientationTools<SpT>::modifyBasisByOrientation(basisGradAtTargetGradEPoints, nonOrientedBasisGradAtTargetGradEPoints, orts, cellBasis);
  }

  auto targetEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetPointsRange());
  auto targetGradEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetDerivPointsRange());

  auto basisGradEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisDerivPointsRange());

  auto refTopologyKey =  Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTopologyKey());

  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(), cellBasis->getAllDofOrdinal());

  CellTools<SpT>::setSubcellParametrization();
  typename CellTools<SpT>::subcellParamViewType  subcellParamFace;
  if(numFaces>0)
    CellTools<SpT>::getSubcellParametrization(subcellParamFace,  faceDim, cellBasis->getBaseCellTopology());

  ordinal_type computedDofsCount = 0;
  for(ordinal_type iv=0; iv<numVertices; ++iv)
    computedDofs(computedDofsCount++) = tagToOrdinal(0, iv, 0);

  const Kokkos::RangePolicy<SpT> policy(0, numCells);
  typedef ComputeBasisCoeffsOnVertices_HGRAD<decltype(basisCoeffs), decltype(tagToOrdinal), decltype(targetEPointsRange),
      decltype(targetAtTargetEPoints), decltype(basisAtTargetEPoints)> functorType;
  Kokkos::parallel_for(policy, functorType(basisCoeffs, tagToOrdinal, targetEPointsRange,
      targetAtTargetEPoints, basisAtTargetEPoints, numVertices));

  for(ordinal_type ie=0; ie<numEdges; ++ie)  {

    ordinal_type edgeCardinality = cellBasis->getDofCount(edgeDim,ie);
    ordinal_type offsetBasis = basisGradEPointsRange(edgeDim, ie).first;
    ordinal_type offsetTarget = targetGradEPointsRange(edgeDim, ie).first;
    ordinal_type numBasisGradEPoints = range_size(basisGradEPointsRange(edgeDim, ie));
    ordinal_type numTargetGradEPoints = range_size(targetGradEPointsRange(edgeDim, ie));
    auto basisGradEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisDerivEvalWeights(edgeDim,ie));
    auto targetGradEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetDerivEvalWeights(edgeDim,ie));

    auto edgeTan = Kokkos::subview(refEdgesTan, ie, Kokkos::ALL());
    auto edgeTanHost = Kokkos::create_mirror_view(edgeTan);
    CellTools<SpT>::getReferenceEdgeTangent(edgeTanHost,ie, cellTopo);
    Kokkos::deep_copy(edgeTan,edgeTanHost);

    ScalarViewType basisTanAtEPoints("basisTanAtEPoints",numCells,edgeCardinality, numBasisGradEPoints);
    ScalarViewType targetGradTanAtTargetGradEPoints("tanBasisAtTargetGradEPoints",numCells, numTargetGradEPoints);
    ScalarViewType wBasisAtBasisGradEPoints("wTanBasisAtBasisGradEPoints",numCells,edgeCardinality, numBasisGradEPoints);
    ScalarViewType wBasisAtTargetGradEPoints("wTanBasisAtTargetGradEPoints",numCells,edgeCardinality, numTargetGradEPoints);
    ScalarViewType negPartialProjGrad("negPartialProjGrad", numCells, numBasisGradEPoints);

    typedef ComputeBasisCoeffsOnEdges_HGRAD<decltype(basisCoeffs), ScalarViewType,  decltype(basisGradEWeights),
        decltype(computedDofs), decltype(tagToOrdinal), decltype(targetGradAtTargetGradEPoints)> functorTypeEdge;

    Kokkos::parallel_for(policy, functorTypeEdge(basisCoeffs, negPartialProjGrad,  basisTanAtEPoints,
        basisGradAtBasisGradEPoints, basisGradEWeights,  wBasisAtBasisGradEPoints,   targetGradEWeights,
        basisGradAtTargetGradEPoints, wBasisAtTargetGradEPoints, computedDofs, tagToOrdinal,
        targetGradTanAtTargetGradEPoints, targetGradAtTargetGradEPoints, refEdgesTan,
        edgeCardinality, offsetBasis,
        offsetTarget, numVertexDofs, edgeDim, dim, ie));

    ScalarViewType edgeMassMat_("edgeMassMat_", numCells, edgeCardinality, edgeCardinality),
        edgeRhsMat_("rhsMat_", numCells, edgeCardinality);

    FunctionSpaceTools<SpT >::integrate(edgeMassMat_, basisTanAtEPoints, wBasisAtBasisGradEPoints);
    FunctionSpaceTools<SpT >::integrate(edgeRhsMat_, targetGradTanAtTargetGradEPoints, wBasisAtTargetGradEPoints);
    FunctionSpaceTools<SpT >::integrate(edgeRhsMat_, negPartialProjGrad, wBasisAtBasisGradEPoints,true);

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
    ScalarViewType t_("t",numCells, edgeCardinality);
    WorkArrayViewType w_("w",numCells, edgeCardinality);

    auto edgeDofs = Kokkos::subview(tagToOrdinal, edgeDim, ie, Kokkos::ALL());

    ElemSystem edgeSystem("edgeSystem", false);
    edgeSystem.solve(basisCoeffs, edgeMassMat_, edgeRhsMat_, t_, w_, edgeDofs, edgeCardinality);

    for(ordinal_type i=0; i<edgeCardinality; ++i)
      computedDofs(computedDofsCount++) = tagToOrdinal(edgeDim, ie, i);
  }

  for(ordinal_type iface=0; iface<numFaces; ++iface) {

    const auto topoKey = refTopologyKey(faceDim,iface);

    ordinal_type faceCardinality = cellBasis->getDofCount(faceDim,iface);

    ordinal_type numGradEPoints = range_size(basisGradEPointsRange(faceDim, iface));
    ordinal_type numTargetGradEPoints = range_size(targetGradEPointsRange(faceDim, iface));
    ordinal_type offsetBasisGrad = basisGradEPointsRange(faceDim, iface).first;
    ordinal_type offsetTargetGrad = targetGradEPointsRange(faceDim, iface).first;
    auto basisGradEWeights =  Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisDerivEvalWeights(faceDim,iface));
    auto targetGradEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetDerivEvalWeights(faceDim,iface));

    ScalarViewType faceBasisGradAtGradEPoints("normaBasisGradAtGradEPoints",numCells,faceCardinality, numGradEPoints,faceDim);
    ScalarViewType wBasisGradAtGradEPoints("wNormalBasisGradAtGradEPoints",numCells,faceCardinality, numGradEPoints,faceDim);
    ScalarViewType wBasisGradBasisAtTargetGradEPoints("wNormalBasisGradAtTargetGradEPoints",numCells,faceCardinality, numTargetGradEPoints,faceDim);
    ScalarViewType targetGradTanAtTargetGradEPoints("targetGradTanAtTargetGradEPoints",numCells, numTargetGradEPoints,faceDim);
    ScalarViewType negPartialProjGrad("mNormalComputedProjection", numCells,numGradEPoints,faceDim);

    typedef ComputeBasisCoeffsOnFaces_HGRAD<decltype(basisCoeffs), ScalarViewType,  decltype(basisGradEWeights),
        decltype(computedDofs), decltype(tagToOrdinal), decltype(orts), decltype(targetGradAtTargetGradEPoints), decltype(subcellParamFace)> functorTypeFace_HGRAD;

    Kokkos::parallel_for(policy, functorTypeFace_HGRAD(basisCoeffs, negPartialProjGrad, faceBasisGradAtGradEPoints,
        basisGradAtBasisGradEPoints, basisGradEWeights, wBasisGradAtGradEPoints, targetGradEWeights,
        basisGradAtTargetGradEPoints,wBasisGradBasisAtTargetGradEPoints, computedDofs, tagToOrdinal,
        orts,targetGradTanAtTargetGradEPoints,targetGradAtTargetGradEPoints,
        subcellParamFace, faceCardinality, offsetBasisGrad,
        offsetTargetGrad, numVertexDofs+numEdgeDofs, numFaces, faceDim,
        dim, iface, topoKey));

    ScalarViewType faceMassMat_("faceMassMat_", numCells, faceCardinality, faceCardinality),
        faceRhsMat_("rhsMat_", numCells, faceCardinality);

    FunctionSpaceTools<SpT >::integrate(faceMassMat_, faceBasisGradAtGradEPoints, wBasisGradAtGradEPoints);

    FunctionSpaceTools<SpT >::integrate(faceRhsMat_, targetGradTanAtTargetGradEPoints, wBasisGradBasisAtTargetGradEPoints);
    FunctionSpaceTools<SpT >::integrate(faceRhsMat_, negPartialProjGrad, wBasisGradAtGradEPoints,true);

    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> faceMassMat("faceMassMat", faceCardinality,faceCardinality);
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> faceRhsMat("faceRhsMat",faceCardinality, 1);

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
    ScalarViewType t_("t",numCells, faceCardinality);
    WorkArrayViewType w_("w",numCells, faceCardinality);

    auto faceDofs = Kokkos::subview(tagToOrdinal, faceDim, iface, Kokkos::ALL());

    ElemSystem faceSystem("faceSystem", false);
    faceSystem.solve(basisCoeffs, faceMassMat_, faceRhsMat_, t_, w_, faceDofs, faceCardinality);

    for(ordinal_type i=0; i<faceCardinality; ++i)
      computedDofs(computedDofsCount++) = tagToOrdinal(faceDim, iface, i);
  }

  ordinal_type numElemDofs = cellBasis->getDofCount(dim,0);
  if(numElemDofs>0) {

    range_type cellTargetGradEPointsRange = targetGradEPointsRange(dim, 0);
    ordinal_type numTargetGradEPoints = range_size(cellTargetGradEPointsRange);
    ordinal_type numGradEPoints = range_size(basisGradEPointsRange(dim, 0));
    ordinal_type offsetBasisGrad = basisGradEPointsRange(dim, 0).first;
    ordinal_type offsetTargetGrad = cellTargetGradEPointsRange.first;
    auto targetGradEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetDerivEvalWeights(dim,0));
    auto basisGradEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisDerivEvalWeights(dim,0));

    ScalarViewType cellBasisGradAtGradEPoints("internalBasisGradAtEPoints",numCells,numElemDofs, numGradEPoints, dim);
    ScalarViewType negPartialProjGrad("negPartialProjGrad", numCells, numGradEPoints, dim);
    ScalarViewType wBasisGradAtGradEPoints("wBasisGradAtGradEPoints",numCells,numElemDofs, numGradEPoints,dim);
    ScalarViewType wBasisGradBasisAtTargetGradEPoints("wBasisGradAtTargetGradEPoints",numCells,numElemDofs, numTargetGradEPoints,dim);

    auto elemDof = Kokkos::subview(tagToOrdinal, dim, 0, Kokkos::ALL());
    typedef ComputeBasisCoeffsOnCells_HGRAD<decltype(basisCoeffs), ScalarViewType,  decltype(basisGradEWeights), decltype(computedDofs), decltype(elemDof)> functorTypeCell_HGRAD;
    Kokkos::parallel_for(policy, functorTypeCell_HGRAD(basisCoeffs, negPartialProjGrad,  cellBasisGradAtGradEPoints,
        basisGradAtBasisGradEPoints, basisGradEWeights,  wBasisGradAtGradEPoints,   targetGradEWeights,
        basisGradAtTargetGradEPoints, wBasisGradBasisAtTargetGradEPoints, computedDofs, elemDof,
        dim, numElemDofs, offsetBasisGrad, offsetTargetGrad, numVertexDofs+numEdgeDofs+numFaceDofs));
    ScalarViewType cellMassMat_("cellMassMat_", numCells, numElemDofs, numElemDofs),
        cellRhsMat_("rhsMat_", numCells, numElemDofs);

    FunctionSpaceTools<SpT >::integrate(cellMassMat_, cellBasisGradAtGradEPoints, wBasisGradAtGradEPoints);
    FunctionSpaceTools<SpT >::integrate(cellRhsMat_, Kokkos::subview(targetGradAtTargetGradEPoints,Kokkos::ALL(),cellTargetGradEPointsRange,Kokkos::ALL()), wBasisGradBasisAtTargetGradEPoints);
    FunctionSpaceTools<SpT >::integrate(cellRhsMat_, negPartialProjGrad, wBasisGradAtGradEPoints, true);

    typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
    ScalarViewType t_("t",numCells, numElemDofs);
    WorkArrayViewType w_("w",numCells, numElemDofs);

    auto cellDofs = Kokkos::subview(tagToOrdinal, dim, 0, Kokkos::ALL());
    ElemSystem cellSystem("cellSystem", true);
    cellSystem.solve(basisCoeffs, cellMassMat_, cellRhsMat_, t_, w_, cellDofs, numElemDofs);
  }
}
}
}

#endif

