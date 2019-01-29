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

template<typename SpT>
template<typename BasisType,
typename ortValueType,       class ...ortProperties>
void
ProjectionTools<SpT>::getHGradEvaluationPoints(typename BasisType::scalarViewType evaluationPoints,
    typename BasisType::scalarViewType extDerivEvaluationPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct,
    const EvalPointsType evalPointType) {
  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> scalarViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type numCells = evaluationPoints.extent(0);
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;

  ordinal_type numVertices = (cellBasis->getDofCount(0, 0) > 0) ? cellTopo.getVertexCount() : 0;
  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;

  Kokkos::View<ordinal_type*> eOrt("eOrt", numEdges), fOrt("fOrt", numFaces);

  if(numVertices>0) {
    //TODO: use lattice to retrieve vertex coordinates.
    scalarViewType dofCoords("dofCoords", cellBasis->getCardinality(), dim);
    cellBasis->getDofCoords(dofCoords);
    for(ordinal_type iv=0; iv<numVertices; ++iv) {
      ordinal_type idof = cellBasis->getDofOrdinal(0, iv, 0);
      for(ordinal_type d=0; d<dim; d++)
        for(ordinal_type ic=0; ic<numCells; ++ic)
          evaluationPoints(ic,iv,d) = dofCoords(idof,d);
    }
  }

  for(ordinal_type ie=0; ie<numEdges; ++ie) {
    range_type edgeGradPointsRange;
    scalarViewType cubPoints;
    if(evalPointType == TARGET) {
      edgeGradPointsRange = projStruct->getTargetDerivPointsRange(edgeDim, ie);
      cubPoints = projStruct->getTargetDerivEvalPoints(edgeDim, ie);
    }
    else {
      edgeGradPointsRange = projStruct->getBasisDerivPointsRange(edgeDim, ie);
      cubPoints = projStruct->getBasisDerivEvalPoints(edgeDim, ie);
    }

    scalarViewType orientedTargetCubPoints("orientedTargetCubPoints", cubPoints.extent(0),edgeDim);

    const auto topoKey = projStruct->getTopologyKey(edgeDim,ie);

    for(ordinal_type ic=0; ic<numCells; ++ic)  {
      orts(ic).getEdgeOrientation(eOrt.data(), numEdges);
      ordinal_type ort = eOrt(ie);
      Impl::OrientationTools::mapToModifiedReference(orientedTargetCubPoints,cubPoints,topoKey,ort);
      CellTools<SpT>::mapToReferenceSubcell(Kokkos::subview(extDerivEvaluationPoints,ic,edgeGradPointsRange,Kokkos::ALL()), orientedTargetCubPoints, edgeDim, ie, cellBasis->getBaseCellTopology());
    }
  }


  for(ordinal_type iface=0; iface<numFaces; ++iface) {

    scalarViewType gradCubPoints;
    range_type faceGradPointsRange;
    if(evalPointType == TARGET) {
      gradCubPoints = projStruct->getTargetDerivEvalPoints(faceDim, iface);
      faceGradPointsRange = projStruct->getTargetDerivPointsRange(faceDim, iface);
    } else {
      gradCubPoints = projStruct->getBasisDerivEvalPoints(faceDim, iface);
      faceGradPointsRange = projStruct->getBasisDerivPointsRange(faceDim, iface);
    }

    scalarViewType faceGradCubPoints("faceGradCubPoints", gradCubPoints.extent(0), faceDim);

    const auto topoKey = projStruct->getTopologyKey(faceDim,iface);
    for(ordinal_type ic=0; ic<numCells; ++ic) {

      orts(ic).getFaceOrientation(fOrt.data(), numFaces);

      ordinal_type ort = fOrt(iface);
      Impl::OrientationTools::mapToModifiedReference(faceGradCubPoints,gradCubPoints,topoKey,ort);
      CellTools<SpT>::mapToReferenceSubcell(Kokkos::subview(extDerivEvaluationPoints,  ic, faceGradPointsRange, Kokkos::ALL()), faceGradCubPoints, faceDim, iface, cellBasis->getBaseCellTopology());
    }
  }

  if(cellBasis->getDofCount(dim,0)>0) {
    range_type cellGradPointsRange;
    scalarViewType gradCubPoints;
    if(evalPointType == TARGET) {
      gradCubPoints = projStruct->getTargetDerivEvalPoints(dim, 0);
      cellGradPointsRange = projStruct->getTargetDerivPointsRange(dim, 0);
    } else {
      gradCubPoints = projStruct->getBasisDerivEvalPoints(dim, 0);
      cellGradPointsRange = projStruct->getBasisDerivPointsRange(dim, 0);
    }
    RealSpaceTools<SpT>::clone(Kokkos::subview(extDerivEvaluationPoints, Kokkos::ALL(), cellGradPointsRange, Kokkos::ALL()), gradCubPoints);
  }
}


template<typename SpT>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<SpT>::getHGradBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtEvalPoints,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetGradAtGradEvalPoints,
    const typename BasisType::scalarViewType evaluationPoints,
    const typename BasisType::scalarViewType extDerivEvaluationPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct){

  typedef typename Kokkos::Impl::is_space<SpT>::host_mirror_space::execution_space host_space_type;
  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> scalarViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type numTotalEvaluationPoints(targetAtEvalPoints.extent(1)),
      numTotalGradEvaluationPoints(targetGradAtGradEvalPoints.extent(1));
  ordinal_type basisCardinality = cellBasis->getCardinality();
  ordinal_type numCells = targetAtEvalPoints.extent(0);
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;

  const std::string& name = cellBasis->getName();

  ordinal_type numVertices = (cellBasis->getDofCount(0, 0) > 0) ? cellTopo.getVertexCount() : 0;
  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;

  Kokkos::View<ordinal_type*> eOrt("eOrt", numEdges);
  Kokkos::View<ordinal_type*> fOrt("fOrt", numFaces);
  scalarViewType refEdgeTan("refEdgeTan",  dim);
  scalarViewType refFaceTangents("refFaceTangents", dim, 2);
  auto refFaceTanU = Kokkos::subview(refFaceTangents, Kokkos::ALL, 0);
  auto refFaceTanV = Kokkos::subview(refFaceTangents, Kokkos::ALL, 1);

  ordinal_type numVertexDofs = numVertices;

  ordinal_type numEdgeDofs(0);
  for(ordinal_type ie=0; ie<numEdges; ++ie)
    numEdgeDofs += cellBasis->getDofCount(edgeDim,ie);

  ordinal_type numFaceDofs(0);
  for(ordinal_type iface=0; iface<numFaces; ++iface)
    numFaceDofs += cellBasis->getDofCount(faceDim,iface);

  Kokkos::View<ordinal_type*> computedDofs("computedDofs",numVertexDofs+numEdgeDofs+numFaceDofs);

  ordinal_type computedDofsCount = 0;

  ordinal_type numTotalCubPoints = projStruct->getNumBasisEvalPoints(),
      numTotalGradCubPoints = projStruct->getNumBasisDerivEvalPoints();
  scalarViewType cubPoints("cubPoints",numCells,numTotalCubPoints, dim);
  scalarViewType gradCubPoints("gradCubPoints",numCells,numTotalGradCubPoints, dim);
  getHGradEvaluationPoints(cubPoints, gradCubPoints, orts, cellBasis, projStruct, BASIS);

  scalarViewType basisAtCubPoints("basisAtCubPoints",numCells,basisCardinality, numTotalCubPoints);
  scalarViewType basisAtTargetCubPoints("basisAtTargetCubPoints",numCells,basisCardinality, numTotalEvaluationPoints);
  {
    scalarViewType nonOrientedBasisAtCubPoints("nonOrientedBasisAtCubPoints",numCells,basisCardinality, numTotalCubPoints);
    scalarViewType nonOrientedBasisAtTargetCubPoints("nonOrientedBasisAtTargetCubPoints",numCells,basisCardinality, numTotalEvaluationPoints);

    for(ordinal_type ic=0; ic<numCells; ++ic) {
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtTargetCubPoints,ic,Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(evaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtCubPoints,ic,Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(cubPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
    }

    OrientationTools<SpT>::modifyBasisByOrientation(basisAtCubPoints, nonOrientedBasisAtCubPoints, orts, cellBasis);
    OrientationTools<SpT>::modifyBasisByOrientation(basisAtTargetCubPoints, nonOrientedBasisAtTargetCubPoints, orts, cellBasis);
  }

  scalarViewType basisGradAtGradCubPoints;
  scalarViewType basisGradAtTargetGradCubPoints;
  if(numTotalGradEvaluationPoints>0) {
    scalarViewType nonOrientedBasisGradAtTargetGradCubPoints, nonOrientedBasisGradAtGradCubPoints;

    basisGradAtGradCubPoints = scalarViewType ("basisGradAtGradCubPoints",numCells,basisCardinality, numTotalGradCubPoints, dim);
    nonOrientedBasisGradAtGradCubPoints = scalarViewType ("nonOrientedBasisGradAtGradCubPoints",numCells,basisCardinality, numTotalGradCubPoints, dim);
    basisGradAtTargetGradCubPoints = scalarViewType("basisGradAtTargetGradCubPoints",numCells,basisCardinality, numTotalGradEvaluationPoints, dim);
    nonOrientedBasisGradAtTargetGradCubPoints = scalarViewType("nonOrientedBasisGradAtTargetGradCubPoints",numCells,basisCardinality, numTotalGradEvaluationPoints, dim);

    for(ordinal_type ic=0; ic<numCells; ++ic) {
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisGradAtGradCubPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(gradCubPoints, ic, Kokkos::ALL(), Kokkos::ALL()),OPERATOR_GRAD);
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisGradAtTargetGradCubPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(extDerivEvaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()),OPERATOR_GRAD);
    }

    OrientationTools<SpT>::modifyBasisByOrientation(basisGradAtGradCubPoints, nonOrientedBasisGradAtGradCubPoints, orts, cellBasis);
    OrientationTools<SpT>::modifyBasisByOrientation(basisGradAtTargetGradCubPoints, nonOrientedBasisGradAtTargetGradCubPoints, orts, cellBasis);
  }

  for(ordinal_type iv=0; iv<numVertices; ++iv) {
    ordinal_type idof = cellBasis->getDofOrdinal(0, iv, 0);
    computedDofs(computedDofsCount++) = idof;
    for(ordinal_type ic=0; ic<numCells; ++ic)
      basisCoeffs(ic,idof) = targetAtEvalPoints(ic,iv);
  }

  for(ordinal_type ie=0; ie<numEdges; ++ie)  {

    ordinal_type edgeCardinality = cellBasis->getDofCount(edgeDim,ie);
    ordinal_type numCubPoints = projStruct->getNumBasisDerivEvalPoints(edgeDim, ie);
    ordinal_type numTargetCubPoints = projStruct->getNumTargetDerivEvalPoints(edgeDim, ie);

    CellTools<SpT>::getReferenceEdgeTangent(refEdgeTan, ie, cellBasis->getBaseCellTopology());

    scalarViewType edgeBasisAtCubPoints("tanBasisAtElemCubPoints",numCells,edgeCardinality, numCubPoints);
    scalarViewType edgeTargetAtTargetCubPoints("tanBasisAtTargetCubPoints",numCells, numTargetCubPoints);
    scalarViewType weightedBasisAtElemCubPoints("weightedTanBasisAtElemCubPoints",numCells,edgeCardinality, numCubPoints);
    scalarViewType weightedBasisAtTargetCubPoints("weightedTanBasisAtTargetCubPoints",numCells,edgeCardinality, numTargetCubPoints);
    scalarViewType mComputedProjection("mComputedProjection", numCells, numCubPoints);

    scalarViewType targetEvalWeights = projStruct->getTargetDerivEvalWeights(edgeDim, ie);
    scalarViewType basisEvalWeights = projStruct->getBasisDerivEvalWeights(edgeDim, ie);

    //Note: we are not considering the jacobian of the orientation map since it is simply a scalar term for the integrals and it does not affect the projection
    ordinal_type offsetBasis = projStruct->getBasisDerivPointsRange(edgeDim, ie).first;
    ordinal_type offsetTarget = projStruct->getTargetDerivPointsRange(edgeDim, ie).first;
    for(ordinal_type j=0; j <edgeCardinality; ++j) {
      ordinal_type jdof = cellBasis->getDofOrdinal(edgeDim, ie, j);
      for(ordinal_type ic=0; ic<numCells; ++ic) {
        for(ordinal_type iq=0; iq <numCubPoints; ++iq) {
          for(ordinal_type d=0; d <dim; ++d)
            edgeBasisAtCubPoints(ic,j,iq) += basisGradAtGradCubPoints(ic,jdof,offsetBasis+iq,d)*refEdgeTan(d);
          weightedBasisAtElemCubPoints(ic,j,iq) = edgeBasisAtCubPoints(ic,j,iq)*basisEvalWeights(iq);
        }
        for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq) {
          for(ordinal_type d=0; d <dim; ++d)
            weightedBasisAtTargetCubPoints(ic,j,iq) += basisGradAtTargetGradCubPoints(ic,jdof,offsetTarget+iq,d)*refEdgeTan(d)*targetEvalWeights(iq);
        }
      }
    }

    for(ordinal_type ic=0; ic<numCells; ++ic)
      for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq)
        for(ordinal_type d=0; d <dim; ++d)
          edgeTargetAtTargetCubPoints(ic,iq) += targetGradAtGradEvalPoints(ic,offsetTarget+iq,d)*refEdgeTan(d);

    for(ordinal_type j=0; j <numVertexDofs; ++j) {
      ordinal_type jdof = computedDofs(j);
      for(ordinal_type ic=0; ic<numCells; ++ic)
        for(ordinal_type iq=0; iq <numCubPoints; ++iq)
          for(ordinal_type d=0; d <dim; ++d)
            mComputedProjection(ic,iq) -=  basisCoeffs(ic,jdof)*basisGradAtGradCubPoints(ic,jdof,offsetBasis+iq,d)*refEdgeTan(d);
    }


    scalarViewType edgeMassMat_("edgeMassMat_", numCells, edgeCardinality, edgeCardinality),
        edgeRhsMat_("rhsMat_", numCells, edgeCardinality);

    scalarViewType cubWeights_("cubWeights_", numCells, 1, basisEvalWeights.extent(0)), targetEvalWeights_("targetEvalWeights", numCells, 1, targetEvalWeights.extent(0));
    RealSpaceTools<SpT>::clone(cubWeights_, basisEvalWeights);
    RealSpaceTools<SpT>::clone(targetEvalWeights_, targetEvalWeights);

    FunctionSpaceTools<SpT >::integrate(edgeMassMat_, edgeBasisAtCubPoints, weightedBasisAtElemCubPoints);
    FunctionSpaceTools<SpT >::integrate(edgeRhsMat_, edgeTargetAtTargetCubPoints, weightedBasisAtTargetCubPoints);
    FunctionSpaceTools<SpT >::integrate(edgeRhsMat_, mComputedProjection, weightedBasisAtElemCubPoints,true);


    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> edgeMassMat("edgeMassMat", edgeCardinality,edgeCardinality);
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> edgeRhsMat("edgeRhsMat",edgeCardinality, 1);

    Teuchos::LAPACK<ordinal_type,funValsValueType> lapack;
    ordinal_type info = 0;
    for(ordinal_type ic=0; ic<numCells; ++ic)  {
      for(ordinal_type i=0; i<edgeCardinality; ++i) {
        edgeRhsMat(i,0) = edgeRhsMat_(ic,i);
        for(ordinal_type j=0; j<edgeCardinality; ++j)
          edgeMassMat(i,j) = edgeMassMat_(ic,i,j);
      }

      lapack.POSV('U', edgeCardinality, 1,
          edgeMassMat.data(),
          edgeMassMat.stride_1(),
          edgeRhsMat.data(),
          edgeRhsMat.stride_1(),
          &info);

      if (info) {
        std::stringstream ss;
        ss << ">>> ERROR (Intrepid::ProjectionTools::getBasisCoeffs): "
            << "LAPACK return with error code: "
            << info;
        INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
      }

      for(ordinal_type i=0; i<edgeCardinality; ++i){
        ordinal_type edge_dof = cellBasis->getDofOrdinal(edgeDim, ie, i);
        basisCoeffs(ic,edge_dof) = edgeRhsMat(i,0);
      }
    }
    for(ordinal_type i=0; i<edgeCardinality; ++i)
      computedDofs(computedDofsCount++) = cellBasis->getDofOrdinal(edgeDim, ie, i);
  }

  scalarViewType ortJacobian("ortJacobian", faceDim, faceDim);

  for(ordinal_type iface=0; iface<numFaces; ++iface) {


    const auto topoKey = projStruct->getTopologyKey(faceDim,iface);


    ordinal_type faceCardinality = cellBasis->getDofCount(faceDim,iface);

    ordinal_type numTargetGradCubPoints = projStruct->getNumTargetDerivEvalPoints(faceDim, iface);
    ordinal_type numGradCubPoints = projStruct->getNumBasisDerivEvalPoints(faceDim, iface);

    CellTools<SpT>::getReferenceFaceTangents(refFaceTanU, refFaceTanV,iface, cellTopo);

    scalarViewType faceBasisGradAtGradCubPoints("normaBasisGradAtGradCubPoints",numCells,faceCardinality, numGradCubPoints,faceDim);
    scalarViewType wBasisGradAtGradCubPoints("weightedNormalBasisGradAtGradCubPoints",numCells,faceCardinality, numGradCubPoints,faceDim);

    scalarViewType faceBasisGradAtTargetGradCubPoints("normalBasisGradAtTargetGradCubPoints",numCells,faceCardinality, numTargetGradCubPoints,faceDim);
    scalarViewType wBasisGradBasisAtTargetGradCubPoints("weightedNormalBasisGradAtTargetGradCubPoints",numCells,faceCardinality, numTargetGradCubPoints,faceDim);

    scalarViewType targetGradAtTargetGradCubPoints("targetGradAtTargetGradCubPoints",numCells, numTargetGradCubPoints,faceDim);
    scalarViewType mComputedProjectionGrad("mNormalComputedProjection", numCells,numGradCubPoints,faceDim);

    ordinal_type offsetBasisGrad = projStruct->getBasisDerivPointsRange(faceDim, iface).first;
    ordinal_type offsetTargetGrad = projStruct->getTargetDerivPointsRange(faceDim, iface).first;
    scalarViewType targetGradCubWeights = projStruct->getTargetDerivEvalWeights(faceDim, iface);
    scalarViewType gradCubWeights = projStruct->getBasisDerivEvalWeights(faceDim, iface);

    //Note: we are not considering the jacobian of the orientation map since it is simply a scalar term for the integrals and it does not affect the projection
    for(ordinal_type ic=0; ic<numCells; ++ic)  {

      orts(ic).getFaceOrientation(fOrt.data(), numFaces);

      ordinal_type ort = fOrt(iface);
      Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, topoKey, ort);
      for(ordinal_type j=0; j <faceCardinality; ++j) {
        ordinal_type jdof = cellBasis->getDofOrdinal(faceDim, iface, j);
        for(ordinal_type itan=0; itan <faceDim; ++itan) {
          for(ordinal_type iq=0; iq <numGradCubPoints; ++iq) {
            for(ordinal_type d=0; d <dim; ++d)
              for(ordinal_type jtan=0; jtan <faceDim; ++jtan)
                faceBasisGradAtGradCubPoints(ic,j,iq,itan) += refFaceTangents(d, jtan)*ortJacobian(jtan,itan)*basisGradAtGradCubPoints(ic,jdof,offsetBasisGrad+iq,d);
            wBasisGradAtGradCubPoints(ic,j,iq,itan) = faceBasisGradAtGradCubPoints(ic,j,iq,itan) * gradCubWeights(iq);
          }
          for(ordinal_type iq=0; iq <numTargetGradCubPoints; ++iq) {
            for(ordinal_type d=0; d <dim; ++d)
              for(ordinal_type jtan=0; jtan <faceDim; ++jtan)
                faceBasisGradAtTargetGradCubPoints(ic,j,iq,itan) += refFaceTangents(d, jtan)*ortJacobian(jtan,itan)*basisGradAtTargetGradCubPoints(ic,jdof,offsetTargetGrad+iq,d);
            wBasisGradBasisAtTargetGradCubPoints(ic,j,iq,itan) = faceBasisGradAtTargetGradCubPoints(ic,j,iq,itan) * targetGradCubWeights(iq);
          }
        }
      }

      for(ordinal_type d=0; d <dim; ++d)
        for(ordinal_type itan=0; itan <faceDim; ++itan)
          for(ordinal_type iq=0; iq <numTargetGradCubPoints; ++iq)
            for(ordinal_type jtan=0; jtan <faceDim; ++jtan)
              targetGradAtTargetGradCubPoints(ic,iq,itan) += refFaceTangents(d, jtan)*ortJacobian(jtan,itan)*targetGradAtGradEvalPoints(ic,offsetTargetGrad+iq,d);

      for(ordinal_type j=0; j <numVertexDofs+numEdgeDofs; ++j) {
        ordinal_type jdof = computedDofs(j);
        for(ordinal_type iq=0; iq <numGradCubPoints; ++iq)
          for(ordinal_type d=0; d <dim; ++d)
            for(ordinal_type itan=0; itan <faceDim; ++itan)
              for(ordinal_type jtan=0; jtan <faceDim; ++jtan)
                mComputedProjectionGrad(ic,iq,itan) -=  basisCoeffs(ic,jdof)*refFaceTangents(d, jtan)*ortJacobian(jtan,itan)*basisGradAtGradCubPoints(ic,jdof,offsetBasisGrad+iq,d);
      }
    }

    scalarViewType faceMassMat_("faceMassMat_", numCells, faceCardinality, faceCardinality),
        faceRhsMat_("rhsMat_", numCells, faceCardinality);

    FunctionSpaceTools<SpT >::integrate(faceMassMat_, faceBasisGradAtGradCubPoints, wBasisGradAtGradCubPoints);

    FunctionSpaceTools<SpT >::integrate(faceRhsMat_, targetGradAtTargetGradCubPoints, wBasisGradBasisAtTargetGradCubPoints);
    FunctionSpaceTools<SpT >::integrate(faceRhsMat_, mComputedProjectionGrad, wBasisGradAtGradCubPoints,true);

    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> faceMassMat("faceMassMat", faceCardinality,faceCardinality);
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> faceRhsMat("faceRhsMat",faceCardinality, 1);

    Teuchos::LAPACK<ordinal_type,funValsValueType> lapack;
    ordinal_type info = 0;
    for(ordinal_type ic=0; ic<numCells; ++ic)  {

      for(ordinal_type i=0; i<faceCardinality; ++i) {
        faceRhsMat(i,0) = faceRhsMat_(ic,i);
        for(ordinal_type j=0; j<faceCardinality; ++j)
          faceMassMat(i,j) = faceMassMat_(ic,i,j);
      }

      lapack.POSV('U', faceCardinality, 1,
          faceMassMat.data(),
          faceMassMat.stride_1(),
          faceRhsMat.data(),
          faceRhsMat.stride_1(),
          &info);

      if (info) {
        std::stringstream ss;
        ss << ">>> ERROR (Intrepid::ProjectionTools::getBasisCoeffs): "
            << "LAPACK return with error code: "
            << info;
        INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
      }

      for(ordinal_type i=0; i<faceCardinality; ++i){
        ordinal_type face_dof = cellBasis->getDofOrdinal(faceDim, iface, i);
        basisCoeffs(ic,face_dof) = faceRhsMat(i,0);
      }
    }

    for(ordinal_type i=0; i<faceCardinality; ++i)
      computedDofs(computedDofsCount++) = cellBasis->getDofOrdinal(faceDim, iface, i);
  }

  ordinal_type numElemDofs = cellBasis->getDofCount(dim,0);
  if(numElemDofs>0) {

    range_type cellGradPointsRange = projStruct->getTargetDerivPointsRange(dim, 0);

    ordinal_type numTargetGradCubPoints = projStruct->getNumTargetDerivEvalPoints(dim,0);
    ordinal_type numGradCubPoints = projStruct->getNumBasisDerivEvalPoints(dim,0);

    scalarViewType internalBasisGradAtGradCubPoints("internalBasisGradAtCubPoints",numCells,numElemDofs, numGradCubPoints, dim);
    scalarViewType internalBasisGradAtTargetGradCubPoints("weightedBasisGradAtGradCubPoints",numCells,numElemDofs, numTargetGradCubPoints,dim);
    scalarViewType mComputedProjectionGrad("mComputedProjectionGrad", numCells, numGradCubPoints, dim);

    scalarViewType targetGradCubWeights = projStruct->getTargetDerivEvalWeights(dim, 0);
    scalarViewType cubGradWeights = projStruct->getBasisDerivEvalWeights(dim, 0);
    ordinal_type offsetBasisGrad = projStruct->getBasisDerivPointsRange(dim, 0).first;
    ordinal_type offsetTargetGrad = projStruct->getTargetDerivPointsRange(dim, 0).first;


    scalarViewType wBasisGradAtGradCubPoints("weightedBasisGradAtGradCubPoints",numCells,numElemDofs, numGradCubPoints,dim);
    scalarViewType wBasisGradBasisAtTargetGradCubPoints("weightedBasisGradAtTargetGradCubPoints",numCells,numElemDofs, numTargetGradCubPoints,dim);
    for(ordinal_type j=0; j <numElemDofs; ++j) {
      ordinal_type idof = cellBasis->getDofOrdinal(dim, 0, j);
      for(ordinal_type ic=0; ic<numCells; ++ic) {
        for(ordinal_type d=0; d <dim; ++d) {
          for(ordinal_type iq=0; iq <numGradCubPoints; ++iq) {
            internalBasisGradAtGradCubPoints(ic,j,iq,d) = basisGradAtGradCubPoints(ic,idof,offsetBasisGrad+iq,d);
            wBasisGradAtGradCubPoints(ic,j,iq,d) = internalBasisGradAtGradCubPoints(ic,j,iq,d) * cubGradWeights(iq);
          }
          for(ordinal_type iq=0; iq <numTargetGradCubPoints; ++iq) {
            internalBasisGradAtTargetGradCubPoints(ic,j,iq,d) = basisGradAtTargetGradCubPoints(ic,idof,offsetTargetGrad+iq,d);
            wBasisGradBasisAtTargetGradCubPoints(ic,j,iq,d )= internalBasisGradAtTargetGradCubPoints(ic,j,iq,d) * targetGradCubWeights(iq);
          }
        }
      }
    }

    for(ordinal_type j=0; j <numVertexDofs+numEdgeDofs+numFaceDofs; ++j) {
      ordinal_type jdof = computedDofs(j);
      for(ordinal_type iq=0; iq <numGradCubPoints; ++iq)
        for(ordinal_type ic=0; ic<numCells; ++ic)  {
          for(ordinal_type d=0; d <dim; ++d) {
            mComputedProjectionGrad(ic,iq,d) -=  basisCoeffs(ic,jdof)*basisGradAtGradCubPoints(ic,jdof,offsetBasisGrad+iq,d);
          }
        }
    }


    scalarViewType cellMassMat_("cellMassMat_", numCells, numElemDofs, numElemDofs),
        cellRhsMat_("rhsMat_", numCells, numElemDofs);

    FunctionSpaceTools<SpT >::integrate(cellMassMat_, internalBasisGradAtGradCubPoints, wBasisGradAtGradCubPoints);
    FunctionSpaceTools<SpT >::integrate(cellRhsMat_, Kokkos::subview(targetGradAtGradEvalPoints,Kokkos::ALL(),cellGradPointsRange,Kokkos::ALL()), wBasisGradBasisAtTargetGradCubPoints);
    FunctionSpaceTools<SpT >::integrate(cellRhsMat_, mComputedProjectionGrad, wBasisGradAtGradCubPoints, true);

    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> cellMassMat("cellMassMat", numElemDofs,numElemDofs);
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> cellRhsMat("cellRhsMat",numElemDofs, 1);

    Teuchos::LAPACK<ordinal_type,funValsValueType> lapack;
    ordinal_type info = 0;
    for(ordinal_type ic=0; ic<numCells; ++ic) {
      for(ordinal_type i=0; i<numElemDofs; ++i) {
        cellRhsMat(i,0) = cellRhsMat_(ic,i);
        for(ordinal_type j=0; j<numElemDofs; ++j)
          cellMassMat(i,j) = cellMassMat_(ic,i,j);
      }

      lapack.POSV('U', numElemDofs, 1,
          cellMassMat.data(),
          cellMassMat.stride_1(),
          cellRhsMat.data(),
          cellRhsMat.stride_1(),
          &info);

      for(ordinal_type i=0; i<numElemDofs; ++i){
        ordinal_type idof = cellBasis->getDofOrdinal(dim, 0, i);
        basisCoeffs(ic,idof) = cellRhsMat(i,0);
      }

      if (info) {
        std::stringstream ss;
        ss << ">>> ERROR (Intrepid::ProjectionTools::getBasisCoeffs): "
            << "LAPACK return with error code: "
            << info;
        INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
      }
    }
  }
}
}
}

#endif

