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

template<typename SpT>
template<typename BasisType,
typename ortValueType,       class ...ortProperties>
void
ProjectionTools<SpT>::getL2EvaluationPoints(typename BasisType::scalarViewType evaluationPoints,
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
    range_type edgePointsRange;
    scalarViewType cubPoints;
    if(evalPointType == TARGET) {
      edgePointsRange = projStruct->getTargetPointsRange(edgeDim, ie);
      cubPoints = projStruct->getTargetEvalPoints(edgeDim, ie);
    }
    else {
      edgePointsRange = projStruct->getBasisPointsRange(edgeDim, ie);
      cubPoints = projStruct->getBasisEvalPoints(edgeDim, ie);
    }

    scalarViewType orientedTargetCubPoints("orientedTargetCubPoints", cubPoints.extent(0),edgeDim);

    const auto topoKey = projStruct->getTopologyKey(edgeDim,ie);

    for(ordinal_type ic=0; ic<numCells; ++ic)  {
      orts(ic).getEdgeOrientation(eOrt.data(), numEdges);
      ordinal_type ort = eOrt(ie);
      Impl::OrientationTools::mapToModifiedReference(orientedTargetCubPoints,cubPoints,topoKey,ort);
      CellTools<SpT>::mapToReferenceSubcell(Kokkos::subview(evaluationPoints,ic,edgePointsRange,Kokkos::ALL()), orientedTargetCubPoints, edgeDim, ie, cellBasis->getBaseCellTopology());
    }
  }


  for(ordinal_type iface=0; iface<numFaces; ++iface) {

    scalarViewType cubPoints;
    range_type facePointsRange;
    if(evalPointType == TARGET) {
      cubPoints = projStruct->getTargetEvalPoints(faceDim, iface);
      facePointsRange = projStruct->getTargetPointsRange(faceDim, iface);
    } else {
      cubPoints = projStruct->getBasisEvalPoints(faceDim, iface);
      facePointsRange = projStruct->getBasisPointsRange(faceDim, iface);
    }

    scalarViewType faceCubPoints("faceCubPoints", cubPoints.extent(0), faceDim);

    const auto topoKey = projStruct->getTopologyKey(faceDim,iface);
    for(ordinal_type ic=0; ic<numCells; ++ic) {

      orts(ic).getFaceOrientation(fOrt.data(), numFaces);

      ordinal_type ort = fOrt(iface);
      Impl::OrientationTools::mapToModifiedReference(faceCubPoints,cubPoints,topoKey,ort);
      CellTools<SpT>::mapToReferenceSubcell(Kokkos::subview(evaluationPoints,  ic, facePointsRange, Kokkos::ALL()), faceCubPoints, faceDim, iface, cellBasis->getBaseCellTopology());
    }
  }

  if(cellBasis->getDofCount(dim,0)>0) {
    range_type cellPointsRange;
    scalarViewType cubPoints;
    if(evalPointType == TARGET) {
      cubPoints = projStruct->getTargetEvalPoints(dim, 0);
      cellPointsRange = projStruct->getTargetPointsRange(dim, 0);
    } else {
      cubPoints = projStruct->getBasisEvalPoints(dim, 0);
      cellPointsRange = projStruct->getBasisPointsRange(dim, 0);
    }
    RealSpaceTools<SpT>::clone(Kokkos::subview(evaluationPoints, Kokkos::ALL(), cellPointsRange, Kokkos::ALL()), cubPoints);
  }
}


template<typename SpT>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<SpT>::getL2BasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtEvalPoints,
    const typename BasisType::scalarViewType evaluationPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct){

  typedef typename Kokkos::Impl::is_space<SpT>::host_mirror_space::execution_space host_space_type;
  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> scalarViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type numTotalEvaluationPoints(targetAtEvalPoints.extent(1));
  ordinal_type basisCardinality = cellBasis->getCardinality();
  ordinal_type numCells = targetAtEvalPoints.extent(0);
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;
  const ordinal_type fieldDim = (targetAtEvalPoints.rank()==2) ? 1 : targetAtEvalPoints.extent(2);

  const std::string& name = cellBasis->getName();

  ordinal_type numVertices = (cellBasis->getDofCount(0, 0) > 0) ? cellTopo.getVertexCount() : 0;
  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;

  Kokkos::View<ordinal_type*> eOrt("eOrt", numEdges);
  Kokkos::View<ordinal_type*> fOrt("fOrt", numFaces);
  scalarViewType refEdgeTan("refEdgeTan",  dim);
  scalarViewType refEdgeNormal("refEdgeNormal",  dim);
  scalarViewType refFaceTangents("refFaceTangents", dim, 2);
  scalarViewType refFaceNormal("refFaceNormal", dim);
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

  ordinal_type numTotalCubPoints = projStruct->getNumBasisEvalPoints();
  scalarViewType cubPoints("cubPoints",numCells,numTotalCubPoints, dim);
  getL2EvaluationPoints(cubPoints, orts, cellBasis, projStruct, BASIS);

  scalarViewType basisAtCubPoints("basisAtCubPoints",numCells,basisCardinality, numTotalCubPoints, fieldDim);
  scalarViewType basisAtTargetCubPoints("basisAtTargetCubPoints",numCells,basisCardinality, numTotalEvaluationPoints, fieldDim);
  {
    if(fieldDim == 1) {
      scalarViewType nonOrientedBasisAtCubPoints("nonOrientedBasisAtCubPoints",numCells,basisCardinality, numTotalCubPoints);
      scalarViewType nonOrientedBasisAtTargetCubPoints("nonOrientedBasisAtTargetCubPoints",numCells,basisCardinality, numTotalEvaluationPoints);
      for(ordinal_type ic=0; ic<numCells; ++ic) {
        cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtTargetCubPoints,ic,Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(evaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
        cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtCubPoints,ic,Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(cubPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
      }
      OrientationTools<SpT>::modifyBasisByOrientation(Kokkos::subview(basisAtCubPoints, Kokkos::ALL(), Kokkos::ALL(),
          Kokkos::ALL(),0), nonOrientedBasisAtCubPoints, orts, cellBasis);
      OrientationTools<SpT>::modifyBasisByOrientation(Kokkos::subview(basisAtTargetCubPoints, Kokkos::ALL(),
          Kokkos::ALL(), Kokkos::ALL(),0), nonOrientedBasisAtTargetCubPoints, orts, cellBasis);

    }
    else {
      scalarViewType nonOrientedBasisAtCubPoints("nonOrientedBasisAtCubPoints",numCells,basisCardinality, numTotalCubPoints,fieldDim);
      scalarViewType nonOrientedBasisAtTargetCubPoints("nonOrientedBasisAtTargetCubPoints",numCells,basisCardinality, numTotalEvaluationPoints,fieldDim);
      for(ordinal_type ic=0; ic<numCells; ++ic) {
        cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtTargetCubPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(evaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
        cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtCubPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(cubPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
      }
      OrientationTools<SpT>::modifyBasisByOrientation(basisAtCubPoints, nonOrientedBasisAtCubPoints, orts, cellBasis);
      OrientationTools<SpT>::modifyBasisByOrientation(basisAtTargetCubPoints, nonOrientedBasisAtTargetCubPoints, orts, cellBasis);
    }
  }

  for(ordinal_type iv=0; iv<numVertices; ++iv) {
    ordinal_type idof = cellBasis->getDofOrdinal(0, iv, 0);
    computedDofs(computedDofsCount++) = idof;
    for(ordinal_type ic=0; ic<numCells; ++ic)
      basisCoeffs(ic,idof) = targetAtEvalPoints(ic,iv);
  }

  bool isHCurlBAsis = name.find("HCURL") != std::string::npos;

  ordinal_type faceDofDim = isHCurlBAsis ? 2 : 1;

  scalarViewType edgeCoeff("edgeCoeff", fieldDim);
  for(ordinal_type ie=0; ie<numEdges; ++ie)  {

    if(fieldDim == 1)
      edgeCoeff(0) = 1;
    else if(isHCurlBAsis) {
      CellTools<SpT>::getReferenceEdgeTangent(refEdgeTan,ie, cellTopo);
      Kokkos::deep_copy(edgeCoeff,refEdgeTan);
    } else {
      CellTools<SpT>::getReferenceSideNormal(refEdgeNormal, ie, cellTopo);
      Kokkos::deep_copy(edgeCoeff,refEdgeNormal);
    }

    ordinal_type edgeCardinality = cellBasis->getDofCount(edgeDim,ie);
    ordinal_type numCubPoints = projStruct->getNumBasisEvalPoints(edgeDim, ie);
    ordinal_type numTargetCubPoints = projStruct->getNumTargetEvalPoints(edgeDim, ie);

    scalarViewType edgeBasisAtCubPoints("tanBasisAtElemCubPoints",numCells,edgeCardinality, numCubPoints);
    scalarViewType edgeTargetAtTargetCubPoints("tanBasisAtTargetCubPoints",numCells, numTargetCubPoints);
    scalarViewType weightedBasisAtElemCubPoints("weightedTanBasisAtElemCubPoints",numCells,edgeCardinality, numCubPoints);
    scalarViewType weightedBasisAtTargetCubPoints("weightedTanBasisAtTargetCubPoints",numCells,edgeCardinality, numTargetCubPoints);
    scalarViewType mComputedProjection("mComputedProjection", numCells, numCubPoints);

    scalarViewType targetEvalWeights = projStruct->getTargetEvalWeights(edgeDim, ie);
    scalarViewType basisEvalWeights = projStruct->getBasisEvalWeights(edgeDim, ie);

    //Note: we are not considering the jacobian of the orientation map since it is simply a scalar term for the integrals and it does not affect the projection
    ordinal_type offsetBasis = projStruct->getBasisPointsRange(edgeDim, ie).first;
    ordinal_type offsetTarget = projStruct->getTargetPointsRange(edgeDim, ie).first;
    for(ordinal_type j=0; j <edgeCardinality; ++j) {
      ordinal_type jdof = cellBasis->getDofOrdinal(edgeDim, ie, j);
      for(ordinal_type ic=0; ic<numCells; ++ic) {
        for(ordinal_type iq=0; iq <numCubPoints; ++iq) {
          for(ordinal_type d=0; d <fieldDim; ++d)
            edgeBasisAtCubPoints(ic,j,iq) += basisAtCubPoints(ic,jdof,offsetBasis+iq,d)*edgeCoeff(d);
          weightedBasisAtElemCubPoints(ic,j,iq) = edgeBasisAtCubPoints(ic,j,iq)*basisEvalWeights(iq);
        }
        for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq) {
          for(ordinal_type d=0; d <fieldDim; ++d)
            weightedBasisAtTargetCubPoints(ic,j,iq) += basisAtTargetCubPoints(ic,jdof,offsetTarget+iq,d)*edgeCoeff(d)*targetEvalWeights(iq);
        }
      }
    }

    for(ordinal_type ic=0; ic<numCells; ++ic)
      for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq)
        for(ordinal_type d=0; d <fieldDim; ++d)
          edgeTargetAtTargetCubPoints(ic,iq) += targetAtEvalPoints(ic,offsetTarget+iq,d)*edgeCoeff(d);

    for(ordinal_type j=0; j <numVertexDofs; ++j) {
      ordinal_type jdof = computedDofs(j);
      for(ordinal_type ic=0; ic<numCells; ++ic)
        for(ordinal_type iq=0; iq <numCubPoints; ++iq)
          for(ordinal_type d=0; d <fieldDim; ++d)
            mComputedProjection(ic,iq) -=  basisCoeffs(ic,jdof)*basisAtCubPoints(ic,jdof,offsetBasis+iq,d)*edgeCoeff(d);
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

  scalarViewType faceCoeff("faceCoeff", fieldDim, faceDofDim);
  for(ordinal_type iface=0; iface<numFaces; ++iface) {


    const auto topoKey = projStruct->getTopologyKey(faceDim,iface);


    ordinal_type faceCardinality = cellBasis->getDofCount(faceDim,iface);

    ordinal_type numTargetCubPoints = projStruct->getNumTargetEvalPoints(faceDim, iface);
    ordinal_type numCubPoints = projStruct->getNumBasisEvalPoints(faceDim, iface);

    if(fieldDim == 1)
      faceCoeff(0,0) = 1;
    else if(isHCurlBAsis) {
      CellTools<SpT>::getReferenceFaceTangents(refFaceTanU, refFaceTanV,iface, cellTopo);
    } else {
      CellTools<SpT>::getReferenceFaceNormal(refFaceNormal, iface, cellTopo);
      for(ordinal_type d=0; d <dim; ++d)
        faceCoeff(d,0) = refFaceNormal(d);
    }

    scalarViewType faceBasisDofAtCubPoints("normaBasisAtCubPoints",numCells,faceCardinality, numCubPoints,faceDofDim);
    scalarViewType wBasisDofAtCubPoints("weightedNormalBasisAtCubPoints",numCells,faceCardinality, numCubPoints,faceDofDim);

    scalarViewType faceBasisAtTargetCubPoints("normalBasisAtTargetCubPoints",numCells,faceCardinality, numTargetCubPoints,faceDofDim);
    scalarViewType wBasisBasisAtTargetCubPoints("weightedNormalBasisAtTargetCubPoints",numCells,faceCardinality, numTargetCubPoints,faceDofDim);

    scalarViewType targetAtTargetCubPoints("targetAtTargetCubPoints",numCells, numTargetCubPoints,faceDofDim);
    scalarViewType mComputedProjection("mNormalComputedProjection", numCells,numCubPoints,faceDofDim);

    ordinal_type offsetBasis = projStruct->getBasisPointsRange(faceDim, iface).first;
    ordinal_type offsetTarget = projStruct->getTargetPointsRange(faceDim, iface).first;
    scalarViewType targetCubWeights = projStruct->getTargetEvalWeights(faceDim, iface);
    scalarViewType CubWeights = projStruct->getBasisEvalWeights(faceDim, iface);

    //Note: we are not considering the jacobian of the orientation map since it is simply a scalar term for the integrals and it does not affect the projection
    for(ordinal_type ic=0; ic<numCells; ++ic)  {

      orts(ic).getFaceOrientation(fOrt.data(), numFaces);

      ordinal_type ort = fOrt(iface);
      if(isHCurlBAsis) {
        Kokkos::deep_copy(faceCoeff,0);
        Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, topoKey, ort);
        for(ordinal_type d=0; d <dim; ++d)
          for(ordinal_type itan=0; itan <faceDim; ++itan)
            for(ordinal_type jtan=0; jtan <faceDim; ++jtan)
              faceCoeff(d,itan) += refFaceTangents(d, jtan)*ortJacobian(jtan,itan);
      }
      for(ordinal_type j=0; j <faceCardinality; ++j) {
        ordinal_type jdof = cellBasis->getDofOrdinal(faceDim, iface, j);
        for(ordinal_type itan=0; itan <faceDofDim; ++itan) {
          for(ordinal_type iq=0; iq <numCubPoints; ++iq) {
            for(ordinal_type d=0; d <fieldDim; ++d)
              faceBasisDofAtCubPoints(ic,j,iq,itan) += faceCoeff(d, itan)*basisAtCubPoints(ic,jdof,offsetBasis+iq,d);
            wBasisDofAtCubPoints(ic,j,iq,itan) = faceBasisDofAtCubPoints(ic,j,iq,itan) * CubWeights(iq);
          }
          for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq) {
            for(ordinal_type d=0; d <fieldDim; ++d)
              faceBasisAtTargetCubPoints(ic,j,iq,itan) += faceCoeff(d, itan)*basisAtTargetCubPoints(ic,jdof,offsetTarget+iq,d);
            wBasisBasisAtTargetCubPoints(ic,j,iq,itan) = faceBasisAtTargetCubPoints(ic,j,iq,itan) * targetCubWeights(iq);
          }
        }
      }

      for(ordinal_type d=0; d <fieldDim; ++d)
        for(ordinal_type itan=0; itan <faceDofDim; ++itan)
          for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq)
            targetAtTargetCubPoints(ic,iq,itan) += faceCoeff(d, itan)*targetAtEvalPoints(ic,offsetTarget+iq,d);

      for(ordinal_type j=0; j <numVertexDofs+numEdgeDofs; ++j) {
        ordinal_type jdof = computedDofs(j);
        for(ordinal_type iq=0; iq <numCubPoints; ++iq)
          for(ordinal_type d=0; d <fieldDim; ++d)
            for(ordinal_type itan=0; itan <faceDofDim; ++itan)
              mComputedProjection(ic,iq,itan) -=  basisCoeffs(ic,jdof)*faceCoeff(d, itan)*basisAtCubPoints(ic,jdof,offsetBasis+iq,d);
      }
    }

    scalarViewType faceMassMat_("faceMassMat_", numCells, faceCardinality, faceCardinality),
        faceRhsMat_("rhsMat_", numCells, faceCardinality);

    FunctionSpaceTools<SpT >::integrate(faceMassMat_, faceBasisDofAtCubPoints, wBasisDofAtCubPoints);
    FunctionSpaceTools<SpT >::integrate(faceRhsMat_, targetAtTargetCubPoints, wBasisBasisAtTargetCubPoints);
    FunctionSpaceTools<SpT >::integrate(faceRhsMat_, mComputedProjection, wBasisDofAtCubPoints,true);

    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> faceMassMat("faceMassMat", faceCardinality,faceCardinality);
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> faceRhsMat("faceRhsMat",faceCardinality, 1);

    Teuchos::LAPACK<ordinal_type,funValsValueType> lapack;
    ordinal_type info = 0;
    for(ordinal_type ic=0; ic<numCells; ++ic)  {

      for(ordinal_type i=0; i<faceCardinality; ++i) {
        faceRhsMat(i,0) = faceRhsMat_(ic,i);
        for(ordinal_type j=0; j<faceCardinality; ++j){
          faceMassMat(i,j) = faceMassMat_(ic,i,j);
        }
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

    range_type cellPointsRange = projStruct->getTargetPointsRange(dim, 0);

    ordinal_type numTargetCubPoints = projStruct->getNumTargetEvalPoints(dim,0);
    ordinal_type numCubPoints = projStruct->getNumBasisEvalPoints(dim,0);

    scalarViewType internalBasisAtCubPoints("internalBasisAtCubPoints",numCells,numElemDofs, numCubPoints, fieldDim);
    scalarViewType mComputedProjection("mComputedProjection", numCells, numCubPoints, fieldDim);

    scalarViewType targetCubWeights = projStruct->getTargetEvalWeights(dim, 0);
    scalarViewType cubWeights = projStruct->getBasisEvalWeights(dim, 0);
    ordinal_type offsetBasis = projStruct->getBasisPointsRange(dim, 0).first;
    ordinal_type offsetTarget = projStruct->getTargetPointsRange(dim, 0).first;


    scalarViewType wBasisAtCubPoints("weightedBasisAtCubPoints",numCells,numElemDofs, numCubPoints,fieldDim);
    scalarViewType wBasisBasisAtTargetCubPoints("weightedBasisAtTargetCubPoints",numCells,numElemDofs, numTargetCubPoints,fieldDim);
    for(ordinal_type j=0; j <numElemDofs; ++j) {
      ordinal_type idof = cellBasis->getDofOrdinal(dim, 0, j);
      for(ordinal_type ic=0; ic<numCells; ++ic) {
        for(ordinal_type d=0; d <fieldDim; ++d) {
          for(ordinal_type iq=0; iq <numCubPoints; ++iq) {
            internalBasisAtCubPoints(ic,j,iq,d) = basisAtCubPoints(ic,idof,offsetBasis+iq,d);
            wBasisAtCubPoints(ic,j,iq,d) = internalBasisAtCubPoints(ic,j,iq,d) * cubWeights(iq);
          }
          for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq) {
            wBasisBasisAtTargetCubPoints(ic,j,iq,d) = basisAtTargetCubPoints(ic,idof,offsetTarget+iq,d)* targetCubWeights(iq);
          }
        }
      }
    }

    for(ordinal_type j=0; j <numVertexDofs+numEdgeDofs+numFaceDofs; ++j) {
      ordinal_type jdof = computedDofs(j);
      for(ordinal_type iq=0; iq <numCubPoints; ++iq)
        for(ordinal_type ic=0; ic<numCells; ++ic)  {
          for(ordinal_type d=0; d <fieldDim; ++d) {
            mComputedProjection(ic,iq,d) -=  basisCoeffs(ic,jdof)*basisAtCubPoints(ic,jdof,offsetBasis+iq,d);
          }
        }
    }


    scalarViewType cellMassMat_("cellMassMat_", numCells, numElemDofs, numElemDofs),
        cellRhsMat_("rhsMat_", numCells, numElemDofs);

    FunctionSpaceTools<SpT >::integrate(cellMassMat_, internalBasisAtCubPoints, wBasisAtCubPoints);
    if(fieldDim==1)
      FunctionSpaceTools<SpT >::integrate(cellRhsMat_, Kokkos::subview(targetAtEvalPoints,Kokkos::ALL(),cellPointsRange,Kokkos::ALL()),
          Kokkos::subview(wBasisBasisAtTargetCubPoints,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),0));
    else
      FunctionSpaceTools<SpT >::integrate(cellRhsMat_, Kokkos::subview(targetAtEvalPoints,Kokkos::ALL(),cellPointsRange,Kokkos::ALL()), wBasisBasisAtTargetCubPoints);
    FunctionSpaceTools<SpT >::integrate(cellRhsMat_, mComputedProjection, wBasisAtCubPoints, true);

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

