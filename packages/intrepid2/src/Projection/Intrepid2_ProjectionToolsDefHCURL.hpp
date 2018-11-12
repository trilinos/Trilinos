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

template<typename SpT>
template<typename BasisType,
typename ortValueType,       class ...ortProperties>
void
ProjectionTools<SpT>::getHCurlEvaluationPoints(typename BasisType::scalarViewType evaluationPoints,
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

  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;

  Kokkos::View<ordinal_type*> eOrt("eOrt", numEdges), fOrt("fOrt", numFaces);

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

    scalarViewType cubPoints;//("cubPoints", numTargetCubPoints, faceDim);
    range_type facePointsRange;
    if(evalPointType == TARGET) {
      cubPoints = projStruct->getTargetEvalPoints(faceDim, iface);
      facePointsRange = projStruct->getTargetPointsRange(faceDim, iface);
    } else {
      cubPoints = projStruct->getBasisEvalPoints(faceDim, iface);
      facePointsRange = projStruct->getBasisPointsRange(faceDim, iface);
    }

    scalarViewType curlCubPoints;//("curlCubPoints", numTargetCurlCubPoints, faceDim);
    range_type faceCurlPointsRange;
    if(evalPointType == TARGET) {
      curlCubPoints = projStruct->getTargetDerivEvalPoints(faceDim, iface);
      faceCurlPointsRange = projStruct->getTargetDerivPointsRange(faceDim, iface);
    } else {
      curlCubPoints = projStruct->getBasisDerivEvalPoints(faceDim, iface);
      faceCurlPointsRange = projStruct->getBasisDerivPointsRange(faceDim, iface);
    }

    scalarViewType faceCubPoints("faceCubPoints", cubPoints.extent(0), faceDim);
    scalarViewType faceCurlCubPoints("faceCurlCubPoints", curlCubPoints.extent(0), faceDim);

    const auto topoKey = projStruct->getTopologyKey(faceDim,iface);
    for(ordinal_type ic=0; ic<numCells; ++ic) {

      orts(ic).getFaceOrientation(fOrt.data(), numFaces);

      ordinal_type ort = fOrt(iface);
      Impl::OrientationTools::mapToModifiedReference(faceCubPoints,cubPoints,topoKey,ort);
      CellTools<SpT>::mapToReferenceSubcell(Kokkos::subview(evaluationPoints, ic, facePointsRange, Kokkos::ALL()), faceCubPoints, faceDim, iface, cellBasis->getBaseCellTopology());

      Impl::OrientationTools::mapToModifiedReference(faceCurlCubPoints,curlCubPoints,topoKey,ort);
      CellTools<SpT>::mapToReferenceSubcell(Kokkos::subview(extDerivEvaluationPoints,  ic, faceCurlPointsRange, Kokkos::ALL()), faceCurlCubPoints, faceDim, iface, cellBasis->getBaseCellTopology());
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

    range_type cellCurlPointsRange;
    scalarViewType curlCubPoints;
    if(evalPointType == TARGET) {
      curlCubPoints = projStruct->getTargetDerivEvalPoints(dim, 0);
      cellCurlPointsRange = projStruct->getTargetDerivPointsRange(dim, 0);
    } else {
      curlCubPoints = projStruct->getBasisDerivEvalPoints(dim, 0);
      cellCurlPointsRange = projStruct->getBasisDerivPointsRange(dim, 0);
    }
    RealSpaceTools<SpT>::clone(Kokkos::subview(extDerivEvaluationPoints, Kokkos::ALL(), cellCurlPointsRange, Kokkos::ALL()), curlCubPoints);
  }
}


template<typename SpT>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<SpT>::getHCurlBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtEvalPoints,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetCurlAtCurlEvalPoints,
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
      numTotalCurlEvaluationPoints(targetCurlAtCurlEvalPoints.extent(1));
  ordinal_type basisCardinality = cellBasis->getCardinality();
  ordinal_type numCells = targetAtEvalPoints.extent(0);
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;
  const ordinal_type derDim = dim == 3 ? dim : 1;

  const std::string& name = cellBasis->getName();

  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;

  Kokkos::View<ordinal_type*> eOrt("eOrt", numEdges);
  Kokkos::View<ordinal_type*> fOrt("fOrt", numFaces);
  scalarViewType refEdgeTan("refEdgeTan",  dim);
  scalarViewType refFaceTangents("refFaceTangents", dim, 2);
  scalarViewType refFaceNormal("refFaceNormal", dim);
  auto refFaceTanU = Kokkos::subview(refFaceTangents, Kokkos::ALL, 0);
  auto refFaceTanV = Kokkos::subview(refFaceTangents, Kokkos::ALL, 1);

  ordinal_type numEdgeDofs(0);
  for(ordinal_type ie=0; ie<numEdges; ++ie)
    numEdgeDofs += cellBasis->getDofCount(edgeDim,ie);

  ordinal_type numFaceDofs(0);
  for(ordinal_type iface=0; iface<numFaces; ++iface)
    numFaceDofs += cellBasis->getDofCount(faceDim,iface);

  Kokkos::View<ordinal_type*> computedDofs("computedDofs",numEdgeDofs+numFaceDofs);

  ordinal_type computedDofsCount = 0;

  ordinal_type numTotalCubPoints = projStruct->getNumBasisEvalPoints(), numTotalCurlCubPoints = projStruct->getNumBasisDerivEvalPoints();
  scalarViewType cubPoints("cubPoints",numCells,numTotalCubPoints, dim);
  scalarViewType curlCubPoints("curlCubPoints",numCells,numTotalCurlCubPoints, dim);
  getHCurlEvaluationPoints(cubPoints, curlCubPoints, orts, cellBasis, projStruct, BASIS);

  scalarViewType basisAtCubPoints("basisAtCubPoints",numCells,basisCardinality, numTotalCubPoints, dim);
  scalarViewType basisAtTargetCubPoints("basisAtTargetCubPoints",numCells,basisCardinality, numTotalEvaluationPoints, dim);
  {
    scalarViewType nonOrientedBasisAtCubPoints("nonOrientedBasisAtCubPoints",numCells,basisCardinality, numTotalCubPoints, dim);
    scalarViewType nonOrientedBasisAtTargetCubPoints("nonOrientedBasisAtTargetCubPoints",numCells,basisCardinality, numTotalEvaluationPoints, dim);
    for(ordinal_type ic=0; ic<numCells; ++ic) {
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtTargetCubPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(evaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtCubPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(cubPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
    }

    OrientationTools<SpT>::modifyBasisByOrientation(basisAtCubPoints, nonOrientedBasisAtCubPoints, orts, cellBasis);
    OrientationTools<SpT>::modifyBasisByOrientation(basisAtTargetCubPoints, nonOrientedBasisAtTargetCubPoints, orts, cellBasis);
  }

  scalarViewType basisCurlAtCurlCubPoints;
  scalarViewType basisCurlAtTargetCurlCubPoints;
  if(numTotalCurlEvaluationPoints>0) {
    scalarViewType nonOrientedBasisCurlAtTargetCurlCubPoints, nonOrientedBasisCurlAtCurlCubPoints;
    if (dim == 3) {
      basisCurlAtCurlCubPoints = scalarViewType ("basisCurlAtCurlCubPoints",numCells,basisCardinality, numTotalCurlCubPoints, dim);
      nonOrientedBasisCurlAtCurlCubPoints = scalarViewType ("nonOrientedBasisCurlAtCurlCubPoints",numCells,basisCardinality, numTotalCurlCubPoints, dim);
      basisCurlAtTargetCurlCubPoints = scalarViewType("basisCurlAtTargetCurlCubPoints",numCells,basisCardinality, numTotalCurlEvaluationPoints, dim);
      nonOrientedBasisCurlAtTargetCurlCubPoints = scalarViewType("nonOrientedBasisCurlAtTargetCurlCubPoints",numCells,basisCardinality, numTotalCurlEvaluationPoints, dim);
    } else {
      basisCurlAtCurlCubPoints = scalarViewType ("basisCurlAtCurlCubPoints",numCells,basisCardinality, numTotalCurlCubPoints);
      nonOrientedBasisCurlAtCurlCubPoints = scalarViewType ("nonOrientedBasisCurlAtCurlCubPoints",numCells,basisCardinality, numTotalCurlCubPoints);
      basisCurlAtTargetCurlCubPoints = scalarViewType("basisCurlAtTargetCurlCubPoints",numCells,basisCardinality, numTotalCurlEvaluationPoints);
      nonOrientedBasisCurlAtTargetCurlCubPoints = scalarViewType("nonOrientedBasisCurlAtTargetCurlCubPoints",numCells,basisCardinality, numTotalCurlEvaluationPoints);
    }
    for(ordinal_type ic=0; ic<numCells; ++ic) {
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisCurlAtCurlCubPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(curlCubPoints, ic, Kokkos::ALL(), Kokkos::ALL()),OPERATOR_CURL);
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisCurlAtTargetCurlCubPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(extDerivEvaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()),OPERATOR_CURL);
    }
    OrientationTools<SpT>::modifyBasisByOrientation(basisCurlAtCurlCubPoints, nonOrientedBasisCurlAtCurlCubPoints, orts, cellBasis);
    OrientationTools<SpT>::modifyBasisByOrientation(basisCurlAtTargetCurlCubPoints, nonOrientedBasisCurlAtTargetCurlCubPoints, orts, cellBasis);
  }

  for(ordinal_type ie=0; ie<numEdges; ++ie)  {

    ordinal_type edgeCardinality = cellBasis->getDofCount(edgeDim,ie);
    ordinal_type numCubPoints = projStruct->getNumBasisEvalPoints(edgeDim, ie);
    ordinal_type numTargetCubPoints = projStruct->getNumTargetEvalPoints(edgeDim, ie);

    CellTools<SpT>::getReferenceEdgeTangent(refEdgeTan, ie, cellBasis->getBaseCellTopology());

    scalarViewType tanBasisAtElemCubPoints("tanBasisAtElemCubPoints",numCells,edgeCardinality, numCubPoints);
    scalarViewType tanBasisAtTargetCubPoints("tanBasisAtTargetCubPoints",numCells,edgeCardinality, numTargetCubPoints);
    scalarViewType weightedTanBasisAtElemCubPoints("weightedTanBasisAtElemCubPoints",numCells,edgeCardinality, numCubPoints);
    scalarViewType weightedTanBasisAtTargetCubPoints("weightedTanBasisAtTargetCubPoints",numCells,edgeCardinality, numTargetCubPoints);
    scalarViewType tanTargetAtTargetCubPoints("normalTargetAtTargetCubPoints",numCells, numTargetCubPoints);

    scalarViewType targetEvalWeights = projStruct->getTargetEvalWeights(edgeDim, ie);
    scalarViewType basisEvalWeights = projStruct->getBasisEvalWeights(edgeDim, ie);

    //Note: we are not considering the jacobian of the orientation map since it is simply a scalar term for the integrals and it does not affect the projection
    ordinal_type offsetBasis = projStruct->getBasisPointsRange(edgeDim, ie).first;
    ordinal_type offsetTarget = projStruct->getTargetPointsRange(edgeDim, ie).first;
    for(ordinal_type ic=0; ic<numCells; ++ic)  {
      for(ordinal_type j=0; j <edgeCardinality; ++j) {
        ordinal_type jdof = cellBasis->getDofOrdinal(edgeDim, ie, j);
        for(ordinal_type iq=0; iq <numCubPoints; ++iq) {
          for(ordinal_type d=0; d <dim; ++d)
            tanBasisAtElemCubPoints(ic,j,iq) += refEdgeTan(d)*basisAtCubPoints(ic,jdof,offsetBasis+iq,d);
          weightedTanBasisAtElemCubPoints(ic,j,iq) = tanBasisAtElemCubPoints(ic,j,iq)*basisEvalWeights(iq);
        }
        for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq) {
          for(ordinal_type d=0; d <dim; ++d)
            tanBasisAtTargetCubPoints(ic,j,iq) += refEdgeTan(d)*basisAtTargetCubPoints(ic,jdof,offsetTarget+iq,d);
          weightedTanBasisAtTargetCubPoints(ic,j,iq) = tanBasisAtTargetCubPoints(ic,j,iq)*targetEvalWeights(iq);
        }
      }
      for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq)
        for(ordinal_type d=0; d <dim; ++d)
          tanTargetAtTargetCubPoints(ic,iq) += refEdgeTan(d)*targetAtEvalPoints(ic,offsetTarget+iq,d);
    }
    scalarViewType edgeMassMat_("edgeMassMat_", numCells, edgeCardinality+1, edgeCardinality+1),
        edgeRhsMat_("rhsMat_", numCells, edgeCardinality+1);

    scalarViewType cubWeights_("cubWeights_", numCells, 1, basisEvalWeights.extent(0)), targetEvalWeights_("targetEvalWeights", numCells, 1, targetEvalWeights.extent(0));
    RealSpaceTools<SpT>::clone(cubWeights_, basisEvalWeights);
    RealSpaceTools<SpT>::clone(targetEvalWeights_, targetEvalWeights);

    range_type range_H(0, edgeCardinality);
    range_type range_B(edgeCardinality, edgeCardinality+1);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(edgeMassMat_,Kokkos::ALL(),range_H,range_H), tanBasisAtElemCubPoints, weightedTanBasisAtElemCubPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(edgeMassMat_,Kokkos::ALL(),range_H,range_B), tanBasisAtElemCubPoints, cubWeights_);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(edgeRhsMat_,Kokkos::ALL(),range_H), tanTargetAtTargetCubPoints, weightedTanBasisAtTargetCubPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(edgeRhsMat_,Kokkos::ALL(),range_B), tanTargetAtTargetCubPoints, targetEvalWeights_);
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> edgeMassMat("edgeMassMat", edgeCardinality+1,edgeCardinality+1);
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> edgeRhsMat("edgeRhsMat",edgeCardinality+1, 1);

    Teuchos::LAPACK<ordinal_type,funValsValueType> lapack;
    ordinal_type info = 0;
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> pivVec("pivVec", edgeCardinality+1, 1);
    for(ordinal_type ic=0; ic<numCells; ++ic)  {
      Kokkos::deep_copy(edgeMassMat,funValsValueType(0));  //LAPACK might overwrite the matrix

      for(ordinal_type i=0; i<edgeCardinality; ++i) {
        edgeRhsMat(i,0) = edgeRhsMat_(ic,i);
        for(ordinal_type j=0; j<edgeCardinality+1; ++j)
          edgeMassMat(i,j) = edgeMassMat_(ic,i,j);
        edgeMassMat(edgeCardinality,i) = edgeMassMat_(ic,i,edgeCardinality);
      }
      edgeRhsMat(edgeCardinality,0) = edgeRhsMat_(ic,edgeCardinality);

      lapack.GESV(edgeCardinality+1, 1,
          edgeMassMat.data(),
          edgeMassMat.stride_1(),
          (ordinal_type*)pivVec.data(),
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

  Basis<host_space_type,scalarType,scalarType> *hgradBasis = NULL;
  for(ordinal_type iface=0; iface<numFaces; ++iface) {

    if(name.find("HEX")!=std::string::npos)
      hgradBasis = new Basis_HGRAD_QUAD_Cn_FEM<host_space_type,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else if(name.find("TET")!=std::string::npos)
      hgradBasis = new Basis_HGRAD_TRI_Cn_FEM<host_space_type,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else  {
      std::stringstream ss;
      ss << ">>> ERROR (Intrepid2::ProjectionTools::getHCurlBasisCoeffs): "
          << "Method not implemented for basis " << name;
      INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
      return;
    }


    ordinal_type numFaceDofs = cellBasis->getDofCount(faceDim,iface);
    ordinal_type numTargetCubPoints = projStruct->getNumTargetEvalPoints(faceDim, iface);

    ordinal_type numTargetCurlCubPoints = projStruct->getNumTargetDerivEvalPoints(faceDim, iface);
    range_type faceCurlPointsRange = projStruct->getTargetDerivPointsRange(faceDim, iface);

    ordinal_type numCubPoints = projStruct->getNumBasisEvalPoints(faceDim, iface);

    scalarViewType hgradBasisGradAtCubPoints("hgradBasisGradAtCubPoints",hgradBasis->getCardinality(), numCubPoints, faceDim);
    scalarViewType hgradBasisGradAtTargetCubPoints("hgradBasisGradAtTargetCubPoints",hgradBasis->getCardinality(), numTargetCubPoints, faceDim);

    ordinal_type internalHgradCardinality = hgradBasis->getDofCount(faceDim,0);
    scalarViewType internalHgradBasisGradAtCubPoints("internalHgradBasisGradAtCubPoints",1, internalHgradCardinality, numCubPoints, faceDim);
    scalarViewType internalHgradBasisGradAtTargetCubPoints("internalHgradBasisGradAtTargetCubPoints",1, internalHgradCardinality, numTargetCubPoints, faceDim);


    CellTools<SpT>::getReferenceFaceNormal(refFaceNormal, iface, cellTopo);
    CellTools<SpT>::getReferenceFaceTangents(refFaceTanU, refFaceTanV,iface, cellTopo);

    hgradBasis->getValues(hgradBasisGradAtCubPoints,projStruct->getBasisEvalPoints(faceDim, iface), OPERATOR_GRAD);
    hgradBasis->getValues(hgradBasisGradAtTargetCubPoints,projStruct->getTargetEvalPoints(faceDim, iface),OPERATOR_GRAD);

    for(ordinal_type j=0; j <internalHgradCardinality; ++j) {
      ordinal_type face_dof = hgradBasis->getDofOrdinal(faceDim, 0, j);
      for(ordinal_type d=0; d <faceDim; ++d) {
        for(ordinal_type iq=0; iq <numCubPoints; ++iq)
          internalHgradBasisGradAtCubPoints(0,j,iq,d)= hgradBasisGradAtCubPoints(face_dof,iq,d);
        for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq)
          internalHgradBasisGradAtTargetCubPoints(0,j,iq,d)= hgradBasisGradAtTargetCubPoints(face_dof,iq,d);
      }
    }

    scalarViewType tanBasisAtElemCubPoints("tanBasisAtElemCubPoints",numCells,numFaceDofs, numCubPoints,dim-1);
    scalarViewType tanBasisAtTargetCubPoints("tanBasisAtTargetCubPoints",numCells,numFaceDofs, numTargetCubPoints,dim-1);
    scalarViewType normalBasisCurlAtElemCubPoints("normaBasisCurlAtElemCubPoints",numCells,numFaceDofs, numCubPoints);
    scalarViewType wNormalBasisCurlAtElemCubPoints("weightedNormalBasisCurlAtElemCubPoints",numCells,numFaceDofs, numCubPoints);

    scalarViewType tanTargetAtTargetCubPoints("tanTargetAtTargetCubPoints",numCells, numTargetCubPoints, dim-1);
    scalarViewType normalTargetCurlAtTargetCubPoints("normalTargetCurlAtTargetCubPoints",numCells, numTargetCurlCubPoints);
    scalarViewType normalBasisCurlAtTargetCurlCubPoints("normalBasisCurlAtTargetCurlCubPoints",numCells,numFaceDofs, numTargetCurlCubPoints);
    scalarViewType wNormalBasisCurlBasisAtTargetCurlCubPoints("weightedNormalBasisCurlAtTargetCurlCubPoints",numCells,numFaceDofs, numTargetCurlCubPoints);

    scalarViewType wHgradBasisGradAtCubPoints("wHgradBasisGradAtCubPoints",1, internalHgradCardinality, numCubPoints, faceDim);
    scalarViewType wHgradBasisGradAtCubPoints_("wHgradBasisGradAtCubPoints_",numCells, internalHgradCardinality, numCubPoints, faceDim);
    scalarViewType wHgradBasisGradAtTargetCubPoints("wHgradBasisGradAtTargetCubPoints",1, internalHgradCardinality, numTargetCubPoints, faceDim);
    scalarViewType wHgradBasisGradAtTargetCubPoints_("wHgradBasisGradAtTargetCubPoints_",numCells, internalHgradCardinality, numTargetCubPoints, faceDim);

    scalarViewType mNormalComputedProjectionCurl("mNormalComputedProjection", numCells,numCubPoints);
    scalarViewType mTanComputedProjection("mTanComputedProjection", numCells,numCubPoints,dim-1);

    scalarViewType targetDerivEvalWeights = projStruct->getTargetDerivEvalWeights(faceDim, iface);
    ordinal_type offsetBasis = projStruct->getBasisPointsRange(faceDim, iface).first;
    ordinal_type offsetBasisCurl = projStruct->getBasisDerivPointsRange(faceDim, iface).first;
    ordinal_type offsetTarget = projStruct->getTargetPointsRange(faceDim, iface).first;
    ordinal_type offsetTargetCurl = projStruct->getTargetDerivPointsRange(faceDim, iface).first;


    //Note: we are not considering the jacobian of the orientation map since it is simply a scalar term for the integrals and it does not affect the projection
    const auto topoKey = projStruct->getTopologyKey(faceDim,iface);
    for(ordinal_type ic=0; ic<numCells; ++ic)  {
      orts(ic).getFaceOrientation(fOrt.data(), numFaces);

      ordinal_type ort = fOrt(iface);

      Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, topoKey, ort);

      for(ordinal_type j=0; j <numFaceDofs; ++j) {
        ordinal_type jdof = cellBasis->getDofOrdinal(faceDim, iface, j);
        for(ordinal_type iq=0; iq <numCubPoints; ++iq)
          for(ordinal_type d=0; d <dim; ++d) {
            normalBasisCurlAtElemCubPoints(ic,j,iq) += refFaceNormal(d)*basisCurlAtCurlCubPoints(ic,jdof,offsetBasisCurl+iq,d);
            for(ordinal_type itan=0; itan <dim-1; ++itan)
              for(ordinal_type jtan=0; jtan <dim-1; ++jtan)
                tanBasisAtElemCubPoints(ic,j,iq,itan) += refFaceTangents(d,jtan)*ortJacobian(jtan,itan)*basisAtCubPoints(ic,jdof,offsetBasis+iq,d);
          }
        for(ordinal_type iq=0; iq <numTargetCurlCubPoints; ++iq) {
          for(ordinal_type d=0; d <dim; ++d)
            normalBasisCurlAtTargetCurlCubPoints(ic,j,iq) += refFaceNormal(d)*basisCurlAtTargetCurlCubPoints(ic,jdof,offsetTargetCurl+iq,d);
          wNormalBasisCurlBasisAtTargetCurlCubPoints(ic,j,iq) = normalBasisCurlAtTargetCurlCubPoints(ic,j,iq)*targetDerivEvalWeights[iq];
        }
      }
      for(ordinal_type j=0; j <numEdgeDofs; ++j) {
        ordinal_type jdof = computedDofs(j);
        for(ordinal_type iq=0; iq <numCubPoints; ++iq)
          for(ordinal_type d=0; d <dim; ++d) {
            mNormalComputedProjectionCurl(ic,iq) -=  refFaceNormal(d)*basisCoeffs(ic,jdof)*basisCurlAtCurlCubPoints(ic,jdof,offsetBasisCurl+iq,d);
            for(ordinal_type itan=0; itan <dim-1; ++itan)
              for(ordinal_type jtan=0; jtan <dim-1; ++jtan)
                mTanComputedProjection(ic,iq,itan) -=  refFaceTangents(d,jtan)*ortJacobian(jtan,itan)*basisCoeffs(ic,jdof)*basisAtCubPoints(ic,jdof,offsetBasis+iq,d);
          }
      }
      for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq)
        for(ordinal_type d=0; d <dim; ++d)
          for(ordinal_type itan=0; itan <dim-1; ++itan)
            for(ordinal_type jtan=0; jtan <dim-1; ++jtan)
              tanTargetAtTargetCubPoints(ic,iq,itan) += refFaceTangents(d,jtan)*ortJacobian(jtan,itan)*targetAtEvalPoints(ic,offsetTarget+iq,d);
      for(ordinal_type iq=0; iq <numTargetCurlCubPoints; ++iq)
        for(ordinal_type d=0; d <dim; ++d)
          normalTargetCurlAtTargetCubPoints(ic,iq) += refFaceNormal(d)*targetCurlAtCurlEvalPoints(ic,offsetTargetCurl+iq,d);
    }

    scalarViewType faceMassMat_("faceMassMat_", numCells, numFaceDofs+internalHgradCardinality, numFaceDofs+internalHgradCardinality),
        faceRhsMat_("rhsMat_", numCells, numFaceDofs+internalHgradCardinality);

    scalarViewType targetCubWeights_("targetCubWeights_", 1, projStruct->getNumTargetEvalPoints(faceDim, iface));
    RealSpaceTools<SpT>::clone(targetCubWeights_, projStruct->getTargetEvalWeights(faceDim, iface));
    scalarViewType cubWeights_("cubWeights_", numCells, 1, numCubPoints);
    RealSpaceTools<SpT>::clone(cubWeights_, projStruct->getBasisEvalWeights(faceDim, iface));
    ArrayTools<SpT>::scalarMultiplyDataField( wNormalBasisCurlAtElemCubPoints, Kokkos::subview(cubWeights_, Kokkos::ALL(),0, Kokkos::ALL()),normalBasisCurlAtElemCubPoints, false);
    ArrayTools<SpT>::scalarMultiplyDataField( wHgradBasisGradAtCubPoints, Kokkos::subview(cubWeights_, 0, Kokkos::ALL(), Kokkos::ALL()),internalHgradBasisGradAtCubPoints, false);
    ArrayTools<SpT>::scalarMultiplyDataField( wHgradBasisGradAtTargetCubPoints, targetCubWeights_, internalHgradBasisGradAtTargetCubPoints , false);

    RealSpaceTools<SpT>::clone(wHgradBasisGradAtCubPoints_,Kokkos::subview(wHgradBasisGradAtCubPoints,0,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()));
    RealSpaceTools<SpT>::clone(wHgradBasisGradAtTargetCubPoints_,Kokkos::subview(wHgradBasisGradAtTargetCubPoints,0,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()));

    range_type range_H(0, numFaceDofs);
    range_type range_B(numFaceDofs, numFaceDofs+internalHgradCardinality);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(faceMassMat_,Kokkos::ALL(),range_H,range_H), normalBasisCurlAtElemCubPoints, wNormalBasisCurlAtElemCubPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(faceMassMat_,Kokkos::ALL(),range_H,range_B), tanBasisAtElemCubPoints, wHgradBasisGradAtCubPoints_);

    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(faceRhsMat_,Kokkos::ALL(),range_H), normalTargetCurlAtTargetCubPoints, wNormalBasisCurlBasisAtTargetCurlCubPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(faceRhsMat_,Kokkos::ALL(),range_H), mNormalComputedProjectionCurl, wNormalBasisCurlAtElemCubPoints,true);

    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(faceRhsMat_,Kokkos::ALL(),range_B), tanTargetAtTargetCubPoints, wHgradBasisGradAtTargetCubPoints_);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(faceRhsMat_,Kokkos::ALL(),range_B), mTanComputedProjection, wHgradBasisGradAtCubPoints_,true);

    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> faceMassMat("faceMassMat", numFaceDofs+internalHgradCardinality,numFaceDofs+internalHgradCardinality);
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> faceRhsMat("faceRhsMat",numFaceDofs+internalHgradCardinality, 1);

    Teuchos::LAPACK<ordinal_type,funValsValueType> lapack;
    ordinal_type info = 0;
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> pivVec("pivVec", numFaceDofs+internalHgradCardinality, 1);
    for(ordinal_type ic=0; ic<numCells; ++ic)  {
      Kokkos::deep_copy(faceMassMat,funValsValueType(0));  //LAPACK might overwrite the matrix

      for(ordinal_type i=0; i<numFaceDofs; ++i) {
        for(ordinal_type j=0; j<numFaceDofs+internalHgradCardinality; ++j)
          faceMassMat(i,j) = faceMassMat_(ic,i,j);
        for(ordinal_type j=0; j<internalHgradCardinality; ++j)
          faceMassMat(numFaceDofs+j,i) = faceMassMat_(ic,i,numFaceDofs+j);
      }

      for(ordinal_type i=0; i<numFaceDofs+internalHgradCardinality; ++i)
        faceRhsMat(i,0) = faceRhsMat_(ic,i);

      lapack.GESV(numFaceDofs+internalHgradCardinality, 1,
          faceMassMat.data(),
          faceMassMat.stride_1(),
          (ordinal_type*)pivVec.data(),
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

      for(ordinal_type i=0; i<numFaceDofs; ++i){
        ordinal_type face_dof = cellBasis->getDofOrdinal(faceDim, iface, i);
        basisCoeffs(ic,face_dof) = faceRhsMat(i,0);
      }
    }

    for(ordinal_type i=0; i<numFaceDofs; ++i)
      computedDofs(computedDofsCount++) = cellBasis->getDofOrdinal(faceDim, iface, i);

    delete hgradBasis;
  }

  ordinal_type numElemDofs = cellBasis->getDofCount(dim,0);
  if(numElemDofs>0) {
    if(name.find("HEX")!=std::string::npos)
      hgradBasis = new Basis_HGRAD_HEX_Cn_FEM<host_space_type,scalarType,scalarType>(cellBasis->getDegree());
    else if(name.find("TET")!=std::string::npos)
      hgradBasis = new Basis_HGRAD_TET_Cn_FEM<host_space_type,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else if(name.find("TRI")!=std::string::npos)
      hgradBasis = new Basis_HGRAD_TRI_Cn_FEM<host_space_type,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else if(name.find("QUAD")!=std::string::npos)
      hgradBasis = new Basis_HGRAD_QUAD_Cn_FEM<host_space_type,scalarType,scalarType>(cellBasis->getDegree(),POINTTYPE_WARPBLEND);
    else  {
      std::stringstream ss;
      ss << ">>> ERROR (Intrepid2::ProjectionTools::getHCurlBasisCoeffs): "
          << "Method not implemented for basis " << name;
      INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
      return;
    }

    range_type cellPointsRange = projStruct->getTargetPointsRange(dim, 0);
    range_type cellCurlPointsRange = projStruct->getTargetDerivPointsRange(dim, 0);

    ordinal_type numTargetCurlCubPoints = projStruct->getNumTargetDerivEvalPoints(dim,0);
    ordinal_type numCubPoints = projStruct->getNumBasisEvalPoints(dim,0);
    ordinal_type numTargetCubPoints = projStruct->getNumTargetEvalPoints(dim,0);

    scalarViewType hgradBasisGradAtCubPoints("hgradBasisGradAtCubPoints",hgradBasis->getCardinality(), numCubPoints, dim);
    scalarViewType hgradBasisGradAtTargetCubPoints("hgradBasisGradAtTargetCubPoints",hgradBasis->getCardinality(), numTargetCubPoints, dim);

    ordinal_type internalHgradCardinality = hgradBasis->getDofCount(dim,0);
    scalarViewType internalHgradBasisGradAtCubPoints("internalHgradBasisGradAtCubPoints",1, internalHgradCardinality, numCubPoints, dim);
    scalarViewType internalHgradBasisGradAtTargetCubPoints("internalHgradBasisGradAtTargetCubPoints",numCells, internalHgradCardinality, numTargetCubPoints, dim);
    scalarViewType wHgradBasisGradAtTargetCubPoints("wHgradBasisGradAtTargetCubPoints",numCells, internalHgradCardinality, numTargetCubPoints, dim);
    scalarViewType wHgradBasisGradAtCubPoints("wHgradBasisGradAtCubPoints",numCells, internalHgradCardinality, numCubPoints, dim);

    scalarViewType targetEvalWeights = projStruct->getTargetEvalWeights(dim, 0);
    scalarViewType basisEvalWeights = projStruct->getBasisEvalWeights(dim, 0);

    hgradBasis->getValues(hgradBasisGradAtCubPoints,projStruct->getBasisEvalPoints(dim, 0), OPERATOR_GRAD);
    hgradBasis->getValues(hgradBasisGradAtTargetCubPoints,projStruct->getTargetEvalPoints(dim, 0),OPERATOR_GRAD);

    for(ordinal_type j=0; j <internalHgradCardinality; ++j) {
      ordinal_type idof = hgradBasis->getDofOrdinal(dim, 0, j);
      for(ordinal_type d=0; d <dim; ++d) {
        for(ordinal_type iq=0; iq <numCubPoints; ++iq) {
          internalHgradBasisGradAtCubPoints(0,j,iq,d)= hgradBasisGradAtCubPoints(idof,iq,d);
          for(ordinal_type ic=0; ic<numCells; ++ic)
            wHgradBasisGradAtCubPoints(ic,j,iq,d) = internalHgradBasisGradAtCubPoints(0,j,iq,d)*basisEvalWeights[iq];
        }
        for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq)
          for(ordinal_type ic=0; ic<numCells; ++ic) {
            internalHgradBasisGradAtTargetCubPoints(ic,j,iq,d)= hgradBasisGradAtTargetCubPoints(idof,iq,d);
            wHgradBasisGradAtTargetCubPoints(ic,j,iq,d) = internalHgradBasisGradAtTargetCubPoints(ic,j,iq,d)*targetEvalWeights[iq];
          }
      }
    }

    scalarViewType internalBasisAtElemcubPoints("basisElemAtcubPoints",numCells,numElemDofs, numCubPoints, dim);
    scalarViewType internalBasisCurlAtElemcubPoints("internalBasisCurlAtElemcubPoints",numCells,numElemDofs, numCubPoints, derDim);
    scalarViewType internalBasisCurlAtTargetCurlCubPoints("weightedBasisCurlAtElemCubPoints",numCells,numElemDofs, numTargetCurlCubPoints,derDim);
    scalarViewType mComputedProjectionCurl("mComputedProjectionCurl", numCells, numCubPoints, derDim);
    scalarViewType mComputedProjection("mComputedProjection", numCells, numCubPoints, dim);
    scalarViewType wBasisCurlAtElemCubPoints("weightedBasisCurlAtElemCubPoints",numCells,numElemDofs, numCubPoints,derDim);
    scalarViewType wBasisCurlBasisAtTargetCurlCubPoints("weightedBasisCurlAtTargetCurlCubPoints",numCells,numElemDofs, numTargetCurlCubPoints,derDim);
    scalarViewType targetDerivEvalWeights = projStruct->getTargetDerivEvalWeights(dim, 0);
    ordinal_type offsetBasis = projStruct->getBasisPointsRange(dim, 0).first;
    ordinal_type offsetBasisCurl = projStruct->getBasisDerivPointsRange(dim, 0).first;
    ordinal_type offsetTargetCurl = projStruct->getTargetDerivPointsRange(dim, 0).first;

    for(ordinal_type j=0; j <numElemDofs; ++j) {
      ordinal_type idof = cellBasis->getDofOrdinal(dim, 0, j);
      for(ordinal_type ic=0; ic<numCells; ++ic) {
        for(ordinal_type d=0; d <dim; ++d)
          for(ordinal_type iq=0; iq <numCubPoints; ++iq)
            internalBasisAtElemcubPoints(ic,j,iq,d)=basisAtCubPoints(ic,idof,offsetBasis+iq,d);

        for(ordinal_type d=0; d <derDim; ++d) {
          for(ordinal_type iq=0; iq <numCubPoints; ++iq) {
            internalBasisCurlAtElemcubPoints(ic,j,iq,d)=basisCurlAtCurlCubPoints(ic,idof,offsetBasisCurl+iq,d);
            wBasisCurlAtElemCubPoints(ic,j,iq,d)=internalBasisCurlAtElemcubPoints(ic,j,iq,d)*basisEvalWeights[iq];
          }
          for(ordinal_type iq=0; iq <numTargetCurlCubPoints; ++iq) {
            internalBasisCurlAtTargetCurlCubPoints(ic,j,iq,d)= basisCurlAtTargetCurlCubPoints(ic,idof,offsetTargetCurl+iq,d);
            wBasisCurlBasisAtTargetCurlCubPoints(ic,j,iq,d) = internalBasisCurlAtTargetCurlCubPoints(ic,j,iq,d)*targetDerivEvalWeights[iq];
          }
        }
      }
    }

    for(ordinal_type j=0; j <numEdgeDofs+numFaceDofs; ++j) {
      ordinal_type jdof = computedDofs(j);
      for(ordinal_type iq=0; iq <numCubPoints; ++iq)
        for(ordinal_type ic=0; ic<numCells; ++ic)  {
          for(ordinal_type d=0; d <derDim; ++d)
            mComputedProjectionCurl(ic,iq,d) -=  basisCoeffs(ic,jdof)*basisCurlAtCurlCubPoints(ic,jdof,offsetBasisCurl+iq,d);
          for(ordinal_type d=0; d <dim; ++d)
            mComputedProjection(ic,iq,d) -=  basisCoeffs(ic,jdof)*basisAtCubPoints(ic,jdof,offsetBasis+iq,d);
        }
    }

    scalarViewType cellMassMat_("cellMassMat_", numCells, numElemDofs+internalHgradCardinality, numElemDofs+internalHgradCardinality),
        cellRhsMat_("rhsMat_", numCells, numElemDofs+internalHgradCardinality);

    range_type range_H(0, numElemDofs);
    range_type range_B(numElemDofs, numElemDofs+internalHgradCardinality);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellMassMat_,Kokkos::ALL(),range_H,range_H), internalBasisCurlAtElemcubPoints, wBasisCurlAtElemCubPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellMassMat_,Kokkos::ALL(),range_H,range_B), internalBasisAtElemcubPoints, wHgradBasisGradAtCubPoints);
    if(dim==3)
      FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_H), Kokkos::subview(targetCurlAtCurlEvalPoints,Kokkos::ALL(),cellCurlPointsRange,Kokkos::ALL()), wBasisCurlBasisAtTargetCurlCubPoints);
    else
      FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_H), Kokkos::subview(targetCurlAtCurlEvalPoints,Kokkos::ALL(),cellCurlPointsRange), Kokkos::subview(wBasisCurlBasisAtTargetCurlCubPoints,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),0));
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_H), mComputedProjectionCurl, wBasisCurlAtElemCubPoints, true);

    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_B), Kokkos::subview(targetAtEvalPoints,Kokkos::ALL(),cellPointsRange,Kokkos::ALL()), wHgradBasisGradAtTargetCubPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(cellRhsMat_,Kokkos::ALL(),range_B), mComputedProjection, wHgradBasisGradAtCubPoints, true);

    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> cellMassMat("cellMassMat", numElemDofs+internalHgradCardinality,numElemDofs+internalHgradCardinality);
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> cellRhsMat("cellRhsMat",numElemDofs+internalHgradCardinality, 1);

    Teuchos::LAPACK<ordinal_type,funValsValueType> lapack;
    ordinal_type info = 0;
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> pivVec("pivVec", numElemDofs+internalHgradCardinality, 1);

    for(ordinal_type ic=0; ic<numCells; ++ic) {
      Kokkos::deep_copy(cellMassMat,funValsValueType(0));  //LAPACK might overwrite the matrix

      for(ordinal_type i=0; i<numElemDofs; ++i) {
        for(ordinal_type j=0; j<numElemDofs+internalHgradCardinality; ++j)
          cellMassMat(i,j) = cellMassMat_(ic,i,j);
        for(ordinal_type j=0; j<internalHgradCardinality; ++j)
          cellMassMat(numElemDofs+j,i) = cellMassMat_(0,i,numElemDofs+j);
      }

      for(ordinal_type i=0; i<numElemDofs+internalHgradCardinality; ++i)
        cellRhsMat(i,0) = cellRhsMat_(ic,i);


      lapack.GESV(numElemDofs+internalHgradCardinality, 1,
          cellMassMat.data(),
          cellMassMat.stride_1(),
          (ordinal_type*)pivVec.data(),
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
    delete hgradBasis;
  }
}
}
}

#endif

