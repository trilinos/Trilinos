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

/** \file   Intrepid2_ProjectionStructDef.hpp
    \brief  Header file for the Intrepid2::Experimental::ProjectionStruct containing definitions.
    \author Created by Mauro Perego
 */
#ifndef __INTREPID2_PROJECTIONSTRUCTDEF_HPP__
#define __INTREPID2_PROJECTIONSTRUCTDEF_HPP__

#include "Intrepid2_Utils.hpp"
#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_Cn_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"
#include "Intrepid2_HCURL_TET_In_FEM.hpp"

#include "Intrepid2_DefaultCubatureFactory.hpp"


namespace Intrepid2 {

namespace Experimental {
template<typename SpT, typename ValueType>
template<typename BasisPtrType>
void ProjectionStruct<SpT,ValueType>::createL2ProjectionStruct(const BasisPtrType cellBasis,
    const ordinal_type targetCubDegree) {
  const auto cellTopo = cellBasis->getBaseCellTopology();
  std::string name = cellBasis->getName();
  ordinal_type dim = cellTopo.getDimension();
  numBasisEvalPoints = 0;
  numBasisDerivEvalPoints = 0;
  numTargetEvalPoints = 0;
  numTargetDerivEvalPoints = 0;
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;

  ordinal_type basisCubDegree = cellBasis->getDegree();
  ordinal_type edgeBasisCubDegree, faceBasisCubDegree;

  ordinal_type numVertices = (cellBasis->getDofCount(0, 0) > 0) ? cellTopo.getVertexCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;
  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;

  numBasisEvalPoints += numVertices;
  numTargetEvalPoints += numVertices;
  view_type dofCoords("dofCoords", cellBasis->getCardinality(), dim);
  for(ordinal_type iv=0; iv<numVertices; ++iv) {
    basisPointsRange[0][iv] = range_type(iv, iv+1);
    basisCubPoints[0][iv] = view_type("basisCubPoints",1,dim);
    targetPointsRange[0][iv] = range_type(iv, iv+1);
    targetCubPoints[0][iv] = view_type("targetCubPoints",1,dim);
    ordinal_type idof = cellBasis->getDofOrdinal(0, iv, 0);
    for(ordinal_type d=0; d<dim; d++) {
      basisCubPoints[0][iv](0,d) = dofCoords(idof,d);
      targetCubPoints[0][iv](0,d) = dofCoords(idof,d);
    }
  }

  if ((name == "Intrepid2_HCURL_QUAD_In_FEM") || (name == "Intrepid2_HCURL_HEX_In_FEM") ||
      (name == "Intrepid2_HCURL_TRI_In_FEM") || (name == "Intrepid2_HCURL_TET_In_FEM")) {
    edgeBasisCubDegree = basisCubDegree - 1;
    faceBasisCubDegree = basisCubDegree;
  }
  else if ((name == "Intrepid2_HDIV_QUAD_In_FEM") || (name == "Intrepid2_HDIV_HEX_In_FEM") ||
      (name == "Intrepid2_HDIV_TRI_In_FEM") || (name == "Intrepid2_HDIV_TET_In_FEM")) {
    edgeBasisCubDegree = basisCubDegree - 1;
    faceBasisCubDegree = basisCubDegree -1;
  }
  else {
    edgeBasisCubDegree = basisCubDegree;
    faceBasisCubDegree = basisCubDegree;
  }

  DefaultCubatureFactory cub_factory;
  for(ordinal_type ie=0; ie<numEdges; ++ie) {
    ordinal_type cub_degree = 2*edgeBasisCubDegree;
    subCellTopologyKey[edgeDim][ie] = cellBasis->getBaseCellTopology().getKey(edgeDim, ie);
    auto edgeBasisCub = cub_factory.create<host_space_type, ValueType, ValueType>(cellBasis->getBaseCellTopology().getKey(edgeDim, ie), cub_degree);
    basisPointsRange[edgeDim][ie] = range_type(numBasisEvalPoints, numBasisEvalPoints+edgeBasisCub->getNumPoints());
    numBasisEvalPoints +=  edgeBasisCub->getNumPoints();
    basisCubPoints[edgeDim][ie] = view_type("basisCubPoints",edgeBasisCub->getNumPoints(),edgeDim);
    basisCubWeights[edgeDim][ie] = view_type("basisCubWeights",edgeBasisCub->getNumPoints());
    edgeBasisCub->getCubature(basisCubPoints[edgeDim][ie], basisCubWeights[edgeDim][ie]);

    cub_degree = edgeBasisCubDegree + targetCubDegree;
    auto edgeTargetCub = cub_factory.create<host_space_type, ValueType, ValueType>(cellBasis->getBaseCellTopology().getKey(edgeDim, ie), cub_degree);
    targetPointsRange[edgeDim][ie] = range_type(numTargetEvalPoints, numTargetEvalPoints+edgeTargetCub->getNumPoints());
    numTargetEvalPoints +=  edgeTargetCub->getNumPoints();
    targetCubPoints[edgeDim][ie] = view_type("targetCubPoints",edgeTargetCub->getNumPoints(),edgeDim);
    targetCubWeights[edgeDim][ie] = view_type("targetCubWeights",edgeTargetCub->getNumPoints());
    edgeTargetCub->getCubature(targetCubPoints[edgeDim][ie], targetCubWeights[edgeDim][ie]);
  }

  for(ordinal_type iface=0; iface<numFaces; ++iface) {
    ordinal_type cub_degree = 2*faceBasisCubDegree;
    subCellTopologyKey[faceDim][iface] = cellBasis->getBaseCellTopology().getKey(faceDim, iface);
    auto faceBasisCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[faceDim][iface], cub_degree);
    basisPointsRange[faceDim][iface] = range_type(numBasisEvalPoints, numBasisEvalPoints+faceBasisCub->getNumPoints());
    numBasisEvalPoints +=  faceBasisCub->getNumPoints();
    basisCubPoints[faceDim][iface] = view_type("basisCubPoints",faceBasisCub->getNumPoints(),faceDim);
    basisCubWeights[faceDim][iface] = view_type("basisCubWeights",faceBasisCub->getNumPoints());
    faceBasisCub->getCubature(basisCubPoints[faceDim][iface], basisCubWeights[faceDim][iface]);

    cub_degree = faceBasisCubDegree + targetCubDegree;
    auto faceTargetCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[faceDim][iface], cub_degree);
    targetPointsRange[faceDim][iface] = range_type(numTargetEvalPoints, numTargetEvalPoints+faceTargetCub->getNumPoints());
    numTargetEvalPoints +=  faceTargetCub->getNumPoints();
    targetCubPoints[faceDim][iface] = view_type("targetCubPoints",faceTargetCub->getNumPoints(),faceDim);
    targetCubWeights[faceDim][iface] = view_type("targetCubWeights",faceTargetCub->getNumPoints());
    faceTargetCub->getCubature(targetCubPoints[faceDim][iface], targetCubWeights[faceDim][iface]);
  }
  subCellTopologyKey[dim][0] = cellBasis->getBaseCellTopology().getBaseKey();
  if(cellBasis->getDofCount(dim,0)>0) {
    ordinal_type cub_degree = 2*basisCubDegree;
    auto elemBasisCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[dim][0], cub_degree);
    basisPointsRange[dim][0] = range_type(numBasisEvalPoints, numBasisEvalPoints+elemBasisCub->getNumPoints());
    numBasisEvalPoints +=  elemBasisCub->getNumPoints();
    basisCubPoints[dim][0] = view_type("basisCubPoints",elemBasisCub->getNumPoints(),dim);
    basisCubWeights[dim][0] = view_type("basisCubWeights",elemBasisCub->getNumPoints());
    elemBasisCub->getCubature(basisCubPoints[dim][0], basisCubWeights[dim][0]);

    cub_degree = basisCubDegree + targetCubDegree;
    auto elemTargetCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[dim][0], cub_degree);
    targetPointsRange[dim][0] = range_type(numTargetEvalPoints, numTargetEvalPoints+elemTargetCub->getNumPoints());
    numTargetEvalPoints +=  elemTargetCub->getNumPoints();
    targetCubPoints[dim][0] = view_type("targetCubPoints",elemTargetCub->getNumPoints(),dim);
    targetCubWeights[dim][0] = view_type("targetCubWeights",elemTargetCub->getNumPoints());
    elemTargetCub->getCubature(targetCubPoints[dim][0], targetCubWeights[dim][0]);
  }
}

template<typename SpT, typename ValueType>
template<typename BasisPtrType>
void ProjectionStruct<SpT,ValueType>::createHGradProjectionStruct(const BasisPtrType cellBasis,
    const ordinal_type targetCubDegree,
    const ordinal_type targetGradCubDegree) {
  const auto cellTopo = cellBasis->getBaseCellTopology();
  std::string name = cellBasis->getName();
  ordinal_type dim = cellTopo.getDimension();
  numBasisEvalPoints = 0;
  numBasisDerivEvalPoints = 0;
  numTargetEvalPoints = 0;
  numTargetDerivEvalPoints = 0;
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;

  ordinal_type basisCubDegree = cellBasis->getDegree();
  ordinal_type edgeBasisCubDegree = basisCubDegree;
  ordinal_type faceBasisCubDegree = basisCubDegree;

  ordinal_type numVertices = (cellBasis->getDofCount(0, 0) > 0) ? cellTopo.getVertexCount() : 0;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount(): 0;
  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;

  numBasisEvalPoints += numVertices;
  numTargetEvalPoints += numVertices;
  view_type dofCoords("dofCoords", cellBasis->getCardinality(), dim);
  for(ordinal_type iv=0; iv<numVertices; ++iv) {
    basisPointsRange[0][iv] = range_type(iv, iv+1);
    basisCubPoints[0][iv] = view_type("basisCubPoints",1,dim);
    targetPointsRange[0][iv] = range_type(iv, iv+1);
    targetCubPoints[0][iv] = view_type("targetCubPoints",1,dim);
    ordinal_type idof = cellBasis->getDofOrdinal(0, iv, 0);
    for(ordinal_type d=0; d<dim; d++) {
      basisCubPoints[0][iv](0,d) = dofCoords(idof,d);
      targetCubPoints[0][iv](0,d) = dofCoords(idof,d);
    }
  }

  DefaultCubatureFactory cub_factory;
  for(ordinal_type ie=0; ie<numEdges; ++ie) {
    ordinal_type cub_degree = 2*edgeBasisCubDegree;
    subCellTopologyKey[edgeDim][ie] = cellBasis->getBaseCellTopology().getKey(edgeDim, ie);
    auto edgeBasisCub = cub_factory.create<host_space_type, ValueType, ValueType>(cellBasis->getBaseCellTopology().getKey(edgeDim, ie), cub_degree);
    basisDerivPointsRange[edgeDim][ie] = range_type(numBasisDerivEvalPoints, numBasisDerivEvalPoints+edgeBasisCub->getNumPoints());
    numBasisDerivEvalPoints +=  edgeBasisCub->getNumPoints();
    basisDerivCubPoints[edgeDim][ie] = view_type("basisDerivCubPoints",edgeBasisCub->getNumPoints(),edgeDim);
    basisDerivCubWeights[edgeDim][ie] = view_type("basisDerivCubWeights",edgeBasisCub->getNumPoints());
    edgeBasisCub->getCubature(basisDerivCubPoints[edgeDim][ie], basisDerivCubWeights[edgeDim][ie]);

    cub_degree = edgeBasisCubDegree + targetGradCubDegree;
    auto edgeTargetCub = cub_factory.create<host_space_type, ValueType, ValueType>(cellBasis->getBaseCellTopology().getKey(edgeDim, ie), cub_degree);
    targetDerivPointsRange[edgeDim][ie] = range_type(numTargetDerivEvalPoints, numTargetDerivEvalPoints+edgeTargetCub->getNumPoints());
    numTargetDerivEvalPoints +=  edgeTargetCub->getNumPoints();
    targetDerivCubPoints[edgeDim][ie] = view_type("targetDerivCubPoints",edgeTargetCub->getNumPoints(),edgeDim);
    targetDerivCubWeights[edgeDim][ie] = view_type("targetDerivCubWeights",edgeTargetCub->getNumPoints());
    edgeTargetCub->getCubature(targetDerivCubPoints[edgeDim][ie], targetDerivCubWeights[edgeDim][ie]);
  }

  for(ordinal_type iface=0; iface<numFaces; ++iface) {
    ordinal_type cub_degree = 2*faceBasisCubDegree;
    subCellTopologyKey[faceDim][iface] = cellBasis->getBaseCellTopology().getKey(faceDim, iface);
    auto faceBasisGradCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[faceDim][iface], cub_degree);
    basisDerivPointsRange[faceDim][iface] = range_type(numBasisDerivEvalPoints, numBasisDerivEvalPoints+faceBasisGradCub->getNumPoints());
    numBasisDerivEvalPoints +=  faceBasisGradCub->getNumPoints();
    basisDerivCubPoints[faceDim][iface] = view_type("basisDerivCubPoints",faceBasisGradCub->getNumPoints(),faceDim);
    basisDerivCubWeights[faceDim][iface] = view_type("basisDerivCubWeights",faceBasisGradCub->getNumPoints());
    faceBasisGradCub->getCubature(basisDerivCubPoints[faceDim][iface], basisDerivCubWeights[faceDim][iface]);

    cub_degree = faceBasisCubDegree + targetGradCubDegree;
    auto faceTargetDerivCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[faceDim][iface], cub_degree);
    targetDerivPointsRange[faceDim][iface] = range_type(numTargetDerivEvalPoints, numTargetDerivEvalPoints+faceTargetDerivCub->getNumPoints());
    numTargetDerivEvalPoints +=  faceTargetDerivCub->getNumPoints();
    targetDerivCubPoints[faceDim][iface] = view_type("targetDerivCubPoints",faceTargetDerivCub->getNumPoints(),faceDim);
    targetDerivCubWeights[faceDim][iface] = view_type("targetDerivCubWeights",faceTargetDerivCub->getNumPoints());
    faceTargetDerivCub->getCubature(targetDerivCubPoints[faceDim][iface], targetDerivCubWeights[faceDim][iface]);
  }
  subCellTopologyKey[dim][0] = cellBasis->getBaseCellTopology().getBaseKey();
  if(cellBasis->getDofCount(dim,0)>0) {
    ordinal_type cub_degree = 2*basisCubDegree;
    auto elemBasisGradCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[dim][0], cub_degree);
    basisDerivPointsRange[dim][0] = range_type(numBasisDerivEvalPoints, numBasisDerivEvalPoints+elemBasisGradCub->getNumPoints());
    numBasisDerivEvalPoints +=  elemBasisGradCub->getNumPoints();
    basisDerivCubPoints[dim][0] = view_type("basisDerivCubPoints",elemBasisGradCub->getNumPoints(),dim);
    basisDerivCubWeights[dim][0] = view_type("basisDerivCubWeights",elemBasisGradCub->getNumPoints());
    elemBasisGradCub->getCubature(basisDerivCubPoints[dim][0], basisDerivCubWeights[dim][0]);

    cub_degree = basisCubDegree + targetGradCubDegree;
    auto elemTargetGradCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[dim][0], cub_degree);
    targetDerivPointsRange[dim][0] = range_type(numTargetDerivEvalPoints, numTargetDerivEvalPoints+elemTargetGradCub->getNumPoints());
    numTargetDerivEvalPoints +=  elemTargetGradCub->getNumPoints();
    targetDerivCubPoints[dim][0] = view_type("targetDerivCubPoints",elemTargetGradCub->getNumPoints(),dim);
    targetDerivCubWeights[dim][0] = view_type("targetDerivCubWeights",elemTargetGradCub->getNumPoints());
    elemTargetGradCub->getCubature(targetDerivCubPoints[dim][0], targetDerivCubWeights[dim][0]);
  }
}

template<typename SpT, typename ValueType>
template<typename BasisPtrType>
void ProjectionStruct<SpT,ValueType>::createHCurlProjectionStruct(const BasisPtrType cellBasis,
    const ordinal_type targetCubDegree,
    const ordinal_type targetCurlCubDegre) {
  const auto cellTopo = cellBasis->getBaseCellTopology();
  std::string name = cellBasis->getName();
  ordinal_type dim = cellTopo.getDimension();
  numBasisEvalPoints = 0;
  numBasisDerivEvalPoints = 0;
  numTargetEvalPoints = 0;
  numTargetDerivEvalPoints = 0;
  const ordinal_type edgeDim = 1;
  const ordinal_type faceDim = 2;

  ordinal_type basisCubDegree = cellBasis->getDegree();
  ordinal_type edgeBasisCubDegree = basisCubDegree - 1;
  ordinal_type faceBasisCubDegree = basisCubDegree;
  ordinal_type numFaces = (cellBasis->getDofCount(2, 0) > 0) ? cellTopo.getFaceCount() : 0;
  ordinal_type numEdges = (cellBasis->getDofCount(1, 0) > 0) ? cellTopo.getEdgeCount() : 0;

  DefaultCubatureFactory cub_factory;
  for(ordinal_type ie=0; ie<numEdges; ++ie) {
    ordinal_type cub_degree = 2*edgeBasisCubDegree;
    subCellTopologyKey[edgeDim][ie] = cellBasis->getBaseCellTopology().getKey(edgeDim, ie);
    auto edgeBasisCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[edgeDim][ie], cub_degree);
    basisPointsRange[edgeDim][ie] = range_type(numBasisEvalPoints, numBasisEvalPoints+edgeBasisCub->getNumPoints());
    numBasisEvalPoints +=  edgeBasisCub->getNumPoints();
    basisCubPoints[edgeDim][ie] = view_type("basisCubPoints",edgeBasisCub->getNumPoints(),edgeDim);
    basisCubWeights[edgeDim][ie] = view_type("basisCubWeights",edgeBasisCub->getNumPoints());
    edgeBasisCub->getCubature(basisCubPoints[edgeDim][ie], basisCubWeights[edgeDim][ie]);

    cub_degree = edgeBasisCubDegree + targetCubDegree;
    auto edgeTargetCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[edgeDim][ie], cub_degree);
    targetPointsRange[edgeDim][ie] = range_type(numTargetEvalPoints, numTargetEvalPoints+edgeTargetCub->getNumPoints());
    numTargetEvalPoints +=  edgeTargetCub->getNumPoints();
    targetCubPoints[edgeDim][ie] = view_type("targetCubPoints",edgeTargetCub->getNumPoints(),edgeDim);
    targetCubWeights[edgeDim][ie] = view_type("targetCubWeights",edgeTargetCub->getNumPoints());
    edgeTargetCub->getCubature(targetCubPoints[edgeDim][ie], targetCubWeights[edgeDim][ie]);
  }

  for(ordinal_type iface=0; iface<numFaces; ++iface) {
    ordinal_type cub_degree = 2*faceBasisCubDegree;
    subCellTopologyKey[faceDim][iface] = cellBasis->getBaseCellTopology().getKey(faceDim, iface);
    auto faceBasisCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[faceDim][iface], cub_degree);
    basisPointsRange[faceDim][iface] = range_type(numBasisEvalPoints, numBasisEvalPoints+faceBasisCub->getNumPoints());
    numBasisEvalPoints +=  faceBasisCub->getNumPoints();
    basisCubPoints[faceDim][iface] = view_type("basisCubPoints",faceBasisCub->getNumPoints(),faceDim);
    basisCubWeights[faceDim][iface] = view_type("basisCubWeights",faceBasisCub->getNumPoints());
    faceBasisCub->getCubature(basisCubPoints[faceDim][iface], basisCubWeights[faceDim][iface]);

    basisDerivPointsRange[faceDim][iface] = range_type(numBasisDerivEvalPoints, numBasisDerivEvalPoints+faceBasisCub->getNumPoints());
    numBasisDerivEvalPoints +=  faceBasisCub->getNumPoints();
    basisDerivCubPoints[faceDim][iface] = view_type("basisDerivCubPoints",faceBasisCub->getNumPoints(),faceDim);
    basisDerivCubWeights[faceDim][iface] = view_type("basisDerivCubWeights",faceBasisCub->getNumPoints());
    faceBasisCub->getCubature(basisDerivCubPoints[faceDim][iface], basisDerivCubWeights[faceDim][iface]);

    cub_degree = faceBasisCubDegree + targetCubDegree;
    auto faceTargetCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[faceDim][iface], cub_degree);
    targetPointsRange[faceDim][iface] = range_type(numTargetEvalPoints, numTargetEvalPoints+faceTargetCub->getNumPoints());
    numTargetEvalPoints +=  faceTargetCub->getNumPoints();
    targetCubPoints[faceDim][iface] = view_type("targetCubPoints",faceTargetCub->getNumPoints(),faceDim);
    targetCubWeights[faceDim][iface] = view_type("targetCubWeights",faceTargetCub->getNumPoints());
    faceTargetCub->getCubature(targetCubPoints[faceDim][iface], targetCubWeights[faceDim][iface]);

    cub_degree = faceBasisCubDegree + targetCurlCubDegre;
    auto faceTargetDerivCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[faceDim][iface], cub_degree);
    targetDerivPointsRange[faceDim][iface] = range_type(numTargetDerivEvalPoints, numTargetDerivEvalPoints+faceTargetDerivCub->getNumPoints());
    numTargetDerivEvalPoints +=  faceTargetDerivCub->getNumPoints();
    targetDerivCubPoints[faceDim][iface] = view_type("targetDerivCubPoints",faceTargetDerivCub->getNumPoints(),faceDim);
    targetDerivCubWeights[faceDim][iface] = view_type("targetDerivCubWeights",faceTargetDerivCub->getNumPoints());
    faceTargetDerivCub->getCubature(targetDerivCubPoints[faceDim][iface], targetDerivCubWeights[faceDim][iface]);
  }

  subCellTopologyKey[dim][0] = cellBasis->getBaseCellTopology().getBaseKey();
  if(cellBasis->getDofCount(dim,0)>0) {
    ordinal_type cub_degree = 2*basisCubDegree;
    auto elemBasisCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[dim][0], cub_degree);
    basisPointsRange[dim][0] = range_type(numBasisEvalPoints, numBasisEvalPoints+elemBasisCub->getNumPoints());
    numBasisEvalPoints +=  elemBasisCub->getNumPoints();
    basisCubPoints[dim][0] = view_type("basisCubPoints",elemBasisCub->getNumPoints(),dim);
    basisCubWeights[dim][0] = view_type("basisCubWeights",elemBasisCub->getNumPoints());
    elemBasisCub->getCubature(basisCubPoints[dim][0], basisCubWeights[dim][0]);

    basisDerivPointsRange[dim][0] = range_type(numBasisDerivEvalPoints, numBasisDerivEvalPoints+elemBasisCub->getNumPoints());
    numBasisDerivEvalPoints +=  elemBasisCub->getNumPoints();
    basisDerivCubPoints[dim][0] = view_type("basisDerivCubPoints",elemBasisCub->getNumPoints(),dim);
    basisDerivCubWeights[dim][0] = view_type("basisDerivCubWeights",elemBasisCub->getNumPoints());
    elemBasisCub->getCubature(basisDerivCubPoints[dim][0], basisDerivCubWeights[dim][0]);

    cub_degree = basisCubDegree + targetCubDegree;
    auto elemTargetCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[dim][0], cub_degree);
    targetPointsRange[dim][0] = range_type(numTargetEvalPoints, numTargetEvalPoints+elemTargetCub->getNumPoints());
    numTargetEvalPoints +=  elemTargetCub->getNumPoints();
    targetCubPoints[dim][0] = view_type("targetCubPoints",elemTargetCub->getNumPoints(),dim);
    targetCubWeights[dim][0] = view_type("targetCubWeights",elemTargetCub->getNumPoints());
    elemTargetCub->getCubature(targetCubPoints[dim][0], targetCubWeights[dim][0]);

    cub_degree = basisCubDegree + targetCurlCubDegre;
    auto elemTargetCurlCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[dim][0], cub_degree);
    targetDerivPointsRange[dim][0] = range_type(numTargetDerivEvalPoints, numTargetDerivEvalPoints+elemTargetCurlCub->getNumPoints());
    numTargetDerivEvalPoints +=  elemTargetCurlCub->getNumPoints();
    targetDerivCubPoints[dim][0] = view_type("targetDerivCubPoints",elemTargetCurlCub->getNumPoints(),dim);
    targetDerivCubWeights[dim][0] = view_type("targetDerivCubWeights",elemTargetCurlCub->getNumPoints());
    elemTargetCurlCub->getCubature(targetDerivCubPoints[dim][0], targetDerivCubWeights[dim][0]);
  }

}

template<typename SpT, typename ValueType>
template<typename BasisPtrType>
void ProjectionStruct<SpT,ValueType>::createHDivProjectionStruct(const BasisPtrType cellBasis,
    const ordinal_type targetCubDegree,
    const ordinal_type targetDivCubDegre) {
  const auto cellTopo = cellBasis->getBaseCellTopology();
  std::string name = cellBasis->getName();
  ordinal_type dim = cellTopo.getDimension();
  numBasisEvalPoints = 0;
  numBasisDerivEvalPoints = 0;
  numTargetEvalPoints = 0;
  numTargetDerivEvalPoints = 0;
  const ordinal_type sideDim = dim - 1;

  ordinal_type basisCubDegree = cellBasis->getDegree();
  ordinal_type sideBasisCubDegree = basisCubDegree - 1;
  ordinal_type numSides = cellTopo.getSideCount()*ordinal_type(cellBasis->getDofCount(sideDim, 0) > 0);

  Basis<host_space_type,ValueType,ValueType> *hcurlBasis = NULL;
  if(name.find("HEX")!=std::string::npos)
    hcurlBasis = new Basis_HCURL_HEX_In_FEM<host_space_type,ValueType,ValueType>(cellBasis->getDegree());
  else if(name.find("TET")!=std::string::npos)
    hcurlBasis = new Basis_HCURL_TET_In_FEM<host_space_type,ValueType,ValueType>(cellBasis->getDegree());
  else if(name.find("QUAD")!=std::string::npos)
    hcurlBasis = new Basis_HGRAD_QUAD_Cn_FEM<host_space_type,ValueType,ValueType>(cellBasis->getDegree());
  else if(name.find("TRI")!=std::string::npos)
    hcurlBasis = new Basis_HGRAD_TRI_Cn_FEM<host_space_type,ValueType,ValueType>(cellBasis->getDegree());
  else  {
    std::stringstream ss;
    ss << ">>> ERROR (Intrepid2::ProjectionTools::createHDivProjectionStruct): "
        << "Method not implemented for basis " << name;
    INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
  }

  bool haveHCurlConstraint = (hcurlBasis != NULL) && (hcurlBasis->getDofCount(dim,0) > 0);
  if(hcurlBasis != NULL) delete hcurlBasis;

  DefaultCubatureFactory cub_factory;

  for(ordinal_type is=0; is<numSides; ++is) {
    ordinal_type cub_degree = 2*sideBasisCubDegree;
    subCellTopologyKey[sideDim][is] = cellBasis->getBaseCellTopology().getKey(sideDim, is);
    auto sideBasisCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[sideDim][is], cub_degree);
    basisPointsRange[sideDim][is] = range_type(numBasisEvalPoints, numBasisEvalPoints+sideBasisCub->getNumPoints());
    numBasisEvalPoints +=  sideBasisCub->getNumPoints();
    basisCubPoints[sideDim][is] = view_type("basisCubPoints",sideBasisCub->getNumPoints(),sideDim);
    basisCubWeights[sideDim][is] = view_type("basisCubWeights",sideBasisCub->getNumPoints());
    sideBasisCub->getCubature(basisCubPoints[sideDim][is], basisCubWeights[sideDim][is]);

    cub_degree = sideBasisCubDegree + targetCubDegree;
    auto sideTargetCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[sideDim][is], cub_degree);
    targetPointsRange[sideDim][is] = range_type(numTargetEvalPoints, numTargetEvalPoints+sideTargetCub->getNumPoints());
    numTargetEvalPoints +=  sideTargetCub->getNumPoints();
    targetCubPoints[sideDim][is] = view_type("targetCubPoints",sideTargetCub->getNumPoints(),sideDim);
    targetCubWeights[sideDim][is] = view_type("targetCubWeights",sideTargetCub->getNumPoints());
    sideTargetCub->getCubature(targetCubPoints[sideDim][is], targetCubWeights[sideDim][is]);
  }

  subCellTopologyKey[dim][0] = cellBasis->getBaseCellTopology().getBaseKey();
  if(cellBasis->getDofCount(dim,0)>0) {
    ordinal_type cub_degree = 2*basisCubDegree - 1;
    auto elemBasisDivCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[dim][0], cub_degree);
    basisDerivPointsRange[dim][0] = range_type(numBasisDerivEvalPoints, numBasisDerivEvalPoints+elemBasisDivCub->getNumPoints());
    numBasisDerivEvalPoints +=  elemBasisDivCub->getNumPoints();
    basisDerivCubPoints[dim][0] = view_type("basisDerivCubPoints",elemBasisDivCub->getNumPoints(),dim);
    basisDerivCubWeights[dim][0] = view_type("basisDerivCubWeights",elemBasisDivCub->getNumPoints());
    elemBasisDivCub->getCubature(basisDerivCubPoints[dim][0], basisDerivCubWeights[dim][0]);

    cub_degree = basisCubDegree - 1 + targetDivCubDegre;
    auto elemTargetDivCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[dim][0], cub_degree);
    targetDerivPointsRange[dim][0] = range_type(numTargetDerivEvalPoints, numTargetDerivEvalPoints+elemTargetDivCub->getNumPoints());
    numTargetDerivEvalPoints +=  elemTargetDivCub->getNumPoints();
    targetDerivCubPoints[dim][0] = view_type("targetDerivCubPoints",elemTargetDivCub->getNumPoints(),dim);
    targetDerivCubWeights[dim][0] = view_type("targetDerivCubWeights",elemTargetDivCub->getNumPoints());
    elemTargetDivCub->getCubature(targetDerivCubPoints[dim][0], targetDerivCubWeights[dim][0]);

    if(haveHCurlConstraint)
    {
      cub_degree = 2*basisCubDegree;
      auto elemBasisCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[dim][0], cub_degree);
      basisPointsRange[dim][0] = range_type(numBasisEvalPoints, numBasisEvalPoints + elemBasisCub->getNumPoints());
      numBasisEvalPoints +=  elemBasisCub->getNumPoints();
      basisCubPoints[dim][0] = view_type("basisCubPoints",elemBasisCub->getNumPoints(),dim);
      basisCubWeights[dim][0] = view_type("basisCubWeights",elemBasisCub->getNumPoints());
      elemBasisCub->getCubature(basisCubPoints[dim][0], basisCubWeights[dim][0]);

      cub_degree = basisCubDegree + targetCubDegree;
      auto elemTargetCub = cub_factory.create<host_space_type, ValueType, ValueType>(subCellTopologyKey[dim][0], cub_degree);
      targetPointsRange[dim][0] = range_type(numTargetEvalPoints, numTargetEvalPoints + elemTargetCub->getNumPoints());
      numTargetEvalPoints +=  elemTargetCub->getNumPoints();
      targetCubPoints[dim][0] = view_type("targetCubPoints",elemTargetCub->getNumPoints(),dim);
      targetCubWeights[dim][0] = view_type("targetCubWeights",elemTargetCub->getNumPoints());
      elemTargetCub->getCubature(targetCubPoints[dim][0], targetCubWeights[dim][0]);
    }
  }
}

template<typename SpT, typename ValueType>
template<typename BasisPtrType>
void ProjectionStruct<SpT,ValueType>::createHVolProjectionStruct(const BasisPtrType cellBasis,
    const ordinal_type targetCubDegree) {
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  numBasisEvalPoints = 0;
  numBasisDerivEvalPoints = 0;
  numTargetEvalPoints = 0;
  numTargetDerivEvalPoints = 0;

  ordinal_type basisCubDegree = cellBasis->getDegree();
  DefaultCubatureFactory cub_factory;
  subCellTopologyKey[dim][0] = cellBasis->getBaseCellTopology().getBaseKey();
  if(cellBasis->getDofCount(dim,0)>0) {

    ordinal_type cub_degree = 2*basisCubDegree;
    auto elemBasisCub = cub_factory.create<host_space_type, ValueType, ValueType>(cellTopo.getBaseKey(), cub_degree);
    basisPointsRange[dim][0] = range_type(0, elemBasisCub->getNumPoints());
    numBasisEvalPoints +=  elemBasisCub->getNumPoints();
    basisCubPoints[dim][0] = view_type("basisCubPoints",elemBasisCub->getNumPoints(),dim);
    basisCubWeights[dim][0] = view_type("basisCubWeights",elemBasisCub->getNumPoints());
    elemBasisCub->getCubature(basisCubPoints[dim][0], basisCubWeights[dim][0]);

    cub_degree = basisCubDegree + targetCubDegree;
    auto elemTargetCub = cub_factory.create<host_space_type, ValueType, ValueType>(cellTopo.getBaseKey(), cub_degree);
    targetPointsRange[dim][0] = range_type(0, elemTargetCub->getNumPoints());
    numTargetEvalPoints +=  elemTargetCub->getNumPoints();
    targetCubPoints[dim][0] = view_type("targetCubPoints",elemTargetCub->getNumPoints(),dim);
    targetCubWeights[dim][0] = view_type("targetCubWeights",elemTargetCub->getNumPoints());
    elemTargetCub->getCubature(targetCubPoints[dim][0], targetCubWeights[dim][0]);
  }
}

}
}
#endif





