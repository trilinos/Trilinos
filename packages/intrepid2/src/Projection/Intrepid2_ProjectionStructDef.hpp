// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ProjectionStructDef.hpp
    \brief  Header file for the Intrepid2::ProjectionStruct containing definitions.
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

template<typename DeviceType, typename ValueType>
template<typename BasisPtrType>
void ProjectionStruct<DeviceType,ValueType>::createL2ProjectionStruct(const BasisPtrType cellBasis,
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

  ordinal_type hasCellDofs = (cellBasis->getDofCount(dim, 0) > 0);

  INTREPID2_TEST_FOR_ABORT( (numVertices > maxSubCellsCount)  || (numFaces > maxSubCellsCount) || (numEdges > maxSubCellsCount),
      ">>> ERROR (Intrepid2::ProjectionStruct:createHDivProjectionStruct, Projections do not support a cell topology with this cub cells count");


  numBasisEvalPoints += numVertices;
  numTargetEvalPoints += numVertices;

  basisPointsRange = range_tag("basisPointsRange", 4,maxSubCellsCount);
  targetPointsRange = range_tag("targetPointsRange", 4,maxSubCellsCount);
  basisDerivPointsRange = range_tag("basisDerivPointsRange", 4,maxSubCellsCount);
  targetDerivPointsRange = range_tag("targetDerivPointsRange", numberSubCellDims,maxSubCellsCount);
  subCellTopologyKey = key_tag("subCellTopologyKey",numberSubCellDims,maxSubCellsCount);

  maxNumBasisEvalPoints = numVertices; maxNumTargetEvalPoints = numVertices;

  // The set of the eval points on the reference vertex contains only point (0.0).
  // Not very useful, updating these for completeness  
  for(ordinal_type iv=0; iv<numVertices; ++iv) {
    basisPointsRange(0,iv) = range_type(iv, iv+1);
    basisCubPoints[0][iv] = host_view_type("basisCubPoints",1,1);
    targetPointsRange(0,iv) = range_type(iv, iv+1);
    targetCubPoints[0][iv] = host_view_type("targetCubPoints",1,1);
  }

  if (cellBasis->getFunctionSpace() == FUNCTION_SPACE_HCURL) {
    edgeBasisCubDegree = basisCubDegree - 1;
    faceBasisCubDegree = basisCubDegree;
  }
  else if (cellBasis->getFunctionSpace() == FUNCTION_SPACE_HDIV) {
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
    subCellTopologyKey(edgeDim,ie) = cellBasis->getBaseCellTopology().getKey(edgeDim, ie);
    auto edgeBasisCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(cellBasis->getBaseCellTopology().getKey(edgeDim, ie), cub_degree, EPolyType::POLYTYPE_MAX, true);
    basisPointsRange(edgeDim,ie) = range_type(numBasisEvalPoints, numBasisEvalPoints+edgeBasisCub->getNumPoints());
    numBasisEvalPoints +=  edgeBasisCub->getNumPoints();
    maxNumBasisEvalPoints = std::max(maxNumBasisEvalPoints, edgeBasisCub->getNumPoints());
    basisCubPoints[edgeDim][ie] = host_view_type("basisCubPoints",edgeBasisCub->getNumPoints(),edgeDim);
    basisCubWeights[edgeDim][ie] = host_view_type("basisCubWeights",edgeBasisCub->getNumPoints());
    edgeBasisCub->getCubature(basisCubPoints[edgeDim][ie], basisCubWeights[edgeDim][ie]);

    cub_degree = edgeBasisCubDegree + targetCubDegree;
    auto edgeTargetCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(cellBasis->getBaseCellTopology().getKey(edgeDim, ie), cub_degree, EPolyType::POLYTYPE_MAX, true);
    targetPointsRange(edgeDim,ie) = range_type(numTargetEvalPoints, numTargetEvalPoints+edgeTargetCub->getNumPoints());
    numTargetEvalPoints +=  edgeTargetCub->getNumPoints();
    maxNumTargetEvalPoints = std::max(maxNumTargetEvalPoints, edgeTargetCub->getNumPoints());
    targetCubPoints[edgeDim][ie] = host_view_type("targetCubPoints",edgeTargetCub->getNumPoints(),edgeDim);
    targetCubWeights[edgeDim][ie] = host_view_type("targetCubWeights",edgeTargetCub->getNumPoints());
    edgeTargetCub->getCubature(targetCubPoints[edgeDim][ie], targetCubWeights[edgeDim][ie]);
  }

  for(ordinal_type iface=0; iface<numFaces; ++iface) {
    ordinal_type cub_degree = 2*faceBasisCubDegree;
    subCellTopologyKey(faceDim,iface) = cellBasis->getBaseCellTopology().getKey(faceDim, iface);
    auto faceBasisCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(faceDim,iface), cub_degree, EPolyType::POLYTYPE_MAX, true);
    basisPointsRange(faceDim,iface) = range_type(numBasisEvalPoints, numBasisEvalPoints+faceBasisCub->getNumPoints());
    numBasisEvalPoints +=  faceBasisCub->getNumPoints();
    maxNumBasisEvalPoints = std::max(maxNumBasisEvalPoints, faceBasisCub->getNumPoints());
    basisCubPoints[faceDim][iface] = host_view_type("basisCubPoints",faceBasisCub->getNumPoints(),faceDim);
    basisCubWeights[faceDim][iface] = host_view_type("basisCubWeights",faceBasisCub->getNumPoints());
    faceBasisCub->getCubature(basisCubPoints[faceDim][iface], basisCubWeights[faceDim][iface]);

    cub_degree = faceBasisCubDegree + targetCubDegree;
    auto faceTargetCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(faceDim,iface), cub_degree, EPolyType::POLYTYPE_MAX, true);
    targetPointsRange(faceDim,iface) = range_type(numTargetEvalPoints, numTargetEvalPoints+faceTargetCub->getNumPoints());
    numTargetEvalPoints +=  faceTargetCub->getNumPoints();
    maxNumTargetEvalPoints = std::max(maxNumTargetEvalPoints, faceTargetCub->getNumPoints());
    targetCubPoints[faceDim][iface] = host_view_type("targetCubPoints",faceTargetCub->getNumPoints(),faceDim);
    targetCubWeights[faceDim][iface] = host_view_type("targetCubWeights",faceTargetCub->getNumPoints());
    faceTargetCub->getCubature(targetCubPoints[faceDim][iface], targetCubWeights[faceDim][iface]);
  }
  subCellTopologyKey(dim,0) = cellBasis->getBaseCellTopology().getBaseKey();
  if(hasCellDofs) {
    ordinal_type cub_degree = 2*basisCubDegree;
    auto elemBasisCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(dim,0), cub_degree, EPolyType::POLYTYPE_MAX, true);
    basisPointsRange(dim,0) = range_type(numBasisEvalPoints, numBasisEvalPoints+elemBasisCub->getNumPoints());
    numBasisEvalPoints +=  elemBasisCub->getNumPoints();
    maxNumBasisEvalPoints = std::max(maxNumBasisEvalPoints, elemBasisCub->getNumPoints());
    basisCubPoints[dim][0] = host_view_type("basisCubPoints",elemBasisCub->getNumPoints(),dim);
    basisCubWeights[dim][0] = host_view_type("basisCubWeights",elemBasisCub->getNumPoints());
    elemBasisCub->getCubature(basisCubPoints[dim][0], basisCubWeights[dim][0]);

    cub_degree = basisCubDegree + targetCubDegree;
    auto elemTargetCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(dim,0), cub_degree, EPolyType::POLYTYPE_MAX, true);
    targetPointsRange(dim,0) = range_type(numTargetEvalPoints, numTargetEvalPoints+elemTargetCub->getNumPoints());
    numTargetEvalPoints +=  elemTargetCub->getNumPoints();
    maxNumTargetEvalPoints = std::max(maxNumTargetEvalPoints, elemTargetCub->getNumPoints());
    targetCubPoints[dim][0] = host_view_type("targetCubPoints",elemTargetCub->getNumPoints(),dim);
    targetCubWeights[dim][0] = host_view_type("targetCubWeights",elemTargetCub->getNumPoints());
    elemTargetCub->getCubature(targetCubPoints[dim][0], targetCubWeights[dim][0]);


  }

  allBasisEPoints = view_type("allBasisPoints", numBasisEvalPoints, dim);
  allTargetEPoints = view_type("allTargetPoints", numTargetEvalPoints, dim);
  allBasisDerivEPoints = view_type("allBasisDerivPoints", numBasisDerivEvalPoints, dim);
  allTargetDerivEPoints = view_type("allTargetDerivPoints", numTargetDerivEvalPoints, dim);

  if(numVertices>0) {
    for(ordinal_type iv=0; iv<numVertices; ++iv) {
      CellTools<DeviceType>::getReferenceVertex(Kokkos::subview(allBasisEPoints, iv, Kokkos::ALL()), cellTopo, iv);
      CellTools<DeviceType>::getReferenceVertex(Kokkos::subview(allTargetEPoints, iv, Kokkos::ALL()), cellTopo, iv);
    }
  }

  if(numEdges>0) {
    auto subcellParamEdge = RefSubcellParametrization<DeviceType>::get(edgeDim, cellTopo.getKey());
    for(ordinal_type ie=0; ie<numEdges; ++ie) {
      auto edgeBasisEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisCubPoints[edgeDim][ie]);
      CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allBasisEPoints,basisPointsRange(edgeDim,ie),Kokkos::ALL()), edgeBasisEPoints, subcellParamEdge, ie);

      auto edgeTargetEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetCubPoints[edgeDim][ie]);
      CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allTargetEPoints,targetPointsRange(edgeDim,ie),Kokkos::ALL()), edgeTargetEPoints, subcellParamEdge, ie);
    }
  }

  if(numFaces>0) {
    auto subcellParamFace = RefSubcellParametrization<DeviceType>::get(faceDim, cellTopo.getKey());
    for(ordinal_type iface=0; iface<numFaces; ++iface) {
      auto faceBasisEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisCubPoints[faceDim][iface]);
      CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allBasisEPoints,basisPointsRange(faceDim,iface),Kokkos::ALL()), faceBasisEPoints, subcellParamFace, iface);

      auto faceTargetEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetCubPoints[faceDim][iface]);
      CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allTargetEPoints,targetPointsRange(faceDim,iface),Kokkos::ALL()), faceTargetEPoints, subcellParamFace, iface);
    }
  }

  if(hasCellDofs) {
    auto cellBasisEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisCubPoints[dim][0]);
    Kokkos::deep_copy(Kokkos::subview(allBasisEPoints,basisPointsRange(dim,0), Kokkos::ALL()), cellBasisEPoints);

    auto cellTargetEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetCubPoints[dim][0]);
    Kokkos::deep_copy(Kokkos::subview(allTargetEPoints,targetPointsRange(dim,0),Kokkos::ALL()), cellTargetEPoints);
  }
}

template<typename DeviceType, typename ValueType>
template<typename BasisPtrType>
void ProjectionStruct<DeviceType,ValueType>::createHGradProjectionStruct(const BasisPtrType cellBasis,
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

  INTREPID2_TEST_FOR_ABORT( (numFaces > maxSubCellsCount) || (numEdges > maxSubCellsCount),
      ">>> ERROR (Intrepid2::ProjectionStruct:createHDivProjectionStruct, Projections do not support a cell topology with this cub cells count");


  ordinal_type hasCellDofs = (cellBasis->getDofCount(dim, 0) > 0);

  maxNumBasisEvalPoints = numVertices; maxNumTargetEvalPoints = numVertices;
  maxNumBasisDerivEvalPoints = 0; maxNumTargetDerivEvalPoints = 0;

  basisPointsRange = range_tag("basisPointsRange", numberSubCellDims,maxSubCellsCount);
  targetPointsRange = range_tag("targetPointsRange", numberSubCellDims,maxSubCellsCount);
  basisDerivPointsRange = range_tag("basisDerivPointsRange", numberSubCellDims,maxSubCellsCount);
  targetDerivPointsRange = range_tag("targetDerivPointsRange", numberSubCellDims,maxSubCellsCount);
  subCellTopologyKey = key_tag("subCellTopologyKey",numberSubCellDims,maxSubCellsCount);

  numBasisEvalPoints += numVertices;
  numTargetEvalPoints += numVertices;
  
  // The set of the eval points on the reference vertex contains only point (0.0).
  // Not very useful, updating these for completeness  
  for(ordinal_type iv=0; iv<numVertices; ++iv) {
    basisPointsRange(0,iv) = range_type(iv, iv+1);
    basisCubPoints[0][iv] = host_view_type("basisCubPoints",1,1);
    targetPointsRange(0,iv) = range_type(iv, iv+1);
    targetCubPoints[0][iv] = host_view_type("targetCubPoints",1,1);
  }

  DefaultCubatureFactory cub_factory;
  for(ordinal_type ie=0; ie<numEdges; ++ie) {
    ordinal_type cub_degree = 2*edgeBasisCubDegree;
    subCellTopologyKey(edgeDim,ie) = cellBasis->getBaseCellTopology().getKey(edgeDim, ie);
    auto edgeBasisCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(cellBasis->getBaseCellTopology().getKey(edgeDim, ie), cub_degree, EPolyType::POLYTYPE_MAX, true);
    basisDerivPointsRange(edgeDim,ie) = range_type(numBasisDerivEvalPoints, numBasisDerivEvalPoints+edgeBasisCub->getNumPoints());
    numBasisDerivEvalPoints +=  edgeBasisCub->getNumPoints();
    maxNumBasisDerivEvalPoints = std::max(maxNumBasisDerivEvalPoints, edgeBasisCub->getNumPoints());
    basisDerivCubPoints[edgeDim][ie] = host_view_type("basisDerivCubPoints",edgeBasisCub->getNumPoints(),edgeDim);
    basisDerivCubWeights[edgeDim][ie] = host_view_type("basisDerivCubWeights",edgeBasisCub->getNumPoints());
    edgeBasisCub->getCubature(basisDerivCubPoints[edgeDim][ie], basisDerivCubWeights[edgeDim][ie]);

    cub_degree = edgeBasisCubDegree + targetGradCubDegree;
    auto edgeTargetCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(cellBasis->getBaseCellTopology().getKey(edgeDim, ie), cub_degree, EPolyType::POLYTYPE_MAX, true);
    targetDerivPointsRange(edgeDim,ie) = range_type(numTargetDerivEvalPoints, numTargetDerivEvalPoints+edgeTargetCub->getNumPoints());
    numTargetDerivEvalPoints +=  edgeTargetCub->getNumPoints();
    maxNumTargetDerivEvalPoints = std::max(maxNumTargetDerivEvalPoints, edgeTargetCub->getNumPoints());
    targetDerivCubPoints[edgeDim][ie] = host_view_type("targetDerivCubPoints",edgeTargetCub->getNumPoints(),edgeDim);
    targetDerivCubWeights[edgeDim][ie] = host_view_type("targetDerivCubWeights",edgeTargetCub->getNumPoints());
    edgeTargetCub->getCubature(targetDerivCubPoints[edgeDim][ie], targetDerivCubWeights[edgeDim][ie]);
  }

  for(ordinal_type iface=0; iface<numFaces; ++iface) {
    ordinal_type cub_degree = 2*faceBasisCubDegree;
    subCellTopologyKey(faceDim,iface) = cellBasis->getBaseCellTopology().getKey(faceDim, iface);
    auto faceBasisGradCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(faceDim,iface), cub_degree, EPolyType::POLYTYPE_MAX, true);
    basisDerivPointsRange(faceDim,iface) = range_type(numBasisDerivEvalPoints, numBasisDerivEvalPoints+faceBasisGradCub->getNumPoints());
    numBasisDerivEvalPoints +=  faceBasisGradCub->getNumPoints();
    maxNumBasisDerivEvalPoints = std::max(maxNumBasisDerivEvalPoints, faceBasisGradCub->getNumPoints());
    basisDerivCubPoints[faceDim][iface] = host_view_type("basisDerivCubPoints",faceBasisGradCub->getNumPoints(),faceDim);
    basisDerivCubWeights[faceDim][iface] = host_view_type("basisDerivCubWeights",faceBasisGradCub->getNumPoints());
    faceBasisGradCub->getCubature(basisDerivCubPoints[faceDim][iface], basisDerivCubWeights[faceDim][iface]);

    cub_degree = faceBasisCubDegree + targetGradCubDegree;
    auto faceTargetDerivCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(faceDim,iface), cub_degree, EPolyType::POLYTYPE_MAX, true);
    targetDerivPointsRange(faceDim,iface) = range_type(numTargetDerivEvalPoints, numTargetDerivEvalPoints+faceTargetDerivCub->getNumPoints());
    numTargetDerivEvalPoints +=  faceTargetDerivCub->getNumPoints();
    maxNumTargetDerivEvalPoints = std::max(maxNumTargetDerivEvalPoints, faceTargetDerivCub->getNumPoints());
    targetDerivCubPoints[faceDim][iface] = host_view_type("targetDerivCubPoints",faceTargetDerivCub->getNumPoints(),faceDim);
    targetDerivCubWeights[faceDim][iface] = host_view_type("targetDerivCubWeights",faceTargetDerivCub->getNumPoints());
    faceTargetDerivCub->getCubature(targetDerivCubPoints[faceDim][iface], targetDerivCubWeights[faceDim][iface]);
  }
  subCellTopologyKey(dim,0) = cellBasis->getBaseCellTopology().getBaseKey();
  if(hasCellDofs) {
    ordinal_type cub_degree = 2*basisCubDegree;
    auto elemBasisGradCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(dim,0), cub_degree, EPolyType::POLYTYPE_MAX, true);
    basisDerivPointsRange(dim,0) = range_type(numBasisDerivEvalPoints, numBasisDerivEvalPoints+elemBasisGradCub->getNumPoints());
    numBasisDerivEvalPoints +=  elemBasisGradCub->getNumPoints();
    maxNumBasisDerivEvalPoints = std::max(maxNumBasisDerivEvalPoints, elemBasisGradCub->getNumPoints());
    basisDerivCubPoints[dim][0] = host_view_type("basisDerivCubPoints",elemBasisGradCub->getNumPoints(),dim);
    basisDerivCubWeights[dim][0] = host_view_type("basisDerivCubWeights",elemBasisGradCub->getNumPoints());
    elemBasisGradCub->getCubature(basisDerivCubPoints[dim][0], basisDerivCubWeights[dim][0]);

    cub_degree = basisCubDegree + targetGradCubDegree;
    auto elemTargetGradCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(dim,0), cub_degree, EPolyType::POLYTYPE_MAX, true);
    targetDerivPointsRange(dim,0) = range_type(numTargetDerivEvalPoints, numTargetDerivEvalPoints+elemTargetGradCub->getNumPoints());
    numTargetDerivEvalPoints +=  elemTargetGradCub->getNumPoints();
    maxNumTargetDerivEvalPoints = std::max(maxNumTargetDerivEvalPoints, elemTargetGradCub->getNumPoints());
    targetDerivCubPoints[dim][0] = host_view_type("targetDerivCubPoints",elemTargetGradCub->getNumPoints(),dim);
    targetDerivCubWeights[dim][0] = host_view_type("targetDerivCubWeights",elemTargetGradCub->getNumPoints());
    elemTargetGradCub->getCubature(targetDerivCubPoints[dim][0], targetDerivCubWeights[dim][0]);
  }

  allBasisEPoints = view_type("allBasisPoints", numBasisEvalPoints, dim);
  allTargetEPoints = view_type("allTargetPoints", numTargetEvalPoints, dim);
  allBasisDerivEPoints = view_type("allBasisDerivPoints", numBasisDerivEvalPoints, dim);
  allTargetDerivEPoints = view_type("allTargetDerivPoints", numTargetDerivEvalPoints, dim);

  if(numVertices>0) {
    for(ordinal_type iv=0; iv<numVertices; ++iv) {
      CellTools<DeviceType>::getReferenceVertex(Kokkos::subview(allBasisEPoints, iv, Kokkos::ALL()), cellTopo, iv);
      CellTools<DeviceType>::getReferenceVertex(Kokkos::subview(allTargetEPoints, iv, Kokkos::ALL()), cellTopo, iv);
    }
  }

  if(numEdges>0) {
    auto subcellParamEdge = RefSubcellParametrization<DeviceType>::get(edgeDim, cellTopo.getKey());
    for(ordinal_type ie=0; ie<numEdges; ++ie) {
      auto edgeBasisDerivEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisDerivCubPoints[edgeDim][ie]);
      CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allBasisDerivEPoints,basisDerivPointsRange(edgeDim,ie),Kokkos::ALL()), edgeBasisDerivEPoints, subcellParamEdge, ie);

      auto edgeTargetDerivEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetDerivCubPoints[edgeDim][ie]);
      CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allTargetDerivEPoints,targetDerivPointsRange(edgeDim,ie),Kokkos::ALL()), edgeTargetDerivEPoints, subcellParamEdge, ie);
    }
  }

  if(numFaces>0) {
    auto subcellParamFace = RefSubcellParametrization<DeviceType>::get(faceDim, cellTopo.getKey());
    for(ordinal_type iface=0; iface<numFaces; ++iface) {
      auto faceBasisDerivEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisDerivCubPoints[faceDim][iface]);
      CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allBasisDerivEPoints,basisDerivPointsRange(faceDim,iface),Kokkos::ALL()), faceBasisDerivEPoints, subcellParamFace, iface);

      auto faceTargetDerivEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetDerivCubPoints[faceDim][iface]);
      CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allTargetDerivEPoints,targetDerivPointsRange(faceDim,iface),Kokkos::ALL()), faceTargetDerivEPoints, subcellParamFace, iface);
    }
  }

  if(hasCellDofs) {
    auto cellBasisDerivEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisDerivCubPoints[dim][0]);
    Kokkos::deep_copy(Kokkos::subview(allBasisDerivEPoints,basisDerivPointsRange(dim,0), Kokkos::ALL()), cellBasisDerivEPoints);

    auto cellTargetDerivEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetDerivCubPoints[dim][0]);
    Kokkos::deep_copy(Kokkos::subview(allTargetDerivEPoints,targetDerivPointsRange(dim,0),Kokkos::ALL()), cellTargetDerivEPoints);
  }
}

template<typename DeviceType, typename ValueType>
template<typename BasisPtrType>
void ProjectionStruct<DeviceType,ValueType>::createHCurlProjectionStruct(const BasisPtrType cellBasis,
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
  ordinal_type hasCellDofs = (cellBasis->getDofCount(dim, 0) > 0);

  maxNumBasisEvalPoints = 0; maxNumTargetEvalPoints = 0;
  maxNumBasisDerivEvalPoints = 0; maxNumTargetDerivEvalPoints = 0;

  basisPointsRange = range_tag("basisPointsRange", numberSubCellDims,maxSubCellsCount);
  targetPointsRange = range_tag("targetPointsRange", numberSubCellDims,maxSubCellsCount);
  basisDerivPointsRange = range_tag("basisDerivPointsRange", numberSubCellDims,maxSubCellsCount);
  targetDerivPointsRange = range_tag("targetDerivPointsRange", numberSubCellDims,maxSubCellsCount);
  subCellTopologyKey = key_tag("subCellTopologyKey",numberSubCellDims,maxSubCellsCount);

  DefaultCubatureFactory cub_factory;
  for(ordinal_type ie=0; ie<numEdges; ++ie) {
    ordinal_type cub_degree = 2*edgeBasisCubDegree;
    subCellTopologyKey(edgeDim,ie) = cellBasis->getBaseCellTopology().getKey(edgeDim, ie);
    auto edgeBasisCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(edgeDim,ie), cub_degree, EPolyType::POLYTYPE_MAX, true);
    basisPointsRange(edgeDim,ie) = range_type(numBasisEvalPoints, numBasisEvalPoints+edgeBasisCub->getNumPoints());
    numBasisEvalPoints +=  edgeBasisCub->getNumPoints();
    maxNumBasisEvalPoints = std::max(maxNumBasisEvalPoints, edgeBasisCub->getNumPoints());
    basisCubPoints[edgeDim][ie] = host_view_type("basisCubPoints",edgeBasisCub->getNumPoints(),edgeDim);
    basisCubWeights[edgeDim][ie] = host_view_type("basisCubWeights",edgeBasisCub->getNumPoints());
    edgeBasisCub->getCubature(basisCubPoints[edgeDim][ie], basisCubWeights[edgeDim][ie]);

    cub_degree = edgeBasisCubDegree + targetCubDegree;
    auto edgeTargetCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(edgeDim,ie), cub_degree, EPolyType::POLYTYPE_MAX, true);
    targetPointsRange(edgeDim,ie) = range_type(numTargetEvalPoints, numTargetEvalPoints+edgeTargetCub->getNumPoints());
    numTargetEvalPoints +=  edgeTargetCub->getNumPoints();
    maxNumTargetEvalPoints = std::max(maxNumTargetEvalPoints, edgeTargetCub->getNumPoints());
    targetCubPoints[edgeDim][ie] = host_view_type("targetCubPoints",edgeTargetCub->getNumPoints(),edgeDim);
    targetCubWeights[edgeDim][ie] = host_view_type("targetCubWeights",edgeTargetCub->getNumPoints());
    edgeTargetCub->getCubature(targetCubPoints[edgeDim][ie], targetCubWeights[edgeDim][ie]);
  }

  for(ordinal_type iface=0; iface<numFaces; ++iface) {
    ordinal_type cub_degree = 2*faceBasisCubDegree;
    subCellTopologyKey(faceDim,iface) = cellBasis->getBaseCellTopology().getKey(faceDim, iface);
    auto faceBasisCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(faceDim,iface), cub_degree, EPolyType::POLYTYPE_MAX, true);
    basisPointsRange(faceDim,iface) = range_type(numBasisEvalPoints, numBasisEvalPoints+faceBasisCub->getNumPoints());
    numBasisEvalPoints +=  faceBasisCub->getNumPoints();
    maxNumBasisEvalPoints = std::max(maxNumBasisEvalPoints, faceBasisCub->getNumPoints());
    basisCubPoints[faceDim][iface] = host_view_type("basisCubPoints",faceBasisCub->getNumPoints(),faceDim);
    basisCubWeights[faceDim][iface] = host_view_type("basisCubWeights",faceBasisCub->getNumPoints());
    faceBasisCub->getCubature(basisCubPoints[faceDim][iface], basisCubWeights[faceDim][iface]);

    auto faceBasisDerivCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(faceDim,iface), cub_degree, EPolyType::POLYTYPE_MAX, true);
    basisDerivPointsRange(faceDim,iface) = range_type(numBasisDerivEvalPoints, numBasisDerivEvalPoints+faceBasisCub->getNumPoints());
    numBasisDerivEvalPoints +=  faceBasisCub->getNumPoints();
    maxNumBasisDerivEvalPoints = std::max(maxNumBasisDerivEvalPoints, faceBasisCub->getNumPoints());
    basisDerivCubPoints[faceDim][iface] = host_view_type("basisDerivCubPoints",faceBasisCub->getNumPoints(),faceDim);
    basisDerivCubWeights[faceDim][iface] = host_view_type("basisDerivCubWeights",faceBasisCub->getNumPoints());
    faceBasisCub->getCubature(basisDerivCubPoints[faceDim][iface], basisDerivCubWeights[faceDim][iface]);

    cub_degree = faceBasisCubDegree + targetCubDegree;
    auto faceTargetCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(faceDim,iface), cub_degree, EPolyType::POLYTYPE_MAX, true);
    targetPointsRange(faceDim,iface) = range_type(numTargetEvalPoints, numTargetEvalPoints+faceTargetCub->getNumPoints());
    numTargetEvalPoints +=  faceTargetCub->getNumPoints();
    maxNumTargetEvalPoints = std::max(maxNumTargetEvalPoints, faceTargetCub->getNumPoints());
    targetCubPoints[faceDim][iface] = host_view_type("targetCubPoints",faceTargetCub->getNumPoints(),faceDim);
    targetCubWeights[faceDim][iface] = host_view_type("targetCubWeights",faceTargetCub->getNumPoints());
    faceTargetCub->getCubature(targetCubPoints[faceDim][iface], targetCubWeights[faceDim][iface]);

    cub_degree = faceBasisCubDegree + targetCurlCubDegre;
    auto faceTargetDerivCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(faceDim,iface), cub_degree, EPolyType::POLYTYPE_MAX, true);
    targetDerivPointsRange(faceDim,iface) = range_type(numTargetDerivEvalPoints, numTargetDerivEvalPoints+faceTargetDerivCub->getNumPoints());
    numTargetDerivEvalPoints +=  faceTargetDerivCub->getNumPoints();
    maxNumTargetDerivEvalPoints = std::max(maxNumTargetDerivEvalPoints, faceTargetDerivCub->getNumPoints());
    targetDerivCubPoints[faceDim][iface] = host_view_type("targetDerivCubPoints",faceTargetDerivCub->getNumPoints(),faceDim);
    targetDerivCubWeights[faceDim][iface] = host_view_type("targetDerivCubWeights",faceTargetDerivCub->getNumPoints());
    faceTargetDerivCub->getCubature(targetDerivCubPoints[faceDim][iface], targetDerivCubWeights[faceDim][iface]);
  }

  subCellTopologyKey(dim,0) = cellBasis->getBaseCellTopology().getBaseKey();
  if(hasCellDofs) {
    ordinal_type cub_degree = 2*basisCubDegree;
    auto elemBasisCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(dim,0), cub_degree, EPolyType::POLYTYPE_MAX, true);
    basisPointsRange(dim,0) = range_type(numBasisEvalPoints, numBasisEvalPoints+elemBasisCub->getNumPoints());
    numBasisEvalPoints +=  elemBasisCub->getNumPoints();
    maxNumBasisEvalPoints = std::max(maxNumBasisEvalPoints, elemBasisCub->getNumPoints());
    basisCubPoints[dim][0] = host_view_type("basisCubPoints",elemBasisCub->getNumPoints(),dim);
    basisCubWeights[dim][0] = host_view_type("basisCubWeights",elemBasisCub->getNumPoints());
    elemBasisCub->getCubature(basisCubPoints[dim][0], basisCubWeights[dim][0]);

    basisDerivPointsRange(dim,0) = range_type(numBasisDerivEvalPoints, numBasisDerivEvalPoints+elemBasisCub->getNumPoints());
    numBasisDerivEvalPoints +=  elemBasisCub->getNumPoints();
    maxNumBasisDerivEvalPoints = std::max(maxNumBasisDerivEvalPoints, elemBasisCub->getNumPoints());
    basisDerivCubPoints[dim][0] = host_view_type("basisDerivCubPoints",elemBasisCub->getNumPoints(),dim);
    basisDerivCubWeights[dim][0] = host_view_type("basisDerivCubWeights",elemBasisCub->getNumPoints());
    elemBasisCub->getCubature(basisDerivCubPoints[dim][0], basisDerivCubWeights[dim][0]);

    cub_degree = basisCubDegree + targetCubDegree;
    auto elemTargetCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(dim,0), cub_degree, EPolyType::POLYTYPE_MAX, true);
    targetPointsRange(dim,0) = range_type(numTargetEvalPoints, numTargetEvalPoints+elemTargetCub->getNumPoints());
    numTargetEvalPoints +=  elemTargetCub->getNumPoints();
    maxNumTargetEvalPoints = std::max(maxNumTargetEvalPoints, elemTargetCub->getNumPoints());
    targetCubPoints[dim][0] = host_view_type("targetCubPoints",elemTargetCub->getNumPoints(),dim);
    targetCubWeights[dim][0] = host_view_type("targetCubWeights",elemTargetCub->getNumPoints());
    elemTargetCub->getCubature(targetCubPoints[dim][0], targetCubWeights[dim][0]);

    cub_degree = basisCubDegree + targetCurlCubDegre;
    auto elemTargetCurlCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(dim,0), cub_degree, EPolyType::POLYTYPE_MAX, true);
    targetDerivPointsRange(dim,0) = range_type(numTargetDerivEvalPoints, numTargetDerivEvalPoints+elemTargetCurlCub->getNumPoints());
    numTargetDerivEvalPoints +=  elemTargetCurlCub->getNumPoints();
    maxNumTargetDerivEvalPoints = std::max(maxNumTargetDerivEvalPoints, elemTargetCurlCub->getNumPoints());
    targetDerivCubPoints[dim][0] = host_view_type("targetDerivCubPoints",elemTargetCurlCub->getNumPoints(),dim);
    targetDerivCubWeights[dim][0] = host_view_type("targetDerivCubWeights",elemTargetCurlCub->getNumPoints());
    elemTargetCurlCub->getCubature(targetDerivCubPoints[dim][0], targetDerivCubWeights[dim][0]);
  }

  allBasisEPoints = view_type("allBasisPoints", numBasisEvalPoints, dim);
  allTargetEPoints = view_type("allTargetPoints", numTargetEvalPoints, dim);
  allBasisDerivEPoints = view_type("allBasisDerivPoints", numBasisDerivEvalPoints, dim);
  allTargetDerivEPoints = view_type("allTargetDerivPoints", numTargetDerivEvalPoints, dim);

  auto subcellParamEdge = RefSubcellParametrization<DeviceType>::get(edgeDim, cellTopo.getKey());
  for(ordinal_type ie=0; ie<numEdges; ++ie) {
    auto edgeBasisEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisCubPoints[edgeDim][ie]);
    CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allBasisEPoints,basisPointsRange(edgeDim,ie),Kokkos::ALL()), edgeBasisEPoints, subcellParamEdge, ie);

    auto edgeTargetEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetCubPoints[edgeDim][ie]);
    CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allTargetEPoints,targetPointsRange(edgeDim,ie),Kokkos::ALL()), edgeTargetEPoints, subcellParamEdge, ie);
  }

  if(numFaces>0) {
    auto subcellParamFace = RefSubcellParametrization<DeviceType>::get(faceDim, cellTopo.getKey());
    for(ordinal_type iface=0; iface<numFaces; ++iface) {
      auto faceBasisEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisCubPoints[faceDim][iface]);
      CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allBasisEPoints,basisPointsRange(faceDim,iface),Kokkos::ALL()), faceBasisEPoints, subcellParamFace, iface);
      
      auto faceBasisDerivEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisDerivCubPoints[faceDim][iface]);
      CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allBasisDerivEPoints,basisDerivPointsRange(faceDim,iface),Kokkos::ALL()), faceBasisDerivEPoints, subcellParamFace, iface);
      
      auto faceTargetEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetCubPoints[faceDim][iface]);
      CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allTargetEPoints,targetPointsRange(faceDim,iface),Kokkos::ALL()), faceTargetEPoints, subcellParamFace, iface);
      
      auto faceTargetDerivEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetDerivCubPoints[faceDim][iface]);
      CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allTargetDerivEPoints,targetDerivPointsRange(faceDim,iface),Kokkos::ALL()), faceTargetDerivEPoints, subcellParamFace, iface);
    }
  }

  if(hasCellDofs) {
    auto cellBasisEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisCubPoints[dim][0]);
    Kokkos::deep_copy(Kokkos::subview(allBasisEPoints,basisPointsRange(dim,0), Kokkos::ALL()), cellBasisEPoints);
    
    auto cellBasisDerivEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisDerivCubPoints[dim][0]);
    Kokkos::deep_copy(Kokkos::subview(allBasisDerivEPoints,basisDerivPointsRange(dim,0), Kokkos::ALL()), cellBasisDerivEPoints);

    auto cellTargetEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetCubPoints[dim][0]);
    Kokkos::deep_copy(Kokkos::subview(allTargetEPoints,targetPointsRange(dim,0),Kokkos::ALL()), cellTargetEPoints);
    
    auto cellTargetDerivEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetDerivCubPoints[dim][0]);
    Kokkos::deep_copy(Kokkos::subview(allTargetDerivEPoints,targetDerivPointsRange(dim,0),Kokkos::ALL()), cellTargetDerivEPoints);
  }
}

template<typename DeviceType, typename ValueType>
template<typename BasisPtrType>
void ProjectionStruct<DeviceType,ValueType>::createHDivProjectionStruct(const BasisPtrType cellBasis,
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
  ordinal_type hasCellDofs = (cellBasis->getDofCount(dim, 0) > 0);

  INTREPID2_TEST_FOR_ABORT( numSides > maxSubCellsCount,
      ">>> ERROR (Intrepid2::ProjectionStruct:createHDivProjectionStruct, Projections do not support a cell topology with so many sides");

  maxNumBasisEvalPoints = 0; maxNumTargetEvalPoints = 0;
  maxNumBasisDerivEvalPoints = 0; maxNumTargetDerivEvalPoints = 0;

  basisPointsRange = range_tag("basisPointsRange", numberSubCellDims,maxSubCellsCount);
  targetPointsRange = range_tag("targetPointsRange", numberSubCellDims,maxSubCellsCount);
  basisDerivPointsRange = range_tag("basisDerivPointsRange", numberSubCellDims,maxSubCellsCount);
  targetDerivPointsRange = range_tag("targetDerivPointsRange", numberSubCellDims,maxSubCellsCount);
  subCellTopologyKey = key_tag("subCellTopologyKey",numberSubCellDims,maxSubCellsCount);

  Basis<HostDeviceType,ValueType,ValueType> *hcurlBasis = NULL;
  if(cellTopo.getKey() == shards::getCellTopologyData<shards::Hexahedron<8> >()->key)
    hcurlBasis = new Basis_HCURL_HEX_In_FEM<HostDeviceType,ValueType,ValueType>(cellBasis->getDegree());
  else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key)
    hcurlBasis = new Basis_HCURL_TET_In_FEM<HostDeviceType,ValueType,ValueType>(cellBasis->getDegree());
  else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Wedge<6> >()->key)
    hcurlBasis = new typename DerivedNodalBasisFamily<HostDeviceType,ValueType,ValueType>::HCURL_WEDGE(cellBasis->getDegree());
// TODO: uncomment the next two lines once H(curl) pyramid implemented
//  else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Pyramid<5> >()->key)
//    hcurlBasis = new typename DerivedNodalBasisFamily<DeviceType,ValueType,ValueType>::HCURL_PYR(cellBasis->getDegree());
  else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Quadrilateral<4> >()->key)
    hcurlBasis = new Basis_HGRAD_QUAD_Cn_FEM<HostDeviceType,ValueType,ValueType>(cellBasis->getDegree());
  else if(cellTopo.getKey() == shards::getCellTopologyData<shards::Triangle<3> >()->key)
    hcurlBasis = new Basis_HGRAD_TRI_Cn_FEM<HostDeviceType,ValueType,ValueType>(cellBasis->getDegree());
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
    subCellTopologyKey(sideDim,is) = cellBasis->getBaseCellTopology().getKey(sideDim, is);
    auto sideBasisCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(sideDim,is), cub_degree, EPolyType::POLYTYPE_MAX, true);
    basisPointsRange(sideDim,is) = range_type(numBasisEvalPoints, numBasisEvalPoints+sideBasisCub->getNumPoints());
    numBasisEvalPoints +=  sideBasisCub->getNumPoints();
    basisCubPoints[sideDim][is] = host_view_type("basisCubPoints",sideBasisCub->getNumPoints(),sideDim);
    basisCubWeights[sideDim][is] = host_view_type("basisCubWeights",sideBasisCub->getNumPoints());
    sideBasisCub->getCubature(basisCubPoints[sideDim][is], basisCubWeights[sideDim][is]);
    maxNumBasisEvalPoints = std::max(maxNumBasisEvalPoints, sideBasisCub->getNumPoints());

    cub_degree = sideBasisCubDegree + targetCubDegree;
    auto sideTargetCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(sideDim,is), cub_degree, EPolyType::POLYTYPE_MAX, true);
    targetPointsRange(sideDim,is) = range_type(numTargetEvalPoints, numTargetEvalPoints+sideTargetCub->getNumPoints());
    numTargetEvalPoints +=  sideTargetCub->getNumPoints();
    targetCubPoints[sideDim][is] = host_view_type("targetCubPoints",sideTargetCub->getNumPoints(),sideDim);
    targetCubWeights[sideDim][is] = host_view_type("targetCubWeights",sideTargetCub->getNumPoints());
    sideTargetCub->getCubature(targetCubPoints[sideDim][is], targetCubWeights[sideDim][is]);
    maxNumTargetEvalPoints = std::max(maxNumTargetEvalPoints, sideTargetCub->getNumPoints());
  }

  subCellTopologyKey(dim,0) = cellBasis->getBaseCellTopology().getBaseKey();
  if(hasCellDofs) {
    ordinal_type cub_degree = 2*basisCubDegree - 1;
    auto elemBasisDivCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(dim,0), cub_degree, EPolyType::POLYTYPE_MAX, true);
    basisDerivPointsRange(dim,0) = range_type(numBasisDerivEvalPoints, numBasisDerivEvalPoints+elemBasisDivCub->getNumPoints());
    numBasisDerivEvalPoints +=  elemBasisDivCub->getNumPoints();
    basisDerivCubPoints[dim][0] = host_view_type("basisDerivCubPoints",elemBasisDivCub->getNumPoints(),dim);
    basisDerivCubWeights[dim][0] = host_view_type("basisDerivCubWeights",elemBasisDivCub->getNumPoints());
    elemBasisDivCub->getCubature(basisDerivCubPoints[dim][0], basisDerivCubWeights[dim][0]);
    maxNumBasisDerivEvalPoints = std::max(maxNumBasisDerivEvalPoints, elemBasisDivCub->getNumPoints());

    cub_degree = basisCubDegree - 1 + targetDivCubDegre;
    auto elemTargetDivCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(dim,0), cub_degree, EPolyType::POLYTYPE_MAX, true);
    targetDerivPointsRange(dim,0) = range_type(numTargetDerivEvalPoints, numTargetDerivEvalPoints+elemTargetDivCub->getNumPoints());
    numTargetDerivEvalPoints +=  elemTargetDivCub->getNumPoints();
    targetDerivCubPoints[dim][0] = host_view_type("targetDerivCubPoints",elemTargetDivCub->getNumPoints(),dim);
    targetDerivCubWeights[dim][0] = host_view_type("targetDerivCubWeights",elemTargetDivCub->getNumPoints());
    elemTargetDivCub->getCubature(targetDerivCubPoints[dim][0], targetDerivCubWeights[dim][0]);
    maxNumTargetDerivEvalPoints = std::max(maxNumTargetDerivEvalPoints, elemTargetDivCub->getNumPoints());

    if(haveHCurlConstraint)
    {
      cub_degree = 2*basisCubDegree;
      auto elemBasisCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(dim,0), cub_degree, EPolyType::POLYTYPE_MAX, true);
      basisPointsRange(dim,0) = range_type(numBasisEvalPoints, numBasisEvalPoints + elemBasisCub->getNumPoints());
      numBasisEvalPoints +=  elemBasisCub->getNumPoints();
      basisCubPoints[dim][0] = host_view_type("basisCubPoints",elemBasisCub->getNumPoints(),dim);
      basisCubWeights[dim][0] = host_view_type("basisCubWeights",elemBasisCub->getNumPoints());
      elemBasisCub->getCubature(basisCubPoints[dim][0], basisCubWeights[dim][0]);
      maxNumBasisEvalPoints = std::max(maxNumBasisEvalPoints, elemBasisCub->getNumPoints());

      cub_degree = basisCubDegree + targetCubDegree;
      auto elemTargetCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(subCellTopologyKey(dim,0), cub_degree, EPolyType::POLYTYPE_MAX, true);
      targetPointsRange(dim,0) = range_type(numTargetEvalPoints, numTargetEvalPoints + elemTargetCub->getNumPoints());
      numTargetEvalPoints +=  elemTargetCub->getNumPoints();
      targetCubPoints[dim][0] = host_view_type("targetCubPoints",elemTargetCub->getNumPoints(),dim);
      targetCubWeights[dim][0] = host_view_type("targetCubWeights",elemTargetCub->getNumPoints());
      elemTargetCub->getCubature(targetCubPoints[dim][0], targetCubWeights[dim][0]);
      maxNumTargetEvalPoints = std::max(maxNumTargetEvalPoints, elemTargetCub->getNumPoints());
    }
  }

  allBasisEPoints = view_type("allBasisPoints", numBasisEvalPoints, dim);
  allTargetEPoints = view_type("allTargetPoints", numTargetEvalPoints, dim);
  allBasisDerivEPoints = view_type("allBasisDerivPoints", numBasisDerivEvalPoints, dim);
  allTargetDerivEPoints = view_type("allTargetDerivPoints", numTargetDerivEvalPoints, dim);


  auto subcellParamSide = RefSubcellParametrization<DeviceType>::get(sideDim, cellTopo.getKey());
  for(ordinal_type is=0; is<numSides; ++is) {
    auto sideBasisEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisCubPoints[sideDim][is]);
    CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allBasisEPoints,basisPointsRange(sideDim,is),Kokkos::ALL()), sideBasisEPoints, subcellParamSide, is);

    auto sideTargetEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetCubPoints[sideDim][is]);
    CellTools<DeviceType>::mapToReferenceSubcell(Kokkos::subview(allTargetEPoints,targetPointsRange(sideDim,is),Kokkos::ALL()), sideTargetEPoints, subcellParamSide, is);
  }

  if(hasCellDofs) {
    if(haveHCurlConstraint) {
      auto cellBasisEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisCubPoints[dim][0]);
      Kokkos::deep_copy(Kokkos::subview(allBasisEPoints,basisPointsRange(dim,0), Kokkos::ALL()), cellBasisEPoints);

      auto cellTargetEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetCubPoints[dim][0]);
      Kokkos::deep_copy(Kokkos::subview(allTargetEPoints,targetPointsRange(dim,0),Kokkos::ALL()), cellTargetEPoints);
    }
    
    auto cellBasisDerivEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisDerivCubPoints[dim][0]);
    Kokkos::deep_copy(Kokkos::subview(allBasisDerivEPoints,basisDerivPointsRange(dim,0), Kokkos::ALL()), cellBasisDerivEPoints);

    auto cellTargetDerivEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetDerivCubPoints[dim][0]);
    Kokkos::deep_copy(Kokkos::subview(allTargetDerivEPoints,targetDerivPointsRange(dim,0),Kokkos::ALL()), cellTargetDerivEPoints);
  }
}

template<typename DeviceType, typename ValueType>
template<typename BasisPtrType>
void ProjectionStruct<DeviceType,ValueType>::createHVolProjectionStruct(const BasisPtrType cellBasis,
    const ordinal_type targetCubDegree) {
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  numBasisEvalPoints = 0;
  numBasisDerivEvalPoints = 0;
  numTargetEvalPoints = 0;
  numTargetDerivEvalPoints = 0;

  basisPointsRange = range_tag("basisPointsRange", 4,maxSubCellsCount);
  targetPointsRange = range_tag("targetPointsRange", 4,maxSubCellsCount);
  basisDerivPointsRange = range_tag("basisDerivPointsRange", 4,maxSubCellsCount);
  targetDerivPointsRange = range_tag("targetDerivPointsRange", 4,maxSubCellsCount);
  subCellTopologyKey = key_tag("subCellTopologyKey",4,maxSubCellsCount);

  ordinal_type basisCubDegree = cellBasis->getDegree();
  DefaultCubatureFactory cub_factory;
  subCellTopologyKey(dim,0) = cellBasis->getBaseCellTopology().getBaseKey();

  maxNumBasisEvalPoints = 0; maxNumTargetEvalPoints =0;

  ordinal_type cub_degree = 2*basisCubDegree;
  auto elemBasisCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(cellTopo.getBaseKey(), cub_degree, EPolyType::POLYTYPE_MAX, true);
  basisPointsRange(dim,0) = range_type(0, elemBasisCub->getNumPoints());
  numBasisEvalPoints +=  elemBasisCub->getNumPoints();
  maxNumBasisEvalPoints = elemBasisCub->getNumPoints();
  basisCubPoints[dim][0] = host_view_type("basisCubPoints",elemBasisCub->getNumPoints(),dim);
  basisCubWeights[dim][0] = host_view_type("basisCubWeights",elemBasisCub->getNumPoints());
  elemBasisCub->getCubature(basisCubPoints[dim][0], basisCubWeights[dim][0]);

  cub_degree = basisCubDegree + targetCubDegree;
  auto elemTargetCub = cub_factory.create<HostDeviceType, ValueType, ValueType>(cellTopo.getBaseKey(), cub_degree, EPolyType::POLYTYPE_MAX, true);
  targetPointsRange(dim,0) = range_type(0, elemTargetCub->getNumPoints());
  numTargetEvalPoints +=  elemTargetCub->getNumPoints();
  maxNumTargetEvalPoints = elemTargetCub->getNumPoints();
  targetCubPoints[dim][0] = host_view_type("targetCubPoints",elemTargetCub->getNumPoints(),dim);
  targetCubWeights[dim][0] = host_view_type("targetCubWeights",elemTargetCub->getNumPoints());
  elemTargetCub->getCubature(targetCubPoints[dim][0], targetCubWeights[dim][0]);

  allBasisEPoints = view_type("allBasisPoints", numBasisEvalPoints, dim);
  auto cellBasisEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),basisCubPoints[dim][0]);
  Kokkos::deep_copy(Kokkos::subview(allBasisEPoints,basisPointsRange(dim,0), Kokkos::ALL()), cellBasisEPoints);
  
  allTargetEPoints = view_type("allTargetPoints", numTargetEvalPoints, dim);
  auto cellTargetEPoints = Kokkos::create_mirror_view_and_copy(MemSpaceType(),targetCubPoints[dim][0]);
  Kokkos::deep_copy(Kokkos::subview(allTargetEPoints,targetPointsRange(dim,0),Kokkos::ALL()), cellTargetEPoints);

  allBasisDerivEPoints = view_type("allBasisDerivPoints", numBasisDerivEvalPoints, dim);
  allTargetDerivEPoints = view_type("allTargetDerivPoints", numTargetDerivEvalPoints, dim);
}

}  // Intrepid2 namespace
#endif





