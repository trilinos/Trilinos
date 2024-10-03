// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_LagrangianInterpolationDef.hpp
    \brief  Header file for the Intrepid2::LagrangianInterpolation containing definitions.
    \author Created by Mauro Perego
 */

#ifndef __INTREPID2_LAGRANGIANINTERPOLATIONDEF_HPP__
#define __INTREPID2_LAGRANGIANINTERPOLATIONDEF_HPP__

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


namespace Intrepid2 {


namespace FunctorsLagrangianTools {
template<typename CoordsViewType,
typename ortViewType,
typename t2oViewType,
typename subcellParamViewType,
typename intViewType,
typename ScalarViewType>
struct computeDofCoords {
  typedef typename ScalarViewType::value_type value_type;

  CoordsViewType dofCoords_;
  const ortViewType orts_;
  const t2oViewType tagToOrdinal_;
  const subcellParamViewType edgeParam_, faceParam_;
  const intViewType edgesInternalDofOrdinals_, facesInternalDofOrdinals_;
  const ScalarViewType edgesInternalDofCoords_;
  const ScalarViewType facesInternalDofCoords_;
  ScalarViewType edgeWorkView_, faceWorkView_;
  const ordinal_type cellDim_, numEdges_, numFaces_;
  const intViewType edgeTopoKey_, numEdgesInternalDofs_;
  const intViewType faceTopoKey_, numFacesInternalDofs_;

  computeDofCoords( CoordsViewType dofCoords,
      const ortViewType orts,
      const t2oViewType tagToOrdinal,
      const subcellParamViewType edgeParam,
      const subcellParamViewType faceParam,
      const ScalarViewType edgesInternalDofCoords,
      const ScalarViewType facesInternalDofCoords,
      const ordinal_type cellDim,
      const ordinal_type numEdges,
      const ordinal_type numFaces,
      const intViewType edgeTopoKey,
      const intViewType numEdgesInternalDofs,
      const intViewType faceTopoKey,
      const intViewType numFacesInternalDofs
  )
  : dofCoords_(dofCoords),
    orts_(orts),
    tagToOrdinal_(tagToOrdinal),
    edgeParam_(edgeParam),
    faceParam_(faceParam),
    edgesInternalDofCoords_(edgesInternalDofCoords),
    facesInternalDofCoords_(facesInternalDofCoords),
    cellDim_(cellDim),
    numEdges_(numEdges),
    numFaces_(numFaces),
    edgeTopoKey_(edgeTopoKey),
    numEdgesInternalDofs_(numEdgesInternalDofs),
    faceTopoKey_(faceTopoKey),
    numFacesInternalDofs_(numFacesInternalDofs)
  {
    if(numEdges > 0)
      edgeWorkView_ = ScalarViewType("edgeWorkView", dofCoords.extent(0), edgesInternalDofCoords.extent(1), cellDim);
    if(numFaces > 0)
      faceWorkView_ = ScalarViewType("faceWorkView", dofCoords.extent(0), facesInternalDofCoords.extent(1), cellDim);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type cell) const {
    typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;


    if(numEdges_ > 0) {
      //compute coordinates associated to edge DoFs
      ordinal_type eOrt[12];
      orts_(cell).getEdgeOrientation(eOrt, numEdges_);
      for (ordinal_type iedge=0; iedge < numEdges_; ++iedge) {
        ordinal_type numInternalDofs = numEdgesInternalDofs_(iedge);
        auto dofRange = range_type(0, numInternalDofs);
        auto edgeInternalDofCoords = Kokkos::subview(edgesInternalDofCoords_, iedge, dofRange, Kokkos::ALL());
        auto cellDofCoordsOriented = Kokkos::subview(edgeWorkView_,cell, dofRange, range_type(0,cellDim_));

        //map edge DoFs coords into parent element coords
        Impl::OrientationTools::mapSubcellCoordsToRefCell(cellDofCoordsOriented, edgeInternalDofCoords, edgeParam_, edgeTopoKey_(iedge), iedge, eOrt[iedge]);

        for (ordinal_type j=0;j<numInternalDofs;++j) {
          auto idof = tagToOrdinal_(1, iedge, j);
          for (ordinal_type d=0;d<cellDim_;++d)
            dofCoords_(cell,idof,d) = cellDofCoordsOriented(j,d);
        }
      }
    }

    if(numFaces_ > 0) {
      //compute coordinates associated to face DoFs
      ordinal_type fOrt[12];
      orts_(cell).getFaceOrientation(fOrt, numFaces_);
      //map face dofs coords into parent element coords
      for (ordinal_type iface=0; iface < numFaces_; ++iface) {
        ordinal_type ort = fOrt[iface];
        ordinal_type numInternalDofs = numFacesInternalDofs_(iface);
        auto dofRange = range_type(0, numInternalDofs);
        auto faceInternalDofCoords = Kokkos::subview(facesInternalDofCoords_, iface, dofRange, Kokkos::ALL());
        auto cellDofCoordsOriented = Kokkos::subview(faceWorkView_,cell, dofRange, range_type(0,cellDim_));

        Impl::OrientationTools::mapSubcellCoordsToRefCell(cellDofCoordsOriented, faceInternalDofCoords, faceParam_, faceTopoKey_(iface), iface, ort);
        for (ordinal_type j=0;j<numInternalDofs;++j) {
          auto idof = tagToOrdinal_(2, iface, j);
          for (ordinal_type d=0;d<cellDim_;++d)
            dofCoords_(cell,idof,d) = cellDofCoordsOriented(j,d);
        }
      }
    }
  }
};
}  // FunctorsLagrangianTools namespace


template<typename DeviceType>
template<typename BasisType,
class ...coordsProperties,
typename ortValueType, class ...ortProperties>
void
LagrangianTools<DeviceType>::getOrientedDofCoords(
    Kokkos::DynRankView<typename BasisType::scalarType, coordsProperties...> dofCoords,
    const BasisType* basis,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts) {

  using HostSpaceType = Kokkos::DefaultHostExecutionSpace;
  using scalarType = typename BasisType::scalarType;
  using ScalarViewType = Kokkos::DynRankView<scalarType, DeviceType>;
  using ScalarViewTypeHost = Kokkos::DynRankView<scalarType, HostSpaceType>;
  using intViewType = Kokkos::DynRankView<ordinal_type, DeviceType>;

  const auto topo = basis->getBaseCellTopology();
  const std::string name(basis->getName());

  ordinal_type numEdges = (basis->getDofCount(1, 0) > 0) ? topo.getEdgeCount() : 0;
  ordinal_type numFaces = (basis->getDofCount(2, 0) > 0) ? topo.getFaceCount() : 0;

  std::vector<Teuchos::RCP<Basis<DeviceType,scalarType,scalarType> > > edgeBases, faceBases;

  for(int i=0;i<numEdges;++i)
    edgeBases.push_back(basis->getSubCellRefBasis(1,i));
  for(int i=0;i<numFaces;++i)
    faceBases.push_back(basis->getSubCellRefBasis(2,i));

  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(typename DeviceType::memory_space(), basis->getAllDofOrdinal());

  const ordinal_type dim = topo.getDimension();

  ScalarViewType refDofCoords("refDofCoords", dofCoords.extent(1), dofCoords.extent(2));
  basis->getDofCoords(refDofCoords);
  RealSpaceTools<DeviceType>::clone(dofCoords,refDofCoords);

  if((numFaces == 0) && (numEdges == 0)) 
    return;

  //*** Pre-compute needed quantities related to edge DoFs that do not depend on the cell ***
  intViewType edgeTopoKey("edgeTopoKey",numEdges);
  intViewType sOrt("eOrt", numEdges);
  intViewType numEdgesInternalDofs("numEdgesInternalDofs", numEdges);
  ScalarViewType  edgesInternalDofCoords;
  intViewType  edgesInternalDofOrdinals;

  ordinal_type maxNumEdgesInternalDofs=0;
  ordinal_type edgeBasisMaxCardinality=0;

  auto hostNumEdgesInternalDofs = Kokkos::create_mirror_view(numEdgesInternalDofs);
  for (ordinal_type iedge=0; iedge < numEdges; ++iedge) {
    ordinal_type numInternalDofs = edgeBases[iedge]->getDofCount(1,0);
    hostNumEdgesInternalDofs(iedge) = numInternalDofs;
    maxNumEdgesInternalDofs = std::max(maxNumEdgesInternalDofs,numInternalDofs);
    ordinal_type edgeBasisCardinality = edgeBases[iedge]->getCardinality();
    edgeBasisMaxCardinality = std::max(edgeBasisMaxCardinality, edgeBasisCardinality);
  }
  Kokkos::deep_copy(numEdgesInternalDofs,hostNumEdgesInternalDofs);

  edgesInternalDofCoords = ScalarViewType("edgeInternalDofCoords", numEdges, maxNumEdgesInternalDofs,1);
  edgesInternalDofOrdinals = intViewType("edgeInternalDofCoords", numEdges, maxNumEdgesInternalDofs);
  auto hostEdgesInternalDofCoords = Kokkos::create_mirror_view(edgesInternalDofCoords);
  auto hostEdgesInternalDofOrdinals = Kokkos::create_mirror_view(edgesInternalDofOrdinals);
  auto hostEdgeTopoKey = Kokkos::create_mirror_view(edgeTopoKey);
  for (ordinal_type iedge=0; iedge < numEdges; ++iedge) {
    auto hostEdgeBasisPtr = edgeBases[iedge]->getHostBasis();
    hostEdgeTopoKey(iedge) = hostEdgeBasisPtr->getBaseCellTopology().getBaseKey();
    ordinal_type edgeBasisCardinality = hostEdgeBasisPtr->getCardinality();
    ScalarViewTypeHost  edgeDofCoords("edgeDofCoords", edgeBasisCardinality, 1);
    hostEdgeBasisPtr->getDofCoords(edgeDofCoords);
    for(ordinal_type i=0; i<hostNumEdgesInternalDofs(iedge); ++i) {
      hostEdgesInternalDofOrdinals(iedge, i) = hostEdgeBasisPtr->getDofOrdinal(1, 0, i);
      hostEdgesInternalDofCoords(iedge, i,0) = edgeDofCoords(hostEdgesInternalDofOrdinals(iedge, i),0);
    }
  }
  Kokkos::deep_copy(edgesInternalDofCoords,hostEdgesInternalDofCoords);
  Kokkos::deep_copy(edgeTopoKey,hostEdgeTopoKey);

  auto edgeParam = RefSubcellParametrization<DeviceType>::get(1, topo.getKey());

  //*** Pre-compute needed quantities related to face DoFs that do not depend on the cell ***
  intViewType faceTopoKey("faceTopoKey",numFaces);
  intViewType fOrt("fOrt", numFaces);
  intViewType numFacesInternalDofs("numFacesInternalDofs", numFaces);
  ScalarViewType  facesInternalDofCoords;
  intViewType  facesInternalDofOrdinals;

  ordinal_type maxNumFacesInternalDofs=0;
  ordinal_type faceBasisMaxCardinality=0;

  auto hostNumFacesInternalDofs = Kokkos::create_mirror_view(numFacesInternalDofs);
  for (ordinal_type iface=0; iface < numFaces; ++iface) {
    ordinal_type numInternalDofs = faceBases[iface]->getDofCount(2,0);
    hostNumFacesInternalDofs(iface) = numInternalDofs;
    maxNumFacesInternalDofs = std::max(maxNumFacesInternalDofs,numInternalDofs);
    ordinal_type faceBasisCardinality = faceBases[iface]->getCardinality();
    faceBasisMaxCardinality = std::max(faceBasisMaxCardinality, faceBasisCardinality);
  }
  Kokkos::deep_copy(numFacesInternalDofs,hostNumFacesInternalDofs);

  facesInternalDofCoords = ScalarViewType("faceInternalDofCoords", numFaces, maxNumFacesInternalDofs, 2);
  facesInternalDofOrdinals = intViewType("faceInternalDofCoords", numFaces, maxNumFacesInternalDofs);

  auto hostFacesInternalDofCoords = Kokkos::create_mirror_view(facesInternalDofCoords);
  auto hostFacesInternalDofOrdinals = Kokkos::create_mirror_view(facesInternalDofOrdinals);
  auto hostFaceTopoKey = Kokkos::create_mirror_view(faceTopoKey);
  for (ordinal_type iface=0; iface < numFaces; ++iface) {
    auto hostFaceBasisPtr = faceBases[iface]->getHostBasis();
    hostFaceTopoKey(iface) = hostFaceBasisPtr->getBaseCellTopology().getBaseKey();
    ordinal_type faceBasisCardinality = hostFaceBasisPtr->getCardinality();
    ScalarViewTypeHost  faceDofCoords("faceDofCoords", faceBasisCardinality, 2);
    hostFaceBasisPtr->getDofCoords(faceDofCoords);
    for(ordinal_type i=0; i<hostNumFacesInternalDofs(iface); ++i) {
      hostFacesInternalDofOrdinals(iface, i) = hostFaceBasisPtr->getDofOrdinal(2, 0, i);
      for(ordinal_type d=0; d <2; ++d)
        hostFacesInternalDofCoords(iface, i,d) = faceDofCoords(hostFacesInternalDofOrdinals(iface, i),d);
    }
  }
  Kokkos::deep_copy(facesInternalDofCoords,hostFacesInternalDofCoords);
  Kokkos::deep_copy(faceTopoKey,hostFaceTopoKey);

  typename RefSubcellParametrization<DeviceType>::ConstViewType faceParam;
  if(dim > 2)
    faceParam = RefSubcellParametrization<DeviceType>::get(2, topo.getKey());


  //*** Loop over cells ***

  const Kokkos::RangePolicy<typename DeviceType::execution_space> policy(0, dofCoords.extent(0));
  using FunctorType = FunctorsLagrangianTools::computeDofCoords<decltype(dofCoords),
      decltype(orts),
      decltype(tagToOrdinal),
      decltype(edgeParam),
      intViewType,
      ScalarViewType>;
  Kokkos::parallel_for(policy,
      FunctorType(dofCoords,
          orts, tagToOrdinal, edgeParam, faceParam,
          edgesInternalDofCoords, facesInternalDofCoords,
          dim, numEdges, numFaces,
          edgeTopoKey, numEdgesInternalDofs,
          faceTopoKey, numFacesInternalDofs));
}


template<typename DeviceType>
template<typename BasisType, 
class ...coeffsProperties,
typename ortValueType, class ...ortProperties>
void
LagrangianTools<DeviceType>::getOrientedDofCoeffs(
    Kokkos::DynRankView<typename BasisType::scalarType, coeffsProperties...> dofCoeffs,
    const BasisType* basis,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts) {

  using ScalarViewType = Kokkos::DynRankView<typename BasisType::scalarType, DeviceType>;
  ScalarViewType refDofCoeffs;
  if(dofCoeffs.rank() == 3) //vector basis
    refDofCoeffs = ScalarViewType("refDofCoeffs", dofCoeffs.extent(1), dofCoeffs.extent(2));
  else //scalar basis
    refDofCoeffs = ScalarViewType("refDofCoeffs",dofCoeffs.extent(1));
  basis->getDofCoeffs(refDofCoeffs); 

  OrientationTools<DeviceType>::modifyBasisByOrientationInverse(dofCoeffs, refDofCoeffs, orts, basis, true);
}


template<typename DeviceType>
template<typename basisCoeffsViewType,
typename funcViewType,
typename BasisType,
typename ortViewType>
void
LagrangianInterpolation<DeviceType>::getBasisCoeffs(basisCoeffsViewType basisCoeffs,
    const funcViewType functionValsAtDofCoords,
    const BasisType* cellBasis,
    const ortViewType  orts){
  using scalarType = typename BasisType::scalarType;
  using ScalarViewType = Kokkos::DynRankView<scalarType, DeviceType>;
  auto dofCoeffs = (functionValsAtDofCoords.rank() == 3) ? 
                   ScalarViewType("dofCoeffs", basisCoeffs.extent(0), cellBasis->getCardinality(), functionValsAtDofCoords.extent(2)) : 
                   ScalarViewType("dofCoeffs", basisCoeffs.extent(0), cellBasis->getCardinality()); 
  auto dofCoeffs0 = Kokkos::subview(dofCoeffs, 0, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
  cellBasis->getDofCoeffs(dofCoeffs0);
  RealSpaceTools<DeviceType>::clone(dofCoeffs,dofCoeffs0);
  Kokkos::DynRankView<typename basisCoeffsViewType::value_type, DeviceType> basisCoeffsRef("basisCoeffsRef", basisCoeffs.extent(0), basisCoeffs.extent(1));
  ArrayTools<DeviceType>::dotMultiplyDataData(basisCoeffsRef,functionValsAtDofCoords,dofCoeffs);
  OrientationTools<DeviceType>::modifyBasisByOrientationInverse(basisCoeffs, basisCoeffsRef, orts, cellBasis, true);
}

} // Intrepid2 namespace

#endif

