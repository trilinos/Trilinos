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

/** \file   Intrepid2_LagrangianInterpolationDef.hpp
    \brief  Header file for the Intrepid2::Experimental::LagrangianInterpolation containing definitions.
    \author Created by Mauro Perego
 */

#ifndef __INTREPID2_LAGRANGIANINTERPOLATIONDEF_HPP__
#define __INTREPID2_LAGRANGIANINTERPOLATIONDEF_HPP__

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


namespace Intrepid2 {
namespace Experimental {



template<typename CoordsViewType,
typename CoeffsViewType,
typename ortViewType,
typename t2oViewType,
typename subcellParamViewType,
typename intViewType,
typename ScalarViewType>
struct computeDofCoordsAndCoeffs {
  typedef typename ScalarViewType::value_type value_type;

  CoordsViewType dofCoords_;
  CoeffsViewType dofCoeffs_;
  const ortViewType orts_;
  const t2oViewType tagToOrdinal_;
  const subcellParamViewType edgeParam_, faceParam_;
  const intViewType edgesInternalDofOrdinals_, facesInternalDofOrdinals_;
  const ScalarViewType edgesInternalDofCoords_, edgeDofCoeffs_;
  const ScalarViewType facesInternalDofCoords_, faceDofCoeffs_;
  ScalarViewType edgeWorkView_, faceWorkView_;
  const ordinal_type cellDim_, numEdges_, numFaces_;
  const intViewType edgeTopoKey_, numEdgesInternalDofs_;
  const intViewType faceTopoKey_, numFacesInternalDofs_;
  const bool isBasisHCURL_, isBasisHDIV_;

  computeDofCoordsAndCoeffs( CoordsViewType dofCoords,
      CoeffsViewType dofCoeffs,
      const ortViewType orts,
      const t2oViewType tagToOrdinal,
      const subcellParamViewType edgeParam,
      const subcellParamViewType faceParam,
      const intViewType edgesInternalDofOrdinals,
      const intViewType facesInternalDofOrdinals,
      const ScalarViewType edgesInternalDofCoords,
      const ScalarViewType edgeDofCoeffs,
      const ScalarViewType facesInternalDofCoords,
      const ScalarViewType faceDofCoeffs,
      const ordinal_type cellDim,
      const ordinal_type numEdges,
      const ordinal_type numFaces,
      const intViewType edgeTopoKey,
      const intViewType numEdgesInternalDofs,
      const intViewType faceTopoKey,
      const intViewType numFacesInternalDofs,
      const bool isBasisHCURL,
      const bool isBasisHDIV
  )
  : dofCoords_(dofCoords),
    dofCoeffs_(dofCoeffs),
    orts_(orts),
    tagToOrdinal_(tagToOrdinal),
    edgeParam_(edgeParam),
    faceParam_(faceParam),
    edgesInternalDofOrdinals_(edgesInternalDofOrdinals),
    facesInternalDofOrdinals_(facesInternalDofOrdinals),
    edgesInternalDofCoords_(edgesInternalDofCoords),
    edgeDofCoeffs_(edgeDofCoeffs),
    facesInternalDofCoords_(facesInternalDofCoords),
    faceDofCoeffs_(faceDofCoeffs),
    cellDim_(cellDim),
    numEdges_(numEdges),
    numFaces_(numFaces),
    edgeTopoKey_(edgeTopoKey),
    numEdgesInternalDofs_(numEdgesInternalDofs),
    faceTopoKey_(faceTopoKey),
    numFacesInternalDofs_(numFacesInternalDofs),
    isBasisHCURL_(isBasisHCURL),
    isBasisHDIV_(isBasisHDIV)
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

      //compute coefficients associated to edge DoFs
      if(isBasisHCURL_) {
        value_type tmp[3];
        ScalarViewType tangents(tmp,1,cellDim_);
        for (ordinal_type iedge=0; iedge < numEdges_; ++iedge) {
          Impl::OrientationTools::getRefSubcellTangents(tangents, edgeParam_,edgeTopoKey_(iedge), iedge, eOrt[iedge]);
          for(ordinal_type j=0; j<numEdgesInternalDofs_(iedge); ++j) {
            auto idof = tagToOrdinal_(1, iedge, j);
            auto jdof = edgesInternalDofOrdinals_(iedge, j);
            for(ordinal_type d=0; d <cellDim_; ++d) {
              dofCoeffs_(cell,idof,d) = 0;
                  dofCoeffs_(cell,idof,d) += tangents(0,d)*edgeDofCoeffs_(iedge,jdof);
            }
          }
        }
      } else if(isBasisHDIV_) {
        value_type tmp[9];
        ScalarViewType tangentsAndNormal(tmp,cellDim_,cellDim_);
        for (ordinal_type iedge=0; iedge < numEdges_; ++iedge) {
          Impl::OrientationTools::getRefSideTangentsAndNormal(tangentsAndNormal, edgeParam_,edgeTopoKey_(iedge), iedge, eOrt[iedge]);
          for(ordinal_type j=0; j<numEdgesInternalDofs_(iedge); ++j) {
            auto idof = tagToOrdinal_(1, iedge, j);
            auto jdof = edgesInternalDofOrdinals_(iedge, j);
            for(ordinal_type d=0; d <cellDim_; ++d)
              dofCoeffs_(cell,idof,d) = tangentsAndNormal(cellDim_-1, d)*edgeDofCoeffs_(iedge,jdof);
          }
        }
      } else {
        for (ordinal_type iedge=0; iedge < numEdges_; ++iedge) {
          for(ordinal_type j=0; j<numEdgesInternalDofs_(iedge); ++j) {
            auto idof = tagToOrdinal_(1, iedge, j);
            auto jdof = edgesInternalDofOrdinals_(iedge, j);
            dofCoeffs_(cell,idof,0) = edgeDofCoeffs_(iedge,jdof);
          }
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
      //compute coefficients associated to face DoFs
      if(isBasisHCURL_) {
        value_type tmp[6];
        ScalarViewType tangents(tmp,2,cellDim_);
        for (ordinal_type iface=0; iface < numFaces_; ++iface) {
          Impl::OrientationTools::getRefSubcellTangents(tangents, faceParam_,faceTopoKey_(iface), iface, fOrt[iface]);
          for(ordinal_type j=0; j<numFacesInternalDofs_(iface); ++j) {
            auto idof = tagToOrdinal_(2, iface, j);
            auto jdof = facesInternalDofOrdinals_(iface, j);
            for(ordinal_type d=0; d <cellDim_; ++d) {
              dofCoeffs_(cell,idof,d) = 0;
              for(ordinal_type k=0; k <2; ++k)
                dofCoeffs_(cell,idof,d) += tangents(k,d)*faceDofCoeffs_(iface,jdof,k);
            }
          }
        }
      } else if(isBasisHDIV_) {
        value_type tmp[9];
        ScalarViewType tangentsAndNormal(tmp,cellDim_,cellDim_);
        for (ordinal_type iface=0; iface < numFaces_; ++iface) {
          Impl::OrientationTools::getRefSideTangentsAndNormal(tangentsAndNormal, faceParam_,faceTopoKey_(iface), iface, fOrt[iface]);
          for(ordinal_type j=0; j<numFacesInternalDofs_(iface); ++j) {
            auto idof = tagToOrdinal_(2, iface, j);
            auto jdof = facesInternalDofOrdinals_(iface, j);
            for(ordinal_type d=0; d <cellDim_; ++d)
              dofCoeffs_(cell,idof,d) = tangentsAndNormal(cellDim_-1,d)*faceDofCoeffs_(iface,jdof);
          }
        }
      } else {

        for (ordinal_type iface=0; iface < numFaces_; ++iface) {
          for(ordinal_type j=0; j<numFacesInternalDofs_(iface); ++j) {
            auto idof = tagToOrdinal_(2, iface, j);
            auto jdof = facesInternalDofOrdinals_(iface, j);
            dofCoeffs_(cell,idof,0) = faceDofCoeffs_(iface,jdof);
          }
        }
      }
    }
  }
};

template<typename DeviceType>
template<typename BasisType,
class ...coordsProperties, class ...coeffsProperties,
typename ortValueType, class ...ortProperties>
void
LagrangianInterpolation<DeviceType>::getDofCoordsAndCoeffs(
    Kokkos::DynRankView<typename BasisType::scalarType, coordsProperties...> dofCoords,
    Kokkos::DynRankView<typename BasisType::scalarType, coeffsProperties...> dofCoeffs,
    const BasisType* basis,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts) {

  using HostSpaceType = typename Kokkos::Impl::is_space<DeviceType>::host_mirror_space::execution_space;
  using scalarType = typename BasisType::scalarType;
  using ScalarViewType = Kokkos::DynRankView<scalarType, DeviceType>;
  using ScalarViewTypeHost = Kokkos::DynRankView<scalarType, HostSpaceType>;
  using intViewType = Kokkos::DynRankView<ordinal_type, DeviceType>;
  using range_type = Kokkos::pair<ordinal_type,ordinal_type>;

  const auto topo = basis->getBaseCellTopology();
  const std::string name(basis->getName());

  bool isBasisHCURL = (basis->getFunctionSpace()==FUNCTION_SPACE_HCURL);
  bool isBasisHDIV = (basis->getFunctionSpace()==FUNCTION_SPACE_HDIV);

  ordinal_type numEdges = (basis->getDofCount(1, 0) > 0) ? topo.getEdgeCount() : 0;
  ordinal_type numFaces = (basis->getDofCount(2, 0) > 0) ? topo.getFaceCount() : 0;

  std::vector<Teuchos::RCP<Basis<DeviceType,scalarType,scalarType> > > edgeBases, faceBases;

  for(int i=0;i<numEdges;++i)
    edgeBases.push_back(basis->getSubCellRefBasis(1,i));
  for(int i=0;i<numFaces;++i)
    faceBases.push_back(basis->getSubCellRefBasis(2,i));

  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(typename DeviceType::memory_space(), basis->getAllDofOrdinal());

  const ordinal_type dim = topo.getDimension();

  const ordinal_type numCells = dofCoeffs.extent(0);

  ScalarViewType refDofCoords("refDofCoords", dofCoords.extent(1), dofCoords.extent(2)), refDofCoeffs;
  basis->getDofCoords(refDofCoords);
  RealSpaceTools<DeviceType>::clone(dofCoords,refDofCoords);

  if(dofCoeffs.rank() == 3) //vector basis
    refDofCoeffs = ScalarViewType("refDofCoeffs", dofCoeffs.extent(1), dofCoeffs.extent(2));
  else //scalar basis
    refDofCoeffs = ScalarViewType("refDofCoeffs",dofCoeffs.extent(1));
  basis->getDofCoeffs(refDofCoeffs);
  RealSpaceTools<DeviceType>::clone(dofCoeffs,refDofCoeffs);

  if((numFaces == 0) && (numEdges == 0)) 
    return;

  //*** Pre-compute needed quantities related to edge DoFs that do not depend on the cell ***
  intViewType edgeTopoKey("edgeTopoKey",numEdges);
  intViewType sOrt("eOrt", numEdges);
  intViewType numEdgesInternalDofs("numEdgesInternalDofs", numEdges);
  ScalarViewType  edgesInternalDofCoords;
  intViewType  edgesInternalDofOrdinals;
  ScalarViewType  edgeDofCoeffs;

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

  edgeDofCoeffs = ScalarViewType("edgeDofCoeffs", numEdges, edgeBasisMaxCardinality);
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
    auto edgeBasisPtr = edgeBases[iedge];
    auto dofRange = range_type(0, edgeBasisPtr->getCardinality());
    edgeBasisPtr->getDofCoeffs(Kokkos::subview(edgeDofCoeffs, iedge, dofRange));
  }
  Kokkos::deep_copy(edgesInternalDofCoords,hostEdgesInternalDofCoords);
  Kokkos::deep_copy(edgesInternalDofOrdinals,hostEdgesInternalDofOrdinals);
  Kokkos::deep_copy(edgeTopoKey,hostEdgeTopoKey);

  auto edgeParam = RefSubcellParametrization<DeviceType>::get(1, topo.getKey());

  //*** Pre-compute needed quantities related to face DoFs that do not depend on the cell ***
  intViewType faceTopoKey("faceTopoKey",numFaces);
  intViewType fOrt("fOrt", numFaces);
  intViewType numFacesInternalDofs("numFacesInternalDofs", numFaces);
  ScalarViewType  facesInternalDofCoords;
  intViewType  facesInternalDofOrdinals;
  ScalarViewType  faceDofCoeffs;

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

  if(isBasisHCURL)
    faceDofCoeffs = ScalarViewType("faceDofCoeffs", numFaces, faceBasisMaxCardinality,2);
  else
    faceDofCoeffs = ScalarViewType("faceDofCoeffs", numFaces, faceBasisMaxCardinality);

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
    auto faceBasisPtr = faceBases[iface];
    auto dofRange = range_type(0, faceBasisPtr->getCardinality());
    faceBasisPtr->getDofCoeffs(Kokkos::subview(faceDofCoeffs, iface, dofRange, Kokkos::ALL()));
  }
  Kokkos::deep_copy(facesInternalDofCoords,hostFacesInternalDofCoords);
  Kokkos::deep_copy(facesInternalDofOrdinals,hostFacesInternalDofOrdinals);
  Kokkos::deep_copy(faceTopoKey,hostFaceTopoKey);

  typename RefSubcellParametrization<DeviceType>::ConstViewType faceParam;
  if(dim > 2)
    faceParam = RefSubcellParametrization<DeviceType>::get(2, topo.getKey());


  //*** Loop over cells ***

  const Kokkos::RangePolicy<typename DeviceType::execution_space> policy(0, numCells);
  typedef computeDofCoordsAndCoeffs
      <decltype(dofCoords),
      decltype(dofCoeffs),
      decltype(orts),
      decltype(tagToOrdinal),
      decltype(edgeParam),
      intViewType,
      ScalarViewType> FunctorType;
  Kokkos::parallel_for(policy,
      FunctorType(dofCoords, dofCoeffs,
          orts, tagToOrdinal, edgeParam, faceParam,
          edgesInternalDofOrdinals, facesInternalDofOrdinals,
          edgesInternalDofCoords, edgeDofCoeffs,
          facesInternalDofCoords, faceDofCoeffs,
          dim, numEdges, numFaces,
          edgeTopoKey, numEdgesInternalDofs,
          faceTopoKey, numFacesInternalDofs,
          isBasisHCURL, isBasisHDIV));
}


template<typename DeviceType>
template<typename basisCoeffsViewType,
typename funcViewType,
typename dofCoeffViewType>
void
LagrangianInterpolation<DeviceType>::getBasisCoeffs(basisCoeffsViewType basisCoeffs,
    const funcViewType functionValsAtDofCoords,
    const dofCoeffViewType dofCoeffs){
  ArrayTools<DeviceType>::dotMultiplyDataData(basisCoeffs,functionValsAtDofCoords,dofCoeffs);
}
}
}

#endif

