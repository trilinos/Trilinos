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



template<typename scalarViewType,
typename ortViewType,
typename t2oViewType,
typename subcellParamViewType,
typename intViewType>
struct F_computeDofCoordsAdCoeff {
  typedef typename scalarViewType::value_type value_type;

  scalarViewType dofCoords, dofCoeffs;
  const ortViewType orts;
  const t2oViewType tagToOrdinal;
  const subcellParamViewType edgeParam, faceParam;
  const intViewType edgeInternalDofOrdinals, facesInternalDofOrdinals;
  const scalarViewType edgeInternalDofCoords, edgeDofCoeffs, ortJacobianEdge, refEdgesTan, refEdgesNormal;
  const scalarViewType facesInternalDofCoords, faceDofCoeffs, ortJacobianFace, refFaceTangents, refFacesNormal;
  scalarViewType edgeWorkView, faceWorkView;
  const ordinal_type cellDim, numEdges, numFaces;
  const ordinal_type edgeTopoKey, numEdgeInternalDofs;
  const intViewType faceTopoKey, numFacesInternalDofs;
  const value_type edgeScale, faceScale;
  const bool isBasisHCURL, isBasisHDIV;

  F_computeDofCoordsAdCoeff( scalarViewType dofCoords_,
      scalarViewType dofCoeffs_,
      const ortViewType orts_,
      const t2oViewType tagToOrdinal_,
      const subcellParamViewType edgeParam_,
      const subcellParamViewType faceParam_,
      const intViewType edgeInternalDofOrdinals_,
      const intViewType facesInternalDofOrdinals_,
      const scalarViewType edgeInternalDofCoords_,
      const scalarViewType edgeDofCoeffs_,
      const scalarViewType refEdgesTan_,
      const scalarViewType refEdgesNormal_,
      const scalarViewType facesInternalDofCoords_,
      const scalarViewType faceDofCoeffs_,
      const scalarViewType refFaceTangents_,
      const scalarViewType refFacesNormal_,
      const ordinal_type cellDim_,
      const ordinal_type numEdges_,
      const ordinal_type numFaces_,
      const ordinal_type edgeTopoKey_,
      const ordinal_type numEdgeInternalDofs_,
      const intViewType faceTopoKey_,
      const intViewType numFacesInternalDofs_,
      const value_type edgeScale_,
      const value_type faceScale_,
      const bool isBasisHCURL_,
      const bool isBasisHDIV_
  )
  : dofCoords(dofCoords_),
    dofCoeffs(dofCoeffs_),
    orts(orts_),
    tagToOrdinal(tagToOrdinal_),
    edgeParam(edgeParam_),
    faceParam(faceParam_),
    edgeInternalDofOrdinals(edgeInternalDofOrdinals_),
    facesInternalDofOrdinals(facesInternalDofOrdinals_),
    edgeInternalDofCoords(edgeInternalDofCoords_),
    edgeDofCoeffs(edgeDofCoeffs_),
    refEdgesTan(refEdgesTan_),
    refEdgesNormal(refEdgesNormal_),
    facesInternalDofCoords(facesInternalDofCoords_),
    faceDofCoeffs(faceDofCoeffs_),
    refFaceTangents(refFaceTangents_),
    refFacesNormal(refFacesNormal_),
    cellDim(cellDim_),
    numEdges(numEdges_),
    numFaces(numFaces_),
    edgeTopoKey(edgeTopoKey_),
    numEdgeInternalDofs(numEdgeInternalDofs_),
    faceTopoKey(faceTopoKey_),
    numFacesInternalDofs(numFacesInternalDofs_),
    edgeScale(edgeScale_),
    faceScale(faceScale_),
    isBasisHCURL(isBasisHCURL_),
    isBasisHDIV(isBasisHDIV_)
  {
    if(numEdges > 0)
      edgeWorkView = scalarViewType("edgeWorkView", dofCoords.extent(0), numEdgeInternalDofs, 1);
    if(numFaces > 0)
      faceWorkView = scalarViewType("faceWorkView", dofCoords.extent(0), facesInternalDofCoords.extent(1), 2);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type cell) const {
    typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;


    if(numEdges > 0) {
      //compute coordinates associated to edge DoFs
      ordinal_type eOrt[12];
      value_type ortJac;
      scalarViewType ortJacobianEdge(&ortJac, 1, 1);
      orts(cell).getEdgeOrientation(eOrt, numEdges);
      auto edgeInternalDofCoordsOriented = Kokkos::subview(edgeWorkView,cell, Kokkos::ALL(), Kokkos::ALL());
      //map edge DoFs coords into parent element coords
      for (ordinal_type iedge=0; iedge < numEdges; ++iedge) {
        Impl::OrientationTools::mapToModifiedReference(edgeInternalDofCoordsOriented,edgeInternalDofCoords,edgeTopoKey, eOrt[iedge]);

        for (ordinal_type j=0;j<numEdgeInternalDofs;++j) {
          const auto u = edgeInternalDofCoordsOriented(j, 0);
          auto idof = tagToOrdinal(1, iedge, j);
          for (ordinal_type d=0;d<cellDim;++d)
            dofCoords(cell,idof,d) = edgeParam(iedge, d, 0) + edgeParam(iedge, d, 1)*u ;
        }
      }

      //compute coefficients associated to edge DoFs
      if(isBasisHCURL) {
        for (ordinal_type iedge=0; iedge < numEdges; ++iedge) {
          Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobianEdge, edgeTopoKey, eOrt[iedge]);
          for(ordinal_type j=0; j<numEdgeInternalDofs; ++j) {
            auto idof = tagToOrdinal(1, iedge, j);
            auto jdof = edgeInternalDofOrdinals(j);
            for(ordinal_type d=0; d <cellDim; ++d) {
              dofCoeffs(cell,idof,d) = 0;
              for(ordinal_type k=0; k <1; ++k)
                for(ordinal_type l=0; l <1; ++l)
                  dofCoeffs(cell,idof,d) += refEdgesTan(iedge,d)*ortJacobianEdge(0,0)*edgeDofCoeffs(jdof)*edgeScale;
            }
          }
        }
      } else if(isBasisHDIV) {
        for (ordinal_type iedge=0; iedge < numEdges; ++iedge) {
          Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobianEdge, edgeTopoKey, eOrt[iedge]);
          auto ortJacobianDet = ortJacobianEdge(0,0);
          for(ordinal_type j=0; j<numEdgeInternalDofs; ++j) {
            auto idof = tagToOrdinal(1, iedge, j);
            auto jdof = edgeInternalDofOrdinals(j);
            for(ordinal_type d=0; d <cellDim; ++d)
              dofCoeffs(cell,idof,d) = refEdgesNormal(iedge, d)*ortJacobianDet*edgeDofCoeffs(jdof)*edgeScale;
          }
        }
      } else {
        for (ordinal_type iedge=0; iedge < numEdges; ++iedge) {
          for(ordinal_type j=0; j<numEdgeInternalDofs; ++j) {
            auto idof = tagToOrdinal(1, iedge, j);
            auto jdof = edgeInternalDofOrdinals(j);
            dofCoeffs(cell,idof,0) = edgeDofCoeffs(jdof);
          }
        }

      }
    }

    if(numFaces > 0) {
      //compute coordinates associated to face DoFs
      ordinal_type fOrt[12];
      value_type ortJac[4];
      scalarViewType ortJacobianFace(ortJac, 2, 2);
      orts(cell).getFaceOrientation(fOrt, numFaces);
      //map face dofs coords into parent element coords
      for (ordinal_type iface=0; iface < numFaces; ++iface) {
        ordinal_type ort = fOrt[iface];
        ordinal_type numInternalDofs = numFacesInternalDofs(iface);
        auto dofRange = range_type(0, numInternalDofs);
        auto faceInternalDofCoords = Kokkos::subview(facesInternalDofCoords, iface, dofRange, Kokkos::ALL());
        auto faceInternalDofCoordsOriented = Kokkos::subview(faceWorkView,cell, dofRange, Kokkos::ALL());
        Impl::OrientationTools::mapToModifiedReference(faceInternalDofCoordsOriented,faceInternalDofCoords,faceTopoKey(iface),ort);

        for (ordinal_type j=0;j<numInternalDofs;++j) {
          const auto u = faceInternalDofCoordsOriented(j, 0);
          const auto v = faceInternalDofCoordsOriented(j, 1);
          auto idof = tagToOrdinal(2, iface, j);
          for (ordinal_type d=0;d<cellDim;++d)
            dofCoords(cell,idof,d) = faceParam(iface, d, 0) + faceParam(iface, d, 1)*u +  faceParam(iface, d, 2)*v;
        }
      }
      //compute coefficients associated to face DoFs
      if(isBasisHCURL) {
        for (ordinal_type iface=0; iface < numFaces; ++iface) {
          Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobianFace, faceTopoKey(iface), fOrt[iface]);
          for(ordinal_type j=0; j<numFacesInternalDofs(iface); ++j) {
            auto idof = tagToOrdinal(2, iface, j);
            auto jdof = facesInternalDofOrdinals(iface, j);
            for(ordinal_type d=0; d <cellDim; ++d) {
              dofCoeffs(cell,idof,d) = 0;
              for(ordinal_type k=0; k <2; ++k)
                for(ordinal_type l=0; l <2; ++l)
                  dofCoeffs(cell,idof,d) += refFaceTangents(iface, d,l)*ortJacobianFace(l,k)*faceDofCoeffs(iface,jdof,k);
            }
          }
        }
      } else if(isBasisHDIV) {
        for (ordinal_type iface=0; iface < numFaces; ++iface) {
          Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobianFace, faceTopoKey(iface), fOrt[iface]);
          auto ortJacobianDet = ortJacobianFace(0,0)*ortJacobianFace(1,1)-ortJacobianFace(1,0)*ortJacobianFace(0,1);
          for(ordinal_type j=0; j<numFacesInternalDofs(iface); ++j) {
            auto idof = tagToOrdinal(2, iface, j);
            auto jdof = facesInternalDofOrdinals(iface, j);
            for(ordinal_type d=0; d <cellDim; ++d)
              dofCoeffs(cell,idof,d) = refFacesNormal(iface,d)*ortJacobianDet*faceDofCoeffs(iface,jdof)*faceScale;
          }
        }
      } else {

        for (ordinal_type iface=0; iface < numFaces; ++iface) {
          for(ordinal_type j=0; j<numFacesInternalDofs(iface); ++j) {
            auto idof = tagToOrdinal(2, iface, j);
            auto jdof = facesInternalDofOrdinals(iface, j);
            dofCoeffs(cell,idof,0) = faceDofCoeffs(iface,jdof);
          }
        }
      }
    }
  }
};

template<typename SpT>
template<typename BasisType,
class ...coordsProperties, class ...coeffsProperties,
typename ortValueType, class ...ortProperties>
void
LagrangianInterpolation<SpT>::getDofCoordsAndCoeffs(
    Kokkos::DynRankView<typename BasisType::scalarType, coordsProperties...> dofCoords,
    Kokkos::DynRankView<typename BasisType::scalarType, coeffsProperties...> dofCoeffs,
    const BasisType* basis,
    EPointType pointType,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts) {

  typedef typename BasisType::scalarType scalarType;
  typedef  Kokkos::DynRankView<scalarType, SpT> scalarViewType;
  typedef  Kokkos::DynRankView<ordinal_type, SpT> intViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

  const auto topo = basis->getBaseCellTopology();
  const std::string name(basis->getName());
  const ordinal_type degree(basis->getDegree());

  bool isBasisHCURL = name.find("HCURL") != std::string::npos;
  bool isBasisHDIV = name.find("HDIV") != std::string::npos;
  bool isBasisTriOrTet = (name.find("TRI") != std::string::npos) || (name.find("TET") != std::string::npos);
  bool isBasisHEXI1 = (name.find("HEX_I1") != std::string::npos);
  bool isBasisTETI1 = (name.find("TET_I1") != std::string::npos);
  bool isBasisI1 = (name.find("I1") != std::string::npos);

  ordinal_type numEdges = (basis->getDofCount(1, 0) > 0) ? topo.getEdgeCount() : 0;
  ordinal_type numFaces = (basis->getDofCount(2, 0) > 0) ? topo.getFaceCount() : 0;

  Teuchos::RCP<Basis<SpT,scalarType,scalarType> > faceBases[6];
  Teuchos::RCP<Basis<SpT,scalarType,scalarType> > edgeBasis;


  if ((name == "Intrepid2_HGRAD_QUAD_Cn_FEM") || (name == "Intrepid2_HGRAD_QUAD_C1_FEM")) {
    edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree, pointType) );
  } else if ((name == "Intrepid2_HGRAD_HEX_Cn_FEM") || (name == "Intrepid2_HGRAD_HEX_C1_FEM")){
    edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree, pointType) );
    faceBases[0] = Teuchos::rcp(new Basis_HGRAD_QUAD_Cn_FEM<SpT,scalarType,scalarType>(degree, pointType) );
    for(ordinal_type i=1; i<6; ++i) faceBases[i]=faceBases[0];
  } else if ((name == "Intrepid2_HGRAD_TRI_Cn_FEM") || (name == "Intrepid2_HGRAD_TRI_C1_FEM")) {
    edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree, pointType) );
  } else if ((name == "Intrepid2_HGRAD_TET_Cn_FEM") || (name == "Intrepid2_HGRAD_TET_C1_FEM")) {
    edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree, pointType) );
    faceBases[0] = Teuchos::rcp(new Basis_HGRAD_TRI_Cn_FEM<SpT,scalarType,scalarType>(degree, pointType) );
    for(ordinal_type i=1; i<4; ++i) faceBases[i]=faceBases[0];
  } else if ((name == "Intrepid2_HCURL_QUAD_In_FEM") || (name == "Intrepid2_HCURL_QUAD_I1_FEM")) {
    edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree-1, POINTTYPE_GAUSS) );
  } else if ((name == "Intrepid2_HCURL_HEX_In_FEM") || (name == "Intrepid2_HCURL_HEX_I1_FEM")) {
    edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree-1, POINTTYPE_GAUSS) );
    faceBases[0] = Teuchos::rcp(new Basis_HCURL_QUAD_In_FEM<SpT,scalarType,scalarType>(degree, pointType) );
    for(ordinal_type i=1; i<6; ++i) faceBases[i]=faceBases[0];
  } else if ((name == "Intrepid2_HCURL_TRI_In_FEM") || (name == "Intrepid2_HCURL_TRI_I1_FEM")) {
    edgeBasis = Teuchos::rcp(new Basis_HVOL_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree-1, pointType) );
  } else if ((name == "Intrepid2_HCURL_TET_In_FEM") || (name == "Intrepid2_HCURL_TET_I1_FEM")) {
    edgeBasis = Teuchos::rcp(new Basis_HVOL_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree-1, pointType) );
    faceBases[0] = Teuchos::rcp(new Basis_HCURL_TRI_In_FEM<SpT,scalarType,scalarType>(degree, pointType) );
    for(ordinal_type i=1; i<4; ++i) faceBases[i]=faceBases[0];
  } else if ((name == "Intrepid2_HDIV_QUAD_In_FEM") || (name == "Intrepid2_HDIV_QUAD_I1_FEM")) {
    edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree-1, POINTTYPE_GAUSS) );
  } else if ((name == "Intrepid2_HDIV_HEX_In_FEM") || (name == "Intrepid2_HDIV_HEX_I1_FEM")) {
    edgeBasis = Teuchos::null;
    faceBases[0] = Teuchos::rcp(new Basis_HGRAD_QUAD_Cn_FEM<SpT,scalarType,scalarType>(degree-1, POINTTYPE_GAUSS) );
    for(ordinal_type i=1; i<6; ++i) faceBases[i]=faceBases[0];
  } else if ((name == "Intrepid2_HDIV_TRI_In_FEM") || (name == "Intrepid2_HDIV_TRI_I1_FEM")) {
    edgeBasis = Teuchos::rcp(new Basis_HVOL_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree-1, pointType) );
  } else if ((name == "Intrepid2_HDIV_TET_In_FEM") || (name == "Intrepid2_HDIV_TET_I1_FEM")) {
    edgeBasis = Teuchos::null;
    faceBases[0] = Teuchos::rcp(new Basis_HVOL_TRI_Cn_FEM<SpT,scalarType,scalarType>(degree-1, pointType) );
    for(ordinal_type i=1; i<4; ++i) faceBases[i]=faceBases[0];
  } else { //HVOL element does not have any face or edge DOFs
    //Throw error when basis is not HVOL.
    INTREPID2_TEST_FOR_ABORT(name.find("HVOL") == std::string::npos,
        ">>> ERROR (Intrepid2::Experimental::LagrangianInterpolation<SpT>::getDofCoordsAndCoeffs): " \
        "method not implemented for this basis function");
  }

  auto tagToOrdinal = Kokkos::create_mirror_view(typename SpT::memory_space(), basis->getAllDofOrdinal());
  Kokkos::deep_copy(tagToOrdinal, basis->getAllDofOrdinal());

  const ordinal_type dim = topo.getDimension();

  const ordinal_type numCells = dofCoeffs.extent(0);

  scalarViewType refDofCoords("refDofCoords", dofCoords.extent(1), dofCoords.extent(2)), refDofCoeffs;
  basis->getDofCoords(refDofCoords);
  RealSpaceTools<SpT>::clone(dofCoords,refDofCoords);

  if(dofCoeffs.rank() == 3) //vector basis
    refDofCoeffs = scalarViewType("refDofCoeffs", dofCoeffs.extent(1), dofCoeffs.extent(2));
  else //scalar basis
    refDofCoeffs = scalarViewType("refDofCoeffs",dofCoeffs.extent(1));
  basis->getDofCoeffs(refDofCoeffs);
  RealSpaceTools<SpT>::clone(dofCoeffs,refDofCoeffs);

  //*** Pre-compute needed quantities related to edge DoFs that do not depend on the cell ***

  ordinal_type edgeTopoKey = Teuchos::nonnull(edgeBasis) ? edgeBasis->getBaseCellTopology().getBaseKey() : 0;
  intViewType eOrt("eOrt", numEdges);
  scalarViewType refEdgesTan("refEdgesTan",  numEdges, dim);
  scalarViewType refEdgesNormal("refEdgesNormal",  numEdges, dim);
  scalarViewType edgeParam;
  ordinal_type edgeBasisCardinality = Teuchos::nonnull(edgeBasis) ? edgeBasis->getCardinality() : ordinal_type(0);
  ordinal_type numEdgeInternalDofs = Teuchos::nonnull(edgeBasis) ? edgeBasis->getDofCount(1,0) : ordinal_type(0);
  scalarViewType edgeDofCoords("edgeDofCoords", edgeBasisCardinality, 1);
  scalarViewType edgeDofCoeffs("edgeDofCoeffs", edgeBasisCardinality);
  scalarViewType edgeInternalDofCoords("edgeInternalDofCoords", numEdgeInternalDofs, 1);
  intViewType edgeInternalDofOrdinals("edgeInternalDofOrdinals", numEdgeInternalDofs);
  //depending on how the reference basis is defined, the edges are scaled differently
  auto edgeScale = (isBasisTriOrTet||isBasisI1) ? 2.0 :1.0;

  if(Teuchos::nonnull(edgeBasis)) {
    edgeBasis->getDofCoords(edgeDofCoords);
    edgeBasis->getDofCoeffs(edgeDofCoeffs);
  }

  for (ordinal_type iedge=0; iedge < numEdges; ++iedge) {
    if(isBasisHCURL) {
      auto edgeTan = Kokkos::subview(refEdgesTan, iedge, Kokkos::ALL);
      auto edgeTanHost = Kokkos::create_mirror_view(edgeTan);
      CellTools<SpT>::getReferenceEdgeTangent(edgeTanHost, iedge, topo);
      Kokkos::deep_copy(edgeTan,edgeTanHost);
    }
    else if(isBasisHDIV) {
      auto edgeNormal = Kokkos::subview(refEdgesNormal, iedge, Kokkos::ALL);
      auto edgeNormalHost = Kokkos::create_mirror_view(edgeNormal);
      CellTools<SpT>::getReferenceSideNormal(edgeNormalHost, iedge, topo);
      Kokkos::deep_copy(edgeNormal,edgeNormalHost);
    }
  }

  //compute DofCoords Oriented
  for(ordinal_type i=0; i<numEdgeInternalDofs; ++i) {
    edgeInternalDofOrdinals(i) = edgeBasis->getDofOrdinal(1, 0, i);
    edgeInternalDofCoords(i,0) = edgeDofCoords(edgeInternalDofOrdinals(i), 0);
    CellTools<SpT>::getSubcellParametrization(edgeParam, 1, topo);
  }


  //*** Pre-compute needed quantities related to face DoFs that do not depend on the cell ***

  intViewType faceTopoKey("faceTopoKey",numFaces);
  intViewType fOrt("fOrt", numFaces);
  scalarViewType refFaceTangents("refFaceTangents", numFaces, dim, 2);
  scalarViewType refFacesNormal("refFacesNormal",  numFaces, dim);
  scalarViewType faceParam;
  intViewType numFacesInternalDofs("numFacesInternalDofs", numFaces);
  scalarViewType  facesInternalDofCoords;
  intViewType  facesInternalDofOrdinals;
  scalarViewType  faceDofCoeffs;
  //depending on how the reference basis is defined, the faces are scaled differently
  auto faceScale = (isBasisHEXI1) ? 4.0 :
      (isBasisTETI1) ? 0.5 : 1.0;


  ordinal_type maxNumFacesInternalDofs=0;
  ordinal_type faceBasisMaxCardinality=0;

  for (ordinal_type iface=0; iface < numFaces; ++iface) {
    ordinal_type numInternalDofs = faceBases[iface]->getDofCount(2,0);
    numFacesInternalDofs(iface) = numInternalDofs;
    maxNumFacesInternalDofs = std::max(maxNumFacesInternalDofs,numInternalDofs);
    ordinal_type faceBasisCardinality = faceBases[iface]->getCardinality();
    faceBasisMaxCardinality = std::max(faceBasisMaxCardinality, faceBasisCardinality);
  }

  facesInternalDofCoords = scalarViewType("faceInternalDofCoords", numFaces, maxNumFacesInternalDofs, 2);
  facesInternalDofOrdinals = intViewType("faceInternalDofCoords", numFaces, maxNumFacesInternalDofs);

  if(isBasisHCURL)
    faceDofCoeffs = scalarViewType("faceDofCoeffs", numFaces, faceBasisMaxCardinality,2);
  else
    faceDofCoeffs = scalarViewType("faceDofCoeffs", numFaces, faceBasisMaxCardinality);

  for (ordinal_type iface=0; iface < numFaces; ++iface) {
    auto faceBasis = faceBases[iface];
    faceTopoKey(iface) = faceBasis->getBaseCellTopology().getBaseKey();
    ordinal_type faceBasisCardinality = faceBasis->getCardinality();
    scalarViewType  faceDofCoords("faceDofCoords", faceBasisCardinality, 2);
    faceBasis->getDofCoords(faceDofCoords);
    for(ordinal_type i=0; i<numFacesInternalDofs(iface); ++i) {
      facesInternalDofOrdinals(iface, i) = faceBasis->getDofOrdinal(2, 0, i);
      for(ordinal_type d=0; d <2; ++d)
        facesInternalDofCoords(iface, i,d) = faceDofCoords(facesInternalDofOrdinals(iface, i),d);
    }

    auto dofRange = range_type(0, faceBasis->getCardinality());
    faceBasis->getDofCoeffs(Kokkos::subview(faceDofCoeffs, iface, dofRange, Kokkos::ALL()));

    if(isBasisHCURL) {
      auto refFaceTanU = Kokkos::subview(refFaceTangents, iface, Kokkos::ALL, 0);
      auto refFaceTanV = Kokkos::subview(refFaceTangents, iface, Kokkos::ALL,1);
      auto refFaceTanUHost = Kokkos::create_mirror_view(refFaceTanU);
      auto refFaceTanVHost = Kokkos::create_mirror_view(refFaceTanV);
      CellTools<SpT>::getReferenceFaceTangents(refFaceTanUHost, refFaceTanVHost, iface, topo);
      Kokkos::deep_copy(refFaceTanU, refFaceTanUHost);
      Kokkos::deep_copy(refFaceTanV, refFaceTanVHost);
    } else if(isBasisHDIV) {
      auto faceNormal = Kokkos::subview(refFacesNormal,iface,Kokkos::ALL());
      auto faceNormalHost = Kokkos::create_mirror_view(faceNormal);
      CellTools<SpT>::getReferenceFaceNormal(faceNormalHost, iface, topo);
      Kokkos::deep_copy(faceNormal, faceNormalHost);
    }
  }

  if(dim > 2)
    CellTools<SpT>::getSubcellParametrization(faceParam, 2, topo);


  //*** Loop over cells ***

  const Kokkos::RangePolicy<SpT> policy(0, numCells);
  typedef F_computeDofCoordsAdCoeff
      <scalarViewType,
      decltype(orts),
      decltype(tagToOrdinal),
      decltype(edgeParam),
      intViewType> FunctorType;
  Kokkos::parallel_for(policy,
      FunctorType(dofCoords, dofCoeffs,
          orts, tagToOrdinal, edgeParam, faceParam,
          edgeInternalDofOrdinals, facesInternalDofOrdinals,
          edgeInternalDofCoords, edgeDofCoeffs, refEdgesTan, refEdgesNormal,
          facesInternalDofCoords, faceDofCoeffs, refFaceTangents, refFacesNormal,
          dim, numEdges, numFaces,
          edgeTopoKey, numEdgeInternalDofs,
          faceTopoKey, numFacesInternalDofs,
          edgeScale, faceScale,
          isBasisHCURL, isBasisHDIV));
}


template<typename SpT>
template<typename basisCoeffsViewType,
typename funcViewType,
typename dofCoeffViewType>
void
LagrangianInterpolation<SpT>::getBasisCoeffs(basisCoeffsViewType basisCoeffs,
    const funcViewType functionValsAtDofCoords,
    const dofCoeffViewType dofCoeffs){
  ArrayTools<SpT>::dotMultiplyDataData(basisCoeffs,functionValsAtDofCoords,dofCoeffs);
}
}
}

#endif

