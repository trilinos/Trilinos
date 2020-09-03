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



template<typename ScalarViewType,
typename ortViewType,
typename t2oViewType,
typename subcellParamViewType,
typename intViewType>
struct computeDofCoordsAndCoeffs {
  typedef typename ScalarViewType::value_type value_type;

  ScalarViewType dofCoords_, dofCoeffs_;
  const ortViewType orts_;
  const t2oViewType tagToOrdinal_;
  const subcellParamViewType edgeParam_, faceParam_;
  const intViewType edgeInternalDofOrdinals_, facesInternalDofOrdinals_;
  const ScalarViewType edgeInternalDofCoords_, edgeDofCoeffs_, ortJacobianEdge_, refEdgesTan_, refEdgesNormal_;
  const ScalarViewType facesInternalDofCoords_, faceDofCoeffs_, ortJacobianFace_, refFaceTangents_, refFacesNormal_;
  ScalarViewType edgeWorkView_, faceWorkView_;
  const ordinal_type cellDim_, numEdges_, numFaces_;
  const ordinal_type edgeTopoKey_, numEdgeInternalDofs_;
  const intViewType faceTopoKey_, numFacesInternalDofs_;
  const value_type edgeScale_, faceScale_;
  const bool isBasisHCURL_, isBasisHDIV_;

  computeDofCoordsAndCoeffs( ScalarViewType dofCoords,
      ScalarViewType dofCoeffs,
      const ortViewType orts,
      const t2oViewType tagToOrdinal,
      const subcellParamViewType edgeParam,
      const subcellParamViewType faceParam,
      const intViewType edgeInternalDofOrdinals,
      const intViewType facesInternalDofOrdinals,
      const ScalarViewType edgeInternalDofCoords,
      const ScalarViewType edgeDofCoeffs,
      const ScalarViewType refEdgesTan,
      const ScalarViewType refEdgesNormal,
      const ScalarViewType facesInternalDofCoords,
      const ScalarViewType faceDofCoeffs,
      const ScalarViewType refFaceTangents,
      const ScalarViewType refFacesNormal,
      const ordinal_type cellDim,
      const ordinal_type numEdges,
      const ordinal_type numFaces,
      const ordinal_type edgeTopoKey,
      const ordinal_type numEdgeInternalDofs,
      const intViewType faceTopoKey,
      const intViewType numFacesInternalDofs,
      const value_type edgeScale,
      const value_type faceScale,
      const bool isBasisHCURL,
      const bool isBasisHDIV
  )
  : dofCoords_(dofCoords),
    dofCoeffs_(dofCoeffs),
    orts_(orts),
    tagToOrdinal_(tagToOrdinal),
    edgeParam_(edgeParam),
    faceParam_(faceParam),
    edgeInternalDofOrdinals_(edgeInternalDofOrdinals),
    facesInternalDofOrdinals_(facesInternalDofOrdinals),
    edgeInternalDofCoords_(edgeInternalDofCoords),
    edgeDofCoeffs_(edgeDofCoeffs),
    refEdgesTan_(refEdgesTan),
    refEdgesNormal_(refEdgesNormal),
    facesInternalDofCoords_(facesInternalDofCoords),
    faceDofCoeffs_(faceDofCoeffs),
    refFaceTangents_(refFaceTangents),
    refFacesNormal_(refFacesNormal),
    cellDim_(cellDim),
    numEdges_(numEdges),
    numFaces_(numFaces),
    edgeTopoKey_(edgeTopoKey),
    numEdgeInternalDofs_(numEdgeInternalDofs),
    faceTopoKey_(faceTopoKey),
    numFacesInternalDofs_(numFacesInternalDofs),
    edgeScale_(edgeScale),
    faceScale_(faceScale),
    isBasisHCURL_(isBasisHCURL),
    isBasisHDIV_(isBasisHDIV)
  {
    if(numEdges > 0)
      edgeWorkView_ = ScalarViewType("edgeWorkView", dofCoords.extent(0), numEdgeInternalDofs, 1);
    if(numFaces > 0)
      faceWorkView_ = ScalarViewType("faceWorkView", dofCoords.extent(0), facesInternalDofCoords.extent(1), 2);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type cell) const {
    typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;


    if(numEdges_ > 0) {
      //compute coordinates associated to edge DoFs
      ordinal_type eOrt[12];
      value_type ortJac;
      ScalarViewType ortJacobianEdge(&ortJac, 1, 1);
      orts_(cell).getEdgeOrientation(eOrt, numEdges_);
      auto edgeInternalDofCoordsOriented = Kokkos::subview(edgeWorkView_,cell, Kokkos::ALL(), Kokkos::ALL());
      //map edge DoFs coords into parent element coords
      for (ordinal_type iedge=0; iedge < numEdges_; ++iedge) {
        Impl::OrientationTools::mapToModifiedReference(edgeInternalDofCoordsOriented,edgeInternalDofCoords_,edgeTopoKey_, eOrt[iedge]);

        for (ordinal_type j=0;j<numEdgeInternalDofs_;++j) {
          const auto u = edgeInternalDofCoordsOriented(j, 0);
          auto idof = tagToOrdinal_(1, iedge, j);
          for (ordinal_type d=0;d<cellDim_;++d)
            dofCoords_(cell,idof,d) = edgeParam_(iedge, d, 0) + edgeParam_(iedge, d, 1)*u ;
        }
      }

      //compute coefficients associated to edge DoFs
      if(isBasisHCURL_) {
        for (ordinal_type iedge=0; iedge < numEdges_; ++iedge) {
          Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobianEdge, edgeTopoKey_, eOrt[iedge]);
          for(ordinal_type j=0; j<numEdgeInternalDofs_; ++j) {
            auto idof = tagToOrdinal_(1, iedge, j);
            auto jdof = edgeInternalDofOrdinals_(j);
            for(ordinal_type d=0; d <cellDim_; ++d) {
              dofCoeffs_(cell,idof,d) = 0;
              for(ordinal_type k=0; k <1; ++k)
                for(ordinal_type l=0; l <1; ++l)
                  dofCoeffs_(cell,idof,d) += refEdgesTan_(iedge,d)*ortJacobianEdge(0,0)*edgeDofCoeffs_(jdof)*edgeScale_;
            }
          }
        }
      } else if(isBasisHDIV_) {
        for (ordinal_type iedge=0; iedge < numEdges_; ++iedge) {
          Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobianEdge, edgeTopoKey_, eOrt[iedge]);
          auto ortJacobianDet = ortJacobianEdge(0,0);
          for(ordinal_type j=0; j<numEdgeInternalDofs_; ++j) {
            auto idof = tagToOrdinal_(1, iedge, j);
            auto jdof = edgeInternalDofOrdinals_(j);
            for(ordinal_type d=0; d <cellDim_; ++d)
              dofCoeffs_(cell,idof,d) = refEdgesNormal_(iedge, d)*ortJacobianDet*edgeDofCoeffs_(jdof)*edgeScale_;
          }
        }
      } else {
        for (ordinal_type iedge=0; iedge < numEdges_; ++iedge) {
          for(ordinal_type j=0; j<numEdgeInternalDofs_; ++j) {
            auto idof = tagToOrdinal_(1, iedge, j);
            auto jdof = edgeInternalDofOrdinals_(j);
            dofCoeffs_(cell,idof,0) = edgeDofCoeffs_(jdof);
          }
        }

      }
    }

    if(numFaces_ > 0) {
      //compute coordinates associated to face DoFs
      ordinal_type fOrt[12];
      value_type ortJac[4];
      ScalarViewType ortJacobianFace(ortJac, 2, 2);
      orts_(cell).getFaceOrientation(fOrt, numFaces_);
      //map face dofs coords into parent element coords
      for (ordinal_type iface=0; iface < numFaces_; ++iface) {
        ordinal_type ort = fOrt[iface];
        ordinal_type numInternalDofs = numFacesInternalDofs_(iface);
        auto dofRange = range_type(0, numInternalDofs);
        auto faceInternalDofCoords = Kokkos::subview(facesInternalDofCoords_, iface, dofRange, Kokkos::ALL());
        auto faceInternalDofCoordsOriented = Kokkos::subview(faceWorkView_,cell, dofRange, Kokkos::ALL());
        Impl::OrientationTools::mapToModifiedReference(faceInternalDofCoordsOriented,faceInternalDofCoords,faceTopoKey_(iface),ort);

        for (ordinal_type j=0;j<numInternalDofs;++j) {
          const auto u = faceInternalDofCoordsOriented(j, 0);
          const auto v = faceInternalDofCoordsOriented(j, 1);
          auto idof = tagToOrdinal_(2, iface, j);
          for (ordinal_type d=0;d<cellDim_;++d)
            dofCoords_(cell,idof,d) = faceParam_(iface, d, 0) + faceParam_(iface, d, 1)*u +  faceParam_(iface, d, 2)*v;
        }
      }
      //compute coefficients associated to face DoFs
      if(isBasisHCURL_) {
        for (ordinal_type iface=0; iface < numFaces_; ++iface) {
          Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobianFace, faceTopoKey_(iface), fOrt[iface]);
          for(ordinal_type j=0; j<numFacesInternalDofs_(iface); ++j) {
            auto idof = tagToOrdinal_(2, iface, j);
            auto jdof = facesInternalDofOrdinals_(iface, j);
            for(ordinal_type d=0; d <cellDim_; ++d) {
              dofCoeffs_(cell,idof,d) = 0;
              for(ordinal_type k=0; k <2; ++k)
                for(ordinal_type l=0; l <2; ++l)
                  dofCoeffs_(cell,idof,d) += refFaceTangents_(iface, d,l)*ortJacobianFace(l,k)*faceDofCoeffs_(iface,jdof,k);
            }
          }
        }
      } else if(isBasisHDIV_) {
        for (ordinal_type iface=0; iface < numFaces_; ++iface) {
          Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobianFace, faceTopoKey_(iface), fOrt[iface]);
          auto ortJacobianDet = ortJacobianFace(0,0)*ortJacobianFace(1,1)-ortJacobianFace(1,0)*ortJacobianFace(0,1);
          for(ordinal_type j=0; j<numFacesInternalDofs_(iface); ++j) {
            auto idof = tagToOrdinal_(2, iface, j);
            auto jdof = facesInternalDofOrdinals_(iface, j);
            for(ordinal_type d=0; d <cellDim_; ++d)
              dofCoeffs_(cell,idof,d) = refFacesNormal_(iface,d)*ortJacobianDet*faceDofCoeffs_(iface,jdof)*faceScale_;
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
  typedef  Kokkos::DynRankView<scalarType, SpT> ScalarViewType;
  typedef  Kokkos::DynRankView<ordinal_type, SpT> intViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

  const auto topo = basis->getBaseCellTopology();
  const std::string name(basis->getName());
  const ordinal_type degree(basis->getDegree());

  bool isBasisHCURL = (basis->getFunctionSpace()==FUNCTION_SPACE_HCURL);
  bool isBasisHDIV = (basis->getFunctionSpace()==FUNCTION_SPACE_HDIV);
  bool isBasisTriOrTet = (basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Triangle<3> >()->key) ||
      (basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key);
  bool isBasisHEXI1 = (name.find("HEX_I1") != std::string::npos);
  bool isBasisTETI1 = (name.find("TET_I1") != std::string::npos);
  bool isBasisI1 = (name.find("I1") != std::string::npos);

  ordinal_type numEdges = (basis->getDofCount(1, 0) > 0) ? topo.getEdgeCount() : 0;
  ordinal_type numFaces = (basis->getDofCount(2, 0) > 0) ? topo.getFaceCount() : 0;

  Teuchos::RCP<Basis<SpT,scalarType,scalarType> > faceBases[6];
  Teuchos::RCP<Basis<SpT,scalarType,scalarType> > edgeBasis;

  bool unsupportedCase = false;
  switch(basis->getFunctionSpace()) {
  case FUNCTION_SPACE_HGRAD:
    if(basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Quadrilateral<4> >()->key) {
      edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree, pointType) );
    }
    else if(basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Hexahedron<8> >()->key){
      faceBases[0] = Teuchos::rcp(new Basis_HGRAD_QUAD_Cn_FEM<SpT,scalarType,scalarType>(degree, pointType) );
      for(ordinal_type i=1; i<6; ++i) faceBases[i]=faceBases[0];
      edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree, pointType) );
    }
    else if(basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Triangle<3> >()->key){
      edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree, pointType) );
    }
    else if(basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key) {
      faceBases[0] = Teuchos::rcp(new Basis_HGRAD_TRI_Cn_FEM<SpT,scalarType,scalarType>(degree, pointType) );
      for(ordinal_type i=1; i<4; ++i) faceBases[i]=faceBases[0];
      edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree, pointType) );
    }
    else
      unsupportedCase = true;
    break;
  case FUNCTION_SPACE_HCURL:
    if(basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Quadrilateral<4> >()->key) {
      edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree-1, POINTTYPE_GAUSS) );
    }
    else if(basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Hexahedron<8> >()->key) {
      faceBases[0] = Teuchos::rcp(new Basis_HCURL_QUAD_In_FEM<SpT,scalarType,scalarType>(degree, pointType) );
      for(ordinal_type i=1; i<6; ++i) faceBases[i]=faceBases[0];
      edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree-1, POINTTYPE_GAUSS) );
    }
    else if(basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Triangle<3> >()->key) {
      edgeBasis = Teuchos::rcp(new Basis_HVOL_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree-1, pointType) );
    }
    else if(basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key) {
      faceBases[0] = Teuchos::rcp(new Basis_HCURL_TRI_In_FEM<SpT,scalarType,scalarType>(degree, pointType) );
      for(ordinal_type i=1; i<4; ++i) faceBases[i]=faceBases[0];
      edgeBasis = Teuchos::rcp(new Basis_HVOL_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree-1, pointType) );
    }
    else
      unsupportedCase = true;
    break;
  case FUNCTION_SPACE_HDIV:
    isBasisHDIV = true;
    if(basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Quadrilateral<4> >()->key) {
      edgeBasis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree-1, POINTTYPE_GAUSS) );
    }
    else if(basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Hexahedron<8> >()->key) {
      faceBases[0] = Teuchos::rcp(new Basis_HGRAD_QUAD_Cn_FEM<SpT,scalarType,scalarType>(degree-1, POINTTYPE_GAUSS) );
      for(ordinal_type i=1; i<6; ++i) faceBases[i]=faceBases[0];
    }
    else if(basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Triangle<3> >()->key) {
      edgeBasis = Teuchos::rcp(new Basis_HVOL_LINE_Cn_FEM<SpT,scalarType,scalarType>(degree-1, pointType) );
    }
    else if(basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key) {
      faceBases[0] = Teuchos::rcp(new Basis_HVOL_TRI_Cn_FEM<SpT,scalarType,scalarType>(degree-1, pointType) );
      for(ordinal_type i=1; i<4; ++i) faceBases[i]=faceBases[0];
    }
    else
      unsupportedCase = true;
    break;
  case FUNCTION_SPACE_HVOL:
    //nothing to do
    break;
  default:
    unsupportedCase = true;
  }

  INTREPID2_TEST_FOR_ABORT(unsupportedCase,
            ">>> ERROR (Intrepid2::Experimental::LagrangianInterpolation<SpT>::getDofCoordsAndCoeffs): " \
            "method not implemented for this basis function");


/*

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
  */

  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(), basis->getAllDofOrdinal());

  const ordinal_type dim = topo.getDimension();

  const ordinal_type numCells = dofCoeffs.extent(0);

  ScalarViewType refDofCoords("refDofCoords", dofCoords.extent(1), dofCoords.extent(2)), refDofCoeffs;
  basis->getDofCoords(refDofCoords);
  RealSpaceTools<SpT>::clone(dofCoords,refDofCoords);

  if(dofCoeffs.rank() == 3) //vector basis
    refDofCoeffs = ScalarViewType("refDofCoeffs", dofCoeffs.extent(1), dofCoeffs.extent(2));
  else //scalar basis
    refDofCoeffs = ScalarViewType("refDofCoeffs",dofCoeffs.extent(1));
  basis->getDofCoeffs(refDofCoeffs);
  RealSpaceTools<SpT>::clone(dofCoeffs,refDofCoeffs);

  //*** Pre-compute needed quantities related to edge DoFs that do not depend on the cell ***

  ordinal_type edgeTopoKey = Teuchos::nonnull(edgeBasis) ? edgeBasis->getBaseCellTopology().getBaseKey() : 0;
  intViewType eOrt("eOrt", numEdges);
  ScalarViewType refEdgesTan("refEdgesTan",  numEdges, dim);
  ScalarViewType refEdgesNormal("refEdgesNormal",  numEdges, dim);
  ScalarViewType edgeParam;
  ordinal_type edgeBasisCardinality = Teuchos::nonnull(edgeBasis) ? edgeBasis->getCardinality() : ordinal_type(0);
  ordinal_type numEdgeInternalDofs = Teuchos::nonnull(edgeBasis) ? edgeBasis->getDofCount(1,0) : ordinal_type(0);
  ScalarViewType edgeDofCoords("edgeDofCoords", edgeBasisCardinality, 1);
  ScalarViewType edgeDofCoeffs("edgeDofCoeffs", edgeBasisCardinality);
  ScalarViewType edgeInternalDofCoords("edgeInternalDofCoords", numEdgeInternalDofs, 1);
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
  ScalarViewType refFaceTangents("refFaceTangents", numFaces, dim, 2);
  ScalarViewType refFacesNormal("refFacesNormal",  numFaces, dim);
  ScalarViewType faceParam;
  intViewType numFacesInternalDofs("numFacesInternalDofs", numFaces);
  ScalarViewType  facesInternalDofCoords;
  intViewType  facesInternalDofOrdinals;
  ScalarViewType  faceDofCoeffs;
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

  facesInternalDofCoords = ScalarViewType("faceInternalDofCoords", numFaces, maxNumFacesInternalDofs, 2);
  facesInternalDofOrdinals = intViewType("faceInternalDofCoords", numFaces, maxNumFacesInternalDofs);

  if(isBasisHCURL)
    faceDofCoeffs = ScalarViewType("faceDofCoeffs", numFaces, faceBasisMaxCardinality,2);
  else
    faceDofCoeffs = ScalarViewType("faceDofCoeffs", numFaces, faceBasisMaxCardinality);

  for (ordinal_type iface=0; iface < numFaces; ++iface) {
    auto faceBasis = faceBases[iface];
    faceTopoKey(iface) = faceBasis->getBaseCellTopology().getBaseKey();
    ordinal_type faceBasisCardinality = faceBasis->getCardinality();
    ScalarViewType  faceDofCoords("faceDofCoords", faceBasisCardinality, 2);
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
  typedef computeDofCoordsAndCoeffs
      <ScalarViewType,
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

