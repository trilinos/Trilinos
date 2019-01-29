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


  auto ordinalToTag = basis->getAllDofTags();
  auto tagToOrdinal = basis->getAllDofOrdinal();

  const ordinal_type dim = basis->getBaseCellTopology().getDimension();

  const ordinal_type numCells = dofCoeffs.extent(0);

  intViewType eOrt("eOrt", numEdges);
  intViewType fOrt("fOrt", numFaces);
  scalarViewType refEdgeTan("refEdgeTan",  dim);
  scalarViewType refFaceTangents("refFaceTangents", dim, 2);
  auto refFaceTanU = Kokkos::subview(refFaceTangents, Kokkos::ALL, 0);
  auto refFaceTanV = Kokkos::subview(refFaceTangents, Kokkos::ALL, 1);
  scalarViewType refNormal("refNormal",  dim);


  scalarViewType refDofCoords("refDofCoords", dofCoords.extent(1), dofCoords.extent(2)), refDofCoeffs;
  basis->getDofCoords(refDofCoords);
  RealSpaceTools<SpT>::clone(dofCoords,refDofCoords);

  if(dofCoeffs.rank() == 3) //vector basis
    refDofCoeffs = scalarViewType("refDofCoeffs", dofCoeffs.extent(1), dofCoeffs.extent(2));
  else //scalar basis
    refDofCoeffs = scalarViewType("refDofCoeffs",dofCoeffs.extent(1));
  basis->getDofCoeffs(refDofCoeffs);
  RealSpaceTools<SpT>::clone(dofCoeffs,refDofCoeffs);


  //compute DofCoords Oriented
  for (ordinal_type iedge=0; iedge < numEdges; ++iedge) {
    scalarViewType ortJacobian("ortJacobian", 1, 1);
    ordinal_type edgeBasisCardinality = edgeBasis->getCardinality();
    ordinal_type numInternalDofs = edgeBasis->getDofCount(1,0);
    scalarViewType edgeDofCoords("edgeDofCoords", edgeBasisCardinality, 1);
    scalarViewType edgeInternalDofCoords("edgeInternalDofCoords", numInternalDofs, 1);
    edgeBasis->getDofCoords(edgeDofCoords);
    for(ordinal_type i=0; i<numInternalDofs; ++i)
      edgeInternalDofCoords(i,0) = edgeDofCoords(edgeBasis->getDofOrdinal(1, 0, i),0);

    scalarViewType  edgeInternalDofCoordsOriented("edgeInternalDofCoordsOriented",  numInternalDofs, 1);
    scalarViewType edgeDofCoordsOriented3d("edgeDofCoordsOriented3d",  numInternalDofs, dim);
    // auto numEdgeDOFs = basis->getDofCount(1,iedge);
    for(ordinal_type i=0; i<numCells; ++i) {
      orts(i).getEdgeOrientation(eOrt.data(), numEdges);
      Impl::OrientationTools::mapToModifiedReference(edgeInternalDofCoordsOriented,edgeInternalDofCoords,edgeBasis->getBaseCellTopology(),eOrt(iedge));
      CellTools<SpT>::mapToReferenceSubcell(edgeDofCoordsOriented3d, edgeInternalDofCoordsOriented, 1, iedge, topo);

      for(ordinal_type j=0; j<numInternalDofs; ++j) {
        auto idof = basis->getDofOrdinal(1, iedge, j);
        for(ordinal_type d=0; d <dim; ++d)
          dofCoords(i,idof,d) = edgeDofCoordsOriented3d(j,d);
      }
    }

    //depending on how the reference basis is defined, the edges are scaled differently
    auto edgeScale = (isBasisTriOrTet||isBasisI1) ? 2.0 :1.0;
    if(isBasisHCURL) {

      scalarViewType edgeDofCoeffs("edgeDofCoeffs", edgeBasisCardinality);
      edgeBasis->getDofCoeffs(edgeDofCoeffs);
      CellTools<SpT>::getReferenceEdgeTangent(refEdgeTan, iedge, topo);
      for(ordinal_type i=0; i<numCells; ++i) {
        orts(i).getEdgeOrientation(eOrt.data(), numEdges);
        Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, edgeBasis->getBaseCellTopology(), eOrt(iedge));
        for(ordinal_type j=0; j<numInternalDofs; ++j) {
          auto idof = basis->getDofOrdinal(1, iedge, j);
          auto jdof = edgeBasis->getDofOrdinal(1, 0, j);
          for(ordinal_type d=0; d <dim; ++d) {
            dofCoeffs(i,idof,d) = 0;
            for(ordinal_type k=0; k <1; ++k)
              for(ordinal_type l=0; l <1; ++l)
                dofCoeffs(i,idof,d) += refEdgeTan(d)*ortJacobian(0,0)*edgeDofCoeffs(jdof)*edgeScale;
          }
        }
      }
    } else if(isBasisHDIV) {
      CellTools<SpT>::getReferenceSideNormal(refNormal, iedge, topo);
      scalarViewType edgeDofCoeffs("edgeDofCoeffs", edgeBasisCardinality);
      edgeBasis->getDofCoeffs(edgeDofCoeffs);
      for(ordinal_type i=0; i<numCells; ++i) {
        orts(i).getEdgeOrientation(eOrt.data(), numEdges);
        Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, edgeBasis->getBaseCellTopology(), eOrt(iedge));
        auto ortJacobianDet = ortJacobian(0,0);
        for(ordinal_type j=0; j<numInternalDofs; ++j) {
          auto idof = basis->getDofOrdinal(1, iedge, j);
          auto jdof = edgeBasis->getDofOrdinal(1, 0, j);
          for(ordinal_type d=0; d <dim; ++d)
            dofCoeffs(i,idof,d) = refNormal(d)*ortJacobianDet*edgeDofCoeffs(jdof)*edgeScale;
        }
      }
    } else {
      scalarViewType edgeDofCoeffs("edgeDofCoeffs", edgeBasisCardinality);
      edgeBasis->getDofCoeffs(edgeDofCoeffs);
      for(ordinal_type i=0; i<numCells; ++i) {
        for(ordinal_type j=0; j<numInternalDofs; ++j) {
          auto idof = basis->getDofOrdinal(1, iedge, j);
          auto jdof = edgeBasis->getDofOrdinal(1, 0, j);
          dofCoeffs(i,idof,0) = edgeDofCoeffs(jdof);
        }
      }
    }
  }

  for (ordinal_type iface=0; iface < numFaces; ++iface) {
    scalarViewType ortJacobian("ortJacobian", 2, 2);
    auto faceBasis = faceBases[iface];
    ordinal_type faceBasisCardinality = faceBasis->getCardinality();
    ordinal_type numInternalDofs = faceBasis->getDofCount(2,0);
    scalarViewType  faceDofCoords("faceDofCoords", faceBasisCardinality, 2);
    scalarViewType  faceInternalDofCoords("faceInternalDofCoords", numInternalDofs, 2);
    faceBasis->getDofCoords(faceDofCoords);
    for(ordinal_type i=0; i<numInternalDofs; ++i)
      for(ordinal_type d=0; d <2; ++d)
        faceInternalDofCoords(i,d) = faceDofCoords(faceBasis->getDofOrdinal(2, 0, i),d);

    scalarViewType  faceInternalDofCoordsOriented("faceInternalDofCoordsOriented",  numInternalDofs, 2);
    scalarViewType  faceDofCoordsOriented3d("faceDofCoordsOriented3d", numInternalDofs, dim);
    for(ordinal_type i=0; i<numCells; ++i) {
      orts(i).getFaceOrientation(fOrt.data(), numFaces);
      ordinal_type ort = fOrt(iface);
      Impl::OrientationTools::mapToModifiedReference(faceInternalDofCoordsOriented,faceInternalDofCoords,faceBasis->getBaseCellTopology(),ort);
      CellTools<SpT>::mapToReferenceSubcell(faceDofCoordsOriented3d, faceInternalDofCoordsOriented, 2, iface, topo);

      for(ordinal_type j=0; j<numInternalDofs; ++j) {
        auto idof = basis->getDofOrdinal(2, iface, j);
        for(ordinal_type d=0; d <dim; ++d)
          dofCoords(i,idof,d) = faceDofCoordsOriented3d(j,d);
      }
    }

    if(isBasisHCURL) {
      scalarViewType faceDofCoeffs("faceDofCoeffs", faceBasisCardinality, 2);
      faceBasis->getDofCoeffs(faceDofCoeffs);
      CellTools<SpT>::getReferenceFaceTangents(refFaceTanU, refFaceTanV, iface, topo);
      for(ordinal_type i=0; i<numCells; ++i) {
        orts(i).getFaceOrientation(fOrt.data(), numFaces);
        Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, faceBasis->getBaseCellTopology(), fOrt(iface));
        for(ordinal_type j=0; j<numInternalDofs; ++j) {
          auto idof = basis->getDofOrdinal(2, iface, j);
          auto jdof = faceBasis->getDofOrdinal(2, 0, j);
          for(ordinal_type d=0; d <dim; ++d) {
            dofCoeffs(i,idof,d) = 0;
            for(ordinal_type k=0; k <2; ++k)
              for(ordinal_type l=0; l <2; ++l)
                dofCoeffs(i,idof,d) += refFaceTangents(d,l)*ortJacobian(l,k)*faceDofCoeffs(jdof,k);
          }
        }
      }
    } else if(isBasisHDIV) {
      //depending on how the reference basis is defined, the faces are scaled differently
      auto faceScale = (isBasisHEXI1) ? 4.0 :
                       (isBasisTETI1) ? 0.5 : 1.0;
      CellTools<SpT>::getReferenceFaceNormal(refNormal, iface, topo);
      scalarViewType faceDofCoeffs("faceDofCoeffs", faceBasisCardinality);
      faceBasis->getDofCoeffs(faceDofCoeffs);
      for(ordinal_type i=0; i<numCells; ++i) {
        orts(i).getFaceOrientation(fOrt.data(), numFaces);
        Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, faceBasis->getBaseCellTopology(), fOrt(iface));
        auto ortJacobianDet = ortJacobian(0,0)*ortJacobian(1,1)-ortJacobian(1,0)*ortJacobian(0,1);
        for(ordinal_type j=0; j<numInternalDofs; ++j) {
          auto idof = basis->getDofOrdinal(2, iface, j);
          auto jdof = faceBasis->getDofOrdinal(2, 0, j);
          for(ordinal_type d=0; d <dim; ++d)
            dofCoeffs(i,idof,d) = refNormal(d)*ortJacobianDet*faceDofCoeffs(jdof)*faceScale;
        }
      }
    } else {
      scalarViewType faceDofCoeffs("faceDofCoeffs", faceBasisCardinality);
      faceBasis->getDofCoeffs(faceDofCoeffs);
      for(ordinal_type i=0; i<numCells; ++i) {
        for(ordinal_type j=0; j<numInternalDofs; ++j) {
          auto idof = basis->getDofOrdinal(2, iface, j);
          auto jdof = faceBasis->getDofOrdinal(2, 0, j);
          dofCoeffs(i,idof,0) = faceDofCoeffs(jdof);
        }
      }
    }
  }

  auto numElemDOFs = basis->getDofCount(dim,0);
  for(ordinal_type j=0; j<numElemDOFs; ++j) {
    auto idof = basis->getDofOrdinal(dim, 0, j);
    for(ordinal_type d=0; d <dim; ++d)
      for(ordinal_type i=1; i <numCells; ++i)
        dofCoords(i,idof,d) = refDofCoords(idof,d);
    for(ordinal_type d=0; d <(ordinal_type)dofCoeffs.extent(2); ++d)
      for(ordinal_type i=1; i <numCells; ++i)
        dofCoeffs(i,idof,d) = refDofCoeffs(idof,d);
  }

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

