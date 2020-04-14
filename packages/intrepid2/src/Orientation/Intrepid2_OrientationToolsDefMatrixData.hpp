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


/** \file   Intrepid2_OrientationToolsDefMatrixData.hpp
    \brief  Definition file for matrix data in the Intrepid2::OrientationTools class.
    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_MATRIX_DATA_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_MATRIX_DATA_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  template<typename SpT>
  template<typename BasisType>
  typename OrientationTools<SpT>::CoeffMatrixDataViewType
  OrientationTools<SpT>::createCoeffMatrixInternal(const BasisType* basis) {
    const ordinal_type order(basis->getDegree());
    const std::string name(basis->getName());
    CoeffMatrixDataViewType matData;

    using ExecutionSpace  = typename BasisType::ExecutionSpace;
    using OutputValueType = typename BasisType::OutputValueType;
    using PointValueType  = typename BasisType::PointValueType;

    using hostExecutionSpace = Kokkos::DefaultHostExecutionSpace;

    auto ordinalToTag = basis->getAllDofTags();
    auto tagToOrdinal = basis->getAllDofOrdinal();

    //
    // High order HGRAD Elements
    //

    //std::cout << "Name: " << name << " " << order << std::endl;

    if((basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Quadrilateral<4> >()->key)
      && (basis->getFunctionSpace() == FUNCTION_SPACE_HGRAD)) {
      if (order >1) {
        const ordinal_type matDim = ordinalToTag(tagToOrdinal(1, 0, 0), 3), numEdges = 4, numOrts = 2;
        matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::"+name,
                                          numEdges,
                                          numOrts,
                                          matDim, 
                                          matDim);

        if(dynamic_cast<const typename NodalBasisFamily<ExecutionSpace, OutputValueType, PointValueType>::HGRAD_QUAD*>(basis)) {
            typename NodalBasisFamily<hostExecutionSpace>::HGRAD_QUAD hostBasis(order);
            init_HGRAD_QUAD(matData, &hostBasis);
        } else if(dynamic_cast<const typename DerivedNodalBasisFamilyModified<ExecutionSpace, OutputValueType, PointValueType>::HGRAD_QUAD*>(basis)) {
            typename DerivedNodalBasisFamilyModified<hostExecutionSpace>::HGRAD_QUAD hostBasis(order);
            init_HGRAD_QUAD(matData, &hostBasis);
        } else if(dynamic_cast<const typename HierarchicalBasisFamily<ExecutionSpace, OutputValueType, PointValueType>::HGRAD_QUAD*>(basis)) {
          typename HierarchicalBasisFamily<hostExecutionSpace>::HGRAD_QUAD hostBasis(order);
          init_HGRAD_QUAD(matData, &hostBasis);
        }
      } else {
        // add dummy
        matData = CoeffMatrixDataViewType();
      }
    }
    else if((basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Hexahedron<8> >()->key)
        && (basis->getFunctionSpace() == FUNCTION_SPACE_HGRAD)) {
      if (order > 1) {
        const ordinal_type matDim = ordinalToTag(tagToOrdinal(2, 0, 0), 3), numSubcells = 18, numOrts = 8;
        matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::Intrepid2_HGRAD_HEX_Cn_FEM",
                                          numSubcells,
                                          numOrts,
                                          matDim, 
                                          matDim);
        
        if(dynamic_cast<const typename NodalBasisFamily<ExecutionSpace, OutputValueType, PointValueType>::HGRAD_HEX*>(basis)) {
            typename NodalBasisFamily<hostExecutionSpace>::HGRAD_HEX hostBasis(order);
            init_HGRAD_HEX(matData, &hostBasis);
        } else if(dynamic_cast<const typename DerivedNodalBasisFamilyModified<ExecutionSpace, OutputValueType, PointValueType>::HGRAD_HEX*>(basis)) {
            typename DerivedNodalBasisFamilyModified<hostExecutionSpace>::HGRAD_HEX hostBasis(order);
            init_HGRAD_HEX(matData, &hostBasis);
        } else if(dynamic_cast<const typename HierarchicalBasisFamily<ExecutionSpace, OutputValueType, PointValueType>::HGRAD_HEX*>(basis)) {
          typename HierarchicalBasisFamily<hostExecutionSpace>::HGRAD_HEX hostBasis(order);
          init_HGRAD_HEX(matData, &hostBasis);
        }
      } else {
        // add dummy
        matData = CoeffMatrixDataViewType();
      }
    }
    else if((basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Triangle<3> >()->key)
        && (basis->getFunctionSpace() == FUNCTION_SPACE_HGRAD)) {
      if (order > 1) {
        const ordinal_type matDim = ordinalToTag(tagToOrdinal(1, 0, 0), 3), numEdges = 3, numOrts = 2;      
        matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::Intrepid2_HGRAD_TRI_Cn_FEM",
                                          numEdges,
                                          numOrts,
                                          matDim, 
                                          matDim);
        
        Basis_HGRAD_TRI_Cn_FEM<hostExecutionSpace> hostBasis(order);
        init_HGRAD_TRI(matData, &hostBasis);
      } else {
        // add dummy
        matData = CoeffMatrixDataViewType();
      }        
    } 
    else if((basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key)
        && (basis->getFunctionSpace() == FUNCTION_SPACE_HGRAD)) {
      const ordinal_type 
        matDim = max(ordinalToTag(tagToOrdinal(1, 0, 0), 3),
                     ordinalToTag(tagToOrdinal(2, 0, 0), 3)),
        numSubcells = 10, numOrts = 6;
      matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::Intrepid2_HGRAD_TET_Cn_FEM",
                                        numSubcells,
                                        numOrts,
                                        matDim, 
                                        matDim);
      
      Basis_HGRAD_TET_Cn_FEM<hostExecutionSpace> hostBasis(order);
      init_HGRAD_TET(matData, &hostBasis);
    } 

    //
    // High order HCURL Elements
    //

    else if((basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Quadrilateral<4> >()->key)
        && (basis->getFunctionSpace() == FUNCTION_SPACE_HCURL)) {
      const ordinal_type matDim = ordinalToTag(tagToOrdinal(1, 0, 0), 3), numEdges = 4, numOrts = 2;
      matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::Intrepid2_HCURL_QUAD_In_FEM",
                                        numEdges,
                                        numOrts,
                                        matDim, 
                                        matDim);

      if (name == "Intrepid2_HCURL_QUAD_I1_FEM") {
        for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId)
          init_EDGE_ELEMENT_I1_FEM(matData, edgeId);
      } else if(dynamic_cast<const typename NodalBasisFamily<ExecutionSpace, OutputValueType, PointValueType>::HCURL_QUAD*>(basis)) {
          typename NodalBasisFamily<hostExecutionSpace>::HCURL_QUAD hostBasis(order);
          init_HCURL_QUAD(matData, &hostBasis);
      } else if(dynamic_cast<const typename DerivedNodalBasisFamilyModified<ExecutionSpace, OutputValueType, PointValueType>::HCURL_QUAD*>(basis)) {
          typename DerivedNodalBasisFamilyModified<hostExecutionSpace>::HCURL_QUAD hostBasis(order);
          init_HCURL_QUAD(matData, &hostBasis);
      } else if(dynamic_cast<const typename HierarchicalBasisFamily<ExecutionSpace, OutputValueType, PointValueType>::HCURL_QUAD*>(basis)) {
        typename HierarchicalBasisFamily<hostExecutionSpace>::HCURL_QUAD hostBasis(order);
        init_HCURL_QUAD(matData, &hostBasis);
      }
    } 
    else if((basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Hexahedron<8> >()->key)
        && (basis->getFunctionSpace() == FUNCTION_SPACE_HCURL)) {
      // if order == 1, there is no face dofs
      const ordinal_type matDim = ordinalToTag(tagToOrdinal((order < 2 ? 1 : 2), 0, 0), 3), numSubcells = (order < 2 ? 12 : 18), numOrts = 8;
      matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::Intrepid2_HCURL_HEX_In_FEM",
                                        numSubcells,
                                        numOrts,
                                        matDim, 
                                        matDim);

      if (name == "Intrepid2_HCURL_HEX_I1_FEM") {
        for (ordinal_type edgeId=0;edgeId<numSubcells;++edgeId)
          init_EDGE_ELEMENT_I1_FEM(matData, edgeId);
      } else if(dynamic_cast<const typename NodalBasisFamily<ExecutionSpace, OutputValueType, PointValueType>::HCURL_HEX*>(basis)) {
          typename NodalBasisFamily<hostExecutionSpace>::HCURL_HEX hostBasis(order);
          init_HCURL_HEX(matData, &hostBasis);
      } else if(dynamic_cast<const typename DerivedNodalBasisFamilyModified<ExecutionSpace, OutputValueType, PointValueType>::HCURL_HEX*>(basis)) {
          typename DerivedNodalBasisFamilyModified<hostExecutionSpace>::HCURL_HEX hostBasis(order);
          init_HCURL_HEX(matData, &hostBasis);
      } else if(dynamic_cast<const typename HierarchicalBasisFamily<ExecutionSpace, OutputValueType, PointValueType>::HCURL_HEX*>(basis)) {
        typename HierarchicalBasisFamily<hostExecutionSpace>::HCURL_HEX hostBasis(order);
        init_HCURL_HEX(matData, &hostBasis);
      }
    } 
    else if((basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Triangle<3> >()->key)
        && (basis->getFunctionSpace() == FUNCTION_SPACE_HCURL)) {
      const ordinal_type matDim = ordinalToTag(tagToOrdinal(1, 0, 0), 3), numEdges = 3, numOrts = 2;
      matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::Intrepid2_HCURL_TRI_In_FEM",
                                        numEdges,
                                        numOrts,
                                        matDim, 
                                        matDim);

      if (name == "Intrepid2_HCURL_TRI_I1_FEM") {
        for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId)
          init_EDGE_ELEMENT_I1_FEM(matData, edgeId);
      } else if (name == "Intrepid2_HCURL_TRI_In_FEM") {
        Basis_HCURL_TRI_In_FEM<hostExecutionSpace> hostBasis(order);
        init_HCURL_TRI(matData, &hostBasis);
      }
    } 
    else if((basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key)
        && (basis->getFunctionSpace() == FUNCTION_SPACE_HCURL)) {
      const ordinal_type matDim = ordinalToTag(tagToOrdinal((order < 2 ? 1 : 2), 0, 0), 3), numSubcells = (order < 2 ? 6 : 10), numOrts = 6;
      matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::Intrepid2_HCURL_TET_In_FEM",
                                        numSubcells,
                                        numOrts,
                                        matDim, 
                                        matDim);

      if (name == "Intrepid2_HCURL_TET_I1_FEM") {
        for (ordinal_type edgeId=0;edgeId<numSubcells;++edgeId)
          init_EDGE_ELEMENT_I1_FEM(matData, edgeId);
      } else if (name == "Intrepid2_HCURL_TET_In_FEM") {
        Basis_HCURL_TET_In_FEM<hostExecutionSpace> hostBasis(order);
        init_HCURL_TET(matData, &hostBasis);
      }
    } 

    //
    // High order HDIV Elements
    //

    else if((basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Quadrilateral<4> >()->key)
        && (basis->getFunctionSpace() == FUNCTION_SPACE_HDIV)) {
            // (name == "Intrepid2_HCURL_QUAD_I1_FEM")) {
      const ordinal_type matDim = ordinalToTag(tagToOrdinal(1, 0, 0), 3), numEdges = 4, numOrts = 2;
      matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::Intrepid2_HCURL_QUAD_In_FEM",
                                        numEdges,
                                        numOrts,
                                        matDim, 
                                        matDim);

      if (name == "Intrepid2_HDIV_QUAD_I1_FEM") {
        for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId)
          init_EDGE_ELEMENT_I1_FEM(matData, edgeId);
      } else if(dynamic_cast<const typename NodalBasisFamily<ExecutionSpace, OutputValueType, PointValueType>::HDIV_QUAD*>(basis)) {
          typename NodalBasisFamily<hostExecutionSpace>::HDIV_QUAD hostBasis(order);
          init_HDIV_QUAD(matData, &hostBasis);
      } else if(dynamic_cast<const typename DerivedNodalBasisFamilyModified<ExecutionSpace, OutputValueType, PointValueType>::HDIV_QUAD*>(basis)) {
          typename DerivedNodalBasisFamilyModified<hostExecutionSpace>::HDIV_QUAD hostBasis(order);
          init_HDIV_QUAD(matData, &hostBasis);
      } else if(dynamic_cast<const typename HierarchicalBasisFamily<ExecutionSpace, OutputValueType, PointValueType>::HDIV_QUAD*>(basis)) {
        typename HierarchicalBasisFamily<hostExecutionSpace>::HDIV_QUAD hostBasis(order);
        init_HDIV_QUAD(matData, &hostBasis);
      }
    } 
    else if((basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Hexahedron<8> >()->key)
        && (basis->getFunctionSpace() == FUNCTION_SPACE_HDIV)) {
      const ordinal_type matDim = ordinalToTag(tagToOrdinal(2, 0, 0), 3), numSubcells = 6, numOrts = 8;
      matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::Intrepid2_HDIV_HEX_In_FEM",
                                        numSubcells,
                                        numOrts,
                                        matDim, 
                                        matDim);

      if (name == "Intrepid2_HDIV_HEX_I1_FEM") {
        for (ordinal_type faceId=0;faceId<numSubcells;++faceId)
          init_QUAD_FACE_ELEMENT_I1_FEM(matData, faceId);
      } else if(dynamic_cast<const typename NodalBasisFamily<ExecutionSpace, OutputValueType, PointValueType>::HDIV_HEX*>(basis)) {
          typename NodalBasisFamily<hostExecutionSpace>::HDIV_HEX hostBasis(order);
          init_HDIV_HEX(matData, &hostBasis);
      } else if(dynamic_cast<const typename DerivedNodalBasisFamilyModified<ExecutionSpace, OutputValueType, PointValueType>::HDIV_HEX*>(basis)) {
          typename DerivedNodalBasisFamilyModified<hostExecutionSpace>::HDIV_HEX hostBasis(order);
          init_HDIV_HEX(matData, &hostBasis);
      } else if(dynamic_cast<const typename HierarchicalBasisFamily<ExecutionSpace, OutputValueType, PointValueType>::HDIV_HEX*>(basis)) {
        typename HierarchicalBasisFamily<hostExecutionSpace>::HDIV_HEX hostBasis(order);
        init_HDIV_HEX(matData, &hostBasis);
      }
    } 
    else if((basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Triangle<3> >()->key)
        && (basis->getFunctionSpace() == FUNCTION_SPACE_HDIV)) {
      const ordinal_type matDim = ordinalToTag(tagToOrdinal(1, 0, 0), 3), numEdges = 3, numOrts = 2;
      matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::Intrepid2_HDIV_TRI_In_FEM",
                                        numEdges,
                                        numOrts,
                                        matDim, 
                                        matDim);
      if (name == "Intrepid2_HDIV_TRI_I1_FEM") {
        for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId)
          init_EDGE_ELEMENT_I1_FEM(matData, edgeId);
      } else if (name == "Intrepid2_HDIV_TRI_In_FEM") {
        Basis_HDIV_TRI_In_FEM<hostExecutionSpace> hostBasis(order);
        init_HDIV_TRI(matData, &hostBasis);
      }
    } 
    else if((basis->getBaseCellTopology().getKey() == shards::getCellTopologyData<shards::Tetrahedron<4> >()->key)
        && (basis->getFunctionSpace() == FUNCTION_SPACE_HDIV)) {
      const ordinal_type matDim = ordinalToTag(tagToOrdinal(2, 0, 0), 3), numSubcells = 4, numOrts = 6;
      matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::Intrepid2_HDIV_TET_In_FEM",
                                        numSubcells,
                                        numOrts,
                                        matDim, 
                                        matDim);
      
      if (name == "Intrepid2_HDIV_TET_I1_FEM") {
        for (ordinal_type faceId=0;faceId<numSubcells;++faceId)
          init_TRI_FACE_ELEMENT_I1_FEM(matData, faceId);
      } else if (name == "Intrepid2_HDIV_TET_In_FEM") {
        Basis_HDIV_TET_In_FEM<hostExecutionSpace> hostBasis(order);
        init_HDIV_TET(matData, &hostBasis);
      }
    } 

    //
    // 3D H(Curl) I1 Elements
    //

    else if (name == "Intrepid2_HCURL_WEDGE_I1_FEM") {
      const ordinal_type matDim = 1, numEdges = 9, numOrts = 2;
      matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::Intrepid2_HCURL_WEDGE_I1_FEM",
                                        numEdges,
                                        numOrts,
                                        matDim, 
                                        matDim);
      for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId)
        init_EDGE_ELEMENT_I1_FEM(matData, edgeId);
    }

    //
    // 3D H(Div) I1 Elements
    //

    else if (name == "Intrepid2_HDIV_WEDGE_I1_FEM") {
      const ordinal_type matDim = 1, numFaces = 5, numOrts = 8;
      matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::Intrepid2_HDIV_WEDGE_I1_FEM",
                                        numFaces,
                                        numOrts,
                                        matDim, 
                                        matDim);
      ordinal_type faceId = 0;
      for ( ;faceId<3;++faceId) 
        init_QUAD_FACE_ELEMENT_I1_FEM(matData, faceId);
      for ( ;faceId<numFaces;++faceId) 
        init_TRI_FACE_ELEMENT_I1_FEM(matData, faceId);
    }
    return matData;
  }

  //
  // Quad elements
  //
  
  template<typename SpT>
  template<typename BasisType>
  void
  OrientationTools<SpT>::
  init_HGRAD_QUAD(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
                       BasisType const *cellBasis) {
    const ordinal_type order(cellBasis->getDegree());
    Basis_HGRAD_LINE_Cn_FEM<typename BasisType::ExecutionSpace> lineBasis(order);
    
    const ordinal_type numEdge = 4, numOrt = 2;
    for (ordinal_type edgeId=0;edgeId<numEdge;++edgeId)
      for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
        auto mat = Kokkos::subview(matData, 
                                   edgeId, edgeOrt,
                                   Kokkos::ALL(), Kokkos::ALL());
        Impl::OrientationTools::getCoeffMatrix_HGRAD(mat,
                                                     lineBasis, *cellBasis,
                                                     edgeId, edgeOrt);
      }
  }

  template<typename SpT>
  template<typename BasisType>
  void
  OrientationTools<SpT>::
  init_HCURL_QUAD(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
                       BasisType const *cellBasis) {
    const ordinal_type order(cellBasis->getDegree());
    Basis_HGRAD_LINE_Cn_FEM<typename BasisType::ExecutionSpace> bubbleBasis(order-1, POINTTYPE_GAUSS);
    
    const ordinal_type numEdge = 4, numOrt = 2;
    for (ordinal_type edgeId=0;edgeId<numEdge;++edgeId)
      for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
        auto mat = Kokkos::subview(matData, 
                                   edgeId, edgeOrt,
                                   Kokkos::ALL(), Kokkos::ALL());
        Impl::OrientationTools::getCoeffMatrix_HCURL(mat,
                                                     bubbleBasis, *cellBasis,
                                                     edgeId, edgeOrt);
      }
  }

  template<typename SpT>
  template<typename BasisType>
  void
  OrientationTools<SpT>::
  init_HDIV_QUAD(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
                          BasisType const *cellBasis) {
    const ordinal_type order(cellBasis->getDegree());
    Basis_HGRAD_LINE_Cn_FEM<typename BasisType::ExecutionSpace> bubbleBasis(order-1, POINTTYPE_GAUSS);
    
    const ordinal_type numEdge = 4, numOrt = 2;
    for (ordinal_type edgeId=0;edgeId<numEdge;++edgeId)
      for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
        auto mat = Kokkos::subview(matData, 
                                   edgeId, edgeOrt,
                                   Kokkos::ALL(), Kokkos::ALL());
        Impl::OrientationTools::getCoeffMatrix_HDIV(mat,
                                                    bubbleBasis, *cellBasis,
                                                    edgeId, edgeOrt);
      }
  }

  //
  // Hexahedral elements
  //

  template<typename SpT>
  template<typename BasisType>
  void
  OrientationTools<SpT>::
  init_HGRAD_HEX(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
                                 BasisType const *cellBasis) {
      const ordinal_type order(cellBasis->getDegree());
      Basis_HGRAD_LINE_Cn_FEM<typename BasisType::ExecutionSpace> lineBasis(order);
      Basis_HGRAD_QUAD_Cn_FEM<typename BasisType::ExecutionSpace> quadBasis(order);
      
      const ordinal_type numEdge = 12, numFace = 6;    
      {
        const ordinal_type numOrt = 2;
        for (ordinal_type edgeId=0;edgeId<numEdge;++edgeId)
          for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
            auto mat = Kokkos::subview(matData, 
                                       edgeId, edgeOrt,
                                       Kokkos::ALL(), Kokkos::ALL());
            Impl::OrientationTools::getCoeffMatrix_HGRAD(mat,
                                                         lineBasis, *cellBasis,
                                                         edgeId, edgeOrt);
          }
      }
      {
        const ordinal_type numOrt = 8;
        for (ordinal_type faceId=0;faceId<numFace;++faceId)
          for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
            auto mat = Kokkos::subview(matData, 
                                       numEdge+faceId, faceOrt,
                                       Kokkos::ALL(), Kokkos::ALL());
            Impl::OrientationTools::getCoeffMatrix_HGRAD(mat,
                                                         quadBasis, *cellBasis,
                                                         faceId, faceOrt);
          }
      }
  }

  template<typename SpT>
  template<typename BasisType>
  void
  OrientationTools<SpT>::
  init_HCURL_HEX(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
      BasisType const *cellBasis) {
    const ordinal_type order(cellBasis->getDegree());
    Basis_HGRAD_LINE_Cn_FEM<typename BasisType::ExecutionSpace> bubbleBasis(order-1, POINTTYPE_GAUSS);
    Basis_HCURL_QUAD_In_FEM<typename BasisType::ExecutionSpace> quadBasis(order);

    const ordinal_type numEdge = 12, numFace = 6;    
    {
      const ordinal_type numOrt = 2;
      for (ordinal_type edgeId=0;edgeId<numEdge;++edgeId)
        for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
          auto mat = Kokkos::subview(matData, 
                                     edgeId, edgeOrt,
                                     Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HCURL(mat,
                                                       bubbleBasis, *cellBasis,
                                                       edgeId, edgeOrt);
        }
    }
    if (order > 1) {
      const ordinal_type numOrt = 8;
      for (ordinal_type faceId=0;faceId<numFace;++faceId)
        for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
          auto mat = Kokkos::subview(matData, 
                                     numEdge+faceId, faceOrt,
                                     Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HCURL(mat,
                                                       quadBasis, *cellBasis,
                                                       faceId, faceOrt);
        }
    }
  }

  template<typename SpT>
  template<typename BasisType>
  void
  OrientationTools<SpT>::
  init_HDIV_HEX(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
      BasisType const *cellBasis) {
    const ordinal_type order(cellBasis->getDegree());
    Basis_HGRAD_QUAD_Cn_FEM<typename BasisType::ExecutionSpace> quadBasis(order-1, POINTTYPE_GAUSS);

    const ordinal_type numFace = 6;    
    {
      const ordinal_type numOrt = 8;
      for (ordinal_type faceId=0;faceId<numFace;++faceId)
        for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
          auto mat = Kokkos::subview(matData, 
                                     faceId, faceOrt,
                                     Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HDIV(mat,
                                                      quadBasis, *cellBasis,
                                                      faceId, faceOrt);
        }
    }
  }
  
  //
  // Triangle elements
  //

  template<typename SpT>
  template<typename BasisType>
  void
  OrientationTools<SpT>::
  init_HGRAD_TRI(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
      BasisType const *cellBasis) {
      const ordinal_type order(cellBasis->getDegree());
      Basis_HGRAD_LINE_Cn_FEM<typename BasisType::ExecutionSpace> lineBasis(order);
      
      const ordinal_type numEdge = 3, numOrt = 2;
      for (ordinal_type edgeId=0;edgeId<numEdge;++edgeId)
        for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
          auto mat = Kokkos::subview(matData, 
                                     edgeId, edgeOrt,
                                     Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HGRAD(mat,
                                                       lineBasis, *cellBasis,
                                                       edgeId, edgeOrt);
        }
  }


  template<typename SpT>
  template<typename BasisType>
  void
  OrientationTools<SpT>::
  init_HCURL_TRI(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
      BasisType const *cellBasis) {
    const ordinal_type order(cellBasis->getDegree());
    Basis_HVOL_LINE_Cn_FEM<typename BasisType::ExecutionSpace> bubbleBasis(order-1);
    
    const ordinal_type numEdge = 3, numOrt = 2;
    for (ordinal_type edgeId=0;edgeId<numEdge;++edgeId)
      for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
        auto mat = Kokkos::subview(matData, 
                                   edgeId, edgeOrt,
                                   Kokkos::ALL(), Kokkos::ALL());
        Impl::OrientationTools::getCoeffMatrix_HCURL(mat,
                                                     bubbleBasis, *cellBasis,
                                                     edgeId, edgeOrt);
      }
  }
  
  template<typename SpT>
  template<typename BasisType>
  void
  OrientationTools<SpT>::
  init_HDIV_TRI(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
      BasisType const *cellBasis) {
    const ordinal_type order(cellBasis->getDegree());
    Basis_HVOL_LINE_Cn_FEM<typename BasisType::ExecutionSpace> bubbleBasis(order-1);
    
    const ordinal_type numEdge = 3, numOrt = 2;
    for (ordinal_type edgeId=0;edgeId<numEdge;++edgeId)
      for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
        auto mat = Kokkos::subview(matData, 
                                   edgeId, edgeOrt,
                                   Kokkos::ALL(), Kokkos::ALL());
        Impl::OrientationTools::getCoeffMatrix_HDIV(mat,
                                                    bubbleBasis, *cellBasis,
                                                    edgeId, edgeOrt);
      }
  }

  //
  // Tetrahedral elements
  //
  
  template<typename SpT>
  template<typename BasisType>
  void
  OrientationTools<SpT>::
  init_HGRAD_TET(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
      BasisType const *cellBasis) {
    const ordinal_type order(cellBasis->getDegree());
      Basis_HGRAD_LINE_Cn_FEM<typename BasisType::ExecutionSpace> lineBasis(order);
      Basis_HGRAD_TRI_Cn_FEM<typename BasisType::ExecutionSpace>  triBasis(order);
      
      const ordinal_type numEdge = 6, numFace = 4;    
      {
        const ordinal_type numOrt = 2;
        for (ordinal_type edgeId=0;edgeId<numEdge;++edgeId)
          for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
            auto mat = Kokkos::subview(matData, 
                                       edgeId, edgeOrt,
                                       Kokkos::ALL(), Kokkos::ALL());
            Impl::OrientationTools::getCoeffMatrix_HGRAD(mat,
                                                         lineBasis, *cellBasis,
                                                         edgeId, edgeOrt);
          }
      }
      if (order > 2) {
        const ordinal_type numOrt = 6;
        for (ordinal_type faceId=0;faceId<numFace;++faceId)
          for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
            auto mat = Kokkos::subview(matData, 
                                       numEdge+faceId, faceOrt,
                                       Kokkos::ALL(), Kokkos::ALL());
            Impl::OrientationTools::getCoeffMatrix_HGRAD(mat,
                                                         triBasis, *cellBasis,
                                                         faceId, faceOrt);
          }
      }
  }

  template<typename SpT>
  template<typename BasisType>
  void
  OrientationTools<SpT>::
  init_HCURL_TET(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
      BasisType const *cellBasis) {
    const ordinal_type order(cellBasis->getDegree());
    Basis_HVOL_LINE_Cn_FEM<typename BasisType::ExecutionSpace> bubbleBasis(order-1);
    Basis_HCURL_TRI_In_FEM<typename BasisType::ExecutionSpace> triBasis(order);
      
    const ordinal_type numEdge = 6, numFace = 4;    
    {
      const ordinal_type numOrt = 2;
      for (ordinal_type edgeId=0;edgeId<numEdge;++edgeId)
        for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
          auto mat = Kokkos::subview(matData, 
                                     edgeId, edgeOrt,
                                     Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HCURL(mat,
                                                       bubbleBasis, *cellBasis,
                                                       edgeId, edgeOrt);
        }
    }
    if (order > 1) {
      const ordinal_type numOrt = 6;
      for (ordinal_type faceId=0;faceId<numFace;++faceId)
        for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
          auto mat = Kokkos::subview(matData, 
                                     numEdge+faceId, faceOrt,
                                     Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HCURL(mat,
                                                       triBasis, *cellBasis,
                                                       faceId, faceOrt);
        }
    }
  }

  template<typename SpT>
  template<typename BasisType>
  void
  OrientationTools<SpT>::
  init_HDIV_TET(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
      BasisType const *cellBasis) {
    const ordinal_type order(cellBasis->getDegree());
    Basis_HVOL_TRI_Cn_FEM<typename BasisType::ExecutionSpace> triBasis(order-1);
    
    const ordinal_type numFace = 4, numOrt = 6;
    for (ordinal_type faceId=0;faceId<numFace;++faceId)
      for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
        auto mat = Kokkos::subview(matData, 
                                   faceId, faceOrt,
                                   Kokkos::ALL(), Kokkos::ALL());
        Impl::OrientationTools::getCoeffMatrix_HDIV(mat,
                                                    triBasis, *cellBasis,
                                                    faceId, faceOrt);
      }
  }
  
  //
  // Lower order I1 elements
  //
  
  template<typename SpT>
  void
  OrientationTools<SpT>::
  init_EDGE_ELEMENT_I1_FEM(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
                           const ordinal_type edgeId) {
    const ordinal_type numOrt = 2;
    const double edgeOrtCoeff[2] = { 1.0, -1.0 };
    for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
      auto mat = Kokkos::subview(matData, 
                                 edgeId, edgeOrt,
                                 Kokkos::ALL(), Kokkos::ALL());
      mat(0,0) = edgeOrtCoeff[edgeOrt];
    }
  }
  
  template<typename SpT>
  void
  OrientationTools<SpT>::
  init_TRI_FACE_ELEMENT_I1_FEM(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
                               const ordinal_type faceId) {
    const ordinal_type numOrt = 6;
    const double faceOrtCoeff[6] = {   1.0,  1.0,  1.0, 
                                      -1.0, -1.0, -1.0 };
    
    for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
      auto mat = Kokkos::subview(matData, 
                                 faceId, faceOrt,
                                 Kokkos::ALL(), Kokkos::ALL());
      mat(0,0) = faceOrtCoeff[faceOrt];
    }
  }

  template<typename SpT>
  void
  OrientationTools<SpT>::
  init_QUAD_FACE_ELEMENT_I1_FEM(typename OrientationTools<SpT>::CoeffMatrixDataViewType matData,
                                const ordinal_type faceId) {
    const ordinal_type numOrt = 8;
    const double faceOrtCoeff[8] = {   1.0,  1.0,  1.0,  1.0, 
                                      -1.0, -1.0, -1.0, -1.0 };
    
    for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
      auto mat = Kokkos::subview(matData, 
                                 faceId, faceOrt,
                                 Kokkos::ALL(), Kokkos::ALL());
      mat(0,0) = faceOrtCoeff[faceOrt];
    }
  }

  template<typename SpT>
  template<typename BasisType>
  typename OrientationTools<SpT>::CoeffMatrixDataViewType
  OrientationTools<SpT>::createCoeffMatrix(const BasisType* basis) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !basis->requireOrientation(), std::invalid_argument,
                                  ">>> ERROR (OrientationTools::createCoeffMatrix): basis does not require orientations." );
#endif
    Kokkos::push_finalize_hook( [=] {
      ortCoeffData=std::map<std::pair<std::string,ordinal_type>, typename OrientationTools<SpT>::CoeffMatrixDataViewType>();
    });

    const std::pair<std::string,ordinal_type> key(basis->getName(), basis->getDegree());
    const auto found = ortCoeffData.find(key);
    
    CoeffMatrixDataViewType matData;
    if (found == ortCoeffData.end()) {
      matData = createCoeffMatrixInternal(basis);
      ortCoeffData.insert(std::make_pair(key, matData));
    } else {
      matData = found->second;
    }
    
    return matData;
  }
  
  template<typename SpT>
  void OrientationTools<SpT>::clearCoeffMatrix() {
    ortCoeffData.clear();
  }
}

#endif


//   template<typename SpT>
//   void
//   OrientationTools<SpT>::CoeffMatrix::import(const OrientationTools<SpT>::DenseMatrix &b,
//                                                 const bool transpose) {
// #ifdef HAVE_INTREPID2_DEBUG
//     INTREPID2_TEST_FOR_ABORT( !( NumRows() == b.NumRows() && NumCols() == b.NumCols() ), std::invalid_argument,
//                                 ">>> ERROR (Intrepid::Orientation::CoeffMatrix::import): "
//                                 "Matrix dimensions are not matched");
// #endif
//     // count size
//     const SpT eps = 1.0e-8;
//     const ordinal_type nrows = (transpose ? b.NumCols() : b.NumRows());
//     const ordinal_type ncols = (transpose ? b.NumRows() : b.NumCols());
//     size_type nnz = b.countNumNonZeros(eps);
//     createInternalArrays(nrows, ncols, nnz);

//     // construct sparse array
//     nnz = 0;
//     for (ordinal_type i=0;i<nrows;++i) {
//       _ap(i) = nnz;
//       for (ordinal_type j=0;j<ncols;++j) {
//         const SpT val  = (transpose ? b.Value(j,i) : b.Value(i,j));
//         const SpT val2 = val*val;

//         // consider it as a nonzero entry
//         if (val2 > eps) {
//           _aj(nnz) = j;
//           _ax(nnz) = val;
//           ++nnz;
//         }
//       }
//     }
//     _ap(nrows) = nnz;
//   }

//   template<typename SpT>
//   std::ostream&
//   OrientationTools<SpT>::CoeffMatrix::showMe(std::ostream &os) const {
//     std::ofstream prec;
//     prec.copyfmt(os);

//     os.precision(3);

//     os << " -- OrientationTools::CoeffMatrix -- " << std::endl
//        << "    # of Rows          = " << _m << std::endl
//        << "    # of Cols          = " << _n << std::endl
//        << std::endl
//        << "    RowPtrArray length = " << _ap.extent(0) << std::endl
//        << "    ColArray    length = " << _aj.extent(0) << std::endl
//        << "    ValueArray  length = " << _ax.extent(0) << std::endl
//        << std::endl;

//     const ordinal_type w = 10;
//     if (_ap.size() && _aj.size() && _ax.size()) {
//       os << std::setw(w) <<  "Row" << "  "
//          << std::setw(w) <<  "Col" << "  "
//          << std::setw(w) <<  "Val" << std::endl;
//       for (ordinal_type i=0;i<_m;++i) {
//         size_type jbegin = _ap[i], jend = _ap[i+1];
//         for (ordinal_type j=jbegin;j<jend;++j) {
//           SpT val = _ax[j];
//           os << std::setw(w) <<      i << "  "
//              << std::setw(w) << _aj[j] << "  "
//              << std::setw(w) <<    val << std::endl;
//         }
//       }
//     }
//     os.copyfmt(prec);

//     return os;
//   }

// template<class Scalar>
// size_type
// OrientationTools<Scalar>::DenseMatrix::countNumNonZeros(const Scalar epsilon) const {
//   size_type nnz = 0;
//   for (ordinal_type j=0;j<NumCols();++j) {
//     for (ordinal_type i=0;i<NumRows();++i) {
//       const Scalar val = Value(i,j);
//       nnz += ((val*val) > epsilon);
//     }
//   }
//   return nnz;
// }

// template<class Scalar>
// std::ostream&
// OrientationTools<Scalar>::DenseMatrix::showMe(std::ostream &os) const {
//   std::ofstream prec;
//   prec.copyfmt(os);

//   os.precision(3);

//   os << " -- OrientationTools::DenseMatrix -- " << std::endl
//      << "    # of Rows              = " << _m << std::endl
//      << "    # of Cols              = " << _n << std::endl
//      << "    Col Stride             = " << _cs << std::endl
//      << "    Row Stride             = " << _rs << std::endl
//      << std::endl
//      << "    ValueArray dimensions  = " << _a.extent(0) << std::endl
//      << std::endl;

//   const ordinal_type w = 10;
//   if (_a.size()) {
//     for (ordinal_type i=0;i<_m;++i) {
//       for (ordinal_type j=0;j<_n;++j) {
//         const Scalar val = this->Value(i,j);
//         os << std::setw(w) << val << "  ";
//       }
//       os << std::endl;
//     }
//   }
//   os.copyfmt(prec);

//   return os;
// }





// template<typename SpT>
// void
// OrientationTools<SpT>::
// initHexahedron(Kokkos::View<CoeffMatrixType***,SpT> MatrixData,
//                const EFunctionSpace space,
//                const ordinal_type order) {
    
//   switch (space) {
//   case FUNCTION_SPACE_HGRAD: {
//     Basis_HGRAD_LINE_Cn_FEM<SpT> lineBasis(order);
//     Basis_HGRAD_QUAD_Cn_FEM<SpT> faceBasis(order);
//     Basis_HGRAD_HEXA_Cn_FEM<SpT> cellBasis(order);
      
//     {
//       const ordinal_type numEdge = 12, numOrt = 2;
//       for (auto edgeId=0;edgeId<numEdge;++edgeId)
//         for (auto edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
//           const auto C = Impl::OrientationTools::getEdgeCoeffMatrix_HGRAD(lineBasis, cellBasis, edgeId, edgeOrt);
//           MatrixData(0, edgeId, edgeOrt) = Kokkos::create_mirror_view(typename SpT::memory_space(), C);
//           Kokkos::deep_copy(MatrixData(0, edgdId, edgeOrt), C);
//         }
//     }
//     {
//       const ordinal_type numFace = 6, numOrt = 8;
//       for (auto faceId=0;faceId<numFace;++faceId)
//         for (auto faceOrt=0;faceOrt<numOrt;++faceOrt) {
//           const auto C = Impl::OrientationTools::getQuadrilateralCoeffMatrix_HGRAD(faceBasis, cellBasis, faceId, faceOrt);
//           MatrixData(1, faceId, faceOrt) = Kokkos::create_mirror_view(typename SpT::memory_space(), C);
//           Kokkos::deep_copy(MatrixData(1, faceId, faceOrt), C);
//         }
//     }
//     break;
//   }            
//   default: {
//     INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
//                                   ">>> ERROR (Intrepid::OrientationTools::initQuadrilateral): Invalid function space.");
//     break;
//   }
//   }
// }

// template<typename SpT>
// void
// OrientationTools<SpT>::
// initTetrahedron(Kokkos::View<CoeffMatrixType***,SpT> MatrixData,
//                 const EFunctionSpace space,
//                 const ordinal_type order) {
    
//   switch (space) {
//   case FUNCTION_SPACE_HGRAD: {
//     Basis_HGRAD_LINE_Cn_FEM<SpT> lineBasis(order);
//     Basis_HGRAD_TRI_Cn_FEM <SpT> faceBasis(order);
//     Basis_HGRAD_TET_Cn_FEM <SpT> cellBasis(order);
      
//     {
//       const ordinal_type numEdge = 6, numOrt = 2;
//       for (auto edgeId=0;edgeId<numEdge;++edgeId)
//         for (auto edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
//           const auto C = Impl::OrientationTools::getEdgeCoeffMatrix_HGRAD(lineBasis, cellBasis, edgeId, edgeOrt);
//           MatrixData(0, edgeId, edgeOrt) = Kokkos::create_mirror_view(typename SpT::memory_space(), C);
//           Kokkos::deep_copy(MatrixData(0, edgdId, edgeOrt), C);
//         }
//     }
//     {
//       const ordinal_type numFace = 4, numOrt = 6;
//       for (auto faceId=0;faceId<numFace;++faceId)
//         for (auto faceOrt=0;faceOrt<numOrt;++faceOrt) {
//           const auto C = Impl::OrientationTools::getTriangleCoeffMatrix_HGRAD(faceBasis, cellBasis, faceId, faceOrt);
//           MatrixData(1, faceId, faceOrt) = Kokkos::create_mirror_view(typename SpT::memory_space(), C);
//           Kokkos::deep_copy(MatrixData(1, faceId, faceOrt), C);
//         }
//     }
//     break;
//   }            
//   default: {
//     INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
//                                   ">>> ERROR (Intrepid::OrientationTools::initQuadrilateral): Invalid function space.");
//     break;
//   }
//   }
// }
