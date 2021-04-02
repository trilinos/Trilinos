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


/** \file   Intrepid2_CellToolsDefNodeInfo.hpp
    \brief  Definition file for node data and subcell functions of the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEF_NODE_INFO_HPP__
#define __INTREPID2_CELLTOOLS_DEF_NODE_INFO_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {
  
  //============================================================================================//          
  //                                                                                            //          
  //                      Reference nodes                                                       //          
  //                                                                                            //          
  //============================================================================================//   

  template<typename DeviceType>
  template<typename cellCenterValueType, class ...cellCenterProperties>
  void 
  CellTools<DeviceType>::
  getReferenceCellCenter( Kokkos::DynRankView<cellCenterValueType,cellCenterProperties...> cellCenter,
                          const shards::CellTopology cell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !hasReferenceCell(cell), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceCellCenter): the specified cell topology does not have a reference cell." );
    
    INTREPID2_TEST_FOR_EXCEPTION( cellCenter.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceCellCenter): cellCenter must have rank 1." );

    INTREPID2_TEST_FOR_EXCEPTION( cellCenter.extent(0) < cell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceCellCenter): cellCenter must have dimension at least as large as cell.getDimension()." );
#endif

    constexpr bool is_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(cellCenter)::memory_space>::accessible;
    static_assert(is_accessible, "CellTools<DeviceType>::getReferenceCellCenter(..): output view's memory space is not compatible with DeviceType");

    const ordinal_type dim = cell.getDimension();

    const auto refCellCenter = RefCellCenter<DeviceType>::get(cell.getKey());

    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,dim),
    KOKKOS_LAMBDA (const int &i) {cellCenter(i) = refCellCenter(i);}
    );
  }


  template<typename DeviceType>
  template<typename cellVertexValueType, class ...cellVertexProperties>
  void
  CellTools<DeviceType>::
  getReferenceVertex(       Kokkos::DynRankView<cellVertexValueType,cellVertexProperties...> cellVertex,
                      const shards::CellTopology cell,
                      const ordinal_type         vertexOrd ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !hasReferenceCell(cell), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceVertex): the specified cell topology does not have a reference cell." );
    
    INTREPID2_TEST_FOR_EXCEPTION( (vertexOrd < 0) || vertexOrd > static_cast<ordinal_type>(cell.getVertexCount()), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceVertex): invalid node ordinal for the specified cell topology." );

    INTREPID2_TEST_FOR_EXCEPTION( cellVertex.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): cellNodes must have rank 1." );

    INTREPID2_TEST_FOR_EXCEPTION( cellVertex.extent(0) < cell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): cellNodes must have dimension at least as large as cell.getDimension()." );
#endif

    constexpr bool is_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(cellVertex)::memory_space>::accessible;
    static_assert(is_accessible, "CellTools<DeviceType>::getReferenceVertex(..): output view's memory space is not compatible with DeviceType");

    const auto refNodes = RefCellNodes<DeviceType>::get(cell.getKey());

    ordinal_type dim = cell.getDimension();
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,dim),
    KOKKOS_LAMBDA (const int &i) {cellVertex(i) = refNodes(vertexOrd,i);}
    );
  }
    
  
  template<typename DeviceType>
  template<typename subcellVertexValueType, class ...subcellVertexProperties>
  void
  CellTools<DeviceType>::
  getReferenceSubcellVertices(       Kokkos::DynRankView<subcellVertexValueType,subcellVertexProperties...> subcellVertices,
                               const ordinal_type         subcellDim,
                               const ordinal_type         subcellOrd,
                               const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !hasReferenceCell(parentCell), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): the specified cell topology does not have a reference cell." );

    INTREPID2_TEST_FOR_EXCEPTION( subcellDim > static_cast<ordinal_type>(parentCell.getDimension()), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): subcell dimension cannot exceed cell dimension." );
    
    INTREPID2_TEST_FOR_EXCEPTION( subcellOrd >= static_cast<ordinal_type>(parentCell.getSubcellCount(subcellDim)), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): subcell ordinal cannot exceed subcell count." );
    
    // Verify subcellVertices rank and dimensions
    INTREPID2_TEST_FOR_EXCEPTION( subcellVertices.rank() != 2, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): subcellVertieces must have rank 2." );
    
    // need to match to node count as it invokes getReferenceSubcellNodes
    INTREPID2_TEST_FOR_EXCEPTION( subcellVertices.extent(0) != parentCell.getNodeCount(subcellDim, subcellOrd), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): subcellVertieces dimension(0) must match to parent cell vertex count." );
    
    INTREPID2_TEST_FOR_EXCEPTION( subcellVertices.extent(1) != parentCell.getDimension(), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): subcellVertieces dimension(1) must match to parent cell dimension." );
#endif 
    getReferenceSubcellNodes(subcellVertices, 
                             subcellDim, 
                             subcellOrd, 
                             parentCell);
  }  
  

  template<typename DeviceType>
  template<typename cellNodeValueType, class ...cellNodeProperties>
  void
  CellTools<DeviceType>::
  getReferenceNode(       Kokkos::DynRankView<cellNodeValueType,cellNodeProperties...> cellNode,
                    const shards::CellTopology  cell,
                    const ordinal_type          nodeOrd ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( nodeOrd >= static_cast<ordinal_type>(cell.getNodeCount()), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): invalid node ordinal for the specified cell topology." );

    INTREPID2_TEST_FOR_EXCEPTION( cellNode.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): cellNode must have rank 1." );

    INTREPID2_TEST_FOR_EXCEPTION( cellNode.extent(0) < cell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): cellNode must have dimension at least as large as cell.getDimension()." );
#endif

    constexpr bool is_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(cellNode)::memory_space>::accessible;
    static_assert(is_accessible, "CellTools<DeviceType>::getReferenceNode(..): output view's memory space is not compatible with DeviceType");

    const auto refNodes = RefCellNodes<DeviceType>::get(cell.getKey());

    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,cell.getDimension()),
    KOKKOS_LAMBDA (const int &i) {cellNode(i) = refNodes(nodeOrd,i);}
    );
  }

  template<typename DeviceType>
  template<typename subcellNodeValueType, class ...subcellNodeProperties>
  void
  CellTools<DeviceType>::
  getReferenceSubcellNodes(       Kokkos::DynRankView<subcellNodeValueType,subcellNodeProperties...> subcellNodes,
                            const ordinal_type         subcellDim,
                            const ordinal_type         subcellOrd,
                            const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !hasReferenceCell(parentCell), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): the specified cell topology does not have a reference cell.");

    INTREPID2_TEST_FOR_EXCEPTION( subcellDim > static_cast<ordinal_type>(parentCell.getDimension()), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): subcell dimension out of range.");
    
    INTREPID2_TEST_FOR_EXCEPTION( subcellOrd >= static_cast<ordinal_type>(parentCell.getSubcellCount(subcellDim)), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): subcell ordinal out of range.");
    
    // Verify subcellNodes rank and dimensions
    INTREPID2_TEST_FOR_EXCEPTION( subcellNodes.rank() != 2, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): subcellNodes must have rank 2.");
      
    INTREPID2_TEST_FOR_EXCEPTION( subcellNodes.extent(0) != parentCell.getNodeCount(subcellDim, subcellOrd), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): subcellNodes dimension (0) must match to node count in cell.");
    
    INTREPID2_TEST_FOR_EXCEPTION( subcellNodes.extent(1) != parentCell.getDimension(), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): subcellNodes dimension (1) must match to cell dimension.");
#endif 
    
    // Find how many cellWorkset does the specified subcell have.
    const auto subcNodeCount = parentCell.getNodeCount(subcellDim, subcellOrd);
    
    // Loop over subcell cellWorkset
    for (size_type subcNodeOrd=0;subcNodeOrd<subcNodeCount;++subcNodeOrd) {      
      // Get the node number relative to the parent reference cell
      const auto cellNodeOrd = parentCell.getNodeMap(subcellDim, subcellOrd, subcNodeOrd);

      auto dst = Kokkos::subdynrankview(subcellNodes, subcNodeOrd, Kokkos::ALL());
      getReferenceNode(dst, parentCell, cellNodeOrd);
    }
  }  

  template<typename DeviceType>
  template<typename refEdgeTangentValueType, class ...refEdgeTangentProperties>
  void
  CellTools<DeviceType>::
  getReferenceEdgeTangent(       Kokkos::DynRankView<refEdgeTangentValueType,refEdgeTangentProperties...> refEdgeTangent,
                           const ordinal_type         edgeOrd,
                           const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 2 &&
                                  parentCell.getDimension() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceEdgeTangent): two or three-dimensional parent cell required");
  
    INTREPID2_TEST_FOR_EXCEPTION( refEdgeTangent.rank() != 1, std::invalid_argument,  
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceEdgeTangent): rank = 1 required for output arrays");
    
    INTREPID2_TEST_FOR_EXCEPTION( refEdgeTangent.extent(0) != parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceEdgeTangent): output array size is required to match space dimension");

    INTREPID2_TEST_FOR_EXCEPTION( edgeOrd <  0 ||
                                  edgeOrd >= static_cast<ordinal_type>(parentCell.getSubcellCount(1)), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceEdgeTangent): edge ordinal out of bounds");

#endif
    constexpr bool is_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(refEdgeTangent)::memory_space>::accessible;
    static_assert(is_accessible, "CellTools<DeviceType>::getReferenceEdgeTangent(..): output view's memory space is not compatible with DeviceType");

    const auto edgeMap = RefSubcellParametrization<DeviceType>::get(1, parentCell.getKey());
  
    // All ref. edge maps have affine coordinate functions: f_dim(u) = C_0(dim) + C_1(dim)*u, 
    //                                     => edge Tangent: -> C_1(*)
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,parentCell.getDimension()),
    KOKKOS_LAMBDA (const int &i) {refEdgeTangent(i) = edgeMap(edgeOrd, i, 1);}
    );
  }


  template<typename DeviceType>
  template<typename refFaceTanValueType, class ...refFaceTanProperties>
  void
  CellTools<DeviceType>::
  getReferenceFaceTangents(       Kokkos::DynRankView<refFaceTanValueType,refFaceTanProperties...> refFaceTanU,
                                  Kokkos::DynRankView<refFaceTanValueType,refFaceTanProperties...> refFaceTanV,
                            const ordinal_type         faceOrd,
                            const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): three-dimensional parent cell required");  
  
    INTREPID2_TEST_FOR_EXCEPTION( faceOrd < 0 || faceOrd >= static_cast<ordinal_type>(parentCell.getSubcellCount(2)), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): face ordinal out of bounds");  
    
    INTREPID2_TEST_FOR_EXCEPTION( refFaceTanU.rank() != 1 || refFaceTanV.rank() != 1, std::invalid_argument,  
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): rank = 1 required for output arrays"); 
  
    INTREPID2_TEST_FOR_EXCEPTION( refFaceTanU.extent(0) != parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): dim0 (spatial dim) must match that of parent cell");  

    INTREPID2_TEST_FOR_EXCEPTION( refFaceTanV.extent(0) != parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): dim0 (spatial dim) must match that of parent cell");  
#endif
    constexpr bool is_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(refFaceTanU)::memory_space>::accessible;
    static_assert(is_accessible, "CellTools<DeviceType>::getReferenceFaceTangents(..): output views' memory spaces are not compatible with DeviceType");


    const auto faceMap = RefSubcellParametrization<DeviceType>::get(2, parentCell.getKey());
  
    /*  All ref. face maps have affine coordinate functions:  f_dim(u,v) = C_0(dim) + C_1(dim)*u + C_2(dim)*v
     *                           `   => Tangent vectors are:  refFaceTanU -> C_1(*);    refFaceTanV -> C_2(*)
     */

    // set refFaceTanU -> C_1(*)
    // set refFaceTanV -> C_2(*)
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,parentCell.getDimension()),
        KOKKOS_LAMBDA (const int &i) {
      refFaceTanU(i) = faceMap(faceOrd, i, 1);
      refFaceTanV(i) = faceMap(faceOrd, i, 2);
      });
  }

  template<typename DeviceType>
  template<typename refSideNormalValueType, class ...refSideNormalProperties>
  void
  CellTools<DeviceType>::
  getReferenceSideNormal(       Kokkos::DynRankView<refSideNormalValueType,refSideNormalProperties...> refSideNormal,
                          const ordinal_type         sideOrd,
                          const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 2 &&
                                  parentCell.getDimension() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSideNormal): two or three-dimensional parent cell required");
  
    // Check side ordinal: by definition side is subcell whose dimension = parentCell.getDimension()-1
    INTREPID2_TEST_FOR_EXCEPTION( sideOrd < 0 || sideOrd >= static_cast<ordinal_type>(parentCell.getSubcellCount(parentCell.getDimension() - 1)), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSideNormal): side ordinal out of bounds");    
#endif 
    constexpr bool is_accessible = Kokkos::Impl::MemorySpaceAccess<
        MemSpaceType,
        typename decltype(refSideNormal)::memory_space>::accessible;
    static_assert(is_accessible, "CellTools<DeviceType>::getReferenceSideNormal(..): output view's memory space is not compatible with DeviceType");

    const auto dim = parentCell.getDimension();
    if (dim == 2) {
      // 2D parent cells: side = 1D subcell (edge), call the edge tangent method and rotate tangents
      getReferenceEdgeTangent(refSideNormal, sideOrd, parentCell);
    
      // rotate t(t1, t2) to get n(t2, -t1) so that (n,t) is positively oriented: det(n1,n2/t1,t2)>0
      Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,1),
              KOKKOS_LAMBDA (const int &) {
        refSideNormalValueType tmp = refSideNormal(0);
        refSideNormal(0) = refSideNormal(1);
        refSideNormal(1) = -tmp;
      });
    } else {
      // 3D parent cell: side = 2D subcell (face), call the face normal method.
      getReferenceFaceNormal(refSideNormal, sideOrd, parentCell);
    }
  }


  template<typename DeviceType>
  template<typename refFaceNormalValueType, class ...refFaceNormalProperties>
  void 
  CellTools<DeviceType>::
  getReferenceFaceNormal(       Kokkos::DynRankView<refFaceNormalValueType,refFaceNormalProperties...> refFaceNormal,
                          const ordinal_type         faceOrd,
                          const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceNormal): three-dimensional parent cell required");  
    
    INTREPID2_TEST_FOR_EXCEPTION( faceOrd < 0 || faceOrd >= static_cast<ordinal_type>(parentCell.getSubcellCount(2)), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceNormal): face ordinal out of bounds");  
    
    INTREPID2_TEST_FOR_EXCEPTION( refFaceNormal.rank() != 1, std::invalid_argument,  
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceNormal): rank = 1 required for output array"); 
  
    INTREPID2_TEST_FOR_EXCEPTION( refFaceNormal.extent(0) != parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceNormal): dim0 (spatial dim) must match that of parent cell");  
#endif
    constexpr bool is_accessible = Kokkos::Impl::MemorySpaceAccess<
        MemSpaceType,
        typename decltype(refFaceNormal)::memory_space>::accessible;
    static_assert(is_accessible, "CellTools<DeviceType>::getReferenceFaceNormal(..): output view's memory space is not compatible with DeviceType");

    // Reference face normal = vector product of reference face tangents. Allocate temp FC storage:
    const auto dim = parentCell.getDimension();
    auto vcprop = Kokkos::common_view_alloc_prop(refFaceNormal);
    using common_value_type = typename decltype(vcprop)::value_type;
    Kokkos::DynRankView< common_value_type, DeviceType > refFaceTanU ( Kokkos::view_alloc("CellTools::getReferenceFaceNormal::refFaceTanU", vcprop), dim );
    Kokkos::DynRankView< common_value_type, DeviceType > refFaceTanV ( Kokkos::view_alloc("CellTools::getReferenceFaceNormal::refFaceTanV", vcprop), dim );
    getReferenceFaceTangents(refFaceTanU, refFaceTanV, faceOrd, parentCell);
  
    RealSpaceTools<DeviceType>::vecprod(refFaceNormal, refFaceTanU, refFaceTanV);
  }

  template<typename DeviceType>
  template<typename edgeTangentValueType,     class ...edgeTangentProperties,
           typename worksetJacobianValueType, class ...worksetJacobianProperties>
  void
  CellTools<DeviceType>::
  getPhysicalEdgeTangents(       Kokkos::DynRankView<edgeTangentValueType,edgeTangentProperties...>         edgeTangents,
                           const Kokkos::DynRankView<worksetJacobianValueType,worksetJacobianProperties...> worksetJacobians,
                           const ordinal_type         worksetEdgeOrd,
                           const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 3 &&
                                  parentCell.getDimension() != 2, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): 2D or 3D parent cell required." );  
  
    // (1) edgeTangents is rank-3 (C,P,D) and D=2, or 3 is required
    INTREPID2_TEST_FOR_EXCEPTION( edgeTangents.rank() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): edgeTangents requires rank 3." );  
    INTREPID2_TEST_FOR_EXCEPTION( edgeTangents.extent(2) != 2 && 
                                  edgeTangents.extent(2) != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): edgeTangents dimension(2) must be 2 or 3." );
 
    // (2) worksetJacobians in rank-4 (C,P,D,D) and D=2, or 3 is required
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.rank() != 4, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): worksetJacobians requires rank 4." );  
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.extent(2) != 2 && 
                                  worksetJacobians.extent(2) != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): worksetJacobians dimension(2) must be 2 or 3." );
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.extent(2) != worksetJacobians.extent(3), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): worksetJacobians dimension(2) and (3) must match each other." );

    // (4) cross-check array dimensions: edgeTangents (C,P,D) vs. worksetJacobians (C,P,D,D)
    for (auto i=0;i<3;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( edgeTangents.extent(i) != worksetJacobians.extent(i), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): edgeTangents dimension (i) does not match to worksetJacobians dimension(i)." );
    }
#endif
    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(edgeTangents)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(worksetJacobians)::memory_space>::accessible;
    static_assert(are_accessible, "CellTools<DeviceType>::getPhysicalEdgeTangents(..): input/output views' memory spaces are not compatible with DeviceType");

  
    // Storage for constant reference edge tangent: rank-1 (D) arrays
    const auto dim = parentCell.getDimension();
    auto vcprop = Kokkos::common_view_alloc_prop(edgeTangents);
    using common_value_type = typename decltype(vcprop)::value_type;
    Kokkos::DynRankView< common_value_type, DeviceType > refEdgeTan ( Kokkos::view_alloc("CellTools::getPhysicalEdgeTangents::refEdgeTan", vcprop), dim);
    getReferenceEdgeTangent(refEdgeTan, worksetEdgeOrd, parentCell);
    
    RealSpaceTools<DeviceType>::matvec(edgeTangents, worksetJacobians, refEdgeTan);
  }


  template<typename DeviceType>
  template<typename faceTanValueType,        class ...faceTanProperties,
           typename worksetJacobianValueType, class ...worksetJacobianProperties>
  void
  CellTools<DeviceType>::
  getPhysicalFaceTangents(       Kokkos::DynRankView<faceTanValueType,faceTanProperties...> faceTanU,
                                 Kokkos::DynRankView<faceTanValueType,faceTanProperties...> faceTanV,
                           const Kokkos::DynRankView<worksetJacobianValueType,worksetJacobianProperties...> worksetJacobians,
                           const ordinal_type         worksetFaceOrd,
                           const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): three-dimensional parent cell required");  
  
    // (1) faceTanU and faceTanV are rank-3 (C,P,D) and D=3 is required
    INTREPID2_TEST_FOR_EXCEPTION( faceTanU.rank() != 3 || 
                                  faceTanV.rank() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): faceTan U,V must have rank 3." );  

    INTREPID2_TEST_FOR_EXCEPTION( faceTanU.extent(2) != 3 ||
                                  faceTanV.extent(2) != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): faceTan U,V dimension (2) must be 3." );  
    
    for (auto i=0;i<3;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( faceTanU.extent(i) != faceTanV.extent(i), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): faceTan U,V dimension (i) must match each other." );  
    }

    // (3) worksetJacobians in rank-4 (C,P,D,D) and D=3 is required
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.rank() != 4, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): worksetJacobians must have rank 4." );  

    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.extent(2) != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): worksetJacobians dimension(2) must be 3." );  

    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.extent(2) != worksetJacobians.extent(3), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): worksetJacobians dimension(2) and dimension(3) must match." );  

    // (4) cross-check array dimensions: faceTanU (C,P,D) vs. worksetJacobians (C,P,D,D)
    for (auto i=0;i<3;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( faceTanU.extent(i) != worksetJacobians.extent(i), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): worksetJacobians dimension(i) and faceTan dimension (i) must match." );  
    }      
#endif
    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(faceTanU)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(worksetJacobians)::memory_space>::accessible;
    static_assert(are_accessible, "CellTools<DeviceType>::getPhysicalFaceTangents(..): input/output views' memory spaces are not compatible with DeviceType");

    // Temp storage for the pair of constant ref. face tangents: rank-1 (D) arrays
    const auto dim = parentCell.getDimension();

    auto vcprop = Kokkos::common_view_alloc_prop(faceTanU);
    using common_value_type = typename decltype(vcprop)::value_type;
    Kokkos::DynRankView< common_value_type, DeviceType > refFaceTanU ( Kokkos::view_alloc("CellTools::getPhysicalFaceTangents::refFaceTanU", vcprop), dim);
    Kokkos::DynRankView< common_value_type, DeviceType > refFaceTanV ( Kokkos::view_alloc("CellTools::getPhysicalFaceTangents::refFaceTanV", vcprop), dim);

    getReferenceFaceTangents(refFaceTanU, refFaceTanV, worksetFaceOrd, parentCell);

    RealSpaceTools<DeviceType>::matvec(faceTanU, worksetJacobians, refFaceTanU);
    RealSpaceTools<DeviceType>::matvec(faceTanV, worksetJacobians, refFaceTanV);
  }

  namespace FunctorCellTools {
  template<typename normalsViewType,
           typename tangentsViewType>
  struct F_edgeNormalsFromTangents {
    normalsViewType edgeNormals_;
    const tangentsViewType edgeTangents_;

    KOKKOS_INLINE_FUNCTION
    F_edgeNormalsFromTangents( normalsViewType  edgeNormals,
                           const tangentsViewType refEdgeTangents)
      : edgeNormals_(edgeNormals), edgeTangents_(refEdgeTangents){};

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type iter) const {
      size_type cell, pt;
      unrollIndex( cell, pt, edgeTangents_.extent(0),
          edgeTangents_.extent(1), iter );
      edgeNormals_(cell,pt,0) =  edgeTangents_(cell,pt,1);
      edgeNormals_(cell,pt,1) = -edgeTangents_(cell,pt,0);
    }
  };
}

  template<typename DeviceType>
  template<typename sideNormalValueType,      class ...sideNormalProperties,
           typename worksetJacobianValueType, class ...worksetJacobianProperties>
  void 
  CellTools<DeviceType>::
  getPhysicalSideNormals(       Kokkos::DynRankView<sideNormalValueType,sideNormalProperties...> sideNormals,
                          const Kokkos::DynRankView<worksetJacobianValueType,worksetJacobianProperties...> worksetJacobians,
                          const ordinal_type         worksetSideOrd,
                          const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 2 && 
                                  parentCell.getDimension() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalSideNormals): two or three-dimensional parent cell required");
  
    // Check side ordinal: by definition side is subcell whose dimension = parentCell.getDimension()-1
    INTREPID2_TEST_FOR_EXCEPTION( worksetSideOrd <  0 ||
                                  worksetSideOrd >= static_cast<ordinal_type>(parentCell.getSubcellCount(parentCell.getDimension() - 1)), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalSideNormals): side ordinal out of bounds");  
#endif  
    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(sideNormals)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(worksetJacobians)::memory_space>::accessible;
    static_assert(are_accessible, "CellTools<DeviceType>::getPhysicalSideNormals(..): input/output views' memory spaces are not compatible with DeviceType");

    const auto dim = parentCell.getDimension();
  
    if (dim == 2) {
      // compute edge tangents and rotate it
      auto vcprop = Kokkos::common_view_alloc_prop(sideNormals);
      using common_value_type = typename decltype(vcprop)::value_type;
      Kokkos::DynRankView< common_value_type, DeviceType > edgeTangents ( Kokkos::view_alloc("CellTools::getPhysicalSideNormals::edgeTan", vcprop),
                                                              sideNormals.extent(0),
                                                              sideNormals.extent(1),
                                                              sideNormals.extent(2));
      getPhysicalEdgeTangents(edgeTangents, worksetJacobians, worksetSideOrd, parentCell);

      //Note: this function has several template parameters and the compiler gets confused if using a lambda function
      using FunctorType = FunctorCellTools::F_edgeNormalsFromTangents<decltype(sideNormals), decltype(edgeTangents)>;
      const auto loopSize = edgeTangents.extent(0)*edgeTangents.extent(1);
          Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
          Kokkos::parallel_for( policy, FunctorType(sideNormals, edgeTangents) );

    } else {
      getPhysicalFaceNormals(sideNormals, worksetJacobians, worksetSideOrd, parentCell);
    }
  }
  

  template<typename DeviceType>
  template<typename faceNormalValueType,      class ...faceNormalProperties,
           typename worksetJacobianValueType, class ...worksetJacobianProperties>
  void
  CellTools<DeviceType>::
  getPhysicalFaceNormals(       Kokkos::DynRankView<faceNormalValueType,faceNormalProperties...> faceNormals,
                          const Kokkos::DynRankView<worksetJacobianValueType,worksetJacobianProperties...> worksetJacobians,
                          const ordinal_type         worksetFaceOrd,
                          const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): three-dimensional parent cell required." );  
    
    // (1) faceNormals is rank-3 (C,P,D) and D=3 is required
    INTREPID2_TEST_FOR_EXCEPTION( faceNormals.rank() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): faceNormals must have a rank 3." );
    INTREPID2_TEST_FOR_EXCEPTION( faceNormals.extent(2) != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): faceNormals dimension (2) must be 3." );
    
    // (3) worksetJacobians in rank-4 (C,P,D,D) and D=3 is required
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.rank() != 4, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): worksetJacobians must have a rank 4." );
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.extent(2) != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): worksetJacobians dimension (2) must be 3." );
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.extent(2) != worksetJacobians.extent(3), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): worksetJacobians dimension (2) must match to dimension (3)." );
  
    // (4) cross-check array dimensions: faceNormals (C,P,D) vs. worksetJacobians (C,P,D,D)
    for (auto i=0;i<3;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( faceNormals.extent(i) != worksetJacobians.extent(i), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): faceNormals dimension (i) must match to worksetJacobians dimension (i)." );
    }        
#endif
    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(faceNormals)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(worksetJacobians)::memory_space>::accessible;
    static_assert(are_accessible, "CellTools<DeviceType>::getPhysicalFaceNormals(..): input/output views' memory spaces are not compatible with DeviceType");

  
    // this should be provided from users
    // Storage for physical face tangents: rank-3 (C,P,D) arrays
    const auto worksetSize = worksetJacobians.extent(0);
    const auto facePtCount = worksetJacobians.extent(1);
    const auto dim = parentCell.getDimension();

    auto vcprop = Kokkos::common_view_alloc_prop(faceNormals);
    using common_value_type = typename decltype(vcprop)::value_type;
    Kokkos::DynRankView< common_value_type, DeviceType > faceTanU ( Kokkos::view_alloc("CellTools::getPhysicalFaceNormals::faceTanU", vcprop), worksetSize, facePtCount, dim);
    Kokkos::DynRankView< common_value_type, DeviceType > faceTanV ( Kokkos::view_alloc("CellTools::getPhysicalFaceNormals::faceTanV", vcprop), worksetSize, facePtCount, dim);

    getPhysicalFaceTangents(faceTanU, faceTanV, 
                            worksetJacobians, 
                            worksetFaceOrd, 
                            parentCell);
  
    RealSpaceTools<DeviceType>::vecprod(faceNormals, faceTanU, faceTanV);
  }
}

#endif
