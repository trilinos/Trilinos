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


/** \file   Intrepid_CellToolsDef.hpp
    \brief  Definition file for the Intrepid2::CellTools class.
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

  template<typename SpT>
  void 
  CellTools<SpT>::
  setReferenceNodeData() {

    {
      // create memory on devices
      refNodeData_.line            = referenceNodeDataViewType("CellTools::ReferenceNodeData::line",       2, 3);
      refNodeData_.line_3          = referenceNodeDataViewType("CellTools::ReferenceNodeData::line_3",     3, 3);
      
      refNodeData_.triangle        = referenceNodeDataViewType("CellTools::ReferenceNodeData::triangle",   3, 3);
      refNodeData_.triangle_4      = referenceNodeDataViewType("CellTools::ReferenceNodeData::triangle_4", 4, 3);    
      refNodeData_.triangle_6      = referenceNodeDataViewType("CellTools::ReferenceNodeData::triangle_6", 6, 3);    
      
      refNodeData_.quadrilateral   = referenceNodeDataViewType("CellTools::ReferenceNodeData::quad",       4, 3);    
      refNodeData_.quadrilateral_8 = referenceNodeDataViewType("CellTools::ReferenceNodeData::quad_8",     8, 3);    
      refNodeData_.quadrilateral_9 = referenceNodeDataViewType("CellTools::ReferenceNodeData::quad_9",     9, 3);
      
      refNodeData_.tetrahedron     = referenceNodeDataViewType("CellTools::ReferenceNodeData::tet",        4, 3);    
      refNodeData_.tetrahedron_8   = referenceNodeDataViewType("CellTools::ReferenceNodeData::tet_8",      8, 3);    
      refNodeData_.tetrahedron_10  = referenceNodeDataViewType("CellTools::ReferenceNodeData::tet_10",    10, 3);
      refNodeData_.tetrahedron_11  = referenceNodeDataViewType("CellTools::ReferenceNodeData::tet_11",    11, 3);    
      
      refNodeData_.hexahedron      = referenceNodeDataViewType("CellTools::ReferenceNodeData::hex",        8, 3);    
      refNodeData_.hexahedron_20   = referenceNodeDataViewType("CellTools::ReferenceNodeData::hex_20",    20, 3);
      refNodeData_.hexahedron_27   = referenceNodeDataViewType("CellTools::ReferenceNodeData::hex_27",    27, 3);    
      
      refNodeData_.pyramid         = referenceNodeDataViewType("CellTools::ReferenceNodeData::pyr",        5, 3);    
      refNodeData_.pyramid_13      = referenceNodeDataViewType("CellTools::ReferenceNodeData::pyr_13",    13, 3);
      refNodeData_.pyramid_14      = referenceNodeDataViewType("CellTools::ReferenceNodeData::pyr_14",    14, 3);    
      
      refNodeData_.wedge           = referenceNodeDataViewType("CellTools::ReferenceNodeData::wedge",      6, 3);
      refNodeData_.wedge_15        = referenceNodeDataViewType("CellTools::ReferenceNodeData::wedge_15",  15, 3);    
      refNodeData_.wedge_18        = referenceNodeDataViewType("CellTools::ReferenceNodeData::wedge_18",  18, 3);    
    }

    {
      // copy static data to devices
      Kokkos::deep_copy(refNodeData_.line,            referenceNodeDataViewHostType(&refNodeDataStatic_.line[0][0],             2, 3));
      Kokkos::deep_copy(refNodeData_.line_3,          referenceNodeDataViewHostType(&refNodeDataStatic_.line_3[0][0],           3, 3));
      
      Kokkos::deep_copy(refNodeData_.triangle,        referenceNodeDataViewHostType(&refNodeDataStatic_.triangle[0][0],         3, 3));
      Kokkos::deep_copy(refNodeData_.triangle_4,      referenceNodeDataViewHostType(&refNodeDataStatic_.triangle_4[0][0],       4, 3));    
      Kokkos::deep_copy(refNodeData_.triangle_6,      referenceNodeDataViewHostType(&refNodeDataStatic_.triangle_6[0][0],       6, 3));    
      
      Kokkos::deep_copy(refNodeData_.quadrilateral,   referenceNodeDataViewHostType(&refNodeDataStatic_.quadrilateral[0][0],    4, 3));    
      Kokkos::deep_copy(refNodeData_.quadrilateral_8, referenceNodeDataViewHostType(&refNodeDataStatic_.quadrilateral_8[0][0],  8, 3));    
      Kokkos::deep_copy(refNodeData_.quadrilateral_9, referenceNodeDataViewHostType(&refNodeDataStatic_.quadrilateral_9[0][0],  9, 3));
      
      Kokkos::deep_copy(refNodeData_.tetrahedron,     referenceNodeDataViewHostType(&refNodeDataStatic_.tetrahedron[0][0],      4, 3));    
      Kokkos::deep_copy(refNodeData_.tetrahedron_8,   referenceNodeDataViewHostType(&refNodeDataStatic_.tetrahedron_8[0][0],    8, 3));    
      Kokkos::deep_copy(refNodeData_.tetrahedron_10,  referenceNodeDataViewHostType(&refNodeDataStatic_.tetrahedron_10[0][0],  10, 3));
      Kokkos::deep_copy(refNodeData_.tetrahedron_11,  referenceNodeDataViewHostType(&refNodeDataStatic_.tetrahedron_11[0][0],  11, 3));    
      
      Kokkos::deep_copy(refNodeData_.hexahedron,      referenceNodeDataViewHostType(&refNodeDataStatic_.hexahedron[0][0],       8, 3));    
      Kokkos::deep_copy(refNodeData_.hexahedron_20,   referenceNodeDataViewHostType(&refNodeDataStatic_.hexahedron_20[0][0],   20, 3));
      Kokkos::deep_copy(refNodeData_.hexahedron_27,   referenceNodeDataViewHostType(&refNodeDataStatic_.hexahedron_27[0][0],   27, 3));    
      
      Kokkos::deep_copy(refNodeData_.pyramid,         referenceNodeDataViewHostType(&refNodeDataStatic_.pyramid[0][0],          5, 3));    
      Kokkos::deep_copy(refNodeData_.pyramid_13,      referenceNodeDataViewHostType(&refNodeDataStatic_.pyramid_13[0][0],      13, 3));
      Kokkos::deep_copy(refNodeData_.pyramid_14,      referenceNodeDataViewHostType(&refNodeDataStatic_.pyramid_14[0][0],      14, 3));    
      
      Kokkos::deep_copy(refNodeData_.wedge,           referenceNodeDataViewHostType(&refNodeDataStatic_.wedge[0][0],            6, 3));
      Kokkos::deep_copy(refNodeData_.wedge_15,        referenceNodeDataViewHostType(&refNodeDataStatic_.wedge_15[0][0],        15, 3));    
      Kokkos::deep_copy(refNodeData_.wedge_18,        referenceNodeDataViewHostType(&refNodeDataStatic_.wedge_18[0][0],        18, 3));    
    }

    isReferenceNodeDataSet_ = true;
  }

  template<typename SpT>
  template<typename cellVertexValueType, class ...cellVertexProperties>
  void
  CellTools<SpT>::
  getReferenceVertex( /**/  Kokkos::DynRankView<cellVertexValueType,cellVertexProperties...> cellVertices,
                      const shards::CellTopology cell,
                      const ordinal_type         vertexOrd ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !hasReferenceCell(cell), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceVertex): the specified cell topology does not have a reference cell." );
    
    INTREPID2_TEST_FOR_EXCEPTION( (vertexOrd < 0) || vertexOrd > cell.getVertexCount(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceVertex): invalid node ordinal for the specified cell topology." );

    INTREPID2_TEST_FOR_EXCEPTION( cellVertices.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): cellNodes must have rank 1." );

    INTREPID2_TEST_FOR_EXCEPTION( cellVertices.dimension(0) < cell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): cellNodes must have dimension bigger than Parameters::MaxDimension." );
#endif
    getReferenceNode(cellVertices, 
                     cell, 
                     vertexOrd);
  }
    
  
  template<typename SpT>
  template<typename subcellVertexValueType, class ...subcellVertexProperties>
  void
  CellTools<SpT>::
  getReferenceSubcellVertices( /**/  Kokkos::DynRankView<subcellVertexValueType,subcellVertexProperties...> subcellVertices,
                               const ordinal_type         subcellDim,
                               const ordinal_type         subcellOrd,
                               const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !hasReferenceCell(parentCell), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): the specified cell topology does not have a reference cell." );

    INTREPID2_TEST_FOR_EXCEPTION( subcellDim > parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): subcell dimension cannot exceed cell dimension." );
    
    INTREPID2_TEST_FOR_EXCEPTION( subcellOrd >= parentCell.getSubcellCount(subcellDim), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): subcell ordinal cannot exceed subcell count." );
    
    // Verify subcellVertices rank and dimensions
    INTREPID2_TEST_FOR_EXCEPTION( subcellVertices.rank() != 2, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): subcellVertieces must have rank 2." );
    
    // need to match to node count as it invokes getReferenceSubcellNodes
    INTREPID2_TEST_FOR_EXCEPTION( subcellVertices.dimension(0) != parentCell.getNodeCount(subcellDim, subcellOrd), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): subcellVertieces dimension(0) must match to parent cell vertex count." );
    
    INTREPID2_TEST_FOR_EXCEPTION( subcellVertices.dimension(1) != parentCell.getDimension(), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): subcellVertieces dimension(1) must match to parent cell dimension." );
#endif 
    getReferenceSubcellNodes(subcellVertices, 
                             subcellDim, 
                             subcellOrd, 
                             parentCell);
  }  
  

  template<typename SpT>
  template<typename cellNodeValueType, class ...cellNodeProperties>
  void
  CellTools<SpT>::
  getReferenceNode( /**/  Kokkos::DynRankView<cellNodeValueType,cellNodeProperties...> cellNodes,
                    const shards::CellTopology  cell,
                    const ordinal_type          nodeOrd ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !hasReferenceCell(cell), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): the specified cell topology does not have a reference cell." );
    
    INTREPID2_TEST_FOR_EXCEPTION( nodeOrd >= cell.getNodeCount(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): invalid node ordinal for the specified cell topology." );

    INTREPID2_TEST_FOR_EXCEPTION( cellNodes.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): cellNodes must have rank 1." );

    INTREPID2_TEST_FOR_EXCEPTION( cellNodes.dimension(0) < cell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): cellNodes must have dimension bigger than Parameters::MaxDimension." );
#endif

    if (!isReferenceNodeDataSet_) 
      setReferenceNodeData();

    //Kokkos::DynRankView<value_type,SpT,Kokkos::MemoryUnmanaged> ref;
    referenceNodeDataViewType ref;

    switch (cell.getKey() ) {
    case shards::Line<2>::key:     
    case shards::ShellLine<2>::key:
    case shards::Beam<2>::key:               ref = refNodeData_.line; break;
    case shards::Line<3>::key:     
    case shards::ShellLine<3>::key:
    case shards::Beam<3>::key:               ref = refNodeData_.line_3; break;
      
    case shards::Triangle<3>::key: 
    case shards::ShellTriangle<3>::key:      ref = refNodeData_.triangle; break;
    case shards::Triangle<4>::key:           ref = refNodeData_.triangle_4; break;
    case shards::Triangle<6>::key:
    case shards::ShellTriangle<6>::key:      ref = refNodeData_.triangle_6; break;
        
    case shards::Quadrilateral<4>::key:
    case shards::ShellQuadrilateral<4>::key: ref = refNodeData_.quadrilateral; break;
    case shards::Quadrilateral<8>::key:
    case shards::ShellQuadrilateral<8>::key: ref = refNodeData_.quadrilateral_8; break;
    case shards::Quadrilateral<9>::key:
    case shards::ShellQuadrilateral<9>::key: ref = refNodeData_.quadrilateral_9; break;

    case shards::Tetrahedron<4>::key:        ref = refNodeData_.tetrahedron; break;
    case shards::Tetrahedron<8>::key:        ref = refNodeData_.tetrahedron_8; break;
    case shards::Tetrahedron<10>::key:       ref = refNodeData_.tetrahedron_10; break;
    case shards::Tetrahedron<11>::key:       ref = refNodeData_.tetrahedron_11; break;

    case shards::Hexahedron<8>::key:         ref = refNodeData_.hexahedron; break;
    case shards::Hexahedron<20>::key:        ref = refNodeData_.hexahedron_20; break;
    case shards::Hexahedron<27>::key:        ref = refNodeData_.hexahedron_27; break;

    case shards::Pyramid<5>::key:            ref = refNodeData_.pyramid; break;
    case shards::Pyramid<13>::key:           ref = refNodeData_.pyramid_13; break;
    case shards::Pyramid<14>::key:           ref = refNodeData_.pyramid_14; break;

    case shards::Wedge<6>::key:              ref = refNodeData_.wedge; break;
    case shards::Wedge<15>::key:             ref = refNodeData_.wedge_15; break;
    case shards::Wedge<18>::key:             ref = refNodeData_.wedge_18; break;

    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::getReferenceNode): invalid cell topology.");
    }
    }
    
    // subview version; this is dangerous that users get control over the static data
    // cellNodes = Kokkos::subdynrankview( ref, nodeOrd, Kokkos::ALL() );

    // let's copy;
    const auto dim = cell.getDimension();
    for (auto i=0;i<dim;++i) 
      cellNodes(i) = ref(nodeOrd, i);
  }

  template<typename SpT>
  template<typename subcellNodeValueType, class ...subcellNodeProperties>
  void
  CellTools<SpT>::
  getReferenceSubcellNodes( /**/  Kokkos::DynRankView<subcellNodeValueType,subcellNodeProperties...> subcellNodes,
                            const ordinal_type         subcellDim,
                            const ordinal_type         subcellOrd,
                            const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !hasReferenceCell(parentCell), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): the specified cell topology does not have a reference cell.");

    INTREPID2_TEST_FOR_EXCEPTION( subcellDim > parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): subcell dimension out of range.");
    
    INTREPID2_TEST_FOR_EXCEPTION( subcellOrd >= parentCell.getSubcellCount(subcellDim), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): subcell ordinal out of range.");
    
    // Verify subcellNodes rank and dimensions
    INTREPID2_TEST_FOR_EXCEPTION( subcellNodes.rank() != 2, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): subcellNodes must have rank 2.");
      
    INTREPID2_TEST_FOR_EXCEPTION( subcellNodes.dimension(0) != parentCell.getNodeCount(subcellDim, subcellOrd), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): subcellNodes dimension (0) must match to node count in cell.");
    
    INTREPID2_TEST_FOR_EXCEPTION( subcellNodes.dimension(1) != parentCell.getDimension(), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): subcellNodes dimension (1) must match to cell dimension.");
#endif 
    
    // Find how many cellWorkset does the specified subcell have.
    const auto subcNodeCount = parentCell.getNodeCount(subcellDim, subcellOrd);
    
    // Loop over subcell cellWorkset
    for (auto subcNodeOrd=0;subcNodeOrd<subcNodeCount;++subcNodeOrd) {      
      // Get the node number relative to the parent reference cell
      const auto cellNodeOrd = parentCell.getNodeMap(subcellDim, subcellOrd, subcNodeOrd);

      auto dst = Kokkos::subdynrankview(subcellNodes, subcNodeOrd, Kokkos::ALL());
      getReferenceNode(dst, parentCell, cellNodeOrd);
    }
  }  

  template<typename SpT>  
  template<typename refEdgeTangentValueType, class ...refEdgeTangentProperties>
  void
  CellTools<SpT>::
  getReferenceEdgeTangent( /**/ Kokkos::DynRankView<refEdgeTangentValueType,refEdgeTangentProperties...> refEdgeTangent,
                           const ordinal_type         edgeOrd,
                           const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 2 &&
                                  parentCell.getDimension() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): two or three-dimensional parent cell required");
  
    INTREPID2_TEST_FOR_EXCEPTION( edgeOrd <  0 ||
                                  edgeOrd >= parentCell.getSubcellCount(1), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): edge ordinal out of bounds");  

    INTREPID2_TEST_FOR_EXCEPTION( refEdgeTangent.rank() != 1, std::invalid_argument,  
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): rank = 1 required for output arrays"); 
    
    INTREPID2_TEST_FOR_EXCEPTION( refEdgeTangent.dimension(0) != parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): output array size is required to match space dimension");  
#endif
    // Edge parametrizations are computed in setSubcellParametrization and stored in rank-3 array 
    // (subcOrd, coordinate, coefficient)
    subcellParamViewType edgeMap;
    getSubcellParametrization(edgeMap, 1, parentCell);
  
    // All ref. edge maps have affine coordinate functions: f_dim(u) = C_0(dim) + C_1(dim)*u, 
    //                                     => edge Tangent: -> C_1(*)
    const auto dim = parentCell.getDimension();
    for (auto i=0;i<dim;++i)
      refEdgeTangent(i) = edgeMap(edgeOrd, i, 1);
  }


  template<typename SpT>
  template<typename refFaceTanUValueType,...refFaceTanUProperties,
           typename refFaceTanVValueType,...refFaceTanVProperties>
  void
  CellTools<SpT>::
  getReferenceFaceTangents( /**/  Kokkos::DynRankView<refFaceTanUValueType,refFaceTanUProperties...> refFaceTanU,
                            /**/  Kokkos::DynRankView<refFaceTanVValueType,refFaceTanVProperties...> refFaceTanV,
                            const ordinal_type         faceOrd,
                            const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): three-dimensional parent cell required");  
  
    INTREPID2_TEST_FOR_EXCEPTION( faceOrd < 0 || faceOrd >= parentCell.getSubcellCount(2), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): face ordinal out of bounds");  
    
    INTREPID2_TEST_FOR_EXCEPTION( uTan.rank() != 1 || vTan.rank() != 1, std::invalid_argument,  
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): rank = 1 required for output arrays"); 
  
    INTREPID2_TEST_FOR_EXCEPTION( uTan.dimension(0) != parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): dim0 (spatial dim) must match that of parent cell");  

    INTREPID2_TEST_FOR_EXCEPTION( vTan.dimension(0) != parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): dim0 (spatial dim) must match that of parent cell");  
#endif
  
    // Face parametrizations are computed in setSubcellParametrization and stored in rank-3 array 
    // (subcOrd, coordinate, coefficient): retrieve this array
    subcellParamViewType faceMap;
    getSubcellParametrization(faceMap, 2, parentCell);
  
    /*  All ref. face maps have affine coordinate functions:  f_dim(u,v) = C_0(dim) + C_1(dim)*u + C_2(dim)*v
     *                           `   => Tangent vectors are:  uTan -> C_1(*);    vTan -> C_2(*)
     */

    // set uTan -> C_1(*)
    // set vTan -> C_2(*)
    const auto dim = parentCell.getDimension();
    for (auto i=0;i<dim;++i) {
      uTan(i) = faceMap(faceOrd, i, 1);
      vTan(i) = faceMap(faceOrd, i, 2);
    }
  }

  template<typename SpT>
  template<typename refSideNormalValueType, class ...refSideNormalProperties>
  void
  CellTools<SpT>::
  getReferenceSideNormal( /**/  Kokkos::DynRankView<refSideNormalValueType,refSideNormalProperties...> refSideNormal,
                          const ordinal_type         sideOrd,
                          const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 2 &&
                                  parentCell.getDimension() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSideNormal): two or three-dimensional parent cell required");
  
    // Check side ordinal: by definition side is subcell whose dimension = parentCell.getDimension()-1
    INTREPID2_TEST_FOR_EXCEPTION( sideOrd < 0 || sideOrd >= parentCell.getSubcellCount(parentCell.getDimension() - 1), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSideNormal): side ordinal out of bounds");    
#endif 

    const auto dim = parentCell.getDimension();
    if (dim == 2) {
      // 2D parent cells: side = 1D subcell (edge), call the edge tangent method and rotate tangents
      getReferenceEdgeTangent(refSideNormal, sideOrd, parentCell);
    
      // rotate t(t1, t2) to get n(t2, -t1) so that (n,t) is positively oriented: det(n1,n2/t1,t2)>0
      refSideNormalValueType tmp = refSideNormal(0);
      refSideNormal(0) = refSideNormal(1);
      refSideNormal(1) = -tmp;
    } else {
      // 3D parent cell: side = 2D subcell (face), call the face normal method.
      getReferenceFaceNormal(refSideNormal, sideOrd, parentCell);
    }
  }


  template<typename SpT>
  template<typename refFaceNormalValueType, class ...refFaceNormalProperties>>
  void 
  CellTools<SpT>::
  getReferenceFaceNormal( /**/  Kokkos::DynRankView<refFaceNormalValueType,refFaceNormalProperties...> refFaceNormal,
                          const ordinal_type         faceOrd,
                          const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceNormal): three-dimensional parent cell required");  
    
    INTREPID2_TEST_FOR_EXCEPTION( faceOrd < 0 || faceOrd >= parentCell.getSubcellCount(2), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceNormal): face ordinal out of bounds");  
    
    INTREPID2_TEST_FOR_EXCEPTION( refFaceNormal.rank() != 1, std::invalid_argument,  
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceNormal): rank = 1 required for output array"); 
  
    INTREPID2_TEST_FOR_EXCEPTION( refFaceNormal.dimension(0) != parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceNormal): dim0 (spatial dim) must match that of parent cell");  
#endif
    
    // Reference face normal = vector product of reference face tangents. Allocate temp FC storage:
    Kokkos::DynRankView<refFaceNormalValueType,SpT>
    FieldContainer<Scalar> uTan(spaceDim);
    FieldContainer<Scalar> vTan(spaceDim);
    getReferenceFaceTangents(uTan, vTan, faceOrd, parentCell);
  
    // Compute the vector product of the reference face tangents:
    RealSpaceTools<Scalar>::vecprod(refFaceNormal, uTan, vTan);
  }

//   template<class Scalar>
//   template<class ArrayEdgeTangent, class ArrayJac>
//   void CellTools<Scalar>::getPhysicalEdgeTangents(ArrayEdgeTangent &            edgeTangents,
//                                                   const ArrayJac &              worksetJacobians,
//                                                   const int &                   worksetEdgeOrd,
//                                                   const shards::CellTopology &  parentCell){
//     index_type worksetSize = static_cast<index_type>(worksetJacobians.dimension(0));
//     index_type edgePtCount = static_cast<index_type>(worksetJacobians.dimension(1)); 
//     int pCellDim    = parentCell.getDimension();
// #ifdef HAVE_INTREPID2_DEBUG
//     std::string errmsg = ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents):";
  
//     INTREPID2_TEST_FOR_EXCEPTION( !( (pCellDim == 3) || (pCellDim == 2) ), std::invalid_argument, 
//                                   ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): 2D or 3D parent cell required");  
  
//     // (1) edgeTangents is rank-3 (C,P,D) and D=2, or 3 is required
//     INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, edgeTangents, 3,3), std::invalid_argument, errmsg);
//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, edgeTangents, 2, 2,3), std::invalid_argument, errmsg);
 
//     // (2) worksetJacobians in rank-4 (C,P,D,D) and D=2, or 3 is required
//     INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, worksetJacobians, 4,4), std::invalid_argument, errmsg);
//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 2, 2,3), std::invalid_argument, errmsg);
//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 3, 2,3), std::invalid_argument, errmsg);
  
//     // (4) cross-check array dimensions: edgeTangents (C,P,D) vs. worksetJacobians (C,P,D,D)
//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, edgeTangents, 0,1,2,2,  worksetJacobians, 0,1,2,3), std::invalid_argument, errmsg);      
  
// #endif
  
//     // Temp storage for constant reference edge tangent: rank-1 (D) arrays
//     FieldContainer<double> refEdgeTan(pCellDim);
//     getReferenceEdgeTangent(refEdgeTan, worksetEdgeOrd, parentCell);
  
//     // Loop over workset faces and edge points
//     for(index_type pCell = 0; pCell < worksetSize; pCell++){
//       for(index_type pt = 0; pt < edgePtCount; pt++){
      
//         // Apply parent cell Jacobian to ref. edge tangent
//         for(int i = 0; i < pCellDim; i++){
//           edgeTangents(pCell, pt, i) = 0.0;
//           for(int j = 0; j < pCellDim; j++){
//             edgeTangents(pCell, pt, i) +=  worksetJacobians(pCell, pt, i, j)*refEdgeTan(j);
//           }// for j
//         }// for i
//       }// for pt
//     }// for pCell
//   }
//   template<class Scalar>
//   template<class ArrayFaceTangentU, class ArrayFaceTangentV, class ArrayJac>
//   void CellTools<Scalar>::getPhysicalFaceTangents(ArrayFaceTangentU &           faceTanU,
//                                                   ArrayFaceTangentV &           faceTanV,
//                                                   const ArrayJac &              worksetJacobians,
//                                                   const int &                   worksetFaceOrd,
//                                                   const shards::CellTopology &  parentCell){
//     index_type worksetSize = static_cast<index_type>(worksetJacobians.dimension(0));
//     index_type facePtCount = static_cast<index_type>(worksetJacobians.dimension(1)); 
//     int pCellDim    = parentCell.getDimension();
// #ifdef HAVE_INTREPID2_DEBUG
//     std::string errmsg = ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents):";

//     INTREPID2_TEST_FOR_EXCEPTION( !(pCellDim == 3), std::invalid_argument, 
//                                   ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): three-dimensional parent cell required");  
  
//     // (1) faceTanU and faceTanV are rank-3 (C,P,D) and D=3 is required
//     INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, faceTanU, 3,3), std::invalid_argument, errmsg);
//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, faceTanU, 2, 3,3), std::invalid_argument, errmsg);

//     INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, faceTanV, 3,3), std::invalid_argument, errmsg);
//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, faceTanV, 2, 3,3), std::invalid_argument, errmsg);

//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, faceTanU,  faceTanV), std::invalid_argument, errmsg);      

//     // (3) worksetJacobians in rank-4 (C,P,D,D) and D=3 is required
//     INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, worksetJacobians, 4,4), std::invalid_argument, errmsg);
//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 2, 3,3), std::invalid_argument, errmsg);
//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 3, 3,3), std::invalid_argument, errmsg);

//     // (4) cross-check array dimensions: faceTanU (C,P,D) vs. worksetJacobians (C,P,D,D)
//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, faceTanU, 0,1,2,2,  worksetJacobians, 0,1,2,3), std::invalid_argument, errmsg);      

// #endif
    
//     // Temp storage for the pair of constant ref. face tangents: rank-1 (D) arrays
//     FieldContainer<double> refFaceTanU(pCellDim);
//     FieldContainer<double> refFaceTanV(pCellDim);
//     getReferenceFaceTangents(refFaceTanU, refFaceTanV, worksetFaceOrd, parentCell);

//     // Loop over workset faces and face points
//     for(index_type pCell = 0; pCell < worksetSize; pCell++){
//       for(index_type pt = 0; pt < facePtCount; pt++){
      
//         // Apply parent cell Jacobian to ref. face tangents
//         for(int dim = 0; dim < pCellDim; dim++){
//           faceTanU(pCell, pt, dim) = 0.0;
//           faceTanV(pCell, pt, dim) = 0.0;
        
//           // Unroll loops: parent cell dimension can only be 3
//           faceTanU(pCell, pt, dim) = \
//             worksetJacobians(pCell, pt, dim, 0)*refFaceTanU(0) + \
//             worksetJacobians(pCell, pt, dim, 1)*refFaceTanU(1) + \
//             worksetJacobians(pCell, pt, dim, 2)*refFaceTanU(2);
//           faceTanV(pCell, pt, dim) = \
//             worksetJacobians(pCell, pt, dim, 0)*refFaceTanV(0) + \
//             worksetJacobians(pCell, pt, dim, 1)*refFaceTanV(1) + \
//             worksetJacobians(pCell, pt, dim, 2)*refFaceTanV(2);
//         }// for dim
//       }// for pt
//     }// for pCell
//   }

//   template<class Scalar>
//   template<class ArraySideNormal, class ArrayJac>
//   void CellTools<Scalar>::getPhysicalSideNormals(ArraySideNormal &             sideNormals,
//                                                  const ArrayJac &              worksetJacobians,
//                                                  const int &                   worksetSideOrd,
//                                                  const shards::CellTopology &  parentCell){
//     index_type worksetSize = static_cast<index_type>(worksetJacobians.dimension(0));
//     index_type sidePtCount = static_cast<index_type>(worksetJacobians.dimension(1));   
//     int spaceDim  = parentCell.getDimension();
// #ifdef HAVE_INTREPID2_DEBUG
//     INTREPID2_TEST_FOR_EXCEPTION( !( (spaceDim == 2) || (spaceDim == 3) ), std::invalid_argument, 
//                                   ">>> ERROR (Intrepid2::CellTools::getPhysicalSideNormals): two or three-dimensional parent cell required");
  
//     // Check side ordinal: by definition side is subcell whose dimension = spaceDim-1
//     INTREPID2_TEST_FOR_EXCEPTION( !( (0 <= worksetSideOrd) && (worksetSideOrd < (int)parentCell.getSubcellCount(spaceDim - 1) ) ), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::getPhysicalSideNormals): side ordinal out of bounds");  
// #endif  
  
//     if(spaceDim == 2){

//       // 2D parent cells: side = 1D subcell (edge), call the edge tangent method and rotate tangents
//       getPhysicalEdgeTangents(sideNormals, worksetJacobians, worksetSideOrd, parentCell);
    
//       // rotate t(t1, t2) to get n(t2, -t1) so that (n,t) is positively oriented: det(n1,n2/t1,t2)>0
//       for(index_type cell = 0; cell < worksetSize; cell++){
//         for(index_type pt = 0; pt < sidePtCount; pt++){
//           Scalar temp = sideNormals(cell, pt, 0);
//           sideNormals(cell, pt, 0) = sideNormals(cell, pt, 1);
//           sideNormals(cell, pt, 1) = -temp;
//         }// for pt
//       }// for cell
//     }
//     else{
//       // 3D parent cell: side = 2D subcell (face), call the face normal method.
//       getPhysicalFaceNormals(sideNormals, worksetJacobians, worksetSideOrd, parentCell);
//     }
//   }
  
  
//   template<class Scalar>
//   template<class ArrayFaceNormal, class ArrayJac>
//   void CellTools<Scalar>::getPhysicalFaceNormals(ArrayFaceNormal &             faceNormals,
//                                                  const ArrayJac &              worksetJacobians,
//                                                  const int &                   worksetFaceOrd,
//                                                  const shards::CellTopology &  parentCell){
//     index_type worksetSize = static_cast<index_type>(worksetJacobians.dimension(0));
//     index_type facePtCount = static_cast<index_type>(worksetJacobians.dimension(1)); 
//     int pCellDim    = parentCell.getDimension();
// #ifdef HAVE_INTREPID2_DEBUG
//     std::string errmsg = ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals):";
  
//     INTREPID2_TEST_FOR_EXCEPTION( !(pCellDim == 3), std::invalid_argument, 
//                                   ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): three-dimensional parent cell required");  
  
//     // (1) faceNormals is rank-3 (C,P,D) and D=3 is required
//     INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, faceNormals, 3,3), std::invalid_argument, errmsg);
//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, faceNormals, 2, 3,3), std::invalid_argument, errmsg);
  
//     // (3) worksetJacobians in rank-4 (C,P,D,D) and D=3 is required
//     INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, worksetJacobians, 4,4), std::invalid_argument, errmsg);
//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 2, 3,3), std::invalid_argument, errmsg);
//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 3, 3,3), std::invalid_argument, errmsg);
  
//     // (4) cross-check array dimensions: faceNormals (C,P,D) vs. worksetJacobians (C,P,D,D)
//     INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, faceNormals, 0,1,2,2,  worksetJacobians, 0,1,2,3), std::invalid_argument, errmsg);        
// #endif
  
//     // Temp storage for physical face tangents: rank-3 (C,P,D) arrays
//     FieldContainer<Scalar> faceTanU(worksetSize, facePtCount, pCellDim);
//     FieldContainer<Scalar> faceTanV(worksetSize, facePtCount, pCellDim);
//     getPhysicalFaceTangents(faceTanU, faceTanV, worksetJacobians, worksetFaceOrd, parentCell);
  
//     // Compute the vector product of the physical face tangents:
//     RealSpaceTools<Scalar>::vecprod(faceNormals, faceTanU, faceTanV);
  
  
//   }


  template<typename SpT>
  bool 
  CellTools<SpT>::
  isReferenceNodeDataSet_ = false;

  template<typename SpT>
  typename CellTools<SpT>::ReferenceNodeData
  CellTools<SpT>::
  refNodeData_ = typename CellTools<SpT>::ReferenceNodeData();

  template<typename SpT>
  const typename CellTools<SpT>::ReferenceNodeDataStatic
  CellTools<SpT>::
  refNodeDataStatic_ = {    
    // line
    { // 2
      {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0} 
    },
    { // 3
      {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 0.0, 0.0}
    },
    // triangle
    { // 3
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0} 
    },
    { // 4
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 1/3, 1/3, 0.0}
    },
    { // 6
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
      { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}
    },
    // quad
    { // 4
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}
    },
    { // 8
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
      { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}
    },
    { // 9
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
      { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}, { 0.0, 0.0, 0.0}
    },
    // tet
    { // 4
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0}
    },
    { // 8
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      { 1/3, 0.0, 1/3}, { 1/3, 1/3, 1/3}, { 1/3, 1/3, 0.0}, { 0.0, 1/3, 1/3} 
    },
    { // 10
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}, { 0.0, 0.0, 0.5}, { 0.5, 0.0, 0.5}, { 0.0, 0.5, 0.5}
    },
    { // 11
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}, { 0.0, 0.0, 0.5}, { 0.5, 0.0, 0.5}, { 0.0, 0.5, 0.5}
    },
    // hex
    { // 8
      {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
      {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}
    },
    { // 20
      {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
      {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0},
      { 0.0,-1.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, {-1.0, 0.0,-1.0}, 
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
      { 0.0,-1.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}, {-1.0, 0.0, 1.0}
    },
    { // 27
      {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
      {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0},
      { 0.0,-1.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, {-1.0, 0.0,-1.0}, 
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
      { 0.0,-1.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}, {-1.0, 0.0, 1.0},
      { 0.0, 0.0, 0.0},
      { 0.0, 0.0,-1.0}, { 0.0, 0.0, 1.0}, {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, {0.0,-1.0, 0.0}, {0.0, 1.0, 0.0} 
    },
    // pyramid
    { // 5
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0}
    },
    { // 13
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0},
      {-0.5,-0.5, 0.5}, { 0.5,-0.5, 0.5}, { 0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}   
    },
    { // 14 
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0},
      {-0.5,-0.5, 0.5}, { 0.5,-0.5, 0.5}, { 0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}, { 0.0, 0.0, 0.0}  
    },
    // wedge
    { // 6
      { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0} 
    },
    { // 15
      { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0},
      { 0.5, 0.0,-1.0}, { 0.5, 0.5,-1.0}, { 0.0, 0.5,-1.0}, { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
      { 0.5, 0.0, 1.0}, { 0.5, 0.5, 1.0}, { 0.0, 0.5, 1.0}
    },
    { // 18
      { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0},
      { 0.5, 0.0,-1.0}, { 0.5, 0.5,-1.0}, { 0.0, 0.5,-1.0}, { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
      { 0.5, 0.0, 1.0}, { 0.5, 0.5, 1.0}, { 0.0, 0.5, 1.0},
      { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}
    }
  };
    
}

#endif
