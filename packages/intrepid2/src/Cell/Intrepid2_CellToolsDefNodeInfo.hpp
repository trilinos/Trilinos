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

    Kokkos::push_finalize_hook( [=] {
      refNodeData_.line            = referenceNodeDataViewType();
      refNodeData_.line_3          = referenceNodeDataViewType();

      refNodeData_.triangle        = referenceNodeDataViewType();
      refNodeData_.triangle_4      = referenceNodeDataViewType();
      refNodeData_.triangle_6      = referenceNodeDataViewType();

      refNodeData_.quadrilateral   = referenceNodeDataViewType();
      refNodeData_.quadrilateral_8 = referenceNodeDataViewType();
      refNodeData_.quadrilateral_9 = referenceNodeDataViewType();

      refNodeData_.tetrahedron     = referenceNodeDataViewType();
      refNodeData_.tetrahedron_8   = referenceNodeDataViewType();
      refNodeData_.tetrahedron_10  = referenceNodeDataViewType();
      refNodeData_.tetrahedron_11  = referenceNodeDataViewType();

      refNodeData_.hexahedron      = referenceNodeDataViewType();
      refNodeData_.hexahedron_20   = referenceNodeDataViewType();
      refNodeData_.hexahedron_27   = referenceNodeDataViewType();

      refNodeData_.pyramid         = referenceNodeDataViewType();
      refNodeData_.pyramid_13      = referenceNodeDataViewType();
      refNodeData_.pyramid_14      = referenceNodeDataViewType();

      refNodeData_.wedge           = referenceNodeDataViewType();
      refNodeData_.wedge_15        = referenceNodeDataViewType();
      refNodeData_.wedge_18        = referenceNodeDataViewType();
    } );

    isReferenceNodeDataSet_ = true;
  }

  template<typename SpT>
  template<typename cellCenterValueType, class ...cellCenterProperties,
           typename cellVertexValueType, class ...cellVertexProperties>
  void 
  CellTools<SpT>::
  getReferenceCellCenter( Kokkos::DynRankView<cellCenterValueType,cellCenterProperties...> cellCenter,
                          Kokkos::DynRankView<cellVertexValueType,cellVertexProperties...> cellVertex,
                          const shards::CellTopology cell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !hasReferenceCell(cell), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceCellCenter): the specified cell topology does not have a reference cell." );
    
    INTREPID2_TEST_FOR_EXCEPTION( cellCenter.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceCellCenter): cellCenter must have rank 1." );

    INTREPID2_TEST_FOR_EXCEPTION( cellCenter.dimension(0) < cell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceCellCenter): cellCenter must have dimension bigger than Parameters::MaxDimension." );

    INTREPID2_TEST_FOR_EXCEPTION( cellVertex.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceCellCenter): cellVertex must have rank 1." );

    INTREPID2_TEST_FOR_EXCEPTION( cellVertex.dimension(0) < cell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceCellCenter): cellVertex must have dimension bigger than Parameters::MaxDimension." );
#endif
    const ordinal_type vertexCount = cell.getVertexCount();
    const ordinal_type dim = cell.getDimension();

    for (ordinal_type i=0;i<dim;++i) {
      cellCenter(i) = 0;
      for (ordinal_type vertOrd=0;vertOrd<vertexCount;++vertOrd) {
        getReferenceVertex(cellVertex, cell, vertOrd); 
        cellCenter(i) += cellVertex(i);
      }
      cellCenter(i) /= vertexCount;
    }
  }


  template<typename SpT>
  template<typename cellVertexValueType, class ...cellVertexProperties>
  void
  CellTools<SpT>::
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

    INTREPID2_TEST_FOR_EXCEPTION( cellVertex.dimension(0) < cell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): cellNodes must have dimension bigger than Parameters::MaxDimension." );
#endif
    getReferenceNode(cellVertex, 
                     cell, 
                     vertexOrd);
  }
    
  
  template<typename SpT>
  template<typename subcellVertexValueType, class ...subcellVertexProperties>
  void
  CellTools<SpT>::
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
  getReferenceNode(       Kokkos::DynRankView<cellNodeValueType,cellNodeProperties...> cellNode,
                    const shards::CellTopology  cell,
                    const ordinal_type          nodeOrd ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !hasReferenceCell(cell), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): the specified cell topology does not have a reference cell." );
    
    INTREPID2_TEST_FOR_EXCEPTION( nodeOrd >= static_cast<ordinal_type>(cell.getNodeCount()), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): invalid node ordinal for the specified cell topology." );

    INTREPID2_TEST_FOR_EXCEPTION( cellNode.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): cellNode must have rank 1." );

    INTREPID2_TEST_FOR_EXCEPTION( cellNode.dimension(0) < cell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): cellNode must have dimension bigger than Parameters::MaxDimension." );
#endif

    if (!isReferenceNodeDataSet_) 
      setReferenceNodeData();

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
    // cellNode = Kokkos::subdynrankview( ref, nodeOrd, Kokkos::ALL() );

    // let's copy;
    const ordinal_type dim = cell.getDimension();

    for (ordinal_type i=0;i<dim;++i) 
      cellNode(i) = ref(nodeOrd, i);
  }

  template<typename SpT>
  template<typename subcellNodeValueType, class ...subcellNodeProperties>
  void
  CellTools<SpT>::
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
      
    INTREPID2_TEST_FOR_EXCEPTION( subcellNodes.dimension(0) != parentCell.getNodeCount(subcellDim, subcellOrd), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): subcellNodes dimension (0) must match to node count in cell.");
    
    INTREPID2_TEST_FOR_EXCEPTION( subcellNodes.dimension(1) != parentCell.getDimension(), std::invalid_argument, 
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

  template<typename SpT>  
  template<typename refEdgeTangentValueType, class ...refEdgeTangentProperties>
  void
  CellTools<SpT>::
  getReferenceEdgeTangent(       Kokkos::DynRankView<refEdgeTangentValueType,refEdgeTangentProperties...> refEdgeTangent,
                           const ordinal_type         edgeOrd,
                           const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 2 &&
                                  parentCell.getDimension() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): two or three-dimensional parent cell required");
  
    INTREPID2_TEST_FOR_EXCEPTION( refEdgeTangent.rank() != 1, std::invalid_argument,  
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): rank = 1 required for output arrays"); 
    
    INTREPID2_TEST_FOR_EXCEPTION( refEdgeTangent.dimension(0) != parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): output array size is required to match space dimension");  

    INTREPID2_TEST_FOR_EXCEPTION( edgeOrd <  0 ||
                                  edgeOrd >= static_cast<ordinal_type>(parentCell.getSubcellCount(1)), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): edge ordinal out of bounds");  

#endif
    // Edge parametrizations are computed in setSubcellParametrization and stored in rank-3 array 
    // (subcOrd, coordinate, coefficient)
    subcellParamViewType edgeMap;
    getSubcellParametrization(edgeMap, 1, parentCell);
  
    // All ref. edge maps have affine coordinate functions: f_dim(u) = C_0(dim) + C_1(dim)*u, 
    //                                     => edge Tangent: -> C_1(*)
    const ordinal_type dim = parentCell.getDimension();
    for (ordinal_type i=0;i<dim;++i)
      refEdgeTangent(i) = edgeMap(edgeOrd, i, 1);
  }


  template<typename SpT>
  template<typename refFaceTanUValueType, class ...refFaceTanUProperties,
           typename refFaceTanVValueType, class ...refFaceTanVProperties>
  void
  CellTools<SpT>::
  getReferenceFaceTangents(       Kokkos::DynRankView<refFaceTanUValueType,refFaceTanUProperties...> refFaceTanU,
                                  Kokkos::DynRankView<refFaceTanVValueType,refFaceTanVProperties...> refFaceTanV,
                            const ordinal_type         faceOrd,
                            const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( parentCell.getDimension() != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): three-dimensional parent cell required");  
  
    INTREPID2_TEST_FOR_EXCEPTION( faceOrd < 0 || faceOrd >= static_cast<ordinal_type>(parentCell.getSubcellCount(2)), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): face ordinal out of bounds");  
    
    INTREPID2_TEST_FOR_EXCEPTION( refFaceTanU.rank() != 1 || refFaceTanV.rank() != 1, std::invalid_argument,  
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): rank = 1 required for output arrays"); 
  
    INTREPID2_TEST_FOR_EXCEPTION( refFaceTanU.dimension(0) != parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): dim0 (spatial dim) must match that of parent cell");  

    INTREPID2_TEST_FOR_EXCEPTION( refFaceTanV.dimension(0) != parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): dim0 (spatial dim) must match that of parent cell");  
#endif
  
    // Face parametrizations are computed in setSubcellParametrization and stored in rank-3 array 
    // (subcOrd, coordinate, coefficient): retrieve this array
    subcellParamViewType faceMap;
    getSubcellParametrization(faceMap, 2, parentCell);
  
    /*  All ref. face maps have affine coordinate functions:  f_dim(u,v) = C_0(dim) + C_1(dim)*u + C_2(dim)*v
     *                           `   => Tangent vectors are:  refFaceTanU -> C_1(*);    refFaceTanV -> C_2(*)
     */

    // set refFaceTanU -> C_1(*)
    // set refFaceTanV -> C_2(*)
    const ordinal_type dim = parentCell.getDimension();
    for (ordinal_type i=0;i<dim;++i) {
      refFaceTanU(i) = faceMap(faceOrd, i, 1);
      refFaceTanV(i) = faceMap(faceOrd, i, 2);
    }
  }

  template<typename SpT>
  template<typename refSideNormalValueType, class ...refSideNormalProperties>
  void
  CellTools<SpT>::
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
  template<typename refFaceNormalValueType, class ...refFaceNormalProperties>
  void 
  CellTools<SpT>::
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
  
    INTREPID2_TEST_FOR_EXCEPTION( refFaceNormal.dimension(0) != parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceNormal): dim0 (spatial dim) must match that of parent cell");  
#endif
    
    // Reference face normal = vector product of reference face tangents. Allocate temp FC storage:
    const auto dim = parentCell.getDimension();
    auto vcprop = Kokkos::common_view_alloc_prop(refFaceNormal);
    using common_value_type = typename decltype(vcprop)::value_type;
    Kokkos::DynRankView< common_value_type, SpT > refFaceTanU ( Kokkos::view_alloc("CellTools::getReferenceFaceNormal::refFaceTanU", vcprop), dim );
    Kokkos::DynRankView< common_value_type, SpT > refFaceTanV ( Kokkos::view_alloc("CellTools::getReferenceFaceNormal::refFaceTanV", vcprop), dim );
    getReferenceFaceTangents(refFaceTanU, refFaceTanV, faceOrd, parentCell);
  
    RealSpaceTools<SpT>::vecprod(refFaceNormal, refFaceTanU, refFaceTanV);
  }

  template<typename SpT>
  template<typename edgeTangentValueType,     class ...edgeTangentProperties,
           typename worksetJacobianValueType, class ...worksetJacobianProperties>
  void
  CellTools<SpT>::
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
    INTREPID2_TEST_FOR_EXCEPTION( edgeTangents.dimension(2) != 2 && 
                                  edgeTangents.dimension(2) != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): edgeTangents dimension(2) must be 2 or 3." );
 
    // (2) worksetJacobians in rank-4 (C,P,D,D) and D=2, or 3 is required
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.rank() != 4, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): worksetJacobians requires rank 4." );  
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.dimension(2) != 2 && 
                                  worksetJacobians.dimension(2) != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): worksetJacobians dimension(2) must be 2 or 3." );
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.dimension(2) != worksetJacobians.dimension(3), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): worksetJacobians dimension(2) and (3) must match each other." );

    // (4) cross-check array dimensions: edgeTangents (C,P,D) vs. worksetJacobians (C,P,D,D)
    for (auto i=0;i<3;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( edgeTangents.dimension(i) != worksetJacobians.dimension(i), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): edgeTangents dimension (i) does not match to worksetJacobians dimension(i)." );
    }
#endif
  
    // Storage for constant reference edge tangent: rank-1 (D) arrays
    const auto dim = parentCell.getDimension();
    auto vcprop = Kokkos::common_view_alloc_prop(edgeTangents);
    using common_value_type = typename decltype(vcprop)::value_type;
    Kokkos::DynRankView< common_value_type, SpT > refEdgeTan ( Kokkos::view_alloc("CellTools::getPhysicalEdgeTangents::refEdgeTan", vcprop), dim);
    getReferenceEdgeTangent(refEdgeTan, worksetEdgeOrd, parentCell);
    
    RealSpaceTools<SpT>::matvec(edgeTangents, worksetJacobians, refEdgeTan);
  }


  template<typename SpT>
  template<typename faceTanUValueType,        class ...faceTanUProperties,
           typename faceTanVValueType,        class ...faceTanVProperties,
           typename worksetJacobianValueType, class ...worksetJacobianProperties>
  void
  CellTools<SpT>::
  getPhysicalFaceTangents(       Kokkos::DynRankView<faceTanUValueType,faceTanUProperties...> faceTanU,
                                 Kokkos::DynRankView<faceTanVValueType,faceTanVProperties...> faceTanV,
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

    INTREPID2_TEST_FOR_EXCEPTION( faceTanU.dimension(2) != 3 ||
                                  faceTanV.dimension(2) != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): faceTan U,V dimension (2) must be 3." );  
    
    for (auto i=0;i<3;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( faceTanU.dimension(i) != faceTanV.dimension(i), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): faceTan U,V dimension (i) must match each other." );  
    }

    // (3) worksetJacobians in rank-4 (C,P,D,D) and D=3 is required
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.rank() != 4, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): worksetJacobians must have rank 4." );  

    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.dimension(2) != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): worksetJacobians dimension(2) must be 3." );  

    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.dimension(2) != worksetJacobians.dimension(3), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): worksetJacobians dimension(2) and dimension(3) must match." );  

    // (4) cross-check array dimensions: faceTanU (C,P,D) vs. worksetJacobians (C,P,D,D)
    for (auto i=0;i<3;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( faceTanU.dimension(i) != worksetJacobians.dimension(i), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): worksetJacobians dimension(i) and faceTan dimension (i) must match." );  
    }      
#endif
    
    // Temp storage for the pair of constant ref. face tangents: rank-1 (D) arrays
    const auto dim = parentCell.getDimension();

    auto vcpropU = Kokkos::common_view_alloc_prop(faceTanU);
    using common_value_typeU = typename decltype(vcpropU)::value_type;
    Kokkos::DynRankView< common_value_typeU, SpT > refFaceTanU ( Kokkos::view_alloc("CellTools::getPhysicalFaceTangents::refFaceTanU", vcpropU), dim);

    auto vcpropV = Kokkos::common_view_alloc_prop(faceTanV);
    using common_value_typeV = typename decltype(vcpropV)::value_type;
    Kokkos::DynRankView< common_value_typeV, SpT > refFaceTanV ( Kokkos::view_alloc("CellTools::getPhysicalFaceTangents::refFaceTanV", vcpropV), dim);

    getReferenceFaceTangents(refFaceTanU, refFaceTanV, worksetFaceOrd, parentCell);

    RealSpaceTools<SpT>::matvec(faceTanU, worksetJacobians, refFaceTanU);    
    RealSpaceTools<SpT>::matvec(faceTanV, worksetJacobians, refFaceTanV);    
  }


  template<typename SpT>
  template<typename sideNormalValueType,      class ...sideNormalProperties,
           typename worksetJacobianValueType, class ...worksetJacobianProperties>
  void 
  CellTools<SpT>::
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
    const auto dim = parentCell.getDimension();
  
    if (dim == 2) {
      // compute edge tangents and rotate it
      auto vcprop = Kokkos::common_view_alloc_prop(sideNormals);
      using common_value_type = typename decltype(vcprop)::value_type;
      Kokkos::DynRankView< common_value_type, SpT > edgeTangents ( Kokkos::view_alloc("CellTools::getPhysicalSideNormals::edgeTan", vcprop),
                                                              sideNormals.dimension(0),
                                                              sideNormals.dimension(1),
                                                              sideNormals.dimension(2));
      getPhysicalEdgeTangents(edgeTangents, worksetJacobians, worksetSideOrd, parentCell);

      Kokkos::DynRankView< common_value_type, SpT > rotation ( Kokkos::view_alloc("CellTools::getPhysicalSideNormals::rotation", vcprop), dim, dim);
      rotation(0,0) =  0; rotation(0,1) =  1;
      rotation(1,0) = -1; rotation(1,1) =  0;

      RealSpaceTools<SpT>::matvec(sideNormals, rotation, edgeTangents);    
    } else {
      getPhysicalFaceNormals(sideNormals, worksetJacobians, worksetSideOrd, parentCell);
    }
  }
  

  template<typename SpT>
  template<typename faceNormalValueType,      class ...faceNormalProperties,
           typename worksetJacobianValueType, class ...worksetJacobianProperties>
  void
  CellTools<SpT>::
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
    INTREPID2_TEST_FOR_EXCEPTION( faceNormals.dimension(2) != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): faceNormals dimension (2) must be 3." );
    
    // (3) worksetJacobians in rank-4 (C,P,D,D) and D=3 is required
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.rank() != 4, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): worksetJacobians must have a rank 4." );
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.dimension(2) != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): worksetJacobians dimension (2) must be 3." );
    INTREPID2_TEST_FOR_EXCEPTION( worksetJacobians.dimension(2) != worksetJacobians.dimension(3), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): worksetJacobians dimension (2) must match to dimension (3)." );
  
    // (4) cross-check array dimensions: faceNormals (C,P,D) vs. worksetJacobians (C,P,D,D)
    for (auto i=0;i<3;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( faceNormals.dimension(i) != worksetJacobians.dimension(i), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): faceNormals dimension (i) must match to worksetJacobians dimension (i)." );
    }        
#endif
  
    // this should be provided from users
    // Storage for physical face tangents: rank-3 (C,P,D) arrays
    const auto worksetSize = worksetJacobians.dimension(0);
    const auto facePtCount = worksetJacobians.dimension(1);
    const auto dim = parentCell.getDimension();

    auto vcprop = Kokkos::common_view_alloc_prop(faceNormals);
    using common_value_type = typename decltype(vcprop)::value_type;
    Kokkos::DynRankView< common_value_type, SpT > faceTanU ( Kokkos::view_alloc("CellTools::getPhysicalFaceNormals::faceTanU", vcprop), worksetSize, facePtCount, dim);
    Kokkos::DynRankView< common_value_type, SpT > faceTanV ( Kokkos::view_alloc("CellTools::getPhysicalFaceNormals::faceTanV", vcprop), worksetSize, facePtCount, dim);

    getPhysicalFaceTangents(faceTanU, faceTanV, 
                            worksetJacobians, 
                            worksetFaceOrd, 
                            parentCell);
  
    RealSpaceTools<SpT>::vecprod(faceNormals, faceTanU, faceTanV);
  }


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
