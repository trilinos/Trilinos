1;2c// @HEADER
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
  
  template<class Scalar>
  const double* CellTools<Scalar>::getReferenceVertex(const shards::CellTopology& cell,
                                                      const int                   vertexOrd){
    
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !(hasReferenceCell(cell) ), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceVertex): the specified cell topology does not have a reference cell.");
    
    INTREPID2_TEST_FOR_EXCEPTION( !( (0 <= vertexOrd) && (vertexOrd < (int)cell.getVertexCount() ) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceVertex): invalid node ordinal for the specified cell topology. ");
#endif
    
    // Simply call getReferenceNode with the base topology of the cell
    return getReferenceNode(cell.getBaseCellTopologyData(), vertexOrd);
  }
    
  
  
  template<class Scalar>
  template<class ArraySubcellVert>
  void CellTools<Scalar>::getReferenceSubcellVertices(ArraySubcellVert &          subcellVertices,
                                                      const int                   subcellDim,
                                                      const int                   subcellOrd,
                                                      const shards::CellTopology& parentCell){
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !(hasReferenceCell(parentCell) ), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): the specified cell topology does not have a reference cell.");

    // subcellDim can equal the cell dimension because the cell itself is a valid subcell! In this case
    // the method will return all cell cellWorkset.
    INTREPID2_TEST_FOR_EXCEPTION( !( (0 <= subcellDim) && (subcellDim <= (int)parentCell.getDimension()) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): subcell dimension out of range.");
    
    INTREPID2_TEST_FOR_EXCEPTION( !( (0 <= subcellOrd) && (subcellOrd < (int)parentCell.getSubcellCount(subcellDim) ) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices): subcell ordinal out of range.");
        
    // Verify subcellVertices rank and dimensions
    {
      std::string errmsg = ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellVertices):";
      INTREPID2_TEST_FOR_EXCEPTION( !( requireRankRange(errmsg, subcellVertices, 2, 2) ), std::invalid_argument, errmsg);
      
      int subcVertexCount = parentCell.getVertexCount(subcellDim, subcellOrd);
      int spaceDim = parentCell.getDimension();
        
      INTREPID2_TEST_FOR_EXCEPTION( !( requireDimensionRange(errmsg, subcellVertices, 0,  subcVertexCount, subcVertexCount) ),
                                    std::invalid_argument, errmsg);
      
      INTREPID2_TEST_FOR_EXCEPTION( !( requireDimensionRange(errmsg, subcellVertices, 1,  spaceDim, spaceDim) ),
                                    std::invalid_argument, errmsg);
    }
#endif 
    
    // Simply call getReferenceNodes with the base topology
    getReferenceSubcellNodes(subcellVertices, subcellDim, subcellOrd, parentCell.getBaseCellTopologyData() );
  }  

  
  
  template<class Scalar>
  const double* CellTools<Scalar>::getReferenceNode(const shards::CellTopology& cell,
                                                    const int                   nodeOrd){
    
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !(hasReferenceCell(cell) ), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): the specified cell topology does not have a reference cell.");

    INTREPID2_TEST_FOR_EXCEPTION( !( (0 <= nodeOrd) && (nodeOrd < (int)cell.getNodeCount() ) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): invalid node ordinal for the specified cell topology. ");
#endif
    
    // Cartesian coordinates of supported reference cell cellWorkset, padded to three-dimensions.
    // Node order follows cell topology definition in Shards
    static const double line[2][3] ={
      {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0} 
    };
    static const double line_3[3][3] = {
      {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0},     
      // Extension node: edge midpoint
      { 0.0, 0.0, 0.0}
    };
    
    
    // Triangle topologies
    static const double triangle[3][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0} 
    };
    static const double triangle_4[4][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
      // Extension node: cell center
      { 1/3, 1/3, 0.0}
    };
    static const double triangle_6[6][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
      // Extension cellWorkset: 3 edge midpoints
      { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}
    };
    
    
    // Quadrilateral topologies
    static const double quadrilateral[4][3] = {
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}
    };
    static const double quadrilateral_8[8][3] = {
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
      // Extension cellWorkset: 4 edge midpoints
      { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}
    };
    static const double quadrilateral_9[9][3] = {
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
      // Extension cellWorkset: 4 edge midpoints + 1 cell center
      { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}, { 0.0, 0.0, 0.0}
    };
    
    
    // Tetrahedron topologies
    static const double tetrahedron[4][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0}
    };
    static const double tetrahedron_8[8][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      // Extension cellWorkset: 4 face centers (do not follow natural face order - see the cell topology!)
      { 1/3, 0.0, 1/3}, { 1/3, 1/3, 1/3}, { 1/3, 1/3, 0.0}, { 0.0, 1/3, 1/3} 
    };
    static const double tetrahedron_10[10][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      // Extension cellWorkset: 6 edge midpoints
      { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}, { 0.0, 0.0, 0.5}, { 0.5, 0.0, 0.5}, { 0.0, 0.5, 0.5}
    };

    static const double tetrahedron_11[10][3] = {
      { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      // Extension cellWorkset: 6 edge midpoints
      { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}, { 0.0, 0.0, 0.5}, { 0.5, 0.0, 0.5}, { 0.0, 0.5, 0.5}
    };

    
    // Hexahedron topologies
    static const double hexahedron[8][3] = {
      {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
      {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}
    };
    static const double hexahedron_20[20][3] = {
      {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
      {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0},
      // Extension cellWorkset: 12 edge midpoints (do not follow natural edge order - see cell topology!)
      { 0.0,-1.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, {-1.0, 0.0,-1.0}, 
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
      { 0.0,-1.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}, {-1.0, 0.0, 1.0}
    };
    static const double hexahedron_27[27][3] = {
      {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
      {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0},
      // Extension cellWorkset: 12 edge midpoints + 1 cell center + 6 face centers  (do not follow natural subcell order!)
      { 0.0,-1.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, {-1.0, 0.0,-1.0}, 
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
      { 0.0,-1.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}, {-1.0, 0.0, 1.0},
      { 0.0, 0.0, 0.0},
      { 0.0, 0.0,-1.0}, { 0.0, 0.0, 1.0}, {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, {0.0,-1.0, 0.0}, {0.0, 1.0, 0.0} 
    };
    
    
    // Pyramid topologies
    static const double pyramid[5][3] = {
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0}
    };
    static const double pyramid_13[13][3] = {
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      // Extension cellWorkset: 8 edge midpoints 
      { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0},
      {-0.5,-0.5, 0.5}, { 0.5,-0.5, 0.5}, { 0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}   
    };
    static const double pyramid_14[14][3] = {
      {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
      // Extension cellWorkset: 8 edge midpoints + quadrilateral face midpoint 
      { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0},
      {-0.5,-0.5, 0.5}, { 0.5,-0.5, 0.5}, { 0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}, { 0.0, 0.0, 0.0}  
    };
    
    
    // Wedge topologies
    static const double wedge[6][3] = {
      { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0} 
    };
    static const double wedge_15[15][3] = {
      { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0},
      // Extension cellWorkset: 9 edge midpoints (do not follow natural edge order - see cell topology!)
      { 0.5, 0.0,-1.0}, { 0.5, 0.5,-1.0}, { 0.0, 0.5,-1.0}, { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
      { 0.5, 0.0, 1.0}, { 0.5, 0.5, 1.0}, { 0.0, 0.5, 1.0}
    };
    static const double wedge_18[18][3] = {
      { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0},
      // Extension cellWorkset: 9 edge midpoints + 3 quad face centers (do not follow natural subcell order - see cell topology!)
      { 0.5, 0.0,-1.0}, { 0.5, 0.5,-1.0}, { 0.0, 0.5,-1.0}, { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
      { 0.5, 0.0, 1.0}, { 0.5, 0.5, 1.0}, { 0.0, 0.5, 1.0},
      { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}
    };
    
    
    switch(cell.getKey() ) {
      
      // Base line topologies
    case shards::Line<2>::key:
    case shards::ShellLine<2>::key:
    case shards::Beam<2>::key:
      return line[nodeOrd];
      break;
        
      // Extended line topologies
    case shards::Line<3>::key:
    case shards::ShellLine<3>::key:
    case shards::Beam<3>::key:
      return line_3[nodeOrd];
      break;
        
        
      // Base triangle topologies
    case shards::Triangle<3>::key:
    case shards::ShellTriangle<3>::key:
      return triangle[nodeOrd];
      break;
        
      // Extened Triangle topologies
    case shards::Triangle<4>::key:
      return triangle_4[nodeOrd];
      break;
    case shards::Triangle<6>::key:
    case shards::ShellTriangle<6>::key:
      return triangle_6[nodeOrd];
      break;
        
        
      // Base Quadrilateral topologies  
    case shards::Quadrilateral<4>::key:
    case shards::ShellQuadrilateral<4>::key:
      return quadrilateral[nodeOrd];
      break;
        
      // Extended Quadrilateral topologies
    case shards::Quadrilateral<8>::key:
    case shards::ShellQuadrilateral<8>::key:
      return quadrilateral_8[nodeOrd];
      break;
    case shards::Quadrilateral<9>::key:
    case shards::ShellQuadrilateral<9>::key:
      return quadrilateral_9[nodeOrd];
      break;
        
        
      // Base Tetrahedron topology
    case shards::Tetrahedron<4>::key:
      return tetrahedron[nodeOrd];
      break;
        
      // Extended Tetrahedron topologies
    case shards::Tetrahedron<8>::key:
      return tetrahedron_8[nodeOrd];
      break;
    case shards::Tetrahedron<10>::key:
      return tetrahedron_10[nodeOrd];
      break;
    case shards::Tetrahedron<11>::key:
      return tetrahedron_11[nodeOrd];
      break;

        
      // Base Hexahedron topology
    case shards::Hexahedron<8>::key:
      return hexahedron[nodeOrd];
      break;
        
      // Extended Hexahedron topologies
    case shards::Hexahedron<20>::key:
      return hexahedron_20[nodeOrd];
      break;
    case shards::Hexahedron<27>::key:
      return hexahedron_27[nodeOrd];
      break;

        
      // Base Pyramid topology  
    case shards::Pyramid<5>::key:
      return pyramid[nodeOrd];
      break;
        
      // Extended pyramid topologies
    case shards::Pyramid<13>::key:
      return pyramid_13[nodeOrd];
      break;
    case shards::Pyramid<14>::key:
      return pyramid_14[nodeOrd];
      break;
      
        
      // Base Wedge topology
    case shards::Wedge<6>::key:
      return wedge[nodeOrd];
      break;
        
      // Extended Wedge topologies
    case shards::Wedge<15>::key:
      return wedge_15[nodeOrd];
      break;
    case shards::Wedge<18>::key:
      return wedge_18[nodeOrd];
      break;
        
    default:
      INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::getReferenceNode): invalid cell topology.");
    }
    // To disable compiler warning, should never be reached
    return line[0];
  }
  
  
  
  template<class Scalar>
  template<class ArraySubcellNode>
  void CellTools<Scalar>::getReferenceSubcellNodes(ArraySubcellNode &          subcellNodes,
                                                   const int                   subcellDim,
                                                   const int                   subcellOrd,
                                                   const shards::CellTopology& parentCell){
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !(hasReferenceCell(parentCell) ), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): the specified cell topology does not have a reference cell.");
    
    // subcellDim can equal the cell dimension because the cell itself is a valid subcell! In this case
    // the method will return all cell cellWorkset.
    INTREPID2_TEST_FOR_EXCEPTION( !( (0 <= subcellDim) && (subcellDim <= (int)parentCell.getDimension()) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): subcell dimension out of range.");
    
    INTREPID2_TEST_FOR_EXCEPTION( !( (0 <= subcellOrd) && (subcellOrd < (int)parentCell.getSubcellCount(subcellDim) ) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes): subcell ordinal out of range.");
    
    // Verify subcellNodes rank and dimensions
    {
      std::string errmsg = ">>> ERROR (Intrepid2::CellTools::getReferenceSubcellNodes):";
      INTREPID2_TEST_FOR_EXCEPTION( !( requireRankRange(errmsg, subcellNodes, 2, 2) ), std::invalid_argument, errmsg);
      
      int subcNodeCount = parentCell.getNodeCount(subcellDim, subcellOrd);
      int spaceDim = parentCell.getDimension();
      
      INTREPID2_TEST_FOR_EXCEPTION( !( requireDimensionRange(errmsg, subcellNodes, 0,  subcNodeCount, subcNodeCount) ),
                                    std::invalid_argument, errmsg);
      
      INTREPID2_TEST_FOR_EXCEPTION( !( requireDimensionRange(errmsg, subcellNodes, 1,  spaceDim, spaceDim) ),
                                    std::invalid_argument, errmsg);
    }
#endif 
    
    // Find how many cellWorkset does the specified subcell have.
    int subcNodeCount = parentCell.getNodeCount(subcellDim, subcellOrd);
    
    // Loop over subcell cellWorkset
    for(int subcNodeOrd = 0; subcNodeOrd < subcNodeCount; subcNodeOrd++){
      
      // Get the node number relative to the parent reference cell
      int cellNodeOrd = parentCell.getNodeMap(subcellDim, subcellOrd, subcNodeOrd);
            
      // Loop over node's Cartesian coordinates
      for(int dim = 0; dim < (int)parentCell.getDimension(); dim++){
        subcellNodes(subcNodeOrd, dim) = CellTools::getReferenceNode(parentCell, cellNodeOrd)[dim];
      }
    }
  }  
   
  template<class Scalar>
  template<class ArrayEdgeTangent>
  void CellTools<Scalar>::getReferenceEdgeTangent(ArrayEdgeTangent &            refEdgeTangent,
                                                  const int &                   edgeOrd,
                                                  const shards::CellTopology &  parentCell){
  
    int spaceDim  = parentCell.getDimension();
  
#ifdef HAVE_INTREPID2_DEBUG
  
    INTREPID2_TEST_FOR_EXCEPTION( !( (spaceDim == 2) || (spaceDim == 3) ), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): two or three-dimensional parent cell required");
  
    INTREPID2_TEST_FOR_EXCEPTION( !( (0 <= edgeOrd) && (edgeOrd < (int)parentCell.getSubcellCount(1) ) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): edge ordinal out of bounds");  
  
    INTREPID2_TEST_FOR_EXCEPTION( !( refEdgeTangent.size() == spaceDim ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): output array size is required to match space dimension");  
#endif
    // Edge parametrizations are computed in setSubcellParametrization and stored in rank-3 array 
    // (subcOrd, coordinate, coefficient)
    const FieldContainer<double>& edgeMap = getSubcellParametrization(1, parentCell);
  
    // All ref. edge maps have affine coordinate functions: f_dim(u) = C_0(dim) + C_1(dim)*u, 
    //                                     => edge Tangent: -> C_1(*)
    refEdgeTangent(0) = edgeMap(edgeOrd, 0, 1);
    refEdgeTangent(1) = edgeMap(edgeOrd, 1, 1);
  
    // Skip last coordinate for 2D parent cells
    if(spaceDim == 3) {
      refEdgeTangent(2) = edgeMap(edgeOrd, 2, 1);  
    }
  }



  template<class Scalar>
  template<class ArrayFaceTangentU, class ArrayFaceTangentV>
  void CellTools<Scalar>::getReferenceFaceTangents(ArrayFaceTangentU &           uTan,
                                                   ArrayFaceTangentV &           vTan,
                                                   const int &                   faceOrd,
                                                   const shards::CellTopology &  parentCell){
  
#ifdef HAVE_INTREPID2_DEBUG
    int spaceDim  = parentCell.getDimension();
    INTREPID2_TEST_FOR_EXCEPTION( !(spaceDim == 3), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): three-dimensional parent cell required");  
  
    INTREPID2_TEST_FOR_EXCEPTION( !( (0 <= faceOrd) && (faceOrd < (int)parentCell.getSubcellCount(2) ) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): face ordinal out of bounds");  

    INTREPID2_TEST_FOR_EXCEPTION( !( (getrank(uTan) == 1)  && (getrank(vTan) == 1) ), std::invalid_argument,  
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): rank = 1 required for output arrays"); 
  
    INTREPID2_TEST_FOR_EXCEPTION( !( uTan.dimension(0) == spaceDim ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): dim0 (spatial dim) must match that of parent cell");  

    INTREPID2_TEST_FOR_EXCEPTION( !( vTan.dimension(0) == spaceDim ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceTangents): dim0 (spatial dim) must match that of parent cell");  
#endif
  
    // Face parametrizations are computed in setSubcellParametrization and stored in rank-3 array 
    // (subcOrd, coordinate, coefficient): retrieve this array
    const FieldContainer<double>& faceMap = getSubcellParametrization(2, parentCell);
  
    /*  All ref. face maps have affine coordinate functions:  f_dim(u,v) = C_0(dim) + C_1(dim)*u + C_2(dim)*v
     *                           `   => Tangent vectors are:  uTan -> C_1(*);    vTan -> C_2(*)
     */
    // set uTan -> C_1(*)
    uTan(0) = faceMap(faceOrd, 0, 1);
    uTan(1) = faceMap(faceOrd, 1, 1);
    uTan(2) = faceMap(faceOrd, 2, 1);
    
    // set vTan -> C_2(*)
    vTan(0) = faceMap(faceOrd, 0, 2);
    vTan(1) = faceMap(faceOrd, 1, 2);
    vTan(2) = faceMap(faceOrd, 2, 2);
  }



  template<class Scalar>
  template<class ArraySideNormal>
  void CellTools<Scalar>::getReferenceSideNormal(ArraySideNormal &             refSideNormal,
                                                 const int &                   sideOrd,
                                                 const shards::CellTopology &  parentCell){
    int spaceDim  = parentCell.getDimension();
 
#ifdef HAVE_INTREPID2_DEBUG
  
    INTREPID2_TEST_FOR_EXCEPTION( !( (spaceDim == 2) || (spaceDim == 3) ), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSideNormal): two or three-dimensional parent cell required");
  
    // Check side ordinal: by definition side is subcell whose dimension = spaceDim-1
    INTREPID2_TEST_FOR_EXCEPTION( !( (0 <= sideOrd) && (sideOrd < (int)parentCell.getSubcellCount(spaceDim - 1) ) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceSideNormal): side ordinal out of bounds");    
#endif 
    if(spaceDim == 2){
    
      // 2D parent cells: side = 1D subcell (edge), call the edge tangent method and rotate tangents
      getReferenceEdgeTangent(refSideNormal, sideOrd, parentCell);
    
      // rotate t(t1, t2) to get n(t2, -t1) so that (n,t) is positively oriented: det(n1,n2/t1,t2)>0
      Scalar temp = refSideNormal(0);
      refSideNormal(0) = refSideNormal(1);
      refSideNormal(1) = -temp;
    }
    else{
      // 3D parent cell: side = 2D subcell (face), call the face normal method.
      getReferenceFaceNormal(refSideNormal, sideOrd, parentCell);
    }
  }
  


  template<class Scalar>
  template<class ArrayFaceNormal>
  void CellTools<Scalar>::getReferenceFaceNormal(ArrayFaceNormal &             refFaceNormal,
                                                 const int &                   faceOrd,
                                                 const shards::CellTopology &  parentCell){
    int spaceDim  = parentCell.getDimension();
#ifdef HAVE_INTREPID2_DEBUG
  
    INTREPID2_TEST_FOR_EXCEPTION( !(spaceDim == 3), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceNormal): three-dimensional parent cell required");  
  
    INTREPID2_TEST_FOR_EXCEPTION( !( (0 <= faceOrd) && (faceOrd < (int)parentCell.getSubcellCount(2) ) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceNormal): face ordinal out of bounds");  
  
    INTREPID2_TEST_FOR_EXCEPTION( !( getrank(refFaceNormal) == 1 ), std::invalid_argument,  
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceNormal): rank = 1 required for output array"); 
    
    INTREPID2_TEST_FOR_EXCEPTION( !( static_cast<index_type>(refFaceNormal.dimension(0)) == static_cast<index_type>(spaceDim) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceFaceNormal): dim0 (spatial dim) must match that of parent cell");  
#endif

    // Reference face normal = vector product of reference face tangents. Allocate temp FC storage:
    FieldContainer<Scalar> uTan(spaceDim);
    FieldContainer<Scalar> vTan(spaceDim);
    getReferenceFaceTangents(uTan, vTan, faceOrd, parentCell);
  
    // Compute the vector product of the reference face tangents:
    RealSpaceTools<Scalar>::vecprod(refFaceNormal, uTan, vTan);
  }

  template<class Scalar>
  template<class ArrayEdgeTangent, class ArrayJac>
  void CellTools<Scalar>::getPhysicalEdgeTangents(ArrayEdgeTangent &            edgeTangents,
                                                  const ArrayJac &              worksetJacobians,
                                                  const int &                   worksetEdgeOrd,
                                                  const shards::CellTopology &  parentCell){
    index_type worksetSize = static_cast<index_type>(worksetJacobians.dimension(0));
    index_type edgePtCount = static_cast<index_type>(worksetJacobians.dimension(1)); 
    int pCellDim    = parentCell.getDimension();
#ifdef HAVE_INTREPID2_DEBUG
    std::string errmsg = ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents):";
  
    INTREPID2_TEST_FOR_EXCEPTION( !( (pCellDim == 3) || (pCellDim == 2) ), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalEdgeTangents): 2D or 3D parent cell required");  
  
    // (1) edgeTangents is rank-3 (C,P,D) and D=2, or 3 is required
    INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, edgeTangents, 3,3), std::invalid_argument, errmsg);
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, edgeTangents, 2, 2,3), std::invalid_argument, errmsg);
 
    // (2) worksetJacobians in rank-4 (C,P,D,D) and D=2, or 3 is required
    INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, worksetJacobians, 4,4), std::invalid_argument, errmsg);
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 2, 2,3), std::invalid_argument, errmsg);
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 3, 2,3), std::invalid_argument, errmsg);
  
    // (4) cross-check array dimensions: edgeTangents (C,P,D) vs. worksetJacobians (C,P,D,D)
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, edgeTangents, 0,1,2,2,  worksetJacobians, 0,1,2,3), std::invalid_argument, errmsg);      
  
#endif
  
    // Temp storage for constant reference edge tangent: rank-1 (D) arrays
    FieldContainer<double> refEdgeTan(pCellDim);
    getReferenceEdgeTangent(refEdgeTan, worksetEdgeOrd, parentCell);
  
    // Loop over workset faces and edge points
    for(index_type pCell = 0; pCell < worksetSize; pCell++){
      for(index_type pt = 0; pt < edgePtCount; pt++){
      
        // Apply parent cell Jacobian to ref. edge tangent
        for(int i = 0; i < pCellDim; i++){
          edgeTangents(pCell, pt, i) = 0.0;
          for(int j = 0; j < pCellDim; j++){
            edgeTangents(pCell, pt, i) +=  worksetJacobians(pCell, pt, i, j)*refEdgeTan(j);
          }// for j
        }// for i
      }// for pt
    }// for pCell
  }
  template<class Scalar>
  template<class ArrayFaceTangentU, class ArrayFaceTangentV, class ArrayJac>
  void CellTools<Scalar>::getPhysicalFaceTangents(ArrayFaceTangentU &           faceTanU,
                                                  ArrayFaceTangentV &           faceTanV,
                                                  const ArrayJac &              worksetJacobians,
                                                  const int &                   worksetFaceOrd,
                                                  const shards::CellTopology &  parentCell){
    index_type worksetSize = static_cast<index_type>(worksetJacobians.dimension(0));
    index_type facePtCount = static_cast<index_type>(worksetJacobians.dimension(1)); 
    int pCellDim    = parentCell.getDimension();
#ifdef HAVE_INTREPID2_DEBUG
    std::string errmsg = ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents):";

    INTREPID2_TEST_FOR_EXCEPTION( !(pCellDim == 3), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceTangents): three-dimensional parent cell required");  
  
    // (1) faceTanU and faceTanV are rank-3 (C,P,D) and D=3 is required
    INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, faceTanU, 3,3), std::invalid_argument, errmsg);
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, faceTanU, 2, 3,3), std::invalid_argument, errmsg);

    INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, faceTanV, 3,3), std::invalid_argument, errmsg);
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, faceTanV, 2, 3,3), std::invalid_argument, errmsg);

    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, faceTanU,  faceTanV), std::invalid_argument, errmsg);      

    // (3) worksetJacobians in rank-4 (C,P,D,D) and D=3 is required
    INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, worksetJacobians, 4,4), std::invalid_argument, errmsg);
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 2, 3,3), std::invalid_argument, errmsg);
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 3, 3,3), std::invalid_argument, errmsg);

    // (4) cross-check array dimensions: faceTanU (C,P,D) vs. worksetJacobians (C,P,D,D)
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, faceTanU, 0,1,2,2,  worksetJacobians, 0,1,2,3), std::invalid_argument, errmsg);      

#endif
    
    // Temp storage for the pair of constant ref. face tangents: rank-1 (D) arrays
    FieldContainer<double> refFaceTanU(pCellDim);
    FieldContainer<double> refFaceTanV(pCellDim);
    getReferenceFaceTangents(refFaceTanU, refFaceTanV, worksetFaceOrd, parentCell);

    // Loop over workset faces and face points
    for(index_type pCell = 0; pCell < worksetSize; pCell++){
      for(index_type pt = 0; pt < facePtCount; pt++){
      
        // Apply parent cell Jacobian to ref. face tangents
        for(int dim = 0; dim < pCellDim; dim++){
          faceTanU(pCell, pt, dim) = 0.0;
          faceTanV(pCell, pt, dim) = 0.0;
        
          // Unroll loops: parent cell dimension can only be 3
          faceTanU(pCell, pt, dim) = \
            worksetJacobians(pCell, pt, dim, 0)*refFaceTanU(0) + \
            worksetJacobians(pCell, pt, dim, 1)*refFaceTanU(1) + \
            worksetJacobians(pCell, pt, dim, 2)*refFaceTanU(2);
          faceTanV(pCell, pt, dim) = \
            worksetJacobians(pCell, pt, dim, 0)*refFaceTanV(0) + \
            worksetJacobians(pCell, pt, dim, 1)*refFaceTanV(1) + \
            worksetJacobians(pCell, pt, dim, 2)*refFaceTanV(2);
        }// for dim
      }// for pt
    }// for pCell
  }

  template<class Scalar>
  template<class ArraySideNormal, class ArrayJac>
  void CellTools<Scalar>::getPhysicalSideNormals(ArraySideNormal &             sideNormals,
                                                 const ArrayJac &              worksetJacobians,
                                                 const int &                   worksetSideOrd,
                                                 const shards::CellTopology &  parentCell){
    index_type worksetSize = static_cast<index_type>(worksetJacobians.dimension(0));
    index_type sidePtCount = static_cast<index_type>(worksetJacobians.dimension(1));   
    int spaceDim  = parentCell.getDimension();
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !( (spaceDim == 2) || (spaceDim == 3) ), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalSideNormals): two or three-dimensional parent cell required");
  
    // Check side ordinal: by definition side is subcell whose dimension = spaceDim-1
    INTREPID2_TEST_FOR_EXCEPTION( !( (0 <= worksetSideOrd) && (worksetSideOrd < (int)parentCell.getSubcellCount(spaceDim - 1) ) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalSideNormals): side ordinal out of bounds");  
#endif  
  
    if(spaceDim == 2){

      // 2D parent cells: side = 1D subcell (edge), call the edge tangent method and rotate tangents
      getPhysicalEdgeTangents(sideNormals, worksetJacobians, worksetSideOrd, parentCell);
    
      // rotate t(t1, t2) to get n(t2, -t1) so that (n,t) is positively oriented: det(n1,n2/t1,t2)>0
      for(index_type cell = 0; cell < worksetSize; cell++){
        for(index_type pt = 0; pt < sidePtCount; pt++){
          Scalar temp = sideNormals(cell, pt, 0);
          sideNormals(cell, pt, 0) = sideNormals(cell, pt, 1);
          sideNormals(cell, pt, 1) = -temp;
        }// for pt
      }// for cell
    }
    else{
      // 3D parent cell: side = 2D subcell (face), call the face normal method.
      getPhysicalFaceNormals(sideNormals, worksetJacobians, worksetSideOrd, parentCell);
    }
  }
  
  
  template<class Scalar>
  template<class ArrayFaceNormal, class ArrayJac>
  void CellTools<Scalar>::getPhysicalFaceNormals(ArrayFaceNormal &             faceNormals,
                                                 const ArrayJac &              worksetJacobians,
                                                 const int &                   worksetFaceOrd,
                                                 const shards::CellTopology &  parentCell){
    index_type worksetSize = static_cast<index_type>(worksetJacobians.dimension(0));
    index_type facePtCount = static_cast<index_type>(worksetJacobians.dimension(1)); 
    int pCellDim    = parentCell.getDimension();
#ifdef HAVE_INTREPID2_DEBUG
    std::string errmsg = ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals):";
  
    INTREPID2_TEST_FOR_EXCEPTION( !(pCellDim == 3), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::getPhysicalFaceNormals): three-dimensional parent cell required");  
  
    // (1) faceNormals is rank-3 (C,P,D) and D=3 is required
    INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, faceNormals, 3,3), std::invalid_argument, errmsg);
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, faceNormals, 2, 3,3), std::invalid_argument, errmsg);
  
    // (3) worksetJacobians in rank-4 (C,P,D,D) and D=3 is required
    INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, worksetJacobians, 4,4), std::invalid_argument, errmsg);
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 2, 3,3), std::invalid_argument, errmsg);
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, worksetJacobians, 3, 3,3), std::invalid_argument, errmsg);
  
    // (4) cross-check array dimensions: faceNormals (C,P,D) vs. worksetJacobians (C,P,D,D)
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, faceNormals, 0,1,2,2,  worksetJacobians, 0,1,2,3), std::invalid_argument, errmsg);        
#endif
  
    // Temp storage for physical face tangents: rank-3 (C,P,D) arrays
    FieldContainer<Scalar> faceTanU(worksetSize, facePtCount, pCellDim);
    FieldContainer<Scalar> faceTanV(worksetSize, facePtCount, pCellDim);
    getPhysicalFaceTangents(faceTanU, faceTanV, worksetJacobians, worksetFaceOrd, parentCell);
  
    // Compute the vector product of the physical face tangents:
    RealSpaceTools<Scalar>::vecprod(faceNormals, faceTanU, faceTanV);
  
  
  }

}

#endif
