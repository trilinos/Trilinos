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

/** \file   Intrepid2_CellTools_Serial.hpp
    \brief  Definition file for the Intrepid2::Impl::CellTools class.
    \author Kyungjoo Kim
*/

#ifndef __INTREPID2_CELLTOOLS_SERIAL_HPP__
#define __INTREPID2_CELLTOOLS_SERIAL_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#include "Shards_CellTopology.hpp"

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_Kernels.hpp"

namespace Intrepid2 {

  namespace Impl {
    
    /**
    \brief See Intrepid2::CellTools
    */
    class CellTools {
    public:
      typedef Kokkos::DynRankView<double,Kokkos::HostSpace> NodeDataHostView;
      typedef Kokkos::DynRankView<const double,Kokkos::HostSpace,Kokkos::MemoryUnmanaged> ConstUnmanagedNodeDataHostView;

      struct ReferenceNodeDataType {
        double line[2][3], line_3[3][3];
        double triangle[3][3], triangle_4[4][3], triangle_6[6][3];
        double quadrilateral[4][3], quadrilateral_8[8][3], quadrilateral_9[9][3];
        double tetrahedron[4][3], tetrahedron_8[8][3], tetrahedron_10[10][3], tetrahedron_11[10][3];
        double hexahedron[8][3], hexahedron_20[20][3], hexahedron_27[27][3];
        double pyramid[5][3], pyramid_13[13][3], pyramid_14[14][3];
        double wedge[6][3], wedge_15[15][3], wedge_18[18][3];
      };
      
      inline 
      static const ReferenceNodeDataType& 
      getRefNodeData() {
        const static ReferenceNodeDataType refNodeData = {
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
        return refNodeData;
      }
      
      template<typename refNodeViewType>
      static 
      void
      getReferenceNode(const refNodeViewType &nodes,
                       const shards::CellTopology  &cell,
                       const ordinal_type nodeOrd ) {
        ConstUnmanagedNodeDataHostView ref;
        switch (cell.getKey() ) {
        case shards::Line<2>::key:
        case shards::ShellLine<2>::key:
        case shards::Beam<2>::key:               ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().line, 2, 3); break;
        case shards::Line<3>::key:
        case shards::ShellLine<3>::key:
        case shards::Beam<3>::key:               ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().line_3, 3, 3); break;

        case shards::Triangle<3>::key:
        case shards::ShellTriangle<3>::key:      ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().triangle, 3, 3); break; 
        case shards::Triangle<4>::key:           ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().triangle_4, 4, 3); break; 
        case shards::Triangle<6>::key:
        case shards::ShellTriangle<6>::key:      ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().triangle_6, 6, 3); break; 

        case shards::Quadrilateral<4>::key:
        case shards::ShellQuadrilateral<4>::key: ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().quadrilateral, 4, 3); break; 
        case shards::Quadrilateral<8>::key:
        case shards::ShellQuadrilateral<8>::key: ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().quadrilateral_8, 8, 3); break; 
        case shards::Quadrilateral<9>::key:
        case shards::ShellQuadrilateral<9>::key: ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().quadrilateral_9, 9, 3); break; 

        case shards::Tetrahedron<4>::key:        ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().tetrahedron, 4, 3); break; 
        case shards::Tetrahedron<8>::key:        ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().tetrahedron_8, 8, 3); break; 
        case shards::Tetrahedron<10>::key:       ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().tetrahedron_10, 10, 3); break; 
        case shards::Tetrahedron<11>::key:       ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().tetrahedron_11, 11, 3); break; 

        case shards::Hexahedron<8>::key:         ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().hexahedron, 8, 3); break; 
        case shards::Hexahedron<20>::key:        ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().hexahedron_20, 20, 3); break; 
        case shards::Hexahedron<27>::key:        ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().hexahedron_27, 27, 3); break; 

        case shards::Pyramid<5>::key:            ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().pyramid, 5, 3); break; 
        case shards::Pyramid<13>::key:           ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().pyramid_13, 13, 3); break; 
        case shards::Pyramid<14>::key:           ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().pyramid_14, 14, 3); break; 

        case shards::Wedge<6>::key:              ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().wedge, 6, 3); break; 
        case shards::Wedge<15>::key:             ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().wedge_15, 15, 3); break; 
        case shards::Wedge<18>::key:             ref = ConstUnmanagedNodeDataHostView((const double*)getRefNodeData().wedge_18, 18, 3); break; 

        default: {
          INTREPID2_TEST_FOR_ABORT( true, "invalid input cell topology.");
          break;
        }
        }

        const ordinal_type D = cell.getDimension();
        for (ordinal_type i=0;i<D;++i)
          nodes(i) = ref(nodeOrd, i);
      }

      struct SubcellParamDataType {
        NodeDataHostView dummy;
        NodeDataHostView lineEdges;  // edge maps for 2d non-standard cells; shell line and beam
        NodeDataHostView triEdges, quadEdges; // edge maps for 2d standard cells
        NodeDataHostView shellTriEdges, shellQuadEdges; // edge maps for 3d non-standard cells; shell tri and quad
        NodeDataHostView tetEdges, hexEdges, pyrEdges, wedgeEdges; // edge maps for 3d standard cells
        NodeDataHostView shellTriFaces, shellQuadFaces; // face maps for 3d non-standard cells
        NodeDataHostView tetFaces, hexFaces, pyrFaces, wedgeFaces; // face maps for 3d standard cells
      };

      inline
      static SubcellParamDataType& getSubcellParamData() { 
        static SubcellParamDataType subcellParamData; 
        Kokkos::push_finalize_hook( [=] {
            subcellParamData.dummy = NodeDataHostView();
            subcellParamData.lineEdges = NodeDataHostView();
            subcellParamData.triEdges = NodeDataHostView();
            subcellParamData.quadEdges = NodeDataHostView();
            subcellParamData.shellTriEdges = NodeDataHostView();
            subcellParamData.shellQuadEdges = NodeDataHostView();
            subcellParamData.tetEdges = NodeDataHostView();
            subcellParamData.hexEdges = NodeDataHostView();
            subcellParamData.pyrEdges = NodeDataHostView();
            subcellParamData.wedgeEdges = NodeDataHostView();
            subcellParamData.shellTriFaces = NodeDataHostView();
            subcellParamData.shellQuadFaces = NodeDataHostView();
            subcellParamData.tetFaces = NodeDataHostView();
            subcellParamData.hexFaces = NodeDataHostView();
            subcellParamData.pyrFaces = NodeDataHostView();
            subcellParamData.wedgeFaces = NodeDataHostView();
          });
        return subcellParamData;
      }

      inline
      static void
      setSubcellParametrization() {
        static bool isSubcellParametrizationSet = false;
        if (!isSubcellParametrizationSet) {
          {
            const auto tet = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >());
            setSubcellParametrization( getSubcellParamData().tetFaces,   2, tet );
            setSubcellParametrization( getSubcellParamData().tetEdges,   1, tet );
          }
          {
            const auto hex = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >());
            setSubcellParametrization( getSubcellParamData().hexFaces,   2, hex );
            setSubcellParametrization( getSubcellParamData().hexEdges,   1, hex );
          }
          {
            const auto pyr = shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<5> >());
            setSubcellParametrization( getSubcellParamData().pyrFaces,   2, pyr );
            setSubcellParametrization( getSubcellParamData().pyrEdges,   1, pyr );
          }
          {
            const auto wedge = shards::CellTopology(shards::getCellTopologyData<shards::Wedge<6> >());
            setSubcellParametrization( getSubcellParamData().wedgeFaces, 2, wedge );
            setSubcellParametrization( getSubcellParamData().wedgeEdges, 1, wedge );
          }
          {
            const auto tri = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >());
            setSubcellParametrization( getSubcellParamData().triEdges,   1, tri );
          }
          {
            const auto quad = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >());
            setSubcellParametrization( getSubcellParamData().quadEdges,  1, quad );
          }
          {
            const auto line = shards::CellTopology(shards::getCellTopologyData<shards::ShellLine<2> >());
            setSubcellParametrization( getSubcellParamData().lineEdges,  1, line );
          }
        }
        isSubcellParametrizationSet = true;
      }

      inline
      static void
      setSubcellParametrization(      NodeDataHostView &subcellParam,
                                const ordinal_type subcellDim,
                                const shards::CellTopology &parentCell ) {
        // get subcellParametrization dimensions: (sc, pcd, coeff)
        const auto sc    = parentCell.getSubcellCount(subcellDim);
        const auto pcd   = parentCell.getDimension();
        const auto coeff = (subcellDim == 1) ? 2 : 3;
        
        // create a view
        subcellParam = NodeDataHostView("subcellParam",
                                        sc, pcd, coeff);
       
        NodeDataHostView v0("v0", 3), v1("v1", 3), v2("v1", 3), v3("v1", 3);
        if (subcellDim == 1) {
          // Edge parametrizations of 2D and 3D cells (shell lines and beams are 2D cells with edges)
          for (size_type subcellOrd=0;subcellOrd<sc;++subcellOrd) {
            // vertexK[0] = x_k; vertexK[1] = y_k; vertexK[2] = z_k; z_k = 0 for 2D cells
            // Note that ShellLine and Beam are 2D cells!
            const auto v0ord = parentCell.getNodeMap(subcellDim, subcellOrd, 0);
            const auto v1ord = parentCell.getNodeMap(subcellDim, subcellOrd, 1);
            
            getReferenceNode(v0, parentCell, v0ord);
            getReferenceNode(v1, parentCell, v1ord);
            
            // x(t) = (x0 + x1)/2 + t*(x1 - x0)/2
            subcellParam(subcellOrd, 0, 0) = (v0[0] + v1[0])/2.0;
            subcellParam(subcellOrd, 0, 1) = (v1[0] - v0[0])/2.0;
            
            // y(t) = (y0 + y1)/2 + t*(y1 - y0)/2
            subcellParam(subcellOrd, 1, 0) = (v0[1] + v1[1])/2.0;
            subcellParam(subcellOrd, 1, 1) = (v1[1] - v0[1])/2.0;
            
            if( pcd == 3 ) {
              // z(t) = (z0 + z1)/2 + t*(z1 - z0)/2
              subcellParam(subcellOrd, 2, 0) = (v0[2] + v1[2])/2.0;
              subcellParam(subcellOrd, 2, 1) = (v1[2] - v0[2])/2.0;
            }
          }
        }
        else if (subcellDim == 2) {
          // Face parametrizations of 3D cells: (shell Tri and Quad are 3D cells with faces)
          // A 3D cell can have both Tri and Quad faces, but because they are affine images of the
          // parametrization domain, 3 coefficients are enough to store them in both cases.
          for (size_type subcellOrd=0;subcellOrd<sc;++subcellOrd) {
            
            switch (parentCell.getKey(subcellDim,subcellOrd)) {
              
            case shards::Triangle<3>::key:
            case shards::Triangle<4>::key:
            case shards::Triangle<6>::key: {
              const auto v0ord = parentCell.getNodeMap(subcellDim, subcellOrd, 0);
              const auto v1ord = parentCell.getNodeMap(subcellDim, subcellOrd, 1);
              const auto v2ord = parentCell.getNodeMap(subcellDim, subcellOrd, 2);
              
              getReferenceNode(v0, parentCell, v0ord);
              getReferenceNode(v1, parentCell, v1ord);
              getReferenceNode(v2, parentCell, v2ord);
              
              // x(u,v) = x0 + (x1 - x0)*u + (x2 - x0)*v
              subcellParam(subcellOrd, 0, 0) = v0[0];
              subcellParam(subcellOrd, 0, 1) = v1[0] - v0[0];
              subcellParam(subcellOrd, 0, 2) = v2[0] - v0[0];
              
              // y(u,v) = y0 + (y1 - y0)*u + (y2 - y0)*v
              subcellParam(subcellOrd, 1, 0) = v0[1];
              subcellParam(subcellOrd, 1, 1) = v1[1] - v0[1];
              subcellParam(subcellOrd, 1, 2) = v2[1] - v0[1];
              
              // z(u,v) = z0 + (z1 - z0)*u + (z2 - z0)*v
              subcellParam(subcellOrd, 2, 0) = v0[2];
              subcellParam(subcellOrd, 2, 1) = v1[2] - v0[2];
              subcellParam(subcellOrd, 2, 2) = v2[2] - v0[2];
              break;
            }
            case shards::Quadrilateral<4>::key:
            case shards::Quadrilateral<8>::key:
            case shards::Quadrilateral<9>::key: {
              const auto v0ord = parentCell.getNodeMap(subcellDim, subcellOrd, 0);
              const auto v1ord = parentCell.getNodeMap(subcellDim, subcellOrd, 1);
              const auto v2ord = parentCell.getNodeMap(subcellDim, subcellOrd, 2);
              const auto v3ord = parentCell.getNodeMap(subcellDim, subcellOrd, 3);
              
              getReferenceNode(v0, parentCell, v0ord);
              getReferenceNode(v1, parentCell, v1ord);
              getReferenceNode(v2, parentCell, v2ord);
              getReferenceNode(v3, parentCell, v3ord);
              
              // x(u,v) = (x0+x1+x2+x3)/4+u*(-x0+x1+x2-x3)/4+v*(-x0-x1+x2+x3)/4+uv*(0=x0-x1+x2-x3)/4
              subcellParam(subcellOrd, 0, 0) = ( v0[0] + v1[0] + v2[0] + v3[0])/4.0;
              subcellParam(subcellOrd, 0, 1) = (-v0[0] + v1[0] + v2[0] - v3[0])/4.0;
              subcellParam(subcellOrd, 0, 2) = (-v0[0] - v1[0] + v2[0] + v3[0])/4.0;
              // y(u,v) = (y0+y1+y2+y3)/4+u*(-y0+y1+y2-y3)/4+v*(-y0-y1+y2+y3)/4+uv*(0=y0-y1+y2-y3)/4
              subcellParam(subcellOrd, 1, 0) = ( v0[1] + v1[1] + v2[1] + v3[1])/4.0;
              subcellParam(subcellOrd, 1, 1) = (-v0[1] + v1[1] + v2[1] - v3[1])/4.0;
              subcellParam(subcellOrd, 1, 2) = (-v0[1] - v1[1] + v2[1] + v3[1])/4.0;
              
              // z(u,v) = (z0+z1+z2+z3)/4+u*(-z0+z1+z2-z3)/4+v*(-z0-z1+z2+z3)/4+uv*(0=z0-z1+z2-z3)/4
              subcellParam(subcellOrd, 2, 0) = ( v0[2] + v1[2] + v2[2] + v3[2])/4.0;
              subcellParam(subcellOrd, 2, 1) = (-v0[2] + v1[2] + v2[2] - v3[2])/4.0;
              subcellParam(subcellOrd, 2, 2) = (-v0[2] - v1[2] + v2[2] + v3[2])/4.0;
              break;
            }
            default: {
              INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                            ">>> ERROR (Intrepid2::CellTools::setSubcellParametrization): parametrization not defined for the specified face topology.");
            }
            }
          }
        }
      }

      struct Serial {

        // output: 
        //   jacobian (D,sD) - jacobian matrix evaluated at a single point 
        // input: 
        //   grads    (N,sD) - hgrad basis grad values evaluated at a single point (C1/C2 element only)
        //   nodes    (N,D) - cell element-to-node connectivity
        template<typename jacobianViewType,
                 typename basisGradViewType,
                 typename nodeViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        computeJacobian(const jacobianViewType  &jacobian, // D,sD
                        const basisGradViewType &grads,    // N,sD
                        const nodeViewType      &nodes) {  // N,D
          const auto N = nodes.extent(0);

          const auto  D = jacobian.extent(0);
          const auto sD = jacobian.extent(1);
          
          INTREPID2_TEST_FOR_ABORT( N != grads.extent(0), "grad dimension_0 does not match to cardinality.");
          INTREPID2_TEST_FOR_ABORT(sD != grads.extent(1), "grad dimension_1 does not match to space dim.");
          INTREPID2_TEST_FOR_ABORT( D != nodes.extent(1), "node dimension_1 does not match to space dim.");

          Kernels::Serial::gemm_trans_notrans(1.0, nodes, grads, 0.0, jacobian);
        }

        // output: 
        //   point (D)   - mapped physical point 
        // input: 
        //   vals  (N)   - hgrad basis values evaluated at a single point (C1/C2 element only)
        //   nodes (N,D) - cell element-to-node connectivity
        template<typename pointViewType,
                 typename basisValViewType,
                 typename nodeViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        mapToPhysicalFrame(const pointViewType    &point,    // D  
                           const basisValViewType &vals,     // N  
                           const nodeViewType     &nodes) {  // N,D 
          const auto N = vals.extent(0);
          const auto D = point.extent(0);

          INTREPID2_TEST_FOR_ABORT(N != nodes.extent(0), "nodes dimension_0 does not match to vals dimension_0.");
          INTREPID2_TEST_FOR_ABORT(D != nodes.extent(1), "node dimension_1 does not match to space dim.");

          Kernels::Serial::gemv_trans(1.0, nodes, vals, 0.0, point);
        }

        // template:
        //   implBasisType - impl basis function type e.g., Impl::Basis_HGRAD_QUAD_C1_FEM
        // output: 
        //   xref (sD)   - point mapped to reference frame (subcell Dim)
        // input: 
        //   xphy  (D)   - point in physical frame
        //   nodes (N,D) - cell element-to-node connectivity
        template<typename implBasisType,
                 typename refPointViewType,
                 typename phyPointViewType,
                 typename nodeViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        mapToReferenceFrame(const refPointViewType &xref, // sD 
                            const phyPointViewType &xphy, // D
                            const nodeViewType &nodes) {  // N,D
          const ordinal_type sD = xref.extent(0);
          const ordinal_type D = xphy.extent(0);
          const ordinal_type N = nodes.extent(0);

          INTREPID2_TEST_FOR_ABORT(sD > D, "subcell dimension is greater than physical cell dimension.");
          INTREPID2_TEST_FOR_ABORT(D != static_cast<ordinal_type>(nodes.extent(1)), "xphy dimension_0 does not match to space dim.");
          
          typedef typename refPointViewType::non_const_value_type value_type;
          
          // I want to use view instead of dynrankview
          // NMAX = 28, MAXDIM = 3
          value_type buf[27*3 + 27 + 9 + 9 + 9 + 9 + 3 + 3] = {}, *ptr = &buf[0];
          Kokkos::DynRankView<value_type,
            Kokkos::Impl::ActiveExecutionMemorySpace,
            Kokkos::MemoryUnmanaged> grads(ptr, N, sD); ptr += N*sD;
          
          Kokkos::DynRankView<value_type,
            Kokkos::Impl::ActiveExecutionMemorySpace,
            Kokkos::MemoryUnmanaged> vals(ptr, N); ptr += N;

          Kokkos::DynRankView<value_type,
            Kokkos::Impl::ActiveExecutionMemorySpace,
            Kokkos::MemoryUnmanaged> jac(ptr, D, sD); ptr += D*sD; 

          Kokkos::DynRankView<value_type,
            Kokkos::Impl::ActiveExecutionMemorySpace,
            Kokkos::MemoryUnmanaged> metric(ptr, sD, sD); ptr += sD*sD;

          Kokkos::DynRankView<value_type,
            Kokkos::Impl::ActiveExecutionMemorySpace,
            Kokkos::MemoryUnmanaged> invmetric(ptr, sD, sD); ptr += sD*sD;

          Kokkos::DynRankView<value_type,
            Kokkos::Impl::ActiveExecutionMemorySpace,
            Kokkos::MemoryUnmanaged> invdf(ptr, sD, D); ptr += sD*D;

          Kokkos::DynRankView<value_type,
            Kokkos::Impl::ActiveExecutionMemorySpace,
            Kokkos::MemoryUnmanaged> xtmp(ptr, sD); ptr += sD;

          Kokkos::DynRankView<value_type,
            Kokkos::Impl::ActiveExecutionMemorySpace,
            Kokkos::MemoryUnmanaged> xold(ptr, sD); ptr += sD;
   
          // set initial guess
          for (ordinal_type j=0;j<D;++j) xold(j) = 0;
       
          const double tol = tolerence();
          for (ordinal_type iter=0;iter<Parameters::MaxNewton;++iter) {
            // xtmp := F(xold);
            implBasisType::template Serial<OPERATOR_VALUE>::getValues(vals, xold);
            mapToPhysicalFrame(xtmp, vals, nodes);

            // DF^{-1}
            implBasisType::template Serial<OPERATOR_GRAD>::getValues(grads, xold);
            CellTools::Serial::computeJacobian(jac, grads, nodes);

            Kernels::Serial::gemm_trans_notrans(1.0, jac, jac, 0.0, metric);
            Kernels::Serial::inverse(invmetric, metric);
            Kernels::Serial::gemm_notrans_trans(1.0, invmetric, jac, 0.0, invdf);

            // Newton
            Kernels::Serial::z_is_axby(xtmp, 1.0, xphy, -1.0, xtmp);  // xtmp := xphy - F(xold);
            Kernels::Serial::gemv_notrans(1.0, invdf, xtmp, 0.0, xref); // xref := DF^{-1}( xphy - F(xold))
            Kernels::Serial::z_is_axby(xref, 1.0, xold,  1.0, xref); // xref += xold
            
            // l2 error
            Kernels::Serial::z_is_axby(xtmp, 1.0, xold, -1.0, xref);

            double err = Kernels::Serial::norm(xtmp, NORM_ONE);

            if (err < tol) 
              break;

            Kernels::Serial::copy(xold, xref);
          }
        }
        
        ///
        /// Function for host only (it requires shards information)
        ///
        
        inline
        static ConstUnmanagedNodeDataHostView
        getSubcellParametrization(const ordinal_type subcell_dim,
                                  const shards::CellTopology parent_cell) {
          // all serial interface assumes that they can be called inside parallel for
          // lazy initilization is not a good idea; init is necessary 
          // setSubcellParametrization();
          Kokkos::DynRankView<const double,Kokkos::HostSpace,Kokkos::MemoryUnmanaged> r_val;          
          if (subcell_dim == 2) {
            switch (parent_cell.getBaseKey()) {
            case shards::Tetrahedron<>::key: r_val = getSubcellParamData().tetFaces; break;
            case shards::Hexahedron<>::key:  r_val = getSubcellParamData().hexFaces; break;
            case shards::Pyramid<>::key:     r_val = getSubcellParamData().pyrFaces; break;
            case shards::Wedge<18>::key:     r_val = getSubcellParamData().wedgeFaces; break;
            }
          } 
          else if (subcell_dim == 1) {
            switch (parent_cell.getBaseKey()) {
            case shards::Tetrahedron<>::key:   r_val = getSubcellParamData().tetEdges; break;
            case shards::Hexahedron<>::key:    r_val = getSubcellParamData().hexEdges; break;
            case shards::Pyramid<>::key:       r_val = getSubcellParamData().pyrEdges; break;
            case shards::Wedge<>::key:         r_val = getSubcellParamData().wedgeEdges; break;
            case shards::Triangle<>::key:      r_val = getSubcellParamData().triEdges; break;
            case shards::Quadrilateral<>::key: r_val = getSubcellParamData().quadEdges; break;
            case shards::Line<>::key:          r_val = getSubcellParamData().lineEdges; break;
            }
          }
          INTREPID2_TEST_FOR_ABORT(r_val.rank() == 0, "subcell param is not set up before.");          
          return r_val;
        }                 

        template<typename refEdgeTanViewType>
        inline
        static void
        getReferenceEdgeTangent(const refEdgeTanViewType &ref_edge_tangent,
                                const ordinal_type edge_ordinal,
                                const shards::CellTopology &parent_cell ) {
          const auto edge_map = getSubcellParametrization(1, parent_cell);

          const ordinal_type D = parent_cell.getDimension();
          for (ordinal_type i=0;i<D;++i)
            ref_edge_tangent(i) = edge_map(edge_ordinal, i, 1);
        }

        template<typename refFaceTanViewType>
        static void
        getReferenceFaceTangent(const refFaceTanViewType &ref_face_tan_u,
                                const refFaceTanViewType &ref_face_tan_v,
                                const ordinal_type face_ordinal,
                                const shards::CellTopology &parent_cell) {
          const auto face_map = getSubcellParametrization(2, parent_cell);
          
          // set refFaceTanU -> C_1(*)
          // set refFaceTanV -> C_2(*)
          const ordinal_type D = parent_cell.getDimension();
          for (ordinal_type i=0;i<D;++i) {
            ref_face_tan_u(i) = face_map(face_ordinal, i, 1);
            ref_face_tan_v(i) = face_map(face_ordinal, i, 2);
          }
        }

        template<typename edgeTangentViewType,
                 typename jacobianViewType>
        inline
        static void
        getPhysicalEdgeTangent(const edgeTangentViewType &edge_tangent, // D
                               const jacobianViewType &jacobian,         // D, D
                               const ordinal_type edge_ordinal,
                               const shards::CellTopology &parent_cell) {
          typedef typename edgeTangentViewType::non_const_value_type value_type;
          const ordinal_type D = parent_cell.getDimension();
          value_type buf[3];
          Kokkos::DynRankView<value_type,
            Kokkos::Impl::ActiveExecutionMemorySpace,
            Kokkos::MemoryUnmanaged> ref_edge_tangent(&buf[0], D); 

          getReferenceEdgeTangent(ref_edge_tangent, edge_ordinal, parent_cell);
          Kernels::Serial::matvec_product(edge_tangent, jacobian, ref_edge_tangent);
        }

        template<typename faceTanViewType,
                 typename jacobianViewType>
        inline
        static void
        getPhysicalFaceTangent(const faceTanViewType &face_tan_u, // D
                               const faceTanViewType &face_tan_v, // D
                               const jacobianViewType &jacobian,       // D, D
                               const ordinal_type face_ordinal,
                               const shards::CellTopology &parent_cell) {
          typedef typename faceTanViewType::non_const_value_type value_type;
          const ordinal_type D = parent_cell.getDimension();
          INTREPID2_TEST_FOR_ABORT(D != 3, "computing face normal requires dimension 3.");
          value_type buf[6];
          Kokkos::DynRankView<value_type,
            Kokkos::Impl::ActiveExecutionMemorySpace,
            Kokkos::MemoryUnmanaged> ref_face_tan_u(&buf[0], D), ref_face_tan_v(&buf[3], D);

          getReferenceFaceTangent(ref_face_tan_u, 
                                  ref_face_tan_v,
                                  face_ordinal,
                                  parent_cell);

          Kernels::Serial::matvec_product_d3(face_tan_u, jacobian, ref_face_tan_u);
          Kernels::Serial::matvec_product_d3(face_tan_v, jacobian, ref_face_tan_v);
        }


        template<typename faceNormalViewType,
                 typename jacobianViewType>
        inline
        static void
        getPhysicalFaceNormal(const faceNormalViewType &face_normal, // D
                              const jacobianViewType &jacobian,       // D, D
                              const ordinal_type face_ordinal,
                              const shards::CellTopology &parent_cell) {
          typedef typename faceNormalViewType::non_const_value_type value_type;
          const ordinal_type D = parent_cell.getDimension();
          INTREPID2_TEST_FOR_ABORT(D != 3, "computing face normal requires dimension 3.");
          value_type buf[6];
          Kokkos::DynRankView<value_type,
            Kokkos::Impl::ActiveExecutionMemorySpace,
            Kokkos::MemoryUnmanaged> face_tan_u(&buf[0], D), face_tan_v(&buf[3], D);

          getPhysicalFaceTangent(face_tan_u, face_tan_v,
                                 jacobian,
                                 face_ordinal,
                                 parent_cell);
          Kernels::Serial::vector_product_d3(face_normal, face_tan_u, face_tan_v);
        }

        template<typename sideNormalViewType,
                 typename jacobianViewType>
        inline
        static void
        getPhysicalSideNormal(const sideNormalViewType &side_normal, // D
                              const jacobianViewType &jacobian,       // D, D
                              const ordinal_type side_ordinal,     
                              const shards::CellTopology &parent_cell ) {
          const ordinal_type D = parent_cell.getDimension();
          typedef typename sideNormalViewType::non_const_value_type value_type;
          switch (D) {
          case 2: {
            value_type buf[3];
            Kokkos::DynRankView<value_type,
                Kokkos::Impl::ActiveExecutionMemorySpace,
                Kokkos::MemoryUnmanaged> edge_tangent(&buf[0], D); 
            getPhysicalEdgeTangent(edge_tangent, jacobian, side_ordinal, parent_cell);
            side_normal(0) =  edge_tangent(1);
            side_normal(1) = -edge_tangent(0);
            break;
          }
          case 3: {
            getPhysicalFaceNormal(side_normal, jacobian, side_ordinal, parent_cell);
            break;
          }
          default: {
            INTREPID2_TEST_FOR_ABORT(true, "cell dimension is out of range.");
            break;
          }
          }
        }
      };
    };
  }
}

#endif

