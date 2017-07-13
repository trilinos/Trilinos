#ifndef INTREPID2_TOYMESH_HPP
#define INTREPID2_TOYMESH_HPP

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

/** \file   Intrepid_ToyMesh.hpp
    \brief  For tests and examples, this toy mesh provides mesh connectivity 
    of high order elements.
    \author Created by Kyungjoo Kim
*/

#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )

namespace Intrepid2 {

  namespace Example {

    constexpr ordinal_type TOYMESH_MAX_NUM_NODES = 100;
    constexpr ordinal_type TOYMESH_NODE_STRIDE = 8;

    class ToyMesh {
    private:
      ordinal_type _nodes   [TOYMESH_MAX_NUM_NODES][TOYMESH_NODE_STRIDE];
      ordinal_type _nids    [TOYMESH_MAX_NUM_NODES];
      ordinal_type _offset  [TOYMESH_MAX_NUM_NODES];
      ordinal_type _boundary[TOYMESH_MAX_NUM_NODES];
      ordinal_type _i;

      void clearCurrentNode() {
        for (ordinal_type j=0;j<TOYMESH_NODE_STRIDE;++j)
          _nodes[_i][j] = -1;
      }
      void sortSubCellVertex(ordinal_type *sorted_verts,
                             const ordinal_type *subcell_verts,
                             const ordinal_type nids) {
        for (ordinal_type j=0;j<nids;++j)
          sorted_verts[j] = subcell_verts[j];
        std::sort(sorted_verts, sorted_verts+nids);
      }
      void addNode(const ordinal_type *subcell_verts, const ordinal_type nids,
                   const ordinal_type offset) {
        clearCurrentNode();
        ordinal_type sorted_verts[TOYMESH_NODE_STRIDE] = {};
        sortSubCellVertex(sorted_verts, subcell_verts, nids);
        for (ordinal_type j=0;j<nids;++j)
          _nodes[_i][j] = sorted_verts[j];
        _offset[_i] = offset;
        _nids[_i] = nids;
        ++_boundary[_i]; // 1 - boundary and interior, 2 - interface

        ++_i;
      }
      bool findNode(ordinal_type &offset, 
                    const ordinal_type *subcell_verts, 
                    const ordinal_type nids, 
                    const bool increase_boundary_flag) {
        ordinal_type sorted_verts[TOYMESH_NODE_STRIDE] = {};
        sortSubCellVertex(sorted_verts, subcell_verts, nids);
        for (ordinal_type i=0;i<_i;++i) {
          ordinal_type diff = std::abs(_nids[i] - nids);
          for (ordinal_type j=0;j<nids;++j)
            diff += std::abs(sorted_verts[j] - _nodes[i][j]);
          if (!diff) {
            offset = _offset[i];
            _boundary[i] += increase_boundary_flag;
            return true;
          }
        }
        return false;
      }
      bool findBoundary(ordinal_type &boundary, const ordinal_type *subcell_verts, const ordinal_type nids) {
        ordinal_type sorted_verts[TOYMESH_NODE_STRIDE] = {};
        sortSubCellVertex(sorted_verts, subcell_verts, nids);
        for (ordinal_type i=0;i<_i;++i) {
          ordinal_type diff = std::abs(_nids[i] - nids);
          for (ordinal_type j=0;j<nids;++j)
            diff += std::abs(sorted_verts[j] - _nodes[i][j]);
          if (!diff) {
            // 1 - boundary or interior (visited once) 2 - interface visited more than twice
            boundary = (_boundary[i] == 1); 
            return true;
          }
        }
        return false;
      }

    public:
      ordinal_type getCurrentNumNodes() const {
        return _i;
      }

      void reset() {
        for (ordinal_type i=0;i<TOYMESH_MAX_NUM_NODES;++i) {
          _nids[i]     = 0;
          _offset[i]   = 0;
          _boundary[i] = 0;
          for (ordinal_type j=0;j<TOYMESH_NODE_STRIDE;++j)
            _nodes[i][j] = -1;
        }
        _i = 0;
      }

      template<typename Scalar,typename ArrayType>
      void getLocalToGlobalMap(ordinal_type (*local2global)[2],
                               ordinal_type &off_global,
                               const Basis<Scalar,ArrayType> &basis,
                               const ordinal_type *element) {
        const ordinal_type local = 0, global = 1;
        const ordinal_type nbf = basis.getCardinality();
        const shards::CellTopology cell = basis.getBaseCellTopology();
        const ordinal_type dim = cell.getDimension();

        ordinal_type cnt = 0, off_element = 0;
        ordinal_type subcell_verts[4], nids;

        const ordinal_type nvert = cell.getVertexCount();
        for (ordinal_type i=0;i<nvert;++i) {
          const ordinal_type ord_vert = (off_element < nbf ? basis.getDofOrdinal(0, i, 0) : 0);
          const ordinal_type dof_vert = (off_element < nbf ? basis.getDofTag(ord_vert)[3] : 0);
      
          local2global[cnt][local] = off_element;
          off_element += dof_vert;
          Orientation::getElementNodeMap(subcell_verts, nids,
                                         cell, element,
                                         0, i);
      
          if (!findNode(local2global[cnt][global], subcell_verts, nids, true)) {
            addNode(subcell_verts, nids, off_global);
            local2global[cnt][global] = off_global;
            off_global += dof_vert;
          }
          ++cnt;
        }
        const ordinal_type nedge = cell.getEdgeCount();
        for (ordinal_type i=0;i<nedge;++i) {
          const ordinal_type ord_edge = (off_element < nbf ? basis.getDofOrdinal(1, i, 0) : 0);
          const ordinal_type dof_edge = (off_element < nbf ? basis.getDofTag(ord_edge)[3] : 0);
      
          local2global[cnt][local] = off_element;
          off_element += dof_edge;
          Orientation::getElementNodeMap(subcell_verts, nids,
                                         cell, element,
                                         1, i);
      
          if (!findNode(local2global[cnt][global], subcell_verts, nids, true)) {
            addNode(subcell_verts, nids, off_global);
            local2global[cnt][global] = off_global;
            off_global += dof_edge;
          }
          ++cnt;
        }
        const ordinal_type nface = cell.getFaceCount();
        for (ordinal_type i=0;i<nface;++i) {
          const ordinal_type ord_face = (off_element < nbf ? basis.getDofOrdinal(2, i, 0) : 0);
          const ordinal_type dof_face = (off_element < nbf ? basis.getDofTag(ord_face)[3] : 0);
      
          local2global[cnt][local] = off_element;
          off_element += dof_face;
          Orientation::getElementNodeMap(subcell_verts, nids,
                                         cell, element,
                                         2, i);
      
          if (!findNode(local2global[cnt][global], subcell_verts, nids, true)) {
            addNode(subcell_verts, nids, off_global);
            local2global[cnt][global] = off_global;
            off_global += dof_face;
          }
          ++cnt;
        }
        {
          const ordinal_type i = 0;
          const ordinal_type ord_intr = (off_element < nbf ? basis.getDofOrdinal(dim, i, 0) : 0);
          const ordinal_type dof_intr = (off_element < nbf ? basis.getDofTag(ord_intr)[3]   : 0);
      
          local2global[cnt][local] = off_element;
          off_element += dof_intr;
          Orientation::getElementNodeMap(subcell_verts, nids,
                                         cell, element,
                                         dim, i);
      
          if (!findNode(local2global[cnt][global], subcell_verts, nids, true)) {
            addNode(subcell_verts, nids, off_global);
            local2global[cnt][global] = off_global;
            off_global += dof_intr;
          }
          ++cnt;
        }
    
        // add the last offset
        local2global[cnt][local] = off_element;
        local2global[cnt][global] = -1; // invalid values
      }

      void getBoundaryFlags(ordinal_type *boundary,
                            const shards::CellTopology cell,
                            const ordinal_type *element) {
        ordinal_type subcell_verts[4], nids;
        const ordinal_type dim = cell.getDimension();
        const ordinal_type nside = cell.getSideCount();
        for (ordinal_type i=0;i<nside;++i) {
          Orientation::getElementNodeMap(subcell_verts, nids,
                                         cell, element,
                                         dim-1, i);
          
          if (!findBoundary(boundary[i], subcell_verts, nids)) {
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error,
                                        ">>> ERROR (Intrepid::HGRAD_TRI_Cn::Test 04): " \
                                        "Side node is not found");
          }
        }
      }

      ToyMesh() {
        reset();
      }

    };
  }

}

#endif

#endif
