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


/** \file   Intrepid2_OrientationDef.hpp
    \brief  Definition file for the Intrepid2::Orientation class.
    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATION_DEF_HPP__
#define __INTREPID2_ORIENTATION_DEF_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  // ------------------------------------------------------------------------------------
  // Orientation
  //
  //
  template<typename cellVertViewType>
  inline
  void
  Orientation::getCellVertexMap(typename cellVertViewType::non_const_value_type *subCellVerts,
                                 ordinal_type &numVerts,
                                 const shards::CellTopology cellTopo,
                                 const cellVertViewType cellVertices,
                                 const ordinal_type subCellDim,
                                 const ordinal_type subCellOrd) {
    static_assert(Kokkos::Impl::MemorySpaceAccess
                  <Kokkos::HostSpace,typename cellVertViewType::device_type::memory_space>::accessible,
                  "host space cannot access cellVertViewType");
    switch (subCellDim) {
    case 0: {
      numVerts = 1;
      subCellVerts[0] = cellVertices(subCellOrd);
      break;
    }
    default: {
      numVerts = cellTopo.getVertexCount(subCellDim, subCellOrd);
      for (ordinal_type i=0;i<numVerts;++i)
        subCellVerts[i] = cellVertices(cellTopo.getNodeMap(subCellDim, subCellOrd, i));
      break;
    }
    }
  }
  
  template<typename subCellVertType>
  inline
  ordinal_type
  Orientation::getOrientation(const subCellVertType subCellVerts[],
                              const ordinal_type numVerts) {
    ordinal_type ort = 0;
    switch (numVerts) {
    case 2: {// edge
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( ( subCellVerts[0] == subCellVerts[1] ), 
                                ">>> ERROR (Intrepid::Orientation::getOrientation): " \
                                "Invalid subCellVerts, same vertex ids are repeated");
#endif
      ort = (subCellVerts[0] > subCellVerts[1]);
      break;
    }
    case 3: {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( ( subCellVerts[0] == subCellVerts[1] ||
                                  subCellVerts[0] == subCellVerts[2] ||
                                  subCellVerts[1] == subCellVerts[2] ), 
                                ">>> ERROR (Intrepid::Orientation::getOrientation): " \
                                "Invalid subCellVerts, same vertex ids are repeated");
#endif
      ordinal_type rotation = 0; // find smallest vertex id
      for (ordinal_type i=1;i<3;++i)
        rotation = ( subCellVerts[i] < subCellVerts[rotation] ? i : rotation );

      const ordinal_type axes[][2] = { {1,2}, {2,0}, {0,1} };
      const ordinal_type flip = (subCellVerts[axes[rotation][0]] > subCellVerts[axes[rotation][1]]);

      ort = flip*3 + rotation;
      break;
    }
    case 4: {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( ( subCellVerts[0] == subCellVerts[1] ||
                                  subCellVerts[0] == subCellVerts[2] ||
                                  subCellVerts[0] == subCellVerts[3] ||
                                  subCellVerts[1] == subCellVerts[2] ||
                                  subCellVerts[1] == subCellVerts[3] ||
                                  subCellVerts[2] == subCellVerts[3] ), 
                                ">>> ERROR (Intrepid::Orientation::getOrientation): " \
                                "Invalid subCellVerts, same vertex ids are repeated");
#endif
      ordinal_type rotation = 0; // find smallest vertex id
      for (ordinal_type i=1;i<4;++i)
        rotation = ( subCellVerts[i] < subCellVerts[rotation] ? i : rotation );

      const ordinal_type axes[][2] = { {1,3}, {2,0}, {3,1}, {0,2} };
      const ordinal_type flip = (subCellVerts[axes[rotation][0]] > subCellVerts[axes[rotation][1]]);

      ort = flip*4 + rotation;
      break;
    }
    default: {
      INTREPID2_TEST_FOR_ABORT( true, 
                                ">>> ERROR (Intrepid::Orientation::getOrientation): " \
                                "Invalid numVerts (2 (edge),3 (triangle) and 4 (quadrilateral) are allowed)");
      break;
    }
    }
    return ort;
  }

  template<typename cellVertViewType>
  inline
  Orientation
  Orientation::getOrientation(const shards::CellTopology cellTopo,
                              const cellVertViewType cellVertices,
                              bool isSide) {
    static_assert(Kokkos::Impl::MemorySpaceAccess
                  <Kokkos::HostSpace,typename cellVertViewType::device_type::memory_space>::accessible,
                  "host space cannot access cellVertViewType");

    Orientation ort;
    auto dim = cellTopo.getDimension();
    const ordinal_type nedge = (isSide && dim==1) ? 1 : cellTopo.getEdgeCount();

    if (nedge > 0) {
      typename cellVertViewType::non_const_value_type vertsSubCell[2];
      ordinal_type orts[12], nvertSubCell;
      for (ordinal_type i=0;i<nedge;++i) {
        Orientation::getCellVertexMap(vertsSubCell,
                                       nvertSubCell,
                                       cellTopo,
                                       cellVertices,
                                       1, i);
        orts[i] = Orientation::getOrientation(vertsSubCell, nvertSubCell);
      }
      ort.setEdgeOrientation(nedge, orts);
    }
    const ordinal_type nface = (isSide && dim==2) ? 1 : cellTopo.getFaceCount();
    if (nface > 0) {
      typename cellVertViewType::non_const_value_type vertsSubCell[4];
      ordinal_type orts[6], nvertSubCell;
      for (ordinal_type i=0;i<nface;++i) {
        Orientation::getCellVertexMap(vertsSubCell,
                                       nvertSubCell,
                                       cellTopo,
                                       cellVertices,
                                       2, i);
        orts[i] = Orientation::getOrientation(vertsSubCell, nvertSubCell);
      }
      ort.setFaceOrientation(nface, orts);
    }
    return ort;
  }

  inline
  ordinal_type 
  Orientation::getEdgeOrdinalOfFace(const ordinal_type subsubcellOrd,
                                    const ordinal_type subcellOrd,
                                    const shards::CellTopology cellTopo) {
    ordinal_type r_val = -1;

    const auto cellBaseKey    = cellTopo.getBaseKey();
    if        (cellBaseKey == shards::Hexahedron<>::key) {
      INTREPID2_TEST_FOR_EXCEPTION( !(subcellOrd < 6) && 
                                    !(subsubcellOrd < 4),
                                    std::logic_error,
                                    "subcell and subsubcell information are not correct" );
      const int quad_to_hex_edges[6][4] = { { 0, 9, 4, 8 },
                                            { 1,10, 5, 9 },
                                            { 2,11, 6,10 },
                                            { 8, 7,11, 3 },
                                            { 3, 2, 1, 0 },
                                            { 4, 5, 6, 7 } };
      r_val = quad_to_hex_edges[subcellOrd][subsubcellOrd];
    } else if (cellBaseKey == shards::Tetrahedron<>::key) {
      INTREPID2_TEST_FOR_EXCEPTION( !(subcellOrd < 4) && 
                                    !(subsubcellOrd < 3),
                                    std::logic_error,
                                    "subcell and subsubcell information are not correct" );
      const ordinal_type tri_to_tet_edges[4][3] = { { 0, 4, 3 },
                                                    { 1, 5, 4 },
                                                    { 3, 5, 2 },
                                                    { 2, 1, 0 } };
      r_val = tri_to_tet_edges[subcellOrd][subsubcellOrd];      
    } else {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    "cellTopo is not supported: try TET and HEX" );
    }
    return r_val;
  }
  
  KOKKOS_INLINE_FUNCTION
  Orientation::Orientation()
    : _edgeOrt(0), _faceOrt(0) {}

  KOKKOS_INLINE_FUNCTION
  bool
  Orientation::isAlignedToReference() const {
    return (_edgeOrt == 0 && _faceOrt == 0);
  }

  KOKKOS_INLINE_FUNCTION
  void
  Orientation::setEdgeOrientation(const ordinal_type numEdge, const ordinal_type edgeOrt[]) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_ABORT( !((numEdge == 1) || (3 <= numEdge && numEdge <= 12 )),
                              ">>> ERROR (Intrepid::Orientation::setEdgeOrientation): " \
                              "Invalid numEdge");
#endif
    _edgeOrt = 0;
    for (ordinal_type i=0;i<numEdge;++i)
      _edgeOrt |= (edgeOrt[i] & 1) << i;
  }

  KOKKOS_INLINE_FUNCTION
  void
  Orientation::getEdgeOrientation(ordinal_type *edgeOrt, const ordinal_type numEdge) const {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_ABORT( !((numEdge == 1) || (3 <= numEdge && numEdge <= 12 )),
                              ">>> ERROR (Intrepid::Orientation::setEdgeOrientation): " \
                              "Invalid numEdge");
#endif
    for (ordinal_type i=0;i<numEdge;++i)
      edgeOrt[i] = (_edgeOrt & (1 << i)) >> i;
  }

  KOKKOS_INLINE_FUNCTION
  void
  Orientation::setFaceOrientation(const ordinal_type numFace, const ordinal_type faceOrt[]) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_ABORT( !((numFace == 1) || (4 <= numFace && numFace <= 6 )),
                              ">>> ERROR (Intrepid::Orientation::setFaceOrientation): "
                              "Invalid numFace");
#endif
    _faceOrt = 0;
    for (ordinal_type i=0;i<numFace;++i) {
      const ordinal_type s = i*3;
      _faceOrt |= (faceOrt[i] & 7) << s;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void
  Orientation::getFaceOrientation(ordinal_type *faceOrt, const ordinal_type numFace) const {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_ABORT( !((numFace == 1) || (4 <= numFace && numFace <= 6 )),
                              ">>> ERROR (Intrepid::Orientation::setEdgeOrientation): "
                              "Invalid numFace");
#endif
    for (ordinal_type i=0;i<numFace;++i) {
      const ordinal_type s = i*3;
      faceOrt[i] = (_faceOrt & (7 << s)) >> s;
    }
  }

  inline std::string Orientation::to_string() const {
    return "Orientation{ face: " + std::to_string(_faceOrt) + "; edge: " + std::to_string(_edgeOrt) + " }";
  }
}

#endif
