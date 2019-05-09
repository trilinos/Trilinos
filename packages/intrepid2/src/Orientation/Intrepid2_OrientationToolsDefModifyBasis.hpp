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


/** \file   Intrepid2_OrientationToolsDefModifyBasis.hpp
    \brief  Definition file for the Intrepid2::OrientationTools class.
    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_MODIFY_BASIS_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_MODIFY_BASIS_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {
  
  template<typename SpT>
  template<typename ptViewType>
  KOKKOS_INLINE_FUNCTION
  bool OrientationTools<SpT>::
  isLeftHandedCell(const ptViewType pts) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_ABORT( pts.rank() != 2,  // npts x ndim
                              ">>> ERROR (Intrepid::OrientationTools::isLeftHandedCell): " \
                              "Point array is supposed to have rank 2.");
#endif
    typedef typename ptViewType::value_type value_type;
    
    const auto dim = pts.extent(1);
    value_type det = 0.0;
    switch (dim) {
    case 2: {
      // need 3 points (origin, x end point, y end point)
      const value_type v[2][2] = { { pts(1,0) - pts(0,0), pts(1,1) - pts(0,1) },
                                   { pts(2,0) - pts(0,0), pts(2,1) - pts(0,1) } };
      
      det = (v[0][0]*v[1][1] - v[1][0]*v[0][1]);
      break;
    }
    case 3: {
      // need 4 points (origin, x end point, y end point, z end point)
      const value_type v[3][3] = { { pts(1,0) - pts(0,0), pts(1,1) - pts(0,1), pts(1,2) - pts(0,2) },
                                   { pts(2,0) - pts(0,0), pts(2,1) - pts(0,1), pts(2,2) - pts(0,2) },
                                   { pts(3,0) - pts(0,0), pts(3,1) - pts(0,1), pts(3,2) - pts(0,2) } };
      
      det = (v[0][0] * v[1][1] * v[2][2] +
             v[0][1] * v[1][2] * v[2][0] +
             v[0][2] * v[1][0] * v[2][1] -
             v[0][2] * v[1][1] * v[2][0] -
             v[0][0] * v[1][2] * v[2][1] -
             v[0][1] * v[1][0] * v[2][2]);
      break;
    }
    default:{
      INTREPID2_TEST_FOR_ABORT( true, 
                                ">>> ERROR (Intrepid::Orientation::isLeftHandedCell): " \
                                "Dimension of points must be 2 or 3");
    }
    }
    return (det < 0.0);
  }
  
  template<typename SpT>
  template<typename elemOrtValueType, class ...elemOrtProperties,
           typename elemNodeValueType, class ...elemNodeProperties>
  void
  OrientationTools<SpT>::
  getOrientation(      Kokkos::DynRankView<elemOrtValueType,elemOrtProperties...> elemOrts,
                 const Kokkos::DynRankView<elemNodeValueType,elemNodeProperties...> elemNodes,
                 const shards::CellTopology cellTopo) {
    // small meta data modification and it uses shards; let's do this on host
    typedef typename Kokkos::Impl::is_space<SpT>::host_mirror_space::execution_space host_space_type;
    auto elemOrtsHost = Kokkos::create_mirror_view(typename host_space_type::memory_space(), elemOrts);
    auto elemNodesHost = Kokkos::create_mirror_view(typename host_space_type::memory_space(), elemNodes);
    
    const ordinal_type numCells = elemNodes.extent(0);
    for (auto cell=0;cell<numCells;++cell) {
      const auto nodes = Kokkos::subview(elemNodesHost, cell, Kokkos::ALL());
      elemOrtsHost(cell) = Orientation::getOrientation(cellTopo, nodes);
    }
    
    Kokkos::deep_copy(elemOrts, elemOrtsHost);
  }

  template<typename SpT>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties,
           typename ortValueType,    class ...ortProperties,
           typename BasisPtrType>
  void
  OrientationTools<SpT>::
  modifyBasisByOrientation(      Kokkos::DynRankView<outputValueType,outputProperties...> output,
                           const Kokkos::DynRankView<inputValueType, inputProperties...>  input,
                           const Kokkos::DynRankView<ortValueType,   ortProperties...>    orts,
                           const BasisPtrType basis ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( input.rank() != output.rank(), std::invalid_argument,
                                    ">>> ERROR (OrientationTools::modifyBasisByOrientation): Input and output rank are not 3.");
      for (size_type i=0;i<input.rank();++i)
        INTREPID2_TEST_FOR_EXCEPTION( input.extent(i) != output.extent(i), std::invalid_argument,
                                      ">>> ERROR (OrientationTools::modifyBasisByOrientation): Input and output dimension does not match.");

      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(input.extent(1)) != basis->getCardinality(), std::invalid_argument,
                                    ">>> ERROR (OrientationTools::modifyBasisByOrientation): Field dimension of input/output does not match to basis cardinality.");
    }
#endif
    typedef typename decltype(input)::non_const_value_type input_value_type;

    if (basis->requireOrientation()) {
      auto ordinalToTag = Kokkos::create_mirror_view(typename SpT::memory_space(), basis->getAllDofTags());
      auto tagToOrdinal = Kokkos::create_mirror_view(typename SpT::memory_space(), basis->getAllDofOrdinal());
      
      Kokkos::deep_copy(ordinalToTag, basis->getAllDofTags());
      Kokkos::deep_copy(tagToOrdinal, basis->getAllDofOrdinal());
      
      const ordinal_type 
        numCells  = output.extent(0),
        //numBasis  = output.extent(1),
        numPoints = output.extent(2),
        dimBasis  = output.extent(3); //returns 1 when output.rank() < 4;
      
      const CoeffMatrixDataViewType matData = createCoeffMatrix(basis);
      const shards::CellTopology cellTopo = basis->getBaseCellTopology();
      
      const ordinal_type 
        numVerts = cellTopo.getVertexCount()*ordinal_type(basis->getDofCount(0, 0) > 0),
        numEdges = cellTopo.getEdgeCount()*ordinal_type(basis->getDofCount(1, 0) > 0),
        numFaces = cellTopo.getFaceCount();
      
      for (auto cell=0;cell<numCells;++cell) {
        auto out = Kokkos::subview(output, cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto in  = Kokkos::subview(input,  cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        
        // vertex copy (no orientation)
        for (ordinal_type vertId=0;vertId<numVerts;++vertId) {
          const ordinal_type i = (static_cast<size_type>(vertId) < tagToOrdinal.extent(1) ? tagToOrdinal(0, vertId, 0) : -1);
          if (i != -1) // if dof does not exist i returns with -1
            for (ordinal_type j=0;j<numPoints;++j)
              for (ordinal_type k=0;k<dimBasis;++k)
                out(i, j, k) = in(i, j, k);
        }
        
        // interior copy
        {
          const ordinal_type cellDim = cellTopo.getDimension();
          const ordinal_type ordIntr = (static_cast<size_type>(cellDim) < tagToOrdinal.extent(0) ? tagToOrdinal(cellDim, 0, 0) : -1);
          if (ordIntr != -1) {
            const ordinal_type ndofIntr = ordinalToTag(ordIntr, 3);
            for (ordinal_type i=0;i<ndofIntr;++i) {
              const ordinal_type ii = tagToOrdinal(cellDim, 0, i);
              for (ordinal_type j=0;j<numPoints;++j)
                for (ordinal_type k=0;k<dimBasis;++k)
                  out(ii, j, k) = in(ii, j, k);
            }
          }
        }
        
        // edge transformation
        ordinal_type existEdgeDofs = 0;
        if (numEdges > 0) {
          ordinal_type ortEdges[12];
          orts(cell).getEdgeOrientation(ortEdges, numEdges);
          
          // apply coeff matrix
          for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId) {
            const ordinal_type ordEdge = (1 < tagToOrdinal.extent(0) ? (static_cast<size_type>(edgeId) < tagToOrdinal.extent(1) ? tagToOrdinal(1, edgeId, 0) : -1) : -1);
            
            if (ordEdge != -1) {
              existEdgeDofs = 1;
              const ordinal_type ndofEdge = ordinalToTag(ordEdge, 3);
              const auto mat = Kokkos::subview(matData, 
                                               edgeId, ortEdges[edgeId], 
                                               Kokkos::ALL(), Kokkos::ALL());
              
              for (ordinal_type j=0;j<numPoints;++j) 
                for (ordinal_type i=0;i<ndofEdge;++i) {
                  const ordinal_type ii = tagToOrdinal(1, edgeId, i);
                  
                  for (ordinal_type k=0;k<dimBasis;++k) {
                    input_value_type temp = 0.0;
                    for (ordinal_type l=0;l<ndofEdge;++l) {
                      const ordinal_type ll = tagToOrdinal(1, edgeId, l);
                      temp += mat(i,l)*in(ll, j, k);
                    }
                    out(ii, j, k) = temp;
                  }
                }
            }
          }
        }
        
        // face transformation
        if (numFaces > 0) {
          ordinal_type ortFaces[12];
          orts(cell).getFaceOrientation(ortFaces, numFaces);
          
          // apply coeff matrix
          for (ordinal_type faceId=0;faceId<numFaces;++faceId) {
            const ordinal_type ordFace = (2 < tagToOrdinal.extent(0) ? (static_cast<size_type>(faceId) < tagToOrdinal.extent(1) ? tagToOrdinal(2, faceId, 0) : -1) : -1);
            
            if (ordFace != -1) {
              const ordinal_type ndofFace = ordinalToTag(ordFace, 3);
              const auto mat = Kokkos::subview(matData, 
                                               numEdges*existEdgeDofs+faceId, ortFaces[faceId], 
                                               Kokkos::ALL(), Kokkos::ALL());
              
              for (ordinal_type j=0;j<numPoints;++j) 
                for (ordinal_type i=0;i<ndofFace;++i) {
                  const ordinal_type ii = tagToOrdinal(2, faceId, i);
                  
                  for (ordinal_type k=0;k<dimBasis;++k) {
                    input_value_type temp = 0.0;
                    for (ordinal_type l=0;l<ndofFace;++l) {
                      const ordinal_type ll = tagToOrdinal(2, faceId, l);
                      temp += mat(i,l)*in(ll, j, k);
                    }
                    out(ii, j, k) = temp;
                  }
                }
            }
          }
        }
        
      }
    } else {
      Kokkos::deep_copy(output, input);      
    }
  }
}

#endif
