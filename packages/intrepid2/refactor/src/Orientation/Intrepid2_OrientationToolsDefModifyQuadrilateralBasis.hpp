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


/** \file   Intrepid_OrientationToolsDef.hpp
    \brief  Definition file for the Intrepid2::OrientationTools class.
    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_MODIFY_QUADRILATERAL_BASIS_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_MODIFY_QUADRILATERAL_BASIS_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {
  
  template<typename SpT>
  template<typename outputValueType,  class ...outputProperties,
           typename inputValueType,   class ...inputProperties,
           typename elemOrtValueType, class ...elemOrtProperties,
           typename quadBasisPtrType>
  void
  OrientationTools<SpT>::  
  getModifiedBasis(/**/  Kokkos::DynRankView<outputValueType,outputProperties...> output,
                                const Kokkos::DynRankView<inputValueType,inputProperties...> input,
                                const Kokkos::DynRankView<elemOrtValueType,elemOrtProperties...> elemOrts,
                                const quadBasisPtrType basis) {
    auto ordinalToTag = Kokkos::create_mirror_view(typename SpT::memory_space(), basis->getAllDofTags());
    auto tagToOrdinal = Kokkos::create_mirror_view(typename SpT::memory_space(), basis->getAllDofOrdinal());
    
    Kokkos::deep_copy(ordinalToTag, basis->getAllDofTags());
    Kokkos::deep_copy(tagToOrdinal, basis->getAllDofOrdinal());
    
    const ordinal_type 
      numCells  = output.dimension(0),
      numBasis  = output.dimension(1),
      numPoints = output.dimension(2),
      dimBasis  = output.dimension(3);
    
    std::pair<std::string,int> key(basis->getName(), basis->getDegree());
    const auto found = ortData.find(key);

    MatDataViewType matData;
    if (found == ortData.end()) {
      matData = createQuadrilateralMatrixData(basis);
      ortData.insert(key, matData);
    } else {
      matData = found->second;
    }

    const shards::CellTopology cellTopo = basis.getBaseCellTopology();

    const ordinal_type 
      numVerts = cellTopo.getVertexCount(), 
      numEdges = cellTopo.getEdgeCount(),
      numFaces = cellTopo.getFaceCount();

    const ordinal_type intrDim = ( numEdges == 0 ? 1 : (numFaces == 0 ? 1 : 2) );

    for (auto cell=0;cell<numCells;++cell) {
      auto out = Kokkos::subview(output, cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
      auto in  = Kokkos::subview(input,  cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

      // vertex copy (no orientation)
      for (auto vertId=0;vertId<numVerts;++vertId) {
        const auto i = tagToOrdinal(0, vertId, 0);
        if (i != -1) // if dof does not exist i returns with -1
          for (auto j=0;j<numPoints;++j)
            for (auto k=0;k<dimBasis;++k)
              out(i, j, k) = in(i, j, k);
      }
      
      // interior copy
      {
        const auto ordIntr = tagToOrdinal(intrDim, 0, 0);
        if (ordIntr != -1) {
          const auto ndofIntr = ordinalToTag(ordIntr, 3);
          for (auto i=0;i<ndofIntr;++i) {
            const auto ii = tagToOrdinal(intrDim, 0, i);
            for (auto j=0;j<numPoints;++j)
              for (auto k=0;k<dimBasis;++k)
                out(ii, j, k) = in(ii, j, k);
          }
        }
      }

      // edge transformation
      if (numEdges > 0) {
        ordinal_type ortEdges[12];
        elemOrts(cell).getEdgeOrientation(ortEdges, numEdges);
        
        // apply coeff matrix
        for (auto edgeId=0;edgeId<numEdges;++edgeId) {
          const auto ordEdge = tagToOrdinal(1, edgeId, 0);
          
          if (ordEdge != -1) {
            const auto ndofEdge = ordinalToTag(ordEdge, 3);
            const auto mat = Kokkos::subview(matData, 
                                             edgeId, ortEdges[edgeId], 
                                             Kokkos::ALL(), Kokkos::ALL());
            
            for (auto j=0;j<numPoints;++j) 
              for (auto i=0;i<ndofEdge;++i) {
                const auto ii = tagToOrdinal(1, edgeId, i);
                
                for (auto k=0;k<dimBasis;++k) {
                  double temp = 0.0;
                  for (auto l=0;l<ndofEdge;++l) {
                    const auto ll = tagToOrdinal(1, edgeId, l);
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
        elemOrts(cell).getFaceOrientation(ortFaces, numFaces);
        
        // apply coeff matrix
        for (auto faceId=0;faceId<numFaces;++faceId) {
          const auto ordFace = tagToOrdinal(2, faceId, 0);
          
          if (ordFace != -1) {
            const auto ndofFace = ordinalToTag(ordFace, 3);
            const auto mat = Kokkos::subview(matData, 
                                             numEdges+faceId, ortFaces[faceId], 
                                             Kokkos::ALL(), Kokkos::ALL());
            
            for (auto j=0;j<numPoints;++j) 
              for (auto i=0;i<ndofFace;++i) {
                const auto ii = tagToOrdinal(2, faceId, i);
                
                for (auto k=0;k<dimBasis;++k) {
                  double temp = 0.0;
                  for (auto l=0;l<ndofFace;++l) {
                    const auto ll = tagToOrdinal(2, faceId, l);
                    temp += mat(i,l)*in(ll, j, k);
                  }
                  out(ii, j, k) = temp;
                }
              }
          }
        }
      }
      
    }
  }

}

#endif
