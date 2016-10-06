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
  template<typename outValueValueType, class ...outValueProperties,
           typename refValueValueType, class ...refValueProperties,
           typename elemOrtValueType,  class ...elemOrtProperties,
           typename quadBasisPtrType>
  void
  OrientationTools<SpT>::  
  getModifiedQuadrilateralBasis(/**/  Kokkos::DynRankView<outValueValueType,outValueProperties...> outValues,
                                const Kokkos::DynRankView<refValueValueType,refValueProperties...> refValues,
                                const Kokkos::DynRankView<elemOrtValueType,elemOrtProperties...> elemOrts,
                                const quadBasisPtrType quadBasis) {
    const ordinal_type order = quadBasis->getDegree();
    
    auto ordinalToTag = Kokkos::create_mirror_view(typename SpT::memory_space(), quadBasis->getAllDofTags());
    auto tagToOrdinal = Kokkos::create_mirror_view(typename SpT::memory_space(), quadBasis->getAllDofOrdinal());

    Kokkos::deep_copy(ordinalToTag, quadBasis->getAllDofTags());
    Kokkos::deep_copy(tagToOrdinal, quadBasis->getAllDofOrdinal());

    const ordinal_type numCells = outValues.dimension(0);
    const ordinal_type ndofBasis = outValues.dimension(1);
    const ordinal_type dofDim = outValues.dimension(2);

    auto matData = Kokkos::subview(quadEdgeData,
                                   static_cast<ordinal_type>(FUNCTION_SPACE_HGRAD), order - 1,
                                   Kokkos::ALL(), Kokkos::ALL(),
                                   Kokkos::ALL(), Kokkos::ALL());

    const ordinal_type numVerts = 4, numEdges = 4;
    for (auto cell=0;cell<numCells;++cell) {
      ordinal_type orts[numEdges];
      elemOrts(cell).getEdgeOrientation(orts, numEdges);

      auto out = Kokkos::subview(outValues, cell, Kokkos::ALL(), Kokkos::ALL());
      auto ref = Kokkos::subview(refValues, cell, Kokkos::ALL(), Kokkos::ALL());

      // vertex copy
      for (auto vertId=0;vertId<numVerts;++vertId) {
        const auto ii = tagToOrdinal(0, vertId, 0);
        for (auto j=0;j<dofDim;++j)
          out(ii, j) = ref(ii, j);
      }

      // interior copy
      {
        const auto ndofIntr = ordinalToTag(tagToOrdinal(2, 0, 0), 3);
        for (auto i=0;i<ndofIntr;++i) {
          const auto ii = tagToOrdinal(2, 0, i);
          for (auto j=0;j<dofDim;++j)
            out(ii, j) = ref(ii, j);
        }
      }

      // apply coeff matrix
      for (auto edgeId=0;edgeId<numEdges;++edgeId) {
        const auto ndofEdge = ordinalToTag(tagToOrdinal(1, edgeId, 0), 3);
        const auto mat = Kokkos::subview(matData, edgeId, orts[edgeId], Kokkos::ALL(), Kokkos::ALL());

        for (auto j=0;j<dofDim;++j) 
          for (auto i=0;i<ndofEdge;++i) {
            const auto ii = tagToOrdinal(1, edgeId, i);

            auto temp = 0.0;
            for (auto k=0;k<ndofEdge;++k) {
              const auto kk = tagToOrdinal(1, edgeId, k);
              temp += mat(i,k)*ref(kk, j);
            }
            out(ii, j) = temp;
          }
      }
    }
  }

}

#endif
