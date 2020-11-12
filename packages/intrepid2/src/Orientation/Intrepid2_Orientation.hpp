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

/** \file   Intrepid2_Orientation.hpp
    \brief  Header file for the Intrepid2::Orientation class.
    \author Created by Kyungjoo Kim
*/

#ifndef __INTREPID2_ORIENTATION_HPP__
#define __INTREPID2_ORIENTATION_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Shards_CellTopology.hpp"

namespace Intrepid2 {

  /**
    \brief Orientation encoding and decoding

     Use input and output as pointer arrays which assumes that input/output
     are located on the stack and contiguous.
  */
  class Orientation {
  private:
    template<typename elemNodeViewType>
    static void getElementNodeMap(typename elemNodeViewType::non_const_value_type *subCellVerts,
                                  ordinal_type &numVerts,
                                  const shards::CellTopology cellTopo,
                                  const elemNodeViewType elemNodes,
                                  const ordinal_type subCellDim,
                                  const ordinal_type subCellOrd);
    
    // orientation is always computed on the right-handed coordinates
    template<typename subCellVertType>
    static ordinal_type getOrientation(const subCellVertType subCellVerts[],
                                       const ordinal_type numVerts);


  public:
    template<typename elemNodeViewType>
    static Orientation getOrientation(const shards::CellTopology cellTopo,
                                      const elemNodeViewType elemNodes);
    
    static ordinal_type getEdgeOrdinalOfFace(const ordinal_type subsubcellOrd,
                                             const ordinal_type subcellOrd,
                                             const shards::CellTopology cellTopo);


    /*
    Function Removed. Use instead Impl::OrientationTools::getRefSubcellTangents
    template<typename refTanViewType>
    static void getReferenceEdgeTangent(const refTanViewType &tanE,
                                        const ordinal_type subcellOrd,
                                        const shards::CellTopology cellTopo,
                                        const ordinal_type ort,
                                        const bool is_normalize = true);
    */
    
    /*
    Function Removed. Use instead Impl::OrientationTools::getRefSubcellTangents
    template<typename refTanViewType>
    static void getReferenceFaceTangents(const refTanViewType &tanU,
                                         const refTanViewType &tanV,
                                         const ordinal_type subcellOrd,
                                         const shards::CellTopology cellTopo,
                                         const ordinal_type ort,
                                         const bool is_normalize = true);
    */

    /*
    Function Removed. Use instead Impl::OrientationTools::getRefSideTangentsAndNormal
    template<typename refNormalViewType>
    static void getReferenceFaceNormal(const refNormalViewType &normalV,
                                       const ordinal_type subcellOrd,
                                       const shards::CellTopology cellTopo,
                                       const ordinal_type ort,
                                       const bool is_normalize = true);
    */

      
  private:
    ordinal_type _edgeOrt, _faceOrt;
    
  public:
    KOKKOS_INLINE_FUNCTION
    Orientation();

    KOKKOS_DEFAULTED_FUNCTION
    Orientation(const Orientation &b) = default;

    KOKKOS_INLINE_FUNCTION
    bool isAlignedToReference() const;

    KOKKOS_INLINE_FUNCTION
    void setEdgeOrientation(const ordinal_type numEdge, const ordinal_type edgeOrt[]);

    KOKKOS_INLINE_FUNCTION
    void getEdgeOrientation(ordinal_type *edgeOrt, const ordinal_type numEdge) const;

    KOKKOS_INLINE_FUNCTION
    void setFaceOrientation(const ordinal_type numFace, const ordinal_type faceOrt[]);

    KOKKOS_INLINE_FUNCTION
    void getFaceOrientation(ordinal_type *faceOrt, const ordinal_type numFace) const;
  };

}

// include templated function definitions
#include "Intrepid2_OrientationDef.hpp"

#endif
