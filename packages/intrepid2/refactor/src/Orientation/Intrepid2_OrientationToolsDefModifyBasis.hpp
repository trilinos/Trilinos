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
    
    const auto dim = pts.dimension(1);
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
  getOrientation(/**/  Kokkos::DynRankView<elemOrtValueType,elemOrtProperties...> elemOrts,
                 const Kokkos::DynRankView<elemNodeValueType,elemNodeProperties...> elemNodes,
                 const shards::CellTopology cellTopo) {
    // small meta data modification and it uses shards; let's do this on host
    typedef typename Kokkos::Impl::is_space<SpT>::host_mirror_space::execution_space host_space_type;
    auto elemOrtsHost = Kokkos::create_mirror_view(typename host_space_type::memory_space(), elemOrts);
    auto elemNodesHost = Kokkos::create_mirror_view(typename host_space_type::memory_space(), elemNodes);
    
    const ordinal_type numCells = elemNodes.dimension(0);
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
  modifyBasisByOrientation(/**/  Kokkos::DynRankView<outputValueType,outputProperties...> output,
                           const Kokkos::DynRankView<inputValueType, inputProperties...>  input,
                           const Kokkos::DynRankView<ortValueType,   ortProperties...>    orts,
                           const BasisPtrType basis ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( input.rank() != output.rank(), std::invalid_argument,
                                    ">>> ERROR (OrientationTools::modifyBasisByOrientation): Input and output rank are not 3.");
      for (ordinal_type i=0;i<input.rank();++i)
        INTREPID2_TEST_FOR_EXCEPTION( input.dimension(i) != output.dimension(i), std::invalid_argument,
                                      ">>> ERROR (OrientationTools::modifyBasisByOrientation): Input and output dimension does not match.");

      INTREPID2_TEST_FOR_EXCEPTION( input.dimension(1) != basis->getCardinality(), std::invalid_argument,
                                    ">>> ERROR (OrientationTools::modifyBasisByOrientation): Field dimension of input/output does not match to basis cardinality.");
      INTREPID2_TEST_FOR_EXCEPTION( input.dimension(3) != basis->getBaseCellTopology().getDimension(), std::invalid_argument,
                                    ">>> ERROR (OrientationTools::modifyBasisByOrientation): Space dimension of input/output does not match to topology dimension.");
    }
#endif
    const auto cellTopo = basis->getBaseCellTopology();
    bool is_ort_applied = false;
    switch ( cellTopo.getKey() ) {
    case shards::Line<2>::key:
      // do nothing
      break;
    case shards::Quadrilateral<4>::key: {
      if (basis->requireOrientation()) {
        OrientationTools<SpT>::getModifiedQuadrilateralBasis(output, input, orts, basis);
        is_ort_applied = true;
      }
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                    ">>> ERROR (OrientationTools::HGRADtransformVALUE): CellTopology is not supported.");
    }
    }

    // if orientation is not applied, copy input to output
    if (!is_ort_applied)
      Kokkos::deep_copy(output, input);
  }


}

#endif
