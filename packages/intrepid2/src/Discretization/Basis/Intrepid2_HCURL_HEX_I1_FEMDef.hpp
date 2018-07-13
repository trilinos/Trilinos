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

/** \file   Intrepid2_HCURL_HEX_I1_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree 1 for H(curl) functions on HEX cells.
    \author Created by P. Bochev, D. Ridzal and K. Peterson.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HCURL_HEX_I1_FEM_DEF_HPP__
#define __INTREPID2_HCURL_HEX_I1_FEM_DEF_HPP__


namespace Intrepid2 {

  // -------------------------------------------------------------------------------------
  namespace Impl {

    template<EOperator opType>
    template<typename outputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HCURL_HEX_I1_FEM::Serial<opType>::
    getValues(       outputViewType output,
               const inputViewType input ) {

      switch (opType) {
      case OPERATOR_VALUE: {
        const auto x = input(0);
        const auto y = input(1);
        const auto z = input(2);

        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        output.access(0, 0) = (1.0 - y)*(1.0 - z)/8.0;
        output.access(0, 1) = 0.0;
        output.access(0, 2) = 0.0;

        output.access(1, 0) = 0.0;
        output.access(1, 1) = (1.0 + x)*(1.0 - z)/8.0;
        output.access(1, 2) = 0.0;

        output.access(2, 0) = -(1.0 + y)*(1.0 - z)/8.0;
        output.access(2, 1) = 0.0;
        output.access(2, 2) = 0.0;

        output.access(3, 0) = 0.0;
        output.access(3, 1) = -(1.0 - x)*(1.0 - z)/8.0;
        output.access(3, 2) = 0.0;

        output.access(4, 0) = (1.0 - y)*(1.0 + z)/8.0;
        output.access(4, 1) = 0.0;
        output.access(4, 2) = 0.0;

        output.access(5, 0) = 0.0;
        output.access(5, 1) = (1.0 + x)*(1.0 + z)/8.0;
        output.access(5, 2) = 0.0;

        output.access(6, 0) = -(1.0 + y)*(1.0 + z)/8.0;
        output.access(6, 1) = 0.0;
        output.access(6, 2) = 0.0;

        output.access(7, 0) = 0.0;
        output.access(7, 1) = -(1.0 - x)*(1.0 + z)/8.0;
        output.access(7, 2) = 0.0;

        output.access(8, 0) = 0.0;
        output.access(8, 1) = 0.0;
        output.access(8, 2) = (1.0 - x)*(1.0 - y)/8.0;

        output.access(9, 0) = 0.0;
        output.access(9, 1) = 0.0;
        output.access(9, 2) = (1.0 + x)*(1.0 - y)/8.0;

        output.access(10, 0) = 0.0;
        output.access(10, 1) = 0.0;
        output.access(10, 2) = (1.0 + x)*(1.0 + y)/8.0;

        output.access(11, 0) = 0.0;
        output.access(11, 1) = 0.0;
        output.access(11, 2) = (1.0 - x)*(1.0 + y)/8.0;
        break;
      }
      case OPERATOR_CURL: {
        const auto x = input(0);
        const auto y = input(1);
        const auto z = input(2);

        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        output.access(0, 0) = 0.0;
        output.access(0, 1) = -(1.0 - y)/8.0;
        output.access(0, 2) = (1.0 - z)/8.0;

        output.access(1, 0) = (1.0 + x)/8.0;
        output.access(1, 1) = 0.0;
        output.access(1, 2) = (1.0 - z)/8.0;

        output.access(2, 0) = 0.0;
        output.access(2, 1) = (1.0 + y)/8.0;
        output.access(2, 2) = (1.0 - z)/8.0;

        output.access(3, 0) = -(1.0 - x)/8.0;
        output.access(3, 1) = 0.0;
        output.access(3, 2) = (1.0 - z)/8.0;

        output.access(4, 0) = 0.0;
        output.access(4, 1) = (1.0 - y)/8.0;
        output.access(4, 2) = (1.0 + z)/8.0;

        output.access(5, 0) = -(1.0 + x)/8.0;
        output.access(5, 1) = 0.0;
        output.access(5, 2) = (1.0 + z)/8.0;

        output.access(6, 0) = 0.0;
        output.access(6, 1) = -(1.0 + y)/8.0;
        output.access(6, 2) = (1.0 + z)/8.0;

        output.access(7, 0) = (1.0 - x)/8.0;
        output.access(7, 1) = 0.0;
        output.access(7, 2) = (1.0 + z)/8.0;

        output.access(8, 0) = -(1.0 - x)/8.0;
        output.access(8, 1) = (1.0 - y)/8.0;
        output.access(8, 2) = 0.0;

        output.access(9, 0) = -(1.0 + x)/8.0;
        output.access(9, 1) = -(1.0 - y)/8.0;
        output.access(9, 2) = 0.0;

        output.access(10, 0) = (1.0 + x)/8.0;
        output.access(10, 1) = -(1.0 + y)/8.0;
        output.access(10, 2) = 0.0;

        output.access(11, 0) = (1.0 - x)/8.0;
        output.access(11, 1) = (1.0 + y)/8.0;
        output.access(11, 2) = 0.0;
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                  opType != OPERATOR_CURL,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_HEX_C1_FEM::Serial::getValues) operator is not supported");
      }
      } //end switch

    }

    template<typename SpT,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void
    Basis_HCURL_HEX_I1_FEM::
    getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const EOperator operatorType ) {
      typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
      typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
      typedef typename ExecSpace<typename inputPointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

      // Number of evaluation points = dim 0 of inputPoints
      const auto loopSize = inputPoints.extent(0);
      Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);

      switch (operatorType) {
      case OPERATOR_VALUE: {
        typedef Functor<outputValueViewType, inputPointViewType, OPERATOR_VALUE> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      case OPERATOR_CURL: {
        typedef Functor<outputValueViewType, inputPointViewType, OPERATOR_CURL> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      case OPERATOR_DIV: {
        INTREPID2_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                                      ">>> ERROR (Basis_HCURL_HEX_I1_FEM::getValues): DIV is invalid operator for HCURL Basis Functions");
        break;
      }

      case OPERATOR_GRAD: {
        INTREPID2_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_GRAD), std::invalid_argument,
                                      ">>> ERROR (Basis_HCURL_HEX_I1_FEM::getValues): GRAD is invalid operator for HCURL Basis Functions");
        break;
      }

      case OPERATOR_D1:
      case OPERATOR_D2:
      case OPERATOR_D3:
      case OPERATOR_D4:
      case OPERATOR_D5:
      case OPERATOR_D6:
      case OPERATOR_D7:
      case OPERATOR_D8:
      case OPERATOR_D9:
      case OPERATOR_D10: {
        INTREPID2_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_D1)    ||
                                      (operatorType == OPERATOR_D2)    ||
                                      (operatorType == OPERATOR_D3)    ||
                                      (operatorType == OPERATOR_D4)    ||
                                      (operatorType == OPERATOR_D5)    ||
                                      (operatorType == OPERATOR_D6)    ||
                                      (operatorType == OPERATOR_D7)    ||
                                      (operatorType == OPERATOR_D8)    ||
                                      (operatorType == OPERATOR_D9)    ||
                                      (operatorType == OPERATOR_D10),
                                      std::invalid_argument,
                                      ">>> ERROR (Basis_HCURL_HEX_I1_FEM::getValues): Invalid operator type");
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( (operatorType != OPERATOR_VALUE) &&
                                      (operatorType != OPERATOR_GRAD)  &&
                                      (operatorType != OPERATOR_CURL)  &&
                                      (operatorType != OPERATOR_DIV)   &&
                                      (operatorType != OPERATOR_D1)    &&
                                      (operatorType != OPERATOR_D2)    &&
                                      (operatorType != OPERATOR_D3)    &&
                                      (operatorType != OPERATOR_D4)    &&
                                      (operatorType != OPERATOR_D5)    &&
                                      (operatorType != OPERATOR_D6)    &&
                                      (operatorType != OPERATOR_D7)    &&
                                      (operatorType != OPERATOR_D8)    &&
                                      (operatorType != OPERATOR_D9)    &&
                                      (operatorType != OPERATOR_D10),
                                      std::invalid_argument,
                                      ">>> ERROR (Basis_HCURL_HEX_I1_FEM::getValues): Invalid operator type");
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------


  template<typename SpT, typename OT, typename PT>
  Basis_HCURL_HEX_I1_FEM<SpT,OT,PT>::
  Basis_HCURL_HEX_I1_FEM() {
    this->basisCardinality_  = 12;
    this->basisDegree_       = 1;
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
    this->basisType_         = BASIS_FEM_DEFAULT;
    this->basisCoordinates_  = COORDINATES_CARTESIAN;

    // initialize tags
    {
      // Basis-dependent intializations
      const ordinal_type tagSize  = 4;        // size of DoF tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

      // An array with local DoF tags assigned to basis functions, in the order of their local enumeration
      ordinal_type tags[48]  = { 1, 0, 0, 1,
                                 1, 1, 0, 1,
                                 1, 2, 0, 1,
                                 1, 3, 0, 1,
                                 1, 4, 0, 1,
                                 1, 5, 0, 1,
                                 1, 6, 0, 1,
                                 1, 7, 0, 1,
                                 1, 8, 0, 1,
                                 1, 9, 0, 1,
                                 1, 10, 0, 1,
                                 1, 11, 0, 1 };

      // when exec space is device, this wrapping relies on uvm.
      ordinal_type_array_1d_host tagView(&tags[0], 48);

      // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
      this->setOrdinalTagData(this->tagToOrdinal_,
                              this->ordinalToTag_,
                              tagView,
                              this->basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
    }

    // dofCoords on host and create its mirror view to device
    Kokkos::DynRankView<typename scalarViewType::value_type,typename SpT::array_layout,Kokkos::HostSpace>
      dofCoords("dofCoordsHost", this->basisCardinality_,this->basisCellTopology_.getDimension());

    dofCoords(0,0)  =  0.0;   dofCoords(0,1)  = -1.0;   dofCoords(0,2)  = -1.0;
    dofCoords(1,0)  =  1.0;   dofCoords(1,1)  =  0.0;   dofCoords(1,2)  = -1.0;
    dofCoords(2,0)  =  0.0;   dofCoords(2,1)  =  1.0;   dofCoords(2,2)  = -1.0;
    dofCoords(3,0)  = -1.0;   dofCoords(3,1)  =  0.0;   dofCoords(3,2)  = -1.0;
    dofCoords(4,0)  =  0.0;   dofCoords(4,1)  = -1.0;   dofCoords(4,2)  =  1.0;
    dofCoords(5,0)  =  1.0;   dofCoords(5,1)  =  0.0;   dofCoords(5,2)  =  1.0;
    dofCoords(6,0)  =  0.0;   dofCoords(6,1)  =  1.0;   dofCoords(6,2)  =  1.0;
    dofCoords(7,0)  = -1.0;   dofCoords(7,1)  =  0.0;   dofCoords(7,2)  =  1.0;
    dofCoords(8,0)  = -1.0;   dofCoords(8,1)  = -1.0;   dofCoords(8,2)  =  0.0;
    dofCoords(9,0)  =  1.0;   dofCoords(9,1)  = -1.0;   dofCoords(9,2)  =  0.0;
    dofCoords(10,0) =  1.0;   dofCoords(10,1) =  1.0;   dofCoords(10,2) =  0.0;
    dofCoords(11,0) = -1.0;   dofCoords(11,1) =  1.0;   dofCoords(11,2) =  0.0;

    this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);


    // dofCoeffs on host and create its mirror view to device
    Kokkos::DynRankView<typename scalarViewType::value_type,typename SpT::array_layout,Kokkos::HostSpace>
      dofCoeffs("dofCoeffsHost", this->basisCardinality_,this->basisCellTopology_.getDimension());

    // for HCURL_HEX_I1 dofCoeffs are the tangents on the hexahedron edges (with tangents magnitude equal to edges' lengths)
    dofCoeffs(0,0)  =  2.0;   dofCoeffs(0,1)  =  0.0;   dofCoeffs(0,2)  =  0.0;
    dofCoeffs(1,0)  =  0.0;   dofCoeffs(1,1)  =  2.0;   dofCoeffs(1,2)  =  0.0;
    dofCoeffs(2,0)  = -2.0;   dofCoeffs(2,1)  =  0.0;   dofCoeffs(2,2)  =  0.0;
    dofCoeffs(3,0)  =  0.0;   dofCoeffs(3,1)  = -2.0;   dofCoeffs(3,2)  =  0.0;
    dofCoeffs(4,0)  =  2.0;   dofCoeffs(4,1)  =  0.0;   dofCoeffs(4,2)  =  0.0;
    dofCoeffs(5,0)  =  0.0;   dofCoeffs(5,1)  =  2.0;   dofCoeffs(5,2)  =  0.0;
    dofCoeffs(6,0)  = -2.0;   dofCoeffs(6,1)  =  0.0;   dofCoeffs(6,2)  =  0.0;
    dofCoeffs(7,0)  =  0.0;   dofCoeffs(7,1)  = -2.0;   dofCoeffs(7,2)  =  0.0;
    dofCoeffs(8,0)  =  0.0;   dofCoeffs(8,1)  =  0.0;   dofCoeffs(8,2)  =  2.0;
    dofCoeffs(9,0)  =  0.0;   dofCoeffs(9,1)  =  0.0;   dofCoeffs(9,2)  =  2.0;
    dofCoeffs(10,0) =  0.0;   dofCoeffs(10,1) =  0.0;   dofCoeffs(10,2) =  2.0;
    dofCoeffs(11,0) =  0.0;   dofCoeffs(11,1) =  0.0;   dofCoeffs(11,2) =  2.0;

    this->dofCoeffs_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoeffs);
    Kokkos::deep_copy(this->dofCoeffs_, dofCoeffs);

  }

}// namespace Intrepid2

#endif
