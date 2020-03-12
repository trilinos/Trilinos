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

/** \file   Intrepid2_HGRAD_HEX_C1_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree 1 for H(grad) functions on HEX cells.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
 */

#ifndef __INTREPID2_HGRAD_HEX_C1_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_HEX_C1_FEM_DEF_HPP__

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------
  namespace Impl {

    template<EOperator opType>
    template<typename OutputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_HEX_C1_FEM::Serial<opType>::
    getValues(       OutputViewType output,
               const inputViewType input ) {
      switch (opType) {
      case OPERATOR_VALUE : {
        const auto x = input(0);
        const auto y = input(1);
        const auto z = input(2);

        // output is a rank-1 array with dimensions (basisCardinality_)
        output.access(0) = (1.0 - x)*(1.0 - y)*(1.0 - z)/8.0;
        output.access(1) = (1.0 + x)*(1.0 - y)*(1.0 - z)/8.0;
        output.access(2) = (1.0 + x)*(1.0 + y)*(1.0 - z)/8.0;
        output.access(3) = (1.0 - x)*(1.0 + y)*(1.0 - z)/8.0;

        output.access(4) = (1.0 - x)*(1.0 - y)*(1.0 + z)/8.0;
        output.access(5) = (1.0 + x)*(1.0 - y)*(1.0 + z)/8.0;
        output.access(6) = (1.0 + x)*(1.0 + y)*(1.0 + z)/8.0;
        output.access(7) = (1.0 - x)*(1.0 + y)*(1.0 + z)/8.0;
        break;
      }
      case OPERATOR_GRAD : {
        const auto x = input(0);
        const auto y = input(1);
        const auto z = input(2);

        // output is a rank-2 array with dimensions (basisCardinality_, spaceDim)
        output.access(0, 0) = -(1.0 - y)*(1.0 - z)/8.0;
        output.access(0, 1) = -(1.0 - x)*(1.0 - z)/8.0;
        output.access(0, 2) = -(1.0 - x)*(1.0 - y)/8.0;

        output.access(1, 0) =  (1.0 - y)*(1.0 - z)/8.0;
        output.access(1, 1) = -(1.0 + x)*(1.0 - z)/8.0;
        output.access(1, 2) = -(1.0 + x)*(1.0 - y)/8.0;

        output.access(2, 0) =  (1.0 + y)*(1.0 - z)/8.0;
        output.access(2, 1) =  (1.0 + x)*(1.0 - z)/8.0;
        output.access(2, 2) = -(1.0 + x)*(1.0 + y)/8.0;

        output.access(3, 0) = -(1.0 + y)*(1.0 - z)/8.0;
        output.access(3, 1) =  (1.0 - x)*(1.0 - z)/8.0;
        output.access(3, 2) = -(1.0 - x)*(1.0 + y)/8.0;

        output.access(4, 0) = -(1.0 - y)*(1.0 + z)/8.0;
        output.access(4, 1) = -(1.0 - x)*(1.0 + z)/8.0;
        output.access(4, 2) =  (1.0 - x)*(1.0 - y)/8.0;

        output.access(5, 0) =  (1.0 - y)*(1.0 + z)/8.0;
        output.access(5, 1) = -(1.0 + x)*(1.0 + z)/8.0;
        output.access(5, 2) =  (1.0 + x)*(1.0 - y)/8.0;

        output.access(6, 0) =  (1.0 + y)*(1.0 + z)/8.0;
        output.access(6, 1) =  (1.0 + x)*(1.0 + z)/8.0;
        output.access(6, 2) =  (1.0 + x)*(1.0 + y)/8.0;

        output.access(7, 0) = -(1.0 + y)*(1.0 + z)/8.0;
        output.access(7, 1) =  (1.0 - x)*(1.0 + z)/8.0;
        output.access(7, 2) =  (1.0 - x)*(1.0 + y)/8.0;
        break;
      }
      case OPERATOR_D2 : {
        const auto x = input(0);
        const auto y = input(1);
        const auto z = input(2);

        // output is a rank-2 array with dimensions (basisCardinality_, D2Cardinality = 6)
        output.access(0, 0) =  0.0;                    // {2, 0, 0}
        output.access(0, 1) =  (1.0 - z)/8.0;          // {1, 1, 0}
        output.access(0, 2) =  (1.0 - y)/8.0;          // {1, 0, 1}
        output.access(0, 3) =  0.0;                    // {0, 2, 0}
        output.access(0, 4) =  (1.0 - x)/8.0;          // {0, 1, 1}
        output.access(0, 5) =  0.0;                    // {0, 0, 2}

        output.access(1, 0) =  0.0;                    // {2, 0, 0}
        output.access(1, 1) = -(1.0 - z)/8.0;          // {1, 1, 0}
        output.access(1, 2) = -(1.0 - y)/8.0;          // {1, 0, 1}
        output.access(1, 3) =  0.0;                    // {0, 2, 0}
        output.access(1, 4) =  (1.0 + x)/8.0;          // {0, 1, 1}
        output.access(1, 5) =  0.0;                    // {0, 0, 2}

        output.access(2, 0) =  0.0;                    // {2, 0, 0}
        output.access(2, 1) =  (1.0 - z)/8.0;          // {1, 1, 0}
        output.access(2, 2) = -(1.0 + y)/8.0;          // {1, 0, 1}
        output.access(2, 3) =  0.0;                    // {0, 2, 0}
        output.access(2, 4) = -(1.0 + x)/8.0;          // {0, 1, 1}
        output.access(2, 5) =  0.0;                    // {0, 0, 2}

        output.access(3, 0) =  0.0;                    // {2, 0, 0}
        output.access(3, 1) = -(1.0 - z)/8.0;          // {1, 1, 0}
        output.access(3, 2) =  (1.0 + y)/8.0;          // {1, 0, 1}
        output.access(3, 3) =  0.0;                    // {0, 2, 0}
        output.access(3, 4) = -(1.0 - x)/8.0;          // {0, 1, 1}
        output.access(3, 5) =  0.0;                    // {0, 0, 2}


        output.access(4, 0) =  0.0;                    // {2, 0, 0}
        output.access(4, 1) =  (1.0 + z)/8.0;          // {1, 1, 0}
        output.access(4, 2) = -(1.0 - y)/8.0;          // {1, 0, 1}
        output.access(4, 3) =  0.0;                    // {0, 2, 0}
        output.access(4, 4) = -(1.0 - x)/8.0;          // {0, 1, 1}
        output.access(4, 5) =  0.0;                    // {0, 0, 2}

        output.access(5, 0) =  0.0;                    // {2, 0, 0}
        output.access(5, 1) = -(1.0 + z)/8.0;          // {1, 1, 0}
        output.access(5, 2) =  (1.0 - y)/8.0;          // {1, 0, 1}
        output.access(5, 3) =  0.0;                    // {0, 2, 0}
        output.access(5, 4) = -(1.0 + x)/8.0;          // {0, 1, 1}
        output.access(5, 5) =  0.0;                    // {0, 0, 2}

        output.access(6, 0) =  0.0;                    // {2, 0, 0}
        output.access(6, 1) =  (1.0 + z)/8.0;          // {1, 1, 0}
        output.access(6, 2) =  (1.0 + y)/8.0;          // {1, 0, 1}
        output.access(6, 3) =  0.0;                    // {0, 2, 0}
        output.access(6, 4) =  (1.0 + x)/8.0;          // {0, 1, 1}
        output.access(6, 5) =  0.0;                    // {0, 0, 2}

        output.access(7, 0) =  0.0;                    // {2, 0, 0}
        output.access(7, 1) = -(1.0 + z)/8.0;          // {1, 1, 0}
        output.access(7, 2) = -(1.0 + y)/8.0;          // {1, 0, 1}
        output.access(7, 3) =  0.0;                    // {0, 2, 0}
        output.access(7, 4) =  (1.0 - x)/8.0;          // {0, 1, 1}
        output.access(7, 5) =  0.0;                    // {0, 0, 2}
        break;
      }
      case OPERATOR_D3:
      {
        // output is a rank-2 array with dimensions (basisCardinality_, D3Cardinality = 10)
        // (1.0 - x)*(1.0 - y)*(1.0 - z)/8.0;
        output.access(0, 0) =  0.0;      // {3, 0, 0}
        output.access(0, 1) =  0.0;      // {2, 1, 0}
        output.access(0, 2) =  0.0;      // {2, 0, 1}
        output.access(0, 3) =  0.0;      // {1, 2, 0}
        output.access(0, 4) =  -1.0/8.0; // {1, 1, 1}
        output.access(0, 5) =  0.0;      // {1, 0, 2}
        output.access(0, 6) =  0.0;      // {0, 3, 0}
        output.access(0, 7) =  0.0;      // {0, 2, 1}
        output.access(0, 8) =  0.0;      // {0, 1, 2}
        output.access(0, 9) =  0.0;      // {0, 0, 3}
        
        // (1.0 + x)*(1.0 - y)*(1.0 - z)/8.0;
        output.access(1, 0) =  0.0;      // {3, 0, 0}
        output.access(1, 1) =  0.0;      // {2, 1, 0}
        output.access(1, 2) =  0.0;      // {2, 0, 1}
        output.access(1, 3) =  0.0;      // {1, 2, 0}
        output.access(1, 4) =  1.0/8.0;  // {1, 1, 1}
        output.access(1, 5) =  0.0;      // {1, 0, 2}
        output.access(1, 6) =  0.0;      // {0, 3, 0}
        output.access(1, 7) =  0.0;      // {0, 2, 1}
        output.access(1, 8) =  0.0;      // {0, 1, 2}
        output.access(1, 9) =  0.0;      // {0, 0, 3}
        
        // (1.0 + x)*(1.0 + y)*(1.0 - z)/8.0;
        output.access(2, 0) =  0.0;      // {3, 0, 0}
        output.access(2, 1) =  0.0;      // {2, 1, 0}
        output.access(2, 2) =  0.0;      // {2, 0, 1}
        output.access(2, 3) =  0.0;      // {1, 2, 0}
        output.access(2, 4) = -1.0/8.0;  // {1, 1, 1}
        output.access(2, 5) =  0.0;      // {1, 0, 2}
        output.access(2, 6) =  0.0;      // {0, 3, 0}
        output.access(2, 7) =  0.0;      // {0, 2, 1}
        output.access(2, 8) =  0.0;      // {0, 1, 2}
        output.access(2, 9) =  0.0;      // {0, 0, 3}
        
        // (1.0 - x)*(1.0 + y)*(1.0 - z)/8.0;
        output.access(3, 0) =  0.0;      // {3, 0, 0}
        output.access(3, 1) =  0.0;      // {2, 1, 0}
        output.access(3, 2) =  0.0;      // {2, 0, 1}
        output.access(3, 3) =  0.0;      // {1, 2, 0}
        output.access(3, 4) =  1.0/8.0;  // {1, 1, 1}
        output.access(3, 5) =  0.0;      // {1, 0, 2}
        output.access(3, 6) =  0.0;      // {0, 3, 0}
        output.access(3, 7) =  0.0;      // {0, 2, 1}
        output.access(3, 8) =  0.0;      // {0, 1, 2}
        output.access(3, 9) =  0.0;      // {0, 0, 3}

        // (1.0 - x)*(1.0 - y)*(1.0 + z)/8.0;
        output.access(4, 0) =  0.0;      // {3, 0, 0}
        output.access(4, 1) =  0.0;      // {2, 1, 0}
        output.access(4, 2) =  0.0;      // {2, 0, 1}
        output.access(4, 3) =  0.0;      // {1, 2, 0}
        output.access(4, 4) =  1.0/8.0;  // {1, 1, 1}
        output.access(4, 5) =  0.0;      // {1, 0, 2}
        output.access(4, 6) =  0.0;      // {0, 3, 0}
        output.access(4, 7) =  0.0;      // {0, 2, 1}
        output.access(4, 8) =  0.0;      // {0, 1, 2}
        output.access(4, 9) =  0.0;      // {0, 0, 3}

        // (1.0 + x)*(1.0 - y)*(1.0 + z)/8.0;
        output.access(5, 0) =  0.0;      // {3, 0, 0}
        output.access(5, 1) =  0.0;      // {2, 1, 0}
        output.access(5, 2) =  0.0;      // {2, 0, 1}
        output.access(5, 3) =  0.0;      // {1, 2, 0}
        output.access(5, 4) = -1.0/8.0;  // {1, 1, 1}
        output.access(5, 5) =  0.0;      // {1, 0, 2}
        output.access(5, 6) =  0.0;      // {0, 3, 0}
        output.access(5, 7) =  0.0;      // {0, 2, 1}
        output.access(5, 8) =  0.0;      // {0, 1, 2}
        output.access(5, 9) =  0.0;      // {0, 0, 3}
        
        // (1.0 + x)*(1.0 + y)*(1.0 + z)/8.0;
        output.access(6, 0) =  0.0;      // {3, 0, 0}
        output.access(6, 1) =  0.0;      // {2, 1, 0}
        output.access(6, 2) =  0.0;      // {2, 0, 1}
        output.access(6, 3) =  0.0;      // {1, 2, 0}
        output.access(6, 4) =  1.0/8.0;  // {1, 1, 1}
        output.access(6, 5) =  0.0;      // {1, 0, 2}
        output.access(6, 6) =  0.0;      // {0, 3, 0}
        output.access(6, 7) =  0.0;      // {0, 2, 1}
        output.access(6, 8) =  0.0;      // {0, 1, 2}
        output.access(6, 9) =  0.0;      // {0, 0, 3}
        
        // (1.0 - x)*(1.0 + y)*(1.0 + z)/8.0;
        output.access(7, 0) =  0.0;      // {3, 0, 0}
        output.access(7, 1) =  0.0;      // {2, 1, 0}
        output.access(7, 2) =  0.0;      // {2, 0, 1}
        output.access(7, 3) =  0.0;      // {1, 2, 0}
        output.access(7, 4) = -1.0/8.0;  // {1, 1, 1}
        output.access(7, 5) =  0.0;      // {1, 0, 2}
        output.access(7, 6) =  0.0;      // {0, 3, 0}
        output.access(7, 7) =  0.0;      // {0, 2, 1}
        output.access(7, 8) =  0.0;      // {0, 1, 2}
        output.access(7, 9) =  0.0;      // {0, 0, 3}
        
        break;
      }
      case OPERATOR_MAX : {
        const ordinal_type jend = output.extent(1);
        const ordinal_type iend = output.extent(0);

        for (ordinal_type j=0;j<jend;++j)
          for (ordinal_type i=0;i<iend;++i)
            output.access(i, j) = 0.0;
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                  opType != OPERATOR_GRAD &&
                                  opType != OPERATOR_CURL &&
                                  opType != OPERATOR_D2 &&
                                  opType != OPERATOR_MAX,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_HEX_C1_FEM::Serial::getValues) operator is not supported");

      }
      }
    }

    template<typename SpT,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void
    Basis_HGRAD_HEX_C1_FEM::
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
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_VALUE> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      case OPERATOR_GRAD:
      case OPERATOR_D1: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_GRAD> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      case OPERATOR_CURL: {
        INTREPID2_TEST_FOR_EXCEPTION( operatorType == OPERATOR_CURL, std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_HEX_C1_FEM): CURL is invalid operator for rank-0 (scalar) functions in 3D");
        break;
      }

      case OPERATOR_DIV: {
        INTREPID2_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_HEX_C1_FEM): DIV is invalid operator for rank-0 (scalar) functions in 3D");
        break;
      }

      case OPERATOR_D2: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D2> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      case OPERATOR_D3:{
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D3> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      case OPERATOR_D4:
      case OPERATOR_D5:
      case OPERATOR_D6:
      case OPERATOR_D7:
      case OPERATOR_D8:
      case OPERATOR_D9:
      case OPERATOR_D10: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_MAX> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( !( Intrepid2::isValidOperator(operatorType) ), std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_HEX_C1_FEM): Invalid operator type");
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------

  template<typename SpT, typename OT, typename PT>
  Basis_HGRAD_HEX_C1_FEM<SpT,OT,PT>::
  Basis_HGRAD_HEX_C1_FEM() {
    this->basisCardinality_  = 8;
    this->basisDegree_       = 1;
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
    this->basisType_         = BASIS_FEM_DEFAULT;
    this->basisCoordinates_  = COORDINATES_CARTESIAN;
    this->functionSpace_     = FUNCTION_SPACE_HGRAD;

    // initialize tags
    {
      // Basis-dependent intializations
      const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

      // An array with local DoF tags assigned to basis functions, in the order of their local enumeration
      ordinal_type tags[32]  = { 0, 0, 0, 1,
                                 0, 1, 0, 1,
                                 0, 2, 0, 1,
                                 0, 3, 0, 1,
                                 0, 4, 0, 1,
                                 0, 5, 0, 1,
                                 0, 6, 0, 1,
                                 0, 7, 0, 1 };

      // host tags
      OrdinalTypeArray1DHost tagView(&tags[0], 32);

      // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
      //OrdinalTypeArray2DHost ordinalToTag;
      //OrdinalTypeArray3DHost tagToOrdinal;
      this->setOrdinalTagData(this->tagToOrdinal_,
                              this->ordinalToTag_,
                              tagView,
                              this->basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);

      //this->tagToOrdinal_ = Kokkos::create_mirror_view(typename SpT::memory_space(), tagToOrdinal);
      //Kokkos::deep_copy(this->tagToOrdinal_, tagToOrdinal);

      //this->ordinalToTag_ = Kokkos::create_mirror_view(typename SpT::memory_space(), ordinalToTag);
      //Kokkos::deep_copy(this->ordinalToTag_, ordinalToTag);
    }

    // dofCoords on host and create its mirror view to device
    Kokkos::DynRankView<typename ScalarViewType::value_type,typename SpT::array_layout,Kokkos::HostSpace>
      dofCoords("dofCoordsHost", this->basisCardinality_,this->basisCellTopology_.getDimension());
    
    dofCoords(0,0) = -1.0;   dofCoords(0,1) = -1.0; dofCoords(0,2) = -1.0;
    dofCoords(1,0) =  1.0;   dofCoords(1,1) = -1.0; dofCoords(1,2) = -1.0;
    dofCoords(2,0) =  1.0;   dofCoords(2,1) =  1.0; dofCoords(2,2) = -1.0;
    dofCoords(3,0) = -1.0;   dofCoords(3,1) =  1.0; dofCoords(3,2) = -1.0;
    dofCoords(4,0) = -1.0;   dofCoords(4,1) = -1.0; dofCoords(4,2) =  1.0;
    dofCoords(5,0) =  1.0;   dofCoords(5,1) = -1.0; dofCoords(5,2) =  1.0;
    dofCoords(6,0) =  1.0;   dofCoords(6,1) =  1.0; dofCoords(6,2) =  1.0;
    dofCoords(7,0) = -1.0;   dofCoords(7,1) =  1.0; dofCoords(7,2) =  1.0;
    
    this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);
  }

}// namespace Intrepid2

#endif

