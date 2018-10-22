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

/** \file   Intrepid2_HGRAD_PYR_C1_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree 1 for H(grad) functions on PYR cells.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_PYR_C1_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_PYR_C1_FEM_DEF_HPP__

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------

  namespace Impl {

    template<EOperator opType>
    template<typename outputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_PYR_C1_FEM::Serial<opType>::
    getValues(       outputViewType output,
               const inputViewType input ) {
      const auto eps = epsilon();

      static_assert(std::is_same<
                    typename outputViewType::value_type,
                    typename inputViewType::value_type>::value,"Input/output view has different value types");

      typedef typename outputViewType::value_type value_type;

      const value_type x = input(0);
      const value_type y = input(1);
      const value_type ztmp = input(2);

      //be sure that the basis functions are defined when z is very close to 1.
      const value_type z = ( (value_type(1.0) - ztmp) < value_type(eps) ? value_type(1.0 - eps) : ztmp );

      switch (opType) {

      case OPERATOR_VALUE: {
        const value_type factor = 0.25/(1.0 - z);

        // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
        output.access(0) = (1.0 - x - z) * (1.0 - y - z) * factor;
        output.access(1) = (1.0 + x - z) * (1.0 - y - z) * factor;
        output.access(2) = (1.0 + x - z) * (1.0 + y - z) * factor;
        output.access(3) = (1.0 - x - z) * (1.0 + y - z) * factor;
        output.access(4) = z;
        break;
      }
      case OPERATOR_GRAD: {
        const value_type factor  = 0.25/(1.0 - z);
        const value_type factor2 = 4.0 * factor * factor;

        // output.accessValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        output.access(0, 0) = (y + z - 1.0) * factor;
        output.access(0, 1) = (x + z - 1.0) * factor;
        output.access(0, 2) = x * y * factor2 - 0.25;

        output.access(1, 0) =  (1.0 - y - z) * factor;
        output.access(1, 1) =  (z - x - 1.0) * factor;
        output.access(1, 2) =  - x*y * factor2 - 0.25;

        output.access(2, 0) =  (1.0 + y - z) * factor;
        output.access(2, 1) =  (1.0 + x - z) * factor;
        output.access(2, 2) =  x * y * factor2 - 0.25;

        output.access(3, 0) =  (z - y - 1.0) * factor;
        output.access(3, 1) =  (1.0 - x - z) * factor;
        output.access(3, 2) =  - x*y * factor2 - 0.25;

        output.access(4, 0) =  0.0;
        output.access(4, 1) =  0.0;
        output.access(4, 2) =  1;
        break;
      }
      case OPERATOR_D2: {
        const value_type factor  = 0.25/(1.0 - z);
        const value_type factor2 = 4.0 * factor * factor;
        const value_type factor3 = 8.0 * factor * factor2;

        // output.accessValues is a rank-3 array with dimensions (basisCardinality_, dim0, D2Cardinality = 6)
        output.access(0, 0) =  0.0;                    // {2, 0, 0}
        output.access(0, 1) =  factor;          	      // {1, 1, 0}
        output.access(0, 2) =  y*factor2;              // {1, 0, 1}
        output.access(0, 3) =  0.0;                    // {0, 2, 0}
        output.access(0, 4) =  x*factor2;              // {0, 1, 1}
        output.access(0, 5) =  x*y*factor3;            // {0, 0, 2}

        output.access(1, 0) =  0.0;                    // {2, 0, 0}
        output.access(1, 1) =  -factor;	              // {1, 1, 0}
        output.access(1, 2) =  -y*factor2; 	      // {1, 0, 1}
        output.access(1, 3) =  0.0;                    // {0, 2, 0}
        output.access(1, 4) =  -x*factor2;             // {0, 1, 1}
        output.access(1, 5) =  -x*y*factor3;           // {0, 0, 2}

        output.access(2, 0) =  0.0;                    // {2, 0, 0}
        output.access(2, 1) =  factor;          	      // {1, 1, 0}
        output.access(2, 2) =  y*factor2;              // {1, 0, 1}
        output.access(2, 3) =  0.0;                    // {0, 2, 0}
        output.access(2, 4) =  x*factor2;       	      // {0, 1, 1}
        output.access(2, 5) =  x*y*factor3;            // {0, 0, 2}

        output.access(3, 0) =  0.0;                    // {2, 0, 0}
        output.access(3, 1) =  -factor;	              // {1, 1, 0}
        output.access(3, 2) =  -y*factor2;             // {1, 0, 1}
        output.access(3, 3) =  0.0;                    // {0, 2, 0}
        output.access(3, 4) =  -x*factor2;	      // {0, 1, 1}
        output.access(3, 5) =  -x*y*factor3;           // {0, 0, 2}

        output.access(4, 0) =  0.0;                    // {2, 0, 0}
        output.access(4, 1) =  0.0;          	      // {1, 1, 0}
        output.access(4, 2) =  0.0;          	      // {1, 0, 1}
        output.access(4, 3) =  0.0;                    // {0, 2, 0}
        output.access(4, 4) =  0.0;                    // {0, 1, 1}
        output.access(4, 5) =  0.0;                    // {0, 0, 2}
        break;
      }
      case OPERATOR_MAX: {
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
                                  opType != OPERATOR_D2 &&
                                  opType != OPERATOR_MAX,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_PYR_C1_FEM::Serial::getValues) operator is not supported");
      }
      }
    }

    template<typename SpT,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void
    Basis_HGRAD_PYR_C1_FEM::
    getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const EOperator operatorType )  {
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
                                      ">>> ERROR (Basis_HGRAD_PYR_C1_FEM): CURL is invalid operator for rank-0 (scalar) functions in 3D");
        break;
      }
      case OPERATOR_DIV: {
        INTREPID2_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_PYR_C1_FEM): DIV is invalid operator for rank-0 (scalar) functions in 3D");
        break;
      }
      case OPERATOR_D2: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D2> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      case OPERATOR_D3:
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
                                      ">>> ERROR (Basis_HGRAD_PYR_C1_FEM): Invalid operator type");
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------

  template<typename SpT, typename OT, typename PT>
  Basis_HGRAD_PYR_C1_FEM<SpT,OT,PT>::
  Basis_HGRAD_PYR_C1_FEM() {
    this->basisCardinality_  = 5;
    this->basisDegree_       = 1;
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<5> >() );
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
      ordinal_type tags[20]  = { 0, 0, 0, 1,
                                 0, 1, 0, 1,
                                 0, 2, 0, 1,
                                 0, 3, 0, 1,
                                 0, 4, 0, 1 };


      // host tags
      ordinal_type_array_1d_host tagView(&tags[0], 20);

      // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
      //ordinal_type_array_2d_host ordinalToTag;
      //ordinal_type_array_3d_host tagToOrdinal;
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
    Kokkos::DynRankView<typename scalarViewType::value_type,typename SpT::array_layout,Kokkos::HostSpace>
      dofCoords("dofCoordsHost", this->basisCardinality_,this->basisCellTopology_.getDimension());

    dofCoords(0,0) = -1.0;  dofCoords(0,1) = -1.0;  dofCoords(0,2) =  0.0;
    dofCoords(1,0) =  1.0;  dofCoords(1,1) = -1.0;  dofCoords(1,2) =  0.0;
    dofCoords(2,0) =  1.0;  dofCoords(2,1) =  1.0;  dofCoords(2,2) =  0.0;
    dofCoords(3,0) = -1.0;  dofCoords(3,1) =  1.0;  dofCoords(3,2) =  0.0;
    dofCoords(4,0) =  0.0;  dofCoords(4,1) =  0.0;  dofCoords(4,2) =  1.0;

    this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);
  }

}

#endif
