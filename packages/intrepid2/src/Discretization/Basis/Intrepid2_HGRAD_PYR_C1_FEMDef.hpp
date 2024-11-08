// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    template<typename OutputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_PYR_C1_FEM::Serial<opType>::
    getValues(       OutputViewType output,
               const inputViewType input ) {
      const auto eps = epsilon();

      static_assert(std::is_same<
                    typename OutputViewType::value_type,
                    typename inputViewType::value_type>::value,"Input/output view has different value types");

      typedef typename OutputViewType::value_type value_type;

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
      default: {
        INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                  opType != OPERATOR_GRAD &&
                                  opType != OPERATOR_D2,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_PYR_C1_FEM::Serial::getValues) operator is not supported");
      }
      }
    }

    template<typename DT,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void
    Basis_HGRAD_PYR_C1_FEM::
    getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const EOperator operatorType )  {
      typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
      typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
      typedef typename ExecSpace<typename inputPointViewType::execution_space,typename DT::execution_space>::ExecSpaceType ExecSpaceType;

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
        INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_PYR_C1_FEM): Operator not implemented yet");
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

  template<typename DT, typename OT, typename PT>
  Basis_HGRAD_PYR_C1_FEM<DT,OT,PT>::
  Basis_HGRAD_PYR_C1_FEM() {
    const ordinal_type spaceDim = 3;
    this->basisCardinality_     = 5;
    this->basisDegree_          = 1;
    this->basisCellTopologyKey_ = shards::Pyramid<5>::key;
    this->basisType_            = BASIS_FEM_DEFAULT;
    this->basisCoordinates_     = COORDINATES_CARTESIAN;
    this->functionSpace_        = FUNCTION_SPACE_HGRAD;

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
      OrdinalTypeArray1DHost tagView(&tags[0], 20);

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
    Kokkos::DynRankView<typename ScalarViewType::value_type,typename DT::execution_space::array_layout,Kokkos::HostSpace>
      dofCoords("dofCoordsHost", this->basisCardinality_,spaceDim);

    dofCoords(0,0) = -1.0;  dofCoords(0,1) = -1.0;  dofCoords(0,2) =  0.0;
    dofCoords(1,0) =  1.0;  dofCoords(1,1) = -1.0;  dofCoords(1,2) =  0.0;
    dofCoords(2,0) =  1.0;  dofCoords(2,1) =  1.0;  dofCoords(2,2) =  0.0;
    dofCoords(3,0) = -1.0;  dofCoords(3,1) =  1.0;  dofCoords(3,2) =  0.0;
    dofCoords(4,0) =  0.0;  dofCoords(4,1) =  0.0;  dofCoords(4,2) =  1.0;

    this->dofCoords_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);
  }

  template<typename DT, typename OT, typename PT>
  void 
  Basis_HGRAD_PYR_C1_FEM<DT,OT,PT>::getScratchSpaceSize(       
                                    ordinal_type& perTeamSpaceSize,
                                    ordinal_type& perThreadSpaceSize,
                              const PointViewType inputPoints,
                              const EOperator operatorType) const {
    perTeamSpaceSize = 0;
    perThreadSpaceSize = 0;
  }

  template<typename DT, typename OT, typename PT>
  KOKKOS_INLINE_FUNCTION
  void 
  Basis_HGRAD_PYR_C1_FEM<DT,OT,PT>::getValues(       
          OutputViewType outputValues,
      const PointViewType  inputPoints,
      const EOperator operatorType,
      const typename Kokkos::TeamPolicy<typename DT::execution_space>::member_type& team_member,
      const typename DT::execution_space::scratch_memory_space & scratchStorage, 
      const ordinal_type subcellDim,
      const ordinal_type subcellOrdinal) const {

      INTREPID2_TEST_FOR_ABORT( !((subcellDim <= 0) && (subcellOrdinal == -1)),
        ">>> ERROR: (Intrepid2::Basis_HGRAD_PYR_C1_FEM::getValues), The capability of selecting subsets of basis functions has not been implemented yet.");

      (void) scratchStorage; //avoid unused variable warning

      const int numPoints = inputPoints.extent(0);

      switch(operatorType) {
        case OPERATOR_VALUE:
          Kokkos::parallel_for (Kokkos::TeamThreadRange (team_member, numPoints), [=] (ordinal_type& pt) {
            auto       output = Kokkos::subview( outputValues, Kokkos::ALL(), pt, Kokkos::ALL() );
            const auto input  = Kokkos::subview( inputPoints,                 pt, Kokkos::ALL() );
            Impl::Basis_HGRAD_PYR_C1_FEM::Serial<OPERATOR_VALUE>::getValues( output, input);
          });
          break;
        case OPERATOR_GRAD:
          Kokkos::parallel_for (Kokkos::TeamThreadRange (team_member, numPoints), [=] (ordinal_type& pt) {
            auto       output = Kokkos::subview( outputValues, Kokkos::ALL(), pt, Kokkos::ALL() );
            const auto input  = Kokkos::subview( inputPoints,                 pt, Kokkos::ALL() );
            Impl::Basis_HGRAD_PYR_C1_FEM::Serial<OPERATOR_GRAD>::getValues( output, input);
          });
          break;
        default: {}
    }
  }
  
}// namespace Intrepid2
#endif
