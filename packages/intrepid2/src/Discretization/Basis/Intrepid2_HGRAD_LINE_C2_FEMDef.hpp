// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_LINE_C2_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree 1 for H(grad) functions on a Line.
 */

#ifndef __INTREPID2_HGRAD_LINE_C2_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_LINE_C2_FEM_DEF_HPP__

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------

  namespace Impl {

    template<EOperator opType>
    template<typename OutputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_LINE_C2_FEM::Serial<opType>::
    getValues(       OutputViewType output,
               const inputViewType input ) {
      switch (opType) {
      case OPERATOR_VALUE : {
        const auto x = input(0);

        output.access(0) = (x - 1.0)*x/2.0;
        output.access(1) = (1.0 + x)*x/2.0;
        output.access(2) = (1.0 + x)*(1.0 - x);
        break;
      }
      case OPERATOR_GRAD : {
        const auto x = input(0);
        
        output.access(0, 0) = x-0.5;
        output.access(1, 0) = x+0.5;
        output.access(2, 0) = -2.0*x;
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
                                  opType != OPERATOR_MAX,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_C2_FEM::Serial::getValues) operator is not supported");

      }
      }
    }

    template<typename DT,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void
    Basis_HGRAD_LINE_C2_FEM::
    getValues( const typename DT::execution_space& space,
                     Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const EOperator operatorType ) {
      typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
      typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
      typedef typename ExecSpace<typename inputPointViewType::execution_space,typename DT::execution_space>::ExecSpaceType ExecSpaceType;

      // Number of evaluation points = dim 0 of inputPoints
      const auto loopSize = inputPoints.extent(0);
      Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(space, 0, loopSize);

      switch (operatorType) {

      case OPERATOR_VALUE: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_VALUE> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      case OPERATOR_GRAD:
      case OPERATOR_DIV:
      case OPERATOR_CURL:
      case OPERATOR_D1: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_GRAD> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      case OPERATOR_D2:
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
        INTREPID2_TEST_FOR_EXCEPTION( !Intrepid2::isValidOperator(operatorType), std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_LINE_C2_FEM): Invalid operator type");
      }
      }
    }



  }

  // -------------------------------------------------------------------------------------

  template<typename DT, typename OT, typename PT>
  Basis_HGRAD_LINE_C2_FEM<DT,OT,PT>::
  Basis_HGRAD_LINE_C2_FEM() {
    const ordinal_type spaceDim = 1;
    this->basisCardinality_     = 3;
    this->basisDegree_          = 2;
    this->basisCellTopologyKey_ = shards::Line<2>::key;
    this->basisType_            = BASIS_FEM_DEFAULT;
    this->basisCoordinates_     = COORDINATES_CARTESIAN;
    this->functionSpace_        = FUNCTION_SPACE_HGRAD;

    // initialize tags
    {
      // Basis-dependent intializations
      const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

      // An array with local DoF tags assigned to basis functions, in the order of their local enumeration
      ordinal_type tags[12]  = { 0, 0, 0, 1,
                                0, 1, 0, 1,
                                1, 0, 0, 1 };

      // host tags
      OrdinalTypeArray1DHost tagView(&tags[0], 12);

      // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
      // tags are constructed on host and sent to devices
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

    dofCoords(0,0) = -1.0;
    dofCoords(1,0) =  1.0;
    dofCoords(2,0) =  0.0;

    this->dofCoords_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);
  }

}

#endif
