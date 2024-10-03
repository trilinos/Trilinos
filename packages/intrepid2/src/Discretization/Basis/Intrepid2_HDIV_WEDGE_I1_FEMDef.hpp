// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HDIV_WEDGE_I1_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree 1 for H(div) functions on WEDGE cells.
    \author Created by P. Bochev, D. Ridzal and K. Peterson.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HDIV_WEDGE_I1_FEM_DEF_HPP__
#define __INTREPID2_HDIV_WEDGE_I1_FEM_DEF_HPP__

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------
  namespace Impl {

    template<EOperator opType>
    template<typename OutputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HDIV_WEDGE_I1_FEM::Serial<opType>::
    getValues(       OutputViewType output,
               const inputViewType input ) {
      switch (opType) {
      case OPERATOR_VALUE: {
        const auto x = input(0);
        const auto y = input(1);
        const auto z = input(2);

        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        output.access(0, 0) = x*2.0;
        output.access(0, 1) = (y - 1.0)*2.0;
        output.access(0, 2) = 0.0;

        output.access(1, 0) = x*2.0;
        output.access(1, 1) = y*2.0;
        output.access(1, 2) = 0.0;

        output.access(2, 0) = (x - 1.0)*2.0;
        output.access(2, 1) = y*2.0;
        output.access(2, 2) = 0.0;

        output.access(3, 0) = 0.0;
        output.access(3, 1) = 0.0;
        output.access(3, 2) = (z - 1.0)/2.0;

        output.access(4, 0) = 0.0;
        output.access(4, 1) = 0.0;
        output.access(4, 2) = (1.0 + z)/2.0;
        break;
      }
      case OPERATOR_DIV: {

        // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
        output.access(0) = 4.0;
        output.access(1) = 4.0;
        output.access(2) = 4.0;
        output.access(3) = 0.5;
        output.access(4) = 0.5;
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                  opType != OPERATOR_DIV,
                                  ">>> ERROR: (Intrepid2::Basis_HDIV_WEDGE_I1_FEM::Serial::getValues) operator is not supported");

      }
      }
    }

    template<typename DT,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void
    Basis_HDIV_WEDGE_I1_FEM::
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
      case OPERATOR_DIV: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_DIV> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( !( Intrepid2::isValidOperator(operatorType) ), std::invalid_argument,
                                      ">>> ERROR (Basis_HDIV_WEDGE_I1_FEM): Invalid operator type");
      }
      }
    }


  }
  // -------------------------------------------------------------------------------------

  template<typename DT, typename OT, typename PT>
  Basis_HDIV_WEDGE_I1_FEM<DT,OT,PT>::
  Basis_HDIV_WEDGE_I1_FEM() {
    const ordinal_type spaceDim = 3;
    this->basisCardinality_     = 5;
    this->basisDegree_          = 1;
    this->basisCellTopologyKey_ = shards::Wedge<6>::key;
    this->basisType_            = BASIS_FEM_DEFAULT;
    this->basisCoordinates_     = COORDINATES_CARTESIAN;
    this->functionSpace_        = FUNCTION_SPACE_HDIV;

    // initialize tags
    {
      // Basis-dependent intializations
      const ordinal_type tagSize  = 4;        // size of DoF tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

      // An array with local DoF tags assigned to basis functions, in the order of their local enumeration
      ordinal_type tags[20]  = { 2, 0, 0, 1,
                                 2, 1, 0, 1,
                                 2, 2, 0, 1,
                                 2, 3, 0, 1,
                                 2, 4, 0, 1 };

      // host tags
      OrdinalTypeArray1DHost tagView(&tags[0], 20);

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
    }

    // dofCoords on host and create its mirror view to device
    Kokkos::DynRankView<typename ScalarViewType::value_type,typename DT::execution_space::array_layout,Kokkos::HostSpace>
      dofCoords("dofCoordsHost", this->basisCardinality_,spaceDim);

    dofCoords(0,0) =  0.5;      dofCoords(0,1) =  0.0;      dofCoords(0,2) =  0.0;
    dofCoords(1,0) =  0.5;      dofCoords(1,1) =  0.5;      dofCoords(1,2) =  0.0;
    dofCoords(2,0) =  0.0;      dofCoords(2,1) =  0.5;      dofCoords(2,2) =  0.0;
    dofCoords(3,0) =  1.0/3.0;  dofCoords(3,1) =  1.0/3.0;  dofCoords(3,2) = -1.0;
    dofCoords(4,0) =  1.0/3.0;  dofCoords(4,1) =  1.0/3.0;  dofCoords(4,2) =  1.0;

    this->dofCoords_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);

    // dofCoords on host and create its mirror view to device
    Kokkos::DynRankView<typename ScalarViewType::value_type,typename DT::execution_space::array_layout,Kokkos::HostSpace>
      dofCoeffs("dofCoeffsHost", this->basisCardinality_,spaceDim);

    // dofCoeffs are normals to edges
    dofCoeffs(0,0) =  0.0;   dofCoeffs(0,1) = -0.5;   dofCoeffs(0,2) =  0.0;
    dofCoeffs(1,0) =  0.5;   dofCoeffs(1,1) =  0.5;   dofCoeffs(1,2) =  0.0;
    dofCoeffs(2,0) = -0.5;   dofCoeffs(2,1) =  0.0;   dofCoeffs(2,2) =  0.0;
    dofCoeffs(3,0) =  0.0;   dofCoeffs(3,1) =  0.0;   dofCoeffs(3,2) = -1.0;
    dofCoeffs(4,0) =  0.0;   dofCoeffs(4,1) =  0.0;   dofCoeffs(4,2) =  1.0;

    this->dofCoeffs_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoeffs);
    Kokkos::deep_copy(this->dofCoeffs_, dofCoeffs);
  }

}// namespace Intrepid2
#endif
