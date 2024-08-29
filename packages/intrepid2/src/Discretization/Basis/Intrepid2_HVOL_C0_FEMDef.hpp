// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file  Intrepid2_HVOL_C0_FEMDef.hpp
    \brief Definition file FEM basis functions of degree 0 for H(vol) functions on all supported topologies.
    \author Created by Kyungjoo Kim
*/

#ifndef __INTREPID2_HVOL_C0_FEM_DEF_HPP__
#define __INTREPID2_HVOL_C0_FEM_DEF_HPP__

namespace Intrepid2 {

  namespace Impl {

    template<EOperator opType>
    template<typename OutputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HVOL_C0_FEM::Serial<opType>::
    getValues(       OutputViewType output,
               const inputViewType /* input */ ) {
      switch (opType) {
      case OPERATOR_VALUE : {
        output.access(0) = 1.0;
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
                                  opType != OPERATOR_MAX,
                                  ">>> ERROR: (Intrepid2::Basis_HVOL_C0_FEM::Serial::getValues) operator is not supported");
      }
      }
    }

    template<typename DT,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void
    Basis_HVOL_C0_FEM::
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
      case OPERATOR_CURL:
      case OPERATOR_DIV:
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
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_MAX> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( !Intrepid2::isValidOperator(operatorType), std::invalid_argument,
                                      ">>> ERROR (Basis_HVOL_C0_FEM): Invalid operator type");
      }
      }
    }
  }

  template<typename DT, typename OT, typename PT>
  Basis_HVOL_C0_FEM<DT,OT,PT>::
  Basis_HVOL_C0_FEM(const shards::CellTopology& cellTopo) {
    const ordinal_type spaceDim = cellTopo.getDimension();

    this->basisCardinality_     = 1;
    this->basisDegree_          = 0;
    this->basisCellTopologyKey_ = cellTopo.getKey();
    this->basisType_            = Intrepid2::BASIS_FEM_DEFAULT;
    this->basisCoordinates_     = Intrepid2::COORDINATES_CARTESIAN;
    this->functionSpace_        = FUNCTION_SPACE_HVOL;
    
    basisName_ = "Intrepid2_HVOL_";
    basisName_ += cellTopo.getName();
    basisName_ += "_C0_FEM";

    // initialize tags
    {
      // Basis-dependent intializations
      const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

      // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration
      ordinal_type tags[4] = { spaceDim, 0, 0, 1 };

      OrdinalTypeArray1DHost tagView(&tags[0], 4);

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
      dofCoords("dofCoordsHost", this->basisCardinality_, spaceDim);

    CellTools<Kokkos::HostSpace>::getReferenceCellCenter(Kokkos::subview(dofCoords, 0, Kokkos::ALL()), cellTopo);

    this->dofCoords_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);
  }

}
#endif
