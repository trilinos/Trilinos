// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_PYR_I2_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree 1 for H(grad) functions on PYR cells.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_PYR_I2_SERENDIPITY_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_PYR_I2_SERENDIPITY_FEM_DEF_HPP__

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------

  namespace Impl {

    template<EOperator opType>
    template<typename OutputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_PYR_I2_FEM::Serial<opType>::
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
        const value_type w = 1.0/(1.0 - z);

        output.access(0) = 0.25 * (-x - y - 1.0)*((1.0-x)*(1.0-y) - z + x*y*z*w);
        output.access(1) = 0.25 * ( x - y - 1.0)*((1.0+x)*(1.0-y) - z - x*y*z*w);
        output.access(2) = 0.25 * ( x + y - 1.0)*((1.0+x)*(1.0+y) - z + x*y*z*w);
        output.access(3) = 0.25 * (-x + y - 1.0)*((1.0-x)*(1.0+y) - z - x*y*z*w);

        output.access(4) =  z * (2.0*z - 1.0);

        output.access(5) = 0.5 * (1.0 + x - z)*(1.0 - x - z)*(1.0 - y - z)*w;
        output.access(6) = 0.5 * (1.0 + y - z)*(1.0 - y - z)*(1.0 + x - z)*w;
        output.access(7) = 0.5 * (1.0 + x - z)*(1.0 - x - z)*(1.0 + y - z)*w;
        output.access(8) = 0.5 * (1.0 + y - z)*(1.0 - y - z)*(1.0 - x - z)*w;

        output.access(9) = z*(1.0 - x - z)*(1.0 - y - z)*w;
        output.access(10) = z*(1.0 + x - z)*(1.0 - y - z)*w;
        output.access(11) = z*(1.0 + x - z)*(1.0 + y - z)*w;
        output.access(12) = z*(1.0 - x - z)*(1.0 + y - z)*w;

        break;
      }
      case OPERATOR_GRAD:
      case OPERATOR_D1: {
        const value_type w  = 1.0/(1.0 - z);

        output.access(0, 0) =  0.25*(-1.0-x-y)*(-1.0+y + y*z*w) - 0.25*((1.0-x)*(1.0-y)-z + x*y*z*w);
        output.access(0, 1) =  0.25*(-1.0-x-y)*(-1.0+x + x*z*w) - 0.25*((1.0-x)*(1.0-y)-z + x*y*z*w);
        output.access(0, 2) =  0.25*(-1.0-x-y)*(-1.0 + x*y*w + x*y*z*w*w);
        
        output.access(1, 0) =  0.25*(-1.0+x-y)*( 1.0-y - y*z*w) + 0.25*((1.0+x)*(1.0-y)-z - x*y*z*w);
        output.access(1, 1) =  0.25*(-1.0+x-y)*(-1.0-x - x*z*w) - 0.25*((1.0+x)*(1.0-y)-z - x*y*z*w);
        output.access(1, 2) =  0.25*(-1.0+x-y)*(-1.0 - x*y*w - x*y*z*w*w);
        
        output.access(2, 0) =  0.25*(-1.0+x+y)*(1.0+y + y*z*w) + 0.25*((1.0+x)*(1.0+y)-z + x*y*z*w);
        output.access(2, 1) =  0.25*(-1.0+x+y)*(1.0+x + x*z*w) + 0.25*((1.0+x)*(1.0+y)-z + x*y*z*w);
        output.access(2, 2) =  0.25*(-1.0+x+y)*(-1.0 + x*y*w + x*y*z*w*w);
        
        output.access(3, 0) =  0.25*(-1.0-x+y)*(-1.0-y - y*z*w) - 0.25*((1.0-x)*(1.0+y)-z - x*y*z*w);
        output.access(3, 1) =  0.25*(-1.0-x+y)*( 1.0-x - x*z*w) + 0.25*((1.0-x)*(1.0+y)-z - x*y*z*w);
        output.access(3, 2) =  0.25*(-1.0-x+y)*(-1.0 - x*y*w - x*y*z*w*w);
        
        output.access(4, 0) =  0.0;
        output.access(4, 1) =  0.0;
        output.access(4, 2) =  -1.0 + 4.0*z;
        
        output.access(5, 0) = -x*w*(1.0-y-z);
        output.access(5, 1) = -0.5*(1.0-x-z)*(1.0+x-z)*w;
        output.access(5, 2) =  0.5*y*x*x*w*w + 0.5*y - 1.0+z;
        
        output.access(6, 0) =  0.5*(1.0-y-z)*(1.0+y-z)*w;
        output.access(6, 1) = -y*w*(1.0+x-z);
        output.access(6, 2) = -0.5*x*y*y*w*w - 0.5*x - 1.0+z; 
        
        output.access(7, 0) = -x*w*(1.0+y-z);
        output.access(7, 1) =  0.5*(1.0-x-z)*(1.0+x-z)*w;
        output.access(7, 2) = -0.5*y*x*x*w*w - 0.5*y - 1.0+z;
        
        output.access(8, 0) = -0.5*(1.0-y-z)*(1.0+y-z)*w;
        output.access(8, 1) = -y*w*(1.0-x-z);
        output.access(8, 2) =  0.5*x*y*y*w*w + 0.5*x - 1.0+z;
        
        output.access(9, 0) = -(1.0-y-z)*z*w;
        output.access(9, 1) = -(1.0-x-z)*z*w;
        output.access(9, 2) =  x*y*w*w + 1.0 - x - y - 2.0*z;
        
        output.access(10,0) =  (1.0-y-z)*z*w;
        output.access(10,1) = -(1.0+x-z)*z*w;
        output.access(10,2) = -x*y*w*w + 1.0 + x - y - 2.0*z;
        
        output.access(11,0) =  (1.0+y-z)*z*w;
        output.access(11,1) =  (1.0+x-z)*z*w;
        output.access(11,2) =  x*y*w*w + 1.0 + x + y - 2.0*z;
        
        output.access(12,0) = -(1.0+y-z)*z*w;
        output.access(12,1) =  (1.0-x-z)*z*w;
        output.access(12,2) = -x*y*w*w + 1.0 - x + y - 2.0*z;

        break;
      }
      case OPERATOR_D2: {
        const value_type w  = 1.0/(1.0 - z);
        
        output.access(0, 0) = -0.5*(-1.0+y+y*z*w);
        output.access(0, 1) = -(-0.25 + 0.5*x + 0.5*y + 0.5*z)*w;
        output.access(0, 2) =  0.25 + (-0.25*y-0.5*x*y-0.25*y*y)*w*w; 
        output.access(0, 3) = -0.5*(-1.0+x+x*z*w);
        output.access(0, 4) =  0.25 + (-0.25*x*x-0.25*x-0.5*x*y)*w*w; 
        output.access(0, 5) =  0.5*x*y*(-1.0-x-y)*w*w*w;
        
        output.access(1, 0) =  0.5*(1.0-y-y*z*w); 
        output.access(1, 1) =-(0.25 + 0.5*x - 0.5*y - 0.5*z)*w;
        output.access(1, 2) = -0.25 + (0.25*y-0.5*x*y+0.25*y*y)*w*w;
        output.access(1, 3) = -0.5*(-1.0-x-x*z*w); 
        output.access(1, 4) =  0.25 + (-0.25*x*x + 0.25*x + 0.5*x*y)*w*w; 
        output.access(1, 5) = -0.5*x*y*(-1.0+x-y)*w*w*w; 
        
        output.access(2, 0) =  0.5*(1.0+y+y*z*w);
        output.access(2, 1) =-(-0.25 - 0.5*x - 0.5*y + 0.5*z)*w; 
        output.access(2, 2) = -0.25 + (-0.25*y+0.5*x*y+0.25*y*y)*w*w;
        output.access(2, 3) =  0.5*(1.0+x+x*z*w); 
        output.access(2, 4) = -0.25 + (0.25*x*x -0.25*x + 0.5*x*y)*w*w;  
        output.access(2, 5) =  0.5*x*y*(-1.0+x+y)*w*w*w;
        
        output.access(3, 0) = -0.5*(-1.0-y-y*z*w);
        output.access(3, 1) =-(0.25 - 0.5*x + 0.5*y - 0.5*z)*w; 
        output.access(3, 2) =  0.25 + (0.25*y+0.5*x*y-0.25*y*y)*w*w; 
        output.access(3, 3) =  0.5*(1.0-x-x*z*w);
        output.access(3, 4) = -0.25 + (0.25*x + 0.25*x*x - 0.5*x*y)*w*w;
        output.access(3, 5) = -0.5*x*y*(-1.0-x+y)*w*w*w;
        
        output.access(4, 0) =  0.0;
        output.access(4, 1) =  0.0;
        output.access(4, 2) =  0.0;
        output.access(4, 3) =  0.0; 
        output.access(4, 4) =  0.0;
        output.access(4, 5) =  4.0; 
        
        output.access(5, 0) = -(1.0-y-z)*w;
        output.access(5, 1) =  x*w; 
        output.access(5, 2) =  x*y*w*w;
        output.access(5, 3) =  0.0; 
        output.access(5, 4) =  0.5*x*x*w*w + 0.5; 
        output.access(5, 5) =  x*x*y*w*w*w + 1.0;
        
        output.access(6, 0) =  0.0;
        output.access(6, 1) = -y*w;
        output.access(6, 2) = -0.5*y*y*w*w - 0.5;
        output.access(6, 3) =-(1.0+x-z)*w; 
        output.access(6, 4) = -x*y*w*w; 
        output.access(6, 5) = -x*y*y*w*w*w + 1.0; 
        
        output.access(7, 0) = -(1.0+y-z)*w;
        output.access(7, 1) = -x*w;
        output.access(7, 2) = -x*y*w*w; 
        output.access(7, 3) =  0.0; 
        output.access(7, 4) = -0.5*x*x*w*w - 0.5;
        output.access(7, 5) = -x*x*y*w*w*w + 1.0; 
        
        output.access(8, 0) =  0.0;
        output.access(8, 1) =  y*w;
        output.access(8, 2) =  0.5*y*y*w*w + 0.5; 
        output.access(8, 3) = -(1.0-x-z)*w; 
        output.access(8, 4) =  x*y*w*w;
        output.access(8, 5) =  x*y*y*w*w*w + 1.0;

        output.access(9, 0) =  0.0;
        output.access(9, 1) =  z*w; 
        output.access(9, 2) =  y*w*w - 1.0;
        output.access(9, 3) =  0.0; 
        output.access(9, 4) =  x*w*w - 1.0; 
        output.access(9, 5) =  2.0*x*y*w*w*w - 2.0; 
         
        output.access(10,0) =  0.0;
        output.access(10,1) = -z*w;
        output.access(10,2) = -y*w*w + 1.0;
        output.access(10,3) =  0.0;
        output.access(10,4) = -x*w*w - 1.0;
        output.access(10,5) = -2.0*x*y*w*w*w - 2.0;
        
        output.access(11,0) =  0.0;
        output.access(11,1) =  z*w; 
        output.access(11,2) =  y*w*w + 1.0;
        output.access(11,3) =  0.0;
        output.access(11,4) =  x*w*w + 1.0; 
        output.access(11,5) =  2.0*x*y*w*w*w - 2.0;      
        
        output.access(12,0) =  0.0;
        output.access(12,1) = -z*w; 
        output.access(12,2) = -y*w*w - 1.0; 
        output.access(12,3) =  0.0;
        output.access(12,4) = -x*w*w + 1.0;    
        output.access(12,5) = -2.0*x*y*w*w*w - 2.0;

        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                  opType != OPERATOR_GRAD &&
                                  opType != OPERATOR_D1 &&
                                  opType != OPERATOR_D2,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_PYR_I2_FEM::Serial::getValues) operator is not supported");
      }
      }
    }

    template<typename DT,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void
    Basis_HGRAD_PYR_I2_FEM::
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
                                      ">>> ERROR (Basis_HGRAD_PYR_I2_FEM): CURL is invalid operator for rank-0 (scalar) functions in 3D");
        break;
      }
      case OPERATOR_DIV: {
        INTREPID2_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_PYR_I2_FEM): DIV is invalid operator for rank-0 (scalar) functions in 3D");
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
                                      ">>> ERROR (Basis_HGRAD_PYR_I2_FEM): Operator not implemented yet");
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( !( Intrepid2::isValidOperator(operatorType) ), std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_PYR_I2_FEM): Invalid operator type");
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------

  template<typename DT, typename OT, typename PT>
  Basis_HGRAD_PYR_I2_FEM<DT,OT,PT>::
  Basis_HGRAD_PYR_I2_FEM() {
    const ordinal_type spaceDim = 3;
    this->basisCardinality_     = 13;
    this->basisDegree_          = 2;
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
      ordinal_type tags[52]  = { 0, 0, 0, 1,
                                 0, 1, 0, 1,
                                 0, 2, 0, 1,
                                 0, 3, 0, 1,
                                 0, 4, 0, 1,
                                 1, 0, 0, 1,
                                 1, 1, 0, 1,
                                 1, 2, 0, 1,
                                 1, 3, 0, 1,
                                 1, 4, 0, 1,
                                 1, 5, 0, 1,
                                 1, 6, 0, 1,
                                 1, 7, 0, 1 };


      // host tags
      OrdinalTypeArray1DHost tagView(&tags[0], 52);

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

    dofCoords(5,0) =   0.0;  dofCoords(5,1) =  -1.0;  dofCoords(5,2) =   0.0;
    dofCoords(6,0) =   1.0;  dofCoords(6,1) =   0.0;  dofCoords(6,2) =   0.0;
    dofCoords(7,0) =   0.0;  dofCoords(7,1) =   1.0;  dofCoords(7,2) =   0.0;
    dofCoords(8,0) =  -1.0;  dofCoords(8,1) =   0.0;  dofCoords(8,2) =   0.0;
    dofCoords(9,0) =  -0.5;  dofCoords(9,1) =  -0.5;  dofCoords(9,2) =   0.5;
    dofCoords(10,0) =  0.5;  dofCoords(10,1) = -0.5;  dofCoords(10,2) =  0.5;
    dofCoords(11,0) =  0.5;  dofCoords(11,1) =  0.5;  dofCoords(11,2) =  0.5;
    dofCoords(12,0) = -0.5;  dofCoords(12,1) =  0.5;  dofCoords(12,2) =  0.5;

    this->dofCoords_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);
  }

}

#endif
