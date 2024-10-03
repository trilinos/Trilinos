// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_TET_COMP12_FEMDef.hpp
    \brief  Definition file for the composite H(grad)-compatible FEM basis
            of degree 1 on Tetrahedron cell with 12 sub-tetrahedrons.
    \author Created by P. Bochev, J. Ostien, K. Peterson and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_TET_COMP12_FEMDEF_HPP__
#define __INTREPID2_HGRAD_TET_COMP12_FEMDEF_HPP__

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------
  namespace Impl {

    // assume that subtets are disjoint and a single point belong to one subtet
    // at the interface, it returns the first one that satisfy the condition
    template<typename pointValueType>
    KOKKOS_INLINE_FUNCTION
    ordinal_type
    Basis_HGRAD_TET_COMP12_FEM::getLocalSubTet( const pointValueType x, 
                                                const pointValueType y, 
                                                const pointValueType z ) {
      
      const pointValueType 
        xyz = x + y + z,
        xy = x + y,
        xz = x + z,
        yz = y + z;
      
      // cycle through each subdomain and push back if the point lies within
      
      // subtet #0 E0 := 0.0 <= r + s + t <= 0.5 && 0.0 <= r <= 0.5 && 0.0 <= s <= 0.5 && 0.0 <= t <= 0.5
      if ( (0.0 <= xyz && xyz <= 0.5) && 
           (0.0 <= x   && x   <= 0.5) && 
           (0.0 <= y   && y   <= 0.5) && 
           (0.0 <= z   && z   <= 0.5) ) 
        return 0;

      // subtet #1 E1 := 0.5 <= r + s + t <= 1.0 && 0.5 <= r <= 1.0 && 0.0 <= s <= 0.5 && 0.0 <= t <= 0.5
      if ( (0.5 <= xyz && xyz <= 1.0) && 
           (0.5 <= x   && x   <= 1.0) && 
           (0.0 <= y   && y   <= 0.5) && 
           (0.0 <= z   && z   <= 0.5) ) 
        return 1;
      
      // subtet #2 E2 := 0.5 <= r + s + t <= 1.0 && 0.0 <= r <= 0.5 && 0.5 <= s <= 1.0 && 0.0 <= t <= 0.5
      if ( (0.5 <= xyz && xyz <= 1.0) && 
           (0.0 <= x   && x   <= 0.5) && 
           (0.5 <= y   && y   <= 1.0) && 
           (0.0 <= z   && z   <= 0.5) ) 
        return 2;
      
      // subtet #3 E3 := 0.5 <= r + s + t <= 1.0 && 0.0 <= r <= 0.5 && 0.0 <= s <= 0.5 && 0.5 <= t <= 1.0
      if ( (0.5 <= xyz && xyz <= 1.0) && 
           (0.0 <= x   && x   <= 0.5) && 
           (0.0 <= y   && y   <= 0.5) && 
           (0.5 <= z   && z   <= 1.0) ) 
        return 3;
      
      // subtet #4 E4 := 0.0 <= s + t <= 0.5 && 0.5 <= r + s <= 1.0 && 0.5 <= r + t <= 1.0 && 0.0 <= r <= 0.5
      if ( (0.0 <= yz && yz <= 0.5) && 
           (0.5 <= xy && xy <= 1.0) && 
           (0.5 <= xz && xz <= 1.0) && 
           (0.0 <= x  && x  <= 0.5) ) 
        return 4;
      
      // subtet #5 E5 := 0.5 <= r + s <= 1.0 && 0.5 <= s + t <= 1.0 && 0.5 <= r + t <= 1.0 && 0.75 <= r + s + t <= 1.0
      if ( (0.5  <= xy  && xy  <= 1.0) && 
           (0.5  <= yz  && yz  <= 1.0) && 
           (0.5  <= xz  && xz  <= 1.0) && 
           (0.75 <= xyz && xyz <= 1.0) ) 
        return 5;
      
      // subtet #6 E6 := 0.5 <= s + t <= 1.0 && 0.0 <= r + s <= 0.5 && 0.5 <= r + t <= 1.0 && 0.0 <= t <= 0.5
      if ( (0.5 <= yz && yz <= 1.0) && 
           (0.0 <= xy && xy <= 0.5) && 
           (0.5 <= xz && xz <= 1.0) && 
           (0.0 <= z  && z  <= 0.5) ) 
        return 6;

      // subtet #7 E7 := 0.0 <= s + t <= 0.5 && 0.0 <= r + s <= 0.5 && 0.5 <= r + t <= 1.0 && 0.0 <= s <= 0.25
      if ( (0.0 <= yz && yz <= 0.5) && 
           (0.0 <= xy && xy <= 0.5) &&
           (0.5 <= xz && xz <= 1.0) && 
           (0.0 <= y  && y  <= 0.25) ) 
        return 7;
      
      // subtet #8 E8 := 0.0 <= r + t <= 0.5 && 0.0 <= s + t <= 0.5 &&  0.5 <= r + s <= 1.0 && 0.0 <= t <= 0.25
      if ( (0.0 <= xz && xz <= 0.5) && 
           (0.0 <= yz && yz <= 0.5) &&
           (0.5 <= xy && xy <= 1.0) &&
           (0.0 <= z  && z  <= 0.25) )
        return 8;
      
      // subtet #9 E9 := 0.0 <= r + t <= 0.5 && 0.5 <= r + s <= 1.0 &&  0.5 <= s + t <= 1.0 && 0.0 <= s <= 0.5
      if ( (0.0 <= xz && xz <= 0.5) && 
           (0.5 <= xy && xy <= 1.0) && 
           (0.5 <= yz && yz <= 1.0) && 
           (0.0 <= y  && y  <= 0.5) ) 
        return 9;
      
      // subtet #10 E10 := 0.0 <= r + t <= 0.5 && 0.5 <= s + t <= 1.0 && 0.0 <= r + s <= 0.5 && 0.0 <= r <= 0.25
      if ( (0.0 <= xz && xz <= 0.5) && 
           (0.5 <= yz && yz <= 1.0) && 
           (0.0 <= xy && xy <= 0.5) && 
           (0.0 <= x  && x  <= 0.25) )
        return 10;
      
      // subtet #11 E11 := 0.5 <= r + s + t <= 0.75 && 0.0 <= r + t <= 0.5 && 0.0 <= s + t <= 0.5 && 0.0 <= r + s <= 0.5
      if ( (0.5 <= xyz && xyz <= 0.75) && 
           (0.0 <= xz  && xz  <= 0.5) &&
           (0.0 <= yz  && yz  <= 0.5) && 
           (0.0 <= xy  && xy  <= 0.5) ) 
        return 11;

      return -1;
    }
    
    template<EOperator opType>
    template<typename outputValueViewType,
             typename inputPointViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_TET_COMP12_FEM::Serial<opType>::
    getValues(       outputValueViewType output,
               const inputPointViewType  input ) {
      switch (opType) {
      case OPERATOR_VALUE: {
        const typename inputPointViewType::value_type r = input(0);
        const typename inputPointViewType::value_type s = input(1);
        const typename inputPointViewType::value_type t = input(2);

        // initialize output 
        for (auto i=0;i<10;++i) 
          output.access(i) = 0.0;
        
        const auto subtet = getLocalSubTet( r, s, t );
        
        // idependent verification shows that a given point will produce
        // the same shape functions for each tet that contains it, so
        // we only need to use the first one returned.
        if (subtet != -1) {
          typename inputPointViewType::value_type aux = 0.0;
          switch (subtet) {
          case 0:
            output.access(0) = 1. - 2. * (r + s + t);
            output.access(4) = 2. * r;
            output.access(6) = 2. * s;
            output.access(7) = 2. * t;
            break;
          case 1:
            output.access(1) = 2. * r - 1.;
            output.access(4) = 2. - 2. * (r + s + t);
            output.access(5) = 2. * s;
            output.access(8) = 2. * t;
            break;
          case 2:
            output.access(2) = 2. * s - 1.;
            output.access(5) = 2. * r;
            output.access(6) = 2. - 2. * (r + s + t);
            output.access(9) = 2. * t;
            break;
          case 3:
            output.access(3) = 2. * t - 1.;
            output.access(7) = 2. - 2. * (r + s + t);
            output.access(8) = 2. * r;
            output.access(9) = 2. * s;
            break;
          case 4:
            output.access(4) = 1. - 2. * (s + t);
            output.access(5) = 2. * (r + s) - 1.;
            output.access(8) = 2. * (r + t) - 1.;
            aux = 2. - 4. * r;
            break;
          case 5:
            output.access(5) = 2. * (r + s) - 1.;
            output.access(8) = 2. * (r + t) - 1.;
            output.access(9) = 2. * (s + t) - 1.;
            aux = 4. - 4. * (r + s + t);
            break;
          case 6:
            output.access(7) = 1. - 2. * (r + s);
            output.access(8) = 2. * (r + t) - 1.;
            output.access(9) = 2. * (s + t) - 1.;
            aux = 2. - 4. * t;
            break;
          case 7:
            output.access(4) = 1. - 2. * (s + t);
            output.access(7) = 1. - 2. * (r + s);
            output.access(8) = 2. * (r + t) - 1.;
            aux = 4. * s;
            break;
          case 8:
            output.access(4) = 1. - 2. * (s + t);
            output.access(5) = 2. * (r + s) - 1.;
            output.access(6) = 1. - 2. * (r + t);
            aux = 4. * t;
            break;
          case 9:
            output.access(5) = 2. * (r + s) - 1.;
            output.access(6) = 1. - 2. * (r + t);
            output.access(9) = 2. * (s + t) - 1.;
            aux = 2. - 4. * s;
            break;
          case 10:
            output.access(6) = 1. - 2. * (r + t);
            output.access(7) = 1. - 2. * (r + s);
            output.access(9) = 2. * (s + t) - 1.;
            aux = 4. * r;
            break;
          case 11:
            output.access(4) = 1. - 2. * (s + t);
            output.access(6) = 1. - 2. * (r + t);
            output.access(7) = 1. - 2. * (r + s);
            aux = 4. * (r + s + t) - 2.;
            break;
          }
          for (auto i=4;i<10;++i)
            output.access(i) += aux/6.0;
        }
        break;
      }
      case OPERATOR_GRAD: {
        const typename inputPointViewType::value_type r = input(0);
        const typename inputPointViewType::value_type s = input(1);
        const typename inputPointViewType::value_type t = input(2);
        
        output.access(0,0) = (-17 + 20*r + 20*s + 20*t)/8.;
        output.access(0,1) = (-17 + 20*r + 20*s + 20*t)/8.;
        output.access(0,2) = (-17 + 20*r + 20*s + 20*t)/8.;
        output.access(1,0) = -0.375 + (5*r)/2.;
        output.access(1,1) = 0.;
        output.access(1,2) = 0.;
        output.access(2,0) = 0.;
        output.access(2,1) = -0.375 + (5*s)/2.;
        output.access(2,2) = 0.;
        output.access(3,0) = 0.;
        output.access(3,1) = 0.;
        output.access(3,2) = -0.375 + (5*t)/2.;
        output.access(4,0) = (-35*(-1 + 2*r + s + t))/12.;
        output.access(4,1) = (-4 - 35*r + 5*s + 10*t)/12.;
        output.access(4,2) = (-4 - 35*r + 10*s + 5*t)/12.;
        output.access(5,0) = (-1 + 5*r + 40*s - 5*t)/12.;
        output.access(5,1) = (-1 + 40*r + 5*s - 5*t)/12.;
        output.access(5,2) = (-5*(-1 + r + s + 2*t))/12.;
        output.access(6,0) = (-4 + 5*r - 35*s + 10*t)/12.;
        output.access(6,1) = (-35*(-1 + r + 2*s + t))/12.;
        output.access(6,2) = (-4 + 10*r - 35*s + 5*t)/12.;
        output.access(7,0) = (-4 + 5*r + 10*s - 35*t)/12.;
        output.access(7,1) = (-4 + 10*r + 5*s - 35*t)/12.;
        output.access(7,2) = (-35*(-1 + r + s + 2*t))/12.;
        output.access(8,0) = (-1 + 5*r - 5*s + 40*t)/12.;
        output.access(8,1) = (-5*(-1 + r + 2*s + t))/12.;
        output.access(8,2) = (-1 + 40*r - 5*s + 5*t)/12.;
        output.access(9,0) = (-5*(-1 + 2*r + s + t))/12.;
        output.access(9,1) = (-1 - 5*r + 5*s + 40*t)/12.;
        output.access(9,2) = (-1 - 5*r + 40*s + 5*t)/12.;
        break;
      }
      case OPERATOR_MAX: {
        const ordinal_type jend = output.extent(1);
        const ordinal_type iend = output.extent(0);

        for (ordinal_type j=0;j<jend;++j)
          for (auto i=0;i<iend;++i)
            output.access(i, j) = 0.0;
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true , 
                                  ">>> ERROR (Basis_HGRAD_TET_COMP12_FEM): Operator type not implemented" );
      }
      }
    }

    template<typename DT, 
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void
    Basis_HGRAD_TET_COMP12_FEM::
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
      case OPERATOR_D1: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_GRAD> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      case OPERATOR_CURL: {
        INTREPID2_TEST_FOR_EXCEPTION( operatorType == OPERATOR_CURL, std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_TET_COMP12_FEM): CURL is invalid operator for rank-0 (scalar) functions in 3D");
        break;
      }
      case OPERATOR_DIV: {
        INTREPID2_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_TET_COMP12_FEM): DIV is invalid operator for rank-0 (scalar) functions in 3D");
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
        INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_TET_COMP12_FEM): Operator type not implemented" );
      }
      }
    }

  }

  // -------------------------------------------------------------------------------------
  template<typename DT, typename OT, typename PT>
  Basis_HGRAD_TET_COMP12_FEM<DT,OT,PT>::
  Basis_HGRAD_TET_COMP12_FEM() {
    const ordinal_type spaceDim = 3;
    this->basisCardinality_     = 10;
    this->basisDegree_          = 1;
    this->basisCellTopologyKey_ = shards::Tetrahedron<4>::key;
    this->basisType_            = BASIS_FEM_DEFAULT;
    this->basisCoordinates_     = COORDINATES_CARTESIAN;
    this->functionSpace_        = FUNCTION_SPACE_HGRAD;

    {
      // Basis-dependent intializations
      const ordinal_type tagSize  = 4;        // size of DoF tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
      
      // An array with local DoF tags assigned to basis functions, in the order of their local enumeration
      ordinal_type tags[]  = { 0, 0, 0, 1,
                               0, 1, 0, 1,
                               0, 2, 0, 1,
                               0, 3, 0, 1,
                               1, 0, 0, 1,
                               1, 1, 0, 1,
                               1, 2, 0, 1,
                               1, 3, 0, 1,
                               1, 4, 0, 1,
                               1, 5, 0, 1 };
      
      // host view
      OrdinalTypeArray1DHost tagView(&tags[0],40);

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
      dofCoords("dofCoordsHost", this->basisCardinality_, spaceDim);

    dofCoords(0,0) = 0.0;   dofCoords(0,1) = 0.0; dofCoords(0,2) = 0.0;
    dofCoords(1,0) = 1.0;   dofCoords(1,1) = 0.0; dofCoords(1,2) = 0.0;
    dofCoords(2,0) = 0.0;   dofCoords(2,1) = 1.0; dofCoords(2,2) = 0.0;
    dofCoords(3,0) = 0.0;   dofCoords(3,1) = 0.0; dofCoords(3,2) = 1.0;
    dofCoords(4,0) = 0.5;   dofCoords(4,1) = 0.0; dofCoords(4,2) = 0.0;
    dofCoords(5,0) = 0.5;   dofCoords(5,1) = 0.5; dofCoords(5,2) = 0.0;
    dofCoords(6,0) = 0.0;   dofCoords(6,1) = 0.5; dofCoords(6,2) = 0.0;
    dofCoords(7,0) = 0.0;   dofCoords(7,1) = 0.0; dofCoords(7,2) = 0.5;
    dofCoords(8,0) = 0.5;   dofCoords(8,1) = 0.0; dofCoords(8,2) = 0.5;
    dofCoords(9,0) = 0.0;   dofCoords(9,1) = 0.5; dofCoords(9,2) = 0.5;

    this->dofCoords_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);
  }
}

#endif
