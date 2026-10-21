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

    // Returns the index (0 - 11) of the subtet containing the point (x,y,z), with the subtet
    // numbering used in Serial<OPERATOR_VALUE>::getValues:
    //
    //   subtets #0..#3  : corner tets of the regular refinement, #k = { lambda_k >= 1/2 }
    //   subtets #4..#11 : the inner octahedron { 0 <= lambda_k <= 1/2 } split into 8 tets
    //                     by the centroid c = (1/4,1/4,1/4); each one is the cone from c
    //                     over one octahedron face, lambda_k = 1/2 or lambda_k = 0
    //
    // where lambda_0 = 1-x-y-z, lambda_1 = x, lambda_2 = y, lambda_3 = z are the barycentric
    // coordinates of the point. A point on an interface may be assigned to either adjacent
    // subtet: the basis functions are continuous, so the result does not depend on the choice.
    //
    // The classification is total: it never fails, including for points on the boundary that
    // round-off puts slightly outside the reference tet. A point outside the tet is assigned
    // to an adjacent subtet, so basis values there are the linear extension from that subtet.
    template<typename pointValueType>
    KOKKOS_INLINE_FUNCTION
    ordinal_type
    Basis_HGRAD_TET_COMP12_FEM::getLocalSubTet( const pointValueType x, 
                                                const pointValueType y, 
                                                const pointValueType z ) {
      // only the point values are needed (not the derivatives, for Fad types)
      typedef typename ScalarTraits<pointValueType>::scalar_type scalar_type;
      const scalar_type
        r = get_scalar_value(x),
        s = get_scalar_value(y),
        t = get_scalar_value(z);

      const scalar_type lambda[4] = { scalar_type(1.0) - r - s - t, r, s, t };

      // corner subtets #0..#3: the point is in subtet #k if lambda_k >= 1/2
      ordinal_type k = 0;
      for (ordinal_type i=1;i<4;++i)
        if (lambda[i] > lambda[k]) k = i;
      if (lambda[k] >= 0.5)
        return k;

      // inner octahedron: the point is in the cone over the octahedron face that is hit first
      // by the ray from the centroid through the point, i.e. the face lambda_k = 1/2 (if
      // lambda_k > 1/4) or lambda_k = 0 (if lambda_k < 1/4) for the k maximizing |lambda_k - 1/4|
      k = 0;
      for (ordinal_type i=1;i<4;++i) {
        const scalar_type di = lambda[i] - 0.25, dk = lambda[k] - 0.25;
        if (di*di > dk*dk) k = i;
      }
      const ordinal_type subTetInnerFace[4] = { 11, 4, 9, 6 }; // face lambda_k = 1/2
      const ordinal_type subTetOuterFace[4] = {  5, 10, 7, 8 }; // face lambda_k = 0
      return (lambda[k] > 0.25) ? subTetInnerFace[k] : subTetOuterFace[k];
    }
    
    template<EOperator OpType>
    template<typename outputValueViewType,
             typename inputPointViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_TET_COMP12_FEM::Serial<OpType>::
    getValues(       outputValueViewType output,
               const inputPointViewType  input ) {

      if constexpr (OpType == OPERATOR_VALUE) {
        const typename inputPointViewType::value_type r = input(0);
        const typename inputPointViewType::value_type s = input(1);
        const typename inputPointViewType::value_type t = input(2);

        // initialize output 
        for (auto i=0;i<10;++i) 
          output(i) = 0.0;
        
        // the basis functions are continuous, so a point on an interface between subtets
        // produces the same values whichever of the adjacent subtets is returned
        const auto subtet = getLocalSubTet( r, s, t );

        typename inputPointViewType::value_type aux = 0.0;
        switch (subtet) {
        case 0:
          output(0) = 1. - 2. * (r + s + t);
          output(4) = 2. * r;
          output(6) = 2. * s;
          output(7) = 2. * t;
          break;
        case 1:
          output(1) = 2. * r - 1.;
          output(4) = 2. - 2. * (r + s + t);
          output(5) = 2. * s;
          output(8) = 2. * t;
          break;
        case 2:
          output(2) = 2. * s - 1.;
          output(5) = 2. * r;
          output(6) = 2. - 2. * (r + s + t);
          output(9) = 2. * t;
          break;
        case 3:
          output(3) = 2. * t - 1.;
          output(7) = 2. - 2. * (r + s + t);
          output(8) = 2. * r;
          output(9) = 2. * s;
          break;
        case 4:
          output(4) = 1. - 2. * (s + t);
          output(5) = 2. * (r + s) - 1.;
          output(8) = 2. * (r + t) - 1.;
          aux = 2. - 4. * r;
          break;
        case 5:
          output(5) = 2. * (r + s) - 1.;
          output(8) = 2. * (r + t) - 1.;
          output(9) = 2. * (s + t) - 1.;
          aux = 4. - 4. * (r + s + t);
          break;
        case 6:
          output(7) = 1. - 2. * (r + s);
          output(8) = 2. * (r + t) - 1.;
          output(9) = 2. * (s + t) - 1.;
          aux = 2. - 4. * t;
          break;
        case 7:
          output(4) = 1. - 2. * (s + t);
          output(7) = 1. - 2. * (r + s);
          output(8) = 2. * (r + t) - 1.;
          aux = 4. * s;
          break;
        case 8:
          output(4) = 1. - 2. * (s + t);
          output(5) = 2. * (r + s) - 1.;
          output(6) = 1. - 2. * (r + t);
          aux = 4. * t;
          break;
        case 9:
          output(5) = 2. * (r + s) - 1.;
          output(6) = 1. - 2. * (r + t);
          output(9) = 2. * (s + t) - 1.;
          aux = 2. - 4. * s;
          break;
        case 10:
          output(6) = 1. - 2. * (r + t);
          output(7) = 1. - 2. * (r + s);
          output(9) = 2. * (s + t) - 1.;
          aux = 4. * r;
          break;
        case 11:
          output(4) = 1. - 2. * (s + t);
          output(6) = 1. - 2. * (r + t);
          output(7) = 1. - 2. * (r + s);
          aux = 4. * (r + s + t) - 2.;
          break;
        }
        for (auto i=4;i<10;++i)
          output(i) += aux/6.0;
      }
      else if constexpr (OpType == OPERATOR_GRAD) {
        const typename inputPointViewType::value_type r = input(0);
        const typename inputPointViewType::value_type s = input(1);
        const typename inputPointViewType::value_type t = input(2);

        // These are NOT the gradients of the piecewise linear functions returned by
        // OPERATOR_VALUE (which are piecewise constant on the 12 subtets), but their L2
        // projection, over the parent tet, onto linear polynomials. See the class documentation.
        output(0,0) = (-17 + 20*r + 20*s + 20*t)/8.;
        output(0,1) = (-17 + 20*r + 20*s + 20*t)/8.;
        output(0,2) = (-17 + 20*r + 20*s + 20*t)/8.;
        output(1,0) = -0.375 + (5*r)/2.;
        output(1,1) = 0.;
        output(1,2) = 0.;
        output(2,0) = 0.;
        output(2,1) = -0.375 + (5*s)/2.;
        output(2,2) = 0.;
        output(3,0) = 0.;
        output(3,1) = 0.;
        output(3,2) = -0.375 + (5*t)/2.;
        output(4,0) = (-35*(-1 + 2*r + s + t))/12.;
        output(4,1) = (-4 - 35*r + 5*s + 10*t)/12.;
        output(4,2) = (-4 - 35*r + 10*s + 5*t)/12.;
        output(5,0) = (-1 + 5*r + 40*s - 5*t)/12.;
        output(5,1) = (-1 + 40*r + 5*s - 5*t)/12.;
        output(5,2) = (-5*(-1 + r + s + 2*t))/12.;
        output(6,0) = (-4 + 5*r - 35*s + 10*t)/12.;
        output(6,1) = (-35*(-1 + r + 2*s + t))/12.;
        output(6,2) = (-4 + 10*r - 35*s + 5*t)/12.;
        output(7,0) = (-4 + 5*r + 10*s - 35*t)/12.;
        output(7,1) = (-4 + 10*r + 5*s - 35*t)/12.;
        output(7,2) = (-35*(-1 + r + s + 2*t))/12.;
        output(8,0) = (-1 + 5*r - 5*s + 40*t)/12.;
        output(8,1) = (-5*(-1 + r + 2*s + t))/12.;
        output(8,2) = (-1 + 40*r - 5*s + 5*t)/12.;
        output(9,0) = (-5*(-1 + 2*r + s + t))/12.;
        output(9,1) = (-1 - 5*r + 5*s + 40*t)/12.;
        output(9,2) = (-1 - 5*r + 40*s + 5*t)/12.;
      }
      else if constexpr (OpType == OPERATOR_MAX) {
        const ordinal_type jend = output.extent(1);
        const ordinal_type iend = output.extent(0);

        for (ordinal_type j=0;j<jend;++j)
          for (auto i=0;i<iend;++i)
            output(i, j) = 0.0;
      }
      else {
        INTREPID2_TEST_FOR_ABORT( true , 
                                  ">>> ERROR (Basis_HGRAD_TET_COMP12_FEM): Operator type not implemented" );
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

  template<typename DT, typename OT, typename PT>
  void 
  Basis_HGRAD_TET_COMP12_FEM<DT,OT,PT>::getScratchSpaceSize(       
                                    ordinal_type& perThreadSpaceSize,
                              const PointViewType inputPoints,
                              const EOperator operatorType) const {
    perThreadSpaceSize = 0;
  }

  template<typename DT, typename OT, typename PT>
  KOKKOS_INLINE_FUNCTION
  void 
  Basis_HGRAD_TET_COMP12_FEM<DT,OT,PT>::getValues(       
          OutputViewType outputValues,
      const PointViewType  inputPoints,
      const EOperator operatorType,
      const typename Kokkos::TeamPolicy<typename DT::execution_space>::member_type& team_member,
      const int threadScratchLevel, 
      const ordinal_type subcellDim,
      const ordinal_type subcellOrdinal) const {

      INTREPID2_TEST_FOR_ABORT( !((subcellDim <= 0) && (subcellOrdinal == -1)),
        ">>> ERROR: (Intrepid2::Basis_HGRAD_TET_COMP12_FEM::getValues), The capability of selecting subsets of basis functions has not been implemented yet.");

      (void) threadScratchLevel; //avoid unused variable warning

      const int numPoints = inputPoints.extent(0);

      switch(operatorType) {
        case OPERATOR_VALUE:
          Kokkos::parallel_for (Kokkos::TeamThreadRange (team_member, numPoints), [=] (const ordinal_type pt) {
            auto       output = Kokkos::subview( outputValues, Kokkos::ALL(), pt, Kokkos::ALL() );
            const auto input  = Kokkos::subview( inputPoints,                 pt, Kokkos::ALL() );
            Impl::Basis_HGRAD_TET_COMP12_FEM::Serial<OPERATOR_VALUE>::getValues( output, input);
          });
          break;
        case OPERATOR_GRAD:
          Kokkos::parallel_for (Kokkos::TeamThreadRange (team_member, numPoints), [=] (const ordinal_type pt) {
            auto       output = Kokkos::subview( outputValues, Kokkos::ALL(), pt, Kokkos::ALL() );
            const auto input  = Kokkos::subview( inputPoints,                 pt, Kokkos::ALL() );
            Impl::Basis_HGRAD_TET_COMP12_FEM::Serial<OPERATOR_GRAD>::getValues( output, input);
          });
          break;
        default: {}
    }
  }
  
}// namespace Intrepid2
#endif
