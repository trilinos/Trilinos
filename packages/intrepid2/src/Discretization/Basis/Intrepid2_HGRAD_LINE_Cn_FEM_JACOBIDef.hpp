// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_LINE_Cn_FEM_JACOBIDef.hpp
    \brief  Definition file for FEM orthogonal basis functions of degree n for H(grad) functions on LINE.
    \author Created by P. Bochev and D. Ridzal and R. Kirby
            Kokkorized by Kyungjoo Kim
*/


#ifndef __INTREPID2_HGRAD_LINE_CN_FEM_JACOBI_DEF_HPP__
#define __INTREPID2_HGRAD_LINE_CN_FEM_JACOBI_DEF_HPP__

namespace Intrepid2 {
  // -------------------------------------------------------------------------------------

  namespace Impl {

    // output (N,P,D)
    // input  (P,D) - assumes that it has a set of points to amortize the function call cost for jacobi polynomial.
    template<EOperator opType>
    template<typename OutputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_LINE_Cn_FEM_JACOBI::Serial<opType>::
    getValues(       OutputViewType output,
               const inputViewType  input,
               const ordinal_type   order,
               const double         alpha,
               const double         beta,
               const ordinal_type   operatorDn ) {
      // cardinality of the evaluation order
      const ordinal_type card = order + 1;
      ordinal_type opDn = operatorDn;

      const auto pts = Kokkos::subview( input, Kokkos::ALL(), 0 );
      const ordinal_type np = input.extent(0);

      switch (opType) {
      case OPERATOR_VALUE: {
        const Kokkos::View<typename inputViewType::value_type*,
          typename inputViewType::memory_space,Kokkos::MemoryUnmanaged> null;
        for (ordinal_type p=0;p<card;++p) {
          auto poly = Kokkos::subview( output, p, Kokkos::ALL() );
          Polylib::Serial::JacobiPolynomial(np, pts, poly, null, p, alpha, beta);
        }
        break;
      }
      case OPERATOR_GRAD:
      case OPERATOR_D1: {
        for (ordinal_type p=0;p<card;++p) {
          auto polyd = Kokkos::subview( output, p, Kokkos::ALL(), 0 );
          Polylib::Serial::JacobiPolynomialDerivative(np, pts, polyd, p, alpha, beta);
        }
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
      case OPERATOR_D10:
        opDn = getOperatorOrder(opType);
      case OPERATOR_Dn: {
        {
          const ordinal_type pend = output.extent(0);
          const ordinal_type iend = output.extent(1);
          const ordinal_type jend = output.extent(2);

          for (ordinal_type p=0;p<pend;++p)
            for (ordinal_type i=0;i<iend;++i)
              for (ordinal_type j=0;j<jend;++j)
                output.access(p, i, j) = 0.0;
        }
        {
          const Kokkos::View<typename inputViewType::value_type*,
            typename inputViewType::memory_space,Kokkos::MemoryUnmanaged> null;

          for (ordinal_type p=opDn;p<card;++p) {
            double scaleFactor = 1.0;
            for (ordinal_type i=1;i<=opDn;++i)
              scaleFactor *= 0.5*(p + alpha + beta + i);

            const auto poly = Kokkos::subview( output, p, Kokkos::ALL(), 0 );
            Polylib::Serial::JacobiPolynomial(np, pts, poly, null, p-opDn, alpha+opDn, beta+opDn);
            for (ordinal_type i=0;i<np;++i)
              poly(i) = scaleFactor*poly(i);
          }
        }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM_JACOBI::Serial::getValues) operator is not supported");
      }
      }
    }

    // -------------------------------------------------------------------------------------

    template<typename DT, ordinal_type numPtsPerEval,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void
    Basis_HGRAD_LINE_Cn_FEM_JACOBI::
    getValues( const typename DT::execution_space& space,
                     Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const ordinal_type order,
               const double alpha,
               const double beta,
               const EOperator operatorType ) {
      typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
      typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
      typedef typename DT::execution_space ExecSpaceType;

      // loopSize corresponds to the # of points
      const auto loopSizeTmp1 = (inputPoints.extent(0)/numPtsPerEval);
      const auto loopSizeTmp2 = (inputPoints.extent(0)%numPtsPerEval != 0);
      const auto loopSize = loopSizeTmp1 + loopSizeTmp2;
      Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(space, 0, loopSize);

      switch (operatorType) {
      case OPERATOR_VALUE: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_VALUE,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints,
                                                  order, alpha, beta) );
        break;
      }
      case OPERATOR_GRAD:
      case OPERATOR_D1: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_GRAD,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints,
                                                  order, alpha, beta) );
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
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_Dn,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints,
                                                  order, alpha, beta,
                                                  getOperatorOrder(operatorType)) );
        break;
      }
      case OPERATOR_DIV:
      case OPERATOR_CURL: {
        INTREPID2_TEST_FOR_EXCEPTION( operatorType == OPERATOR_DIV ||
                                      operatorType == OPERATOR_CURL, std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_LINE_Cn_FEM_JACOBI): invalid operator type (div and curl).");
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( !Intrepid2::isValidOperator(operatorType), std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_LINE_Cn_FEM_JACOBI): invalid operator type");
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------

  template<typename DT, typename OT, typename PT>
  Basis_HGRAD_LINE_Cn_FEM_JACOBI<DT,OT,PT>::
  Basis_HGRAD_LINE_Cn_FEM_JACOBI( const ordinal_type order,
                                  const double alpha,
                                  const double beta ) {
    this->basisCardinality_     = order+1;
    this->basisDegree_          = order;
    this->basisCellTopologyKey_ = shards::Line<>::key;
    this->basisType_            = BASIS_FEM_HIERARCHICAL;
    this->basisCoordinates_     = COORDINATES_CARTESIAN;
    this->functionSpace_        = FUNCTION_SPACE_HGRAD;

    // jacobi
    this->alpha_ = alpha;
    this->beta_  = beta;

    // initialize tags
    {
      // Basis-dependent intializations
      const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

      ordinal_type tags[Parameters::MaxOrder+1][4];
      const ordinal_type card = this->basisCardinality_;
      for (ordinal_type i=0;i<card;++i) {
        tags[i][0] = 1;     // these are all "internal" i.e. "volume" DoFs
        tags[i][1] = 0;     // there is only one line
        tags[i][2] = i;     // local DoF id
        tags[i][3] = card;  // total number of DoFs
      }

      OrdinalTypeArray1DHost tagView(&tags[0][0], card*4);

      // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
      // tags are constructed on host
      this->setOrdinalTagData(this->tagToOrdinal_,
                              this->ordinalToTag_,
                              tagView,
                              this->basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
    }

    // dof coords is not applicable to hierarchical functions
  }


}

#endif
