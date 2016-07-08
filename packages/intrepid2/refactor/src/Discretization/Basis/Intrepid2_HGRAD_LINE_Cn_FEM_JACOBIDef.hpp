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

/** \file   Intrepid_HGRAD_LINE_Cn_FEM_JACOBIDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(grad) orthogonal on LINE.
    \author Created by P. Bochev and D. Ridzal and R. Kirby
            Kokkorized by Kyungjoo Kim
*/


#ifndef __INTREPID2_HGRAD_LINE_CN_FEM_JACOBI_DEF_HPP__
#define __INTREPID2_HGRAD_LINE_CN_FEM_JACOBI_DEF_HPP__

namespace Intrepid2 {
  // -------------------------------------------------------------------------------------

  template<typename SpT, typename OT, typename PT>
  template<EOperator opType>
  template<typename outputValueValueType, class ...outputValueProperties,
           typename inputPointValueType,  class ...inputPointProperties>
  KOKKOS_INLINE_FUNCTION
  void
  Basis_HGRAD_LINE_Cn_FEM_JACOBI<SpT,OT,PT>::Serial<opType>::
  getValues( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> output,
             const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  input,
             const ordinal_type order,
             const double alpha,
             const double beta,
             const ordinal_type opDn ) {
    switch (opType) {
    case OPERATOR_VALUE: {
      const auto np = input.dimension(0);
      const Kokkos::View<ValueType*,Kokkos::Impl::ActiveExecutionMemorySpace,Kokkos::MemoryUnmanaged> null;

      Polylib::JacobiPolynomial(np, input, output, null, order, alpha, beta);
      break;
    }
    case OPERATOR_GRAD: {
      const auto np = input.dimension(0);
      const auto polyd = Kokkos::subdynrankview( output, Kokkos::ALL(), 0 );

      Polylib::JacobiPolynomialDerivative(np, input, polyd, order, alpha, beta);      
      break;
    }
    case OPERATOR_MAX: {
      if (order < opDn) {
        const auto jend = output.dimension(1);
        const auto iend = output.dimension(0);
        
        for (auto j=0;j<jend;++j)
          for (auto i=0;i<iend;++i)
            output(i, j) = 0.0;
      } else {
        double scaleFactor = 1.0;
        for (auto i=1;i<=opDn;++i) 
          scaleFactor *= 0.5*(order + alpha + beta + i);

        const auto np = input.dimension(0);
        const auto poly = Kokkos::subdynrankview( output, Kokkos::ALL(), 0 );        
        const Kokkos::View<ValueType*,Kokkos::Impl::ActiveExecutionMemorySpace,Kokkos::MemoryUnmanaged> null;

        Polylib::JacobiPolynomial(np, input, poly, null, order-opDn, alpha+opDn, beta+opDn);
        for (auto i=0;i<np;++i) 
          poly(i) = scaleFactor*poly(i);
      }
      break;
    }
    default: {
      INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                opType != OPERATOR_GRAD &&
                                opType != OPERATOR_MAX,
                                ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM_JACOBI::Serial::getValues) operator is not supported");
    }
    }
  }


  // -------------------------------------------------------------------------------------

  template<typename SpT, typename OT, typename PT>
  Basis_HGRAD_LINE_Cn_FEM_JACOBI<SpT,OT,PT>::
  Basis_HGRAD_LINE_Cn_FEM_JACOBI( const ordinal_type order, 
                                  const double alpha, 
                                  const double beta ) {
    this->basisCardinality_  = order+1;
    this->basisDegree_       = order;    
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
    this->basisType_         = BASIS_FEM_HIERARCHICAL;
    this->basisCoordinates_  = COORDINATES_CARTESIAN;

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
      
      ordinal_type tags[MaxOrder+1][4];
      const auto card = this->basisCardinality_;
      for (auto i=0;i<card;++i) {
        tags[i][0] = 1;     // these are all "internal" i.e. "volume" DoFs
        tags[i][1] = 0;     // there is only one line
        tags[i][2] = i;     // local DoF id 
        tags[i][3] = card;  // total number of DoFs 
      }
     
      ordinal_type_array_1d_host tagView(&tag[0][0], card*4);
 
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

  template<typename SpT, typename OT, typename PT>
  template<typename outputValueValueType, class ...outputValueProperties,
           typename inputPointValueType,  class ...inputPointProperties>
  void 
  Basis_HGRAD_LINE_Cn_FEM_JACOBI<SpT,OT,PT>::Internal::
  getValues(  /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
              const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
              const EOperator operatorType ) const {
#ifdef HAVE_INTREPID2_DEBUG
    Intrepid2::getValues_HGRAD_Args(outputValues,
                                    inputPoints,
                                    operatorType,
                                    obj_->getBaseCellTopology(),
                                    obj_->getCardinality() );
#endif
    
    typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
    typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
    typedef typename ExecSpace<typename inputPointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

    // loopSize corresponds to cardinality
    const auto loopSize = outputValues.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);

    switch (operatorType) {
    case OPERATOR_VALUE: {
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_VALUE> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, 
                                                obj_->alpha_, obj_->beta_) );
      break;
    }
    case OPERATOR_GRAD:
    case OPERATOR_D1: {
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_GRAD> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, 
                                                obj_->alpha_, obj_->beta_) );
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
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, 
                                                obj_->alpha_, obj_->beta_, 
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

#endif
