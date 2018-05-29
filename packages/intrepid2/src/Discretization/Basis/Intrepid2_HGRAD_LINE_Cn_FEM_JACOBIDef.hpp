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
    template<typename outputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_LINE_Cn_FEM_JACOBI::Serial<opType>::
    getValues(       outputViewType output,
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
    
    template<typename SpT, ordinal_type numPtsPerEval,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void 
    Basis_HGRAD_LINE_Cn_FEM_JACOBI::
    getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const ordinal_type order,
               const double alpha,
               const double beta,
               const EOperator operatorType ) {
      typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
      typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
      typedef typename ExecSpace<typename inputPointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

      // loopSize corresponds to the # of points
      const auto loopSizeTmp1 = (inputPoints.extent(0)/numPtsPerEval);
      const auto loopSizeTmp2 = (inputPoints.extent(0)%numPtsPerEval != 0);
      const auto loopSize = loopSizeTmp1 + loopSizeTmp2;
      Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
      
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

  template<typename SpT, typename OT, typename PT>
  Basis_HGRAD_LINE_Cn_FEM_JACOBI<SpT,OT,PT>::
  Basis_HGRAD_LINE_Cn_FEM_JACOBI( const ordinal_type order, 
                                  const double alpha, 
                                  const double beta ) {
    this->basisCardinality_  = order+1;
    this->basisDegree_       = order;    
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
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
      
      ordinal_type tags[Parameters::MaxOrder+1][4];
      const ordinal_type card = this->basisCardinality_;
      for (ordinal_type i=0;i<card;++i) {
        tags[i][0] = 1;     // these are all "internal" i.e. "volume" DoFs
        tags[i][1] = 0;     // there is only one line
        tags[i][2] = i;     // local DoF id 
        tags[i][3] = card;  // total number of DoFs 
      }
     
      ordinal_type_array_1d_host tagView(&tags[0][0], card*4);
 
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
