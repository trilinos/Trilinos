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

/** \file   Intrepid_HGRAD_TRI_Cn_FEM_ORTHDef.hpp
    \brief  Definition file for FEM orthogonal basis functions of arbitrary degree
    for H(grad) functions on TRI.
    \author Created by R. Kirby
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_TRI_CN_FEM_ORTH_DEF_HPP__
#define __INTREPID2_HGRAD_TRI_CN_FEM_ORTH_DEF_HPP__


namespace Intrepid2 {
  // -------------------------------------------------------------------------------------

  namespace Impl {

    template<>
    template<ordinal_type maxOrder,
             ordinal_type maxNumPts,
             typename outputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void OrthPolynomial<0>::generate( /**/  outputViewType output,
                                      const inputViewType input,
                                      const ordinal_type order ) {
      typedef outputViewType::value_type value_type;

      const ordinal_type
        npts = input.dimension(0);

      const auto z = input;

      // each point needs to be transformed from Pavel's element
      // z(i,0) --> (2.0 * z(i,0) - 1.0)
      // z(i,1) --> (2.0 * z(i,1) - 1.0)

      auto idx = [](const ordinal_type p,
                    const ordinal_type q) {
        return (p+q)*(p+q+1)/2+q;
      };

      auto jrc = [](const value_type alpha,
                    const value_type beta ,
                    const ordinal_type n ,
                    /**/  value_type  &an,
                    /**/  value_type  &bn,
                    /**/  value_type  &cn) {
        an = ( (2.0 * n + 1.0 + alpha + beta) * ( 2.0 * n + 2.0 + alpha + beta ) /
               (2.0 * ( n + 1 ) * ( n + 1 + alpha + beta ) ) );
        bn = ( (alpha*alpha-beta*beta)*(2.0*n+1.0+alpha+beta) /
               (2.0*(n+1.0)*(2.0*n+alpha+beta)*(n+1.0+alpha+beta) ) );
        cn = ( (n+alpha)*(n+beta)*(2.0*n+2.0+alpha+beta) /
               ( (n+1.0)*(n+1.0+alpha+beta)*(2.0*n+alpha+beta) ) );
      }

      // set D^{0,0} = 1.0
      {
        const ordinal_type loc = idx(0,0);
        for (ordinal_type i=0;i<npts;++i)
          output(loc, i) = 1.0 + z(i,0) - z(i,0) + z(i,1) - z(i,1);
      }

      if (p > 0) {
        value_type f1[maxNumPts],f2[maxNumPts],f3[maxNumPts];

        for (ordinal_type i=0;i<npts;++i) {
          f1[i] = 0.5 * (1.0+2.0*(2.0*z(i,0)-1.0)+(2.0*z(i,1)-1.0));
          f2[i] = 0.5 * (1.0-(2.0*z(i,1)-1.0));
          f3[i] = f2[i] * f2[i];
        }

        // set D^{1,0} = f1
        {
          const ordinal_type loc = idx(1,0);
          for (ordinal_type i=0;i<npts;++i)
            output(loc, i) = f1[i];
        }

        // recurrence in p
        for (ordinal_type p=1;p<order;p++) {
          const ordinal_type
            loc = idx(p,0),
            loc_p1 = idx(p+1,0),
            loc_m1 = idx(p-1,0);

          const value_type
            a = (2.0*p+1.0)/(1.0+p),
            b = p / (p+1.0);

          for (ordinal_type i=0;i<npts;++i)
            output(loc_p1,i) = ( a * f1[i] * output(loc,i) -
                                 b * f3[i] * output(loc_m1,i) );
        }

        // D^{p,1}
        for (ordinal_type p=0;p<order;++p) {
          const ordinal_type
            loc_p_0 = idx(p,0),
            loc_p_1 = idx(p,1);

          for (ordinal_type i=0;i<npts;++i)
            output(loc_p_1,i) = output(loc_p_0,i)*0.5*(1.0+2.0*p+(3.0+2.0*p)*(2.0*z(i,1)-1.0));
        }


        // recurrence in q
        for (ordinal_type p=0;p<order-1;++p)
          for (ordinal_type q=1;q<order-p;++q) {
            const ordinal_type
              loc_p_qp1 = idx(p,q+1),
              loc_p_q = idx(p,q),
              loc_p_qm1 = idx(p,q-1);

            value_type a,b,c;
            jrc((value_type)(2*p+1),(value_type)0,q,a,b,c);
            for (ordinal_type i=0;i<npts;++i)
              output(loc_p_qp1,i) = ( (a*(2.0*z(i,1)-1.0)+b)*outputValues(loc_p_q,i)
                                      - c*outputValues(loc_p_qm1,i) );
          }
      }
      
      // orthogonalize
      for (ordinal_type p=0;p<=order;++p)
        for (ordinal_type q=0;q<=order-p;++q)
          for (ordinal_type i=0;i<npts;++i)
            output(idx(p,q),i) *= sqrt( (p+0.5)*(p+q+1.0));
    }
    
    template<>
    template<ordinal_type maxOrder,
             ordinal_type maxNumPts,
             typename outputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void OrthPolynomial<1>::generate( /**/  outputViewType output,
                                      const inputViewType input,
                                      const ordinal_type p ) {
      typedef outputViewType::value_type value_type;
      typedef Sacado::Fad::SFad<value_type,2> fad_type;

      constexpr ordinal_type maxCard = (maxOrder+1)*(maxOrder+2)/2;

      const ordinal_type
        npts = input.dimension(0),
        card = output.dimension(0);

      // use stack buffer
      fad_type inBuf[maxNumPts][2], outBuf[maxCard][maxNumPts];

      Kokkos::View<fad_type**, Kokkos::Impl::ActiveExecutionMemorySpace> in(&inBuf[0][0],            npts, 2);
      Kokkos::View<fad_type***,Kokkos::Impl::ActiveExecutionMemorySpace> out(&outBuf[0][0][0], card, npts);

      for (ordinal_type i=0;i<npts;++i)
        for (ordinal_type j=0;j<2;++j) {
          in(i,j) = Sacado::Fad::SFad<value_type,2>( input(i,j) );
          in(i,j).diff(j,2);
        }

      OrthPolynomial<0>::generate<maxOrder,maxNumPts>(out, in, p);

      for (ordinal_type i=0;i<card;++i)
        for (ordinal_type j=0;j<npts;++j)
          for (ordinal_type k=0;k<2;++k)
            output(i,j,k) = out(i,j).dx(k);
    }

    // when n >= 2, use recursion
    template<ordinal_type n>
    template<ordinal_type maxOrder,
             ordinal_type maxNumPts,
             typename outputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void OrthPolynomial<n>::generate( /**/  outputViewType output,
                                      const inputViewType input,
                                      const ordinal_type p ) {
      typedef outputViewType::value_type value_type;
      typedef Sacado::Fad::SFad<value_type,2> fad_type;

      constexpr ordinal_type maxCard = (maxOrder+1)*(maxOrder+2)/2;

      const ordinal_type
        npts = input.dimension(0),
        card = output.dimension(0);

      // use stack buffer
      fad_type inBuf[maxNumPts][2], outBuf[maxCard][maxNumPts][n+1];

      Kokkos::View<fad_type**, Kokkos::Impl::ActiveExecutionMemorySpace> in(&inBuf[0][0],            npts, 2);
      Kokkos::View<fad_type***,Kokkos::Impl::ActiveExecutionMemorySpace> out(&outBuf[0][0][0], card, npts, n+1);

      for (ordinal_type i=0;i<npts;++i)
        for (ordinal_type j=0;j<2;++j) {
          in(i,j) = Sacado::Fad::SFad<value_type,2>( input(i,j) );
          in(i,j).diff(j,2);
        }

      if (n == 2)
        OrthPolynomial<1>::generate<maxOrder,maxNumPts>(out, in, p);
      else
        OrthPolynomial<n-1>::generate<maxOrder,maxNumPts>(out, in, p);

      for (ordinal_type i=0;i<card;++i)
        for (ordinal_type j=0;j<npts;++j) {
          output(i,j,0) = out(i,j,0).dx(0);
          for (ordinal_type k=0;k<n;++k)
            output(i,j,k+1) = out(i,j,k).dx(1);
        }
    }



    template<EOperator opType>
    template<typename outputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_TRI_Cn_FEM_ORTH::Serial<opType>::
    getValues( /**/  outputViewType output,
               const inputViewType  input,
               const ordinal_type   order,
               const ordinal_type   operatorDn ) {
      ordinal_type opDn = operatorDn;

      switch (opType) {
      case OPERATOR_VALUE: {
        OrthPolynomial<0>( output, order, input );
        break;
      }
      case OPERATOR_GRAD:
      case OPERATOR_D1: {
        OrthPolynomial<1>( output, order, input );
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
        OrthPolynomial( output, order, input, opDn );
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH::Serial::getValues) operator is not supported");
      }
      }
    }

    // -------------------------------------------------------------------------------------

    template<typename SpT, ordinal_type numPtsPerEval,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void
    Basis_HGRAD_TRI_Cn_FEM_ORTH::
    getValues( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const ordinal_type order,
               const EOperator operatorType ) {
      typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
      typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
      typedef typename ExecSpace<typename inputPointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

      // loopSize corresponds to the # of points
      const auto loopSizeTmp1 = (inputPoints.dimension(0)/numPtsPerEval);
      const auto loopSizeTmp2 = (inputPoints.dimension(0)%numPtsPerEval != 0);
      const auto loopSize = loopSizeTmp1 + loopSizeTmp2;
      Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);

      switch (operatorType) {
      case OPERATOR_VALUE: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_VALUE,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints,
                                                  order) );
        break;
      }
      case OPERATOR_GRAD:
      case OPERATOR_D1: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_GRAD,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints,
                                                  order) );
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
                                                  order, getOperatorOrder(operatorType)) );
        break;
      }
      case OPERATOR_DIV:
      case OPERATOR_CURL: {
        INTREPID2_TEST_FOR_EXCEPTION( operatorType == OPERATOR_DIV ||
                                      operatorType == OPERATOR_CURL, std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_TRI_Cn_FEM_ORTH): invalid operator type (div and curl).");
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( !Intrepid2::isValidOperator(operatorType), std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_TRI_Cn_FEM_ORTH): invalid operator type");
      }
      }
    }
  }


  // -------------------------------------------------------------------------------------

  template<typename SpT, typename OT, typename PT>
  Basis_HGRAD_TRI_Cn_FEM_ORTH<SpT,OT,PT>::
  Basis_HGRAD_TRI_Cn_FEM_ORTH( const ordinal_type order ) {
    this->basisCardinality_  = (order+1)*(order+2)/2;
    this->basisDegree_       = order;
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >() );
    this->basisType_         = BASIS_FEM_HIERARCHICAL;
    this->basisCoordinates_  = COORDINATES_CARTESIAN;

    // initialize tags
    {
      // Basis-dependent initializations
      const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

      constexpr ordinal_type maxCard = (Parameters::MaxOrder+1)*(Parameters::MaxOrder+2)/2;
      ordinal_type tags[maxCard][4];
      const ordinal_type card = this->basisCardinality_;
      for (ordinal_type i=0;i<card;++i) {
        tags[i][0] = 2;     // these are all "internal" i.e. "volume" DoFs
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
