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
    Kokkorized by Kyungjoo Kim and Mauro Perego
*/

#ifndef __INTREPID2_HGRAD_TRI_CN_FEM_ORTH_DEF_HPP__
#define __INTREPID2_HGRAD_TRI_CN_FEM_ORTH_DEF_HPP__

namespace Intrepid2 {
  // -------------------------------------------------------------------------------------

  namespace Impl {

    template<typename outputViewType,
             typename inputViewType,
             bool hasDeriv>
    KOKKOS_INLINE_FUNCTION
    void OrthPolynomialTri<outputViewType,inputViewType,hasDeriv,0>::generate(
                                                                              outputViewType output,
                                                                              const inputViewType input,
                                                                              const ordinal_type order ) {

      constexpr ordinal_type spaceDim = 2;
      constexpr ordinal_type maxNumPts = Parameters::MaxNumPtsPerBasisEval;

      typedef typename outputViewType::value_type value_type;

      auto output0 = (hasDeriv) ? Kokkos::subview(output,  Kokkos::ALL(), Kokkos::ALL(),0) : Kokkos::subview(output,  Kokkos::ALL(), Kokkos::ALL());

      const ordinal_type
        npts = input.dimension(0);

      const auto z = input;

      // each point needs to be transformed from Pavel's element
      // z(i,0) --> (2.0 * z(i,0) - 1.0)
      // z(i,1) --> (2.0 * z(i,1) - 1.0)
  

      // set D^{0,0} = 1.0
      {
        const ordinal_type loc = Intrepid2::getPnEnumeration<spaceDim>(0,0);
        for (ordinal_type i=0;i<npts;++i) {
          output0(loc, i) = 1.0;
          if(hasDeriv) {
            output(loc,i,1) = 0;
            output(loc,i,2) = 0;
          }
        }
      }

      if (order > 0) {
        value_type f1[maxNumPts],f2[maxNumPts], df2_1[maxNumPts];;
        value_type df1_0, df1_1;

        for (ordinal_type i=0;i<npts;++i) {
          f1[i] = 0.5 * (1.0+2.0*(2.0*z(i,0)-1.0)+(2.0*z(i,1)-1.0));   // \eta_1 * (1 - \eta_2)/2
          f2[i] = std::pow(z(i,1)-1,2);  //( (1 - \eta_2)/2 )^2
          if(hasDeriv) {
            df1_0 = 2.0;
            df1_1 = 1.0;
            df2_1[i] = 2.0*(z(i,1)-1);
          }
        }

        // set D^{1,0} = f1
        {
          const ordinal_type loc = Intrepid2::getPnEnumeration<spaceDim>(1,0);
          for (ordinal_type i=0;i<npts;++i) {
            output0(loc, i) = f1[i];
            if(hasDeriv) {
              output(loc,i,1) = df1_0;
              output(loc,i,2) = df1_1;
            }
          }
        }

        // recurrence in p
        for (ordinal_type p=1;p<order;p++) {
          const ordinal_type
            loc = Intrepid2::getPnEnumeration<spaceDim>(p,0),
            loc_p1 = Intrepid2::getPnEnumeration<spaceDim>(p+1,0),
            loc_m1 = Intrepid2::getPnEnumeration<spaceDim>(p-1,0);

          const value_type
            a = (2.0*p+1.0)/(1.0+p),
            b = p / (p+1.0);

          for (ordinal_type i=0;i<npts;++i) {
            output0(loc_p1,i) = ( a * f1[i] * output0(loc,i) -
                                  b * f2[i] * output0(loc_m1,i) );
            if(hasDeriv) {
              output(loc_p1,i,1) =  a * (f1[i] * output(loc,i,1) + df1_0 * output0(loc,i))  -
                b * f2[i] * output(loc_m1,i,1) ;
              output(loc_p1,i,2) =  a * (f1[i] * output(loc,i,2) + df1_1 * output0(loc,i))  -
                b * (df2_1[i] * output0(loc_m1,i) + f2[i] * output(loc_m1,i,2)) ;
            }
          }
        }

        // D^{p,1}
        for (ordinal_type p=0;p<order;++p) {
          const ordinal_type
            loc_p_0 = Intrepid2::getPnEnumeration<spaceDim>(p,0),
            loc_p_1 = Intrepid2::getPnEnumeration<spaceDim>(p,1);

          for (ordinal_type i=0;i<npts;++i) {
            output0(loc_p_1,i) = output0(loc_p_0,i)*0.5*(1.0+2.0*p+(3.0+2.0*p)*(2.0*z(i,1)-1.0));
            if(hasDeriv) {
              output(loc_p_1,i,1) = output(loc_p_0,i,1)*0.5*(1.0+2.0*p+(3.0+2.0*p)*(2.0*z(i,1)-1.0));
              output(loc_p_1,i,2) = output(loc_p_0,i,2)*0.5*(1.0+2.0*p+(3.0+2.0*p)*(2.0*z(i,1)-1.0)) + output0(loc_p_0,i)*(3.0+2.0*p);
            }
          }
        }


        // recurrence in q
        for (ordinal_type p=0;p<order-1;++p)
          for (ordinal_type q=1;q<order-p;++q) {
            const ordinal_type
              loc_p_qp1 = Intrepid2::getPnEnumeration<spaceDim>(p,q+1),
              loc_p_q = Intrepid2::getPnEnumeration<spaceDim>(p,q),
              loc_p_qm1 = Intrepid2::getPnEnumeration<spaceDim>(p,q-1);

            value_type a,b,c;
            Intrepid2::getJacobyRecurrenceCoeffs(a,b,c, 2*p+1,0,q);
            for (ordinal_type i=0;i<npts;++i) {
              output0(loc_p_qp1,i) =  (a*(2.0*z(i,1)-1.0)+b)*output0(loc_p_q,i)
                - c*output0(loc_p_qm1,i) ;
              if(hasDeriv) {
                output(loc_p_qp1,i,1) =  (a*(2.0*z(i,1)-1.0)+b)*output(loc_p_q,i,1)
                  - c*output(loc_p_qm1,i,1) ;
                output(loc_p_qp1,i,2) =  (a*(2.0*z(i,1)-1.0)+b)*output(loc_p_q,i,2) +2*a*output0(loc_p_q,i)
                  - c*output(loc_p_qm1,i,2) ;
              }
            }
          }
      }

      // orthogonalize
      for (ordinal_type p=0;p<=order;++p)
        for (ordinal_type q=0;q<=order-p;++q)
          for (ordinal_type i=0;i<npts;++i) {
            output0(Intrepid2::getPnEnumeration<spaceDim>(p,q),i) *= std::sqrt( (p+0.5)*(p+q+1.0));
            if(hasDeriv) {
              output(Intrepid2::getPnEnumeration<spaceDim>(p,q),i,1) *= std::sqrt( (p+0.5)*(p+q+1.0));
              output(Intrepid2::getPnEnumeration<spaceDim>(p,q),i,2) *= std::sqrt( (p+0.5)*(p+q+1.0));
            }
          }
    }

    template<typename outputViewType,
             typename inputViewType,
             bool hasDeriv>
    KOKKOS_INLINE_FUNCTION
    void OrthPolynomialTri<outputViewType,inputViewType,hasDeriv,1>::generate(
                                                                              outputViewType output,
                                                                              const inputViewType input,
                                                                              const ordinal_type order ) {
      typedef typename outputViewType::value_type value_type;
      typedef typename outputViewType::pointer_type pointer_type;
      typedef typename Kokkos::DynRankView<value_type, Kokkos::Impl::ActiveExecutionMemorySpace> outViewType;

      constexpr ordinal_type spaceDim = 2;
      constexpr ordinal_type maxCard = Intrepid2::getPnCardinality<spaceDim, Parameters::MaxOrder>();
      const ordinal_type
        npts = input.dimension(0),
        card = output.dimension(0);

      // use stack buffer
      value_type outBuf[maxCard][Parameters::MaxNumPtsPerBasisEval][spaceDim+1]; //spaceDim for derivatives, 1 for value
      // outViewType out = Kokkos::createDynRankViewWithType<outViewType>(output, (pointer_type)(&outBuf[0][0][0]), card, npts, spaceDim+1); // createDynRankViewWithType not KOKKOS_INLINE_FUNCTION

      typedef typename Kokkos::View<value_type***, Kokkos::Impl::ActiveExecutionMemorySpace> outHackViewType; // Issue trying to do this directly with DynRankView - no matching ctor
      outHackViewType hack_view( (pointer_type)&outBuf[0][0][0], card, npts, spaceDim+1); // As a hack, wrapped with a View then wrapped the View with a DynRankView
      outViewType out(hack_view);

      OrthPolynomialTri<outViewType,inputViewType,hasDeriv,0>::generate(out, input, order);
      for (ordinal_type i=0;i<card;++i)
        for (ordinal_type j=0;j<npts;++j)
          for (ordinal_type k=0;k<spaceDim;++k)
            output(i,j,k) = out(i,j,k+1);
    }

    // when n >= 2, use recursion
    template<typename outputViewType,
             typename inputViewType,
             bool hasDeriv,
             ordinal_type n>
    KOKKOS_INLINE_FUNCTION
    void OrthPolynomialTri<outputViewType,inputViewType,hasDeriv,n>::generate(
                                                                              outputViewType output,
                                                                              const inputViewType input,
                                                                              const ordinal_type order ) {

      constexpr ordinal_type spaceDim = 2;
      constexpr ordinal_type maxCard = Intrepid2::getPnCardinality<spaceDim, Parameters::MaxOrder>();
  
      typedef typename outputViewType::value_type value_type;
      typedef Sacado::Fad::SFad<value_type,spaceDim> fad_type;

      const ordinal_type
        npts = input.dimension(0),
        card = output.dimension(0);

      // use stack buffer
      fad_type inBuf[Parameters::MaxNumPtsPerBasisEval][spaceDim],
        outBuf[maxCard][Parameters::MaxNumPtsPerBasisEval][n];

      typedef typename Kokkos::View<fad_type***, Kokkos::Impl::ActiveExecutionMemorySpace> outViewType;
      typedef typename Kokkos::View<fad_type**, Kokkos::Impl::ActiveExecutionMemorySpace> inViewType;

      inViewType in((value_type*)(&inBuf[0][0]), npts, spaceDim);
      outViewType out((value_type*)(&outBuf[0][0][0]), card, npts, n);

      for (ordinal_type i=0;i<npts;++i)
        for (ordinal_type j=0;j<spaceDim;++j) {
          in(i,j) = input(i,j);
          in(i,j).diff(j,spaceDim);
        }

      OrthPolynomialTri<outViewType,inViewType,hasDeriv,n-1>::generate(out, in, order);
      for (ordinal_type i=0;i<card;++i)
        for (ordinal_type j=0;j<npts;++j) {
          for (ordinal_type i_dx = 0; i_dx <= n; ++i_dx) {
            ordinal_type i_dy =  n-i_dx;
            ordinal_type i_Dn = i_dy;
            if(i_dx > 0) {
              //n=2:  (f_x)_x, (f_y)_x
              //n=3:  (f_xx)_x, (f_xy)_x, (f_yy)_x
              ordinal_type i_Dnm1 = i_dy;
              output(i,j,i_Dn) = out(i,j,i_Dnm1).dx(0);
            }
            else {
              //n=2:  (f_y)_y, (f_z)_y
              //n=3:  (f_yy)_y
              ordinal_type i_Dnm1 = i_dy-1;
              output(i,j,i_Dn) = out(i,j,i_Dnm1).dx(1);
            }
          }
        }
    }


    template<EOperator opType>
    template<typename outputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_TRI_Cn_FEM_ORTH::Serial<opType>::
    getValues( outputViewType output,
               const inputViewType  input,
               const ordinal_type   order) {
      switch (opType) {
      case OPERATOR_VALUE: {
        OrthPolynomialTri<outputViewType,inputViewType,false,0>::generate( output, input, order );
        break;
      }
      case OPERATOR_GRAD:
      case OPERATOR_D1:
        {
          OrthPolynomialTri<outputViewType,inputViewType,true,1>::generate( output, input, order );
          break;
        }
        /*
          case OPERATOR_D2: {
          OrthPolynomialTri<outputViewType,inputViewType,true,2>::generate( output, input, order );
          break;
          }
          case OPERATOR_D3: {
          OrthPolynomialTri<outputViewType,inputViewType,true,3>::generate( output, input, order );
          break;
          }
          case OPERATOR_D4: {
          OrthPolynomialTri<outputViewType,inputViewType,true,4>::generate( output, input, order );
          break;
          }
          case OPERATOR_D5: {
          OrthPolynomialTri<outputViewType,inputViewType,true,5>::generate( output, input, order );
          break;
          }
          case OPERATOR_D6: {
          OrthPolynomialTri<outputViewType,inputViewType,true,6>::generate( output, input, order );
          break;
          }
          case OPERATOR_D7: {
          OrthPolynomialTri<outputViewType,inputViewType,true,7>::generate( output, input, order );
          break;
          }
          case OPERATOR_D8: {
          OrthPolynomialTri<outputViewType,inputViewType,true,8>::generate( output, input, order );
          break;
          }
          case OPERATOR_D9: {
          OrthPolynomialTri<outputViewType,inputViewType,true,9>::generate( output, input, order );
          break;
          }
          case OPERATOR_D10: {
          OrthPolynomialTri<outputViewType,inputViewType,true,10>::generate( output, input, order );
          break;
          }
        */
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
    getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
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
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, order) );
        break;
      }
      case OPERATOR_GRAD:
      case OPERATOR_D1: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D1,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, order) );
        break;
      }
      case OPERATOR_CURL: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_CURL,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, order) );
        break;
      }
      case OPERATOR_D2:{
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D2,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, order) );
        break;
      }
      case OPERATOR_D3: {
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D3,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, order) );
        break;
      }
        /*
          case OPERATOR_D4: {
          typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D4,numPtsPerEval> FunctorType;
          Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, order) );
          break;
          }
          case OPERATOR_D5: {
          typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D5,numPtsPerEval> FunctorType;
          Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, order) );
          break;
          }
          case OPERATOR_D6: {
          typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D6,numPtsPerEval> FunctorType;
          Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, order) );
          break;
          }
          case OPERATOR_D7: {
          typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D7,numPtsPerEval> FunctorType;
          Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, order) );
          break;
          }
          case OPERATOR_D8: {
          typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D8,numPtsPerEval> FunctorType;
          Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, order) );
          break;
          }
          case OPERATOR_D9: {
          typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D9,numPtsPerEval> FunctorType;
          Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, order) );
          break;
          }
          case OPERATOR_D10: {
          typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D10,numPtsPerEval> FunctorType;
          Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, order) );
          break;
          }
        */
      case OPERATOR_DIV: {
        INTREPID2_TEST_FOR_EXCEPTION( operatorType == OPERATOR_DIV ||
                                      operatorType == OPERATOR_CURL, std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_TRI_Cn_FEM_ORTH): invalid operator type div.");
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

    constexpr ordinal_type spaceDim = 2;

    this->basisCardinality_  = Intrepid2::getPnCardinality<spaceDim>(order);
    this->basisDegree_       = order;
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >() );
    this->basisType_         = BASIS_FEM_HIERARCHICAL;
    this->basisCoordinates_  = COORDINATES_CARTESIAN;

    // initialize tags
    {
      // Basis-dependent initializations
      constexpr ordinal_type tagSize  = 4;    // size of DoF tag, i.e., number of fields in the tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

      constexpr ordinal_type maxCard = Intrepid2::getPnCardinality<spaceDim, Parameters::MaxOrder>();
      ordinal_type tags[maxCard][tagSize];
      const ordinal_type card = this->basisCardinality_;
      for (ordinal_type i=0;i<card;++i) {
        tags[i][0] = 2;     // these are all "internal" i.e. "volume" DoFs
        tags[i][1] = 0;     // there is only one line
        tags[i][2] = i;     // local DoF id
        tags[i][3] = card;  // total number of DoFs
      }

      ordinal_type_array_1d_host tagView(&tags[0][0], card*tagSize);

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
