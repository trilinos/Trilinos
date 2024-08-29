// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_TET_Cn_FEM_ORTHDef.hpp
    \brief  Definition file for FEM orthogonal basis functions of degree n
    for H(grad) functions on TET cells.
    \author Created by R. Kirby
    Kokkorized by Kyungjoo Kim and Mauro Perego
 */
#ifndef __INTREPID2_HGRAD_TET_CN_FEM_ORTH_DEF_HPP__
#define __INTREPID2_HGRAD_TET_CN_FEM_ORTH_DEF_HPP__

namespace Intrepid2 {
// -------------------------------------------------------------------------------------

namespace Impl {

template<typename OutputViewType,
typename inputViewType,
typename workViewType,
bool hasDeriv>
KOKKOS_INLINE_FUNCTION
void OrthPolynomialTet<OutputViewType,inputViewType,workViewType,hasDeriv,0>::generate(
    OutputViewType output,
    const inputViewType input,
    workViewType  /*work*/,
    const ordinal_type order ) {

  constexpr ordinal_type spaceDim = 3;
  constexpr ordinal_type maxNumPts = Parameters::MaxNumPtsPerBasisEval;

  typedef typename OutputViewType::value_type value_type;

  auto output0 = (hasDeriv) ? Kokkos::subview(output,  Kokkos::ALL(), Kokkos::ALL(),0) : Kokkos::subview(output,  Kokkos::ALL(), Kokkos::ALL());

  const ordinal_type
  npts = input.extent(0);

  const auto z = input;

  // each point needs to be transformed from Pavel's element
  // z(i,0) --> (2.0 * z(i,0) - 1.0)
  // z(i,1) --> (2.0 * z(i,1) - 1.0)

  // set D^{0,0,0} = 1.0
  {
    const ordinal_type loc = Intrepid2::getPnEnumeration<spaceDim>(0,0,0);
    for (ordinal_type i=0;i<npts;++i) {
      output0(loc, i) = 1.0;
      if(hasDeriv) {
        output.access(loc,i,1) = 0;
        output.access(loc,i,2) = 0;
        output.access(loc,i,3) = 0;
      }
    }
  }

  if (order > 0) {
    value_type f1[maxNumPts]={},f2[maxNumPts]={}, f3[maxNumPts]={}, f4[maxNumPts]={},f5[maxNumPts]={},
        df2_1[maxNumPts]={}, df2_2[maxNumPts]={}, df5_2[maxNumPts]={};
    value_type df1_0, df1_1, df1_2, df3_1, df3_2, df4_2;

    for (int i=0;i<npts;i++) {
      f1[i] = 2.0*z(i,0) + z(i,1) + z(i,2)- 1.0;
      f2[i] = pow(z(i,1) + z(i,2) -1.0, 2);
      f3[i] = 2.0*z(i,1) + z(i,2) -1.0;
      f4[i] = 1.0 - z(i,2);
      f5[i] = f4[i] * f4[i];
      if(hasDeriv) {
        df2_1[i] = 2*(z(i,1) + z(i,2) -1.0);
        df2_2[i] = df2_1[i];
        df5_2[i] = -2*f4[i];
      }
    }

    df1_0 = 2.0;
    df1_1 = 1.0;
    df1_2 = 1.0;
    df3_1 = 2;
    df3_2 = 1;
    df4_2 = -1;

    // set D^{1,0,0} = f1
        {
      const ordinal_type loc = Intrepid2::getPnEnumeration<spaceDim>(1,0,0);
      for (ordinal_type i=0;i<npts;++i) {
        output0(loc, i) = f1[i];
        if(hasDeriv) {
          output.access(loc,i,1) = df1_0;
          output.access(loc,i,2) = df1_1;
          output.access(loc,i,3) = df1_2;
        }
      }
        }

        // recurrence in p
        // D^{p+1,0,0}
        for (ordinal_type p=1;p<order;p++) {
          const ordinal_type
          loc = Intrepid2::getPnEnumeration<spaceDim>(p,0,0),
          loc_p1 = Intrepid2::getPnEnumeration<spaceDim>(p+1,0,0),
          loc_m1 = Intrepid2::getPnEnumeration<spaceDim>(p-1,0,0);

          const value_type
          a = (2.0*p+1.0)/(1.0+p),
          b = p / (p+1.0);

          for (ordinal_type i=0;i<npts;++i) {
            output0(loc_p1,i) = ( a * f1[i] * output0(loc,i) -
                b * f2[i] * output0(loc_m1,i) );
            if(hasDeriv) {
              output.access(loc_p1,i,1) =  a * (f1[i] * output.access(loc,i,1) + df1_0 * output0(loc,i))  -
                  b * f2[i] * output.access(loc_m1,i,1) ;
              output.access(loc_p1,i,2) =  a * (f1[i] * output.access(loc,i,2) + df1_1 * output0(loc,i))  -
                  b * (df2_1[i] * output0(loc_m1,i) + f2[i] * output.access(loc_m1,i,2)) ;
              output.access(loc_p1,i,3) =  a * (f1[i] * output.access(loc,i,3) + df1_2 * output0(loc,i))  -
                  b * (df2_2[i] * output0(loc_m1,i) + f2[i] * output.access(loc_m1,i,3)) ;
            }
          }
        }

        // D^{p,1,0}
        for (ordinal_type p=0;p<order;++p) {
          const ordinal_type
          loc_p_0 = Intrepid2::getPnEnumeration<spaceDim>(p,0,0),
          loc_p_1 = Intrepid2::getPnEnumeration<spaceDim>(p,1,0);

          for (ordinal_type i=0;i<npts;++i) {
            const value_type coeff = (2.0 * p + 3.0) * z(i,1) + z(i,2)-1.0;
            output0(loc_p_1,i) = output0(loc_p_0,i)*coeff;
            if(hasDeriv) {
              output.access(loc_p_1,i,1) = output.access(loc_p_0,i,1)*coeff;
              output.access(loc_p_1,i,2) = output.access(loc_p_0,i,2)*coeff + output0(loc_p_0,i)*(2.0 * p + 3.0);
              output.access(loc_p_1,i,3) = output.access(loc_p_0,i,3)*coeff + output0(loc_p_0,i);
            }
          }
        }


        // recurrence in q
        // D^{p,q+1,0}
        for (ordinal_type p=0;p<order-1;++p)
          for (ordinal_type q=1;q<order-p;++q) {
            const ordinal_type
            loc_p_qp1 = Intrepid2::getPnEnumeration<spaceDim>(p,q+1,0),
            loc_p_q = Intrepid2::getPnEnumeration<spaceDim>(p,q,0),
            loc_p_qm1 = Intrepid2::getPnEnumeration<spaceDim>(p,q-1,0);

            value_type a,b,c, coeff, dcoeff_1, dcoeff_2;
            Intrepid2::getJacobyRecurrenceCoeffs(a,b,c, 2*p+1,0,q);
            for (ordinal_type i=0;i<npts;++i) {
              coeff = a * f3[i] + b * f4[i];
              dcoeff_1 = a * df3_1;
              dcoeff_2 = a * df3_2 + b * df4_2;
              output0(loc_p_qp1,i) =  coeff * output0(loc_p_q,i)
              - c* f5[i] * output0(loc_p_qm1,i) ;
              if(hasDeriv) {
                output.access(loc_p_qp1,i,1) =  coeff * output.access(loc_p_q,i,1) +
                    - c * f5[i] * output.access(loc_p_qm1,i,1) ;
                output.access(loc_p_qp1,i,2) =  coeff * output.access(loc_p_q,i,2)  + dcoeff_1 * output0(loc_p_q,i)
                - c * f5[i] * output.access(loc_p_qm1,i,2) ;
                output.access(loc_p_qp1,i,3) =  coeff * output.access(loc_p_q,i,3)  + dcoeff_2 * output0(loc_p_q,i)
                - c * f5[i] * output.access(loc_p_qm1,i,3) - c * df5_2[i] * output0(loc_p_qm1,i);
              }
            }
          }


        // D^{p,q,1}
        for (ordinal_type p=0;p<order;++p)
          for (ordinal_type q=0;q<order-p;++q) {

            const ordinal_type
            loc_p_q_0 = Intrepid2::getPnEnumeration<spaceDim>(p,q,0),
            loc_p_q_1 = Intrepid2::getPnEnumeration<spaceDim>(p,q,1);

            for (ordinal_type i=0;i<npts;++i) {
              const value_type coeff =  2.0 * ( 2.0 + q + p ) * z(i,2) - 1.0;
              output0(loc_p_q_1,i) = output0(loc_p_q_0,i) * coeff;
              if(hasDeriv) {
                output.access(loc_p_q_1,i,1) = output.access(loc_p_q_0,i,1) * coeff;
                output.access(loc_p_q_1,i,2) = output.access(loc_p_q_0,i,2) * coeff;
                output.access(loc_p_q_1,i,3) = output.access(loc_p_q_0,i,3) * coeff + 2.0 * ( 2.0 + q + p ) * output0(loc_p_q_0,i);
              }
            }
          }

        // general r recurrence
        // D^{p,q,r+1}
        for (ordinal_type p=0;p<order-1;++p)
          for (ordinal_type q=0;q<order-p-1;++q)
            for (ordinal_type r=1;r<order-p-q;++r) {
              const ordinal_type
              loc_p_q_rp1 = Intrepid2::getPnEnumeration<spaceDim>(p,q,r+1),
              loc_p_q_r = Intrepid2::getPnEnumeration<spaceDim>(p,q,r),
              loc_p_q_rm1 = Intrepid2::getPnEnumeration<spaceDim>(p,q,r-1);

              value_type a,b,c, coeff;
              Intrepid2::getJacobyRecurrenceCoeffs(a,b,c, 2*p+2*q+2,0,r);
              for (ordinal_type i=0;i<npts;++i) {
                coeff = 2.0 * a * z(i,2) - a + b;
                output0(loc_p_q_rp1,i) =  coeff * output0(loc_p_q_r,i) - c * output0(loc_p_q_rm1,i) ;
                if(hasDeriv) {
                  output.access(loc_p_q_rp1,i,1) =  coeff * output.access(loc_p_q_r,i,1) - c * output.access(loc_p_q_rm1,i,1);
                  output.access(loc_p_q_rp1,i,2) =  coeff * output.access(loc_p_q_r,i,2) - c * output.access(loc_p_q_rm1,i,2);
                  output.access(loc_p_q_rp1,i,3) =  coeff * output.access(loc_p_q_r,i,3) + 2 * a * output0(loc_p_q_r,i) - c * output.access(loc_p_q_rm1,i,3);
                }
              }
            }

  }

  // orthogonalize
  for (ordinal_type p=0;p<=order;++p)
    for (ordinal_type q=0;q<=order-p;++q)
      for (ordinal_type r=0;r<=order-p-q;++r) {
        ordinal_type loc_p_q_r = Intrepid2::getPnEnumeration<spaceDim>(p,q,r);
        value_type scal = std::sqrt( (p+0.5)*(p+q+1.0)*(p+q+r+1.5) );
        for (ordinal_type i=0;i<npts;++i) {
          output0(loc_p_q_r,i) *= scal;
          if(hasDeriv) {
            output.access(loc_p_q_r,i,1) *= scal;
            output.access(loc_p_q_r,i,2) *= scal;
            output.access(loc_p_q_r,i,3) *= scal;
          }
        }
      }
}

template<typename OutputViewType,
typename inputViewType,
typename workViewType,
bool hasDeriv>
KOKKOS_INLINE_FUNCTION
void OrthPolynomialTet<OutputViewType,inputViewType,workViewType,hasDeriv,1>::generate(
    OutputViewType output,
    const inputViewType  input,
    workViewType   work,
    const ordinal_type   order ) {
  constexpr ordinal_type spaceDim = 3;
  const ordinal_type
  npts = input.extent(0),
  card = output.extent(0);

  workViewType dummyView;
  OrthPolynomialTet<workViewType,inputViewType,workViewType,hasDeriv,0>::generate(work, input, dummyView, order);
  for (ordinal_type i=0;i<card;++i)
    for (ordinal_type j=0;j<npts;++j)
      for (ordinal_type k=0;k<spaceDim;++k)
        output.access(i,j,k) = work(i,j,k+1);
}


// when n >= 2, use recursion
template<typename OutputViewType,
typename inputViewType,
typename workViewType,
bool hasDeriv,
ordinal_type n>
KOKKOS_INLINE_FUNCTION
void OrthPolynomialTet<OutputViewType,inputViewType,workViewType,hasDeriv,n>::generate(
    OutputViewType /* output */,
    const inputViewType  /* input */,
    workViewType   /* work */,
    const ordinal_type   /* order */ ) {
#if 0 //#ifdef HAVE_INTREPID2_SACADO

constexpr ordinal_type spaceDim = 3;
constexpr ordinal_type maxCard = Intrepid2::getPnCardinality<spaceDim, Parameters::MaxOrder>();

typedef typename OutputViewType::value_type value_type;
typedef Sacado::Fad::SFad<value_type,spaceDim> fad_type;

const ordinal_type
npts = input.extent(0),
card = output.extent(0);

// use stack buffer
fad_type inBuf[Parameters::MaxNumPtsPerBasisEval][spaceDim];
fad_type outBuf[maxCard][Parameters::MaxNumPtsPerBasisEval][n*(n+1)/2];

typedef typename inputViewType::memory_space memory_space;
typedef typename inputViewType::memory_space memory_space;
typedef typename Kokkos::View<fad_type***, memory_space> outViewType;
typedef typename Kokkos::View<fad_type**, memory_space> inViewType;
auto vcprop = Kokkos::common_view_alloc_prop(input);

inViewType in(Kokkos::view_wrap((value_type*)&inBuf[0][0], vcprop), npts, spaceDim);
outViewType out(Kokkos::view_wrap((value_type*)&outBuf[0][0][0], vcprop), card, npts, n*(n+1)/2);

for (ordinal_type i=0;i<npts;++i)
  for (ordinal_type j=0;j<spaceDim;++j) {
    in(i,j) = input(i,j);
    in(i,j).diff(j,spaceDim);
  }

typedef typename Kokkos::DynRankView<fad_type, memory_space> outViewType_;
outViewType_ workView;
if (n==2) {
  //char outBuf[bufSize*sizeof(typename inViewType::value_type)];
  fad_type outBuf[maxCard][Parameters::MaxNumPtsPerBasisEval][spaceDim+1];
  auto vcprop = Kokkos::common_view_alloc_prop(in);
  workView = outViewType_( Kokkos::view_wrap((value_type*)&outBuf[0][0][0], vcprop), card, npts, spaceDim+1);
}
OrthPolynomialTet<outViewType,inViewType,outViewType_,hasDeriv,n-1>::generate(out, in, workView, order);

for (ordinal_type i=0;i<card;++i)
  for (ordinal_type j=0;j<npts;++j) {
    for (ordinal_type i_dx = 0; i_dx <= n; ++i_dx)
      for (ordinal_type i_dy = 0; i_dy <= n-i_dx; ++i_dy) {
        ordinal_type i_dz = n - i_dx - i_dy;
        ordinal_type i_Dn = Intrepid2::getDkEnumeration<spaceDim>(i_dx,i_dy,i_dz);
        if(i_dx > 0) {
          //n=2:  (f_x)_x, (f_y)_x, (f_z)_x
          //n=3:  (f_xx)_x, (f_xy)_x, (f_xz)_x, (f_yy)_x, (f_yz)_x, (f_zz)_x
          ordinal_type i_Dnm1 = Intrepid2::getDkEnumeration<spaceDim>(i_dx-1, i_dy, i_dz);
          output.access(i,j,i_Dn) = out(i,j,i_Dnm1).dx(0);
        }
        else if (i_dy > 0) {
          //n=2:  (f_y)_y, (f_z)_y
          //n=3:  (f_yy)_y, (f_yz)_y, (f_zz)_y
          ordinal_type i_Dnm1 = Intrepid2::getDkEnumeration<spaceDim>(i_dx, i_dy-1, i_dz);
          output.access(i,j,i_Dn) = out(i,j,i_Dnm1).dx(1);
        }
        else  {
          //n=2: (f_z)_z;
          //n=3: (f_zz)_z
          ordinal_type i_Dnm1 = Intrepid2::getDkEnumeration<spaceDim>(i_dx, i_dy, i_dz-1);
          output.access(i,j,i_Dn) = out(i,j,i_Dnm1).dx(2);
        }
      }
  }
#else
INTREPID2_TEST_FOR_ABORT( true,
    ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH::OrthPolynomialTri) Computing of second and higher-order derivatives is not currently supported");
#endif
}


template<EOperator opType>
template<typename OutputViewType,
typename inputViewType,
typename workViewType>
KOKKOS_INLINE_FUNCTION
void
Basis_HGRAD_TET_Cn_FEM_ORTH::Serial<opType>::
getValues( OutputViewType output,
    const inputViewType  input,
    workViewType   work,
    const ordinal_type   order) {
  switch (opType) {
  case OPERATOR_VALUE: {
    OrthPolynomialTet<OutputViewType,inputViewType,workViewType,false,0>::generate( output, input, work, order );
    break;
  }
  case OPERATOR_GRAD:
  case OPERATOR_D1: {
    OrthPolynomialTet<OutputViewType,inputViewType,workViewType,true,1>::generate( output, input, work, order );
    break;
  }
  case OPERATOR_D2: {
    OrthPolynomialTet<OutputViewType,inputViewType,workViewType,true,2>::generate( output, input, work, order );
    break;
  }
  default: {
    INTREPID2_TEST_FOR_ABORT( true,
        ">>> ERROR: (Intrepid2::Basis_HGRAD_TET_Cn_FEM_ORTH::Serial::getValues) operator is not supported");
  }
  }
}

// -------------------------------------------------------------------------------------

template<typename DT, ordinal_type numPtsPerEval,
typename outputValueValueType, class ...outputValueProperties,
typename inputPointValueType,  class ...inputPointProperties>
void
Basis_HGRAD_TET_Cn_FEM_ORTH::
getValues( 
    const typename DT::execution_space& space,
          Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
    const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
    const ordinal_type order,
    const EOperator operatorType ) {
  typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
  typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
  typedef typename DT::execution_space ExecSpaceType;

  // loopSize corresponds to the # of points
  const auto loopSizeTmp1 = (inputPoints.extent(0)/numPtsPerEval);
  const auto loopSizeTmp2 = (inputPoints.extent(0)%numPtsPerEval != 0);
  const auto loopSize = loopSizeTmp1 + loopSizeTmp2;
  Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(space, 0, loopSize);

  typedef typename inputPointViewType::value_type inputPointType;
  const ordinal_type cardinality = outputValues.extent(0);
  const ordinal_type spaceDim = 3;

  auto vcprop = Kokkos::common_view_alloc_prop(inputPoints);
  typedef typename Kokkos::DynRankView< inputPointType, typename inputPointViewType::memory_space> workViewType;

  switch (operatorType) {
  case OPERATOR_VALUE: {
    workViewType  dummyWorkView;
    typedef Functor<outputValueViewType,inputPointViewType,workViewType,OPERATOR_VALUE,numPtsPerEval> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, dummyWorkView, order) );
    break;
  }
  case OPERATOR_GRAD:
  case OPERATOR_D1: {
    workViewType  work(Kokkos::view_alloc(space, "Basis_HGRAD_TET_In_FEM_ORTH::getValues::work", vcprop), cardinality, inputPoints.extent(0), spaceDim+1);
    typedef Functor<outputValueViewType,inputPointViewType,workViewType,OPERATOR_D1,numPtsPerEval> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, work, order) );
    break;
  }
  case OPERATOR_D2:{
    workViewType  dummyWorkView;
    typedef Functor<outputValueViewType,inputPointViewType,workViewType,OPERATOR_D2,numPtsPerEval> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints ,dummyWorkView, order) );
    break;
  }
  default: {
    INTREPID2_TEST_FOR_EXCEPTION( !Intrepid2::isValidOperator(operatorType), std::invalid_argument,
        ">>> ERROR (Basis_HGRAD_TET_Cn_FEM_ORTH): invalid operator type");
  }
  }
}
}


// -------------------------------------------------------------------------------------
template<typename DT, typename OT, typename PT>
Basis_HGRAD_TET_Cn_FEM_ORTH<DT,OT,PT>::
Basis_HGRAD_TET_Cn_FEM_ORTH( const ordinal_type order ) {

  constexpr ordinal_type spaceDim = 3;
  this->basisCardinality_     = Intrepid2::getPnCardinality<spaceDim>(order);
  this->basisDegree_          = order;
  this->basisCellTopologyKey_ = shards::Tetrahedron<4>::key;
  this->basisType_            = BASIS_FEM_HIERARCHICAL;
  this->basisCoordinates_     = COORDINATES_CARTESIAN;
  this->functionSpace_        = FUNCTION_SPACE_HGRAD;

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

    OrdinalTypeArray1DHost tagView(&tags[0][0], card*tagSize);

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
