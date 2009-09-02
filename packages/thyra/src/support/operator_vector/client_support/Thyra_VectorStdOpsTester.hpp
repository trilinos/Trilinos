// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_VECTOR_STD_OPS_TESTER_HPP
#define THYRA_VECTOR_STD_OPS_TESTER_HPP

#include "Thyra_VectorStdOpsTesterDecl.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_arrayArg.hpp"

//#define THYRA_VECTOR_STD_OPS_TESTER_DUMP

#ifdef THYRA_VECTOR_STD_OPS_TESTER_DUMP
#  include "RTOpPack_SPMD_apply_op_decl.hpp"
#  include "Thyra_SpmdVectorBase.hpp"
#endif // THYRA_VECTOR_STD_OPS_TESTER_DUMP

namespace Thyra {

// VectorStdOpsTesterComparable (using partial specialization to only do tests in some cases)

template <bool isComparable, class Scalar>
class VectorStdOpsTesterComparable {
public:
  static bool checkComparableStdOps(
    const VectorSpaceBase<Scalar>                                 &vecSpc
    ,VectorBase<Scalar>                                           *z
    ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &error_tol
    ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &warning_tol
    ,std::ostream                                                 *out
    ,const bool                                                   &dumpAll
    )
    {
      return Teuchos::ScalarTraits<Scalar>::ThisShouldNotCompile();
    }
};

template <class Scalar>
class VectorStdOpsTesterComparable<false,Scalar> {
public:
  static bool checkComparableStdOps(
    const VectorSpaceBase<Scalar>                                 &vecSpc
    ,VectorBase<Scalar>                                           *z
    ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &error_tol
    ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &warning_tol
    ,std::ostream                                                 *out
    ,const bool                                                   &dumpAll
    )
    {
      if(out) *out << "\nThis scalar type does not support comparable operations so we can not test min(), max() and other such functions.\n";
      return true;
    }
};

template <class Scalar>
class VectorStdOpsTesterComparable<true,Scalar> {
public:
  static bool checkComparableStdOps(
    const VectorSpaceBase<Scalar>                                 &vecSpc
    ,VectorBase<Scalar>                                           *z
    ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &error_tol
    ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &warning_tol
    ,std::ostream                                                 *out
    ,const bool                                                   &dumpAll
    )
    {
      typedef Teuchos::ScalarTraits<Scalar> ST;

      bool success = true, result;
      
      if(out) *out << "\nTesting comparable operations ...\n";
      
      const Scalar scalarSmall(1e-5), scalarMedium(2.0), scalarLarge(100.0);
      if(out) *out << "\nassign(&*z,"<<scalarMedium<<");\n";
      assign(&*z,Scalar(scalarMedium));
      if(out && dumpAll) *out << "\nz =\n" << *z;
      if(out) *out << "\nset_ele(0,"<<scalarSmall<<",&*z);\n";
      set_ele(0,scalarSmall,&*z);
      if(out && dumpAll) *out << "\nz =\n" << *z;
      if(out) *out << "\nset_ele(1,"<<scalarLarge<<",&*z);\n";
      set_ele(1,scalarLarge,&*z);
      if(out && dumpAll) *out << "\nz =\n" << *z;
      if(out) *out << "\nset_ele(vecSpc.dim()-2,"<<scalarSmall<<",&*z);\n";
      set_ele(vecSpc.dim()-2,scalarSmall,&*z);
      if(out && dumpAll) *out << "\nz =\n" << *z;
      if(out) *out << "\nset_ele(vecSpc.dim()-1,"<<scalarLarge<<",&*z);\n";
      set_ele(vecSpc.dim()-1,scalarLarge,&*z);
      if(out && dumpAll) *out << "\nz =\n" << *z;

      Scalar minEle; Index minIndex;
      Scalar maxEle; Index maxIndex;

      if(!testRelErr<Scalar>(
           "min(*z)",min(*z),"scalarSmall",scalarSmall
           ,"error_tol",error_tol,"warning_tol",warning_tol,out
           )
        ) success=false;

      if(out) *out << "\nmin(*z,&minEle,&minIndex);\n";
      minEle = ST::zero(); minIndex = 0;
      min(*z,&minEle,&minIndex);
      if(!testRelErr<Scalar>(
           "minEle",minEle,"scalarSmall",scalarSmall
           ,"error_tol",error_tol,"warning_tol",warning_tol,out
           )
        ) success=false;
      result = minIndex == 0;
      if(out) *out << "\nminIndex = " << minIndex << " == 0 ? " << passfail(result) << std::endl;
      if(!result) success = false;

      if(out) *out << "\nminGreaterThanBound(*z,"<<scalarMedium<<",&minEle,&minIndex);\n";
      minEle = ST::zero(); minIndex = 0;
      minGreaterThanBound(*z,scalarMedium,&minEle,&minIndex);
      if(!testRelErr<Scalar>(
           "minEle",minEle,"scalarLarge",scalarLarge
           ,"error_tol",error_tol,"warning_tol",warning_tol,out
           )
        ) success=false;
      result = minIndex == 1;
      if(out) *out << "\nminIndex = " << minIndex << " == 1 ? " << passfail(result) << std::endl;
      if(!result) success = false;

      if(out) *out << "\nminGreaterThanBound(*z,"<<scalarLarge<<",&minEle,&minIndex);\n";
      minEle = ST::zero(); minIndex = 0;
      minGreaterThanBound(*z,scalarLarge,&minEle,&minIndex);
      result = minIndex < 0;
      if(out) *out << "\nminIndex = " << minIndex << " < 0 ? " << passfail(result) << std::endl;
      if(!result) success = false;
    
      if(!testRelErr<Scalar>(
           "max(*z)",max(*z),"scalarLarge",scalarLarge
           ,"error_tol",error_tol,"warning_tol",warning_tol,out)
        ) success=false;

      if(out) *out << "\nmax(*z,&maxEle,&maxIndex);\n";
      maxEle = ST::zero(); maxIndex = 0;
      max(*z,&maxEle,&maxIndex);
      if(!testRelErr<Scalar>(
           "maxEle",maxEle,"scalarLarge",scalarLarge
           ,"error_tol",error_tol,"warning_tol",warning_tol,out)
        ) success=false;
      result = maxIndex == 1;
      if(out) *out << "\nmaxIndex = " << maxIndex << " == 1 ? " << passfail(result) << std::endl;
      if(!result) success = false;

      if(out) *out << "\nmaxLessThanBound(*z,"<<scalarMedium<<",&maxEle,&maxIndex);\n";
      maxEle = ST::zero(); maxIndex = 0;
      maxLessThanBound(*z,scalarMedium,&maxEle,&maxIndex);
      if(!testRelErr<Scalar>(
           "maxEle",maxEle,"scalarSmall",scalarSmall
           ,"error_tol",error_tol,"warning_tol",warning_tol,out)
        ) success=false;
      result = maxIndex == 0;
      if(out) *out << "\nmaxIndex = " << maxIndex << " == 0 ? " << passfail(result) << std::endl;
      if(!result) success = false;

      if(out) *out << "\nmaxLessThanBound(*z,"<<scalarSmall<<",&maxEle,&maxIndex);\n";
      maxEle = ST::zero(); maxIndex = 0;
      maxLessThanBound(*z,scalarSmall,&maxEle,&maxIndex);
      result = ( maxIndex < 0 );
      if(out) *out << "\nmaxIndex = " << maxIndex << " < 0 ? " << passfail(result) << std::endl;
      if(!result) success = false;
      
      return success;
    }
};

// VectorStdOpsTester


template <class Scalar>
VectorStdOpsTester<Scalar>::VectorStdOpsTester(
  const ScalarMag    &warning_tol
  ,const ScalarMag   &error_tol
  )
  :warning_tol_(warning_tol)
  ,error_tol_(error_tol)
{}

template <class Scalar>
bool VectorStdOpsTester<Scalar>::checkStdOps(
  const VectorSpaceBase<Scalar>    &vecSpc
  ,std::ostream                    *out
  ,const bool                      &dumpAll
  )
{
  using Teuchos::arrayArg;
  using Teuchos::outArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  if(out)
    *out << "\n*** Entering VectorStdOpsTester<"<<ST::name()<<">::checkStdOps(...) ...\n"
         << "using a \'" << vecSpc.description() << "\' object ...\n";

  bool success = true;
  if(out) *out << "\nvecSpc.dim() = " << vecSpc.dim() << std::endl;

  if(out) *out << "\nCreating vectors v1, v2, v3, v4, x and z ...\n";
  Teuchos::RCP<VectorBase<Scalar> >
    v1 = createMember(vecSpc),
    v2 = createMember(vecSpc),
    v3 = createMember(vecSpc),
    v4 = createMember(vecSpc),
    x  = createMember(vecSpc),
    z  = createMember(vecSpc);

  if(out) *out << "\nassign(&*v1,-2.0);\n";
  assign(&*v1,Scalar(-2.0));
  if(out) *out << "\nassign(&*v2,-3.0);\n";
  assign(&*v2,Scalar(-3.0));
  if(out) *out << "\nassign(&*v3,-4.0);\n";
  assign(&*v3,Scalar(-4.0));

  if(!testRelErr<Scalar>(
       "norm_inf(*v1)",norm_inf(*v1),"2.0",Scalar(2.0)
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out)
    ) success=false;

  if(!testRelErr<Scalar>(
       "norm_2(*v1)",norm_2(*v1),"2.0*sqrt(vecSpc.dim())",Scalar(2.0)*ST::squareroot(vecSpc.dim())
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out)
    ) success=false;

  if(!testRelErr<Scalar>(
       "norm_1(*v1)",norm_1(*v1),"2.0*vecSpc.dim()",Scalar(2.0)*Scalar(vecSpc.dim())
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out)
    ) success=false;
  
  if(out) *out << "\nabs(&*z,*v1);\n";
  abs(&*z,*v1);
#ifdef THYRA_VECTOR_STD_OPS_TESTER_DUMP
  SpmdVectorBase<Scalar>::show_dump = true;
  RTOpPack::show_spmd_apply_op_dump = true;
#endif // THYRA_VECTOR_STD_OPS_TESTER_DUMP
  if(!testRelErr<Scalar>(
       "sum(*z)",sum(*z),"2.0*vecSpc.dim()",Scalar(2.0)*Scalar(vecSpc.dim())
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out)
    ) success=false;
#ifdef THYRA_VECTOR_STD_OPS_TESTER_DUMP
  SpmdVectorBase<Scalar>::show_dump = false;
  RTOpPack::show_spmd_apply_op_dump = false;
#endif // THYRA_VECTOR_STD_OPS_TESTER_DUMP
  
  if(out) *out << "\nreciprocal(&*z,*v1);\n";
  reciprocal(&*z,*v1);
  if(!testRelErr<Scalar>(
       "sum(*z)",sum(*z),"-0.5*vecSpc.dim()",Scalar(-0.5)*Scalar(vecSpc.dim())
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out)
    ) success=false;

  if(out) *out << "\nlinear_combination(2,{0.5,0.25},{&*v1,&*v2},0.0,&*z);\n";
  linear_combination(2,arrayArg<Scalar>(0.5,0.25)(),arrayArg<const VectorBase<Scalar>*>(&*v1,&*v2)(),Scalar(0.0),&*z);
  if(!testRelErr<Scalar>(
       "sum(*z)",sum(*z),"(-0.5*2.0-0.25*3.0)*vecSpc.dim()",Scalar(-0.5*2.0-0.25*3.0)*Scalar(vecSpc.dim())
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out)
    ) success=false;

  if(out) *out << "\nassign(&*z,2.0);\n";
  assign(&*z,Scalar(2.0));
  if(out) *out << "\nlinear_combination(3,{0.5,0.25,0.125},{&*v1,&*v2,&*v2},0.5,&*z);\n";
  linear_combination(3,arrayArg<Scalar>(0.5,0.25,0.125)(),arrayArg<const VectorBase<Scalar>*>(&*v1,&*v2,&*v3)(),Scalar(0.5),&*z);
  if(!testRelErr<Scalar>(
       "sum(*z)",sum(*z)
       ,"(0.5*2.0-0.5*2.0-0.25*3.0-0.125*4.0)*vecSpc.dim()",Scalar(0.5*2.0-0.5*2.0-0.25*3.0-0.125*4.0)*Scalar(vecSpc.dim())
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out
       )
    ) success=false;

  if(out) *out << "\nassign(&*z,2.0);\n";
  assign(&*z,Scalar(2.0));
  if(!testRelErr<Scalar>(
       "norm_2(*z,*v2)",norm_2(*z,*v2)
       ,"sqrt(2.0*3.0*3.0*vecSpc.dim())",ST::magnitude(ST::squareroot(Scalar(2.0*3.0*3.0)*Scalar(vecSpc.dim())))
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out
       )
    ) success=false;

  if(!VectorStdOpsTesterComparable<ST::isComparable,Scalar>::checkComparableStdOps(
       vecSpc,&*z,error_tol(),warning_tol(),out,dumpAll)
    ) success=false;

  // ToDo: Add tests for *all* standard operators!
  
  Scalar alpha;
  Scalar beta;

  // Test Vt_S
  if(out) *out << "\nTesting Vt_S(&*z,alpha) ...\n";
  v1  = createMember(vecSpc);
  v2  = createMember(vecSpc);
/* RAB: I commented this out since it was not setting
 * (real,imag) but just the real part.
 */
/*
  if (ST::isComplex)
    alpha = (Scalar(1.2345), Scalar(1.2345));
  else
    alpha = Scalar(1.2345);
*/
  alpha = Scalar(1.2345);
  seed_randomize<Scalar>(12345);
  randomize(Scalar(-ST::one()),ST::one(),&*v1);
  V_V(&*v2,*v1);
  Vt_S(&*v1,alpha);
  Scalar norm_alpha_v1 = norm_2(*v1);
  Scalar alpha_norm_v1 = ST::magnitude(alpha)*norm_2(*v2);
  if(!testMaxErr<Scalar>(
       "norm_alpha_v1 - alpha_norm_v1",ST::magnitude(norm_alpha_v1-alpha_norm_v1)
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out
       )
    ) success=false;
  
  // Test V_StV
  if(out) *out << "\nTesting V_StV(&*z,alpha,*v) ...\n";
  v1  = createMember(vecSpc);
  v2  = createMember(vecSpc);
  z   = createMember(vecSpc);
  alpha = Scalar(-1.2345);
  seed_randomize<Scalar>(12345);
  randomize(Scalar(-ST::one()),ST::one(),&*v1);
  V_StV(&*v2,alpha,*v1);
  Vt_S(&*v1,alpha);
  V_V(&*z,*v1);
  Vp_V(&*z,*v2,Scalar(-ST::one()));
  if(!testMaxErr<Scalar>(
       "norm_2(*z)",norm_2(*z)
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out
       )
    ) success=false;


  // Test Vp_StV
  if(out) *out << "\nTesting Vp_StV(&*z,alpha,*v) ...\n";
  v1  = createMember(vecSpc);
  v2  = createMember(vecSpc);
  v3  = createMember(vecSpc);
  z   = createMember(vecSpc);
  alpha = Scalar(-1.2345);
  seed_randomize<Scalar>(12345);
  randomize(Scalar(-ST::one()),ST::one(),&*v1); // v1 = rand
  randomize(Scalar(-ST::one()),ST::one(),&*v2); // v2 = rand
  V_V(&*v3,*v1); // v3 = v1
  Vp_StV(&*v1,alpha,*v2); // v1 += alpha*v2
  V_StV(&*z,alpha,*v2); // z = alpha*v2
  Vp_V(&*z,*v3); // z += v3
  V_V(&*v3,*v1); // v3 = v1
  Vp_V(&*v3,*z,Scalar(-ST::one())); // v3 -= z
  if(!testMaxErr<Scalar>(
       "norm_2(*v3)",norm_2(*v3)
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out
       )
    ) success=false;
  
  // Test ele_wise_prod
  if(out) *out << "\nTesting ele_wise_prod(alpha,*v1, *v2, &*z) ...\n";
  v1  = createMember(vecSpc);
  v2  = createMember(vecSpc);
  v3  = createMember(vecSpc);
  z   = createMember(vecSpc);
  alpha = Scalar(-1.2345);
  seed_randomize<Scalar>(12345);
  randomize(Scalar(-ST::one()),ST::one(),&*v1); // v1 = rand
  randomize(Scalar(-ST::one()),ST::one(),&*v2); // v2 = rand
  randomize(Scalar(-ST::one()),ST::one(),&*v3); // v3 = rand
  V_V(&*v4, *v1); // v4 = v1
  V_V(&*z, *v2); // z = v2
  ele_wise_prod(alpha, *v2, *v3, &*v1); // v1 += alpha * v2 * v3
  ele_wise_prod_update(alpha, *v3, &*z); // z *= alpha * v3
  Vp_V(&*z, *v4); // z += v4
  V_V(&*v2, *v1); // v2 = v1
  Vp_V(&*v2, *z, Scalar(-ST::one())); // v2 -= z
  if(!testMaxErr<Scalar>(
       "norm_2(*v2)",norm_2(*v2)
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out
       )
    ) success=false;

  // Test Vt_StV
  if(out) *out << "\nTesting Vt_StV(&*z, alpha, *v) ...\n";
  v1  = createMember(vecSpc);
  v2  = createMember(vecSpc);
  v3  = createMember(vecSpc);
  z   = createMember(vecSpc);
  alpha = Scalar(-1.2345);
  seed_randomize<Scalar>(12345);
  randomize(Scalar(-ST::one()),ST::one(),&*v1); // v1 = rand
  randomize(Scalar(-ST::one()),ST::one(),&*v2); // v2 = rand
  V_V(&*v3,*v1); // v3 = v1
  Vt_StV(&*v1,alpha,*v2); // v1 *= alpha*v2
  V_S(&*z,ST::zero()); // z = 0
  Vp_StVtV(&*z,alpha,*v3,*v2); // z += alpha*v3*v2
  V_V(&*v2,*v1); // v2 = v1
  Vp_V(&*v2,*z,Scalar(-ST::one())); // v2 -= z
  if(!testMaxErr<Scalar>(
       "norm_2(*v2)",norm_2(*v2)
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out
       )
    ) success=false;

  // Test V_StVpV
  if(out) *out << "\nTesting V_StVpV(&*z,alpha,*v1,*v2) ...\n";
  v1 = createMember(vecSpc);
  v2 = createMember(vecSpc);
  v3 = createMember(vecSpc);
  x  = createMember(vecSpc);
  z  = createMember(vecSpc);
  alpha = Scalar(1.2345);
  seed_randomize<Scalar>(12345);
  randomize(Scalar(Scalar(-10)*ST::one()),Scalar(Scalar(10)*ST::one()),&*v1);
  randomize(Scalar(Scalar(-10)*ST::one()),Scalar(Scalar(10)*ST::one()),&*v2);
  V_StVpV(&*v3,alpha,*v1,*v2);
  V_V(&*z,*v1);
  Vp_V(&*z,*v2,alpha);
  V_V(&*x,*v3);
  Vp_V(&*x,*z,Scalar(-ST::one()));
  if(!testMaxErr<Scalar>(
       "norm_2(*x)",norm_2(*x)
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out
       )
    ) success=false;

  // Test V_VpStV
  if(out) *out << "\nTesting V_VpStV(&*z,*v1,alpha,*v2) ...\n";
  {
    v1 = createMember(vecSpc);
    v2 = createMember(vecSpc);
    v3 = createMember(vecSpc);
    x  = createMember(vecSpc);
    z  = createMember(vecSpc);
    alpha = Scalar(1.2345);
    seed_randomize<Scalar>(12345);
    randomize(Scalar(Scalar(-10)*ST::one()),Scalar(Scalar(10)*ST::one()),&*v1);
    randomize(Scalar(Scalar(-10)*ST::one()),Scalar(Scalar(10)*ST::one()),&*v2);
    V_VpStV(outArg(*v3), *v1, alpha, *v2);
    V_V(&*z, *v1);
    Vp_StV(&*z, alpha, *v2);
    V_VmV(outArg(*x), *z, *v3);
    if(!testMaxErr<Scalar>(
         "norm_2(*x)",norm_2(*x)
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),out
        )
      ) success=false;
  }
  // Test V_StVpStV
  if(out) *out << "\nTesting V_StVpStV(&*z,alpha,*v1,beta,*v2) ...\n";
  v1 = createMember(vecSpc);
  v2 = createMember(vecSpc);
  v3 = createMember(vecSpc);
  x  = createMember(vecSpc);
  z  = createMember(vecSpc);
  alpha = Scalar(1.2345);
  beta = Scalar(5.4321);
  seed_randomize<Scalar>(12345);
  randomize(Scalar(Scalar(-10)*ST::one()),Scalar(Scalar(10)*ST::one()),&*v1);
  randomize(Scalar(Scalar(-10)*ST::one()),Scalar(Scalar(10)*ST::one()),&*v2);
  V_StVpStV(&*v3,alpha,*v1,beta,*v2);
  V_StV(&*z,alpha,*v1);
  Vp_StV(&*z,beta,*v2);
  V_V(&*x,*v3);
  Vp_V(&*x,*z,Scalar(-ST::one()));
  if(!testMaxErr<Scalar>(
       "norm_2(*x)",norm_2(*x)
       ,"10*error_tol",ScalarMag(ScalarMag(10)*error_tol()),"warning_tol",warning_tol(),out
       )
    ) success=false;

  // Test Vp_V
  if(out) *out << "\nTesting Vp_V(&*v1,*v2,beta) ...\n";
  v1 = createMember(vecSpc);
  v2 = createMember(vecSpc);
  v3 = createMember(vecSpc);
  x  = createMember(vecSpc);
  z  = createMember(vecSpc);
  alpha = Scalar(-2.0);
  beta = Scalar(10.0);
  V_S(&*v1,alpha);
  seed_randomize<Scalar>(12345);
  randomize(Scalar(Scalar(-10)*ST::one()),Scalar(Scalar(10)*ST::one()),&*v2);
  Vp_V(&*v1,*v2,beta); 
  V_S(&*v3,alpha);
  V_StVpV(&*z,beta,*v3,*v2);
  V_StVpV(&*x,Scalar(-ST::one()),*z,*v1);
  if(!testMaxErr<Scalar>(
       "norm_2(*x)",norm_2(*x)
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out
       )
    ) success=false;
  
  // Test Vp_V
  if(out) *out << "\nTesting Vp_V(&*v1,*v2) ...\n";
  v1 = createMember(vecSpc);
  v2 = createMember(vecSpc);
  v3 = createMember(vecSpc);
  x  = createMember(vecSpc);
  z  = createMember(vecSpc);
  alpha = Scalar(-2.0);
  V_S(&*v1,alpha);
  seed_randomize<Scalar>(12345);
  randomize(Scalar(Scalar(-10)*ST::one()),Scalar(Scalar(10)*ST::one()),&*v2);
  Vp_V(&*v1,*v2); 
  V_S(&*v3,alpha);
  V_StVpV(&*z,ST::one(),*v3,*v2);
  V_StVpV(&*x,Scalar(-ST::one()),*z,*v1);
  if(!testMaxErr<Scalar>(
       "norm_2(*x)",norm_2(*x)
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out
       )
    ) success=false;

  // Test V_S
  if(out) *out << "\nTesting V_S(&*v1,alpha) ...\n";
  v1 = createMember(vecSpc);
  v2 = createMember(vecSpc);
  z  = createMember(vecSpc);
  alpha = Scalar(1.2345);
  assign(&*v1,alpha);
  V_S(&*v2,alpha);
  V_StVpV(&*z,Scalar(-ST::one()),*v1,*v2);
  if(!testMaxErr<Scalar>(
       "norm_2(*z)",norm_2(*z)
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out
       )
    ) success=false;

  
  // Test V_V
  if(out) *out << "\nTesting V_V(&*v1,*v2) ...\n";
  v1 = createMember(vecSpc);
  v2 = createMember(vecSpc);
  z  = createMember(vecSpc);
  seed_randomize<Scalar>(12345);
  randomize(Scalar(Scalar(-10)*ST::one()),Scalar(Scalar(10)*ST::one()),&*v1);
  V_V(&*v2,*v1);
  V_StVpV(&*z,Scalar(-ST::one()),*v1,*v2);
  if(!testMaxErr<Scalar>(
       "norm_2(*z)",norm_2(*z)
       ,"error_tol",error_tol(),"warning_tol",warning_tol(),out
       )
    ) success=false;
  


  if(out) *out
    << "\n*** Leaving VectorStdOpsTester<"<<ST::name()<<">::checkStdOps(...) ...\n";

  return success;

}

} // namespace Thyra

#endif // THYRA_VECTOR_STD_OPS_TESTER_HPP
