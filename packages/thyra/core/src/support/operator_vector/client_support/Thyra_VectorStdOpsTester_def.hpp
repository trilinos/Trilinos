// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_VECTOR_STD_OPS_TESTER_HPP
#define THYRA_VECTOR_STD_OPS_TESTER_HPP

#include "Thyra_VectorStdOpsTester_decl.hpp"
#include "Thyra_TestingTools.hpp"
#include "RTOpPack_TOpSetAssendingValues.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_Assert.hpp"

//#define THYRA_VECTOR_STD_OPS_TESTER_DUMP

#ifdef THYRA_VECTOR_STD_OPS_TESTER_DUMP
#  include "RTOpPack_SPMD_apply_op.hpp"
#endif // THYRA_VECTOR_STD_OPS_TESTER_DUMP


namespace Thyra {


// VectorStdOpsTesterComparable (using partial specialization to only do tests in some cases)


template <bool isComparable, class Scalar>
class VectorStdOpsTesterComparable {
public:
  static bool checkComparableStdOps(
    const VectorSpaceBase<Scalar>                                 &vecSpc
    ,const Ptr<VectorBase<Scalar> >                               &z
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
    ,const Ptr<VectorBase<Scalar> >                               &z
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
    ,const Ptr<VectorBase<Scalar> >                               &z
    ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &error_tol
    ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &warning_tol
    ,std::ostream                                                 *out
    ,const bool                                                   &dumpAll
    )
    {
      typedef Teuchos::ScalarTraits<Scalar> ST;
      using Teuchos::outArg;

      bool success = true, result;
      
      if(out) *out << "\nTesting comparable operations ...\n";
      
      const Scalar scalarSmall(1e-5), scalarMedium(2.0), scalarLarge(100.0);
      if(out) *out << "\nassign(z.ptr(),"<<scalarMedium<<");\n";
      assign(z.ptr(),Scalar(scalarMedium));
      if(out && dumpAll) *out << "\nz =\n" << *z;
      if(out) *out << "\nset_ele(0,"<<scalarSmall<<",z.ptr());\n";
      set_ele(0,scalarSmall,z.ptr());
      if(out && dumpAll) *out << "\nz =\n" << *z;
      if(out) *out << "\nset_ele(1,"<<scalarLarge<<",z.ptr());\n";
      set_ele(1,scalarLarge,z.ptr());
      if(out && dumpAll) *out << "\nz =\n" << *z;
      if(out) *out << "\nset_ele(vecSpc.dim()-2,"<<scalarSmall<<",z.ptr());\n";
      set_ele(vecSpc.dim()-2,scalarSmall,z.ptr());
      if(out && dumpAll) *out << "\nz =\n" << *z;
      if(out) *out << "\nset_ele(vecSpc.dim()-1,"<<scalarLarge<<",z.ptr());\n";
      set_ele(vecSpc.dim()-1,scalarLarge,z.ptr());
      if(out && dumpAll) *out << "\nz =\n" << *z;

      Scalar minEle; Ordinal minIndex;
      Scalar maxEle; Ordinal maxIndex;

      if(!testRelErr<Scalar>(
           "min(*z)",min(*z),"scalarSmall",scalarSmall
           ,"error_tol",error_tol,"warning_tol",warning_tol,out
           )
        ) success=false;

      if(out) *out << "\nmin(*z,&minEle,&minIndex);\n";
      minEle = ST::zero(); minIndex = 0;
      min(*z, outArg(minEle), outArg(minIndex));
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
      minGreaterThanBound(*z, scalarMedium, outArg(minEle), outArg(minIndex));
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
      minGreaterThanBound(*z,scalarLarge, outArg(minEle), outArg(minIndex));
      result = minIndex < 0;
      if(out) *out << "\nminIndex = " << minIndex << " < 0 ? " << passfail(result) << std::endl;
      if(!result) success = false;
    
      if(!testRelErr<Scalar>(
           "max(*z)",max(*z),"scalarLarge",scalarLarge
           ,"error_tol",error_tol,"warning_tol",warning_tol,out)
        ) success=false;

      if(out) *out << "\nmax(*z,&maxEle,&maxIndex);\n";
      maxEle = ST::zero(); maxIndex = 0;
      max(*z, outArg(maxEle), outArg(maxIndex));
      if(!testRelErr<Scalar>(
           "maxEle",maxEle,"scalarLarge",scalarLarge
           ,"error_tol",error_tol,"warning_tol",warning_tol,out)
        ) success=false;
      result = maxIndex == 1;
      if(out) *out << "\nmaxIndex = " << maxIndex << " == 1 ? " << passfail(result) << std::endl;
      if(!result) success = false;

      if(out) *out << "\nmaxLessThanBound(*z,"<<scalarMedium<<",&maxEle,&maxIndex);\n";
      maxEle = ST::zero(); maxIndex = 0;
      maxLessThanBound(*z, scalarMedium, outArg(maxEle), outArg(maxIndex));
      if(!testRelErr<Scalar>(
           "maxEle",maxEle,"scalarSmall",scalarSmall
           ,"error_tol",error_tol,"warning_tol",warning_tol,out)
        ) success=false;
      result = maxIndex == 0;
      if(out) *out << "\nmaxIndex = " << maxIndex << " == 0 ? " << passfail(result) << std::endl;
      if(!result) success = false;

      if(out) *out << "\nmaxLessThanBound(*z,"<<scalarSmall<<",&maxEle,&maxIndex);\n";
      maxEle = ST::zero(); maxIndex = 0;
      maxLessThanBound(*z, scalarSmall, outArg(maxEle), outArg(maxIndex));
      result = ( maxIndex < 0 );
      if(out) *out << "\nmaxIndex = " << maxIndex << " < 0 ? " << passfail(result) << std::endl;
      if(!result) success = false;
      
      return success;
    }
};


// Other helpers


template<class Scalar>
void setEleTestCase( const Ptr<VectorBase<Scalar> > &z, const Ordinal i, int &tc,
  std::ostream &out, bool &success)
{
  using Teuchos::as;
  out << "\n"<<tc<<") set_ele(z, "<<i<<");\n";
  ++tc;
  {
    typedef ScalarTraits<Scalar> ST;
    const Scalar val_i = as<Scalar>(i+1);
    assign<Scalar>(z, ST::zero());
    set_ele(i, val_i, z);
    TEUCHOS_TEST_EQUALITY_CONST(get_ele(*z, i), val_i, out, success);
    TEUCHOS_TEST_EQUALITY_CONST(sum(*z), val_i, out, success);
  }
}


// VectorStdOpsTester


template <class Scalar>
VectorStdOpsTester<Scalar>::VectorStdOpsTester(
  const ScalarMag &warning_tol_in,
  const ScalarMag &error_tol_in
  )
  :warning_tol_(warning_tol_in),
   error_tol_(error_tol_in)
{}


template <class Scalar>
bool VectorStdOpsTester<Scalar>::checkStdOps(
  const VectorSpaceBase<Scalar> &vecSpc,
  std::ostream *out_out,
  const bool &dumpAll
  )
{
  using Teuchos::as;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::outArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEUCHOS_ASSERT(out_out);
  std::ostream &out = *out_out;

  out << "\n*** Entering VectorStdOpsTester<"<<ST::name()<<">::checkStdOps(...) ...\n"
      << "using a \'" << vecSpc.description() << "\' object ...\n";

  bool success = true;
  out << "\nvecSpc.dim() = " << vecSpc.dim() << std::endl;

  const Ordinal n = vecSpc.dim();

  TEST_FOR_EXCEPTION( n < 4, std::logic_error,
    "Error: n = "<<n<<" must be least 4 or greater to"
    " run Thyra::VectorStdOpsTester::checkStdOps(...)!" );

  const Scalar
    two = as<Scalar>(2.0),
    three = as<Scalar>(3.0),
    four = as<Scalar>(4.0);

  int tc = 0;

  out << "\nCreating vectors v1, v2, v3, v4, x and z ...\n";
  Teuchos::RCP<VectorBase<Scalar> >
    v1 = createMember(vecSpc),
    v2 = createMember(vecSpc),
    v3 = createMember(vecSpc),
    v4 = createMember(vecSpc),
    y = createMember(vecSpc),
    x = createMember(vecSpc),
    z = createMember(vecSpc);

  out << "\nassign(v1.ptr(), -2.0);\n";
  assign<Scalar>(v1.ptr(), -two);
  out << "\nassign(v2.ptr(), -3.0);\n";
  assign<Scalar>(v2.ptr(), -three);
  out << "\nassign(v3.ptr(), -4.0);\n";
  assign<Scalar>(v3.ptr(), -four);
  out << "\ny[i] = i+1\n";
  {
    RTOpPack::TOpSetAssendingValues<Scalar> setAssendOp(ST::zero());
    applyOp<Scalar>( setAssendOp,
      ArrayView<const Ptr<const VectorBase<Scalar> > >(null),
      tuple<Ptr<VectorBase<Scalar> > >(y.ptr())(),
      null );
  }

  // sum
  out << "\n"<<tc<<") sum(*y);\n";
  ++tc;
  TEUCHOS_TEST_ASSERT(
    testRelErr<Scalar>(
      "sum(*y)", sum(*y),
      "0.5*(n+1)*n", as<Scalar>(0.5*(n+1)*n),
      "error_tol", error_tol(), "warning_tol", warning_tol(), &out),
    out, success);

  // norm_inf
  out << "\n"<<tc<<") nom_inf(*v1);\n";
  ++tc;
  TEUCHOS_TEST_ASSERT(
    testRelErr<Scalar>(
      "norm_inf(*v1)", norm_inf(*v1),
      "2.0", two,
      "error_tol", error_tol(), "warning_tol", warning_tol(), &out),
    out, success);

  // norm_2
  out << "\n"<<tc<<") norm_2(*v1);\n";
  ++tc;
  TEUCHOS_TEST_ASSERT(
    testRelErr<Scalar>(
      "norm_2(*v1)", norm_2(*v1),
      "2.0*sqrt(vecSpc.dim())", as<Scalar>(2.0)*ST::squareroot(vecSpc.dim()),
      "error_tol", error_tol(), "warning_tol",warning_tol(), &out),
    out, success);

  // norm_1
  out << "\n"<<tc<<") norm_1(*v1);\n";
  ++tc;
  TEUCHOS_TEST_ASSERT(
    testRelErr<Scalar>(
      "norm_1(*v1)" ,norm_1(*v1),
      "2.0*vecSpc.dim()", as<Scalar>(2.0)*as<Scalar>(vecSpc.dim()),
      "error_tol", error_tol(), "warning_tol", warning_tol(), &out),
    out, success);
  
  // abs
  out << "\n"<<tc<<") abs(z.ptr(),*v1);\n";
  ++tc;
  {
    abs(z.ptr(), *v1);
    if(!testRelErr<Scalar>(
         "sum(*z)",sum(*z),"2.0*vecSpc.dim()",as<Scalar>(2.0)*as<Scalar>(vecSpc.dim())
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out)
      ) success=false;
  }

  // get_ele

  out << "\n"<<tc<<") val = get_ele(y, 0);\n";
  ++tc;
  {
    const Scalar val = get_ele<Scalar>(*y, 0);
    TEUCHOS_TEST_EQUALITY_CONST( val, as<Scalar>(1), out, success );
  }

  out << "\n"<<tc<<") val = get_ele<Scalar>(*y, 1);\n";
  ++tc;
  {
    const Scalar val = get_ele<Scalar>(*y, 1);
    TEUCHOS_TEST_EQUALITY_CONST( val, as<Scalar>(2), out, success );
  }

  out << "\n"<<tc<<") val = get_ele<Scalar>(*y, n-2);\n";
  ++tc;
  {
    const Scalar val = get_ele<Scalar>(*y, n-2);
    TEUCHOS_TEST_EQUALITY_CONST( val, as<Scalar>(n-1), out, success );
  }

  out << "\n"<<tc<<") val = get_ele<Scalar>(*y, n-1);\n";
  ++tc;
  {
    const Scalar val = get_ele<Scalar>(*y, n-1);
    TEUCHOS_TEST_EQUALITY_CONST( val, as<Scalar>(n), out, success );
  }

#ifdef THYRA_DEBUG

  out << "\n"<<tc<<") get_ele<Scalar>(*y, -1);\n";
  ++tc;
  {
    TEUCHOS_TEST_THROW(get_ele<Scalar>(*y, -1), std::out_of_range, out, success );
  }

  out << "\n"<<tc<<") get_ele<Scalar>(*y, n);\n";
  ++tc;
  {
    TEUCHOS_TEST_THROW(get_ele<Scalar>(*y, n), std::out_of_range, out, success );
  }

#endif // THYRA_DEBUG

  // set_ele

  setEleTestCase<Scalar>(z.ptr(), 0, tc, out, success);

  setEleTestCase<Scalar>(z.ptr(), 1, tc, out, success);

  setEleTestCase<Scalar>(z.ptr(), n-2, tc, out, success);

  setEleTestCase<Scalar>(z.ptr(), n-1, tc, out, success);

#ifdef THYRA_DEBUG

  TEUCHOS_TEST_THROW(set_ele(-1, two, z.ptr()),
    std::out_of_range, out, success);

  TEUCHOS_TEST_THROW(set_ele(n, two, z.ptr()),
    std::out_of_range, out, success);

#endif // THYRA_DEBUG
    
  // reciprocal
  out << "\n"<<tc<<") reciprocal(z.ptr(),*v1);\n";
  ++tc;
  {
    reciprocal(z.ptr(), *v1);
    if(!testRelErr<Scalar>(
         "sum(*z)",sum(*z),"-0.5*vecSpc.dim()",as<Scalar>(-0.5)*as<Scalar>(vecSpc.dim())
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out)
      ) success=false;
  }
    
  // linear_combination

  out << "\n"<<tc<<") linear_combination(2,{0.5,0.25},{v1.ptr(),v2.ptr()},0.0,z.ptr());\n";
  ++tc;
  {
    linear_combination<Scalar>(
      tuple<Scalar>(0.5, 0.25)(),
      tuple<Ptr<const VectorBase<Scalar> > >(v1.ptr(), v2.ptr())(),
      ST::zero(),
      z.ptr());
    TEUCHOS_TEST_ASSERT(
      testRelErr<Scalar>(
        "sum(*z)", sum(*z),
        "(-0.5*2.0-0.25*3.0)*vecSpc.dim()", as<Scalar>((-0.5 * 2.0 - 0.25 * 3.0) *vecSpc.dim()),
        "error_tol", error_tol(), "warning_tol", warning_tol(), &out),
      out, success);
  }

  out << "\nassign(z.ptr(), 2.0);\n";
  ++tc;
  assign(z.ptr(), as<Scalar>(2.0));

  out << "\n"<<tc<<") linear_combination(3,{0.5,0.25,0.125},{v1.ptr(),v2.ptr(),v2.ptr()},0.5,z.ptr());\n";
  ++tc;
  {
    linear_combination<Scalar>(
      tuple<Scalar>(0.5, 0.25, 0.125)(),
      tuple<Ptr<const VectorBase<Scalar> > >(v1.ptr(), v2.ptr(),  v3.ptr())(),
      as<Scalar>(0.5),
      z.ptr());
    if(!testRelErr<Scalar>(
         "sum(*z)", sum(*z),
         "(0.5*2.0-0.5*2.0-0.25*3.0-0.125*4.0)*vecSpc.dim()",
         as<Scalar>(0.5*2.0-0.5*2.0-0.25*3.0-0.125*4.0)*as<Scalar>(vecSpc.dim()),
         "error_tol", error_tol(), "warning_tol", warning_tol(), &out
        )
      ) success=false;
  }
  
  // assgin
  out << "\n"<<tc<<") assign(z.ptr(),2.0);\n";
  ++tc;
  {
    assign(z.ptr(),as<Scalar>(2.0));
    if(!testRelErr<Scalar>(
         "norm_2(*z,*v2)",norm_2(*z,*v2)
         ,"sqrt(2.0*3.0*3.0*vecSpc.dim())",ST::magnitude(ST::squareroot(as<Scalar>(2.0*3.0*3.0)*as<Scalar>(vecSpc.dim())))
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out
        )
      ) success=false;

    if(!VectorStdOpsTesterComparable<ST::isComparable,Scalar>::checkComparableStdOps(
         vecSpc,z.ptr(),error_tol(),warning_tol(),&out,dumpAll)
      ) success=false;
  }

  // Test Vt_S
  out << "\n"<<tc<<") Testing Vt_S(z.ptr(),alpha) ...\n";
  ++tc;
  {
#ifdef THYRA_VECTOR_STD_OPS_TESTER_DUMP
    RTOpPack::show_spmd_apply_op_dump = true;
#endif // THYRA_VECTOR_STD_OPS_TESTER_DUMP
    v1  = createMember(vecSpc);
    v2  = createMember(vecSpc);
    const Scalar alpha = as<Scalar>(1.2345);
    seed_randomize<Scalar>(12345);
    randomize(as<Scalar>(-ST::one()),ST::one(),v1.ptr());
    V_V(v2.ptr(),*v1);
    Vt_S(v1.ptr(), alpha);
    const Scalar norm_alpha_v1 = norm_2(*v1);
    //out << "norm_alpha_v1 = " << norm_alpha_v1 << "\n";
    const Scalar mag_alpha = ST::magnitude(alpha);
    //out << "mag_alpha = " << mag_alpha << "\n";
    const Scalar norm_2_v2 = norm_2(*v2);
    //out << "norm_2_v2 = " << norm_2_v2 << "\n";
    const Scalar alpha_norm_v1 = mag_alpha * norm_2_v2;
    //out << "alpha_norm_v1 = " << alpha_norm_v1 << "\n";
    if(!testMaxErr<Scalar>(
         "norm_alpha_v1 - alpha_norm_v1",ST::magnitude(norm_alpha_v1-alpha_norm_v1)
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out
        )
      ) success=false;
#ifdef THYRA_VECTOR_STD_OPS_TESTER_DUMP
    RTOpPack::show_spmd_apply_op_dump = false;
#endif // THYRA_VECTOR_STD_OPS_TESTER_DUMP
  }
  
  // Test V_StV
  out << "\n"<<tc<<") Testing V_StV(z.ptr(),alpha,*v) ...\n";
  ++tc;
  {
    v1  = createMember(vecSpc);
    v2  = createMember(vecSpc);
    z   = createMember(vecSpc);
    const Scalar alpha = as<Scalar>(-1.2345);
    seed_randomize<Scalar>(12345);
    randomize(as<Scalar>(-ST::one()),ST::one(),v1.ptr());
    V_StV(v2.ptr(),alpha,*v1);
    Vt_S(v1.ptr(),alpha);
    V_V(z.ptr(),*v1);
    Vp_V(z.ptr(),*v2,as<Scalar>(-ST::one()));
    if(!testMaxErr<Scalar>(
         "norm_2(*z)",norm_2(*z)
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out
        )
      ) success=false;
  }

  // Test Vp_StV
  out << "\n"<<tc<<") Testing Vp_StV(z.ptr(),alpha,*v) ...\n";
  ++tc;
  {
    v1  = createMember(vecSpc);
    v2  = createMember(vecSpc);
    v3  = createMember(vecSpc);
    z   = createMember(vecSpc);
    const Scalar alpha = as<Scalar>(-1.2345);
    seed_randomize<Scalar>(12345);
    randomize(as<Scalar>(-ST::one()),ST::one(),v1.ptr()); // v1 = rand
    randomize(as<Scalar>(-ST::one()),ST::one(),v2.ptr()); // v2 = rand
    V_V(v3.ptr(),*v1); // v3 = v1
    Vp_StV(v1.ptr(),alpha,*v2); // v1 += alpha*v2
    V_StV(z.ptr(),alpha,*v2); // z = alpha*v2
    Vp_V(z.ptr(),*v3); // z += v3
    V_V(v3.ptr(),*v1); // v3 = v1
    Vp_V(v3.ptr(),*z,as<Scalar>(-ST::one())); // v3 -= z
    if(!testMaxErr<Scalar>(
         "norm_2(*v3)",norm_2(*v3)
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out
        )
      ) success=false;
  }
  
  // Test ele_wise_prod
  out << "\n"<<tc<<") Testing ele_wise_prod(alpha,*v1, *v2, z.ptr()) ...\n";
  ++tc;
  {
    v1  = createMember(vecSpc);
    v2  = createMember(vecSpc);
    v3  = createMember(vecSpc);
    z   = createMember(vecSpc);
    const Scalar alpha = as<Scalar>(-1.2345);
    seed_randomize<Scalar>(12345);
    randomize(as<Scalar>(-ST::one()),ST::one(),v1.ptr()); // v1 = rand
    randomize(as<Scalar>(-ST::one()),ST::one(),v2.ptr()); // v2 = rand
    randomize(as<Scalar>(-ST::one()),ST::one(),v3.ptr()); // v3 = rand
    V_V(v4.ptr(), *v1); // v4 = v1
    V_V(z.ptr(), *v2); // z = v2
    ele_wise_prod(alpha, *v2, *v3, v1.ptr()); // v1 += alpha * v2 * v3
    ele_wise_prod_update(alpha, *v3, z.ptr()); // z *= alpha * v3
    Vp_V(z.ptr(), *v4); // z += v4
    V_V(v2.ptr(), *v1); // v2 = v1
    Vp_V(v2.ptr(), *z, as<Scalar>(-ST::one())); // v2 -= z
    if(!testMaxErr<Scalar>(
         "norm_2(*v2)",norm_2(*v2)
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out
        )
      ) success=false;
  }
  
  // Test ele_wise_scale
  out << "\n"<<tc<<") Testing ele_wise_scale(*v1, z.ptr()) ...\n";
  ++tc;
  {
    v1  = createMember(vecSpc);
    z   = createMember(vecSpc);
    V_S(v1.ptr(), as<Scalar>(2.0));
    V_S(z.ptr(), as<Scalar>(3.0));
    ele_wise_scale( *v1, z.ptr() );
    if (!testRelErr(
        "norm_2(*z)", norm_2(*z),
        "ST::squareroot(n*sqr(3.0*2.0))", ST::squareroot(n*36.0),
        "error_tol", error_tol(),
        "warning_tol", warning_tol(),
        &out
        )
      ) success=false;
  }

  // Test Vt_StV
  out << "\n"<<tc<<") Testing Vt_StV(z.ptr(), alpha, *v) ...\n";
  ++tc;
  {
    v1  = createMember(vecSpc);
    v2  = createMember(vecSpc);
    v3  = createMember(vecSpc);
    z   = createMember(vecSpc);
    const Scalar alpha = as<Scalar>(-1.2345);
    seed_randomize<Scalar>(12345);
    randomize(as<Scalar>(-ST::one()),ST::one(),v1.ptr()); // v1 = rand
    randomize(as<Scalar>(-ST::one()),ST::one(),v2.ptr()); // v2 = rand
    V_V(v3.ptr(),*v1); // v3 = v1
    Vt_StV(v1.ptr(),alpha,*v2); // v1 *= alpha*v2
    V_S(z.ptr(),ST::zero()); // z = 0
    Vp_StVtV(z.ptr(),alpha,*v3,*v2); // z += alpha*v3*v2
    V_V(v2.ptr(),*v1); // v2 = v1
    Vp_V(v2.ptr(),*z,as<Scalar>(-ST::one())); // v2 -= z
    if(!testMaxErr<Scalar>(
         "norm_2(*v2)",norm_2(*v2)
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out
        )
      ) success=false;
  }

  // Test V_StVpV
  out << "\n"<<tc<<") Testing V_StVpV(z.ptr(),alpha,*v1,*v2) ...\n";
  ++tc;
  {
    v1 = createMember(vecSpc);
    v2 = createMember(vecSpc);
    v3 = createMember(vecSpc);
    x  = createMember(vecSpc);
    z  = createMember(vecSpc);
    const Scalar alpha = as<Scalar>(1.2345);
    seed_randomize<Scalar>(12345);
    randomize(as<Scalar>(as<Scalar>(-10)*ST::one()),as<Scalar>(as<Scalar>(10)*ST::one()),v1.ptr());
    randomize(as<Scalar>(as<Scalar>(-10)*ST::one()),as<Scalar>(as<Scalar>(10)*ST::one()),v2.ptr());
    V_StVpV(v3.ptr(),alpha,*v1,*v2);
    V_V(z.ptr(),*v1);
    Vp_V(z.ptr(),*v2,alpha);
    V_V(x.ptr(),*v3);
    Vp_V(x.ptr(),*z,as<Scalar>(-ST::one()));
    if(!testMaxErr<Scalar>(
         "norm_2(*x)",norm_2(*x)
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out
        )
      ) success=false;
  }

  // Test V_VpStV
  out << "\n"<<tc<<") Testing V_VpStV(z.ptr(),*v1,alpha,*v2) ...\n";
  ++tc;
  {
    v1 = createMember(vecSpc);
    v2 = createMember(vecSpc);
    v3 = createMember(vecSpc);
    x  = createMember(vecSpc);
    z  = createMember(vecSpc);
    const Scalar alpha = as<Scalar>(1.2345);
    seed_randomize<Scalar>(12345);
    randomize<Scalar>(as<Scalar>(as<Scalar>(-10)*ST::one()),
      as<Scalar>(as<Scalar>(10)*ST::one()), v1.ptr());
    randomize<Scalar>(as<Scalar>(as<Scalar>(-10)*ST::one()),
      as<Scalar>(as<Scalar>(10)*ST::one()), v2.ptr());
    V_VpStV(outArg(*v3), *v1, alpha, *v2);
    V_V(z.ptr(), *v1);
    Vp_StV(z.ptr(), alpha, *v2);
    V_VmV(outArg(*x), *z, *v3);
    if(!testMaxErr<Scalar>(
         "norm_2(*x)",norm_2(*x)
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out
        )
      ) success=false;
  }

  // Test V_StVpStV
  out << "\n"<<tc<<") Testing V_StVpStV(z.ptr(),alpha,*v1,beta,*v2) ...\n";
  ++tc;
  {
    v1 = createMember(vecSpc);
    v2 = createMember(vecSpc);
    v3 = createMember(vecSpc);
    x  = createMember(vecSpc);
    z  = createMember(vecSpc);
    const Scalar alpha = as<Scalar>(1.2345);
    const Scalar beta = as<Scalar>(5.4321);
    seed_randomize<Scalar>(12345);
    randomize(as<Scalar>(as<Scalar>(-10)*ST::one()),as<Scalar>(as<Scalar>(10)*ST::one()),v1.ptr());
    randomize(as<Scalar>(as<Scalar>(-10)*ST::one()),as<Scalar>(as<Scalar>(10)*ST::one()),v2.ptr());
    V_StVpStV(v3.ptr(),alpha,*v1,beta,*v2);
    V_StV(z.ptr(),alpha,*v1);
    Vp_StV(z.ptr(),beta,*v2);
    V_V(x.ptr(),*v3);
    Vp_V(x.ptr(),*z,as<Scalar>(-ST::one()));
    if(!testMaxErr<Scalar>(
         "norm_2(*x)",norm_2(*x)
         ,"10*error_tol",ScalarMag(ScalarMag(10)*error_tol()),"warning_tol",warning_tol(),&out
        )
      ) success=false;
  }

  // Test Vp_V
  out << "\n"<<tc<<") Testing Vp_V(v1.ptr(),*v2,beta) ...\n";
  ++tc;
  {
    v1 = createMember(vecSpc);
    v2 = createMember(vecSpc);
    v3 = createMember(vecSpc);
    x  = createMember(vecSpc);
    z  = createMember(vecSpc);
    const Scalar alpha = as<Scalar>(-2.0);
    const Scalar beta = as<Scalar>(10.0);
    V_S(v1.ptr(),alpha);
    seed_randomize<Scalar>(12345);
    randomize(as<Scalar>(as<Scalar>(-10)*ST::one()),as<Scalar>(as<Scalar>(10)*ST::one()),v2.ptr());
    Vp_V(v1.ptr(),*v2,beta); 
    V_S(v3.ptr(),alpha);
    V_StVpV(z.ptr(),beta,*v3,*v2);
    V_StVpV(x.ptr(),as<Scalar>(-ST::one()),*z,*v1);
    if(!testMaxErr<Scalar>(
         "norm_2(*x)",norm_2(*x)
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out
        )
      ) success=false;
  }
  
  // Test Vp_V
  out << "\n"<<tc<<") Testing Vp_V(v1.ptr(),*v2) ...\n";
  ++tc;
  {
    v1 = createMember(vecSpc);
    v2 = createMember(vecSpc);
    v3 = createMember(vecSpc);
    x  = createMember(vecSpc);
    z  = createMember(vecSpc);
    const Scalar alpha = as<Scalar>(-2.0);
    V_S(v1.ptr(),alpha);
    seed_randomize<Scalar>(12345);
    randomize(as<Scalar>(as<Scalar>(-10)*ST::one()),as<Scalar>(as<Scalar>(10)*ST::one()),v2.ptr());
    Vp_V(v1.ptr(),*v2); 
    V_S(v3.ptr(),alpha);
    V_StVpV(z.ptr(),ST::one(),*v3,*v2);
    V_StVpV(x.ptr(),as<Scalar>(-ST::one()),*z,*v1);
    if(!testMaxErr<Scalar>(
         "norm_2(*x)",norm_2(*x)
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out
        )
      ) success=false;
  }

  // Test V_S
  out << "\n"<<tc<<") Testing V_S(v1.ptr(),alpha) ...\n";
  ++tc;
  {
    v1 = createMember(vecSpc);
    v2 = createMember(vecSpc);
    z  = createMember(vecSpc);
    const Scalar alpha = as<Scalar>(1.2345);
    assign(v1.ptr(),alpha);
    V_S(v2.ptr(),alpha);
    V_StVpV(z.ptr(),as<Scalar>(-ST::one()),*v1,*v2);
    if(!testMaxErr<Scalar>(
         "norm_2(*z)",norm_2(*z)
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out
        )
      ) success=false;
  }
  
  // Test V_V
  out << "\n"<<tc<<") Testing V_V(v1.ptr(),*v2) ...\n";
  ++tc;
  {
    v1 = createMember(vecSpc);
    v2 = createMember(vecSpc);
    z  = createMember(vecSpc);
    seed_randomize<Scalar>(12345);
    randomize(as<Scalar>(as<Scalar>(-10)*ST::one()),as<Scalar>(as<Scalar>(10)*ST::one()),v1.ptr());
    V_V(v2.ptr(),*v1);
    V_StVpV(z.ptr(),as<Scalar>(-ST::one()),*v1,*v2);
    if(!testMaxErr<Scalar>(
         "norm_2(*z)",norm_2(*z)
         ,"error_tol",error_tol(),"warning_tol",warning_tol(),&out
        )
      ) success=false;
  }

  // ToDo: Add tests for *all* standard operators!

  out << "\n*** Leaving VectorStdOpsTester<"<<ST::name()<<">::checkStdOps(...) ...\n";

  return success;


}


} // namespace Thyra


#endif // THYRA_VECTOR_STD_OPS_TESTER_HPP
