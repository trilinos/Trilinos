// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_STD_OPS_TESTER_HPP
#define THYRA_MULTI_VECTOR_STD_OPS_TESTER_HPP

#include "Thyra_MultiVectorStdOpsTester_decl.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"

namespace Thyra {

// MultiVectorStdOpsTester

template <class Scalar>
MultiVectorStdOpsTester<Scalar>::MultiVectorStdOpsTester(
  const ScalarMag    &warning_tol_in
  ,const ScalarMag   &error_tol_in
  ,const int         num_mv_cols_in
  )
  :warning_tol_(warning_tol_in)
  ,error_tol_(error_tol_in)
  ,num_mv_cols_(num_mv_cols_in)
{}

template <class Scalar>
bool MultiVectorStdOpsTester<Scalar>::checkStdOps(
  const VectorSpaceBase<Scalar>    &vecSpc
  ,std::ostream                    *out
  ,const bool                      &/* dumpAll */
  )
{
  using Teuchos::as;
  using Teuchos::ptr;
  using Teuchos::tuple;
  using ST = Teuchos::ScalarTraits<Scalar>;

  if(out)
    *out << "\n*** Entering MultiVectorStdOpsTester<"<<ST::name()<<">::checkStdOps(...) ...\n"
         << "using a \'" << vecSpc.description() << "\' object ...\n";

  bool success = true;
  if(out) *out << "\nvecSpc.dim() = " << vecSpc.dim() << std::endl;

  const Ordinal n = vecSpc.dim();

  const Scalar
    two = as<Scalar>(2.0),
    three = as<Scalar>(3.0),
    four = as<Scalar>(4.0);

  int tc = 0;

  if(out) *out << "\nCreating MultiVectorBase objects V1, V2, and V3 ...\n";
  Teuchos::RCP<MultiVectorBase<Scalar> >
    V1 = createMembers(vecSpc,num_mv_cols()),
    V2 = createMembers(vecSpc,num_mv_cols()),
    V3 = createMembers(vecSpc,num_mv_cols());

  if(out) *out << "\nassign(V1.ptr(),-2.0);\n";
  assign(V1.ptr(),Scalar(-two));

  Teuchos::Array<Scalar> scalars1(num_mv_cols());
  Teuchos::Array<Scalar> scalars2(num_mv_cols());
  Teuchos::Array<ScalarMag> mags1(num_mv_cols());
  Teuchos::Array<ScalarMag> mags2(num_mv_cols());

  // sums
  if(out) *out << "\n"<<tc<<") sums(*V1);\n";
  ++tc;
  {
    sums(*V1, scalars1());
    for (Ordinal i = 0; i < scalars2.size(); ++i)
      scalars2[i] = -two*as<Scalar>(vecSpc.dim());
    if (!testRelErrors<Scalar, Scalar, ScalarMag>(
           "sums(*V1)", scalars1(),
           "-2.0*n", scalars2(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  // norms_1
  if(out) *out << "\n"<<tc<<") norms_1(*V1);\n";
  ++tc;
  {
    norms_1(*V1, mags1());
    for (Ordinal i = 0; i < mags2.size(); ++i)
      mags2[i] = ST::magnitude(two)*as<ScalarMag>(vecSpc.dim());
    if (!testRelErrors<ScalarMag, ScalarMag, ScalarMag>(
           "norms_1(*V1)", mags1(),
           "2.0*n", mags2(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  // norms_2
  if(out) *out << "\n"<<tc<<") norms_2(*V1);\n";
  ++tc;
  {
    norms_2(*V1, mags1());
    for (Ordinal i = 0; i < mags2.size(); ++i)
      mags2[i] = ST::magnitude(two * ST::squareroot(as<Scalar>(n)));
    if (!testRelErrors<ScalarMag, ScalarMag, ScalarMag>(
          "norms_2(*V1)", mags1(),
          "2.0*sqrt(n)", mags2(),
          "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  // norms_inf
  if(out) *out << "\n"<<tc<<") norms_inf(*V1);\n";
  ++tc;
  {
    norms_inf(*V1, mags1());
    for (Ordinal i = 0; i < mags2.size(); ++i)
      mags2[i] = ST::magnitude(two);
    if (!testRelErrors<ScalarMag, ScalarMag, ScalarMag>(
           "norms_inf(*V1)", mags1(),
           "2.0", mags2(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  // assign(scalar)
  if(out) *out << "\n"<<tc<<") assign(V2.ptr(), alpha);\n";
  ++tc;
  {
    assign(V2.ptr(), three);
    norms_2(*V2, mags1());
    for (Ordinal i = 0; i < mags2.size(); ++i)
      mags2[i] = ST::magnitude(three * ST::squareroot(as<Scalar>(n)));
    if (!testRelErrors<ScalarMag, ScalarMag, ScalarMag>(
           "norms_2(*V2)", mags1(),
           "3.0*sqrt(n)", mags2(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  // assign(MV)
  if(out) *out << "\n"<<tc<<") assign(V2.ptr(), *V1);\n";
  ++tc;
  assign(V2.ptr(), *V1);
  norms_2(*V1, mags1());
  norms_2(*V2, mags2());
  if (!testRelErrors<ScalarMag, ScalarMag, ScalarMag>(
         "norms_2(*V1)", mags1(),
         "norms_2(*V2)", mags2(),
         "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
       )
     ) success = false;

  // scale
  if(out) *out << "\n"<<tc<<") scale(alpha,V2.ptr());\n";
  ++tc;
  {
    Scalar alpha = as<Scalar>(1.2345);
    assign(V2.ptr(), *V1);
    scale(alpha, V2.ptr());
    norms_2(*V1, mags1());
    for (Ordinal i = 0; i < mags1.size(); ++i)
      mags1[i] *= ST::magnitude(alpha);
    norms_2(*V2, mags2());
    if (!testRelErrors<ScalarMag, ScalarMag, ScalarMag>(
           "norms_2(alpha*V1)", mags1(),
           "alpha*norms_2(V1)", mags2(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  // scaleUpdate
  if(out) *out << "\n"<<tc<<") scaleUpdate(a,V1,V2.ptr());\n";
  ++tc;
  {
    Teuchos::RCP<VectorBase<Scalar> > a = createMember(vecSpc);
    assign(a.ptr(), two);
    assign(V2.ptr(), four);
    scaleUpdate(*a, *V1, V2.ptr()); // V2(i,j) = 2.0*(-2.0) + 4.0
    norms_2(*V2, mags1());
    if (!testMaxErrors<Scalar>(
           "norms_2(*V2)", mags1(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  // update(alpha, U, V.ptr())
  if(out) *out << "\n"<<tc<<") update(a,V1,V2.ptr());\n";
  ++tc;
  {
    Scalar alpha = as<Scalar>(1.2345);
    assign(V2.ptr(), three);
    assign(V3.ptr(), *V1);
    scale(alpha, V3.ptr());
    Vp_V(V3.ptr(), *V2);
    update(alpha, *V1, V2.ptr());
    norms_2(*V2, mags1());
    norms_2(*V3, mags2());
    if (!testRelErrors<ScalarMag, ScalarMag, ScalarMag>(
           "norms_2(*V2)", mags1(),
           "norms_2(*V3)", mags2(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  // update(alpha, beta, U, V.ptr())
  if(out) *out << "\n"<<tc<<") update(alpha,beta,*V1,V2.ptr());\n";
  ++tc;
  {
    Teuchos::Array<Scalar> alpha(num_mv_cols());
    for (Ordinal i = 0; i < alpha.size(); ++i)
      alpha[i] = as<Scalar>(i+1);
    Scalar beta = as<Scalar>(1.2345);
    assign(V2.ptr(), three);
    assign(V3.ptr(), *V1);
    scale(beta, V3.ptr());
    for (Ordinal i = 0; i < alpha.size(); ++i)
      scale(alpha[i], V3->col(i).ptr());
    Vp_V(V3.ptr(), *V2);
    Teuchos::ArrayView<const Scalar> alphaView = alpha();
    update(alphaView, beta, *V1, V2.ptr());
    norms_2(*V2, mags1());
    norms_2(*V3, mags2());
    if (!testRelErrors<ScalarMag, ScalarMag, ScalarMag>(
           "norms_2(*V2)", mags1(),
           "norms_2(*V3)", mags2(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  // update(U, alpha, beta, V.ptr())
  if(out) *out << "\n"<<tc<<") update(*V1,alpha,beta,V2.ptr());\n";
  ++tc;
  {
    Teuchos::Array<Scalar> alpha(num_mv_cols());
    for (Ordinal i = 0; i < alpha.size(); ++i)
      alpha[i] = as<Scalar>(i+1);
    Scalar beta = as<Scalar>(1.2345);
    assign(V2.ptr(), three);
    assign(V3.ptr(), *V2);
    scale(beta, V3.ptr());
    for (Ordinal i = 0; i < alpha.size(); ++i)
      scale(alpha[i], V3->col(i).ptr());
    Vp_V(V3.ptr(), *V1);
    Teuchos::ArrayView<const Scalar> alphaView = alpha();
    update(*V1, alphaView, beta, V2.ptr());
    norms_2(*V2, mags1());
    norms_2(*V3, mags2());
    if (!testRelErrors<ScalarMag, ScalarMag, ScalarMag>(
           "norms_2(*V2)", mags1(),
           "norms_2(*V3)", mags2(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  // linear_combination
  if(out) *out << "\n"<<tc<<") linear_combination({alpha,beta,gamma},{V1.ptr(),V2.ptr(),V3.ptr()},0.0,V4.ptr());\n";
  ++tc;
  {
    Scalar alpha = two, beta = -three, gamma = three;
    Teuchos::RCP<MultiVectorBase<Scalar> > V4 = createMembers(vecSpc,num_mv_cols());
    assign(V2.ptr(), two);
    assign(V3.ptr(), four);
    linear_combination<Scalar>(
      tuple<Scalar>(alpha, beta, gamma),
      tuple<Ptr<const MultiVectorBase<Scalar> > >(V1.ptr(), V2.ptr(), V3.ptr()),
      as<Scalar>(0.0),
      V4.ptr()); // V4(i,j) = 2.0(-2.0) + (-3.0)(2.0) + 3.0(4.0)
    norms_2(*V4, mags1());
    for (Ordinal i = 0; i < mags2.size(); ++i)
      mags2[i] = ST::magnitude(two * ST::squareroot(as<Scalar>(n)));
    if (!testRelErrors<ScalarMag, ScalarMag, ScalarMag>(
           "norms_2(*V4)", mags1(),
           "2.0*sqrt(n)", mags2(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  if(out) *out << "\n"<<tc<<") linear_combination({alpha,beta,gamma},{V1.ptr(),V2.ptr(),V3.ptr()},0.5,V4.ptr());\n";
  ++tc;
  {
    Scalar alpha = two, beta = -three, gamma = three;
    Teuchos::RCP<MultiVectorBase<Scalar> > V4 = createMembers(vecSpc,num_mv_cols());
    assign(V2.ptr(), two);
    assign(V3.ptr(), four);
    assign(V4.ptr(), -four);
    linear_combination<Scalar>(
      tuple<Scalar>(alpha, beta, gamma),
      tuple<Ptr<const MultiVectorBase<Scalar> > >(V1.ptr(), V2.ptr(), V3.ptr()),
      as<Scalar>(0.5),
      V4.ptr()); // V4(i,j) = 0.5(-4.0) + 2.0(-2.0) + (-3.0)(2.0) + 3.0(4.0)
    norms_2(*V4, mags1());
    if (!testMaxErrors<Scalar>(
           "norms_2(*V4)", mags1(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  // Vt_S
  if(out) *out << "\n"<<tc<<") Vt_S(V1.ptr(),alpha);\n";
  ++tc;
  {
    Scalar alpha = as<Scalar>(1.2345);
    assign(V2.ptr(), *V1);
    Vt_S(V2.ptr(), alpha);
    norms_2(*V1, mags1());
    for (Ordinal i = 0; i < mags1.size(); ++i)
      mags1[i] *= ST::magnitude(alpha);
    norms_2(*V2, mags2());
    if (!testRelErrors<ScalarMag, ScalarMag, ScalarMag>(
           "norms_2(alpha*V1)", mags1(),
           "alpha*norms_2(V1)", mags2(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  // Vp_S
  if(out) *out << "\n"<<tc<<") Vp_S(V2.ptr(),alpha);\n";
  ++tc;
  {
    assign(V2.ptr(), *V1);
    Vp_S(V2.ptr(), two); // V2(i,j) = -2.0 + 2.0
    norms_2(*V2, mags1());
    if (!testMaxErrors<Scalar>(
           "norms_2(V2)", mags1(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }

  // Vp_V
  if(out) *out << "\n"<<tc<<") Vp_V(V2.ptr(),*V1);\n";
  ++tc;
  {
    assign(V2.ptr(), two);
    Vp_V(V2.ptr(), *V1); // V2(i,j) = -2.0 + 2.0
    norms_2(*V2, mags1());
    if (!testMaxErrors<Scalar>(
           "norms_2(V2)", mags1(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }


  // V_VpV
  if(out) *out << "\n"<<tc<<") V_VpV(V3.ptr(),*V1,*V2);\n";
  ++tc;
  {
    assign(V2.ptr(), two);
    V_VpV(V3.ptr(), *V1, *V2); // V3(i,j) = -2.0 + 2.0
    norms_2(*V3, mags1());
    if (!testMaxErrors<Scalar>(
           "norms_2(V3)", mags1(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }


  // V_VmV
  if(out) *out << "\n"<<tc<<") V_VmV(V3.ptr(),*V1,*V2);\n";
  ++tc;
  {
    assign(V2.ptr(), -two);
    V_VmV(V3.ptr(), *V1, *V2); // V3(i,j) = -2.0 - (-2.0)
    norms_2(*V3, mags1());
    if (!testMaxErrors<Scalar>(
           "norms_2(V3)", mags1(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }


  // V_StVpV
  if(out) *out << "\n"<<tc<<") V_StVpV(V3.ptr(),alpha,*V1,*V2);\n";
  ++tc;
  {
    Scalar alpha = as<Scalar>(1.2345);
    assign(V2.ptr(), three);
    V_StVpV(V3.ptr(), alpha, *V1, *V2);
    scale(alpha, V1.ptr());
    Vp_V(V2.ptr(), *V1);
    V_VmV(V3.ptr(), *V2, *V3);
    norms_2(*V3, mags1());
    if (!testMaxErrors<Scalar>(
           "norms_2(V3)", mags1(),
           "error_tol", error_tol(), "warning_tol", warning_tol(), ptr(out)
         )
       ) success = false;
  }


  if(out) *out
    << "\n*** Leaving MultiVectorStdOpsTester<"<<ST::name()<<">::checkStdOps(...) ...\n";

  return success;

}

} // namespace Thyra

#endif // THYRA_MULTI_VECTOR_STD_OPS_TESTER_HPP
