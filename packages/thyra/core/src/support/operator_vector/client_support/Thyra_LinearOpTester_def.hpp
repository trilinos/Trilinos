// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_LINEAR_OP_TESTER_DEF_HPP
#define THYRA_LINEAR_OP_TESTER_DEF_HPP

#include "Thyra_LinearOpTester_decl.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_describeLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_UniversalMultiVectorRandomizer.hpp"
#include "Teuchos_TestingHelpers.hpp"

namespace Thyra {

template <class Scalar>
class SymmetricLinearOpTester {
 public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  static void checkSymmetry(
      const LinearOpBase<Scalar> &op,
      const Ptr<MultiVectorRandomizerBase<Scalar> > &dRand,
      Teuchos::FancyOStream &testOut,
      const int num_rhs,
      const int num_random_vectors,
      const Teuchos::EVerbosityLevel verbLevel,
      const bool dump_all,
      const ScalarMag &symmetry_error_tol,
      const ScalarMag &symmetry_warning_tol,
      const Ptr<bool> &these_results) {
    using std::endl;

    bool result;
    using Teuchos::OSTab;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar half                          = Scalar(0.4) * ST::one();
    RCP<const VectorSpaceBase<Scalar> > domain = op.domain();

    testOut << endl
            << "op.domain()->isCompatible(*op.range()) == true : ";
    result = op.domain()->isCompatible(*op.range());
    if (!result) *these_results = false;
    testOut << passfail(result) << endl;

    if (result) {
      testOut
          << endl
          << "Checking that the operator is symmetric as:\n"
          << endl
          << "  <0.5*op*v2,v1> == <v2,0.5*op*v1>"
          << endl
          << "   \\_______/            \\_______/"
          << endl
          << "      v4                    v3"
          << endl
          << ""
          << endl
          << "         <v4,v1> == <v2,v3>"
          << endl
          << std::flush;

      for (int rand_vec_i = 1; rand_vec_i <= num_random_vectors; ++rand_vec_i) {
        testOut << endl
                << "Random vector tests = " << rand_vec_i << endl;

        OSTab tab(testOut);

        if (dump_all) testOut << endl
                              << "v1 = randomize(-1,+1); ...\n";
        RCP<MultiVectorBase<Scalar> > v1 = createMembers(domain, num_rhs);
        dRand->randomize(v1.ptr());
        if (dump_all) testOut << endl
                              << "v1 =\n"
                              << describe(*v1, verbLevel);

        if (dump_all) testOut << endl
                              << "v2 = randomize(-1,+1); ...\n";
        RCP<MultiVectorBase<Scalar> > v2 = createMembers(domain, num_rhs);
        dRand->randomize(v2.ptr());
        if (dump_all) testOut << endl
                              << "v2 =\n"
                              << describe(*v2, verbLevel);

        if (dump_all) testOut << endl
                              << "v3 = 0.5*op*v1 ...\n";
        RCP<MultiVectorBase<Scalar> > v3 = createMembers(domain, num_rhs);
        apply(op, NOTRANS, *v1, v3.ptr(), half);
        if (dump_all) testOut << endl
                              << "v3 =\n"
                              << describe(*v3, verbLevel);

        if (dump_all) testOut << endl
                              << "v4 = 0.5*op*v2 ...\n";
        RCP<MultiVectorBase<Scalar> > v4 = createMembers(domain, num_rhs);
        apply(op, NOTRANS, *v2, v4.ptr(), half);
        if (dump_all) testOut << endl
                              << "v4 =\n"
                              << describe(*v4, verbLevel);

        Array<Scalar> prod1(num_rhs), prod2(num_rhs);
        domain->scalarProds(*v4, *v1, prod1());
        domain->scalarProds(*v2, *v3, prod2());

        result = testRelErrors<Scalar, Scalar, ScalarMag>(
            "<v4,v1>", prod1(),
            "<v2,v3>", prod2(),
            "symmetry_error_tol()", symmetry_error_tol,
            "symmetry_warning_tol()", symmetry_warning_tol,
            inOutArg(testOut));
        if (!result) *these_results = false;
      }
    } else {
      testOut << endl
              << "Range and domain spaces are different, skipping check!\n";
    }
  }
};

//
// LinearOpTester
//

template <class Scalar>
LinearOpTester<Scalar>::LinearOpTester()
  : check_linear_properties_(true)
  , linear_properties_warning_tol_(-1.0)
  , linear_properties_error_tol_(-1.0)
  , check_adjoint_(true)
  , adjoint_warning_tol_(-1.0)
  , adjoint_error_tol_(-1.0)
  , check_for_symmetry_(false)
  , symmetry_warning_tol_(-1.0)
  , symmetry_error_tol_(-1.0)
  , num_random_vectors_(1)
  , show_all_tests_(false)
  , dump_all_(false)
  , num_rhs_(1) {
  setDefaultTols();
}

template <class Scalar>
void LinearOpTester<Scalar>::enable_all_tests(const bool enable_all_tests_in) {
  check_linear_properties_ = enable_all_tests_in;
  check_adjoint_           = enable_all_tests_in;
  check_for_symmetry_      = enable_all_tests_in;
}

template <class Scalar>
void LinearOpTester<Scalar>::set_all_warning_tol(const ScalarMag warning_tol_in) {
  linear_properties_warning_tol_ = warning_tol_in;
  adjoint_warning_tol_           = warning_tol_in;
  symmetry_warning_tol_          = warning_tol_in;
}

template <class Scalar>
void LinearOpTester<Scalar>::set_all_error_tol(const ScalarMag error_tol_in) {
  linear_properties_error_tol_ = error_tol_in;
  adjoint_error_tol_           = error_tol_in;
  symmetry_error_tol_          = error_tol_in;
}

template <class Scalar>
bool LinearOpTester<Scalar>::check(
    const LinearOpBase<Scalar> &op,
    const Ptr<MultiVectorRandomizerBase<Scalar> > &rangeRandomizer,
    const Ptr<MultiVectorRandomizerBase<Scalar> > &domainRandomizer,
    const Ptr<Teuchos::FancyOStream> &out_inout) const {
  using std::endl;
  using Teuchos::as;
  using Teuchos::fancyOStream;
  using Teuchos::FancyOStream;
  using Teuchos::inoutArg;
  using Teuchos::OSTab;
  using Teuchos::outArg;
  using Teuchos::rcp;
  using Teuchos::rcpFromPtr;
  using Teuchos::rcpFromRef;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  bool success          = true, result;
  const int loc_num_rhs = this->num_rhs();
  const Scalar r_one    = ST::one();
  const Scalar d_one    = ST::one();
  const Scalar r_half   = as<Scalar>(0.5) * r_one;
  const Scalar d_half   = as<Scalar>(0.5) * d_one;

  RCP<FancyOStream> out;
  if (!is_null(out_inout))
    out = Teuchos::rcpFromPtr(out_inout);
  else
    out = Teuchos::fancyOStream(rcp(new Teuchos::oblackholestream));

  const Teuchos::EVerbosityLevel verbLevel =
      (dump_all() ? Teuchos::VERB_EXTREME : Teuchos::VERB_MEDIUM);

  OSTab tab2(out, 1, "THYRA");

  // ToDo 04/28/2005:
  // * Test the MultiVectorBase apply() function and output to the VectorBase apply() function!

  *out << endl
       << "*** Entering LinearOpTester<" << ST::name() << "," << ST::name() << ">::check(op,...) ...\n";
  if (show_all_tests()) {
    *out << endl
         << "describe op:\n"
         << Teuchos::describe(op, verbLevel);
    /*
      if(op.applyTransposeSupports(CONJ_ELE) && verbLevel==Teuchos::VERB_EXTREME) {
      *out << endl << "describe adjoint op:\n";
      describeLinearOp(*adjoint(Teuchos::rcp(&op,false)),*out,verbLevel);
      }
    */
  } else {
    *out << endl
         << "describe op: " << Teuchos::describe(op, Teuchos::VERB_LOW);
  }

  RCP<MultiVectorRandomizerBase<Scalar> > rRand;
  if (!is_null(rangeRandomizer))
    rRand = rcpFromPtr(rangeRandomizer);
  else
    rRand = universalMultiVectorRandomizer<Scalar>();
  RCP<MultiVectorRandomizerBase<Scalar> > dRand;
  if (!is_null(domainRandomizer))
    dRand = rcpFromPtr(domainRandomizer);
  else
    dRand = universalMultiVectorRandomizer<Scalar>();

  *out << endl
       << "Checking the domain and range spaces ... ";

  RCP<const VectorSpaceBase<Scalar> > range  = op.range();
  RCP<const VectorSpaceBase<Scalar> > domain = op.domain();

  {
    TestResultsPrinter testResultsPrinter(out, show_all_tests());
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();

    bool these_results = true;

    *testOut << endl
             << "op.domain().get() != NULL ? ";
    result = domain.get() != NULL;
    if (!result) these_results = false;
    *testOut << passfail(result) << endl;

    *testOut << endl
             << "op.range().get() != NULL ? ";
    result = range.get() != NULL;
    if (!result) these_results = false;
    *testOut << passfail(result) << endl;

    testResultsPrinter.printTestResults(these_results, inoutArg(success));
  }

  if (check_linear_properties()) {
    *out << endl
         << "this->check_linear_properties()==true:"
         << "Checking the linear properties of the forward linear operator ... ";

    TestResultsPrinter testResultsPrinter(out, show_all_tests());
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();

    bool these_results = true;

    TEUCHOS_TEST_EQUALITY_CONST(op.opSupported(NOTRANS), true, *testOut, these_results);

    if (result) {
      *testOut
          << endl
          << "Checking that the forward operator is truly linear:\n"
          << endl
          << "  0.5*op*(v1 + v2) == 0.5*op*v1 + 0.5*op*v2"
          << endl
          << "          \\_____/         \\___/"
          << endl
          << "             v3            v5"
          << endl
          << "  \\_____________/     \\___________________/"
          << endl
          << "         v4                    v5"
          << endl
          << ""
          << endl
          << "           sum(v4) == sum(v5)"
          << endl
          << std::flush;

      for (int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i) {
        *testOut << endl
                 << "Random vector tests = " << rand_vec_i << endl;

        OSTab tab3(testOut);

        *testOut << endl
                 << "v1 = randomize(-1,+1); ...\n";
        RCP<MultiVectorBase<Scalar> > v1 = createMembers(domain, loc_num_rhs);
        dRand->randomize(v1.ptr());
        if (dump_all()) *testOut << endl
                                 << "v1 =\n"
                                 << describe(*v1, verbLevel);

        *testOut << endl
                 << "v2 = randomize(-1,+1); ...\n";
        RCP<MultiVectorBase<Scalar> > v2 = createMembers(domain, loc_num_rhs);
        dRand->randomize(v2.ptr());
        if (dump_all()) *testOut << endl
                                 << "v2 =\n"
                                 << describe(*v2, verbLevel);

        *testOut << endl
                 << "v3 = v1 + v2 ...\n";
        RCP<MultiVectorBase<Scalar> > v3 = createMembers(domain, loc_num_rhs);
        V_VpV(v3.ptr(), *v1, *v2);
        if (dump_all()) *testOut << endl
                                 << "v3 =\n"
                                 << describe(*v3, verbLevel);

        *testOut << endl
                 << "v4 = 0.5*op*v3 ...\n";
        RCP<MultiVectorBase<Scalar> > v4 = createMembers(range, loc_num_rhs);
        apply(op, NOTRANS, *v3, v4.ptr(), r_half);
        if (dump_all()) *testOut << endl
                                 << "v4 =\n"
                                 << describe(*v4, verbLevel);

        *testOut << endl
                 << "v5 = op*v1 ...\n";
        RCP<MultiVectorBase<Scalar> > v5 = createMembers(range, loc_num_rhs);
        apply(op, NOTRANS, *v1, v5.ptr());
        if (dump_all()) *testOut << endl
                                 << "v5 =\n"
                                 << describe(*v5, verbLevel);

        *testOut << endl
                 << "v5 = 0.5*op*v2 + 0.5*v5 ...\n";
        apply(op, NOTRANS, *v2, v5.ptr(), r_half, r_half);
        if (dump_all()) *testOut << endl
                                 << "v5 =\n"
                                 << describe(*v5, verbLevel);

        Array<Scalar> sum_v4(loc_num_rhs), sum_v5(loc_num_rhs);
        sums(*v4, sum_v4());
        sums(*v5, sum_v5());

        result = testRelErrors<Scalar, Scalar, ScalarMag>(
            "sum(v4)", sum_v4(),
            "sum(v5)", sum_v5(),
            "linear_properties_error_tol()", linear_properties_error_tol(),
            "linear_properties_warning_tol()", linear_properties_warning_tol(),
            testOut.ptr());
        if (!result) these_results = false;
      }
    } else {
      *testOut << endl
               << "Forward operator not supported, skipping check!\n";
    }

    testResultsPrinter.printTestResults(these_results, inoutArg(success));

  } else {
    *out << endl
         << "this->check_linear_properties()==false: Skipping the check of the linear properties of the forward operator!\n";
  }

  if (check_linear_properties() && check_adjoint()) {
    *out << endl
         << "(this->check_linear_properties()&&this->check_adjoint())==true:"
         << " Checking the linear properties of the adjoint operator ... ";

    TestResultsPrinter testResultsPrinter(out, show_all_tests());
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();

    bool these_results = true;

    TEUCHOS_TEST_EQUALITY_CONST(op.opSupported(CONJTRANS), true, *testOut, these_results);

    if (result) {
      *testOut
          << endl
          << "Checking that the adjoint operator is truly linear:\n"
          << endl
          << "  0.5*op'*(v1 + v2) == 0.5*op'*v1 + 0.5*op'*v2"
          << endl
          << "           \\_____/         \\____/"
          << endl
          << "              v3             v5"
          << endl
          << "  \\_______________/    \\_____________________/"
          << endl
          << "         v4                      v5"
          << endl
          << ""
          << endl
          << "           sum(v4) == sum(v5)"
          << endl
          << std::flush;

      for (int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i) {
        *testOut << endl
                 << "Random vector tests = " << rand_vec_i << endl;

        OSTab tab(testOut);

        *testOut << endl
                 << "v1 = randomize(-1,+1); ...\n";
        RCP<MultiVectorBase<Scalar> > v1 = createMembers(range, loc_num_rhs);
        rRand->randomize(v1.ptr());
        if (dump_all()) *testOut << endl
                                 << "v1 =\n"
                                 << describe(*v1, verbLevel);

        *testOut << endl
                 << "v2 = randomize(-1,+1); ...\n";
        RCP<MultiVectorBase<Scalar> > v2 = createMembers(range, loc_num_rhs);
        rRand->randomize(v2.ptr());
        if (dump_all()) *testOut << endl
                                 << "v2 =\n"
                                 << describe(*v2, verbLevel);

        *testOut << endl
                 << "v3 = v1 + v2 ...\n";
        RCP<MultiVectorBase<Scalar> > v3 = createMembers(range, loc_num_rhs);
        V_VpV(v3.ptr(), *v1, *v2);
        if (dump_all()) *testOut << endl
                                 << "v3 =\n"
                                 << describe(*v3, verbLevel);

        *testOut << endl
                 << "v4 = 0.5*op'*v3 ...\n";
        RCP<MultiVectorBase<Scalar> > v4 = createMembers(domain, loc_num_rhs);
        apply(op, CONJTRANS, *v3, v4.ptr(), d_half);
        if (dump_all()) *testOut << endl
                                 << "v4 =\n"
                                 << describe(*v4, verbLevel);

        *testOut << endl
                 << "v5 = op'*v1 ...\n";
        RCP<MultiVectorBase<Scalar> > v5 = createMembers(domain, loc_num_rhs);
        apply(op, CONJTRANS, *v1, v5.ptr());
        if (dump_all()) *testOut << endl
                                 << "v5 =\n"
                                 << describe(*v5, verbLevel);

        *testOut << endl
                 << "v5 = 0.5*op'*v2 + 0.5*v5 ...\n";
        apply(op, CONJTRANS, *v2, v5.ptr(), d_half, d_half);
        if (dump_all()) *testOut << endl
                                 << "v5 =\n"
                                 << describe(*v5, verbLevel);

        Array<Scalar> sum_v4(loc_num_rhs), sum_v5(loc_num_rhs);
        sums(*v4, sum_v4());
        sums(*v5, sum_v5());

        result = testRelErrors<Scalar, Scalar, ScalarMag>(
            "sum(v4)", sum_v4(),
            "sum(v5)", sum_v5(),
            "linear_properties_error_tol()", linear_properties_error_tol(),
            "linear_properties_warning_tol()", linear_properties_warning_tol(),
            testOut.ptr());
        if (!result) these_results = false;
      }
    } else {
      *testOut << endl
               << "Adjoint operator not supported, skipping check!\n";
    }

    testResultsPrinter.printTestResults(these_results, inoutArg(success));

  } else {
    *out << endl
         << "(this->check_linear_properties()&&this->check_adjoint())==false: Skipping the check of the linear properties of the adjoint operator!\n";
  }

  if (check_adjoint()) {
    *out << endl
         << "this->check_adjoint()==true:"
         << " Checking the agreement of the adjoint and forward operators ... ";

    TestResultsPrinter testResultsPrinter(out, show_all_tests());
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();

    bool these_results = true;

    TEUCHOS_TEST_EQUALITY_CONST(op.opSupported(CONJTRANS), true, *testOut, these_results);

    if (result) {
      *testOut
          << endl
          << "Checking that the adjoint agrees with the non-adjoint operator as:\n"
          << endl
          << "  <0.5*op'*v2,v1> == <v2,0.5*op*v1>"
          << endl
          << "   \\________/            \\_______/"
          << endl
          << "       v4                   v3"
          << endl
          << ""
          << endl
          << "         <v4,v1>  == <v2,v3>"
          << endl
          << std::flush;

      for (int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i) {
        *testOut << endl
                 << "Random vector tests = " << rand_vec_i << endl;

        OSTab tab(testOut);

        *testOut << endl
                 << "v1 = randomize(-1,+1); ...\n";
        RCP<MultiVectorBase<Scalar> > v1 = createMembers(domain, loc_num_rhs);
        dRand->randomize(v1.ptr());
        if (dump_all()) *testOut << endl
                                 << "v1 =\n"
                                 << describe(*v1, verbLevel);

        *testOut << endl
                 << "v2 = randomize(-1,+1); ...\n";
        RCP<MultiVectorBase<Scalar> > v2 = createMembers(range, loc_num_rhs);
        rRand->randomize(v2.ptr());
        if (dump_all()) *testOut << endl
                                 << "v2 =\n"
                                 << describe(*v2, verbLevel);

        *testOut << endl
                 << "v3 = 0.5*op*v1 ...\n";
        RCP<MultiVectorBase<Scalar> > v3 = createMembers(range, loc_num_rhs);
        apply(op, NOTRANS, *v1, v3.ptr(), r_half);
        if (dump_all()) *testOut << endl
                                 << "v3 =\n"
                                 << describe(*v3, verbLevel);

        *testOut << endl
                 << "v4 = 0.5*op'*v2 ...\n";
        RCP<MultiVectorBase<Scalar> > v4 = createMembers(domain, loc_num_rhs);
        apply(op, CONJTRANS, *v2, v4.ptr(), d_half);
        if (dump_all()) *testOut << endl
                                 << "v4 =\n"
                                 << describe(*v4, verbLevel);

        Array<Scalar> prod_v4_v1(loc_num_rhs);
        domain->scalarProds(*v4, *v1, prod_v4_v1());
        Array<Scalar> prod_v2_v3(loc_num_rhs);
        range->scalarProds(*v2, *v3, prod_v2_v3());

        result = testRelErrors<Scalar, Scalar, ScalarMag>(
            "<v4,v1>", prod_v4_v1(),
            "<v2,v3>", prod_v2_v3(),
            "adjoint_error_tol()", adjoint_error_tol(),
            "adjoint_warning_tol()", adjoint_warning_tol(),
            testOut.ptr());
        if (!result) these_results = false;
      }
    } else {
      *testOut << endl
               << "Adjoint operator not supported, skipping check!\n";
    }

    testResultsPrinter.printTestResults(these_results, inoutArg(success));

  } else {
    *out << endl
         << "this->check_adjoint()==false:"
         << " Skipping check for the agreement of the adjoint and forward operators!\n";
  }

  if (check_for_symmetry()) {
    *out << endl
         << "this->check_for_symmetry()==true: Performing check of symmetry ... ";

    TestResultsPrinter testResultsPrinter(out, show_all_tests());
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();

    bool these_results = true;

    SymmetricLinearOpTester<Scalar>::checkSymmetry(
        op, dRand.ptr(), *testOut, loc_num_rhs, num_random_vectors(), verbLevel, dump_all(),
        symmetry_error_tol(), symmetry_warning_tol(),
        outArg(these_results));

    testResultsPrinter.printTestResults(these_results, inoutArg(success));

  } else {
    *out << endl
         << "this->check_for_symmetry()==false: Skipping check of symmetry ...\n";
  }

  if (success)
    *out << endl
         << "Congratulations, this LinearOpBase object seems to check out!\n";
  else
    *out << endl
         << "Oh no, at least one of the tests performed with this LinearOpBase object failed (see above failures)!\n";
  *out << endl
       << "*** Leaving LinearOpTester<" << ST::name() << "," << ST::name() << ">::check(...)\n";

  return success;
}

template <class Scalar>
bool LinearOpTester<Scalar>::check(
    const LinearOpBase<Scalar> &op,
    const Ptr<Teuchos::FancyOStream> &out) const {
  using Teuchos::null;
  return check(op, null, null, out);
}

template <class Scalar>
bool LinearOpTester<Scalar>::compare(
    const LinearOpBase<Scalar> &op1,
    const LinearOpBase<Scalar> &op2,
    const Ptr<MultiVectorRandomizerBase<Scalar> > &domainRandomizer,
    const Ptr<Teuchos::FancyOStream> &out_arg) const {
  using std::endl;
  using Teuchos::FancyOStream;
  using Teuchos::inoutArg;
  using Teuchos::OSTab;
  using Teuchos::rcpFromPtr;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  bool success                             = true, result;
  const int loc_num_rhs                    = this->num_rhs();
  const Scalar r_half                      = Scalar(0.5) * ST::one();
  const RCP<FancyOStream> out              = rcpFromPtr(out_arg);
  const Teuchos::EVerbosityLevel verbLevel = (dump_all() ? Teuchos::VERB_EXTREME : Teuchos::VERB_MEDIUM);

  OSTab tab(out, 1, "THYRA");

  if (out.get()) {
    *out
        << endl
        << "*** Entering LinearOpTester<" << ST::name() << "," << ST::name() << ">::compare(op1,op2,...) ...\n";
    if (show_all_tests())
      *out << endl
           << "describe op1:\n"
           << Teuchos::describe(op1, verbLevel);
    else
      *out << endl
           << "describe op1: " << op1.description() << endl;
    if (show_all_tests())
      *out << endl
           << "describe op2:\n"
           << Teuchos::describe(op2, verbLevel);
    else
      *out << endl
           << "describe op2: " << op2.description() << endl;
  }

  RCP<MultiVectorRandomizerBase<Scalar> > dRand;
  if (nonnull(domainRandomizer))
    dRand = rcpFromPtr(domainRandomizer);
  else
    dRand = universalMultiVectorRandomizer<Scalar>();

  RCP<const VectorSpaceBase<Scalar> > range  = op1.range();
  RCP<const VectorSpaceBase<Scalar> > domain = op1.domain();

  if (out.get()) *out << endl
                      << "Checking that range and domain spaces are compatible ... ";

  {
    TestResultsPrinter testResultsPrinter(out, show_all_tests());
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();

    bool these_results = true;

    *testOut << endl
             << "op1.domain()->isCompatible(*op2.domain()) ? ";
    result = op1.domain()->isCompatible(*op2.domain());
    if (!result) these_results = false;
    *testOut << passfail(result) << endl;

    *testOut << endl
             << "op1.range()->isCompatible(*op2.range()) ? ";
    result = op1.range()->isCompatible(*op2.range());
    if (!result) these_results = false;
    *testOut << passfail(result) << endl;

    testResultsPrinter.printTestResults(these_results, inoutArg(success));
  }

  if (!success) {
    if (out.get()) *out << endl
                        << "Skipping further checks since operators are not compatible!\n";
    return success;
  }

  if (out.get()) *out << endl
                      << "Checking that op1 == op2 ... ";

  {
    TestResultsPrinter testResultsPrinter(out, show_all_tests());
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();

    bool these_results = true;

    *testOut
        << endl
        << "Checking that op1 and op2 produce the same results:\n"
        << endl
        << "  0.5*op1*v1 == 0.5*op2*v1"
        << endl
        << "  \\________/    \\________/"
        << endl
        << "      v2            v3"
        << endl
        << ""
        << endl
        << "   norm(v2-v3) ~= 0"
        // << endl << "   |sum(v2)| == |sum(v3)|"
        << endl
        << std::flush;

    for (int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i) {
      *testOut << endl
               << "Random vector tests = " << rand_vec_i << endl;

      OSTab tab2(testOut);

      if (dump_all()) *testOut << endl
                               << "v1 = randomize(-1,+1); ...\n";
      RCP<MultiVectorBase<Scalar> > v1 = createMembers(domain, loc_num_rhs);
      dRand->randomize(v1.ptr());
      if (dump_all()) *testOut << endl
                               << "v1 =\n"
                               << *v1;

      if (dump_all()) *testOut << endl
                               << "v2 = 0.5*op1*v1 ...\n";
      RCP<MultiVectorBase<Scalar> > v2 = createMembers(range, loc_num_rhs);
      apply(op1, NOTRANS, *v1, v2.ptr(), r_half);
      if (dump_all()) *testOut << endl
                               << "v2 =\n"
                               << *v2;

      if (dump_all()) *testOut << endl
                               << "v3 = 0.5*op2*v1 ...\n";
      RCP<MultiVectorBase<Scalar> > v3 = createMembers(range, loc_num_rhs);
      apply(op2, NOTRANS, *v1, v3.ptr(), r_half);
      if (dump_all()) *testOut << endl
                               << "v3 =\n"
                               << *v3;

      // check error in each column
      for (int col_id = 0; col_id < v1->domain()->dim(); col_id++) {
        std::stringstream ss;
        ss << ".col[" << col_id << "]";

        result = Thyra::testRelNormDiffErr(
            "v2" + ss.str(), *v2->col(col_id),
            "v3" + ss.str(), *v3->col(col_id),
            "linear_properties_error_tol()", linear_properties_error_tol(),
            "linear_properties_warning_tol()", linear_properties_warning_tol(),
            &*testOut);
        if (!result) these_results = false;
      }
      /*
      Array<Scalar> sum_v2(loc_num_rhs), sum_v3(loc_num_rhs);
      sums(*v2,&sum_v2[0]);
      sums(*v3,&sum_v3[0]);

      result = testRelErrors(
        loc_num_rhs
        ,"sum(v2)", &sum_v2[0]
        ,"sum(v3)", &sum_v3[0]
        ,"linear_properties_error_tol()", linear_properties_error_tol()
        ,"linear_properties_warning_tol()", linear_properties_warning_tol()
        ,inOutArg(testOut)
        );
      */
      if (!result) these_results = false;
    }

    testResultsPrinter.printTestResults(these_results, inoutArg(success));
  }

  if (out.get()) {
    if (success)
      *out << endl
           << "Congratulations, these two LinearOpBase objects seem to be the same!\n";
    else
      *out << endl
           << "Oh no, these two LinearOpBase objects seem to be different (see above failures)!\n";
    *out << endl
         << "*** Leaving LinearOpTester<" << ST::name() << "," << ST::name() << ">::compare(...)\n";
  }

  return success;
}

template <class Scalar>
bool LinearOpTester<Scalar>::compare(
    const LinearOpBase<Scalar> &op1,
    const LinearOpBase<Scalar> &op2,
    const Ptr<Teuchos::FancyOStream> &out_arg) const {
  return compare(op1, op2, Teuchos::null, out_arg);
}

// private

// Nonmember helper
template <class ScalarMag>
inline void setDefaultTol(const ScalarMag &defaultDefaultTol,
                          ScalarMag &defaultTol) {
  typedef ScalarTraits<ScalarMag> SMT;
  if (defaultTol < SMT::zero()) {
    defaultTol = defaultDefaultTol;
  }
}

template <class Scalar>
void LinearOpTester<Scalar>::setDefaultTols() {
  typedef ScalarTraits<ScalarMag> SMT;
  const ScalarMag defaultWarningTol = 1e+2 * SMT::eps();
  const ScalarMag defaultErrorTol   = 1e+4 * SMT::eps();
  setDefaultTol(defaultWarningTol, linear_properties_warning_tol_);
  setDefaultTol(defaultErrorTol, linear_properties_error_tol_);
  setDefaultTol(defaultWarningTol, adjoint_warning_tol_);
  setDefaultTol(defaultErrorTol, adjoint_error_tol_);
  setDefaultTol(defaultWarningTol, symmetry_warning_tol_);
  setDefaultTol(defaultErrorTol, symmetry_error_tol_);
}

}  // namespace Thyra

#endif  // THYRA_LINEAR_OP_TESTER_DEF_HPP
