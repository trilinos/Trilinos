// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

// For debugging
#include "RTOpPack_SPMD_apply_op.hpp"


/** \brief Main test driver function for testing various composite linear operator classes
 */
template <class Scalar>
bool run_composite_linear_ops_tests(
  const Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > comm,
  const int n,
  const bool useSpmd,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType &tol,
  const bool dumpAll,
  Teuchos::FancyOStream *out_arg
  )
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType    ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> STM;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::dyn_cast;
  using Teuchos::OSTab;
  using Thyra::relErr;
  using Thyra::passfail;

  RCP<Teuchos::FancyOStream>
    out = rcp(new Teuchos::FancyOStream(rcp(out_arg,false)));

  const Teuchos::EVerbosityLevel
    verbLevel = dumpAll?Teuchos::VERB_EXTREME:Teuchos::VERB_HIGH;

  if (nonnull(out)) *out
    << "\n*** Entering run_composite_linear_ops_tests<"<<ST::name()<<">(...) ...\n";

  bool success = true, result;

  const ScalarMag warning_tol = ScalarMag(1e-2)*tol, error_tol = tol;
  Thyra::LinearOpTester<Scalar> linearOpTester;
  linearOpTester.linear_properties_warning_tol(warning_tol);
  linearOpTester.linear_properties_error_tol(error_tol);
  linearOpTester.adjoint_warning_tol(warning_tol);
  linearOpTester.adjoint_error_tol(error_tol);
  linearOpTester.dump_all(dumpAll);
  Thyra::LinearOpTester<Scalar> symLinearOpTester(linearOpTester);
  symLinearOpTester.check_for_symmetry(true);
  symLinearOpTester.symmetry_warning_tol(STM::squareroot(warning_tol));
  symLinearOpTester.symmetry_error_tol(STM::squareroot(error_tol));

  RCP<const Thyra::VectorSpaceBase<Scalar> > space;
  if(useSpmd) space = Thyra::defaultSpmdVectorSpace<Scalar>(comm,n,-1);
  else space = Thyra::defaultSpmdVectorSpace<Scalar>(n);
  if (nonnull(out)) *out
    << "\nUsing a basic vector space described as " << describe(*space,verbLevel) << " ...\n";

  if (nonnull(out)) *out << "\nCreating random n x (n/2) multi-vector origA ...\n";
  RCP<Thyra::MultiVectorBase<Scalar> >
    mvOrigA = createMembers(space,n/2,"origA");
  Thyra::seed_randomize<Scalar>(0);
  //RTOpPack::show_spmd_apply_op_dump = true;
  Thyra::randomize( as<Scalar>(as<Scalar>(-1)*ST::one()), as<Scalar>(as<Scalar>(+1)*ST::one()),
    mvOrigA.ptr() );
  RCP<const Thyra::LinearOpBase<Scalar> >
    origA = mvOrigA;
  if (nonnull(out)) *out << "\norigA =\n" << describe(*origA,verbLevel);
  //RTOpPack::show_spmd_apply_op_dump = false;

  if (nonnull(out)) *out << "\nTesting origA ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*origA, out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out
    << "\nCreating implicit scaled linear operator A1 = scale(0.5,origA) ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A1 = scale(as<Scalar>(0.5),origA);
  if (nonnull(out)) *out << "\nA1 =\n" << describe(*A1,verbLevel);

  if (nonnull(out)) *out << "\nTesting A1 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A1,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nTesting that A1.getOp() == origA ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.compare(
    *dyn_cast<const Thyra::DefaultScaledAdjointLinearOp<Scalar> >(*A1).getOp(),
    *origA,out.ptr());
  if(!result) success = false;

  {

    if (nonnull(out)) *out
      << "\nUnwrapping origA to get non-persisting pointer to origA_1, scalar and transp ...\n";
    Scalar  scalar;
    Thyra::EOpTransp transp;
    const Thyra::LinearOpBase<Scalar> *origA_1 = NULL;
    unwrap( *origA, &scalar, &transp, &origA_1 );
    TEUCHOS_TEST_FOR_EXCEPT( origA_1 == NULL );

    if (nonnull(out)) *out << "\nscalar = " << scalar << " == 1 ? ";
    result = (scalar == ST::one());
    if(!result) success = false;
    if (nonnull(out))	*out << passfail(result) << std::endl;

    if (nonnull(out)) *out << "\ntransp = " << toString(transp) << " == NOTRANS ? ";
    result = (transp == Thyra::NOTRANS);
    if(!result) success = false;
    if (nonnull(out))	*out << passfail(result) << std::endl;

    if (nonnull(out)) *out << "\nTesting that origA_1 == origA ...\n";
    Thyra::seed_randomize<Scalar>(0);
    result = linearOpTester.compare(*origA_1,*origA,out.ptr());
    if(!result) success = false;

  }

  {

    if (nonnull(out)) *out << "\nUnwrapping A1 to get non-persisting pointer to origA_2 ...\n";
    Scalar  scalar;
    Thyra::EOpTransp transp;
    const Thyra::LinearOpBase<Scalar> *origA_2 = NULL;
    unwrap( *A1, &scalar, &transp, &origA_2 );
    TEUCHOS_TEST_FOR_EXCEPT( origA_2 == NULL );

    if (nonnull(out)) *out << "\nscalar = " << scalar << " == 0.5 ? ";
    result = (scalar == as<Scalar>(0.5));
    if(!result) success = false;
    if (nonnull(out))	*out << passfail(result) << std::endl;

    if (nonnull(out)) *out << "\ntransp = " << toString(transp) << " == NOTRANS ? ";
    result = (transp == Thyra::NOTRANS);
    if(!result) success = false;
    if (nonnull(out))	*out << passfail(result) << std::endl;

    if (nonnull(out)) *out << "\nTesting that origA_2 == origA ...\n";
    Thyra::seed_randomize<Scalar>(0);
    result = linearOpTester.compare(*origA_2,*origA,out.ptr());
    if(!result) success = false;

  }

  if (nonnull(out)) *out << "\nCreating implicit scaled linear operator A2 = adjoint(A1) ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A2 = adjoint(A1);
  if (nonnull(out)) *out << "\nA2 =\n" << describe(*A2,verbLevel);

  if (nonnull(out)) *out << "\nTesting A2 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A2,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nTesting that A2.getOp() == A1 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.compare(*dyn_cast<const Thyra::DefaultScaledAdjointLinearOp<Scalar> >(*A2).getOp(),*A1,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nCreating implicit scaled, adjoined linear operator A3 = adjoint(scale(2.0,(A2)) ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A3 = adjoint(scale(as<Scalar>(2.0),A2));
  if (nonnull(out)) *out << "\nA3 =\n" << describe(*A3,verbLevel);

  if (nonnull(out)) *out << "\nTesting A3 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A3,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nTesting that A3 == origA ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.compare(*A3,*origA,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nCalling all of the rest of the functions for non-const just to test them ...\n";
  RCP<Thyra::LinearOpBase<Scalar> >
    A4 = nonconstScale(
      as<Scalar>(0.25)
      ,nonconstAdjoint(
        nonconstTranspose(
          nonconstAdjoint(
            nonconstScaleAndAdjoint(
              as<Scalar>(4.0)
              ,Thyra::TRANS
              ,Teuchos::rcp_const_cast<Thyra::LinearOpBase<Scalar> >(origA)
              )
            )
          )
        )
      );
  if(!ST::isComplex) A4 = nonconstTranspose(nonconstAdjoint(A4)); // Should result in CONJ
  if (nonnull(out)) *out << "\nA4 =\n" << describe(*A4,verbLevel);

  if (nonnull(out)) *out << "\nTesting A4 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A4,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nCalling all of the rest of the functions for const just to test them ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A5 = scale(
      as<Scalar>(0.25)
      ,adjoint(
        transpose(
          adjoint(
            scaleAndAdjoint(
              as<Scalar>(4.0)
              ,Thyra::TRANS
              ,origA
              )
            )
          )
        )
      );
  if(!ST::isComplex) A5 = transpose(adjoint(A5)); // Should result in CONJ
  if (nonnull(out)) *out << "\nA5 =\n" << describe(*A5,verbLevel);

  if (nonnull(out)) *out << "\nTesting A5 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A5,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nCreating a multiplied operator A6 = origA^H*A1 ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A6 = multiply(adjoint(origA),A1);
  if (nonnull(out)) *out << "\nA6 =\n" << describe(*A6,verbLevel);

  if (nonnull(out)) *out << "\nTesting A6 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A6,out.ptr());
  if(!result) success = false;
  // Note that testing the symmetry above helps to check the transpose mode
  // against the non-transpose mode!

#ifdef TEUCHOS_DEBUG
  if (nonnull(out)) *out << "\nCreating an invalid multiplied operator A6b = origA*origA (should throw an exception) ...\n\n";
  try {
    RCP<const Thyra::LinearOpBase<Scalar> >
      A6b = multiply(origA,origA);
    result = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,result)
  if (nonnull(out))
    *out << "\nCaught expected exception : " << (result?"failed\n":"passed\n");
  if(result) success = false;
#endif // TEUCHOS_DEBUG

  if (nonnull(out)) *out << "\nCreating a non-const multiplied operator A7 = origA^H*A1 ...\n";
  RCP<Thyra::LinearOpBase<Scalar> >
    A7 = nonconstMultiply(
      rcp_const_cast<Thyra::LinearOpBase<Scalar> >(adjoint(origA))
      ,rcp_const_cast<Thyra::LinearOpBase<Scalar> >(A1)
      );
  if (nonnull(out)) *out << "\nA7 =\n" << describe(*A7,verbLevel);

  if (nonnull(out)) *out << "\nTesting A7 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A7,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nCreating an added operator A8 = origA + A1 ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A8 = add(origA,A1);
  if (nonnull(out)) *out << "\nA8 =\n" << describe(*A8,verbLevel);

  if (nonnull(out)) *out << "\nTesting A8 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A8,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nCreating a symmetric subtracted operator A8b = A6 + adjoint(origA)*origA ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A8b = subtract(A6,multiply(adjoint(origA),origA));
  if (nonnull(out)) *out << "\nA8b =\n" << describe(*A8b,verbLevel);

  if (nonnull(out)) *out << "\nTesting A8b ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A8b,out.ptr());
  if(!result) success = false;

#ifdef TEUCHOS_DEBUG
  if (nonnull(out)) *out << "\nCreating an invalid added operator A8c = origA + adjoint(origA) (should throw an exception) ...\n\n";
  try {
    RCP<const Thyra::LinearOpBase<Scalar> >
      A8c = add(origA,adjoint(origA));
    result = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,result)
  if (nonnull(out))
    *out << "\nCaught expected exception : " << (result?"failed\n":"passed\n");
  if(result) success = false;
#endif // TEUCHOS_DEBUG

  RCP<const Thyra::LinearOpBase<Scalar> >
    nullOp = null;

  if (nonnull(out)) *out << "\nCreating a blocked 2x2 linear operator A9 = [ A6, A1^H; A1, null ] ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A9 = Thyra::block2x2<Scalar>(
      A6,  adjoint(A1)
      ,A1, nullOp
      );
  if (nonnull(out)) *out << "\nA9 =\n" << describe(*A9,verbLevel);

  if (nonnull(out)) *out << "\nTesting A9 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A9,out.ptr());
  if(!result) success = false;
  // Note that testing the symmetry above helps to check the transpose mode
  // against the non-transpose mode!

  if (nonnull(out)) *out << "\nCreating a blocked 2x2 linear operator A9_a = [ A6, A1^H; A1, null ] using pre-formed range and domain product spaces ...\n";
  RCP<Thyra::PhysicallyBlockedLinearOpBase<Scalar> >
    A9_a = rcp(new Thyra::DefaultBlockedLinearOp<Scalar>());
  A9_a->beginBlockFill(
    rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<Scalar> >(A9,true)->productRange()
    ,rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<Scalar> >(A9,true)->productDomain()
    );
  A9_a->setBlock(0,0,A6);
  A9_a->setBlock(0,1,adjoint(A1));
  A9_a->setBlock(1,0,A1);
  A9_a->endBlockFill();
  if (nonnull(out)) *out << "\nA9_a =\n" << describe(*A9_a,verbLevel);

  if (nonnull(out)) *out << "\nTesting A9_a ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A9_a,out.ptr());
  if(!result) success = false;
  // Note that testing the symmetry above helps to check the transpose mode
  // against the non-transpose mode!

  if (nonnull(out)) *out << "\nComparing A9 == A9_a ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.compare(*A9,*A9_a,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nCreating a blocked 2x2 linear operator A9_b = [ A6, A1^H; A1, null ] using flexible fill ...\n";
  RCP<Thyra::PhysicallyBlockedLinearOpBase<Scalar> >
    A9_b = rcp(new Thyra::DefaultBlockedLinearOp<Scalar>());
  A9_b->beginBlockFill();
  A9_b->setBlock(0,0,A6);
  A9_b->setBlock(0,1,adjoint(A1));
  A9_b->setBlock(1,0,A1);
  A9_b->endBlockFill();
  if (nonnull(out)) *out << "\nA9_b =\n" << describe(*A9_b,verbLevel);

  if (nonnull(out)) *out << "\nTesting A9_b ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A9_b,out.ptr());
  if(!result) success = false;
  // Note that testing the symmetry above helps to check the transpose mode
  // against the non-transpose mode!

  if (nonnull(out)) *out << "\nComparing A9 == A9_b ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.compare(*A9,*A9_b,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nCreating a blocked 2x2 linear operator A9a = [ null, A1^H; A1, null ] ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A9a = Thyra::block2x2<Scalar>(
      nullOp,  adjoint(A1),
      A1,      nullOp
      );
  if (nonnull(out)) *out << "\nA9a =\n" << describe(*A9a,verbLevel);

  if (nonnull(out)) *out << "\nTesting A9a ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A9a,out.ptr());
  if(!result) success = false;
  // Note that testing the symmetry above helps to check the transpose mode
  // against the non-transpose mode!

#ifdef TEUCHOS_DEBUG
  if (nonnull(out)) *out << "\nCreating an invalid blocked 2x2 operator A9b = [ A6, A1^H; A1, A1 ] (should throw an exception) ...\n\n";
  try {
    RCP<const Thyra::LinearOpBase<Scalar> >
      A9b = Thyra::block2x2<Scalar>(
        A6,  adjoint(A1),
        A1,  A1
        );
    result = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,result)
  if (nonnull(out))
    *out << "\nCaught expected exception : " << (result?"failed\n":"passed\n");
  if(result) success = false;
#endif // TEUCHOS_DEBUG

#ifdef TEUCHOS_DEBUG
  if (nonnull(out)) *out << "\nCreating an invalid blocked 2x2 operator A9c = [ A1, A1 ; null, null ] (should throw an exception) ...\n\n";
  try {
    RCP<const Thyra::LinearOpBase<Scalar> >
      A9c = Thyra::block2x2<Scalar>(
        A1,       A1,
        nullOp,   nullOp
        );
    result = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,result)
  if (nonnull(out))
    *out << "\nCaught expected exception : " << (result?"failed\n":"passed\n");
  if(result) success = false;
#endif // TEUCHOS_DEBUG

#ifdef TEUCHOS_DEBUG
  if (nonnull(out)) *out << "\nCreating an invalid blocked 2x2 operator A9d = [ A1, null; A1, null ] (should throw an exception) ...\n\n";
  try {
    RCP<const Thyra::LinearOpBase<Scalar> >
      A9d = Thyra::block2x2<Scalar>(
        A1,  nullOp,
        A1,  nullOp
        );
    result = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,result)
  if (nonnull(out))
    *out << "\nCaught expected exception : " << (result?"failed\n":"passed\n");
  if(result) success = false;
#endif // TEUCHOS_DEBUG

  if (nonnull(out)) *out << "\nCreating a blocked 2x1 linear operator A10 = [ A6; A1 ] ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A10 = Thyra::block2x1<Scalar>(
      A6,
      A1
      );
  if (nonnull(out)) *out << "\nA10 =\n" << describe(*A10,verbLevel);

  if (nonnull(out)) *out << "\nTesting A10 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A10,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nCreating a blocked 1x2 linear operator A11 = [ A9, A10 ] ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A11 = Thyra::block1x2<Scalar>( A9, A10 );
  if (nonnull(out)) *out << "\nA11 =\n" << describe(*A11,verbLevel);

  if (nonnull(out)) *out << "\nTesting A11 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A11,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nCreating a zero linear operator A12 = 0 (range and domain spaces of origA) ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A12 = Thyra::zero(origA->range(),origA->domain());
  if (nonnull(out)) *out << "\nA12 =\n" << describe(*A12,verbLevel);

  if (nonnull(out)) *out << "\nTesting A12 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A12,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nCreating a blocked 2x2 linear operator A13 = [ zero, A1^H; A1, zero ] ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A13 = Thyra::block2x2<Scalar>(
      Thyra::zero(A1->domain(),A1->domain()),   adjoint(A1),
      A1,                                       Thyra::zero(A1->range(),A1->range())
      );
  if (nonnull(out)) *out << "\nA13 =\n" << describe(*A13,verbLevel);

  if (nonnull(out)) *out << "\nComparing A9a == A13 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.compare(*A9a,*A13,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\nCreating a zero linear operator A14 = I (range space of origA) ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A14 = Thyra::identity(origA->range());
  if (nonnull(out)) *out << "\nA14 =\n" << describe(*A14,verbLevel);

  if (nonnull(out)) *out << "\nTesting A14 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A14,out.ptr());
  if(!result) success = false;

  if (nonnull(out)) *out << "\n*** Leaving run_composite_linear_ops_tests<"<<ST::name()<<">(...) ...\n";

  return success;

} // end run_composite_linear_ops_tests() [Doxygen looks for this!]

int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {
    const Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> >
      comm = Teuchos::DefaultComm<Thyra::Ordinal>::getComm();

    //
    // Read options from command-line
    //

    int n         = 4;
    bool useSpmd   = true;
    bool dumpAll  = false;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "local-dim", &n, "Local number of elements in each constituent vector." );
    clp.setOption( "use-spmd", "use-serial", &useSpmd, "Determines if MPI or serial vector space is used." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    //
    // Run the tests
    //

#if defined(HAVE_TEUCHOS_INST_FLOAT) && defined(HAVE_TEUCHOS_BLASFLOAT)
    if( !run_composite_linear_ops_tests<float>(comm,n,useSpmd,float(1e-4),dumpAll,verbose?&*out:NULL) ) success = false;
#endif
    if( !run_composite_linear_ops_tests<double>(comm,n,useSpmd,double(1e-12),dumpAll,verbose?&*out:NULL) ) success = false;
#if defined(HAVE_TEUCHOS_INST_COMPLEX_FLOAT) && defined(HAVE_TEUCHOS_INST_FLOAT) && defined(HAVE_TEUCHOS_BLASFLOAT)
    if( !run_composite_linear_ops_tests<std::complex<float> >(comm,n,useSpmd,float(1e-4),dumpAll,verbose?&*out:NULL) ) success = false;
#endif
#if defined(HAVE_TEUCHOS_INST_COMPLEX_DOUBLE)
    if( !run_composite_linear_ops_tests<std::complex<double> >(comm,n,useSpmd,double(1e-12),dumpAll,verbose?&*out:NULL) ) success = false;
#endif
#if defined(HAVE_TEUCHOS_GNU_MP) && !defined(USE_MPI) // mpf_class can not be used with MPI yet!
    if( !run_composite_linear_ops_tests<mpf_class>(comm,n,useSpmd,mpf_class(1e-12),dumpAll,verbose?&*out:NULL) ) success = false;
#endif

    if( verbose ) {
      if(success) *out << "\nAll of the tests seem to have run successfully!\n";
      else        *out << "\nOh no! at least one of the tests failed!\n";
    }

  } // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)

  return success ? EXIT_SUCCESS : EXIT_FAILURE;

} // end main() [Doxygen looks for this!]
