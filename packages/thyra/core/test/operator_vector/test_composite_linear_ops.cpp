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
  const Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> >   comm
  ,const int                                                       n
  ,const bool                                                      useSpmd
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType     &tol
  ,const bool                                                      dumpAll
  ,Teuchos::FancyOStream                                           *out_arg
  )
{

  using Thyra::relErr;
  using Thyra::passfail;
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

  RCP<Teuchos::FancyOStream>
    out = rcp(new Teuchos::FancyOStream(rcp(out_arg,false)));

  const Teuchos::EVerbosityLevel
    verbLevel = dumpAll?Teuchos::VERB_EXTREME:Teuchos::VERB_HIGH;

  if(out.get()) *out << "\n*** Entering run_composite_linear_ops_tests<"<<ST::name()<<">(...) ...\n";

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
  if(out.get()) *out << "\nUsing a basic vector space described as " << describe(*space,verbLevel) << " ...\n";
  
  if(out.get()) *out << "\nCreating random n x (n/2) multi-vector origA ...\n";
  RCP<Thyra::MultiVectorBase<Scalar> >
    mvOrigA = createMembers(space,n/2,"origA");
  Thyra::seed_randomize<Scalar>(0);
  //RTOpPack::show_spmd_apply_op_dump = true;
  Thyra::randomize( Scalar(Scalar(-1)*ST::one()), Scalar(Scalar(+1)*ST::one()), &*mvOrigA );
  RCP<const Thyra::LinearOpBase<Scalar> >
    origA = mvOrigA;
  if(out.get()) *out << "\norigA =\n" << describe(*origA,verbLevel);
  //RTOpPack::show_spmd_apply_op_dump = false;

  if(out.get()) *out << "\nTesting origA ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*origA,out.get());
  if(!result) success = false;
  
  if(out.get()) *out << "\nCreating implicit scaled linear operator A1 = scale(0.5,origA) ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A1 = scale(Scalar(0.5),origA);
  if(out.get()) *out << "\nA1 =\n" << describe(*A1,verbLevel);

  if(out.get()) *out << "\nTesting A1 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A1,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nTesting that A1.getOp() == origA ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.compare(*dyn_cast<const Thyra::DefaultScaledAdjointLinearOp<Scalar> >(*A1).getOp(),*origA,out.get());
  if(!result) success = false;

  {

    if(out.get()) *out << "\nUnwrapping origA to get non-persisting pointer to origA_1, scalar and transp ...\n";
    Scalar  scalar;
    Thyra::EOpTransp transp;
    const Thyra::LinearOpBase<Scalar> *origA_1 = NULL;
    unwrap( *origA, &scalar, &transp, &origA_1 );
    TEST_FOR_EXCEPT( origA_1 == NULL );

    if(out.get()) *out << "\nscalar = " << scalar << " == 1 ? ";
    result = (scalar == ST::one());
    if(!result) success = false;
    if(out.get())	*out << passfail(result) << std::endl;

    if(out.get()) *out << "\ntransp = " << toString(transp) << " == NOTRANS ? ";
    result = (transp == Thyra::NOTRANS);
    if(!result) success = false;
    if(out.get())	*out << passfail(result) << std::endl;
    
    if(out.get()) *out << "\nTesting that origA_1 == origA ...\n";
    Thyra::seed_randomize<Scalar>(0);
    result = linearOpTester.compare(*origA_1,*origA,out.get());
    if(!result) success = false;
    
  }

  {

    if(out.get()) *out << "\nUnwrapping A1 to get non-persisting pointer to origA_2 ...\n";
    Scalar  scalar;
    Thyra::EOpTransp transp;
    const Thyra::LinearOpBase<Scalar> *origA_2 = NULL;
    unwrap( *A1, &scalar, &transp, &origA_2 );
    TEST_FOR_EXCEPT( origA_2 == NULL );

    if(out.get()) *out << "\nscalar = " << scalar << " == 0.5 ? ";
    result = (scalar == Scalar(0.5));
    if(!result) success = false;
    if(out.get())	*out << passfail(result) << std::endl;

    if(out.get()) *out << "\ntransp = " << toString(transp) << " == NOTRANS ? ";
    result = (transp == Thyra::NOTRANS);
    if(!result) success = false;
    if(out.get())	*out << passfail(result) << std::endl;

    if(out.get()) *out << "\nTesting that origA_2 == origA ...\n";
    Thyra::seed_randomize<Scalar>(0);
    result = linearOpTester.compare(*origA_2,*origA,out.get());
    if(!result) success = false;

  }
  
  if(out.get()) *out << "\nCreating implicit scaled linear operator A2 = adjoint(A1) ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A2 = adjoint(A1);
  if(out.get()) *out << "\nA2 =\n" << describe(*A2,verbLevel);

  if(out.get()) *out << "\nTesting A2 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A2,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nTesting that A2.getOp() == A1 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.compare(*dyn_cast<const Thyra::DefaultScaledAdjointLinearOp<Scalar> >(*A2).getOp(),*A1,out.get());
  if(!result) success = false;
  
  if(out.get()) *out << "\nCreating implicit scaled, adjoined linear operator A3 = adjoint(scale(2.0,(A2)) ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A3 = adjoint(scale(Scalar(2.0),A2));
  if(out.get()) *out << "\nA3 =\n" << describe(*A3,verbLevel);

  if(out.get()) *out << "\nTesting A3 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A3,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nTesting that A3 == origA ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.compare(*A3,*origA,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nCalling all of the rest of the functions for non-const just to test them ...\n";
  RCP<Thyra::LinearOpBase<Scalar> >
    A4 = nonconstScale(
      Scalar(0.25)
      ,nonconstAdjoint(
        nonconstTranspose(
          nonconstAdjoint(
            nonconstScaleAndAdjoint(
              Scalar(4.0)
              ,Thyra::TRANS
              ,Teuchos::rcp_const_cast<Thyra::LinearOpBase<Scalar> >(origA)
              )
            )
          )
        )
      );
  if(!ST::isComplex) A4 = nonconstTranspose(nonconstAdjoint(A4)); // Should result in CONJ
  if(out.get()) *out << "\nA4 =\n" << describe(*A4,verbLevel);

  if(out.get()) *out << "\nTesting A4 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A4,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nCalling all of the rest of the functions for const just to test them ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A5 = scale(
      Scalar(0.25)
      ,adjoint(
        transpose(
          adjoint(
            scaleAndAdjoint(
              Scalar(4.0)
              ,Thyra::TRANS
              ,origA
              )
            )
          )
        )
      );
  if(!ST::isComplex) A5 = transpose(adjoint(A5)); // Should result in CONJ
  if(out.get()) *out << "\nA5 =\n" << describe(*A5,verbLevel);

  if(out.get()) *out << "\nTesting A5 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A5,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nCreating a multiplied operator A6 = origA^H*A1 ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A6 = multiply(adjoint(origA),A1);
  if(out.get()) *out << "\nA6 =\n" << describe(*A6,verbLevel);

  if(out.get()) *out << "\nTesting A6 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A6,out.get());
  if(!result) success = false;
  // Note that testing the symmetry above helps to check the transpose mode
  // against the non-transpose mode!

#ifdef TEUCHOS_DEBUG
  if(out.get()) *out << "\nCreating an invalid multiplied operator A6b = origA*origA (should throw an exception) ...\n\n";
  try {
    RCP<const Thyra::LinearOpBase<Scalar> >
      A6b = multiply(origA,origA);
    result = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,result)
  if(out.get())
    *out << "\nCaught expected exception : " << (result?"failed\n":"passed\n");
  if(result) success = false;
#endif // TEUCHOS_DEBUG

  if(out.get()) *out << "\nCreating a non-const multiplied operator A7 = origA^H*A1 ...\n";
  RCP<Thyra::LinearOpBase<Scalar> >
    A7 = nonconstMultiply(
      rcp_const_cast<Thyra::LinearOpBase<Scalar> >(adjoint(origA))
      ,rcp_const_cast<Thyra::LinearOpBase<Scalar> >(A1)
      );
  if(out.get()) *out << "\nA7 =\n" << describe(*A7,verbLevel);

  if(out.get()) *out << "\nTesting A7 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A7,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nCreating an added operator A8 = origA + A1 ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A8 = add(origA,A1);
  if(out.get()) *out << "\nA8 =\n" << describe(*A8,verbLevel);

  if(out.get()) *out << "\nTesting A8 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A8,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nCreating a symmetric subtracted operator A8b = A6 + adjoint(origA)*origA ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A8b = subtract(A6,multiply(adjoint(origA),origA));
  if(out.get()) *out << "\nA8b =\n" << describe(*A8b,verbLevel);

  if(out.get()) *out << "\nTesting A8b ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A8b,out.get());
  if(!result) success = false;

#ifdef TEUCHOS_DEBUG
  if(out.get()) *out << "\nCreating an invalid added operator A8c = origA + adjoint(origA) (should throw an exception) ...\n\n";
  try {
    RCP<const Thyra::LinearOpBase<Scalar> >
      A8c = add(origA,adjoint(origA));
    result = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,result)
  if(out.get())
    *out << "\nCaught expected exception : " << (result?"failed\n":"passed\n");
  if(result) success = false;
#endif // TEUCHOS_DEBUG

  RCP<const Thyra::LinearOpBase<Scalar> >
    nullOp = null;

  if(out.get()) *out << "\nCreating a blocked 2x2 linear operator A9 = [ A6, A1^H; A1, null ] ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A9 = Thyra::block2x2<Scalar>(
      A6,  adjoint(A1)
      ,A1, nullOp
      );
  if(out.get()) *out << "\nA9 =\n" << describe(*A9,verbLevel);
  
  if(out.get()) *out << "\nTesting A9 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A9,out.get());
  if(!result) success = false;
  // Note that testing the symmetry above helps to check the transpose mode
  // against the non-transpose mode!

  if(out.get()) *out << "\nCreating a blocked 2x2 linear operator A9_a = [ A6, A1^H; A1, null ] using pre-formed range and domain product spaces ...\n";
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
  if(out.get()) *out << "\nA9_a =\n" << describe(*A9_a,verbLevel);
  
  if(out.get()) *out << "\nTesting A9_a ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A9_a,out.get());
  if(!result) success = false;
  // Note that testing the symmetry above helps to check the transpose mode
  // against the non-transpose mode!
  
  if(out.get()) *out << "\nComparing A9 == A9_a ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.compare(*A9,*A9_a,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nCreating a blocked 2x2 linear operator A9_b = [ A6, A1^H; A1, null ] using flexible fill ...\n";
  RCP<Thyra::PhysicallyBlockedLinearOpBase<Scalar> >
    A9_b = rcp(new Thyra::DefaultBlockedLinearOp<Scalar>());
  A9_b->beginBlockFill();
  A9_b->setBlock(0,0,A6);
  A9_b->setBlock(0,1,adjoint(A1));
  A9_b->setBlock(1,0,A1);
  A9_b->endBlockFill();
  if(out.get()) *out << "\nA9_b =\n" << describe(*A9_b,verbLevel);
  
  if(out.get()) *out << "\nTesting A9_b ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A9_b,out.get());
  if(!result) success = false;
  // Note that testing the symmetry above helps to check the transpose mode
  // against the non-transpose mode!
  
  if(out.get()) *out << "\nComparing A9 == A9_b ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.compare(*A9,*A9_b,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nCreating a blocked 2x2 linear operator A9a = [ null, A1^H; A1, null ] ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A9a = Thyra::block2x2<Scalar>(
      nullOp,  adjoint(A1),
      A1,      nullOp
      );
  if(out.get()) *out << "\nA9a =\n" << describe(*A9a,verbLevel);
  
  if(out.get()) *out << "\nTesting A9a ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A9a,out.get());
  if(!result) success = false;
  // Note that testing the symmetry above helps to check the transpose mode
  // against the non-transpose mode!
  
#ifdef TEUCHOS_DEBUG
  if(out.get()) *out << "\nCreating an invalid blocked 2x2 operator A9b = [ A6, A1^H; A1, A1 ] (should throw an exception) ...\n\n";
  try {
    RCP<const Thyra::LinearOpBase<Scalar> >
      A9b = Thyra::block2x2<Scalar>(
        A6,  adjoint(A1),
        A1,  A1
        );
    result = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,result)
  if(out.get())
    *out << "\nCaught expected exception : " << (result?"failed\n":"passed\n");
  if(result) success = false;
#endif // TEUCHOS_DEBUG

#ifdef TEUCHOS_DEBUG
  if(out.get()) *out << "\nCreating an invalid blocked 2x2 operator A9c = [ A1, A1 ; null, null ] (should throw an exception) ...\n\n";
  try {
    RCP<const Thyra::LinearOpBase<Scalar> >
      A9c = Thyra::block2x2<Scalar>(
        A1,       A1,
        nullOp,   nullOp
        );
    result = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,result)
  if(out.get())
    *out << "\nCaught expected exception : " << (result?"failed\n":"passed\n");
  if(result) success = false;
#endif // TEUCHOS_DEBUG

#ifdef TEUCHOS_DEBUG
  if(out.get()) *out << "\nCreating an invalid blocked 2x2 operator A9d = [ A1, null; A1, null ] (should throw an exception) ...\n\n";
  try {
    RCP<const Thyra::LinearOpBase<Scalar> >
      A9d = Thyra::block2x2<Scalar>(
        A1,  nullOp,
        A1,  nullOp
        );
    result = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,result)
  if(out.get())
    *out << "\nCaught expected exception : " << (result?"failed\n":"passed\n");
  if(result) success = false;
#endif // TEUCHOS_DEBUG

  if(out.get()) *out << "\nCreating a blocked 2x1 linear operator A10 = [ A6; A1 ] ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A10 = Thyra::block2x1<Scalar>(
      A6,
      A1
      );
  if(out.get()) *out << "\nA10 =\n" << describe(*A10,verbLevel);
  
  if(out.get()) *out << "\nTesting A10 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A10,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nCreating a blocked 1x2 linear operator A11 = [ A9, A10 ] ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A11 = Thyra::block1x2<Scalar>( A9, A10 );
  if(out.get()) *out << "\nA11 =\n" << describe(*A11,verbLevel);
  
  if(out.get()) *out << "\nTesting A11 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A11,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nCreating a zero linear operator A12 = 0 (range and domain spaces of origA) ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A12 = Thyra::zero(origA->range(),origA->domain());
  if(out.get()) *out << "\nA12 =\n" << describe(*A12,verbLevel);

  if(out.get()) *out << "\nTesting A12 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.check(*A12,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nCreating a blocked 2x2 linear operator A13 = [ zero, A1^H; A1, zero ] ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A13 = Thyra::block2x2<Scalar>(
      Thyra::zero(A1->domain(),A1->domain()),   adjoint(A1),
      A1,                                       Thyra::zero(A1->range(),A1->range())
      );
  if(out.get()) *out << "\nA13 =\n" << describe(*A13,verbLevel);
  
  if(out.get()) *out << "\nComparing A9a == A13 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = linearOpTester.compare(*A9a,*A13,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nCreating a zero linear operator A14 = I (range space of origA) ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> >
    A14 = Thyra::identity(origA->range());
  if(out.get()) *out << "\nA14 =\n" << describe(*A14,verbLevel);
  
  if(out.get()) *out << "\nTesting A14 ...\n";
  Thyra::seed_randomize<Scalar>(0);
  result = symLinearOpTester.check(*A14,out.get());
  if(!result) success = false;
  
  if(out.get()) *out << "\n*** Leaving run_composite_linear_ops_tests<"<<ST::name()<<">(...) ...\n";

  return success;

} // end run_composite_linear_ops_tests() [Doxygen looks for this!]

int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  const Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> >
    comm = Teuchos::DefaultComm<Thyra::Ordinal>::getComm();

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  try {

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

#ifdef HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_composite_linear_ops_tests<float>(comm,n,useSpmd,float(1e-4),dumpAll,verbose?&*out:NULL) ) success = false;
#endif // HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_composite_linear_ops_tests<double>(comm,n,useSpmd,double(1e-12),dumpAll,verbose?&*out:NULL) ) success = false;
#if defined(HAVE_THYRA_COMPLEX)
#ifdef HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_composite_linear_ops_tests<std::complex<float> >(comm,n,useSpmd,float(1e-4),dumpAll,verbose?&*out:NULL) ) success = false;
#endif // HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_composite_linear_ops_tests<std::complex<double> >(comm,n,useSpmd,double(1e-12),dumpAll,verbose?&*out:NULL) ) success = false;
#endif
#if defined(HAVE_TEUCHOS_GNU_MP) && !defined(USE_MPI) // mpf_class can not be used with MPI yet!
    if( !run_composite_linear_ops_tests<mpf_class>(comm,n,useSpmd,mpf_class(1e-12),dumpAll,verbose?&*out:NULL) ) success = false;
#endif

  } // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)

  if( verbose ) {
    if(success) *out << "\nAll of the tests seem to have run successfully!\n";
    else        *out << "\nOh no! at least one of the tests failed!\n";	
  }
  
  return success ? 0 : 1;

} // end main() [Doxygen looks for this!]
