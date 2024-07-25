// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Bring in all of the operator/vector ANA client support software
#include "Thyra_OperatorVectorClientSupport.hpp"

// Just use a default SPMD space for this example
#include "Thyra_DefaultSpmdVectorSpace.hpp"

// Some utilities from Teuchos
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_GlobalMPISession.hpp"


/** This is an example of how to use the implicitly composed linear operator
 *  subclasses to build a complex implicit linear operator.
 *
 * The linear operator that we are trying to build is:
 
 \verbatim

   M = [ M00, M01 ]
       [ M10, M11 ]

 \endverbatim

 * where:

 \verbatim

   M00 = [ gama*B*A + C,  E + F ] ^H
         [ J^H * A,       I     ]
      
       = ( n0, n1 ) x ( n0, n1 )

   M01 = beta * [ Q ]
                [ K ]

       = ( n0, n1 ) x n2

   M10 = [ L * N^H,  eta*P ]

       = n2 x ( n0, n1 )

   M11 = D - Q^H*Q

 \endverbatim

 * Above, <tt>I1</tt> and <tt>I2</tt> are the identity operators.
 *
 * This system of linear operators is build around three basic vector spaces
 * <tt>space0</tt>, <tt>space1</tt> and <tt>space2</tt> with the dimensions
 * <tt>n0</tt>, <tt>n1</tt>, and <tt>n2</tt> respectively.  For this example,
 * we will just create random Multi-vectors for each of the component linear
 * operators.
 */
template<class Scalar>
int exampleImplicitlyComposedLinearOperators(
  const int n0,
  const int n1,
  const int n2,
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel,
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType errorTol,
  const bool testAdjoint
  )
{

  // Using and other declarations
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::OSTab;
  using Thyra::VectorSpaceBase;
  using Thyra::VectorBase;
  using Thyra::MultiVectorBase;
  using Thyra::LinearOpBase;
  using Thyra::defaultSpmdVectorSpace;
  using Thyra::randomize;
  using Thyra::identity;
  using Thyra::diagonal;
  using Thyra::multiply;
  using Thyra::add;
  using Thyra::subtract;
  using Thyra::scale;
  using Thyra::adjoint;
  using Thyra::block1x2;
  using Thyra::block2x2;
  using Thyra::block2x2;

  out << "\n***"
      << "\n*** Demonstrating building linear operators for scalar type "
      << ST::name()
      << "\n***\n";

  OSTab tab(out);

  //
  // A) Set up the basic objects and other inputs to build the implicitly
  // composed linear operators.
  //
  
  // Create serial vector spaces in this case
  const RCP<const VectorSpaceBase<Scalar> >
    space0 = defaultSpmdVectorSpace<Scalar>(n0),
    space1 = defaultSpmdVectorSpace<Scalar>(n1),
    space2 = defaultSpmdVectorSpace<Scalar>(n2);

  // Create the component linear operators first as multi-vectors
  const RCP<MultiVectorBase<Scalar> >
    mvA = createMembers(space2, n0, "A"),
    mvB = createMembers(space0, n2, "B"),
    mvC = createMembers(space0, n0, "C"),
    mvE = createMembers(space0, n1, "E"),
    mvF = createMembers(space0, n1, "F"),
    mvJ = createMembers(space2, n1, "J"),
    mvK = createMembers(space1, n2, "K"),
    mvL = createMembers(space2, n1, "L"),
    mvN = createMembers(space0, n1, "N"),
    mvP = createMembers(space2, n1, "P"),
    mvQ = createMembers(space0, n2, "Q");

  // Create the vector diagonal for D
  const RCP<VectorBase<Scalar> > d = createMember(space2);

  // Get the constants
  const Scalar
    one = 1.0,
    beta = 2.0,
    gamma = 3.0,
    eta = 4.0;

  // Randomize the values in the Multi-Vector
  randomize( -one, +one, mvA.ptr() );
  randomize( -one, +one, mvB.ptr() );
  randomize( -one, +one, mvC.ptr() );
  randomize( -one, +one, d.ptr() );
  randomize( -one, +one, mvE.ptr() );
  randomize( -one, +one, mvF.ptr() );
  randomize( -one, +one, mvJ.ptr() );
  randomize( -one, +one, mvK.ptr() );
  randomize( -one, +one, mvL.ptr() );
  randomize( -one, +one, mvN.ptr() );
  randomize( -one, +one, mvP.ptr() );
  randomize( -one, +one, mvQ.ptr() );

  // Get the linear operator forms of the basic component linear operators
  const RCP<const LinearOpBase<Scalar> >
    A = mvA,
    B = mvB,
    C = mvC,
    E = mvE,
    F = mvF,
    J = mvJ,
    K = mvK,
    L = mvL,
    N = mvN,
    P = mvP,
    Q = mvQ;

  out << describe(*A, verbLevel);
  out << describe(*B, verbLevel);
  out << describe(*C, verbLevel);
  out << describe(*E, verbLevel);
  out << describe(*F, verbLevel);
  out << describe(*J, verbLevel);
  out << describe(*K, verbLevel);
  out << describe(*L, verbLevel);
  out << describe(*N, verbLevel);
  out << describe(*P, verbLevel);
  out << describe(*Q, verbLevel);

  //
  // B) Create the composed linear operators
  //

  // I
  const RCP<const LinearOpBase<Scalar> > I = identity(space1, "I");

  // D = diag(d)
  const RCP<const LinearOpBase<Scalar> > D = diagonal(d, "D");

  // M00 = [ gama*B*A + C,  E + F ] ^H
  //       [ J^H * A,       I     ]
  const RCP<const LinearOpBase<Scalar> > M00 =
    adjoint(
      block2x2(
        add( scale(gamma,multiply(B,A)), C ),  add( E, F ),
        multiply(adjoint(J),A),                I
        ),
      "M00"
      );

  out << "\nM00 = " << describe(*M00, verbLevel);

  // M01 = beta * [ Q ]
  //              [ K ]
  const RCP<const LinearOpBase<Scalar> > M01 =
    scale(
      beta,
      block2x1( Q, K ),
      "M01"
      );

  out << "\nM01 = "  << describe(*M01, verbLevel);
            
  // M10 = [ L * N^H,  eta*P ]
  const RCP<const LinearOpBase<Scalar> > M10 =
    block1x2(
      multiply(L,adjoint(N)),  scale(eta,P),
      "M10"
      );

  out << "\nM10 = " << describe(*M10, verbLevel);

  // M11 = D - Q^H*Q
  const RCP<const LinearOpBase<Scalar> > M11 =
    subtract( D, multiply(adjoint(Q),Q), "M11" );

  out << "\nM11 = "  << describe(*M11, verbLevel);
  

  // M = [ M00, M01 ]
  //     [ M10, M11 ]
  const RCP<const LinearOpBase<Scalar> > M =
    block2x2(
      M00, M01,
      M10, M11,
      "M"
      );

  out << "\nM = " << describe(*M, verbLevel);

  
  //
  // C) Test the final composed operator
  //

  Thyra::LinearOpTester<Scalar> linearOpTester;
  linearOpTester.set_all_error_tol(errorTol);
  linearOpTester.check_adjoint(testAdjoint);
  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH))
    linearOpTester.show_all_tests(true);
  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME))
    linearOpTester.dump_all(true);

  const bool result = linearOpTester.check(*M, Teuchos::inOutArg(out));

  return result;

}


//
// Main driver program
//
int main( int argc, char *argv[] )
{
  
  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool result;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  // Above is needed to run in an MPI build with some MPI implementations

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read in command-line options
    //

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    int n0 = 2;
    clp.setOption( "n0", &n0 );

    int n1 = 3;
    clp.setOption( "n1", &n1 );

    int n2 = 4;
    clp.setOption( "n2", &n2 );

    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_MEDIUM;
    setVerbosityLevelOption( "verb-level", &verbLevel,
      "Top-level verbosity level.  By default, this gets deincremented as you go deeper into numerical objects.",
      &clp );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

#if defined(HAVE_TEUCHOS_INST_FLOAT)
    // Run using float
    result = exampleImplicitlyComposedLinearOperators<float>(
      n0, n1, n2, *out, verbLevel, 1e-5, true );
    if (!result) success = false;
#endif

    // Run using double
    result = exampleImplicitlyComposedLinearOperators<double>(
      n0, n1, n2, *out, verbLevel, 1e-12, true );
    if (!result) success = false;

#if defined(HAVE_TEUCHOS_INST_COMPLEX_FLOAT) && defined(HAVE_TEUCHOS_INST_FLOAT)
    // Run using std::complex<float>
    result = exampleImplicitlyComposedLinearOperators<std::complex<float> >(
      n0, n1, n2, *out, verbLevel, 1e-5, false );
    if (!result) success = false;
    // 2007/11/02: rabartl: Above, I skip the test of the adjoint since it
    // fails but a lot.  On my machine, the relative error is:
    //    rel_err((-3.00939,-0.836347),(-0.275689,1.45244)) = 1.14148.
    // Since this works just fine for the next complex<double> case, I am
    // going to just skip this test.
#endif

#if defined(HAVE_TEUCHOS_INST_COMPLEX_DOUBLE)
    // Run using std::complex<double>
    result = exampleImplicitlyComposedLinearOperators<std::complex<double> >(
      n0, n1, n2, *out, verbLevel, 1e-12, true );
    if (!result) success = false;
#endif

#ifdef HAVE_TEUCHOS_GNU_MP

    // Run using mpf_class
    result = exampleImplicitlyComposedLinearOperators<mpf_class>(
      n0, n1, n2, *out, verbLevel, 1e-20, true );
    if (!result) success = false;

#endif // HAVE_TEUCHOS_GNU_MP

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)
    
  if(success)
    *out << "\nCongratulations! All of the tests checked out!\n";
  else
    *out << "\nOh no! At least one of the tests failed!\n";
  
  return success ? 0 : 1;

} // end main()
