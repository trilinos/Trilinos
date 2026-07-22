// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraVector.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"


namespace {

// Helper function to compute a single norm for a vector
template<class Scalar, class LO, class GO, class Node>
double tpetraNorm2( const Tpetra::Vector<Scalar,LO,GO,Node> &v )
{
  Teuchos::Tuple<double,1> norm;
  norm[0] = -1.0;
  v.norm2(norm());
  return norm[0];
}

} // namespace


bool readAndSolveLinearSystem(int argc, char* argv[])
{

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::CommandLineProcessor;
  using Teuchos::EVerbosityLevel;
  typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;

  using Scalar = double;
  using LocalOrdinal = Tpetra::Map<>::local_ordinal_type;
  using GlobalOrdinal = Tpetra::Map<>::global_ordinal_type;
  using NodeType = Tpetra::Map<>::node_type;
  using CrsMatrix_t = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, NodeType>;
  using Vector_t = Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, NodeType>;
  using ST = Teuchos::ScalarTraits<Scalar>;

  bool success = true;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  //
  // A) Program setup code
  //

  //
  // Read options from command-line
  //

  std::string matrixFile = "";
  std::string extraParamsFile = "";
  double tol = 1e-5;
  bool onlyPrintOptions = false;
  bool printXmlFormat = false;
  bool showDoc = true;
  EVerbosityLevel verbLevel = Teuchos::VERB_LOW;

  auto validVerbLevels = Teuchos::getValidVerbLevels();
  auto validVerbLevelsNamesRawStrings = Teuchos::getValidVerbLevelsNamesRawStrings();

  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

  CommandLineProcessor  clp(false); // Don't throw exceptions

  // Set up command-line options for the linear solver that will be used!
  linearSolverBuilder.setupCLP(&clp);

  clp.setOption( "matrix-file", &matrixFile,
    "Defines the matrix and perhaps the RHS and LHS for a linear system to be solved." );
  clp.setOption( "tol", &tol,
    "Tolerance to check against the scaled residual norm." );
  clp.setOption( "extra-params-file", &extraParamsFile,
    "File containing extra parameters in XML format.");
  clp.setOption( "only-print-options", "continue-after-printing-options",
    &onlyPrintOptions,
    "Only print options and stop or continue on" );
  clp.setOption( "print-xml-format", "print-readable-format", &printXmlFormat,
    "Print the valid options in XML format or in readable format." );
  clp.setOption( "show-doc", "hide-doc", &showDoc,
    "Print the valid options with or without documentation." );
  clp.setOption("verb-level", &verbLevel, validVerbLevels.size(),
    validVerbLevels.getRawPtr(), validVerbLevelsNamesRawStrings.getRawPtr(),
    "The verbosity level." );

  clp.setDocString(
    "Simple example for using TrilinosLinearSolver interface.\n"
    "\n"
    "To print out just the valid options use --only-print-options with"
    " --print-xml-format or --print-readable-format.\n"
    "\n"
    "To solve a linear system read from a file use --matrix-file=\"<matrix>.mtx\""
    " with options read from an XML file using"
    " --linear-solver-params-file=\"<paramList>.xml\"\n"
    );

  CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
  if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

  //
  // Print out the valid options and stop if asked
  //

  if(onlyPrintOptions) {
    if(printXmlFormat)
      Teuchos::writeParameterListToXmlOStream(
        *linearSolverBuilder.getValidParameters(), *out );
    else
      linearSolverBuilder.getValidParameters()->print(
        *out,PLPrintOptions().indent(2).showTypes(true).showDoc(showDoc) );
    return 0;
  }

  //
  // B) Set Tpetra-specific code that sets up the linear system to be solved
  //
  // While the below code reads in the Tpetra objects from a file, you can
  // setup the Tpetra objects any way you would like.
  //

  *out << "\nReading in Tpetra matrix A from the file"
       << " \'"<<matrixFile<<"\' ...\n";

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  RCP<const CrsMatrix_t> tpetra_A =
    Tpetra::MatrixMarket::Reader<CrsMatrix_t>::readSparseFile(matrixFile, comm);

  *out << "\nGenerate a random RHS random vector ...\n";
  auto tpetra_b = rcp(new Vector_t(tpetra_A->getRangeMap()));
  tpetra_b->randomize();

  *out << "\nGenerate LHS of zeros ...\n";
  auto tpetra_x = rcp(new Vector_t(tpetra_A->getDomainMap()));
  tpetra_x->putScalar(0.0); // Initial guess is critical!

  *out << "\nPrinting statistics of the Tpetra linear system ...\n";

  *out
    << "\n  Tpetra::CrsMatrix tpetra_A of dimension "
    << tpetra_A->getRangeMap()->getGlobalNumElements()
    << " x " << tpetra_A->getDomainMap()->getGlobalNumElements() << "\n"
    << "  ||tpetra_A||_F = " << tpetra_A->getFrobeniusNorm() << "\n"
    << "  ||tpetra_b||2 = " << tpetraNorm2(*tpetra_b) << "\n"
    << "  ||tpetra_x||2 = " << tpetraNorm2(*tpetra_x) << "\n";

  //
  // C) The "Glue" code that takes Tpetra objects and wraps them as Thyra
  // objects
  //
  // This next set of code wraps the Tpetra objects that define the linear
  // system to be solved as Thyra objects so that they can be passed to the
  // linear solver.
  //

  using TpetraVectorSpace_t =
    const Thyra::TpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, NodeType>;

  RCP<const Thyra::LinearOpBase<double>> A =
    Thyra::createConstLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,NodeType>(tpetra_A);
  RCP<Thyra::VectorBase<Scalar>> x =
    Thyra::tpetraVector(rcp_dynamic_cast<const TpetraVectorSpace_t>(A->domain()),
      tpetra_x);
  RCP<const Thyra::VectorBase<Scalar>> b =
    Thyra::tpetraVector(rcp_dynamic_cast<const TpetraVectorSpace_t>(A->range()),
      tpetra_b);

  //
  // D) Thyra-specific code for solving the linear system
  //

  // Read in the solver parameters from the parameters file and/or from the
  // command line.  This was setup by the command-line options set by the
  // setupCLP(...) function above.
  linearSolverBuilder.readParameters(out.get());

  // Augment parameters if appropriate
  if(extraParamsFile.length()) {
    Teuchos::updateParametersFromXmlFile( "./"+extraParamsFile,
      linearSolverBuilder.getNonconstParameterList().ptr() );
  }

  // Create a linear solver factory given information read from the
  // parameter list.
  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > lowsFactory =
    linearSolverBuilder.createLinearSolveStrategy("");
  *out << "\nlowsFactory: " << describe(*lowsFactory, verbLevel);

  // Setup output stream and the verbosity level
  lowsFactory->setOStream(out);
  lowsFactory->setVerbLevel(verbLevel);

  // Create a linear solver based on the forward operator A
  RCP<Thyra::LinearOpWithSolveBase<Scalar> > A_lows =
    Thyra::linearOpWithSolve(*lowsFactory, A);
  *out << "\nA: " << describe(*A, verbLevel);
  *out << "\nA_lows: " << describe(*A_lows, verbLevel);

  // Solve the linear system (note: the initial guess in 'x' is critical)
  Thyra::SolveStatus<Scalar> status =
    Thyra::solve<Scalar>(*A_lows, Thyra::NOTRANS, *b, x.ptr());
  *out << "\nSolve status:\n" << status;

  // Write the linear solver parameters after they were read
  linearSolverBuilder.writeParamsFile(*lowsFactory);

  //
  // E) Post process the solution and check the error
  //
  // Note that the below code is based only on the Tpetra objects themselves
  // and does not in any way depend or interact with any Thyra-based objects.
  // The point is that most users of Thyra can largely gloss over the fact
  // that Thyra is really being used for anything.
  //

  // Wipe out the Thyra wrapper for x to guarantee that the solution will be
  // written back to tpetra_x!  At the time of this writting this is not
  // really needed but the behavior may change at some point so this is a
  // good idea.
  x = Teuchos::null;

  *out << "\nSolution ||tpetra_x||2 = " << tpetraNorm2(*tpetra_x) << "\n";

  *out << "\nTesting the solution error ||b-A*x||/||b||:\n";

  // r = b - A*x
  auto tpetra_r = rcp(new Vector_t(tpetra_b->getMap()));
  tpetra_r->assign(*tpetra_b);
  tpetra_A->apply(*tpetra_x, *tpetra_r, Teuchos::NO_TRANS, -ST::one(), ST::one());

  const double
    nrm_r = tpetraNorm2(*tpetra_r),
    nrm_b = tpetraNorm2(*tpetra_b),
    rel_err = ( nrm_r / nrm_b );
  const bool
    passed = (rel_err <= tol);

  *out
    << "||b-A*x||/||b|| = " << nrm_r << "/" << nrm_b << " = " << rel_err
    << " < tol = " << tol << " ? " << ( passed ? "passed" : "failed" ) << "\n";

  if(!passed) success = false;

  return success;

}


int main(int argc, char* argv[])
{

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  bool success = false;

  try {
    success = readAndSolveLinearSystem(argc, argv);
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success)

  if(success) *out << "\nCongratulations! All of the tests checked out!\n";
  else        *out << "\nOh no! At least one of the tests failed!\n";

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );

}
