// @HEADER
// ***********************************************************************
//
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "MatrixMarket_Tpetra.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryExamples.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"

int main(int argc, char* argv[])
{

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::OSTab;
  using Teuchos::ParameterList;
  using Teuchos::getParameter;
  typedef double Scalar;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type NO;

  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool verbose = true;


  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);
  {

    //
    // Read options from command-line
    //

    std::string     inputFile              = "";
    std::string     extraParams            = "";
    bool            printUnusedParams      = false;
    RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    int MyPID = rank(*comm);
    verbose = ( out && (MyPID==0) ); 

    CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.addOutputSetupOptions(true);
    clp.setOption( "input-file", &inputFile, "Input file [Required].", true );
    clp.setOption( "extra-params", &extraParams, "Extra parameters overriding the parameters read in from --input-file");
    clp.setOption( "print-unused-params", "no-print-unused-params", &printUnusedParams, "Whether to print all unused parameters at the end.");
    clp.setDocString(
      "Testing program for Trilinos Tpetra linear solvers access through Thyra."
      );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    RCP<ParameterList> paramList = rcp(new Teuchos::ParameterList());
    if(verbose) *out << "\nReading parameters from XML file \""<<inputFile<<"\" ...\n";
    Teuchos::updateParametersFromXmlFile(inputFile, Teuchos::inOutArg(*paramList));
    if(extraParams.length()) {
      if(verbose) *out << "\nAppending extra parameters from the XML string \""<<extraParams<<"\" ...\n";
      Teuchos::updateParametersFromXmlString(extraParams, Teuchos::inOutArg(*paramList));
    }

    if(verbose) {
      *out << "\nEchoing input parameters ...\n";
      paramList->print(*out,1,true,false);
    }

    // Create list of valid parameter sublists
    Teuchos::ParameterList validParamList("test_tpetra_stratimikos_solver");
    validParamList.set("Matrix File","fileName");
    validParamList.sublist("Linear Solver Builder").disableRecursiveValidation();

    if(verbose) *out << "\nValidating top-level input parameters ...\n";
    paramList->validateParametersAndSetDefaults(validParamList);

    const std::string
      &matrixFile = getParameter<std::string>(*paramList,"Matrix File");
    RCP<ParameterList>
      solverBuilderSL  = sublist(paramList,"Linear Solver Builder",true);

    if(verbose) *out << "\nReading in an tpetra matrix A from the file \'"<<matrixFile<<"\' ...\n";

    RCP<const Tpetra::CrsMatrix<Scalar, LO, GO, NO>> tpetra_A = Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<Scalar, LO, GO, NO>>::readSparseFile( matrixFile, comm );

    RCP<Tpetra::Vector<Scalar,LO,GO,NO>> tpetra_x0 = rcp(new Tpetra::Vector<Scalar,LO,GO,NO>(tpetra_A->getDomainMap()));
    RCP<Tpetra::Vector<Scalar,LO,GO,NO>> tpetra_b = rcp(new Tpetra::Vector<Scalar,LO,GO,NO>(tpetra_A->getRangeMap()));
    tpetra_x0->putScalar(1.0);
    tpetra_A->apply(*tpetra_x0, *tpetra_b);
    tpetra_x0->putScalar(0.0);

    RCP<const Thyra::LinearOpBase<double> >
      A = Thyra::createConstLinearOp(Teuchos::rcp_dynamic_cast<const Tpetra::Operator<Scalar,LO,GO,NO>>(tpetra_A));

   RCP<Thyra::VectorBase<double> >
      x0 = Thyra::createVector( tpetra_x0);
   RCP<const Thyra::VectorBase<double> >
      b = Thyra::createVector( tpetra_b); 

    if(verbose) *out << "\nCreating a Stratimikos::DefaultLinearSolverBuilder object ...\n";

    RCP<Thyra::LinearSolverBuilderBase<double> >
      linearSolverBuilder = rcp(new Stratimikos::DefaultLinearSolverBuilder);

    linearSolverBuilder->setParameterList(solverBuilderSL);

    if(verbose) *out << "\nCreating the LinearOpWithSolveFactoryBase object lowsFactory ...\n";
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
      lowsFactory = createLinearSolveStrategy(*linearSolverBuilder);

    if(verbose) *out << "\nChecking the LOWSB interface ...\n";
    RCP<Thyra::LinearOpWithSolveBase<Scalar> >
      lowsA = Thyra::linearOpWithSolve<Scalar>(*lowsFactory, A);

   // Solve the linear system 
   Thyra::SolveStatus<double> status; 
   status = Thyra::solve<double>(*lowsA, Thyra::NOTRANS, *b, x0.ptr());

    success = (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED ? true : false);

    if(verbose && printUnusedParams) {
      *out << "\nPrinting the parameter list (showing what was used) ...\n";
      paramList->print(*out,1,true,true);
    }

  } //End Tpetra scope
  }
  catch( const std::exception &excpt ) {
    std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
    success = false;
  }

  if (verbose) {
    if(success)  *out << "\nCongratulations! The test passed!\n";
    else         *out << "\nOh no! At least one of the solves failed!\n";
  }

   return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
