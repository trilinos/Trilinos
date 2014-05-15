//
//
// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Stratimikos_Config.h"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "simpleStratimikosSolve.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraLinearOp.hpp"
//  #include "Thyra_LinearOpWithSolveFactoryExamples.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Teuchos_ParameterList.hpp"

namespace Teuchos { class ParameterList; }

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif

bool TestSingleStratimikosSolver(
  Teuchos::ParameterList                  *paramList_inout
  ,const bool                             dumpAll
  ,Teuchos::FancyOStream                  *out
  )
{

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::OSTab;
  using Teuchos::ParameterList;
  using Teuchos::getParameter;
  bool result, success = true;

  try {

    TEUCHOS_TEST_FOR_EXCEPT(!paramList_inout);

    RCP<ParameterList>
      paramList = rcp(paramList_inout,false);

#if 0
    if(out) {
      *out << "\nEchoing input parameters ...\n";
      paramList->print(*out,1,true,false);
    }
#endif

    // Create list of valid parameter sublists
    Teuchos::ParameterList validParamList("TestSingleStratimikosSolver");
    validParamList.set("Matrix File","fileName");
    validParamList.sublist("Linear Solver Builder");
    
    Teuchos::ParameterList &StratimikosList = paramList_inout->sublist("Linear Solver Builder");

#if 0 

    this fails because validParamList has only been set one level deep - see 5-10 lines above

    if(out) *out << "\nValidating top-level input parameters ...\n";
    StratimikosList.validateParameters(validParamList.sublist("Linear Solver Builder"),0);
#endif

#if 0
    paramList->validateParameters(validParamList,0);

    if(out) *out << "\nValidating top-level input parameters ...\n";
    paramList->validateParameters(StratimikosList,0);
#endif

    const std::string
      &matrixFile = getParameter<std::string>(*paramList,"Matrix File");

    if(out) *out << "\nReading in an epetra matrix A from the file \'"<<matrixFile<<"\' ...\n";
  
#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
    Teuchos::RCP<Epetra_CrsMatrix> epetra_A;
    EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A );

    const int num_vecs = 1 ; 
    const Epetra_Map& A_domainmap = epetra_A->DomainMap() ; 
    Epetra_MultiVector X( A_domainmap, num_vecs ) ; 
    const Epetra_Map& A_rangemap = epetra_A->RangeMap() ; 
    Epetra_MultiVector B( A_rangemap, num_vecs ) ; 
    B.Random(); 

    int status = simpleStratimikosSolve( *epetra_A, B, &X, &StratimikosList );
    assert( status==0 ) ; 


    int nrows =epetra_A->NumGlobalRows();
    Epetra_MultiVector Residual( B ) ;
    epetra_A->Apply( X, Residual );
    Residual.Update( 1.0, B, -1.0 );
    double ResidualNorm;
    double BNorm;
    Residual.Norm2( &ResidualNorm );
    B.Norm2( &BNorm );

    const Teuchos::ParameterList& LST_PL =  StratimikosList.sublist("Linear Solver Types"); 
    const Teuchos::ParameterList& Amesos_PL =  LST_PL.sublist("Amesos"); 
    std::string SolverType = Amesos_PL.get<std::string>("Solver Type");
    assert( SolverType == "Lapack" ) ; // for now - remove later 
    const Teuchos::ParameterList& AmesosSettings_PL =  Amesos_PL.sublist("Amesos Settings"); 
    const Teuchos::ParameterList& SolverType_PL =  AmesosSettings_PL.sublist(SolverType); 
    double AddToDiag = -13.0;
    if (  SolverType == "Lapack" ) {
      AddToDiag = SolverType_PL.get<double>("AddToDiag");
    }
    assert( AddToDiag >= 0.0 ) ; 
    double MaxError = 1e-15*BNorm + 10.0 * AddToDiag ; 
    double MinError = 0.02 * ( 1e-15*BNorm + AddToDiag ); 
    success = ResidualNorm < nrows * MaxError && 
      ResidualNorm > MinError ; 

#if 0
    std::cout << __FILE__ << "::" << __LINE__ 
	 << " ResidualNorm = " << ResidualNorm  
	 << " MinError = " << MinError 
	 << " MaxError = " << MaxError << std::endl ; 
    std::cout << " B = " ; 
    B.Print( std::cout ) ; 
    std::cout << " X = " ; 
    X.Print( std::cout ) ; 
    std::cout << " epetra_A = " ; 
    epetra_A->Print( std::cout ) ; 
    std::cout << " success = " << success << " ResidualNorm = " <<  ResidualNorm << std::endl ; 
#endif

    if(false && out) {
      *out << "\nPrinting the parameter list (showing what was used) ...\n";
      paramList->print(*out,1,true,true);
    }
    
  }
  catch( const std::exception &excpt ) {
    std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
    success = false;
  }
  
  return success;
  
}

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char* argv[])
{

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  
  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool verbose = true;

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read options from command-line
    //
    
    std::string     inputFile              = "stratimikos_amesos_lapack.xml";
    std::string     extraParams            = "";
    bool            dumpAll                = false;

    CommandLineProcessor  clp(false); // Don't throw exceptions
    std::string default_inputFile ;
    sprintf( &default_inputFile[0], "Input file [%s]", &inputFile[0] );
    clp.setOption( "input-file", &inputFile, &default_inputFile[0], false );
    clp.setOption( "extra-params", &extraParams, "Extra parameters overriding the parameters read in from --input-file");
    clp.setDocString(
      "Testing program for Trilinos (and non-Trilinos) linear solvers access through Thyra."
      );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    Teuchos::ParameterList paramList;
    if(verbose) *out << "\nReading parameters from XML file \""<<inputFile<<"\" ...\n";
    Teuchos::updateParametersFromXmlFile(inputFile,&paramList);
    if(extraParams.length()) {
      if(verbose) *out << "\nAppending extra parameters from the XML string \""<<extraParams<<"\" ...\n";
      Teuchos::updateParametersFromXmlString(extraParams,&paramList);
    }
    
    success
      = TestSingleStratimikosSolver(
        &paramList,dumpAll,verbose?&*out:0
        );
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success)
  
  if (verbose) {
    if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
    else         *out << "\nOh no! At least one of the tests failed!\n";
  }

  return ( success ? 0 : 1 );
}
