// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA.H"
#include "LOCA_LAPACK.H"

#include "LOCA_Eigensolver_AbstractStrategy.H"
#include "LOCA_Parameter_SublistParser.H"

#include "ChanProblemInterface.H"
#include "NOX_TestCompare.H"

#include "Teuchos_GlobalMPISession.hpp"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpi_session(&argc, &argv);

  int n = 100;
  double alpha = 0.0;
  double beta = 0.0;
  double scale = 1.0;
  int ierr = 0;
  int nev = 10;
  int narn = 20;
  double arntol = 1.0e-12;

  alpha = alpha / scale;

  try {

    bool verbose = false;
    // Check for verbose output
    if (argc>1)
      if (argv[1][0]=='-' && argv[1][1]=='v') 
	verbose = true;

    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");

    // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
    Teuchos::ParameterList& aList = stepperList.sublist("Eigensolver");
    aList.set("Method", "Anasazi");
    aList.set("Operator", "Jacobian Inverse");
    aList.set("Block Size", 1);
    aList.set("Num Blocks", narn);
    aList.set("Num Eigenvalues", nev);
    aList.set("Convergence Tolerance", arntol);
    aList.set("Step Size", 1);
    aList.set("Maximum Restarts",2);
    aList.set("Sorting Order","LM");
    if (verbose)
      aList.set("Debug Level",
		Anasazi::Errors + 
		Anasazi::Warnings +
		Anasazi::FinalSummary);
    else
      aList.set("Debug Level", Anasazi::Errors);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");

    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    if (verbose)
       nlPrintParams.set("Output Information", 
				  NOX::Utils::Error +
				  NOX::Utils::Details +
				  NOX::Utils::OuterIteration + 
				  NOX::Utils::InnerIteration + 
				  NOX::Utils::Warning +
				  NOX::Utils::TestDetails + 
				  NOX::Utils::StepperIteration +
				  NOX::Utils::StepperDetails);
     else
       nlPrintParams.set("Output Information", NOX::Utils::Error);

    // Create LAPACK factory
    Teuchos::RCP<LOCA::Abstract::Factory> lapackFactory =
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, lapackFactory);

    // Create parsed parameter list
    Teuchos::RCP<LOCA::Parameter::SublistParser> parsedParams = 
      Teuchos::rcp(new LOCA::Parameter::SublistParser(globalData));
    parsedParams->parseSublists(paramList);

    // Set up the problem interface
    ChanProblemInterface chan(globalData, n, alpha, beta, scale);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("scale",scale);
  
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    LOCA::LAPACK::Group grp(globalData, chan);
    
    grp.setParams(p);

    grp.computeF();
    grp.computeJacobian();

    // Create Anasazi eigensolver
    Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy> anasaziStrategy
      = globalData->locaFactory->createEigensolverStrategy(
				     parsedParams, 
				     parsedParams->getSublist("Eigensolver"));

    Teuchos::RCP< std::vector<double> > anasazi_evals_r;
    Teuchos::RCP< std::vector<double> > anasazi_evals_i;
    Teuchos::RCP< NOX::Abstract::MultiVector > anasazi_evecs_r;
    Teuchos::RCP< NOX::Abstract::MultiVector > anasazi_evecs_i;
    NOX::Abstract::Group::ReturnType anasaziStatus = 
      anasaziStrategy->computeEigenvalues(grp, 
					  anasazi_evals_r, 
					  anasazi_evals_i,
					  anasazi_evecs_r,
					  anasazi_evecs_i);

    if (anasaziStatus != NOX::Abstract::Group::Ok)
      ++ierr;

    // Change strategy to DGGEV
    aList.set("Method", "DGGEV");
    aList.set("Sorting Order","SM");

    // Create DGGEV eigensolver
    Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy> dggevStrategy
      = globalData->locaFactory->createEigensolverStrategy(
				      parsedParams,
				      parsedParams->getSublist("Eigensolver"));

    Teuchos::RCP< std::vector<double> > dggev_evals_r;
    Teuchos::RCP< std::vector<double> > dggev_evals_i;
    Teuchos::RCP< NOX::Abstract::MultiVector > dggev_evecs_r;
    Teuchos::RCP< NOX::Abstract::MultiVector > dggev_evecs_i;
    NOX::Abstract::Group::ReturnType dggevStatus = 
      dggevStrategy->computeEigenvalues(grp, 
					dggev_evals_r, 
					dggev_evals_i,
					dggev_evecs_r,
					dggev_evecs_i);

    if (dggevStatus != NOX::Abstract::Group::Ok)
      ++ierr;

    // Check some statistics on the solution
    NOX::TestCompare testCompare(globalData->locaUtils->out(), 
				 *(globalData->locaUtils));

    if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
      globalData->locaUtils->out() 
	<< std::endl 
	<< "***** Checking solution statistics *****" 
	<< std::endl;

    // Check eigenvalues
    for (int i=0; i<nev; i++) {
      stringstream sstr;
      sstr << "Eigenvalue " << i;
      ierr += testCompare.testValue((*anasazi_evals_r)[i], 
				    (*dggev_evals_r)[i], arntol*1e3,
				    sstr.str(),
				    NOX::TestCompare::Relative);
    }

    LOCA::destroyGlobalData(globalData);

  }

 catch (std::exception& e) {
    cout << e.what() << endl;
    ierr = 1;
  } 
  catch (const char *s) {
    cout << s << endl;
    ierr = 1;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
    ierr = 1;
  }

   if (ierr == 0)
     cout << "All tests passed!" << endl;
   else
     cout << ierr << " test(s) failed!" << endl;

  return ierr;
}
