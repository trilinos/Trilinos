// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA.H"
#include "LOCA_LAPACK.H"

#include "LOCA_GlobalData.H"
#include "LOCA_Utils.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Factory.H"
#include "LOCA_Eigensolver_AbstractStrategy.H"
#include "LOCA_Parameter_SublistParser.H"

#include "ChanProblemInterface.H"
#include "NOX_TestCompare.H"

int main(int argc, char *argv[])
{
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

    // Set up the problem interface
    ChanProblemInterface chan(n, alpha, beta, scale);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("scale",scale);
  
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    LOCA::LAPACK::Group grp(chan);
    
    grp.setParams(p);

    grp.computeF();
    grp.computeJacobian();

    // Create parameter list
    Teuchos::RefCountPtr<NOX::Parameter::List> paramList = 
      Teuchos::rcp(new NOX::Parameter::List);

    // Create LOCA sublist
    NOX::Parameter::List& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    NOX::Parameter::List& stepperList = locaParamsList.sublist("Stepper");

    // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
    NOX::Parameter::List& aList = stepperList.sublist("Eigensolver");
    aList.setParameter("Method", "Anasazi");
    aList.setParameter("Operator", "Jacobian Inverse");
    aList.setParameter("Block Size", 1);
    aList.setParameter("Arnoldi Size", narn);
    aList.setParameter("NEV", nev);
    aList.setParameter("Tol", arntol);
    aList.setParameter("Convergence Check", 1);
    aList.setParameter("Restarts",2);
    aList.setParameter("Sorting Order","LM");
    aList.setParameter("Debug Level",0);

    // Set the LOCA Utilities
    NOX::Parameter::List& locaUtilsList = locaParamsList.sublist("Utilities");
    if (verbose) {
      locaUtilsList.setParameter("Output Information", 
				 LOCA::Utils::Error + 
				 LOCA::Utils::Warning +
				 LOCA::Utils::StepperIteration +
				 LOCA::Utils::StepperDetails +
				 LOCA::Utils::Solver +
				 LOCA::Utils::Parameters +
				 LOCA::Utils::SolverDetails);
    }
    else
      locaUtilsList.setParameter("Output Information", LOCA::Utils::Error);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    NOX::Parameter::List& nlParams = paramList->sublist("NOX");

    NOX::Parameter::List& nlPrintParams = nlParams.sublist("Printing");
    if (verbose)
       nlPrintParams.setParameter("Output Information", 
				  NOX::Utils::Error +
				  NOX::Utils::Details +
				  NOX::Utils::OuterIteration + 
				  NOX::Utils::InnerIteration + 
				  NOX::Utils::Warning +
				  NOX::Utils::TestDetails);
     else
       nlPrintParams.setParameter("Output Information", NOX::Utils::Error);

    // Create printing utils
    Teuchos::RefCountPtr<LOCA::Utils> locaUtils =
      Teuchos::rcp(new LOCA::Utils);
    locaUtils->setUtils(*paramList);

    // Create error check
    Teuchos::RefCountPtr<LOCA::ErrorCheck> locaErrorCheck = 
      Teuchos::rcp(new LOCA::ErrorCheck);

    Teuchos::RefCountPtr<LOCA::Factory> locaFactory;

    // Create global data object
    Teuchos::RefCountPtr<LOCA::GlobalData> locaGlobalData =
      Teuchos::rcp(new LOCA::GlobalData(locaUtils, 
					locaErrorCheck, 
					locaFactory));

    // Create parsed parameter list
    Teuchos::RefCountPtr<LOCA::Parameter::SublistParser> parsedParams = 
      Teuchos::rcp(new LOCA::Parameter::SublistParser(locaGlobalData));
    parsedParams->parseSublists(paramList);

    // Create LAPACK factory
    Teuchos::RefCountPtr<LOCA::Abstract::Factory> lapackFactory =
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create factory
    locaFactory = Teuchos::rcp(new LOCA::Factory(locaGlobalData, 
						 lapackFactory));

    // Create Anasazi eigensolver
    Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy> anasaziStrategy
      = locaFactory->createEigensolverStrategy(
				     parsedParams, 
				     parsedParams->getSublist("Eigensolver"));

    Teuchos::RefCountPtr< std::vector<double> > anasazi_evals_r;
    Teuchos::RefCountPtr< std::vector<double> > anasazi_evals_i;
    Teuchos::RefCountPtr< NOX::Abstract::MultiVector > anasazi_evecs_r;
    Teuchos::RefCountPtr< NOX::Abstract::MultiVector > anasazi_evecs_i;
    NOX::Abstract::Group::ReturnType anasaziStatus = 
      anasaziStrategy->computeEigenvalues(grp, 
					  anasazi_evals_r, 
					  anasazi_evals_i,
					  anasazi_evecs_r,
					  anasazi_evecs_i);

    if (anasaziStatus != NOX::Abstract::Group::Ok)
      ++ierr;

    // Change strategy to DGGEV
    aList.setParameter("Method", "DGGEV");
    aList.setParameter("Sorting Order","SM");

    // Create DGGEV eigensolver
    Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy> dggevStrategy
      = locaFactory->createEigensolverStrategy(
				      parsedParams,
				      parsedParams->getSublist("Eigensolver"));

    Teuchos::RefCountPtr< std::vector<double> > dggev_evals_r;
    Teuchos::RefCountPtr< std::vector<double> > dggev_evals_i;
    Teuchos::RefCountPtr< NOX::Abstract::MultiVector > dggev_evecs_r;
    Teuchos::RefCountPtr< NOX::Abstract::MultiVector > dggev_evecs_i;
    NOX::Abstract::Group::ReturnType dggevStatus = 
      dggevStrategy->computeEigenvalues(grp, 
					dggev_evals_r, 
					dggev_evals_i,
					dggev_evecs_r,
					dggev_evecs_i);

    if (dggevStatus != NOX::Abstract::Group::Ok)
      ++ierr;

    // Check some statistics on the solution
    NOX::Utils utils(nlPrintParams);
    NOX::TestCompare testCompare(cout, utils);

    if (utils.isPrintProcessAndType(NOX::Utils::TestDetails))
      cout << endl << "***** Checking solution statistics *****" << endl;

    // Check eigenvalues
    for (int i=0; i<nev; i++) {
      stringstream sstr;
      sstr << "Eigenvalue " << i;
      ierr += testCompare.testValue((*anasazi_evals_r)[i], 
				    (*dggev_evals_r)[i], arntol*1e3,
				    sstr.str(),
				    NOX::TestCompare::Relative);
    }

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
