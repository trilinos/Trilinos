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

#include "ChanProblemInterface.H"
#include "NOX_TestCompare.H"

int main(int argc, char *argv[])
{
  int n = 100;
  double alpha = 0.0;
  double beta = 0.0;
  double scale = 1.0;
  int ierr = 0;

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
    NOX::Parameter::List paramList;

    // Create LOCA sublist
    NOX::Parameter::List& locaParamsList = paramList.sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    NOX::Parameter::List& stepperList = locaParamsList.sublist("Stepper");

    // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
    NOX::Parameter::List& aList = stepperList.sublist("Eigensolver");
    aList.setParameter("Method", "Anasazi");
    aList.setParameter("Operator", "Jacobian Inverse");
    aList.setParameter("Block Size", 1);
    aList.setParameter("Arnoldi Size", 10);
    aList.setParameter("NEV", 3);
    aList.setParameter("Tol", 2.0e-7);
    aList.setParameter("Convergence Check", 1);
    aList.setParameter("Restarts",2);
    aList.setParameter("Sorting Order","LM");
    aList.setParameter("Debug Level",1);

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
    NOX::Parameter::List& nlParams = paramList.sublist("NOX");
    nlParams.setParameter("Nonlinear Solver", "Line Search Based");

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
    locaUtils->setUtils(paramList);

    // Create error check
    Teuchos::RefCountPtr<LOCA::ErrorCheck> locaErrorCheck = 
      Teuchos::rcp(new LOCA::ErrorCheck);

    Teuchos::RefCountPtr<LOCA::Factory> locaFactory;

    // Create global data object
    Teuchos::RefCountPtr<LOCA::GlobalData> locaGlobalData =
      Teuchos::rcp(new LOCA::GlobalData(locaUtils, 
					locaErrorCheck, 
					locaFactory));

    // Create factory
    locaFactory = Teuchos::rcp(new LOCA::Factory(locaGlobalData, 
						 Teuchos::rcp(&paramList,
							      false)));

    // Creat eigensolver
    Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy> anasaziStrategy
      = locaFactory->createEigensolver();

    NOX::Abstract::Group::ReturnType anasaziStatus = 
      anasaziStrategy->computeEigenvalues(grp);

    if (anasaziStatus != NOX::Abstract::Group::Ok)
      ierr = 1;

    // Check some statistics on the solution
    NOX::Utils utils(nlPrintParams);
    NOX::TestCompare testCompare(cout, utils);

    if (utils.isPrintProcessAndType(NOX::Utils::TestDetails))
      cout << endl << "***** Checking solutions statistics *****" << endl;

    if (ierr == 0)
      cout << "All tests passed!" << endl;
    else
      cout << ierr << " test(s) failed!" << endl;
  }

  catch (string& s) {
    cout << s << endl;
  }
  catch (char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

  return ierr;
}
