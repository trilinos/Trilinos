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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "LOCA.H"
#include "LOCA_LAPACK.H"
#include "ThermoSyphonProblemInterface.H"

int main()
{
  int m_max = 512;
  int maxNewtonIters = 5;

  try {

    // Set up the problem interface
    ThermoSyphonProblemInterface therm("therm/in","therm/oldfile",m_max);
    LOCA::ParameterVector p = therm.getParams();

    p.print(cout);

    // Get dimension of problem
    int m = therm.getSize();
  
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    LOCA::LAPACK::Group grp(therm, m, m, 3*m_max+3, 3*m_max+3);
    grp.setParams(p);

    // Create parameter list
    NOX::Parameter::List nlParams;
    nlParams.setParameter("Nonlinear Solver", "Line Search Based");

    // Create the line search parameters sublist
    NOX::Parameter::List& lineSearchParameters = 
      nlParams.sublist("Line Search");

    // Set the line search method
    lineSearchParameters.setParameter("Method","Full Step");

    // Create the newton and  linear solver parameters sublist
    NOX::Parameter::List& directionParameters = 
      nlParams.sublist("Direction");
    NOX::Parameter::List& newtonParameters = 
      directionParameters.sublist("Newton");
    NOX::Parameter::List& linearSolverParameters = 
      newtonParameters.sublist("Linear Solver");

    NOX::Parameter::List& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.setParameter("Output Information", 
			  NOX::Utils::Details +
			  NOX::Utils::OuterIteration + 
			  NOX::Utils::InnerIteration + 
			  NOX::Utils::Warning);

    // Set up the status tests
    NOX::StatusTest::NormF statusTestA(1.0e-8);
    NOX::StatusTest::MaxIters statusTestB(maxNewtonIters);
    NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR, statusTestA, 
				 statusTestB);

    // Create the solver
    NOX::Solver::Manager solver(grp, combo, nlParams); 

    // Solve the nonlinear system
    NOX::StatusTest::StatusType status = solver.solve();

    // Print the answer
    cout << "\n" << "-- Parameter List From Solver --" << "\n";
    solver.getParameterList().print(cout);
  
    // Get the answer
    grp = solver.getSolutionGroup();

    // Print the answer
    cout << "\n" << "-- Final Solution From Solver --" << "\n";
    grp.printSolution(0.0);
  }

  catch (string& s) {
    cout << s << endl;
  }
  catch (char *s) {
    cout << s << endl;
  }
  catch (exception& e) {
    cout << e.what() << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

  return 0;
}
