// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*  Example of solving a problem with bound and inequality constraints
 *
 */

#define OPTIMIZATION_PROBLEM_REFACTOR 

#include "ROL_OptimizationSolver.hpp"

#include "ROL_RandomVector.hpp"
#include "ROL_StdObjective.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"


/* OBJECTIVE FUNCTION */

template<class Real> 
class ObjectiveRosenbrock : public ROL::StdObjective<Real> {
public:
  ObjectiveRosenbrock(void) {}

  Real value(const std::vector<Real> &x, Real &tol) {
    const Real one(1), alpha(100);
    Real val = alpha * std::pow(std::pow(x[0], 2) - x[1], 2)
               + std::pow(x[0] - one, 2);
    return val;
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    const Real one(1), two(2), alpha(100);
    g[0] = two*alpha*(std::pow(x[0], 2) - x[1]) * two*x[0] + two*(x[0]-one);
    g[1] = -two*alpha*(std::pow(x[0], 2) - x[1]);
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    const Real two(2), three(3), alpha(100);
    Real h11 = two*two*three*alpha*std::pow(x[0], 2) - two*two*alpha*x[1] + two;
    Real h12 = -two*two*alpha*x[0];
    Real h21 = h12;
    Real h22 = two*alpha;
    hv[0] = h11*v[0] + h12*v[1];
    hv[1] = h21*v[0] + h22*v[1];
  }

}; // class ObjectiveRosenbrock


int main(int argc, char *argv[]) {

   

  typedef double RealT;
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);


  int errorFlag   = 0;

  try {
    ROL::ParameterList parlist;
    parlist.sublist("General").sublist("Secant").set("Use as Hessian",false);
    parlist.sublist("Step").set("Type","Trust Region");
    parlist.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");

    ROL::Ptr<std::vector<RealT> > x_ptr  = ROL::makePtr<std::vector<RealT>>(2);
    ROL::Ptr<ROL::Vector<RealT> > x      = ROL::makePtr<ROL::StdVector<RealT>>(x_ptr);
    (*x_ptr)[0] = static_cast<RealT>(-3);
    (*x_ptr)[1] = static_cast<RealT>(-4);

    ROL::Ptr<ROL::Objective<RealT> > obj = ROL::makePtr<ObjectiveRosenbrock<RealT>>();

    ROL::OptimizationProblem<RealT> problem( obj, x );
    problem.check(*outStream);

    ROL::OptimizationSolver<RealT> solver( problem, parlist );
    solver.solve(*outStream); 

    *outStream << "x_opt = [" << (*x_ptr)[0] << ", " << (*x_ptr)[1] << "]" << std::endl;
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;


  return 0;
}





