// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_02.cpp
    \brief Shows how to solve an optimal control problem constrained by
           unsteady Burgers' equation with the SimOpt interface and
           hand-coded time integration loops.
*/

#include "example_02.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_OptimizationSolver.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  // *** Example body.

  try {
    // Initialize full objective function.
    int nx      = 80;    // Set spatial discretization.
    int nt      = 80;    // Set temporal discretization.
    RealT T     = 1.0;   // Set end time.
    RealT alpha = 5e-2;  // Set penalty parameter.
    RealT nu    = 1e-2;  // Set viscosity parameter.
    Objective_BurgersControl<RealT> obj(alpha,nx,nt,T);

    // Initialize equality constraints
    Constraint_BurgersControl<RealT> con(nx, nt, T, nu);

    // Initialize control, state, constraint and Lagrange multiplier vectors.
    ROL::Ptr<std::vector<RealT> > z_ptr = ROL::makePtr<std::vector<RealT>>((nx+2)*(nt+1), 1.0);
    ROL::StdVector<RealT> z(z_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zp = ROL::makePtrFromRef(z);

    ROL::Ptr<std::vector<RealT> > u_ptr = ROL::makePtr<std::vector<RealT>>(nx*nt, 1.0);
    ROL::StdVector<RealT> u(u_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up = ROL::makePtrFromRef(u);

    ROL::Ptr<std::vector<RealT> > c_ptr = ROL::makePtr<std::vector<RealT>>(nx*nt, 1.0);
    ROL::StdVector<RealT> c(c_ptr);

    ROL::Ptr<std::vector<RealT> > l_ptr = ROL::makePtr<std::vector<RealT>>(nx*nt, 1.0);
    ROL::StdVector<RealT> l(l_ptr);
    ROL::Ptr<ROL::Vector<RealT> > lp = ROL::makePtrFromRef(l);

    // Initialize reduced objective function.
    ROL::Ptr<ROL::Objective_SimOpt<RealT> >  pobj = ROL::makePtrFromRef(obj);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pcon = ROL::makePtrFromRef(con);
    ROL::Reduced_Objective_SimOpt<RealT> robj(pobj,pcon,up,zp,lp);
    ROL::Ptr<ROL::Objective<RealT> > probj = ROL::makePtrFromRef(robj);

    // Get input parameter list.
    auto parlist = ROL::getParametersFromXmlFile( "example_02.xml" );

    // Run optimization.

    // Solve in reduced space.
    z.zero();
    ROL::OptimizationProblem<RealT> problem_rs(probj, zp);
    problem_rs.check(*outStream);
    ROL::OptimizationSolver<RealT> solver_rs(problem_rs, *parlist);
    std::clock_t timer_rs = std::clock();
    solver_rs.solve(*outStream);
    *outStream << "Reduced-space formulation required " << (std::clock()-timer_rs)/(RealT)CLOCKS_PER_SEC
               << " seconds.\n";
    ROL::Ptr<ROL::Vector<RealT> > zTR = z.clone();
    zTR->set(z);

    // Solve in full space.
    ROL::Vector_SimOpt<RealT> x(up,zp);  // Build full-space vector.
    ROL::Ptr<ROL::Vector<RealT> > xp = ROL::makePtrFromRef(x);
    ROL::Elementwise::Fill<RealT> setFunc(0.25);
    x.applyUnary(setFunc);
    ROL::OptimizationProblem<RealT> problem_fs(pobj, xp, pcon, lp);
    problem_fs.check(*outStream);
    ROL::OptimizationSolver<RealT> solver_fs(problem_fs, *parlist);
    std::clock_t timer_fs = std::clock();
    solver_fs.solve(*outStream);
    *outStream << "Full-space formulation required " << (std::clock()-timer_fs)/(RealT)CLOCKS_PER_SEC
               << " seconds.\n";

    // Compute error between solutions
    ROL::Ptr<ROL::Vector<RealT> > err = z.clone();
    err->set(*zTR); err->axpy(-1.,z);
    errorFlag += (err->norm() > 1.e-4) ? 1 : 0;
    if (errorFlag) {
      *outStream << "\n\nControl error = " << err->norm() << "\n";
    }

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

}

