// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_04.cpp
    \brief Test default Constraint_SimOpt solve.
*/

#include "ROL_Stream.hpp"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>

#include "test_04.hpp"

typedef double RealT;
typedef H1VectorPrimal<RealT> PrimalStateVector;
typedef H1VectorDual<RealT> DualStateVector;
typedef L2VectorPrimal<RealT> PrimalControlVector;
typedef L2VectorDual<RealT> DualControlVector;
typedef H1VectorDual<RealT> PrimalConstraintVector;
typedef H1VectorPrimal<RealT> DualConstraintVector;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  bool print = (iprint>0);
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (print)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {
    /*************************************************************************/
    /************* INITIALIZE BURGERS FEM CLASS ******************************/
    /*************************************************************************/
    int nx      = 512;   // Set spatial discretization.
    RealT nl    = 1.0;   // Nonlinearity parameter (1 = Burgers, 0 = linear).
    RealT cH1   = 1.0;   // Scale for derivative term in H1 norm.
    RealT cL2   = 0.0;   // Scale for mass term in H1 norm.
    ROL::Ptr<BurgersFEM<RealT> > fem
      = ROL::makePtr<BurgersFEM<RealT>>(nx,nl,cH1,cL2);
    fem->test_inverse_mass(*outStream);
    fem->test_inverse_H1(*outStream);
    /*************************************************************************/
    /************* INITIALIZE SIMOPT CONSTRAINT ******************************/
    /*************************************************************************/
    bool hess = true;
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<Constraint_BurgersControl<RealT>>(fem,hess);
    /*************************************************************************/
    /************* INITIALIZE VECTOR STORAGE *********************************/
    /*************************************************************************/
    // INITIALIZE CONTROL VECTORS
    ROL::Ptr<std::vector<RealT> > z_ptr
      = ROL::makePtr<std::vector<RealT>>(nx+2, 0.0);
    ROL::Ptr<ROL::Vector<RealT> > zp
      = ROL::makePtr<PrimalControlVector>(z_ptr,fem);
    // INITIALIZE STATE VECTORS
    ROL::Ptr<std::vector<RealT> > u_ptr
      = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    ROL::Ptr<ROL::Vector<RealT> > up
      = ROL::makePtr<PrimalStateVector>(u_ptr,fem);
    // INITIALIZE CONSTRAINT VECTORS
    ROL::Ptr<std::vector<RealT> > c_ptr
      = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    ROL::Ptr<ROL::Vector<RealT> > cp
      = ROL::makePtr<PrimalConstraintVector>(c_ptr,fem);
    /*************************************************************************/
    /************* CHECK DERIVATIVES AND CONSISTENCY *************************/
    /*************************************************************************/
    RealT rnorm(0), cnorm(0);
    ROL::ParameterList list;
    RealT tol = std::sqrt(ROL::ROL_EPSILON<RealT>());
    list.sublist("SimOpt").sublist("Solve").set("Output Iteration History",print);
    list.sublist("SimOpt").sublist("Solve").set("Step Tolerance",ROL::ROL_EPSILON<RealT>());
    // Newton
    list.sublist("SimOpt").sublist("Solve").set("Absolute Residual Tolerance",tol);
    list.sublist("SimOpt").sublist("Solve").set("Solver Type",0);
    con->setSolveParameters(list);
    u_ptr->assign(nx, 1.0);
    con->solve(*cp,*up,*zp,tol);
    rnorm = cp->norm();
    con->value(*cp,*up,*zp,tol);
    cnorm = cp->norm();
    errorFlag += ((cnorm > tol) ? 1 : 0) + ((rnorm > tol) ? 1 : 0);
    *outStream << std::scientific << std::setprecision(8);
    *outStream << std::endl;
    *outStream << "Test SimOpt solve at feasible (u,z):" << std::endl;
    *outStream << "  Solver Residual = " << rnorm << std::endl;
    *outStream << "       ||c(u,z)|| = " << cnorm;
    *outStream << std::endl << std::endl;
    // Levenberg-Marquardt
    list.sublist("SimOpt").sublist("Solve").set("Absolute Residual Tolerance",1e-4*tol);
    list.sublist("SimOpt").sublist("Solve").set("Solver Type",1);
    con->setSolveParameters(list);
    u_ptr->assign(nx, 1.0);
    con->solve(*cp,*up,*zp,tol);
    rnorm = cp->norm();
    con->value(*cp,*up,*zp,tol);
    cnorm = cp->norm();
    errorFlag += ((cnorm > tol) ? 1 : 0) + ((rnorm > tol) ? 1 : 0);
    *outStream << std::scientific << std::setprecision(8);
    *outStream << std::endl;
    *outStream << "Test SimOpt solve at feasible (u,z):" << std::endl;
    *outStream << "  Solver Residual = " << rnorm << std::endl;
    *outStream << "       ||c(u,z)|| = " << cnorm;
    *outStream << std::endl << std::endl;
    // Composite Step
    list.sublist("SimOpt").sublist("Solve").set("Absolute Residual Tolerance",tol);
    list.sublist("SimOpt").sublist("Solve").set("Solver Type",2);
    con->setSolveParameters(list);
    u_ptr->assign(nx, 1.0);
    con->solve(*cp,*up,*zp,tol);
    rnorm = cp->norm();
    con->value(*cp,*up,*zp,tol);
    cnorm = cp->norm();
    tol *= 100.0;  // Probably will not need this when we have composite step
    errorFlag += ((cnorm > tol) ? 1 : 0) + ((rnorm > tol) ? 1 : 0);
    *outStream << std::scientific << std::setprecision(8);
    *outStream << std::endl;
    *outStream << "Test SimOpt solve at feasible (u,z):" << std::endl;
    *outStream << "  Solver Residual = " << rnorm << std::endl;
    *outStream << "       ||c(u,z)|| = " << cnorm;
    *outStream << std::endl << std::endl;
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
