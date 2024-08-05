// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_02.cpp
    \brief Shows how to solve a steady Burgers' optimal control problem using
    the SimOpt interface.  We solve the control problem using Composite
           Step and trust regions.
*/

#include "ROL_TypeU_TrustRegionAlgorithm.hpp"
#include "ROL_TypeE_CompositeStepAlgorithm.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Stream.hpp"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <algorithm>

#include "example_02.hpp"

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
    int nx      = 256;   // Set spatial discretization.
    RealT alpha = 1.e-3; // Set penalty parameter.
    RealT nu    = 1e-2;  // Viscosity parameter.
    Objective_BurgersControl<RealT> obj(alpha,nx);
    // Initialize equality constraints
    Constraint_BurgersControl<RealT> con(nx,nu);
    ROL::ParameterList list;
    list.sublist("SimOpt").sublist("Solve").set("Absolute Residual Tolerance",1.e2*ROL::ROL_EPSILON<RealT>());
    con.setSolveParameters(list);
    // Initialize iteration vectors.
    ROL::Ptr<std::vector<RealT> > z_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2, 1.0);
    ROL::Ptr<std::vector<RealT> > gz_ptr = ROL::makePtr<std::vector<RealT>>(nx+2, 1.0);
    ROL::Ptr<std::vector<RealT> > yz_ptr = ROL::makePtr<std::vector<RealT>>(nx+2, 1.0);
    for (int i=0; i<nx+2; i++) {
      (*z_ptr)[i]  = (RealT)rand()/(RealT)RAND_MAX;
      (*yz_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    ROL::StdVector<RealT> z(z_ptr);
    ROL::StdVector<RealT> gz(gz_ptr);
    ROL::StdVector<RealT> yz(yz_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zp  = ROL::makePtrFromRef(z);
    ROL::Ptr<ROL::Vector<RealT> > gzp = ROL::makePtrFromRef(z);
    ROL::Ptr<ROL::Vector<RealT> > yzp = ROL::makePtrFromRef(yz);

    ROL::Ptr<std::vector<RealT> > u_ptr  = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    ROL::Ptr<std::vector<RealT> > gu_ptr = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    ROL::Ptr<std::vector<RealT> > yu_ptr = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    for (int i=0; i<nx; i++) {
      (*u_ptr)[i]  = (RealT)rand()/(RealT)RAND_MAX;
      (*yu_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    ROL::StdVector<RealT> u(u_ptr);
    ROL::StdVector<RealT> gu(gu_ptr);
    ROL::StdVector<RealT> yu(yu_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up  = ROL::makePtrFromRef(u);
    ROL::Ptr<ROL::Vector<RealT> > gup = ROL::makePtrFromRef(gu);
    ROL::Ptr<ROL::Vector<RealT> > yup = ROL::makePtrFromRef(yu);

    ROL::Ptr<std::vector<RealT> > c_ptr = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    ROL::Ptr<std::vector<RealT> > l_ptr = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    ROL::StdVector<RealT> c(c_ptr);
    ROL::StdVector<RealT> l(l_ptr);

    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> g(gup,gzp);
    ROL::Vector_SimOpt<RealT> y(yup,yzp);

    // Check derivatives.
    obj.checkGradient(x,x,y,true,*outStream);
    obj.checkHessVec(x,x,y,true,*outStream);
    con.checkApplyJacobian(x,y,c,true,*outStream);
    con.checkApplyAdjointJacobian(x,yu,c,x,true,*outStream);
    con.checkApplyAdjointHessian(x,yu,y,x,true,*outStream);

    // Initialize reduced objective function.
    ROL::Ptr<std::vector<RealT> > p_ptr  = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    ROL::StdVector<RealT> p(p_ptr);
    ROL::Ptr<ROL::Vector<RealT> > pp              = ROL::makePtrFromRef(p);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > pobj  = ROL::makePtrFromRef(obj);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pcon = ROL::makePtrFromRef(con);
    ROL::Reduced_Objective_SimOpt<RealT> robj(pobj,pcon,up,zp,pp);
    // Check derivatives.
    robj.checkGradient(z,z,yz,true,*outStream);
    robj.checkHessVec(z,z,yz,true,*outStream);

    // Get parameter list.
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-14);
    parlist->sublist("Status Test").set("Constraint Tolerance",1.e-14);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-16);
    parlist->sublist("Status Test").set("Iteration Limit",1000);

    // Run equality-constrained optimization.
    RealT zerotol = std::sqrt(ROL::ROL_EPSILON<RealT>());
    z.zero();
    con.solve(c,u,z,zerotol);
    c.zero(); l.zero();
    {
      // Define algorithm.
      ROL::TypeE::CompositeStepAlgorithm<RealT> algo(*parlist);
      // Run Algorithm
      algo.run(x, obj, con, l, *outStream);
    }
    ROL::Ptr<ROL::Vector<RealT> > zCS = z.clone();
    zCS->set(z);

    // Run unconstrained optimization.
    z.zero();
    {
      // Define algorithm.
      ROL::TypeU::TrustRegionAlgorithm<RealT> algo(*parlist);
      // Run Algorithm
      algo.run(z, z.dual(), robj, *outStream);
    }

    // Check solutions.
    ROL::Ptr<ROL::Vector<RealT> > err = z.clone();
    err->set(*zCS); err->axpy(-1.,z);
    errorFlag += ((err->norm()) > 1.e-8) ? 1 : 0;
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

