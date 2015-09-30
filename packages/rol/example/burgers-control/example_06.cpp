// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  example_03.cpp
    \brief Shows how to solve a steady Burgers' optimal control problem using
           full-space methods.
*/

#include "example_06.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Example body.

  try {
    // Initialize finite element class.
    int nx      = 128;   // Set spatial discretization.
    RealT alpha = 1.e-3; // Set penalty parameter.
    RealT nu    = 1e-2;  // Viscosity parameter.
    Teuchos::RCP<BurgersFEM<RealT> > fem
      = Teuchos::rcp(new BurgersFEM<RealT>(nx,nu));
    fem->test_inverse_mass(*outStream);
    // Initialize objective function.
    Teuchos::RCP<std::vector<RealT> > ud_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<ROL::Vector<RealT> > ud = Teuchos::rcp(new L2Vector<RealT>(ud_rcp,fem));
    Objective_BurgersControl<RealT> obj(fem,ud,alpha);
    // Initialize equality constraints
    EqualityConstraint_BurgersControl<RealT> con(fem);
    // Initialize equality constraints
    std::vector<RealT> lo1(nx, 0.0), hi1(nx, 1.0);
    //std::vector<RealT> lo1(nx, -1.e6), hi1(nx, 1.e6);
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd1 = Teuchos::rcp(new L2BoundConstraint<RealT>(lo1,hi1,fem));
    //bnd1.deactivate();
    //std::vector<RealT> lo2(nx+2, -0.1*ROL::ROL_OVERFLOW), hi2(nx+2, 0.1*ROL::ROL_OVERFLOW);
    std::vector<RealT> lo2(nx+2,0.0), hi2(nx+2,2.0);
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd2 = Teuchos::rcp(new L2BoundConstraint<RealT>(lo2,hi2,fem));
    //bnd2.deactivate();
    ROL::BoundConstraint_SimOpt<RealT> bnd(bnd1,bnd2);
    // Initialize iteration vectors.
    Teuchos::RCP<std::vector<RealT> > z_rcp  = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    Teuchos::RCP<std::vector<RealT> > gz_rcp = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    Teuchos::RCP<std::vector<RealT> > yz_rcp = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    for (int i=0; i<nx+2; i++) {
      (*z_rcp)[i]  = (RealT)rand()/(RealT)RAND_MAX;
      (*yz_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    L2Vector<RealT> z(z_rcp,fem);
    L2Vector<RealT> gz(gz_rcp,fem);
    L2Vector<RealT> yz(yz_rcp,fem);
    Teuchos::RCP<ROL::Vector<RealT> > zp  = Teuchos::rcp(&z,false);
    Teuchos::RCP<ROL::Vector<RealT> > gzp = Teuchos::rcp(&z,false);
    Teuchos::RCP<ROL::Vector<RealT> > yzp = Teuchos::rcp(&yz,false);

    Teuchos::RCP<std::vector<RealT> > u_rcp  = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<std::vector<RealT> > gu_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<std::vector<RealT> > yu_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    for (int i=0; i<nx; i++) {
      (*u_rcp)[i]  = (RealT)rand()/(RealT)RAND_MAX;
      (*yu_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    L2Vector<RealT> u(u_rcp,fem);
    L2Vector<RealT> gu(gu_rcp,fem);
    L2Vector<RealT> yu(yu_rcp,fem);
    Teuchos::RCP<ROL::Vector<RealT> > up  = Teuchos::rcp(&u,false);
    Teuchos::RCP<ROL::Vector<RealT> > gup = Teuchos::rcp(&u,false);
    Teuchos::RCP<ROL::Vector<RealT> > yup = Teuchos::rcp(&yu,false);

    Teuchos::RCP<std::vector<RealT> > c_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<std::vector<RealT> > l_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    L2Vector<RealT> c(c_rcp,fem);
    L2Vector<RealT> l(l_rcp,fem);
    for (int i=0; i<nx; i++) {
      (*l_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }

    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> g(gup,gzp);
    ROL::Vector_SimOpt<RealT> y(yup,yzp);

    // Check derivatives.
    obj.checkGradient(x,x,y,true,*outStream);
    obj.checkHessVec(x,x,y,true,*outStream);

    con.checkApplyJacobian(x,y,c,true,*outStream);
    con.checkApplyAdjointHessian(x,yu,y,x,true,*outStream);

    con.checkSolve(u,z,c,true,*outStream);
    con.checkAdjointConsistencyJacobian_1(c,yu,u,z,true,*outStream);
    con.checkAdjointConsistencyJacobian_2(c,yz,u,z,true,*outStream);
    con.checkInverseJacobian_1(yu,l,u,z,true,*outStream);
    con.checkInverseAdjointJacobian_1(yu,l,u,z,true,*outStream);
    *outStream << "\n";

    ROL::AugmentedLagrangian<RealT> augLag(obj,con,x,c);
    augLag.updateMultipliers(l,1.0);
    augLag.checkGradient(x, y, true, *outStream);
    augLag.checkHessVec(x, y, true, *outStream);

    // Optimization 
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, Teuchos::Ptr<Teuchos::ParameterList>(&*parlist) );
    // Define status test.
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-6);
    parlist->sublist("Status Test").set("Constraint Tolerance",1.e-6);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-10);
    parlist->sublist("Status Test").set("Iteration Limit",100);
    ROL::StatusTestSQP<RealT> status(*parlist);
    // Define step.
    ROL::AugmentedLagrangianStep<RealT> step(*parlist);
    // Define algorithm.
    ROL::DefaultAlgorithm<RealT> algo(step,status,false);
    // Run Algorithm
    RealT zerotol = 0.0;
    //z.scale(50.0);
    con.solve(u,z,zerotol);
    obj.gradient_1(gu,u,z,zerotol);
    gu.scale(-1.0);
    con.applyInverseAdjointJacobian_1(l,gu,u,z,zerotol);
    gu.zero();
    c.zero();
    algo.run(x, g, l, c, obj, con, bnd, true, *outStream);

//    for ( int i = 0; i < nx+2; i++ ) {
//      std::cout << std::scientific << std::setprecision(10);
//      std::cout << std::setw(20) << std::left << (*z_rcp)[i];
//      if ( i == 0 ) {
//        std::cout << std::setw(20) << std::left << 1.0;
//      }
//      if ( i != 0 && i != nx+1 ) {
//        std::cout << std::setw(20) << std::left << (*u_rcp)[i-1];
//      }
//      if ( i == nx+1 ) {
//        std::cout << std::setw(20) << std::left << 0.0;
//      }
//      std::cout << "\n";
//    }
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
