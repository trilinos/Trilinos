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

/*! \file  example_02.cpp
    \brief Shows how to solve a steady Burgers' optimal control problem using
	   the SimOpt interface.  We solve the control problem using Composite
           Step and trust regions.
*/

#include "example_02.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag = 0;

  // *** Example body.

  try {
    // Initialize full objective function.
    int nx      = 256;   // Set spatial discretization.
    RealT alpha = 1.e-3; // Set penalty parameter.
    RealT nu    = 1e-2;  // Viscosity parameter.
    Objective_BurgersControl<RealT> obj(alpha,nx);
    // Initialize equality constraints
    EqualityConstraint_BurgersControl<RealT> con(nx,nu);
    Teuchos::ParameterList list;
    list.sublist("SimOpt").sublist("Solve").set("Residual Tolerance",1.e2*ROL::ROL_EPSILON<RealT>());
    con.setSolveParameters(list);
    // Initialize iteration vectors.
    Teuchos::RCP<std::vector<RealT> > z_rcp  = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    Teuchos::RCP<std::vector<RealT> > gz_rcp = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    Teuchos::RCP<std::vector<RealT> > yz_rcp = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    for (int i=0; i<nx+2; i++) {
      (*z_rcp)[i]  = (RealT)rand()/(RealT)RAND_MAX;
      (*yz_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    ROL::StdVector<RealT> z(z_rcp);
    ROL::StdVector<RealT> gz(gz_rcp);
    ROL::StdVector<RealT> yz(yz_rcp);
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
    ROL::StdVector<RealT> u(u_rcp);
    ROL::StdVector<RealT> gu(gu_rcp);
    ROL::StdVector<RealT> yu(yu_rcp);
    Teuchos::RCP<ROL::Vector<RealT> > up  = Teuchos::rcp(&u,false);
    Teuchos::RCP<ROL::Vector<RealT> > gup = Teuchos::rcp(&u,false);
    Teuchos::RCP<ROL::Vector<RealT> > yup = Teuchos::rcp(&yu,false);

    Teuchos::RCP<std::vector<RealT> > c_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<std::vector<RealT> > l_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    ROL::StdVector<RealT> c(c_rcp);
    ROL::StdVector<RealT> l(l_rcp);

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
    Teuchos::RCP<std::vector<RealT> > p_rcp  = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    ROL::StdVector<RealT> p(p_rcp);
    Teuchos::RCP<ROL::Vector<RealT> > pp  = Teuchos::rcp(&p,false);
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > pobj = Teuchos::rcp(&obj,false);
    Teuchos::RCP<ROL::EqualityConstraint_SimOpt<RealT> > pcon = Teuchos::rcp(&con,false);
    ROL::Reduced_Objective_SimOpt<RealT> robj(pobj,pcon,up,pp);
    // Check derivatives.
    robj.checkGradient(z,z,yz,true,*outStream);
    robj.checkHessVec(z,z,yz,true,*outStream);

    // Get parameter list.
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-14);
    parlist->sublist("Status Test").set("Constraint Tolerance",1.e-14);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-16);
    parlist->sublist("Status Test").set("Iteration Limit",1000);
    // Declare ROL algorithm pointer.
    Teuchos::RCP<ROL::Algorithm<RealT> > algo;

    // Run optimization with Composite Step.
    algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Composite Step",*parlist,false));
    RealT zerotol = std::sqrt(ROL::ROL_EPSILON<RealT>());
    z.zero();
    con.solve(c,u,z,zerotol);
    c.zero(); l.zero();
    algo->run(x, g, l, c, obj, con, true, *outStream);
    Teuchos::RCP<ROL::Vector<RealT> > zCS = z.clone();
    zCS->set(z);

    // Run Optimization with Trust-Region algorithm.
    algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Trust Region",*parlist,false));
    z.zero();
    algo->run(z,robj,true,*outStream);

    // Check solutions.
    Teuchos::RCP<ROL::Vector<RealT> > err = z.clone();
    err->set(*zCS); err->axpy(-1.,z);
    errorFlag += ((err->norm()) > 1.e-8) ? 1 : 0;
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

