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
    \brief Shows how to solve an optimal control problem constrained by 
           unsteady Burgers' equation with the SimOpt interface.
*/

#include "example_03.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
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
    // Initialize iteration vectors.
    ROL::Ptr<std::vector<RealT> > z_ptr  = ROL::makePtr<std::vector<RealT>>((nx+2)*(nt+1), 1.0);
    ROL::Ptr<std::vector<RealT> > gz_ptr = ROL::makePtr<std::vector<RealT>>((nx+2)*(nt+1), 1.0);
    ROL::Ptr<std::vector<RealT> > yz_ptr = ROL::makePtr<std::vector<RealT>>((nx+2)*(nt+1), 1.0);
    for (int i=0; i<(nx+2)*(nt+1); i++) {
      (*z_ptr)[i]  = (RealT)rand()/(RealT)RAND_MAX;
      (*yz_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    ROL::StdVector<RealT> z(z_ptr);
    ROL::StdVector<RealT> gz(gz_ptr);
    ROL::StdVector<RealT> yz(yz_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zp  = ROL::makePtrFromRef(z);
    ROL::Ptr<ROL::Vector<RealT> > gzp = ROL::makePtrFromRef(gz);
    ROL::Ptr<ROL::Vector<RealT> > yzp = ROL::makePtrFromRef(yz);

    ROL::Ptr<std::vector<RealT> > u_ptr  = ROL::makePtr<std::vector<RealT>>(nx*nt, 1.0);
    ROL::Ptr<std::vector<RealT> > gu_ptr = ROL::makePtr<std::vector<RealT>>(nx*nt, 1.0);
    ROL::Ptr<std::vector<RealT> > yu_ptr = ROL::makePtr<std::vector<RealT>>(nx*nt, 1.0);
    for (int i=0; i<nx*nt; i++) {
      (*u_ptr)[i]  = (RealT)rand()/(RealT)RAND_MAX;
      (*yu_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    ROL::StdVector<RealT> u(u_ptr);
    ROL::StdVector<RealT> gu(gu_ptr);
    ROL::StdVector<RealT> yu(yu_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up  = ROL::makePtrFromRef(u);
    ROL::Ptr<ROL::Vector<RealT> > gup = ROL::makePtrFromRef(gu);
    ROL::Ptr<ROL::Vector<RealT> > yup = ROL::makePtrFromRef(yu);

    ROL::Ptr<std::vector<RealT> > c_ptr = ROL::makePtr<std::vector<RealT>>(nx*nt, 1.0);
    ROL::Ptr<std::vector<RealT> > l_ptr = ROL::makePtr<std::vector<RealT>>(nx*nt, 1.0);
    ROL::StdVector<RealT> c(c_ptr);
    ROL::StdVector<RealT> l(l_ptr);

    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> g(gup,gzp);
    ROL::Vector_SimOpt<RealT> y(yup,yzp);
    // Check derivatives.
    obj.checkGradient(x,x,y,true,*outStream);
    obj.checkHessVec(x,x,y,true,*outStream);
    con.checkApplyJacobian(x,y,c,true,*outStream);
    //con.checkApplyAdjointJacobian(x,yu,c,x,true,*outStream);
    con.checkApplyAdjointHessian(x,yu,y,x,true,*outStream);
    // Check consistency of Jacobians and adjoint Jacobians.
    con.checkAdjointConsistencyJacobian_1(c,yu,u,z,true,*outStream);
    con.checkAdjointConsistencyJacobian_2(c,yz,u,z,true,*outStream);
    // Check consistency of solves.
    con.checkSolve(u,z,c,true,*outStream);
    con.checkInverseJacobian_1(c,yu,u,z,true,*outStream);
    con.checkInverseAdjointJacobian_1(yu,c,u,z,true,*outStream);

    // Initialize reduced objective function.
    ROL::Ptr<std::vector<RealT> > p_ptr  = ROL::makePtr<std::vector<RealT>>(nx*nt, 1.0);
    ROL::StdVector<RealT> p(p_ptr);
    ROL::Ptr<ROL::Vector<RealT> > pp              = ROL::makePtrFromRef(p);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > pobj  = ROL::makePtrFromRef(obj);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pcon = ROL::makePtrFromRef(con);
    ROL::Reduced_Objective_SimOpt<RealT> robj(pobj,pcon,up,zp,pp);
    // Check derivatives.
    robj.checkGradient(z,z,yz,true,*outStream);
    robj.checkHessVec(z,z,yz,true,*outStream);
    // Get input parameter list.
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-10);
    parlist->sublist("Status Test").set("Constraint Tolerance",1.e-10);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-16);
    parlist->sublist("Status Test").set("Iteration Limit",100);
    // Build Algorithm pointer.
    ROL::Ptr<ROL::Algorithm<RealT> > algo;

    // Solve using trust regions.
    algo = ROL::makePtr<ROL::Algorithm<RealT>>("Trust Region",*parlist,false);
    z.zero();
    std::clock_t timer_tr = std::clock();
    algo->run(z,robj,true,*outStream);
    *outStream << "Trust-Region Newton required " << (std::clock()-timer_tr)/(RealT)CLOCKS_PER_SEC
               << " seconds.\n";
    ROL::Ptr<ROL::Vector<RealT> > zTR = z.clone();
    zTR->set(z);

    // Solve using a composite step method.
    algo = ROL::makePtr<ROL::Algorithm<RealT>>("Composite Step",*parlist,false);
    x.zero();
    ROL::Elementwise::Fill<RealT> setFunc(0.25);
    x.applyUnary(setFunc);
    std::clock_t timer_cs = std::clock();
    algo->run(x,g,l,c,obj,con,true,*outStream);
    *outStream << "Composite Step required " << (std::clock()-timer_cs)/(RealT)CLOCKS_PER_SEC
               << " seconds.\n";

    // Compute error between solutions
    ROL::Ptr<ROL::Vector<RealT> > err = z.clone();
    err->set(*zTR); err->axpy(-1.,z);
    errorFlag += (err->norm() > 1.e-4) ? 1 : 0;
    if (errorFlag) {
      *outStream << "\n\nControl error = " << err->norm() << "\n";
    }

//    std::ofstream control;
//    control.open("control.txt");
//    for (int t = 0; t < nt+1; t++) {
//      for (int n = 0; n < nx+2; n++) {
//        control << (RealT)t/(RealT)nt       << "  " 
//                << (RealT)n/((RealT)(nx+1)) << "  " 
//                << (*z_ptr)[t*(nx+2)+n]     << "\n";
//      }
//    } 
//    control.close();
//
//    std::ofstream state;
//    state.open("state.txt");
//    for (int t = 0; t < nt; t++) {
//      for (int n = 0; n < nx; n++) {
//        state << (RealT)(t+1)/(RealT)nt       << "  " 
//              << (RealT)(n+1)/((RealT)(nx+1)) << "  " 
//              << (*u_ptr)[t*nx+n]             << "\n";
//      }
//    } 
//    state.close();
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

