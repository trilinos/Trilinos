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

/*! \file  example_01.cpp
    \brief Shows how to solve a golf putting control problem.
*/

#include "example_01.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::SharedPointer<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makeSharedFromRef(std::cout);
  else
    outStream = ROL::makeSharedFromRef(bhs);

  int errorFlag = 0;

  // *** Example body.

  try {
    // Read in parameterlist
    std::string filename = "input.xml";
    ROL::SharedPointer<Teuchos::ParameterList> parlist = ROL::makeShared<Teuchos::ParameterList>();
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Physical parameters
    RealT g  = 9.834;   // Acceleration due to gravity (m/s^2)
    RealT m  = 0.04593; // Mass of golf ball (kg)
    RealT mu = 0.07;    // Coefficient of friction: range (0.06,0.20)
    // User input parameters
    RealT x0 = parlist->sublist("Problem").get("Initial x-location",  1.0); // Initial x-location of ball
    RealT y0 = parlist->sublist("Problem").get("Initial y-location",  2.0); // Initial y-location of ball
    RealT xn = parlist->sublist("Problem").get("Hole x-location",     1.0); // x-location of hole
    RealT yn = parlist->sublist("Problem").get("Hole y-location",    -2.0); // y-location of hole
    RealT Rg = parlist->sublist("Problem").get("Green radius",        4.0); // Green radius (m)
    RealT sf = parlist->sublist("Problem").get("Final speed",         0.1); // Final speed
    int    n = parlist->sublist("Problem").get("Number of time steps", 50); // Number of time steps

    // Initialize iteration vectors.
    ROL::SharedPointer<std::vector<RealT> > u_ptr    = ROL::makeShared<std::vector<RealT>>(9*n+10);
    ROL::SharedPointer<std::vector<RealT> > ulb_ptr  = ROL::makeShared<std::vector<RealT>>(9*n+10);
    ROL::SharedPointer<std::vector<RealT> > uub_ptr  = ROL::makeShared<std::vector<RealT>>(9*n+10);
    ROL::SharedPointer<std::vector<RealT> > z_ptr    = ROL::makeShared<std::vector<RealT>>(2);
    ROL::SharedPointer<std::vector<RealT> > zlb_ptr  = ROL::makeShared<std::vector<RealT>>(2);
    ROL::SharedPointer<std::vector<RealT> > zub_ptr  = ROL::makeShared<std::vector<RealT>>(2);
    ROL::SharedPointer<std::vector<RealT> > emul_ptr = ROL::makeShared<std::vector<RealT>>(9*n+10);
    ROL::SharedPointer<std::vector<RealT> > imul_ptr = ROL::makeShared<std::vector<RealT>>(n+1);
    ROL::SharedPointer<std::vector<RealT> > ilb_ptr  = ROL::makeShared<std::vector<RealT>>(n+1);
    ROL::SharedPointer<std::vector<RealT> > iub_ptr  = ROL::makeShared<std::vector<RealT>>(n+1);
    ROL::SharedPointer<ROL::Vector<RealT> > up       = ROL::makeShared<ROL::StdVector<RealT>>(u_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > ulbp     = ROL::makeShared<ROL::StdVector<RealT>>(ulb_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > uubp     = ROL::makeShared<ROL::StdVector<RealT>>(uub_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > zp       = ROL::makeShared<ROL::StdVector<RealT>>(z_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > zlbp     = ROL::makeShared<ROL::StdVector<RealT>>(zlb_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > zubp     = ROL::makeShared<ROL::StdVector<RealT>>(zub_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > emul     = ROL::makeShared<ROL::StdVector<RealT>>(emul_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > imul     = ROL::makeShared<ROL::StdVector<RealT>>(imul_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > ilbp     = ROL::makeShared<ROL::StdVector<RealT>>(ilb_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > iubp     = ROL::makeShared<ROL::StdVector<RealT>>(iub_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > x        = ROL::makeShared<ROL::Vector_SimOpt<RealT>>(up,zp);

    // Fill initial state
    RealT N(n), one(1), T(1.5);
    for (int i = 0; i < n+1; ++i) {
      RealT I(i);
      (*u_ptr)[i]         = I/N * x0 + (one-I/N) * xn;
      (*u_ptr)[n+1+i]     = I/N * y0 + (one-I/N) * yn;
      (*u_ptr)[3*n+3+i]   = (xn - x0)/T;
      (*u_ptr)[4*n+4+i]   = (yn - y0)/T;
    }
    (*u_ptr)[9*n+9] = T;
    (*z_ptr)[0] = (xn - x0)/T;
    (*z_ptr)[1] = (yn - y0)/T;

    // Fill bounds
    for (int i = 0; i < 9*n+10; ++i) {
      (*ulb_ptr)[i] = ROL::ROL_NINF<RealT>();
      (*uub_ptr)[i] = ROL::ROL_INF<RealT>();
    }
    (*ulb_ptr)[9*n+9] = static_cast<RealT>(0);
    for (int i = 0; i < 2; ++i) {
      (*zlb_ptr)[i] = ROL::ROL_NINF<RealT>();
      (*zub_ptr)[i] = ROL::ROL_INF<RealT>();
    }
    for (int i = 0; i < n+1; ++i) {
      (*ilb_ptr)[i] = static_cast<RealT>(0);
      (*iub_ptr)[i] = Rg*Rg;
    }

    // Initialize bound constraints
    ROL::SharedPointer<ROL::Bounds<RealT> > ubnd
      = ROL::makeShared<ROL::Bounds<RealT>>(ulbp,uubp);
    ROL::SharedPointer<ROL::Bounds<RealT> > zbnd
      = ROL::makeShared<ROL::Bounds<RealT>>(zlbp,zubp);
    zbnd->deactivate();
    ROL::SharedPointer<ROL::BoundConstraint<RealT> > xbnd
      = ROL::makeShared<ROL::BoundConstraint_SimOpt<RealT>>(ubnd,zbnd);
    ROL::SharedPointer<ROL::Bounds<RealT> > ibnd
      = ROL::makeShared<ROL::Bounds<RealT>>(ilbp,iubp);

    // Dynamic constraints
    ROL::SharedPointer<PuttingConstraint<RealT> > econ
      = ROL::makeShared<PuttingConstraint<RealT>>(g,m,x0,y0,xn,yn,mu);

    // Green constraints
    ROL::SharedPointer<GreenConstraint<RealT> > icon
      = ROL::makeShared<GreenConstraint<RealT>>();

    // Final speed objective
    RealT target = sf*sf;
    ROL::SharedPointer<PuttingObjective<RealT> > obj
      = ROL::makeShared<PuttingObjective<RealT>>(target);

    // Initialize optimization problem
    ROL::SharedPointer<ROL::OptimizationProblem<RealT> > problem;
    bool useReduced = parlist->sublist("Problem").get("Use reduced space",false);
    if (!useReduced) {
      bool useGreenCon = parlist->sublist("Problem").get("Use green constraint",false);
      if (!useGreenCon) {
        problem = ROL::makeShared<ROL::OptimizationProblem<RealT>>(obj,x,xbnd,econ,emul);
      }
      else {
        problem = ROL::makeShared<ROL::OptimizationProblem<RealT>>(obj,x,xbnd,econ,emul,icon,imul,ibnd);
      }
    }
    else {
      econ->setSolveParameters(*parlist);
      ROL::SharedPointer<ROL::SimController<RealT> > stateStore
        = ROL::makeShared<ROL::SimController<RealT>>();
      ROL::SharedPointer<ROL::Objective<RealT> > robj
        = ROL::makeShared<ROL::Reduced_Objective_SimOpt<RealT>>(obj,econ,stateStore,up,zp,emul,true,false);
      problem = ROL::makeShared<ROL::OptimizationProblem<RealT>>(robj,zp);
    }

    // Check derivatives
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if (checkDeriv) {
      problem->check(*outStream);
    }

    // Initialize optimization solver
    ROL::OptimizationSolver<RealT> solver(*problem,*parlist);

    // Solve putting control problem
    solver.solve(*outStream);

    // Print optimal control (initial velocity vector)
    *outStream << "Initial x-velocity: " << (*z_ptr)[0]     << std::endl;
    *outStream << "Initial y-velocity: " << (*z_ptr)[1]     << std::endl;
    *outStream << "Final time: "         << (*u_ptr)[9*n+9] << std::endl;

    // Print optimal trajectory
    bool printTrajectory = parlist->sublist("Problem").get("Print trajectory",true);
    if (printTrajectory) {
      std::ofstream xFile, vFile, aFile;
      xFile.open("position.txt");
      vFile.open("velocity.txt");
      aFile.open("acceleration.txt");
      for (int i = 0; i < n+1; ++i) {
        xFile << std::scientific << std::setprecision(8) << std::setw(12) << std::left;
        xFile << (*u_ptr)[i]       << "  " << (*u_ptr)[n+1+i]   << "  " << (*u_ptr)[2*n+2+i] << std::endl;
        vFile << std::scientific << std::setprecision(8) << std::setw(12) << std::left;
        vFile << (*u_ptr)[3*n+3+i] << "  " << (*u_ptr)[4*n+4+i] << "  " << (*u_ptr)[5*n+5+i] << std::endl;
        aFile << std::scientific << std::setprecision(8) << std::setw(12) << std::left;
        aFile << (*u_ptr)[6*n+6+i] << "  " << (*u_ptr)[7*n+7+i] << "  " << (*u_ptr)[8*n+8+i] << std::endl;
      }
      xFile.close();
      vFile.close();
      aFile.close();
    }

/*
    RealT tol(1e-8);
    std::ofstream jacFile;
    jacFile.open("jacobian.txt");
    ROL::SharedPointer<ROL::Vector<RealT> > jv = emul->clone();
    for (int i = 0; i < 9*n+10; ++i) {
      jacFile << std::scientific << std::setprecision(8) << std::setw(12) << std::left;
      for (int j = 0; j < 9*n+10; ++j) {
        econ->applyJacobian_1(*jv,*(up->basis(i)),*up,*zp,tol);
        jacFile << jv->dot(*(jv->basis(j))) << "  ";
      }
      jacFile << std::endl;
    }
    jacFile.close();

    ROL::SharedPointer<ROL::Vector<RealT> > ajv = up->clone();
    std::ofstream ajacFile;
    ajacFile.open("adjoint-jacobian.txt");
    for (int j = 0; j < 9*n+10; ++j) {
      ajacFile << std::scientific << std::setprecision(8) << std::setw(12) << std::left;
      for (int i = 0; i < 9*n+10; ++i) {
        econ->applyAdjointJacobian_1(*ajv,*(emul->basis(j)),*up,*zp,tol);
        ajacFile << ajv->dot(*(ajv->basis(i))) << "  ";
      }
      ajacFile << std::endl;
    }
    ajacFile.close();
*/

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
