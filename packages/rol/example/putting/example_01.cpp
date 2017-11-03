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
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag = 0;

  // *** Example body.

  try {
    // Read in parameterlist
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
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
    Teuchos::RCP<std::vector<RealT> > u_rcp    = Teuchos::rcp(new std::vector<RealT>(9*n+10));
    Teuchos::RCP<std::vector<RealT> > ulb_rcp  = Teuchos::rcp(new std::vector<RealT>(9*n+10));
    Teuchos::RCP<std::vector<RealT> > uub_rcp  = Teuchos::rcp(new std::vector<RealT>(9*n+10));
    Teuchos::RCP<std::vector<RealT> > z_rcp    = Teuchos::rcp(new std::vector<RealT>(2));
    Teuchos::RCP<std::vector<RealT> > zlb_rcp  = Teuchos::rcp(new std::vector<RealT>(2));
    Teuchos::RCP<std::vector<RealT> > zub_rcp  = Teuchos::rcp(new std::vector<RealT>(2));
    Teuchos::RCP<std::vector<RealT> > emul_rcp = Teuchos::rcp(new std::vector<RealT>(9*n+10));
    Teuchos::RCP<std::vector<RealT> > imul_rcp = Teuchos::rcp(new std::vector<RealT>(n+1));
    Teuchos::RCP<std::vector<RealT> > ilb_rcp  = Teuchos::rcp(new std::vector<RealT>(n+1));
    Teuchos::RCP<std::vector<RealT> > iub_rcp  = Teuchos::rcp(new std::vector<RealT>(n+1));
    Teuchos::RCP<ROL::Vector<RealT> > up       = Teuchos::rcp(new ROL::StdVector<RealT>(u_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > ulbp     = Teuchos::rcp(new ROL::StdVector<RealT>(ulb_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > uubp     = Teuchos::rcp(new ROL::StdVector<RealT>(uub_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > zp       = Teuchos::rcp(new ROL::StdVector<RealT>(z_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > zlbp     = Teuchos::rcp(new ROL::StdVector<RealT>(zlb_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > zubp     = Teuchos::rcp(new ROL::StdVector<RealT>(zub_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > emul     = Teuchos::rcp(new ROL::StdVector<RealT>(emul_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > imul     = Teuchos::rcp(new ROL::StdVector<RealT>(imul_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > ilbp     = Teuchos::rcp(new ROL::StdVector<RealT>(ilb_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > iubp     = Teuchos::rcp(new ROL::StdVector<RealT>(iub_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > x        = Teuchos::rcp(new ROL::Vector_SimOpt<RealT>(up,zp));

    // Fill initial state
    RealT N(n), one(1), T(1.5);
    for (int i = 0; i < n+1; ++i) {
      RealT I(i);
      (*u_rcp)[i]         = I/N * x0 + (one-I/N) * xn;
      (*u_rcp)[n+1+i]     = I/N * y0 + (one-I/N) * yn;
      (*u_rcp)[3*n+3+i]   = (xn - x0)/T;
      (*u_rcp)[4*n+4+i]   = (yn - y0)/T;
    }
    (*u_rcp)[9*n+9] = T;
    (*z_rcp)[0] = (xn - x0)/T;
    (*z_rcp)[1] = (yn - y0)/T;

    // Fill bounds
    for (int i = 0; i < 9*n+10; ++i) {
      (*ulb_rcp)[i] = ROL::ROL_NINF<RealT>();
      (*uub_rcp)[i] = ROL::ROL_INF<RealT>();
    }
    (*ulb_rcp)[9*n+9] = static_cast<RealT>(0);
    for (int i = 0; i < 2; ++i) {
      (*zlb_rcp)[i] = ROL::ROL_NINF<RealT>();
      (*zub_rcp)[i] = ROL::ROL_INF<RealT>();
    }
    for (int i = 0; i < n+1; ++i) {
      (*ilb_rcp)[i] = static_cast<RealT>(0);
      (*iub_rcp)[i] = Rg*Rg;
    }

    // Initialize bound constraints
    Teuchos::RCP<ROL::Bounds<RealT> > ubnd
      = Teuchos::rcp(new ROL::Bounds<RealT>(ulbp,uubp));
    Teuchos::RCP<ROL::Bounds<RealT> > zbnd
      = Teuchos::rcp(new ROL::Bounds<RealT>(zlbp,zubp));
    zbnd->deactivate();
    Teuchos::RCP<ROL::BoundConstraint<RealT> > xbnd
      = Teuchos::rcp(new ROL::BoundConstraint_SimOpt<RealT>(ubnd,zbnd));
    Teuchos::RCP<ROL::Bounds<RealT> > ibnd
      = Teuchos::rcp(new ROL::Bounds<RealT>(ilbp,iubp));

    // Dynamic constraints
    Teuchos::RCP<PuttingConstraint<RealT> > econ
      = Teuchos::rcp(new PuttingConstraint<RealT>(g,m,x0,y0,xn,yn,mu));

    // Green constraints
    Teuchos::RCP<GreenConstraint<RealT> > icon
      = Teuchos::rcp(new GreenConstraint<RealT>());

    // Final speed objective
    RealT target = sf*sf;
    Teuchos::RCP<PuttingObjective<RealT> > obj
      = Teuchos::rcp(new PuttingObjective<RealT>(target));

    // Initialize optimization problem
    Teuchos::RCP<ROL::OptimizationProblem<RealT> > problem;
    bool useReduced = parlist->sublist("Problem").get("Use reduced space",false);
    if (!useReduced) {
      bool useGreenCon = parlist->sublist("Problem").get("Use green constraint",false);
      if (!useGreenCon) {
        problem = Teuchos::rcp(new ROL::OptimizationProblem<RealT>(obj,x,xbnd,econ,emul));
      }
      else {
        problem = Teuchos::rcp(new ROL::OptimizationProblem<RealT>(obj,x,xbnd,econ,emul,icon,imul,ibnd));
      }
    }
    else {
      econ->setSolveParameters(*parlist);
      Teuchos::RCP<ROL::SimController<RealT> > stateStore
        = Teuchos::rcp(new ROL::SimController<RealT>());
      Teuchos::RCP<ROL::Objective<RealT> > robj
        = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(obj,econ,stateStore,up,zp,emul,true,false));
      problem = Teuchos::rcp(new ROL::OptimizationProblem<RealT>(robj,zp));
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
    *outStream << "Initial x-velocity: " << (*z_rcp)[0]     << std::endl;
    *outStream << "Initial y-velocity: " << (*z_rcp)[1]     << std::endl;
    *outStream << "Final time: "         << (*u_rcp)[9*n+9] << std::endl;

    // Print optimal trajectory
    bool printTrajectory = parlist->sublist("Problem").get("Print trajectory",true);
    if (printTrajectory) {
      std::ofstream xFile, vFile, aFile;
      xFile.open("position.txt");
      vFile.open("velocity.txt");
      aFile.open("acceleration.txt");
      for (int i = 0; i < n+1; ++i) {
        xFile << std::scientific << std::setprecision(8) << std::setw(12) << std::left;
        xFile << (*u_rcp)[i]       << "  " << (*u_rcp)[n+1+i]   << "  " << (*u_rcp)[2*n+2+i] << std::endl;
        vFile << std::scientific << std::setprecision(8) << std::setw(12) << std::left;
        vFile << (*u_rcp)[3*n+3+i] << "  " << (*u_rcp)[4*n+4+i] << "  " << (*u_rcp)[5*n+5+i] << std::endl;
        aFile << std::scientific << std::setprecision(8) << std::setw(12) << std::left;
        aFile << (*u_rcp)[6*n+6+i] << "  " << (*u_rcp)[7*n+7+i] << "  " << (*u_rcp)[8*n+8+i] << std::endl;
      }
      xFile.close();
      vFile.close();
      aFile.close();
    }

/*
    RealT tol(1e-8);
    std::ofstream jacFile;
    jacFile.open("jacobian.txt");
    Teuchos::RCP<ROL::Vector<RealT> > jv = emul->clone();
    for (int i = 0; i < 9*n+10; ++i) {
      jacFile << std::scientific << std::setprecision(8) << std::setw(12) << std::left;
      for (int j = 0; j < 9*n+10; ++j) {
        econ->applyJacobian_1(*jv,*(up->basis(i)),*up,*zp,tol);
        jacFile << jv->dot(*(jv->basis(j))) << "  ";
      }
      jacFile << std::endl;
    }
    jacFile.close();

    Teuchos::RCP<ROL::Vector<RealT> > ajv = up->clone();
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
