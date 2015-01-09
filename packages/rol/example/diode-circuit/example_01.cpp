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
    \brief 
*/

#include "ROL_DiodeCircuit.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_Algorithm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "ROL_BoundConstraint.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "ROL_TrustRegionStep.hpp"

#include <iostream>
#include <fstream>
#include <math.h>

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
    
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, Teuchos::Ptr<Teuchos::ParameterList>(&*parlist) );

    RealT V_th      = parlist->get("Thermal Voltage", 0.02585);
    RealT lo_Vsrc   = parlist->get("Source Voltage Lower Bound", 0.0);
    RealT up_Vsrc   = parlist->get("Source Voltage Upper Bound", 1.0);
    RealT step_Vsrc = parlist->get("Source Voltage Step", 1.e-2);

    RealT true_Is = parlist->get("True Saturation Current", 1.e-12);
    RealT true_Rs = parlist->get("True Saturation Resistance", 0.25);
    RealT init_Is = parlist->get("Initial Saturation Current", 1.e-12);
    RealT init_Rs = parlist->get("Initial Saturation Resistance", 0.25);
    RealT lo_Is   = parlist->get("Saturation Current Lower Bound", 1.e-16);
    RealT up_Is   = parlist->get("Saturation Current Upper Bound", 1.e-1);
    RealT lo_Rs   = parlist->get("Saturation Resistance Lower Bound", 1.e-2);
    RealT up_Rs   = parlist->get("Saturation Resistance Upper Bound", 30.0);

    bool use_lambertw   = parlist->get("Use Analytical Solution",true); 
    bool use_scale      = parlist->get("Use Scaling For Epsilon-Active Sets",true);
    bool use_linesearch = parlist->get("Use Line Search", true);
    bool datatype       = parlist->get("Get Data From Input File",false);
    bool use_adjoint    = parlist->get("Use Adjoint Gradient Computation",false);
    int  use_hessvec    = parlist->get("Use Hessian-Vector Implementation",1); // 0 - FD, 1 - exact, 2 - GN
    bool plot           = parlist->get("Generate Plot Data",false);
    RealT noise         = parlist->get("Measurement Noise",0.0);

    Teuchos::RCP< ROL::ZOO::Objective_DiodeCircuit<RealT> > obj;
            
    if(datatype){
      // Get objective with data from file
      std::ifstream input_file("diode_forTimur.cir.dat");
      obj = Teuchos::rcp( new ROL::ZOO::Objective_DiodeCircuit<RealT> (V_th,input_file,use_lambertw,noise,use_adjoint,use_hessvec) );
    }
    else{
      // Generate data and get objective
      obj = Teuchos::rcp( new ROL::ZOO::Objective_DiodeCircuit<RealT> (V_th,lo_Vsrc,up_Vsrc,step_Vsrc,true_Is,true_Rs,use_lambertw,noise,use_adjoint,use_hessvec) );
    }

    
    // Generate data for plot
    if(plot){
      (*obj).generate_plot(std::min(true_Is,init_Is),std::max(true_Is,init_Is),fabs((true_Is-init_Is)/100),
			   std::min(true_Rs,init_Rs),std::max(true_Rs,init_Rs),fabs((true_Rs-init_Rs)/100));
    }
    
    //ROL::getTestObjectives<RealT>(obj,con,x0,z,prob);
    
    int dim = 2; // Set problem dimension. Must be even.

    // Define Step
    Teuchos::RCP< ROL::Step<RealT> > step;
    if(use_linesearch){
      step = Teuchos::rcp( new ROL::LineSearchStep<RealT> (*parlist) );
    }
    else{
      step = Teuchos::rcp( new ROL::TrustRegionStep<RealT> (*parlist) );
    }
    
    // Define Status Test                                                      
    RealT gtol = parlist->get("Gradient Tolerance",1.e-6);
    RealT stol = parlist->get("Step Tolerance",1.e-12);
    int maxit  = parlist->get("Maximum Number of Iterations",100);
    ROL::StatusTest<RealT> status(gtol,stol,maxit);

    // Define Algorithm
    ROL::DefaultAlgorithm<RealT> algo(*step,status,false);

    // Iteration Vector
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    // Set Initial Guess
    (*x_rcp)[0] = init_Is; /// Is
    (*x_rcp)[1] = init_Rs;    /// Rs
    ROL::StdVector<RealT> x(x_rcp);

    RealT tol = 1.e-12;
    Teuchos::RCP<std::vector<RealT> > g0_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );;
    ROL::StdVector<RealT> g0p(g0_rcp);
    (*obj).gradient(g0p,x,tol);
    *outStream << std::scientific <<  "Initial gradient = " << (*g0_rcp)[0] << " " << (*g0_rcp)[1] << "\n";
    *outStream << std::scientific << "Norm of Gradient = " << g0p.norm() << "\n";

    // Define scaling for epsilon-active sets (used in inequality constraints)
    RealT scale;
    if(use_scale){ scale = 1.0e-2/g0p.norm();}
    else{ scale = 1.0;}
     *outStream << std::scientific << "Scaling: " << scale << "\n";

    /// Define constraints on Is and Rs
    ROL::ZOO::BoundConstraint_DiodeCircuit<RealT> con(scale,lo_Is,up_Is,lo_Rs,up_Rs);

    /*--------------------------------------------------------------------------------------------
    // Gradient and Hessian check
    // direction for gradient check
    Teuchos::RCP<std::vector<RealT> > d_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    RealT left = 0.0, right = 1.0;
    RealT Is_scale = pow(10,int(log10(init_Is)));
    RealT Rs_scale = pow(10,int(log10(init_Rs)));
    (*d_rcp)[0] = Is_scale*(( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left);
    (*d_rcp)[1] = Rs_scale*(( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left);
    ROL::StdVector<RealT> d(d_rcp);
    // check gradient and Hessian-vector computation using finite differences
    (*obj).checkGradient(x, d, true);
    (*obj).checkHessVec(x, d, true);
    ---------------------------------------------------------------------------------------------*/

    // Run Algorithm
    std::vector<std::string> output = algo.run(x, *obj, con, false);
    for ( unsigned i = 0; i < output.size(); i++ ) {
      std::cout << output[i];
    }
    
    Teuchos::RCP<std::vector<RealT> > gf_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    ROL::StdVector<RealT> gfp(gf_rcp);
    (*obj).gradient(gfp,x,tol);
     *outStream << std::scientific << "Final gradient = " << (*gf_rcp)[0] << " " << (*gf_rcp)[1] << "\n";
     *outStream << std::scientific << "Norm of Gradient = " << gfp.norm() << "\n";
    
    // Get True Solution
    Teuchos::RCP<std::vector<RealT> > xtrue_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    (*xtrue_rcp)[0] = true_Is;
    (*xtrue_rcp)[1] = true_Rs;
    ROL::StdVector<RealT> xtrue(xtrue_rcp);
    
    // Print
    *outStream << "Solution:" << "\n";
    *outStream << "x[0] = " <<  std::scientific  << (*x_rcp)[0] <<  "\n";
    *outStream << "x[1] = " <<  std::scientific  << (*x_rcp)[1] <<  "\n";
    
    // Compute Error
    x.axpy(-1.0, xtrue);
    RealT abserr = x.norm();
    RealT relerr = abserr/xtrue.norm();
    *outStream << std::scientific << "\n   Absolute Error: " << abserr;
    *outStream << std::scientific << "\n   Relative Error: " << relerr << "\n";
    if ( relerr > sqrt(ROL::ROL_EPSILON) ) {
      errorFlag += 1;
    }
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

