// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief 
*/

#include "ROL_DiodeCircuit.hpp"
#include "ROL_StepFactory.hpp"
#include "ROL_StatusTestFactory.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_ParameterList.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {
    
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );

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
    bool datatype       = parlist->get("Get Data From Input File",false);
    bool use_adjoint    = parlist->get("Use Adjoint Gradient Computation",false);
    int  use_hessvec    = parlist->get("Use Hessian-Vector Implementation",1); // 0 - FD, 1 - exact, 2 - GN
    bool plot           = parlist->get("Generate Plot Data",false);
    RealT noise         = parlist->get("Measurement Noise",0.0);

    ROL::Ptr< ROL::ZOO::Objective_DiodeCircuit<RealT> > obj;
            
    if(datatype){
      // Get objective with data from file
      std::ifstream input_file("diode_forTimur.cir.dat");
      obj = ROL::makePtr<ROL::ZOO::Objective_DiodeCircuit<RealT>>(V_th,input_file,use_lambertw,noise,use_adjoint,use_hessvec);
    }
    else{
      // Generate data and get objective
      obj = ROL::makePtr<ROL::ZOO::Objective_DiodeCircuit<RealT>>(V_th,lo_Vsrc,up_Vsrc,step_Vsrc,true_Is,true_Rs,use_lambertw,noise,use_adjoint,use_hessvec);
    }

    
    // Generate data for plot
    if(plot){
      (*obj).generate_plot(std::min(true_Is,init_Is),std::max(true_Is,init_Is),std::abs((true_Is-init_Is)/100),
			   std::min(true_Rs,init_Rs),std::max(true_Rs,init_Rs),std::abs((true_Rs-init_Rs)/100));
    }
    
    int dim = 2;

    // Define Algorithm
    std::string stepname(parlist->get("Step Type", "Line Search"));
    ROL::StepFactory<RealT> stepFactory;
    ROL::Ptr<ROL::Step<RealT>> step = stepFactory.getStep(stepname,*parlist);
    ROL::StatusTestFactory<RealT> statusFactory;
    ROL::Ptr<ROL::StatusTest<RealT>> status = statusFactory.getStatusTest(stepname,*parlist);
    ROL::Algorithm<RealT> algo(step,status,false);

    // Iteration Vector
    ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    // Set Initial Guess
    (*x_ptr)[0] = init_Is; /// Is
    (*x_ptr)[1] = init_Rs; /// Rs
    // Scaling Vector
    ROL::Ptr<std::vector<RealT> > scaling_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    (*scaling_ptr)[0] = 1e24; /// Is
    (*scaling_ptr)[1] = 1e01; /// Rs
    ROL::PrimalScaledStdVector<RealT> x(x_ptr,scaling_ptr);

    RealT tol = 1.e-12;
    ROL::Ptr<std::vector<RealT> > g0_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::DualScaledStdVector<RealT> g0p(g0_ptr,scaling_ptr);
    (*obj).gradient(g0p,x,tol);
    *outStream << std::scientific <<  "Initial gradient = " << (*g0_ptr)[0] << " " << (*g0_ptr)[1] << "\n";
    *outStream << std::scientific << "Norm of Gradient = " << g0p.norm() << "\n";

    // Define scaling for epsilon-active sets (used in inequality constraints)
    RealT scale;
    if(use_scale){ scale = 1.0e-2/g0p.norm();}
    else{ scale = 1.0;}
     *outStream << std::scientific << "Scaling: " << scale << "\n";

    /// Define constraints on Is and Rs.
    // Bound vectors.
    ROL::Ptr<std::vector<RealT> > IsRs_lower_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    (*IsRs_lower_ptr)[0] = lo_Is; /// Is lower bound
    (*IsRs_lower_ptr)[1] = lo_Rs; /// Rs lower bound
    ROL::Ptr<std::vector<RealT> > IsRs_upper_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    (*IsRs_upper_ptr)[0] = up_Is; /// Is upper bound
    (*IsRs_upper_ptr)[1] = up_Rs; /// Rs upper bound
    ROL::Ptr<ROL::PrimalScaledStdVector<RealT> > lo_IsRs = ROL::makePtr<ROL::PrimalScaledStdVector<RealT>>(IsRs_lower_ptr, scaling_ptr);
    ROL::Ptr<ROL::PrimalScaledStdVector<RealT> > up_IsRs = ROL::makePtr<ROL::PrimalScaledStdVector<RealT>>(IsRs_upper_ptr, scaling_ptr);
    // Bound constraint.
    ROL::Bounds<RealT> con2(lo_IsRs, up_IsRs, scale);

    // Gradient and Hessian check
    // direction for gradient check
    ROL::Ptr<std::vector<RealT> > d_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    RealT left = 0.0, right = 1.0;
    RealT Is_scale = pow(10,int(log10(init_Is)));
    RealT Rs_scale = pow(10,int(log10(init_Rs)));
    (*d_ptr)[0] = Is_scale*(( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left);
    (*d_ptr)[1] = Rs_scale*(( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left);
    ROL::PrimalScaledStdVector<RealT> d(d_ptr,scaling_ptr);
    // check gradient and Hessian-vector computation using finite differences
    (*obj).checkGradient(x, g0p, d, true, *outStream);
    (*obj).checkHessVec(x, g0p, d, true, *outStream);

    // Run Algorithm
    algo.run(x, *obj, con2, true, *outStream);
    
    ROL::Ptr<std::vector<RealT> > gf_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::DualScaledStdVector<RealT> gfp(gf_ptr,scaling_ptr);
    (*obj).gradient(gfp,x,tol);
     *outStream << std::scientific << "Final gradient = " << (*gf_ptr)[0] << " " << (*gf_ptr)[1] << "\n";
     *outStream << std::scientific << "Norm of Gradient = " << gfp.norm() << "\n";
    
    // Get True Solution
    ROL::Ptr<std::vector<RealT> > xtrue_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    (*xtrue_ptr)[0] = true_Is;
    (*xtrue_ptr)[1] = true_Rs;
    ROL::PrimalScaledStdVector<RealT> xtrue(xtrue_ptr,scaling_ptr);
    
    // Print
    *outStream << "Solution:" << "\n";
    *outStream << "x[0] = " <<  std::scientific  << (*x_ptr)[0] <<  "\n";
    *outStream << "x[1] = " <<  std::scientific  << (*x_ptr)[1] <<  "\n";
    
    // Compute Error
    x.axpy(-1.0, xtrue);
    RealT abserr = x.norm();
    RealT relerr = abserr/xtrue.norm();
    *outStream << std::scientific << "\n   Absolute Error: " << abserr;
    *outStream << std::scientific << "\n   Relative Error: " << relerr << "\n";
    if ( relerr > sqrt(ROL::ROL_EPSILON<RealT>()) ) {
      errorFlag += 1;
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

