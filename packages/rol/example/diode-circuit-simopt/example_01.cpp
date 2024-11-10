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

#include "example_01.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_ConstraintStatusTest.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_ParameterList.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

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

    // bool use_lambertw   = parlist->get("Use Analytical Solution",true); 
    bool use_scale      = parlist->get("Use Scaling For Epsilon-Active Sets",true);
    bool use_sqp        = parlist->get("Use SQP", true);
    // bool use_linesearch = parlist->get("Use Line Search", true);
    // bool datatype       = parlist->get("Get Data From Input File",false);
    // bool use_adjoint    = parlist->get("Use Adjoint Gradient Computation",false);
    // int  use_hessvec    = parlist->get("Use Hessian-Vector Implementation",1); // 0 - FD, 1 - exact, 2 - GN
    // bool plot           = parlist->get("Generate Plot Data",false);
    // RealT noise         = parlist->get("Measurement Noise",0.0);

    
    Constraint_DiodeCircuit<RealT> con(V_th,lo_Vsrc,up_Vsrc,step_Vsrc);

    RealT alpha = 1.e-4; // regularization parameter (unused)
    int ns = 101; // number of Vsrc components
    int nz = 2; // number of optimization variables

    Objective_DiodeCircuit<RealT> obj(alpha,ns,nz);
    
    // Initialize iteration vectors.
    ROL::Ptr<std::vector<RealT> > z_ptr    = ROL::makePtr<std::vector<RealT>>(nz, 0.0);
    ROL::Ptr<std::vector<RealT> > yz_ptr   = ROL::makePtr<std::vector<RealT>>(nz, 0.0);
    ROL::Ptr<std::vector<RealT> > soln_ptr = ROL::makePtr<std::vector<RealT>>(nz, 0.0);
    (*z_ptr)[0]     = init_Is;
    (*z_ptr)[1]     = init_Rs;
    (*yz_ptr)[0]    = init_Is;
    (*yz_ptr)[1]    = init_Rs;
    (*soln_ptr)[0]  = true_Is;
    (*soln_ptr)[1]  = true_Rs;
    ROL::StdVector<RealT> z(z_ptr);
    ROL::StdVector<RealT> yz(yz_ptr);
    ROL::StdVector<RealT> soln(soln_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zp  = ROL::makePtrFromRef(z);
    ROL::Ptr<ROL::Vector<RealT> > yzp = ROL::makePtrFromRef(yz);

    ROL::Ptr<std::vector<RealT> > u_ptr  = ROL::makePtr<std::vector<RealT>>(ns, 0.0);
    ROL::Ptr<std::vector<RealT> > yu_ptr = ROL::makePtr<std::vector<RealT>>(ns, 0.0);
    std::ifstream input_file("measurements.dat");
    RealT temp, temp_scale;
    for (int i=0; i<ns; i++) {
      input_file >> temp;
      input_file >> temp;
      temp_scale = pow(10,int(log10(temp)));
      (*u_ptr)[i] = temp_scale*(RealT)rand()/(RealT)RAND_MAX;
      (*yu_ptr)[i] = temp_scale*(RealT)rand()/(RealT)RAND_MAX;
    }
    input_file.close();
    ROL::StdVector<RealT> u(u_ptr);
    ROL::StdVector<RealT> yu(yu_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up  = ROL::makePtrFromRef(u);
    ROL::Ptr<ROL::Vector<RealT> > yup = ROL::makePtrFromRef(yu);

    ROL::Ptr<std::vector<RealT> > jv_ptr  = ROL::makePtr<std::vector<RealT>>(ns, 1.0);
    ROL::StdVector<RealT> jv(jv_ptr);
    ROL::Ptr<ROL::Vector<RealT> > jvp = ROL::makePtrFromRef(jv);

    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> y(yup,yzp);

    // Check derivatives
    obj.checkGradient(x,x,y,true,*outStream);
    obj.checkHessVec(x,x,y,true,*outStream);

    con.checkApplyJacobian(x,y,jv,true,*outStream);
    con.checkApplyAdjointJacobian(x,yu,jv,x,true,*outStream);
    con.checkApplyAdjointHessian(x,yu,y,x,true,*outStream);
    // Check consistency of Jacobians and adjoint Jacobians.
    con.checkAdjointConsistencyJacobian_1(jv,yu,u,z,true,*outStream);
    con.checkAdjointConsistencyJacobian_2(jv,yz,u,z,true,*outStream);
    // Check consistency of solves.
    con.checkSolve(u,z,jv,true,*outStream);
    con.checkInverseJacobian_1(jv,yu,u,z,true,*outStream);
    con.checkInverseAdjointJacobian_1(yu,jv,u,z,true,*outStream);

    // Initialize reduced objective function.
    ROL::Ptr<std::vector<RealT> > p_ptr  = ROL::makePtr<std::vector<RealT>>(ns, 0.0);
    ROL::StdVector<RealT> p(p_ptr);
    ROL::Ptr<ROL::Vector<RealT> > pp  = ROL::makePtrFromRef(p);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > pobj  = ROL::makePtrFromRef(obj);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pcon = ROL::makePtrFromRef(con);
    ROL::Reduced_Objective_SimOpt<RealT> robj(pobj,pcon,up,zp,pp);
    // Check derivatives.
    *outStream << "Derivatives of reduced objective" << std::endl;
    robj.checkGradient(z,z,yz,true,*outStream);
    robj.checkHessVec(z,z,yz,true,*outStream);
    
    // Bound constraints
    RealT tol = 1.e-12;
    ROL::Ptr<std::vector<RealT> > g0_ptr = ROL::makePtr<std::vector<RealT>>(nz, 0.0);;
    ROL::StdVector<RealT> g0p(g0_ptr);
    robj.gradient(g0p,z,tol);
    *outStream << std::scientific <<  "Initial gradient = " << (*g0_ptr)[0] << " " << (*g0_ptr)[1] << "\n";
    *outStream << std::scientific << "Norm of Gradient = " << g0p.norm() << "\n";

    // Define scaling for epsilon-active sets (used in inequality constraints)
    RealT scale;
    if(use_scale){ scale = 1.0e-2/g0p.norm();}
    else{ scale = 1.0;}
    *outStream << std::scientific << "Scaling: " << scale << "\n";

    /// Define constraints on Is and Rs
    BoundConstraint_DiodeCircuit<RealT> bcon(scale,lo_Is,up_Is,lo_Rs,up_Rs);
    //bcon.deactivate();
    
    // Optimization 
    *outStream << "\n Initial guess " << (*z_ptr)[0] << " " << (*z_ptr)[1] << std::endl;
      
    if (!use_sqp){    
      // Trust Region
      ROL::Ptr<ROL::Step<RealT>>
        step = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
      ROL::Ptr<ROL::StatusTest<RealT>>
        status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
      ROL::Algorithm<RealT> algo_tr(step,status,false);
      std::clock_t timer_tr = std::clock();
      algo_tr.run(z,robj,bcon,true,*outStream);
      *outStream << "\n Solution " << (*z_ptr)[0] << " " << (*z_ptr)[1] << "\n" << std::endl;
      *outStream << "Trust-Region required " << (std::clock()-timer_tr)/(RealT)CLOCKS_PER_SEC
                 << " seconds.\n";
    }
    else{
      // SQP.
      //ROL::Ptr<std::vector<RealT> > gz_ptr = ROL::makePtr<std::vector<RealT>>(nz, 0.0);
      //ROL::StdVector<RealT> gz(gz_ptr);
      //ROL::Ptr<ROL::Vector<RealT> > gzp = &gz,false;
      ROL::Ptr<std::vector<RealT> > gu_ptr = ROL::makePtr<std::vector<RealT>>(ns, 0.0);
      ROL::StdVector<RealT> gu(gu_ptr);
      ROL::Ptr<ROL::Vector<RealT> > gup = ROL::makePtrFromRef(gu);
      //ROL::Vector_SimOpt<RealT> g(gup,gzp);
      ROL::Vector_SimOpt<RealT> g(gup,zp);
      ROL::Ptr<std::vector<RealT> > c_ptr = ROL::makePtr<std::vector<RealT>>(ns, 0.0);
      ROL::Ptr<std::vector<RealT> > l_ptr = ROL::makePtr<std::vector<RealT>>(ns, 0.0);
      ROL::StdVector<RealT> c(c_ptr);
      ROL::StdVector<RealT> l(l_ptr);
      
      ROL::Ptr<ROL::Step<RealT>>
        step = ROL::makePtr<ROL::CompositeStep<RealT>>(*parlist);
      ROL::Ptr<ROL::StatusTest<RealT>>
        status = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(*parlist);
      ROL::Algorithm<RealT> algo_cs(step,status,false);
      //x.zero();
      std::clock_t timer_cs = std::clock();
      algo_cs.run(x,g,l,c,obj,con,true,*outStream);
      *outStream << "\n Solution " << (*z_ptr)[0] << " " << (*z_ptr)[1] << "\n" << std::endl;
      *outStream << "Composite Step required " << (std::clock()-timer_cs)/(RealT)CLOCKS_PER_SEC
		 << " seconds.\n";
    }
    soln.axpy(-1.0, z);
    *outStream << "Norm of error: " << soln.norm() << std::endl;
    if (soln.norm() > 1e4*ROL::ROL_EPSILON<RealT>()) {
      errorFlag = 1;
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

