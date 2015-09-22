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

/** \file
    \brief An example equality constrained problem combining ROL and Sacado 
           This is the same problem as found in rol/examples/simple-eq-constr
           with the objective gradient, objective Hessian direction, constraint 
           Jacobian direction, constraint adjoint Jacobian direction, and
           constraint adjoint Hessian direction computed via automatic 
           differentiation with Sacado.  

    \author Created by G. von Winckel
**/

#include <iostream>


#include "ROL_Sacado_Objective_SimOpt.hpp"
#include "ROL_Sacado_EqualityConstraint_SimOpt.hpp"
#include "ROL_Vector_SimOpt.hpp"

#include "ROL_LineSearchStep.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_CompositeStepSQP.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "example_03.hpp"

using namespace ROL;

typedef double RealT;

int main(int argc, char **argv)
{


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

      Teuchos::ParameterList parlist;
      parlist.set("Nominal SQP Optimality Solver Tolerance", 1.e-2);
      ROL::CompositeStepSQP<RealT> step(parlist);

      int ni = 12;     // Number of interpolation points
      int nq = 20;     // Number of quadrature points

      RealT xl = 0.0;  // Left end of domain
      RealT xr = 5.0;  // Right end of domain

      Teuchos::RCP<NodalBasis<RealT> > basisp = Teuchos::rcp(new NodalBasis<RealT> (ni,nq) );

      // simulation vector (state u) and optimization vector (control z)
      Teuchos::RCP<std::vector<RealT> > u_rcp = Teuchos::rcp(new std::vector<RealT>(ni,0));
      Teuchos::RCP<std::vector<RealT> > z_rcp = Teuchos::rcp(new std::vector<RealT>(ni,0));

      // change of simulation vector (state u) and optimization vector (control z)
      Teuchos::RCP<std::vector<RealT> > yu_rcp = Teuchos::rcp(new std::vector<RealT>(ni,0));
      Teuchos::RCP<std::vector<RealT> > yz_rcp = Teuchos::rcp(new std::vector<RealT>(ni,0));

      // Gradients 
      Teuchos::RCP<std::vector<RealT> > gu_rcp = Teuchos::rcp(new std::vector<RealT>(ni,0));
      Teuchos::RCP<std::vector<RealT> > gz_rcp = Teuchos::rcp(new std::vector<RealT>(ni,0));
      
      // Constraint and Lagrange multiplier value
      Teuchos::RCP<std::vector<RealT> > c_rcp = Teuchos::rcp(new std::vector<RealT>(ni,1.0));
      Teuchos::RCP<std::vector<RealT> > l_rcp = Teuchos::rcp(new std::vector<RealT>(ni,1.0));
      
        

      for(int i=0; i<ni; ++i) {
          (*u_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
          (*z_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
          (*yu_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
          (*yz_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;

      }
      

      StdVector<RealT> u(u_rcp);
      StdVector<RealT> z(z_rcp);
      StdVector<RealT> yu(yu_rcp);
      StdVector<RealT> yz(yz_rcp);
      StdVector<RealT> gu(gu_rcp);
      StdVector<RealT> gz(gz_rcp);
      StdVector<RealT> c(c_rcp);
      StdVector<RealT> l(l_rcp);

      Teuchos::RCP<Vector<RealT> > zp  = Teuchos::rcp(&z,false);
      Teuchos::RCP<Vector<RealT> > up  = Teuchos::rcp(&u,false);
      Teuchos::RCP<Vector<RealT> > yzp  = Teuchos::rcp(&yz,false);
      Teuchos::RCP<Vector<RealT> > yup  = Teuchos::rcp(&yu,false);
      Teuchos::RCP<Vector<RealT> > gzp  = Teuchos::rcp(&gz,false);
      Teuchos::RCP<Vector<RealT> > gup  = Teuchos::rcp(&gu,false);

      // Set tracking target 
      Teuchos::RCP<std::vector<RealT> > u_targ_rcp = Teuchos::rcp(new std::vector<RealT>(ni,0));

      for(int i=0;i<ni;++i) {
          (*u_targ_rcp)[i] = tanh(2.0*((*basisp->xip_)[i]-2.5));            
      }
     
      // Regularization parameter
      RealT gamma = 1e-4; 

      StdVector<RealT> u_targ(u_targ_rcp);

      BoundaryValueProblem<RealT> bvp(xl,xr,basisp); 
      QuadraticTracking<RealT> quadtrack(gamma,u_targ,basisp);

      // Instantiate constraint and objective objects
      BVP_Constraint<RealT,BoundaryValueProblem> constr(bvp);
      Sacado_Objective_SimOpt<RealT,QuadraticTracking> obj(quadtrack); 


      Vector_SimOpt<RealT> x(up,zp);
      Vector_SimOpt<RealT> g(gup,gzp);
      Vector_SimOpt<RealT> y(yup,yzp); 
 
      // Check derivatives
      obj.checkGradient(x,x,y,true,*outStream);
      obj.checkHessVec(x,x,y,true,*outStream);
      constr.checkApplyJacobian(x,y,c,true,*outStream); 
      constr.checkApplyAdjointJacobian(x,yu,c,x,true,*outStream); 
      constr.checkApplyAdjointHessian(x,yu,y,x,true,*outStream);
      try {
        constr.checkInverseJacobian_1(c,yu,u,z,true,*outStream);
      } catch (const std::logic_error &e) {
        *outStream << e.what();
      } 
      try {
        constr.checkInverseAdjointJacobian_1(c,yu,u,z,true,*outStream);
      } catch (const std::logic_error &e) {
        *outStream << e.what();
      } 
      try {
        constr.checkSolve(u,z,c,true,*outStream);
      } catch (const std::logic_error &e) {
        *outStream << e.what();
      } 
 
      // Define Status Test
      RealT gtol  = 1e-12;  // norm of gradient tolerance
      RealT ctol  = 1e-12;  // norm of constraint tolerance
      RealT stol  = 1e-18;  // norm of step tolerance
      int   maxit = 1000;    // maximum number of iterations
      ROL::StatusTestSQP<RealT> status(gtol, ctol, stol, maxit);    

      // Define Algorithm
      ROL::DefaultAlgorithm<RealT> algo(step, status, false);
      algo.run(x,g,l,c,obj,constr,true,*outStream);

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
