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

#include "ROL_Sacado_Objective.hpp"
#include "ROL_Sacado_Constraint.hpp"

#include "ROL_Algorithm.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_ConstraintStatusTest.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_ParameterList.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "example_02.hpp"

using namespace ROL;

typedef double RealT;

int main(int argc, char **argv)
{


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

    // Run derivative checks, etc.
    int dim = 5;
    int nc = 3;

    ROL::Ptr< Sacado_Objective<RealT,Example_Objective> > obj = 
        ROL::makePtr<Sacado_Objective<RealT,Example_Objective>>();

    ROL::Ptr< Sacado_Constraint<RealT,Example_Constraint > > constr =
        ROL::makePtr<Sacado_Constraint<RealT,Example_Constraint >>(nc);

    ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);

    ROL::Ptr<std::vector<RealT> > sol_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::StdVector<RealT> x(x_ptr);      // Iteration vector.
    ROL::StdVector<RealT> sol(sol_ptr);  // Reference solution vector.

    // Get initial guess
    (*x_ptr)[0] = -1.8;
    (*x_ptr)[1] = 1.7;
    (*x_ptr)[2] = 1.9;
    (*x_ptr)[3] = -0.8;
    (*x_ptr)[4] = -0.8;

    // Get solution
    (*sol_ptr)[0] = -1.717143570394391e+00;
    (*sol_ptr)[1] =  1.595709690183565e+00;
    (*sol_ptr)[2] =  1.827245752927178e+00;
    (*sol_ptr)[3] = -7.636430781841294e-01;
    (*sol_ptr)[4] = -7.636430781841294e-01;

    RealT left = -1e0, right = 1e0;
    ROL::Ptr<std::vector<RealT> > xtest_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > g_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > d_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > v_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > vc_ptr = ROL::makePtr<std::vector<RealT>>(nc, 0.0);
    ROL::Ptr<std::vector<RealT> > vl_ptr = ROL::makePtr<std::vector<RealT>>(nc, 0.0);
    ROL::StdVector<RealT> xtest(xtest_ptr);
    ROL::StdVector<RealT> g(g_ptr);
    ROL::StdVector<RealT> d(d_ptr);
    ROL::StdVector<RealT> v(v_ptr);
    ROL::StdVector<RealT> vc(vc_ptr);
    ROL::StdVector<RealT> vl(vl_ptr);

    // set xtest, d, v
    for (int i=0; i<dim; i++) {
      (*xtest_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*d_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*v_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }
    // set vc, vl
    for (int i=0; i<nc; i++) {
      (*vc_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*vl_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }

    obj->checkGradient(xtest, d, true, *outStream);                             *outStream << "\n"; 
    obj->checkHessVec(xtest, v, true, *outStream);                              *outStream << "\n";
    obj->checkHessSym(xtest, d, v, true, *outStream);                           *outStream << "\n";
    constr->checkApplyJacobian(xtest, v, vc, true, *outStream);                 *outStream << "\n";
    constr->checkApplyAdjointJacobian(xtest, vl, vc, xtest, true, *outStream);  *outStream << "\n";
    constr->checkApplyAdjointHessian(xtest, vl, d, xtest, true, *outStream);    *outStream << "\n";

    ROL::Ptr<std::vector<RealT> > v1_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > v2_ptr = ROL::makePtr<std::vector<RealT>>(nc, 0.0);
    ROL::StdVector<RealT> v1(v1_ptr);
    ROL::StdVector<RealT> v2(v2_ptr);
    RealT augtol = 1e-8;
    constr->solveAugmentedSystem(v1, v2, d, vc, xtest, augtol);
    
    // Define algorithm.
    std::string paramfile = "parameters.xml";
    auto parlist = ROL::getParametersFromXmlFile(paramfile);
    ROL::Algorithm<RealT> algo("Composite Step", *parlist);

    // Run algorithm.
    vl.zero();
    algo.run(x, g, vl, vc, *obj, *constr, true, *outStream);

    // Compute Error
    *outStream << "\nReference solution x_r =\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[0] << "\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[1] << "\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[2] << "\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[3] << "\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[4] << "\n";
    *outStream << "\nOptimal solution x =\n";
    *outStream << std::scientific << "  " << (*x_ptr)[0] << "\n";
    *outStream << std::scientific << "  " << (*x_ptr)[1] << "\n";
    *outStream << std::scientific << "  " << (*x_ptr)[2] << "\n";
    *outStream << std::scientific << "  " << (*x_ptr)[3] << "\n";
    *outStream << std::scientific << "  " << (*x_ptr)[4] << "\n";
    x.axpy(-1.0, sol);
    RealT abserr = x.norm();
    RealT relerr = abserr/sol.norm();
    *outStream << std::scientific << "\n   Absolute Error: " << abserr;
    *outStream << std::scientific << "\n   Relative Error: " << relerr << "\n";
    if ( relerr > sqrt(ROL::ROL_EPSILON<RealT>()) ) {
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
