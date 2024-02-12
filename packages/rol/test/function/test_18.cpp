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

/*! \file  test_18.cpp
    \brief Test ChainRuleObjective class

    Compose the objective \f$f(y)\f$ where \f$y = g(x)\f$ with 

    \f[
       y_1 &= g_1(x) &= x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 - 10 \\
       y_2 &= g_2(x) &= x_2*x_3 - 5*x_4*x_5                        \\
       y_3 &= g_3(x) &= x_1^3 + x_2^3 + 1
    \f]

    and

    \f[
      f(\mathbf{y}) = \mathbf{y}^\top\mathbf{y} 
                    + \frac{1}{4}(\mathbf{k}^\top \mathbf{y})^2 
                    + \frac{1}{16}(\mathbf{k}^\top \mathbf{y})^4 
    \f]
    Where \f$\mathbf{k}=(1,\cdots,n)\f$
*/

#include "ROL_RandomVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_ChainRuleObjective.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_SimpleEqConstrained.hpp"
#include "ROL_Zakharov.hpp"


int main(int argc, char *argv[]) {

  using RealT        = double;
  using VectorT      = ROL::StdVector<RealT>;
  using ObjectiveT   = ROL::ZOO::Objective_Zakharov<RealT>;
  using ConstraintT  = ROL::ZOO::EqualityConstraint_SimpleEqConstrained<RealT>;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  // Save the format state of the original std::cout.
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(std::cout);

//  RealT errtol = std::sqrt(ROL::ROL_THRESHOLD<RealT>());

  int errorFlag  = 0;

  // *** Test body.

  try {

    uint x_dim = 5; // Constraint domain space dimension
    uint y_dim = 3; // Constraint range space dimension and objective domain space dimension


    // Make a Zakharov objective function f(y)
    auto k_ptr = ROL::makePtr<VectorT>(y_dim);
    auto& k = *k_ptr;
    k[0] = 1;
    k[1] = 2;
    k[2] = 3;
    
    auto x = VectorT(x_dim);
    auto l = VectorT(y_dim);

    auto obj_ptr = ROL::makePtr<ObjectiveT>(k_ptr);
    auto con_ptr = ROL::makePtr<ConstraintT>();
    
    VectorT v(x_dim), g(x_dim), hv(x_dim), u(x_dim);

    auto obj = ROL::ChainRuleObjective<RealT>(obj_ptr,con_ptr,x,l);

    ROL::RandomizeVector(x);
    ROL::RandomizeVector(v);
    ROL::RandomizeVector(u);

    RealT tol = std::sqrt(ROL::ROL_EPSILON<RealT>());

    auto result_1 = obj.checkGradient(x,v,true,*outStream,7,4);

    bool gradient_passed = false;

    for( auto& row : result_1 ) {
      if(row[3] < tol) {
        gradient_passed = true;
        break;
      }
    }

    errorFlag += (!gradient_passed);

    auto result_2 = obj.checkHessVec(x,hv,v,true,*outStream,7,4);

    bool hessVec_passed = false;

    for( auto& row : result_2 ) {
      if(row[3] < tol) {
        hessVec_passed = true;
        break;
      }
    }

    errorFlag += (!hessVec_passed) << 1;

    auto result_3 = obj.checkHessSym(x,hv,v,u,true,*outStream);
    auto hessSym_passed = (result_3[2] < tol);
    
    errorFlag += (!hessSym_passed) << 2;

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

