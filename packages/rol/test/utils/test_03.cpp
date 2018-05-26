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

/*! \file  test_03.cpp
    \brief Test of FiniteDifference class for performing derivative checks
*/

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_ValidateFunction.hpp"

#include "ROL_SimpleEqConstrained.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_Objective_CheckInterface.hpp"
#include "ROL_Constraint_CheckInterface.hpp"


int main(int argc, char *argv[]) {

  using namespace ROL;
  using RealT = double;
  using std::vector;
  using std::bind;
  using std::mem_fn;
  using std::cref;
  using namespace std::placeholders;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int iprint     = argc - 1;
  Ptr<std::ostream> os;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    os = makePtrFromRef(std::cout);
  else
    os = makePtrFromRef(bhs);

  int errorFlag = 0;

  try {

    // Retrieve objective, constraint, iteration vector, solution vector.
    ZOO::getSimpleEqConstrained <RealT> SEC;

    auto obj_ptr = SEC.getObjective();
    auto con_ptr = SEC.getEqualityConstraint();
    auto& obj = *obj_ptr;
    auto& con = *con_ptr;

    auto x = SEC.getInitialGuess();
    auto l = SEC.getEqualityMultiplier();
    auto d = x->clone();
    auto v = x->clone();
    auto g = x->dual().clone();
    auto c = l->dual().clone();

    RandomizeVector( *d );
    RandomizeVector( *v );
    RandomizeVector( *l );

    auto obj_check = make_check(obj);
    auto obj_up    = obj_check.update(); 
    auto obj_val   = obj_check.value();    
    auto obj_grad  = obj_check.gradient(); 
    auto obj_hess  = obj_check.hessVec();  

    auto con_check = make_check(con);
    auto con_up    = con_check.update();
    auto con_val   = con_check.value();
    auto con_jac   = con_check.jacobian();
    auto con_ajac  = con_check.adjointJacobian();
    auto con_ajacl = fix_direction(con_ajac,*l);
    auto con_hess  = con_check.adjointHessian( cref(*l) );

    ValidateFunction<RealT> validator( 1, 13, 20, 11, true, *os );

    // Validate Objective Function

    validator.derivative_check( obj_val,  obj_grad, obj_up, *g, *v, *x, "grad'*dir" );
    validator.derivative_check( obj_grad, obj_hess, obj_up, *g, *v, *x, "norm(Hess*vec)" );
    validator.symmetry_check( obj_hess, obj_up, *d, *v, *x, "objective Hessian", "H" );

    validator.derivative_check( con_val,  con_jac,  con_up, *c, *v, *x, "norm(Jac*vec)" ); 
    validator.derivative_check( con_ajacl, con_hess, con_up, *d, *v, *x, "norm(adj(H)(u,v))" ); 
    validator.adjoint_consistency_check( con_jac, con_ajac, con_up, *v, *l, *x, "Jacobian", "J");


  }
  catch (std::logic_error err) {
    *os << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

