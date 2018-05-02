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

#include "ROL_FiniteDifference.hpp"

#include "ROL_SimpleEqConstrained.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_RandomVector.hpp"

#include <functional>




int main(int argc, char *argv[]) {

  using namespace ROL;
  using RealT = double;
  using V = Vector<RealT>;
  using std::vector;
  using std::bind;
  using std::reference_wrapper;
  using namespace std::placeholders;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int iprint     = argc - 1;
  Ptr<std::ostream> os;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    os = makePtrFromRef(std::cout);
  else
    os = makePtrFromRef(bhs);

  int errorFlag = 0;

  RealT tol = std::sqrt(ROL_EPSILON<RealT>());

  try {

    int dim = 5;
    int nc = 3;
        
    // Retrieve objective, constraint, iteration vector, solution vector.
    ZOO::getSimpleEqConstrained <RealT> SEC;

    auto obj_ptr = SEC.getObjective();
    auto con_ptr = SEC.getEqualityConstraint();
    auto& obj = *obj_ptr;
    auto& con = *con_ptr;

    auto x_ptr = SEC.getInitialGuess();
    auto l_ptr = SEC.getEqualityMultiplier();
    auto d_ptr = x_ptr->clone();
    auto v_ptr = x_ptr->clone();
    auto g_ptr = x_ptr->dual().clone();
    auto c_ptr = l_ptr->dual().clone();

    RandomizeVector( *d_ptr );
    RandomizeVector( *v_ptr );
    RandomizeVector( *l_ptr );

    reference_wrapper<V> x_ref(*x_ptr);
    reference_wrapper<V> d_ref(*d_ptr);
    reference_wrapper<V> v_ref(*v_ptr);
    reference_wrapper<V> l_ref(*l_ptr);
    reference_wrapper<V> g_ref(*g_ptr);
    reference_wrapper<V> c_ref(*c_ptr);
  
    auto obj_up   = bind( &Objective<RealT>::update, &obj, _1, true, 0 );
    auto obj_val  = bind( &Objective<RealT>::value, &obj, _1, tol );
    auto obj_grad = bind( &Objective<RealT>::gradient, &obj, _1, _2, tol );
    auto obj_hess = bind( &Objective<RealT>::hessVec, &obj, _1, _2, _3, tol );

    auto con_up   = bind( &Constraint<RealT>::update, &con, _1, true, 0 );
    auto con_val  = bind( &Constraint<RealT>::value, &con, _1, _2, tol );
    auto con_jac  = bind( &Constraint<RealT>::applyJacobian, &con, _1, _2, _3, tol );

    // Extra work is needed here because applyAdjointJacobian is overloaded 
//    auto con_ajac = bind( &Constraint<RealT>::applyAdjointJacobian, &con, _1, l_ref, _2, tol );
//    auto con_hess = bind( &Constraint<RealT>::applyAdjointHessian, &con, _1, l_ref, _2, _3, tol );
 
    FiniteDifference<RealT> fd( 1, 13, 20, 11, true, *os );
 
    fd.scalar_check( obj_val,  obj_grad, obj_up, g_ref, v_ref, x_ref, "grad'*dir" );
    fd.vector_check( obj_grad, obj_hess, obj_up, g_ref, v_ref, x_ref, "norm(Hess*vec)" );
 
    fd.vector_check( con_val,  con_jac,  con_up, c_ref, v_ref, x_ref, "norm(Jac*vec)" ); 
//    fd.vector_check( con_ajac, con_hess, con_up, d_ref, v_ref, x_ref, "norm(adj(H)(u,v))" ); 

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

