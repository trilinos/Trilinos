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

/*! \file  test_13.cpp
    \brief Validate BinaryConstraint class

    Solve a bound constrained problem where a solution is sought on 
    the surface of the bounding set by converting the bound constraint
    into an equality constraint. 

*/

#include "ROL_OptimizationSolver.hpp"

#include "ROL_BinaryConstraint.hpp"
#include "ROL_DiagonalOperator.hpp"
#include "ROL_QuadraticObjective.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Bounds.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"




// Creates f(x) = <x,Dx>/2 - <x,b> where D is a diagonal matrix
// If D_{ii}>0 for all i, then the minimizer is the solution to 
// the linear system x_i=b_i/d_i
template<class Real> 
ROL::Ptr<ROL::Objective<Real>> 
createDiagonalQuadraticObjective( const ROL::Vector<Real> &a, 
                                  const ROL::Vector<Real> &b ) {
  
  auto op = ROL::makePtr<ROL::DiagonalOperator<Real>>(a);
  auto vec = b.clone();
  vec->set(b);
  auto obj = ROL::makePtr<ROL::QuadraticObjective<Real>>(op,vec);
  return obj;
}

typedef double RealT;

int main(int argc, char *argv[]) {
  
//  typedef ROL::Vector<RealT>    V;
  typedef ROL::StdVector<RealT> SV;

   

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a
  // (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try { 
 
    auto parlist = ROL::getParametersFromXmlFile("binary_constraint.xml");

    // Penalty parameter
    RealT gamma = 1.0;

    RealT INF  = ROL::ROL_INF<RealT>();
    RealT NINF = ROL::ROL_NINF<RealT>();
  
    /*-----  Optimization Vector -----*/

    // Initial guess
    auto x0_ptr = ROL::makePtr<std::vector<RealT>>(4);  
    auto x0 = ROL::makePtr<SV>(x0_ptr);
    ROL::RandomizeVector(*x0);

    auto x = x0->clone(); x->set(*x0);

    /*----- Objective Function -----*/

    // Diagonal quadratic objective scaling vector
    auto d_ptr = ROL::makePtr<std::vector<RealT>>(4);  
    auto d = ROL::makePtr<SV>(d_ptr);

    // Quadratic objective offset vector
    auto b_ptr = ROL::makePtr<std::vector<RealT>>(4);  
    auto b = ROL::makePtr<SV>(b_ptr);

    // Set values for objective
    (*b_ptr)[0] = 1.0;  (*b_ptr)[1] = 1.0;  
    (*b_ptr)[2] = 1.0;  (*b_ptr)[3] = 1.0;   

    (*d_ptr)[0] = 1.0;  (*d_ptr)[1] = 2.0;  
    (*d_ptr)[2] = 4.0;  (*d_ptr)[3] = 8.0;

    auto obj = createDiagonalQuadraticObjective( *d, *b ); 

    // Unconstrained minimizer: x = [ 1.0, 0.5, 0.25, 0.125 ]


    /*----- Bound Constraint -----*/

    // Lower bound vector
    auto xl_ptr = ROL::makePtr<std::vector<RealT>>(4);  
    auto xl = ROL::makePtr<SV>(xl_ptr);

    // Upper bound vector
    auto xu_ptr = ROL::makePtr<std::vector<RealT>>(4);  
    auto xu = ROL::makePtr<SV>(xu_ptr);
    
    // Set bounds
    (*xl_ptr)[0] = 0.0;   (*xl_ptr)[1] = 0.0;  
    (*xl_ptr)[2] = NINF;  (*xl_ptr)[3] = NINF;

    (*xu_ptr)[0] = 1.0;   (*xu_ptr)[1] = INF;  
    (*xu_ptr)[2] = 1.0;   (*xu_ptr)[3] = INF;

//    ROL::BoundConstraint<RealT> bnd(xl,xu);
    auto bnd = ROL::makePtr<ROL::Bounds<RealT>>(xl,xu);

    /*---- Constraint and Lagrange Multiplier -----*/

    auto con = ROL::makePtr<ROL::BinaryConstraint<RealT>>( bnd, gamma );

    // Lagrange multiplier
    auto l = x->dual().clone();
    ROL::Elementwise::Fill<RealT> FillWithOnes(1.0);
    l->applyUnary( ROL::Elementwise::Fill<RealT>(1.0) );

    // Constrained minimizer set X = { [ 0, 0, 1, 0.125 ], [ 1, 0, 1, 0.125 ] }
   
    // Create Optimization problems
    ROL::OptimizationProblem<RealT> problem_E( obj, x, ROL::nullPtr, con, l ); 
    ROL::OptimizationProblem<RealT> problem_EB( obj, x, bnd, con, l ); 
   
    // Perform linear algebra and finite difference checks
    problem_E.check( *outStream );
    

    // Solve using Composite Step where the bound is not enforced and
    // equality constraints are satisfied asymptotically
    parlist->sublist("Step").set("Type","Composite Step");
    ROL::OptimizationSolver<RealT> solver_cs( problem_E, *parlist );
    solver_cs.solve( *outStream );    
    *outStream << "\n\nFinal optimization vector:";
    x->print(*outStream);

    // Reset optimization vector and Lagrange multiplier to initial values 
    x->set(*x0); l->applyUnary(FillWithOnes);
  

    // Solve using Augmented Lagrangian where the bound is enforced explicitly 
    // and equality constraints are enforced through penalization
    parlist->sublist("Step").set("Type","Augmented Lagrangian");
    ROL::OptimizationSolver<RealT> solver_al( problem_EB, *parlist );
    solver_al.solve( *outStream );    
    *outStream << "\n\nFinal optimization vector:";
    x->print(*outStream);
    

    // Reset optimization vector and Lagrange multiplier to initial values 
    x->set(*x0); l->applyUnary(FillWithOnes);


    // Solve using Moreau-Yosida where the bound is enforced through penalization 
    // and equality constraints are satisfied asymptotically
    parlist->sublist("Step").set("Type","Moreau-Yosida Penalty");
    ROL::OptimizationSolver<RealT> solver_my( problem_EB, *parlist );
    solver_my.solve( *outStream );    
    *outStream << "\n\nFinal optimization vector:";
    x->print(*outStream);
    

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

