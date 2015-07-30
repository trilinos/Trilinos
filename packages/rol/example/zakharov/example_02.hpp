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

#ifndef ROL_ZAKHAROV_EXAMPLE_02_H
#define ROL_ZAKHAROV_EXAMPLE_02_H


#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <vector>


#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StdVector.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

template<class Real>
class UnconstrainedTest {

  typedef ROL::AlgorithmState<Real>   State;
  typedef ROL::DefaultAlgorithm<Real> Algorithm;
  typedef ROL::LineSearchStep<Real>   Step;
  typedef ROL::Objective<Real>        Obj;
  typedef ROL::Vector<Real>           V;
  typedef ROL::StatusTest<Real>       Status;
  
  typedef std::string                 String;
  typedef std::vector<String>         StrList;  
  typedef std::map<String,StrList>    OP;
  typedef StrList::iterator           Iter;    

  typedef Teuchos::ParameterList      ParList;
  typedef Teuchos::RCP<Obj>           PObj;
  typedef Teuchos::RCP<V>             PV;
  typedef Teuchos::RCP<const State>   PState;


  public:

    UnconstrainedTest( PObj obj, PV x0, const ParList &parlist, 
                       const int maxit ) : 
      obj_(obj), x0_(x0), x_(x0->clone()), parlist_(parlist), 
      maxit_(maxit), status_(1e-12,1e-14,maxit) {

      // Various descent types
      StrList descent{"Steepest Descent","Nonlinear CG", "Quasi-Newton Method", 
                      "Newton's Method","Newton-Krylov"};

      // Various NCG types 
      StrList ncg{"Hestenes-Stiefel","Fletcher-Reeves","Daniel (uses Hessian)", 
                  "Polak-Ribiere","Fletcher Conjugate Descent","Liu-Storey",
                  "Dai-Yuan","Hagar-Zhang","Oren-Luenberger"};

      // Various Line Search types
      StrList ls{"Iteration Scaling","Path-Based Target Level","Backtracking",
                 "Bisection","Golden Section","Cubic Interpolation","Brents"};

      // Various curvature condition types
      StrList curve{"Wolfe Conditions","Strong Wolfe Conditions",
                    "Generalized Wolfe Conditions","Goldstein Conditions",
                    "Approximate Wolfe Conditions","Null Curvature Condition"};
 
      options_.insert( OP::value_type("Descent Type",descent) );
      options_.insert( OP::value_type("Nonlinear CG Type",ncg) );
      options_.insert( OP::value_type("Linesearch Type",ls) );
      options_.insert( OP::value_type("Linesearch Curvature Condition",curve) );


    }

    void run(const String &opt1, const String &opt2) {
      
      StrList items1( options_[opt1] );
      StrList items2( options_[opt2] );

      for(Iter it1 = items1.begin(); it1 != items1.end(); ++it1) {
        for(Iter it2 = items2.begin(); it2 != items2.end(); ++it2) {

          // Reset optimization vector to initial value
          x_->set(*x0_);      

          // Update the parameters to use
          parlist_.set(opt1,*it1);
          parlist_.set(opt2,*it2);

          // Create the step
          Step step(parlist_);

          // Create the algorithm using step and status
          Algorithm algo(step,status_,false);

          // Run the algorithm
          StrList output = algo.run(*x_,*obj_,false);

          PState state = algo.getState();         

          std::cout << *it1 << ", " << *it2 << ": " << state->gnorm << std::endl;            
        }  
      } 

    }

  private:

    PObj    obj_;      // Pointer to an objective  
    PV      x0_;       // Initial optimization vector
    PV      x_;        // Optimization vector     
    ParList parlist_;  // Parameters 
    int     maxit_;    // Maximum number of iterations
    Status  status_;   // Status test object
    OP      options_;  // Map of possible parameters to iterate over
    
};


#endif // ROL_ZAKHAROV_EXAMPLE_02_H
