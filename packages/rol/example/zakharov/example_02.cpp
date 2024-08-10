// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_02.cpp
    \brief Test the OptimizationProblem interface on Step by taking
           alternating Linesearch and Trust-Region steps on the 
           Zahkarov function
*/

#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Zakharov.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  using namespace Teuchos;

  typedef std::vector<RealT>          vector;  
  typedef ROL::Vector<RealT>          V;      // Abstract vector
  typedef ROL::StdVector<RealT>       SV;     // Concrete vector containing std::vector data

  GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  auto outStream = ROL::makeStreamPtr( std::cout, argc > 1 );

  int errorFlag  = 0;

  // *** Example body.
 
  try {

    int dim = 10; // Set problem dimension. 

    std::string paramfile = "parameters.xml";
    auto parlist = ROL::getParametersFromXmlFile( paramfile );

    ROL::Ptr<vector> x_ptr = ROL::makePtr<vector>(dim, 1.0);
    ROL::Ptr<vector> k_ptr = ROL::makePtr<vector>(dim, 0.0);

    ROL::Ptr<V> x = ROL::makePtr<SV>(x_ptr);  // Optimization vector
    ROL::Ptr<V> k = ROL::makePtr<SV>(k_ptr);  // Vector appearing in Zakharov objective

    ROL::Ptr<V> s = x->clone();            // Step vector

    for( int i=0; i<dim; ++i ) {
      (*k_ptr)[i] = i+1.0;
    }
    
    ROL::Ptr<ROL::Objective<RealT>> obj = ROL::makePtr<ROL::ZOO::Objective_Zakharov<RealT>>(k);
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd = ROL::makePtr<ROL::BoundConstraint<RealT>>();
    bnd->deactivate();
    
    ROL::OptimizationProblem<RealT> opt(obj,x);
    ROL::AlgorithmState<RealT> state;
   
    // Allocate iterate vector in algorithm state 
    state.iterateVec = x->clone();
    state.iterateVec->set(*x);
    state.minIterVec = x->clone();
 
    ROL::LineSearchStep<RealT>  ls(*parlist);
    ROL::TrustRegionStep<RealT> tr(*parlist);

    ls.initialize( *opt.getSolutionVector(),
                    opt.getSolutionVector()->dual(),
                   *opt.getObjective(),
                   *bnd, state );
    tr.initialize( *opt.getSolutionVector(),
                    opt.getSolutionVector()->dual(),
                   *opt.getObjective(),
                   *bnd, state );

    for( int iter = 0; iter<10; ++iter ) {
      ls.compute( *s,
                  *opt.getSolutionVector(), 
                  *opt.getObjective(),
                  *bnd, state ); 
      ls.update( *opt.getSolutionVector(),
                 *s,
                 *opt.getObjective(),
                 *bnd, state );

      state.minIterVec->set(*x);
      state.minIter = state.iter;
      state.minValue = state.value;      

      *outStream << "LS fval = " << state.minValue << std::endl; 

      tr.compute( *s,
                  *opt.getSolutionVector(), 
                  *opt.getObjective(),
                  *bnd, state ); 
      tr.update( *opt.getSolutionVector(),
                 *s,
                 *opt.getObjective(),
                 *bnd, state );

      state.minIterVec->set(*x);
      state.minIter = state.iter;
      state.minValue = state.value;      

      *outStream << "TR fval = " << state.minValue << std::endl; 
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



