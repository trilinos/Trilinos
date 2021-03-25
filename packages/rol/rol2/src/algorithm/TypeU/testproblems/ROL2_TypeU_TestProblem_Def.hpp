#pragma once
#ifndef ROL2_TYPEU_TESTPROBLEM_DEF_HPP
#define ROL2_TYPEU_TESTPROBLEM_DEF_HPP

namespace ROL2 {
namespace TypeU {

template<class Real>
void TestProblem<Real>::get( Ptr<OptimizationProblem<Real>>& problem,
                             Ptr<Vector<Real>>&              x0,
                             std::vector<Ptr<Vector<Real>>>& xsols ) const {
  if( x0.is_nullPtr() ) x0 = getInitialGuess();
  x0->set( *getInitialGuess() );
  
  xsols.resize(getNumSolutions());
  for( int i=0; i<getNumSolutions(); ++i ) {
    if( x[i].is_nullPtr() ) x[i] = getSolution(i)->clone();
    x[i]->set(*getSolution(i));
  }
  
  problem = makePtr<OptimizationProblem<Real>>( getObjective(), x0 );
                          
}

} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_TESTPROBLEM_DEF_HPP

