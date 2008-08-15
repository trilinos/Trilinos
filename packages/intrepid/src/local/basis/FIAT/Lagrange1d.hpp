#ifndef LAGRANGE1D_HPP
#define LAGRANGE1D_HPP


#include <vector>
#include <iostream>
#include "Teuchos_TestForException.hpp"

using std::vector;
using std::cout;
using std::endl;

namespace Lagrange
{
  // function for equispaced points on [a,b]
  template<typename ScalarType>
  void equispacedPoints( const int &n,
			 const ScalarType &a,
			 const ScalarType &b,
			 vector<ScalarType> &pts );

  // function for getting divided difference coefficients
  template<typename ScalarType>
  void dividedDifferences( const vector<ScalarType> &abscissa ,
			   const vector<ScalarType> &ordinates ,
			   vector<ScalarType> &coefficients );

  // third, need a function for evaluating polynomials in divided difference
  // format.  Two different types since I might do ad w.r.t. x, but not
  // the coefficients/pts.
  template<typename ScalarType1,
	   typename ScalarType2>
  ScalarType1 evaluateDividedDifferencePoly( const vector<ScalarType1> &abscissa ,
					     const vector<ScalarType1> &coefficients ,
					     ScalarType2 &x );
  

  // need a class for getting all Lagrange polynomials
  template<typename ScalarType>
  class Lagrange
  {
  public:
    Lagrange( const vector<ScalarType> &pts );
    ~Lagrange() {}
    
    template<typename ScalarType2>
    ScalarType2 eval( int i , ScalarType2 & x) {
      return evaluateDividedDifferencePoly( pts_ , coeffs_[i] , x );
    }

    int getDegree() const { return pts_.size() - 1; }

  private:
    vector<ScalarType> pts_;
    vector<vector<ScalarType> > coeffs_;
  };

}

#include "Lagrange1dDef.hpp"

#endif
