#include <vector>

using std::vector;

namespace Lagrange
{
  // function for equispaced points on [a,b]
  template<typename ScalarType>
  void equispacedPoints( int n,
			 ScalarType a,
			 ScalarType b,
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
  ScalarType evaluateDivideDifferencePoly( const vector<ScalarType> &abscissa ,
					   const vector<ScalarType> &coefficients ,
					   const ScalarType2 x );


  // need a class for getting all Lagrange polynomials
  template<typename ScalarType>
  class Lagrange
  {
  public:
    Lagrange( const vector<ScalarType> &pts );
    ~Lagrange() {}
    
    template<typename ScalarType2>
    ScalarType2 eval( int i , ScalarType2 & x);

    int getDegree() const { return pts_.size() - 1; }

  private:
    vector<ScalarType> pts_;
    vector<vector<ScalarType> > coeffs_;
  };

}
