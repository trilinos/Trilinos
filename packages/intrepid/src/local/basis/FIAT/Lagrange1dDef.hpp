namespace Lagrange
{
  template<typename ScalarType>
  void equispacedPoints( const int & n,
			 const ScalarType & a,
			 const ScalarType & b,
			 vector<ScalarType> &pts )
  {
    TEST_FOR_EXCEPTION( n < 1 ,
			std::invalid_argument,
			">>> ERROR (Lagrange::equispacedPoints): n must be at least 1" );
    
    TEST_FOR_EXCEPTION( pts.size() != (unsigned)(n+1) ,
			std::invalid_argument,
			">>> ERROR (Lagrange::equispacedPoints): pts array improperly sized" );

    const ScalarType h = (b - a) / n;

    for (unsigned i=0;i<= (unsigned) n;i++) {
      pts[i] = a + i * h;
    }

    return;
  }
  
  template<typename ScalarType>
  void dividedDifferences( const vector<ScalarType> &x ,
			   const vector<ScalarType> &f ,
			   vector<ScalarType> &a )
  {
    TEST_FOR_EXCEPTION( x.size() != f.size() 
			|| x.size() != a.size() ,
			std::invalid_argument ,
			">>> ERROR (Lagrange::dividedDifferences): arrays not properly sized" );

    const int np = x.size();
    const int n = np - 1;
    vector<ScalarType> t( x.size() );
    for (int i=0;i<=n;i++) {
      t[i] = f[i];
      for (int j=i-1;j>=0;j--) {
	t[j] = ( t[j+1]-t[j])/(x[i]-x[j]);
      }
      a[i]=t[0];
    }
    
    return;
  }
  

  template<typename ScalarType1,
	   typename ScalarType2>
  ScalarType1 evaluateDividedDifferencePoly( const vector<ScalarType1> &abscissa ,
					     const vector<ScalarType1> &coefficients ,
					     ScalarType2 &x )
  {
    TEST_FOR_EXCEPTION( abscissa.size() != coefficients.size() ,
			std::invalid_argument ,
			">>>ERROR (Lagrange::evalDividedDifferencePoly) : input arrays not same size" );
    const unsigned np = coefficients.size();
    ScalarType1 val = coefficients[np-1];
    for (int i=abscissa.size()-2;i>=0;--i) {
      val = val * ( x - abscissa[i] ) + coefficients[i];
    }

    return val;
  }

  template<typename ScalarType>
  Lagrange<ScalarType>::Lagrange( const vector<ScalarType> &pts ) :
    pts_( pts ), coeffs_( pts.size() ) 
  {
    vector<ScalarType> vals(pts.size());

    for (unsigned i=0;i<pts.size();i++) {
      vals[i] = 0.0;
    }
    
    for (unsigned i=0;i<pts.size();i++) {
      vals[i] = 1.0;
      coeffs_[i].resize(pts.size());
      dividedDifferences( pts , vals , coeffs_[i] );
      vals[i] = 0.0;
    }
  }


}
