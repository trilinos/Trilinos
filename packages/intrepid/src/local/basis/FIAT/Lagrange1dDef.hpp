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
    
    TEST_FOR_EXCEPTION( pts.size() != n+1 ,
			std::invalid_argument,
			">>> ERROR (Lagrange::equispacedPoints): pts array improperly sized" );

    const ScalarType h = (b - a) / n;

    for (int i=0;i<=n;i++) {
      pts[i] = a + i * h;
    }

    return;
  }
  
  template<typename ScalarType>
  void dividedDifferences( const vector<ScalarType> &abscissa ,
			   const vector<ScalarType> &ordinates ,
			   vector<ScalarType> &coefficients )
  {
    TEST_FOR_EXCEPTION( abscissa.size() != ordinates.size() 
			|| abscissa.size() != coefficients.size() ,
			std::invalid_argument ,
			">>> ERROR (Lagrange::dividedDifferences): arrays not properly sized" );

    const int np = abscissa.size();
    const int order = np - 1;
    vector<ScalarType> work( abscissa.size() );
    for (int i=0;i<np;i++) {
      work[i] = ordinates[i];
    }
    for (int stage=0;stage<=order;stage++) {
      for (int i=0;i<=order-stage;i++) {
	work[i] = (work[i+stage]-work[i]) / (abscissa[i+stage] - abscissa[i]);
      }
      coefficients[stage] = work[0];
    }
    
    return;

  }

  template<typename ScalarType>
  ScalarType evalDividedDifferencePoly( const vector<ScalarType> &abscissa ,
					const vector<ScalarType> &coefficients ,
					ScalarType x )
  {
    TEST_FOR_EXCEPTION( abscissa.size() != coefficients.size() ,
			std::invalid_argument ,
			">>>ERROR (Lagrange::evalDividedDifferencePoly) : input arrays not same size" );
    const int np = coefficients.size();
    ScalarType val = coefficients[np-1];
    for (int i=np-2;i>=0;i--) {
      val = val * ( x - abscissa[i] ) + coefficients[i];
    }

    return val;
  }

  template<typename ScalarType>
  Lagrange<ScalarType>::Lagrange( const vector<ScalarType> &pts ) :
    pts_( pts ), coeffs_( pts.size() ) 
  {
    vector<ScalarType> vals(pts.size());

    for (int i=0;i<pts.size();i++) {
      vals[i] = 0.0;
    }
    
    for (int i=0;i<pts.size();i++) {
      vals[i] = 1.0;
      coeffs_[i].resize(pts.size());
      dividedDifferences( pts , vals , coeffs_[i] );
      vals[i] = 0.0;
    }
  }

  template<typename ScalarType>
  template<typename ScalarType2>
  Lagrange<ScalarType>::eval<ScalarType2>( int i , ScalarType2 &x )
  {
    return evaluateDividedDifferencePoly( pts_ , coeffs_[i] , x );
  }

}
