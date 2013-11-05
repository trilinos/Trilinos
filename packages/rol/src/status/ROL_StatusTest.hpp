//@HEADER
// ***********************************************************************
//
//                     Rapid Optimization Library
//
// Questions? Contact:    Drew Kouri (dpkouri@sandia.gov)
//                      Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ROL_STATUSTEST_H
#define ROL_STATUSTEST_H

/** \class ROL::StatusTest
    \brief Provides an interface to check status of optimization algorithms.
*/


namespace ROL {


template <class Real>
class StatusTest {
private:

  Real gtol_;
  Real stol_;
  int  max_iter_;

public:

  virtual ~StatusTest() {}

  StatusTest( Real gtol = 1.e-6, Real stol = 1.e-12, int max_iter = 100 ) :  
    gtol_(gtol), stol_(stol), max_iter_(max_iter) {}

  /** \brief Check algorithm status.
  */
  virtual bool check( AlgorithmState<Real> &state ) {
     if ( state.gnorm > gtol_ && state.snorm > stol_ && state.iter < max_iter_ ) {
       return true;
     }
     else {
       return false;
     }
  }

}; // class StatusTest

} // namespace ROL

#endif
