// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Objective.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_VectorWorkspace.hpp"
#include "ROL_Stream.hpp"

#include <iomanip>

using RealT = double;
using size_type = typename std::vector<RealT>::size_type;

namespace details {

using namespace ROL;

// An alternative, more expensive way to compute the dot product!
template<typename Real> 
class PolarizationIdentity {
private:

  mutable VectorWorkspace<Real> workspace_;

public:

  Real dot( const Vector<Real>& x, const Vector<Real>& y ) {

    auto w = workspace_.copy(x);
    auto z = workspace_.copy(x);
    w->plus(y);
    z->axpy(-1.0,y);
    auto result = 0.25*( w->dot(*w) - z->dot(*z) );
    return result;
  } // w and z should go out of scope - decrement ref counts

  void status( std::ostream& os ) const { workspace_.status(os); }
};




} // namespace details

using details::PolarizationIdentity;


int main( int argc, char* argv[] ) {

  using namespace ROL;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  auto outStream = makeStreamPtr( std::cout, argc > 1 );
  
  int errorFlag  = 0;
  RealT errtol = std::sqrt(ROL_EPSILON<RealT>());

  try {

    PolarizationIdentity<RealT> polar;

    size_type N = 20;

    auto xp = makePtr<std::vector<RealT>>(N);
    auto x = makePtr<StdVector<RealT>>(xp);
    auto y = x->clone();

    auto xx = PartitionedVector<RealT>::create({x,x});
    auto yy = PartitionedVector<RealT>::create({y,y});

    RandomizeVector( *x );
    RandomizeVector( *y );

    auto result = 0.5*polar.dot(*x,*y);
    result     += 0.5*polar.dot(*y,*x);

    auto x_dot_y = x->dot(*y);
    errorFlag += ( std::abs( x_dot_y - result ) > errtol );

    *outStream << std::setprecision(16) <<  x_dot_y << std::endl;
    *outStream << std::setprecision(16) <<  result  << std::endl;
 

    polar.dot(*xx,*yy);

    polar.status(*outStream);
    
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

