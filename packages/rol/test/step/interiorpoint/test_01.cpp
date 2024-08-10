// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_InteriorPointPenalty.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_ParameterList.hpp"

#include <iomanip>


/*! \file test_01.cpp
    \brief Verify that the interior point log-barrier penalized objective
           passes gradient and Hessian checks.
*/

template<class Real>
class NullObjective : public ROL::Objective<Real> {
  typedef ROL::Vector<Real>  V;
public:
  Real value( const V &x, Real &tol ) {
    return Real(0.0);
  }
  void gradient( V &g, const V &x, Real &tol ) {
    g.zero();
  }
  void hessVec( V &hv, const V &v, const V &x, Real &tol ) {
    hv.zero();
  }
};

template<class Real>
void printVector( const ROL::Vector<Real> &x, std::ostream &outStream ) {
  ROL::Ptr<const std::vector<Real> > xp = 
    dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();

  for( size_t i=0; i<xp->size(); ++i ) {
    outStream << (*xp)[i] << std::endl;
  }
}




typedef double RealT;

int main(int argc, char *argv[]) {

  typedef std::vector<RealT>          vector;

  typedef ROL::Vector<RealT>          V;
  typedef ROL::StdVector<RealT>       SV;
  typedef ROL::Objective<RealT>       OBJ;
  typedef ROL::BoundConstraint<RealT> BND;

  typedef ROL::ParameterList      PL;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs;
  if( iprint > 0 ) 
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try {

    PL parlist;
    PL &iplist = parlist.sublist("Step").sublist("Primal Dual Interior Point");
    PL &lblist = iplist.sublist("Barrier Objective");

    RealT mu = 1.e2;
    RealT kappaD = 1.e-4;
    bool useLinearDamping = true;

    lblist.set("Use Linear Damping", useLinearDamping);
    lblist.set("Linear Damping Coefficient",kappaD);
    lblist.set("Initial Barrier Parameter",mu);

    RealT ninf = ROL::ROL_NINF<RealT>();
    RealT inf  = ROL::ROL_INF<RealT>();

    int dim = 4;
    int numTestVectors = 19;
 
    ROL::Ptr<vector> x_ptr  = ROL::makePtr<vector>(dim, 0.0); 
    ROL::Ptr<vector> d_ptr  = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> v_ptr  = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> l_ptr  = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> u_ptr  = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> e0_ptr = ROL::makePtr<vector>(dim, 0.0); // First canonical vector

    (*e0_ptr)[0] = 1.0;

    SV e0(e0_ptr);

    // Lower Bounds         // Upper Bounds
    (*l_ptr)[0] = ninf;     (*u_ptr)[0] = 5.0;
    (*l_ptr)[1] = ninf;     (*u_ptr)[1] = inf;
    (*l_ptr)[2] = -5.0;     (*u_ptr)[2] = inf;
    (*l_ptr)[3] = -5.0;     (*u_ptr)[3] = 5.0;

    RealT left = -1.0;  RealT right = 1.0;

    RealT xmax = 4.99;

    ROL::Ptr<V> x  = ROL::makePtr<SV>( x_ptr );
    ROL::Ptr<V> d  = ROL::makePtr<SV>( d_ptr );
    ROL::Ptr<V> v  = ROL::makePtr<SV>( v_ptr );
    ROL::Ptr<V> l  = ROL::makePtr<SV>( l_ptr );
    ROL::Ptr<V> u  = ROL::makePtr<SV>( u_ptr );

    ROL::Ptr<const V> maskL, maskU;

    ROL::RandomizeVector(*d,left,right);
    ROL::RandomizeVector(*v,left,right);

    std::vector<RealT>   values(numTestVectors);        // Computed objective value for each
    std::vector<RealT>   exact_values(numTestVectors);

    std::vector<ROL::Ptr<V> > x_test;

    for(int i=0; i<numTestVectors; ++i) {
      x_test.push_back(x->clone());
      RealT t = static_cast<RealT>(i)/static_cast<RealT>(numTestVectors-1);
      RealT fillValue = xmax*(2.0*t-1.0);
      x_test[i]->applyUnary(ROL::Elementwise::Fill<RealT>(fillValue));
    }

    ROL::Ptr<OBJ> obj = ROL::makePtr<NullObjective<RealT>>();
    ROL::Ptr<BND> bnd = ROL::makePtr<ROL::Bounds<RealT>>(l,u);

    ROL::InteriorPointPenalty<RealT> ipobj(obj,bnd,parlist);

    maskL = ipobj.getLowerMask();
    maskU = ipobj.getUpperMask();

    ROL::Ptr<const std::vector<RealT> > maskL_ptr = dynamic_cast<const SV&>(*maskL).getVector();
    ROL::Ptr<const std::vector<RealT> > maskU_ptr = dynamic_cast<const SV&>(*maskU).getVector();

    *outStream << "\nLower bound vector" << std::endl;
    printVector(*l,*outStream);

    *outStream << "\nUpper bound vector" << std::endl;
    printVector(*u,*outStream);

    *outStream << "\nLower mask vector" << std::endl;
    printVector(*maskL, *outStream);

    *outStream << "\nUpper mask vector" << std::endl;
    printVector(*maskU, *outStream);

    *outStream << "\nChecking Objective value" << std::endl;

    RealT tol = std::sqrt(ROL::ROL_EPSILON<RealT>());
    *outStream   << std::setw(16) << "x[i], i=0,1,2,3"
                 << std::setw(20) << "Computed Objective"
                 << std::setw(20) << "Exact Objective" << std::endl;

    RealT valueError(0.0);

    for(int i=0; i<numTestVectors; ++i) {
      values[i] = ipobj.value(*(x_test[i]),tol);

      exact_values[i] = 0;

      // Extract the value from the test vector that is in every element
      RealT xval = x_test[i]->dot(e0);


      for(int j=0; j<dim; ++j) {
        if( (*maskL_ptr)[j] ) {
          RealT diff = xval-(*l_ptr)[j];
          exact_values[i] -= mu*std::log(diff);

          if( useLinearDamping && !(*maskU_ptr)[j] ) {
            exact_values[i] += mu*kappaD*diff;
          }

        }
        if( (*maskU_ptr)[j] ) {
          RealT diff = (*u_ptr)[j]-xval;
          exact_values[i] -= mu*std::log(diff);

          if(useLinearDamping && !(*maskL_ptr)[j] ) {
            exact_values[i] += mu*kappaD*diff;
          }

        }
      } // end loop over elements

      *outStream << std::setw(16) << xval
                 << std::setw(20) << values[i]
                 << std::setw(20) << exact_values[i] << std::endl;
      RealT valDiff = exact_values[i] - values[i];
      valueError += valDiff*valDiff;
    } // end loop over vectors

    if(valueError>ROL::ROL_EPSILON<RealT>()) {
      errorFlag++;
    }

    *outStream << "\nPerforming finite difference checks" << std::endl;

    ipobj.checkGradient(*x,*v,true,*outStream);       *outStream << std::endl;
    ipobj.checkHessVec(*x,*d,true,*outStream);        *outStream << std::endl;
    ipobj.checkHessSym(*x,*d,*v,true,*outStream);     *outStream << std::endl;

  }
  catch (std::logic_error& err) {
    *outStream << err.what() << std::endl;
    errorFlag = -1000;
  }

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED" << std::endl;
  else
    std::cout << "End Result: TEST PASSED" << std::endl;

  return 0;
}
