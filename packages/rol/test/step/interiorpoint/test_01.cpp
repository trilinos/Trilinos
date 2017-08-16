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
// @HEADER


#include "ROL_InteriorPointPenalty.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Bounds.hpp"

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
  Teuchos::RCP<const std::vector<Real> > xp = 
    Teuchos::dyn_cast<const ROL::StdVector<Real> >(x).getVector();

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

  typedef Teuchos::ParameterList      PL;

  using Teuchos::RCP; using Teuchos::rcp;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int iprint = argc - 1;
  RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs;
  if( iprint > 0 ) 
    outStream = rcp(&std::cout,false);
  else
    outStream = rcp(&bhs,false);

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
 
    RCP<vector> x_rcp  = rcp( new vector(dim, 0.0) ); 
    RCP<vector> d_rcp  = rcp( new vector(dim, 0.0) );
    RCP<vector> v_rcp  = rcp( new vector(dim, 0.0) );
    RCP<vector> l_rcp  = rcp( new vector(dim, 0.0) );
    RCP<vector> u_rcp  = rcp( new vector(dim, 0.0) );
    RCP<vector> e0_rcp = rcp( new vector(dim, 0.0) ); // First canonical vector

    (*e0_rcp)[0] = 1.0;

    SV e0(e0_rcp);

    // Lower Bounds         // Upper Bounds
    (*l_rcp)[0] = ninf;     (*u_rcp)[0] = 5.0;
    (*l_rcp)[1] = ninf;     (*u_rcp)[1] = inf;
    (*l_rcp)[2] = -5.0;     (*u_rcp)[2] = inf;
    (*l_rcp)[3] = -5.0;     (*u_rcp)[3] = 5.0;

    RealT left = -1.0;  RealT right = 1.0;

    RealT xmax = 4.99; 

    RCP<V> x  = rcp( new SV( x_rcp ) );
    RCP<V> d  = rcp( new SV( d_rcp ) );
    RCP<V> v  = rcp( new SV( v_rcp ) );
    RCP<V> l  = rcp( new SV( l_rcp ) );
    RCP<V> u  = rcp( new SV( u_rcp ) );

    RCP<const V> maskL, maskU;

    ROL::RandomizeVector(*d,left,right);
    ROL::RandomizeVector(*v,left,right);

    std::vector<RealT>   values(numTestVectors);        // Computed objective value for each 
    std::vector<RealT>   exact_values(numTestVectors);  

    std::vector<RCP<V> > x_test;

    for(int i=0; i<numTestVectors; ++i) {
      x_test.push_back(x->clone());
      RealT t = static_cast<RealT>(i)/static_cast<RealT>(numTestVectors-1);
      RealT fillValue = xmax*(2.0*t-1.0);
      x_test[i]->applyUnary(ROL::Elementwise::Fill<RealT>(fillValue));
    }

    RCP<OBJ> obj = rcp( new NullObjective<RealT> );
    RCP<BND> bnd = rcp( new ROL::Bounds<RealT>(l,u) );

    ROL::InteriorPointPenalty<RealT> ipobj(obj,bnd,parlist);

    maskL = ipobj.getLowerMask();
    maskU = ipobj.getUpperMask();

    RCP<const std::vector<RealT> > maskL_rcp = Teuchos::dyn_cast<const SV>(*maskL).getVector();
    RCP<const std::vector<RealT> > maskU_rcp = Teuchos::dyn_cast<const SV>(*maskU).getVector();

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
        if( (*maskL_rcp)[j] ) {
          RealT diff = xval-(*l_rcp)[j];
          exact_values[i] -= mu*std::log(diff);

          if( useLinearDamping && !(*maskU_rcp)[j] ) {
            exact_values[i] += mu*kappaD*diff;
          }

        }
        if( (*maskU_rcp)[j] ) {
          RealT diff = (*u_rcp)[j]-xval;
          exact_values[i] -= mu*std::log(diff);

          if(useLinearDamping && !(*maskL_rcp)[j] ) {
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
  catch (std::logic_error err) {
    *outStream << err.what() << std::endl;
    errorFlag = -1000;
  }

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED" << std::endl;
  else
    std::cout << "End Result: TEST PASSED" << std::endl;

  return 0;
}

