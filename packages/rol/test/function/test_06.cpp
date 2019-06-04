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

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_StdObjective.hpp"

typedef double RealT;

template<class Real> 
class ObjectiveFunctionTest06 : public ROL::StdObjective<Real> {
public:
  Real value( const std::vector<Real> &x, Real &tol ) {
    Real quad(0), lin(0);
    unsigned size = x.size();
    std::vector<Real> p(size+2,1);
    for ( unsigned i = 0; i < size; i++ ) {
      quad += x[i]*x[i]; 
      lin  += x[i]*p[i+1];
    }
    return std::exp(p[0])*quad + lin + p[size+1];
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    unsigned size = x.size();
    std::vector<Real> p(size+2,1);
    const Real two(2);
    for ( unsigned i = 0; i < size; i++ ) {
      g[i] = two*std::exp(p[0])*x[i] + p[i+1];
    }
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    unsigned size = x.size();
    std::vector<Real> p(size+2,1);
    const Real two(2);
    for ( unsigned i = 0; i < size; i++ ) {
      hv[i] = two*std::exp(p[0])*v[i]; 
    } 
  }
};

void setRandomVector(std::vector<RealT> &x) {
  unsigned dim = x.size();
  for ( unsigned i = 0; i < dim; i++ ) {
    x[i] = (RealT)rand()/(RealT)RAND_MAX;
  }
}

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build vectors
    unsigned dim = 4;
    ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT> > x = ROL::makePtr<ROL::StdVector<RealT>>(x_ptr);
    setRandomVector(*x_ptr);
    ROL::Ptr<std::vector<RealT> > d_ptr = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT> > d = ROL::makePtr<ROL::StdVector<RealT>>(d_ptr);
    setRandomVector(*d_ptr);
    // Build objective function
    ROL::Ptr<ROL::StdObjective<RealT> > obj =
      ROL::makePtr<ObjectiveFunctionTest06<RealT>>();
    // Test parametrized objective functions
    *outStream << "Check Derivatives of StdObjective\n";
    obj->checkGradient(*x,*d,true,*outStream);
    obj->checkHessVec(*x,*d,true,*outStream);
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
