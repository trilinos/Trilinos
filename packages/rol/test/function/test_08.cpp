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

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_CompositeObjective_SimOpt.hpp"

typedef double RealT;

template<class Real> 
class ObjectiveFunctionTest08_1 : public ROL::Objective_SimOpt<Real> {
public:
  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    Real half(0.5), quadu(0), quadz(0);
    unsigned usize = up->size();
    for ( unsigned i = 0; i < usize; i++ ) {
      quadu += (*up)[i]*(*up)[i]; 
    }
    unsigned zsize = zp->size();
    for ( unsigned i = 0; i < zsize; i++ ) {
      quadz += (*zp)[i]*(*zp)[i]; 
    }
    return half*(quadu + quadz);
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > gp
      = dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    gp->assign(up->begin(),up->end());
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > gp
      = dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    gp->assign(zp->begin(),zp->end());
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > hvp
      = dynamic_cast<ROL::StdVector<Real>&>(hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    hvp->assign(vp->begin(),vp->end());
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > hvp
      = dynamic_cast<ROL::StdVector<Real>&>(hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    hvp->assign(vp->begin(),vp->end());
  }
};

template<class Real> 
class ObjectiveFunctionTest08_2 : public ROL::Objective_SimOpt<Real> {
public:
  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    Real linu(0), linz(0);
    unsigned usize = up->size();
    for ( unsigned i = 0; i < usize; i++ ) {
      linu += (*up)[i];
    }
    unsigned zsize = zp->size();
    for ( unsigned i = 0; i < zsize; i++ ) {
      linz += (*zp)[i];
    }
    return linu + linz;
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > gp
      = dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    gp->assign(up->size(),1);
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > gp
      = dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    gp->assign(zp->size(),1);
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }
};

template<class Real> 
class ObjectiveFunctionTest08_scalarize : public ROL::StdObjective<Real> {
public:
  Real value( const std::vector<Real> &x, Real &tol ) {
    return std::log(x[0]) * std::exp(x[1]);
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    g[0] = std::exp(x[1])/x[0];
    g[1] = std::exp(x[1]) * std::log(x[0]);
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    Real H11 = -std::exp(x[1])/(x[0]*x[0]);
    Real H12 = std::exp(x[1])/x[0];
    Real H21 = std::exp(x[1])/x[0];
    Real H22 = std::exp(x[1]) * std::log(x[0]);
    hv[0] = H11*v[0] + H12*v[1];
    hv[1] = H21*v[0] + H22*v[1];
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
  Teuchos::oblackholestream bhs; // outputs nothing
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
    unsigned dim = 2;
    ROL::Ptr<std::vector<RealT> > u_ptr = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT> > u = ROL::makePtr<ROL::StdVector<RealT>>(u_ptr);
    setRandomVector(*u_ptr);
    ROL::Ptr<std::vector<RealT> > z_ptr = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT> > z = ROL::makePtr<ROL::StdVector<RealT>>(z_ptr);
    setRandomVector(*z_ptr);
    ROL::Ptr<ROL::Vector<RealT> > x = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(u,z);
    ROL::Ptr<std::vector<RealT> > du_ptr = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT> > du = ROL::makePtr<ROL::StdVector<RealT>>(du_ptr);
    setRandomVector(*du_ptr);
    ROL::Ptr<std::vector<RealT> > dz_ptr = ROL::makePtr<std::vector<RealT>>(dim,0.0);
    ROL::Ptr<ROL::Vector<RealT> > dz = ROL::makePtr<ROL::StdVector<RealT>>(dz_ptr);
    setRandomVector(*dz_ptr);
    ROL::Ptr<ROL::Vector<RealT> > d = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(du,dz);
    // Build objective function
    std::vector<ROL::Ptr<ROL::Objective_SimOpt<RealT> > > vec_obj(2,ROL::nullPtr);
    vec_obj[0] = ROL::makePtr<ObjectiveFunctionTest08_1<RealT>>();
    vec_obj[1] = ROL::makePtr<ObjectiveFunctionTest08_2<RealT>>();
    ROL::Ptr<ROL::StdObjective<RealT> > obj_scalarize
      = ROL::makePtr<ObjectiveFunctionTest08_scalarize<RealT>>();
    ROL::Ptr<ROL::CompositeObjective_SimOpt<RealT> > obj
      = ROL::makePtr<ROL::CompositeObjective_SimOpt<RealT>>(vec_obj,obj_scalarize);
    // Test parametrized objective functions
    *outStream << "Check Derivatives of CompositeObjective_SimOpt\n";
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
