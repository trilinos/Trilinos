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

/*! \file  example_01.cpp
    \brief Shows how to minimize a quadratic functional
*/

#include "ROL_Algorithm.hpp"
#include "ROL_ConjugateGradients.hpp"

#include "ROL_LinearOperator.hpp"
#include "ROL_QuadraticObjective.hpp"
#include "ROL_StdVector.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <iostream>

typedef double RealT;

template<class Real>
class matrix : public ROL::LinearOperator<Real> {
private:
  const int dim_;

public:
  matrix(const int dim) : dim_(dim) {}

  void apply( ROL::Vector<Real> &Hv, const ROL::Vector<Real> &v, Real &tol ) const {
    ROL::Ptr<std::vector<Real> > Hvp =
      dynamic_cast<ROL::StdVector<Real>&>(Hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    for (int i = 0; i < dim_; i++) {
      (*Hvp)[i] = 2.0*(*vp)[i];
      if ( i > 0 ) {
        (*Hvp)[i] -= (*vp)[i-1];
      }
      if ( i < dim_-1 ) {
        (*Hvp)[i] -= (*vp)[i+1];
      }
    }
  }
};

template<class Real>
class precond : public ROL::LinearOperator<Real> {
private:
  const int dim_;

public:
  precond(const int dim) : dim_(dim) {}

  void apply( ROL::Vector<Real> &Hv, const ROL::Vector<Real> &v, Real &tol ) const {
    ROL::Ptr<std::vector<Real> > Hvp =
      dynamic_cast<ROL::StdVector<Real>&>(Hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    for (int i = 0; i < dim_; i++) {
      (*Hvp)[i] = 2.0*(*vp)[i];
    }
  }

  void applyInverse( ROL::Vector<Real> &Hv, const ROL::Vector<Real> &v, Real &tol ) const {
    ROL::Ptr<std::vector<Real> > Hvp =
      dynamic_cast<ROL::StdVector<Real>&>(Hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    for (int i = 0; i < dim_; i++) {
      (*Hvp)[i] = 0.5*(*vp)[i];
    }
  }
};

int main(int argc, char *argv[]) {
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

  // *** Example body.
 
  try {

    // Set up problem data
    int dim = 10; // Set problem dimension. 
    ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > g_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > d_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > v_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    for (int i=0; i<dim; i++) {
      (*x_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*g_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*d_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*v_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    ROL::Ptr<ROL::Vector<RealT> > x = ROL::makePtr<ROL::StdVector<RealT>>(x_ptr);
    ROL::Ptr<ROL::Vector<RealT> > g = ROL::makePtr<ROL::StdVector<RealT>>(g_ptr);
    ROL::Ptr<ROL::Vector<RealT> > d = ROL::makePtr<ROL::StdVector<RealT>>(d_ptr);
    ROL::Ptr<ROL::Vector<RealT> > v = ROL::makePtr<ROL::StdVector<RealT>>(v_ptr);
    ROL::Ptr<ROL::LinearOperator<RealT> > op = ROL::makePtr<matrix<RealT>>(dim);
    ROL::Ptr<ROL::Objective<RealT> > obj = ROL::makePtr<ROL::QuadraticObjective<RealT>>(op,g);

   // Define algorithm
    Teuchos::RCP<Teuchos::ParameterList> parlist
      = Teuchos::rcp( new Teuchos::ParameterList() );
    std::string paramfile = "input.xml";
    Teuchos::updateParametersFromXmlFile(paramfile,parlist.ptr());
    ROL::Algorithm<RealT> algo("Trust-Region",*parlist);

    // Test objective
    obj->checkGradient(*x, *d, true, *outStream);
    *outStream << "\n"; 
    obj->checkHessVec(*x, *v, true, *outStream);
    *outStream << "\n";
    obj->checkHessSym(*x, *d, *v, true, *outStream);
    *outStream << "\n";

    // Run algorithm
    algo.run(*x, *obj, true, *outStream);

    // Solve using Krylov
    RealT absTol = 1.e-4, relTol = 1.e-2;
    int iter = 2*dim+1;
    ROL::Ptr<ROL::Krylov<RealT> > kv
      = ROL::makePtr<ROL::ConjugateGradients<RealT>>(absTol,relTol,iter,false);
    ROL::Ptr<ROL::LinearOperator<RealT> > M = ROL::makePtr<precond<RealT>>(dim);
    int iterCG = 0, flagCG = 0;
    kv->run(*v,*op,*g,*M,iterCG,flagCG);
    v->scale(-1.0);

    // Check error
    d->set(*x);
    d->axpy(-1.0,*v);
    RealT abserr = d->norm();
    *outStream << std::scientific << "\n   Absolute Error: " << abserr << std::endl;
    if ( abserr > std::sqrt(ROL::ROL_EPSILON<RealT>()) ) {
      errorFlag += 1;
    }
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

