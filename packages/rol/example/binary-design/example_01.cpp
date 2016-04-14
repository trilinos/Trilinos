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
    \brief Shows how to minimize a function with binary (0/1) constraints.
*/

#include "ROL_Algorithm.hpp"
#include "ROL_StdVector.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <iostream>

typedef double RealT;

template<class Real>
class BinaryDesignObjective : public ROL::Objective<Real> {
private:
  int nvars_;
  Real alpha_;

public:
  BinaryDesignObjective(int nvars, Real alpha) : nvars_(nvars), alpha_(alpha) {}

  RealT value(const ROL::Vector<Real> &x, Real &tol) {
    Teuchos::RCP<const std::vector<Real> > ex
        = Teuchos::dyn_cast<const ROL::StdVector<Real> >(x).getVector();

    RealT val(0);
    for (int i=0; i<nvars_; ++i) {
      val += (*ex)[i] + alpha_ * (*ex)[i] * (*ex)[i];
    }
    return val;
  }

  void gradient(ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real &tol) {
    Teuchos::RCP<const std::vector<Real> > ex
        = Teuchos::dyn_cast<const ROL::StdVector<Real> >(x).getVector();
    Teuchos::RCP<std::vector<Real> > eg
        = Teuchos::dyn_cast<ROL::StdVector<Real> >(g).getVector();

    Real two(2);
    for (int i=0; i<nvars_; ++i) {
      (*eg)[i] = 1 + alpha_*two*(*ex)[i];
    }
  }

  void hessVec(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol) {
    Teuchos::RCP<const std::vector<Real> > ex
        = Teuchos::dyn_cast<const ROL::StdVector<Real> >(x).getVector();
    Teuchos::RCP<const std::vector<Real> > ev
        = Teuchos::dyn_cast<const ROL::StdVector<Real> >(v).getVector();
    Teuchos::RCP<std::vector<Real> > ehv
        = Teuchos::dyn_cast<ROL::StdVector<Real> >(hv).getVector();

    Real two(2);
    for (int i=0; i<nvars_; ++i) {
      (*ehv)[i] = alpha_*two*(*ev)[i];
    }
  }

};

template<class Real>
class BinaryDesignEqualityConstraint : public ROL::EqualityConstraint<Real> {
private:
  int nvars_;
  Real vol_;

public:
  BinaryDesignEqualityConstraint(const int &nvars, const Real &vol) : nvars_(nvars), vol_(vol) {}

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &x, Real &tol) {
    Teuchos::RCP<const std::vector<Real> > ex
        = Teuchos::dyn_cast<const ROL::StdVector<Real> >(x).getVector();
    Teuchos::RCP<std::vector<Real> > ec
        = Teuchos::dyn_cast<ROL::StdVector<Real> >(c).getVector();

    Real one(1);
    for (int i=0; i<nvars_; ++i) {
      (*ec)[i] = (*ex)[i] * ((*ex)[i] - one);
    }
    (*ec)[nvars_] = -vol_;
    for (int i=0; i<nvars_; ++i) {
      (*ec)[nvars_] += (*ex)[i];
    }
  }

  void applyJacobian(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol) {
    Teuchos::RCP<const std::vector<Real> > ex
        = Teuchos::dyn_cast<const ROL::StdVector<Real> >(x).getVector();
    Teuchos::RCP<const std::vector<Real> > ev
        = Teuchos::dyn_cast<const ROL::StdVector<Real> >(v).getVector();
    Teuchos::RCP<std::vector<Real> > ejv
        = Teuchos::dyn_cast<ROL::StdVector<Real> >(jv).getVector();

    Real one(1), two(2);
    for (int i=0; i<nvars_; ++i) {
      (*ejv)[i] = (two*(*ex)[i]-one) * (*ev)[i];
    }
    (*ejv)[nvars_] = 0;
    for (int i=0; i<nvars_; ++i) {
      (*ejv)[nvars_] += (*ev)[i];
    }
  }

  void applyAdjointJacobian(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol) {
    Teuchos::RCP<const std::vector<Real> > ex
        = Teuchos::dyn_cast<const ROL::StdVector<Real> >(x).getVector();
    Teuchos::RCP<const std::vector<Real> > ev
        = Teuchos::dyn_cast<const ROL::StdVector<Real> >(v).getVector();
    Teuchos::RCP<std::vector<Real> > eajv
        = Teuchos::dyn_cast<ROL::StdVector<Real> >(ajv).getVector();

    Real one(1), two(2);
    for (int i=0; i<nvars_; ++i) {
      (*eajv)[i] = (two*(*ex)[i]-one) * (*ev)[i] + (*ev)[nvars_];
    }
  }

  void applyAdjointHessian(ROL::Vector<Real> &ahuv, const ROL::Vector<Real> &u, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol) {
    Teuchos::RCP<const std::vector<Real> > ex
        = Teuchos::dyn_cast<const ROL::StdVector<Real> >(x).getVector();
    Teuchos::RCP<const std::vector<Real> > eu
        = Teuchos::dyn_cast<const ROL::StdVector<Real> >(u).getVector();
    Teuchos::RCP<const std::vector<Real> > ev
        = Teuchos::dyn_cast<const ROL::StdVector<Real> >(v).getVector();
    Teuchos::RCP<std::vector<Real> > eahuv
        = Teuchos::dyn_cast<ROL::StdVector<Real> >(ahuv).getVector();

    Real two(2);
    for (int i=0; i<nvars_; ++i) {
      (*eahuv)[i] = two * (*eu)[i] * (*ev)[i];
    }
  }

};

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Example body.
 
  try {

    // Set up problem data
    int   dim   = 10; // Set problem dimension. 
    RealT vol   = 2;  // Set desired volume. 
    RealT alpha = 1;  // Set quadratic penalty. 
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT>(dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > g_rcp = Teuchos::rcp( new std::vector<RealT>(dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > d_rcp = Teuchos::rcp( new std::vector<RealT>(dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > v_rcp = Teuchos::rcp( new std::vector<RealT>(dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > jv_rcp = Teuchos::rcp( new std::vector<RealT>(dim+1, 0.0) );
    Teuchos::RCP<std::vector<RealT> > ajv_rcp = Teuchos::rcp( new std::vector<RealT>(dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > c_rcp = Teuchos::rcp( new std::vector<RealT>(dim+1, 0.0) );
    for (int i=0; i<dim; i++) {
      (*x_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*g_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*d_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*v_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*ajv_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    for (int i=0; i<dim+1; i++) {
      (*jv_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*c_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    Teuchos::RCP<ROL::Vector<RealT> > x = Teuchos::rcp( new ROL::StdVector<RealT>(x_rcp) );
    Teuchos::RCP<ROL::Vector<RealT> > g = Teuchos::rcp( new ROL::StdVector<RealT>(g_rcp) );
    Teuchos::RCP<ROL::Vector<RealT> > d = Teuchos::rcp( new ROL::StdVector<RealT>(d_rcp) );
    Teuchos::RCP<ROL::Vector<RealT> > v = Teuchos::rcp( new ROL::StdVector<RealT>(v_rcp) );
    Teuchos::RCP<ROL::Vector<RealT> > jv = Teuchos::rcp( new ROL::StdVector<RealT>(jv_rcp) );
    Teuchos::RCP<ROL::Vector<RealT> > ajv = Teuchos::rcp( new ROL::StdVector<RealT>(ajv_rcp) );
    Teuchos::RCP<ROL::Vector<RealT> > c = Teuchos::rcp( new ROL::StdVector<RealT>(c_rcp) );
    Teuchos::RCP<ROL::Objective<RealT> > obj = Teuchos::rcp(new BinaryDesignObjective<RealT>(dim, alpha));
    Teuchos::RCP<ROL::EqualityConstraint<RealT> > con = Teuchos::rcp(new BinaryDesignEqualityConstraint<RealT>(dim, vol));

   // Define algorithm
    Teuchos::RCP<Teuchos::ParameterList> parlist
      = Teuchos::rcp(new Teuchos::ParameterList());
    std::string paramfile = "input.xml";
    Teuchos::updateParametersFromXmlFile(paramfile,parlist.ptr());
    ROL::Algorithm<RealT> algo("Composite Step",*parlist);

    // Test objective
    obj->checkGradient(*x, *d, true, *outStream);
    *outStream << "\n"; 
    obj->checkHessVec(*x, *v, true, *outStream);
    *outStream << "\n";
    obj->checkHessSym(*x, *d, *v, true, *outStream);
    *outStream << "\n";
    // Test constraint.
    con->checkApplyJacobian(*x, *v, *jv, true, *outStream);
    con->checkAdjointConsistencyJacobian(*ajv, *v, *x, true, *outStream);
    con->checkApplyAdjointHessian(*x, *ajv, *v, *v, true, *outStream);

    // Run algorithm
    for (int i=0; i<dim; ++i) {
      (*x_rcp)[i] = 1.234*(i<2);
    }
    for (int i=0; i<dim+1; ++i) {
      (*c_rcp)[i] = 0.0;
    }
    algo.run(*x, *c, *obj, *con, true, *outStream);
    *outStream << "x = [";
    for (int i=0; i<dim; ++i) {
      *outStream << (*x_rcp)[i] << "  ";
    }
    *outStream << "]\n";

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

