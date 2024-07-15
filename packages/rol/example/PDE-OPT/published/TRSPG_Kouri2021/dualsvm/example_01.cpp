// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to minimize a function with binary (0/1) constraints.
*/

#include "ROL_Solver.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_Stream.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <iostream>
#include <random>

typedef double RealT;

template<class Real>
class DualSVMObjective : public ROL::StdObjective<Real> {
private:
  const int nrows_;
  std::vector<std::vector<Real>> K_;
  int fval_, grad_, hess_;

  Real kernel(const std::vector<Real> &x, const std::vector<Real> &y) const {
    typename std::vector<Real>::size_type size = x.size();
    assert(size==y.size());
    Real val(0);
    for (typename std::vector<Real>::size_type i = 0; i < size; ++i) val += x[i]*y[i];
    return val;
  }

public:
  DualSVMObjective(const std::vector<std::vector<Real>> &X, const std::vector<Real> &y)
    : nrows_(X.size()), fval_(0), grad_(0), hess_(0) {
    K_.resize(nrows_);
    for (int i = 0; i < nrows_; ++i) {
      K_[i].resize(nrows_);
      for (int j = 0; j < nrows_; ++j) {
        K_[i][j] = y[i]*y[j]*kernel(X[i],X[j]);
      }
    }
  }

  void summarize(std::ostream &stream) const {
    stream << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << "  LassoObjective::summarize" << std::endl;
    stream << "    Number of calls to value:           " << fval_ << std::endl;
    stream << "    Number of calls to gradient:        " << grad_ << std::endl;
    stream << "    Number of calls to hessVec:         " << hess_ << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << std::endl;
  }

  void reset(void) {
    fval_ = 0; grad_ = 0; hess_ = 0;
  }

  Real value(const std::vector<Real> &x, Real &tol) {
    fval_++;
    Real val(0);
    for (int i = 0; i < nrows_; ++i) {
      val -= x[i];
      for (int j = 0; j < nrows_; ++j) {
        val += static_cast<Real>(0.5)*x[i]*K_[i][j]*x[j];
      }
    }
    return val;
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    grad_++;
    g.assign(x.size(),static_cast<Real>(-1));
    for (int i = 0; i < nrows_; ++i) {
      for (int j = 0; j < nrows_; ++j) {
        g[i] += K_[i][j]*x[j];
      }
    }
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    hess_++;
    hv.assign(x.size(),static_cast<Real>(0));
    for (int i = 0; i < nrows_; ++i) {
      for (int j = 0; j < nrows_; ++j) {
        hv[i] += K_[i][j]*v[j];
      }
    }
  }

};

template<typename Real>
class DualSVMConstraint : public ROL::StdConstraint<Real> {
private:
  const std::vector<Real> y_;
  const int size_;

public:
  DualSVMConstraint(const std::vector<Real> &y) : y_(y), size_(y.size()) {}

  void value(std::vector<Real> &c, const std::vector<Real> &x, Real &tol) {
    c[0] = static_cast<Real>(0);
    for (int i = 0; i < size_; ++i) c[0] += x[i]*y_[i]; 
  }

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    jv[0] = static_cast<Real>(0);
    for (int i = 0; i < size_; ++i) jv[0] += v[i]*y_[i]; 
  }

  void applyAdjointJacobian(std::vector<Real> &ajv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    ajv.assign(size_,static_cast<Real>(0));
    for (int i = 0; i < size_; ++i) ajv[i] = v[0]*y_[i]; 
  }

  void applyAdjointHessian(std::vector<Real> &ahuv, const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    ahuv.assign(size_,static_cast<Real>(0));
  }
};

int main(int argc, char *argv[]) {
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

  // *** Example body.
 
  try {

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Set up problem data
    int nInstances = 50;
    int nVars      = 100;
    unsigned seed  = 111;
    std::default_random_engine generator(seed);
    std::normal_distribution<RealT> dist(0.0,1.0);
    std::vector<std::vector<RealT>> X(nInstances);
    for (int i = 0; i < nInstances; ++i) {
      for (int j = 0; j < nVars; ++j) {
        X[i].push_back(dist(generator));
      }
    }
    std::vector<RealT> w(nVars);
    for (int j = 0; j < nVars; ++j) {
      w[j] = dist(generator);
    }
    RealT yval(0);
    std::vector<RealT> y(nInstances,0);
    for (int i = 0; i < nInstances; ++i) {
      yval = dist(generator);
      for (int j = 0; j < nVars; ++j) {
        yval += X[i][j]*w[j];
      }
      y[i] = (yval < static_cast<RealT>(0) ? static_cast<RealT>(-1)
              : (yval > static_cast<RealT>(0) ? static_cast<RealT>(1)
              : static_cast<RealT>(0)));
    }
    ROL::Ptr<DualSVMObjective<RealT>>  obj = ROL::makePtr<DualSVMObjective<RealT>>(X,y);
    ROL::Ptr<DualSVMConstraint<RealT>> con = ROL::makePtr<DualSVMConstraint<RealT>>(y);
    ROL::Ptr<ROL::StdVector<RealT>>    mul = ROL::makePtr<ROL::StdVector<RealT>>(1,0);
    ROL::Ptr<ROL::StdVector<RealT>>      x = ROL::makePtr<ROL::StdVector<RealT>>(nInstances,0);
    ROL::Ptr<ROL::StdVector<RealT>>      l = ROL::makePtr<ROL::StdVector<RealT>>(nInstances,0);
    ROL::Ptr<ROL::Bounds<RealT>>       bnd = ROL::makePtr<ROL::Bounds<RealT>>(*l,true);
    ROL::Ptr<ROL::Problem<RealT>> prob
      = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    prob->addBoundConstraint(bnd);
    prob->addLinearConstraint("Dual SVM",con,mul);
    prob->setProjectionAlgorithm(*parlist);
    prob->finalize(false,true,*outStream);
    prob->check(true,*outStream);

    ROL::Ptr<ROL::Solver<RealT>> solver;
    std::vector<RealT> vals(6);

    // SPG
    obj->reset();
    x->zero();
    parlist->sublist("Step").set("Type","Spectral Gradient");
    Teuchos::Time spgTimer("SPG Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    spgTimer.stop();
    *outStream << "Total optimization time = " << spgTimer.totalElapsedTime() << " seconds.\n";
    obj->summarize(*outStream);
    vals[0] = solver->getAlgorithmState()->value;

    // PQN
    obj->reset();
    x->zero();
    parlist->sublist("Step").set("Type","Line Search");
    parlist->sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type","Quasi-Newton Method");
    Teuchos::Time pqnTimer("PQN Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    pqnTimer.stop();
    *outStream << "Total optimization time = " << pqnTimer.totalElapsedTime() << " seconds.\n";
    obj->summarize(*outStream);
    vals[1] = solver->getAlgorithmState()->value;

    // LMTR
    obj->reset();
    x->zero();
    parlist->sublist("Step").set("Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","Lin-More");
    Teuchos::Time lmtrTimer("LMTR Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    lmtrTimer.stop();
    *outStream << "Total optimization time = " << lmtrTimer.totalElapsedTime() << " seconds.\n";
    obj->summarize(*outStream);
    vals[2] = solver->getAlgorithmState()->value;

    // TRSPG
    obj->reset();
    x->zero();
    parlist->sublist("Step").set("Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","SPG");
    Teuchos::Time trspgTimer("TRSPG Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    trspgTimer.stop();
    *outStream << "Total optimization time = " << trspgTimer.totalElapsedTime() << " seconds.\n";
    obj->summarize(*outStream);
    vals[3] = solver->getAlgorithmState()->value;

    // AL-LMTR
    obj->reset();
    x->zero();
    mul->zero();
    prob->edit();
    prob->finalize(true,true,*outStream);
    parlist->sublist("Step").set("Type","Augmented Lagrangian");
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Step Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","Lin-More");
    Teuchos::Time allmtrTimer("AL-LMTR Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    allmtrTimer.stop();
    *outStream << "Total optimization time = " << allmtrTimer.totalElapsedTime() << " seconds.\n";
    obj->summarize(*outStream);
    vals[4] = solver->getAlgorithmState()->value;

    // AL-LMTR
    obj->reset();
    x->zero();
    mul->zero();
    parlist->sublist("Step").set("Type","Augmented Lagrangian");
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Step Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","SPG");
    Teuchos::Time altrspgTimer("AL-TRSPG Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    altrspgTimer.stop();
    *outStream << "Total optimization time = " << altrspgTimer.totalElapsedTime() << " seconds.\n";

    int cnt(0);
    for (int i = 0; i < nInstances; ++i) {
      if ((*x->getVector())[i] >= static_cast<RealT>(1e-4)) cnt++;
    }
    *outStream << std::endl << "Number of Support Vectors: " << cnt << std::endl << std::endl;
    obj->summarize(*outStream);
    vals[5] = solver->getAlgorithmState()->value;

    RealT minval = std::min(vals[0],std::min(vals[1],std::min(vals[2],vals[3])));
    *outStream << std::scientific << std::setprecision(16);
    *outStream << "SPG      " << (vals[0]-minval)/minval << std::endl;
    *outStream << "PQN      " << (vals[1]-minval)/minval << std::endl;
    *outStream << "LMTR     " << (vals[2]-minval)/minval << std::endl;
    *outStream << "TRSPG    " << (vals[3]-minval)/minval << std::endl;
    *outStream << "AL-LMTR  " << (vals[4]-minval)/minval << std::endl;
    *outStream << "AL-TRSPG " << (vals[5]-minval)/minval << std::endl;
    *outStream << std::endl;

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


