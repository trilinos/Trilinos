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
#include "ROL_LinearAlgebra.hpp"
#include "ROL_BLAS.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <iostream>
#include <fstream>
#include <random>

typedef double RealT;

template<class Real>
class LassoObjective : public ROL::StdObjective<Real> {
private:
  const ROL::Ptr<ROL::BLAS<int,Real>> blas_;
  const ROL::LA::Matrix<Real> X_;
  const std::vector<Real> y_;
  std::vector<Real> res_;
  const int nrows_, ncols_;
  const Real wt_;
  int fval_, grad_, hess_;
  bool computeRes_;

  void computeResidual(const std::vector<Real> &x) {
    if (computeRes_) {
      const Real one(1);
      std::vector<Real> d(ncols_);
      for (int j = 0; j < ncols_; ++j) d[j] = (x[j]-x[ncols_+j]);
      res_.assign(y_.begin(),y_.end());
      applyK(res_,d,one,-one,false);
      computeRes_ = false;
    }
  }

  void applyK(std::vector<Real> &Kx, const std::vector<Real> &x, Real alpha, Real beta, bool trans) const {
    //const Real zero(0), one(1);
    if (trans)
      blas_->GEMV(   Teuchos::TRANS,nrows_,ncols_,alpha,X_.values(),nrows_,&x[0],1,beta,&Kx[0],1);
    else
      blas_->GEMV(Teuchos::NO_TRANS,nrows_,ncols_,alpha,X_.values(),nrows_,&x[0],1,beta,&Kx[0],1);
  }

public:
  LassoObjective(const ROL::LA::Matrix<Real> &X, const std::vector<Real> &y)
    : blas_(ROL::makePtr<ROL::BLAS<int,Real>>()),
      X_(X), y_(y), nrows_(X.numRows()), ncols_(X.numCols()),
      wt_(static_cast<Real>(1)/static_cast<Real>(nrows_)),
      fval_(0), grad_(0), hess_(0), computeRes_(true) {
    res_.resize(nrows_);
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

  void update(const std::vector<Real> &x, ROL::UpdateType type, int iter = -1) {
    if (type != ROL::UpdateType::Revert && type != ROL::UpdateType::Accept)
      computeRes_ = true;
  }

  Real value(const std::vector<Real> &x, Real &tol) {
    fval_++;
    computeResidual(x);
    return wt_*static_cast<Real>(0.5)*blas_->DOT(nrows_,&res_[0],1,&res_[0],1);
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    grad_++;
    const Real zero(0), one(1);
    computeResidual(x);
    std::vector<Real> Kr(ncols_);
    applyK(Kr,res_,one,zero,true);
    for (int j = 0; j < ncols_; ++j) {
      g[j]        =  wt_*Kr[j];
      g[ncols_+j] = -wt_*Kr[j];
    }
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    hess_++;
    const Real zero(0), one(1);
    std::vector<Real> d(ncols_), res(nrows_), Kr(ncols_);
    for (int j = 0; j < ncols_; ++j) d[j] = (v[j]-v[ncols_+j]);
    res_.assign(y_.begin(),y_.end());
    applyK(res,d,one,zero,false);
    applyK(Kr,res,one,zero,true);
    for (int j = 0; j < ncols_; ++j) {
      hv[j]        =  wt_*Kr[j];
      hv[ncols_+j] = -wt_*Kr[j];
    }
  }

};

template<typename Real>
class LassoConstraint : public ROL::StdConstraint<Real> {
private:
  const Real t_;
  const std::vector<Real> d_;
  const int size_;

public:
  LassoConstraint(Real t, const std::vector<Real> &d) : t_(t), d_(d), size_(d.size()) {}

  void value(std::vector<Real> &c, const std::vector<Real> &x, Real &tol) {
    c[0] = -t_;
    for (int j = 0; j < size_; ++j) c[0] += (x[j]+x[size_+j])/d_[j];
  }

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    jv[0] = static_cast<Real>(0);
    for (int j = 0; j < size_; ++j) jv[0] += (v[j]+v[size_+j])/d_[j];
  }

  void applyAdjointJacobian(std::vector<Real> &ajv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    for (int j = 0; j < size_; ++j) {
      ajv[j]       = v[0]/d_[j];
      ajv[size_+j] = v[0]/d_[j];
    }
  }

  void applyAdjointHessian(std::vector<Real> &ahuv, const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    ahuv.assign(x.size(),static_cast<Real>(0));
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
    bool useRandomData = false;
    int nfactors, nVars;
    ROL::LA::Matrix<RealT> X;
    std::vector<RealT> y, d;
    if (useRandomData) {
      unsigned seed = 111;
      nfactors = 500000;
      nVars    = 500;
      X.reshape(nfactors,nVars);
      y.resize(nfactors,0.0);
      d.resize(nVars,1.0);
      std::default_random_engine generator(seed);
      std::normal_distribution<RealT> dist(0.0,1.0);
      for (int i = 0; i < nfactors; ++i) {
        for (int j = 0; j < nVars; ++j) X(i,j) = dist(generator);
      }
      std::vector<RealT> w(nVars);
      for (int j = 0; j < nVars; ++j) {
        RealT u = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
        w[j] = (u > static_cast<RealT>(0.5) ? dist(generator) : static_cast<RealT>(0));
      }
      for (int i = 0; i < nfactors; ++i) {
        y[i] = dist(generator);
        for (int j = 0; j < nVars; ++j) y[i] += X(i,j)*w[j];
      }
    }
    else {
      nfactors = 20640;
      nVars    = 9;
      X.reshape(nfactors,nVars);
      y.resize(nfactors,0.0);
      d.resize(nVars,0.0);
      std::ifstream file;
      file.open("lasso-data.txt");
      RealT medVal(0), medInc(0), medAge(0), totRoom(0), totBed(0), pop(0), house(0), lon(0), lat(0);
      for (int i = 0; i < nfactors; ++i) {
        file >> medVal;
        file >> medInc;
        file >> medAge;
        file >> totRoom;
        file >> totBed;
        file >> pop;
        file >> house;
        file >> lon;
        file >> lat;
        y[i] = std::log(medVal);
        X(i,0) = static_cast<RealT>(1);
        X(i,1) = medInc;
        X(i,2) = std::pow(medInc,2);
        X(i,3) = std::pow(medInc,3);
        X(i,4) = std::log(medAge);
        X(i,5) = std::log(totRoom)-std::log(pop);
        X(i,6) = std::log(totBed)-std::log(pop);
        X(i,7) = std::log(pop)-std::log(house);
        X(i,8) = std::log(house);
        for (int j = 0; j < nVars; ++j) d[j] = std::max(std::abs(X(i,j)),d[j]);
      }
      file.close();
      for (int i = 0; i < nfactors; ++i) {
        for (int j = 0; j < nVars; ++j) X(i,j) /= d[j];
      }
    }

    RealT alpha(1);
    ROL::Ptr<LassoObjective<RealT>>  obj = ROL::makePtr<LassoObjective<RealT>>(X,y);
    ROL::Ptr<LassoConstraint<RealT>> con = ROL::makePtr<LassoConstraint<RealT>>(alpha,d);
    ROL::Ptr<ROL::StdVector<RealT>>  mul = ROL::makePtr<ROL::StdVector<RealT>>(1,0);
    ROL::Ptr<ROL::StdVector<RealT>>  res = ROL::makePtr<ROL::StdVector<RealT>>(1,0);
    ROL::Ptr<ROL::StdVector<RealT>>    u = ROL::makePtr<ROL::StdVector<RealT>>(1,0);
    ROL::Ptr<ROL::Bounds<RealT>>    ibnd = ROL::makePtr<ROL::Bounds<RealT>>(*u,false);
    ROL::Ptr<ROL::StdVector<RealT>>    x = ROL::makePtr<ROL::StdVector<RealT>>(2*nVars,0);
    ROL::Ptr<ROL::StdVector<RealT>>    l = ROL::makePtr<ROL::StdVector<RealT>>(2*nVars,0);
    ROL::Ptr<ROL::Bounds<RealT>>     bnd = ROL::makePtr<ROL::Bounds<RealT>>(*l,true);
    ROL::Ptr<ROL::Problem<RealT>>   prob = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    prob->addBoundConstraint(bnd);
    prob->addLinearConstraint("L1",con,mul,ibnd,res);
    prob->setProjectionAlgorithm(*parlist);
    prob->finalize(false,true,*outStream);
    prob->check(true,*outStream);

    ROL::Ptr<ROL::Solver<RealT>> solver;
    std::vector<RealT> vals(6);

    // SPG
    obj->reset();
    x->zero();
    prob->edit();
    prob->finalize(false,true,*outStream);
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
    prob->edit();
    prob->finalize(false,true,*outStream);
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
    prob->edit();
    prob->finalize(false,true,*outStream);
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
    prob->edit();
    prob->finalize(false,true,*outStream);
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

    // AL-TRSPG
    obj->reset();
    x->zero();
    mul->zero();
    prob->edit();
    prob->finalize(true,true,*outStream);
    parlist->sublist("Step").set("Type","Augmented Lagrangian");
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Step Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","SPG");
    Teuchos::Time altrspgTimer("AL-TRSPG Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    altrspgTimer.stop();
    *outStream << "Total optimization time = " << altrspgTimer.totalElapsedTime() << " seconds.\n";
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


