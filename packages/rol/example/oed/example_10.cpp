// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test StatusTest input mechanism in OptimizationSolver.
*/

#include "ROL_StdConstraint.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_UserInputGenerator.hpp"
#include "ROL_Solver.hpp"
#include "ROL_Stream.hpp"

#include "ROL_OED_Factory.hpp"
#include "ROL_OED_StdMomentOperator.hpp"
#include "ROL_OED_GreedyAlgorithm.hpp"
#include "ROL_OED_PrintDesign.hpp"

#include "ROL_GlobalMPISession.hpp"

#include <iostream>

template<typename Real>
class HelmholtzModel : public ROL::StdConstraint<Real> {
private:
  const unsigned nx_;
  const Real dx_;
  ROL::LA::Matrix<Real> data_;

  void generateData(Real k, Real c) {
    const unsigned N(2*nx_);
    const Real sd0(1.0/dx_), sd(2.0/dx_), so(-1.0/dx_);    // Stiffness matrix
    const Real md0(dx_/6.0), md(dx_*2.0/3.0), mo(dx_/6.0); // Mass matrix
    const Real k2(k*k), kc(k*c), one(1);
    ROL::LA::Matrix<Real> K(N,N);
    data_.reshape(N,4);
    data_.putScalar(0.0);
    for (unsigned i=0u; i<nx_; ++i) {
      if (i==0u) {
        K(i,i)           = sd0-k2*md0;
        K(i,i+1)         = so-k2*mo;
        K(i,i+nx_)       = kc*md0;
        K(i,i+nx_+1)     = kc*mo;
        K(i+nx_,i)       = -K(i,i+nx_);
        K(i+nx_,i+1)     = -K(i,i+nx_+1);
        K(i+nx_,i+nx_)   = K(i,i);
        K(i+nx_,i+nx_+1) = K(i,i+1);
      }
      else if (i==nx_-1u) {
        K(i,i)           = sd0-k2*md0;
        K(i,i-1)         = so-k2*mo;
        K(i,i+nx_)       = kc*md0;
        K(i,i+nx_-1)     = kc*mo;
        K(i+nx_,i)       = -K(i,i+nx_);
        K(i+nx_,i-1)     = -K(i,i+nx_-1);
        K(i+nx_,i+nx_)   = K(i,i);
        K(i+nx_,i+nx_-1) = K(i,i-1);
      }
      else {
        K(i,i)           = sd-k2*md;
        K(i,i-1)         = so-k2*mo;
        K(i,i+1)         = K(i,i-1);
        K(i,i+nx_)       = kc*md;
        K(i,i+nx_-1)     = kc*mo;
        K(i,i+nx_+1)     = K(i,i+nx_-1);
        K(i+nx_,i)       = -K(i,i+nx_);
        K(i+nx_,i-1)     = -K(i,i+nx_-1);
        K(i+nx_,i+1)     = -K(i,i+nx_+1);
        K(i+nx_,i+nx_)   = K(i,i);
        K(i+nx_,i+nx_-1) = K(i,i-1);
        K(i+nx_,i+nx_+1) = K(i,i+1);
      }
    }
    //K.print(std::cout);
    data_(0,0) = one;
    data_(nx_-1,1) = one;
    data_(nx_,2) = one;
    data_(2*nx_-1,3) = one;
    ROL::LAPACK<int,Real> lapack;
    int info;
    std::vector<int> ipiv(N);
    lapack.GETRF(N,N,K.values(),N,&ipiv[0],&info);
    lapack.GETRS('N',N,4,K.values(),N,&ipiv[0],data_.values(),N,&info);
  }

  void processParameter(unsigned& ind, Real& frac0, Real& frac1, Real x) const {
    const Real tol(std::sqrt(ROL::ROL_EPSILON<Real>())*dx_);
    if (std::abs(x-static_cast<Real>(0))<tol) {
      ind = 0u;
      frac0 = static_cast<Real>(1);
      frac1 = static_cast<Real>(0);
    }
    else if (std::abs(x-static_cast<Real>(1))<tol) {
      ind = nx_-2u;
      frac0 = static_cast<Real>(0);
      frac1 = static_cast<Real>(1);
    }
    else {
      const Real xi(x/dx_);
      ind = static_cast<unsigned>(std::floor(xi));
      frac0 = (xi-static_cast<Real>(ind));
      frac1 = (static_cast<Real>(ind+1u)-xi);
    }
  }

public:
  HelmholtzModel(unsigned nx = 65, Real k = Real(2.0*M_PI), Real c = Real(0.2))
    : nx_(nx), dx_(1.0/(static_cast<Real>(nx)-1.0)) {
    generateData(k,c);
  }

  void printState(std::string filename, const std::vector<Real>& theta, unsigned nx=65) {
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    const Real dx(1.0/(static_cast<Real>(nx)-1));
    std::vector<Real> u(2);
    Real x;
    std::ofstream file;
    file.open(filename);
    file << std::scientific << std::setprecision(15);
    for (unsigned i=0u; i<nx; ++i) {
      x = static_cast<Real>(i)*dx;
      ROL::StdConstraint<Real>::setParameter({x});
      value(u,theta,tol);
      file << std::right << std::setw(25) << x;
      file << std::right << std::setw(25) << u[0];
      file << std::right << std::setw(25) << u[1];
      file << std::endl;
    }
    file.close();
  }

  void value(std::vector<Real>& c, const std::vector<Real> &theta, Real &tol) override {
    //const unsigned ind(ROL::Constraint<Real>::getParameter()[0]);
    const Real x(ROL::Constraint<Real>::getParameter()[0]);
    unsigned ind;
    Real frac0, frac1;
    processParameter(ind,frac0,frac1,x);
    c[0] = static_cast<Real>(0);
    c[1] = static_cast<Real>(0);
    for (unsigned i=0u; i<4u; ++i) {
      c[0] += (frac0*data_(    ind,i)+frac1*data_(    ind+1,i))*theta[i];
      c[1] += (frac0*data_(nx_+ind,i)+frac1*data_(nx_+ind+1,i))*theta[i];
    }
  }

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real>& v, const std::vector<Real> &theta, Real &tol) override {
    //const unsigned ind(ROL::Constraint<Real>::getParameter()[0]);
    const Real x(ROL::Constraint<Real>::getParameter()[0]);
    unsigned ind;
    Real frac0, frac1;
    processParameter(ind,frac0,frac1,x);
    jv[0] = static_cast<Real>(0);
    jv[1] = static_cast<Real>(0);
    for (unsigned i=0u; i<4u; ++i) {
      jv[0] += (frac0*data_(    ind,i)+frac1*data_(    ind+1,i))*v[i];
      jv[1] += (frac0*data_(nx_+ind,i)+frac1*data_(nx_+ind+1,i))*v[i];
    }
  }

  void applyAdjointJacobian(std::vector<Real> &ajv, const std::vector<Real>& v, const std::vector<Real> &theta, Real &tol) override {
    //const unsigned ind(ROL::Constraint<Real>::getParameter()[0]);
    const Real x(ROL::Constraint<Real>::getParameter()[0]);
    unsigned ind;
    Real frac0, frac1;
    processParameter(ind,frac0,frac1,x);
    for (unsigned i=0u; i<4u; ++i)
      ajv[i] = (frac0*data_(    ind,i)+frac1*data_(    ind+1,i))*v[0]
              +(frac0*data_(nx_+ind,i)+frac1*data_(nx_+ind+1,i))*v[1];
  }
};

template<typename Real>
class HelmholtzNoise : public ROL::OED::Noise<Real> {
private:
  const Real alpha_;

public:
  HelmholtzNoise(Real alpha = Real(2)) : alpha_(alpha) {}

  bool isHomoscedastic() const override { return false; }
  bool isCorrelated() const override { return true; }

  Real evaluate(const std::vector<Real> &x) const override {
    return std::exp(-alpha_ * std::abs(x[0]-static_cast<Real>(0.5))) / alpha_;
  }
  void apply(ROL::Vector<Real>& Sx, const ROL::Vector<Real>& x, const std::vector<Real>& param) const override {
    auto& Sxd = *static_cast<ROL::StdVector<Real>&>(Sx).getVector();
    auto& xd = *static_cast<const ROL::StdVector<Real>&>(x).getVector();
    const Real s = evaluate(param), d(2.0/6.0), o(1.0/6.0);
    Sxd[0] = s*s*(d * xd[0] + o * xd[1]);
    Sxd[1] = s*s*(o * xd[0] + d * xd[1]);
  }
  void applyInverse(ROL::Vector<Real>& Sx, const ROL::Vector<Real>& x, const std::vector<Real>& param) const override {
    auto& Sxd = *static_cast<ROL::StdVector<Real>&>(Sx).getVector();
    auto& xd = *static_cast<const ROL::StdVector<Real>&>(x).getVector();
    const Real s = evaluate(param), d0(2.0/6.0), o0(1.0/6.0);
    const Real det(d0*d0-o0*o0), d(d0/det), o(-o0/det);
    Sxd[0] = (d * xd[0] + o * xd[1])/(s*s);
    Sxd[1] = (o * xd[0] + d * xd[1])/(s*s);
  }
};

template<typename Real>
class RegularizationOperator : public ROL::LinearOperator<Real> {
private:
  const Real alpha_;

public:
  RegularizationOperator(Real alpha = Real(1)) : alpha_(alpha) {}

  void apply(ROL::Vector<Real> &Px, const ROL::Vector<Real> &x, Real &tol) const override {
    Px.set(x);
    Px.scale(alpha_);
  }

  void applyInverse(ROL::Vector<Real> &Px, const ROL::Vector<Real> &x, Real &tol) const override {
    Px.set(x);
    Px.scale(static_cast<Real>(1)/alpha_);
  }

  void applyAdjoint(ROL::Vector<Real> &Px, const ROL::Vector<Real> &x, Real &tol) const override {
    apply(Px,x,tol);
  }

  void applyAdjointInverse(ROL::Vector<Real> &Px, const ROL::Vector<Real> &x, Real &tol) const override {
    applyInverse(Px,x,tol);
  }
};

typedef double RealT;

int main(int argc, char *argv[]) {

  ROL::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Test body.

  try {
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );

    // Setup parameter vector and polynomial model
    const unsigned nx = static_cast<unsigned>(parlist->sublist("Problem").get("Number of DoFs", 65));
    const RealT alpha = parlist->sublist("Problem").get("Noise Decay Rate", 5.0);
    auto theta = ROL::makePtr<ROL::StdVector<RealT>>(4,1);
    auto obs   = ROL::makePtr<ROL::StdVector<RealT>>(2,0);
    auto model = ROL::makePtr<HelmholtzModel<RealT>>(nx);
    auto noise = ROL::makePtr<HelmholtzNoise<RealT>>(alpha);

    // Setup experiment sample generator
    const int nsamp = parlist->sublist("Problem").get("Number of Samples", 100);
    std::ofstream ptfile, wtfile;
    ptfile.open("points.txt");
    wtfile.open("weights.txt");
    for (int i = 0; i < nsamp; ++i) {
      ptfile << static_cast<RealT>(i+1) / static_cast<RealT>(nsamp+1) << std::endl;
      wtfile << static_cast<RealT>(1) / static_cast<RealT>(nsamp) << std::endl;
    }
    ptfile.close();
    wtfile.close();
    auto bman = ROL::makePtr<ROL::BatchManager<RealT>>();
    auto sampler = ROL::makePtr<ROL::UserInputGenerator<RealT>>("points.txt","weights.txt",nsamp,1,bman);

    // Setup factory
    bool homNoise = true;
    std::string regType = "Least Squares";
    std::string ocType = parlist->sublist("OED").get("Optimality Type","A");
    auto type = ROL::OED::StringToRegressionType(regType);
    auto M = ROL::makePtr<ROL::OED::StdMomentOperator<RealT>>(type,(homNoise ? ROL::nullPtr : noise));
    bool addTik = parlist->sublist("Problem").get("Use Tikhonov",false);
    if (addTik) {
      RealT beta  = parlist->sublist("Problem").get("Tikhonov Parameter",1e-4);
      auto P = ROL::makePtr<RegularizationOperator<RealT>>(beta);
      M->setPerturbation(P);
    }
    auto factory = ROL::makePtr<ROL::OED::Factory<RealT>>(model,sampler,theta,obs,M,*parlist);
    if (parlist->sublist("Problem").get("Use Budget Constraint",false)) {
      auto cost = factory->createDesignVector();
      cost->setScalar(static_cast<RealT>(1));
      RealT budget = parlist->sublist("Problem").get("Budget",5.0);
      bool useBudgetEquality = parlist->sublist("Problem").get("Use Budget Equality Constraint",false);
      factory->setBudgetConstraint(cost,budget,useBudgetEquality);
    }
    if (parlist->sublist("Problem").get("Use Probability Scaling",false)) {
      auto prob = factory->createDesignVector();
      auto prob0 = ROL::dynamicPtrCast<ROL::StdVector<RealT>>(prob)->getVector();
      for (unsigned i = 0u; i < prob0->size(); ++i)
        (*prob0)[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      factory->setProbabilityVector(prob);
    }
    //if (ocType == "A" || ocType == "I")
    //  parlist->sublist("General").sublist("Polyhedral Projection").set("Type","Brents");
    //else
    //  parlist->sublist("General").sublist("Polyhedral Projection").set("Type","Dai-Fletcher");

    auto x = theta->clone(); x->randomize(-1.0,1.0);
    auto v = theta->clone(); v->randomize(-1.0,1.0);
    auto jv = obs->clone(); obs->randomize(-1.0,1.0);
    auto w = obs->dual().clone(); w->randomize(-1.0,1.0);
    model->setParameter(sampler->getMyPoint(0));
    model->checkApplyJacobian(*x,*v,*jv,true,*outStream);
    model->checkAdjointConsistencyJacobian(*w,*v,*x,true,*outStream);
    
    // Generate optimization problem
    auto problem = factory->get(*parlist,sampler);
    problem->setProjectionAlgorithm(*parlist);
    problem->finalize(false,true,*outStream);
    auto test = factory->getDesign()->clone();
    test->randomize(1,2);
    problem->check(true,*outStream,test,0.1);

    auto obj  = problem->getObjective();
    auto x0   = factory->getDesign()->clone();
    auto v0   = factory->getDesign()->clone();
    auto xtmp = factory->getDesign()->clone();
    auto g    = xtmp->dual().clone();
    x0->randomize(0.0,2.0);
    v0->randomize(-1.0,1.0);
    RealT tol(std::sqrt(ROL::ROL_EPSILON<RealT>()));
    RealT t(0.1), a(0.125), ia(8), c(ia), tmp, dd;
    RealT err(ROL::ROL_INF<RealT>()), err0(err), minerr(err);
    const unsigned Nmax(10);
    unsigned N(1);
    std::vector<RealT> vals(Nmax);
    obj->update(*x0,ROL::UpdateType::Temp);
    RealT val = obj->value(*x0,tol);
    obj->gradient(*g,*x0,tol);
    RealT dd0 = g->apply(*v0);
    xtmp->set(*x0); xtmp->axpy(t,*v0);
    obj->update(*xtmp,ROL::UpdateType::Temp);
    vals[0] = (obj->value(*xtmp,tol)-val)/t;
    dd = vals[0];
    const RealT atol(0e-4), rtol(0e-2), two(2);
    *outStream << std::endl;
    *outStream << std::right
               << std::setw(20) << "Step size"
               << std::setw(20) << "grad'*dir"
               << std::setw(20) << "FD approx"
               << std::setw(20) << "abs error"
               << std::endl
               << std::setw(20) << "---------"
               << std::setw(20) << "---------"
               << std::setw(20) << "---------"
               << std::setw(20) << "---------"
               << std::endl;
    *outStream << std::scientific << std::setprecision(11) << std::right
               << std::setw(20) << t
               << std::setw(20) << dd0
               << std::setw(20) << vals[0]
               << std::setw(20) << std::abs(vals[0]-dd0)
               << std::endl;
    *outStream << std::scientific << std::setprecision(11) << std::right
               << std::setw(20) << "Extrapolation"
               << std::setw(20) << dd0
               << std::setw(20) << dd
               << std::setw(20) << std::abs(dd-dd0)
               << std::endl;
    while(N <= Nmax) {
      t *= a;
      c = ia;
      xtmp->set(*x0); xtmp->axpy(t,*v0);
      obj->update(*xtmp,ROL::UpdateType::Temp);
      vals[N] = (obj->value(*xtmp,tol)-val)/t;
      *outStream << std::scientific << std::setprecision(11) << std::right
                 << std::setw(20) << t
                 << std::setw(20) << dd0
                 << std::setw(20) << vals[N]
                 << std::setw(20) << std::abs(vals[N]-dd0)
                 << std::endl;
      N++;
      minerr = ROL::ROL_INF<RealT>();
      for(unsigned i=N-1u; i>0u; --i) {
        tmp = vals[i-1u];
        vals[i-1u] = vals[i] + (vals[i]-tmp)/(c-1);
	err0 = std::abs(vals[i-1u]-tmp);
	minerr = std::min(minerr,err0);
	if (err0 < err) {
          dd = vals[i-1u];
	  err = err0;
	}
	c *= ia;
      }
      *outStream << std::scientific << std::setprecision(11) << std::right
                 << std::setw(20) << "Extrapolation"
                 << std::setw(20) << dd0
                 << std::setw(20) << dd
                 << std::setw(20) << std::abs(dd-dd0)
                 << std::endl;
      if (minerr > two*err || err <= std::min(atol,rtol*std::abs(dd))) break;
    }

  }
  catch (std::logic_error& err) {
    *outStream << err.what() << std::endl;
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED" << std::endl;
  else
    std::cout << "End Result: TEST PASSED" << std::endl;

  return 0;

}
