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

#include "ROL_StdObjective.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_UserInputGenerator.hpp"
#include "ROL_Stream.hpp"

#include "ROL_OED_BinaryDesignAlgorithm.hpp"
#include "ROL_OED_StdMomentOperator.hpp"

#include "ROL_GlobalMPISession.hpp"

#include <iostream>

template<typename Real>
class PolynomialModel : public ROL::StdObjective<Real> {
private:
  const unsigned deg_;

public:
  PolynomialModel(unsigned deg = 2) : deg_(deg) {}

  Real value(const std::vector<Real> &theta, Real &tol) override {
    Real val(0), xpow(1);
    for (unsigned i = 0u; i <= deg_; ++i) {
      val  += theta[i] * xpow;
      xpow *= ROL::Objective<Real>::getParameter()[0];
    }
    return val;
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &theta, Real &tol) override {
    Real xpow(1);
    for (unsigned i = 0u; i <= deg_; ++i) {
      g[i]  = xpow;
      xpow *= ROL::Objective<Real>::getParameter()[0];
    }
  }
};

template<typename Real>
class PolynomialNoise : public ROL::OED::Noise<Real> {
private:
  const Real alpha_;

public:
  PolynomialNoise(Real alpha = Real(1)) : alpha_(alpha) {}

  Real evaluate(const std::vector<Real> &x) const override {
    return std::exp(alpha_ * std::abs(x[0])) / std::exp(alpha_);
  }

  bool isHomoscedastic() const override { return false; }
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
    std::string filename = "input_11.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );

    // Setup parameter vector and polynomial model
    const RealT alpha = parlist->sublist("Problem").get("Noise Decay Rate", 5.0);
    const int deg = parlist->sublist("Problem").get("Polynomial Degree", 5);
    auto theta = ROL::makePtr<ROL::StdVector<RealT>>(deg+1,1);
    auto model = ROL::makePtr<PolynomialModel<RealT>>(deg);
    auto noise = ROL::makePtr<PolynomialNoise<RealT>>(alpha);

    // Setup experiment sample generator
    const RealT lb = parlist->sublist("Problem").get("X Lower Bound", -1.0);
    const RealT ub = parlist->sublist("Problem").get("X Upper Bound",  1.0);
    const int nsamp = parlist->sublist("Problem").get("Number of Samples", 100);
    std::ofstream ptfile, wtfile;
    ptfile.open("points.txt");
    wtfile.open("weights.txt");
    for (int i = 0; i < nsamp; ++i) {
      ptfile << std::scientific << std::setprecision(16);
      wtfile << std::scientific << std::setprecision(16);
      ptfile << (ub - lb) * static_cast<RealT>(i) / static_cast<RealT>(nsamp - 1) + lb << std::endl;
      wtfile << static_cast<RealT>(1) / static_cast<RealT>(nsamp)                      << std::endl;
    }
    ptfile.close();
    wtfile.close();
    auto bman = ROL::makePtr<ROL::BatchManager<RealT>>();
    auto sampler = ROL::makePtr<ROL::UserInputGenerator<RealT>>("points.txt","weights.txt",nsamp,1,bman);

    // Setup binary design algorithm
    bool homNoise = true;
    std::string regType = "Least Squares";
    std::string ocType = parlist->sublist("OED").get("Optimality Type","A");
    auto type = ROL::OED::StringToRegressionType(regType);
    auto M = ROL::makePtr<ROL::OED::StdMomentOperator<RealT>>(type,(homNoise ? ROL::nullPtr : noise));
    bool addTik = parlist->sublist("Problem").get("Use Tikhonov",false);
    if (addTik) {
      const RealT beta = parlist->sublist("Problem").get("Tikhonov Parameter",1e-4);
      auto P = ROL::makePtr<RegularizationOperator<RealT>>(beta);
      M->setPerturbation(P);
    }
    bool useBudgetEquality = parlist->sublist("Problem").get("Use Budget Equality Constraint",true);
    auto bdalg = ROL::makePtr<ROL::OED::BinaryDesignAlgorithm<RealT>>(model,sampler,theta,M,*parlist);
    auto cost  = bdalg->createCostVector();
    cost->setScalar(static_cast<RealT>(1));
    RealT budget = parlist->sublist("Problem").get("Budget",5.0);
    bdalg->setBudgetConstraint(cost,budget,useBudgetEquality);
    bdalg->run(*parlist,sampler,*outStream);
    
    //auto factory = ROL::makePtr<ROL::OED::Factory<RealT>>(model,sampler,theta,M,*parlist);
    //auto cost = factory->createDesignVector();
    //cost->setScalar(static_cast<RealT>(1));
    //RealT budget = parlist->sublist("Problem").get("Budget",5.0);
    //factory->setBudgetConstraint(cost,budget,useBudgetEquality);
    //// Generate optimization problem
    //auto problem0 = factory->get(*parlist,sampler);
    //problem0->setProjectionAlgorithm(*parlist);
    //problem0->finalize(false,true,*outStream);
    //auto test = factory->getDesign()->clone();
    //test->randomize(1,2);
    //problem0->check(true,*outStream,test,0.1);

    //// Get problem components
    //auto crit  = problem0->getObjective();
    //auto x     = problem0->getPrimalOptimizationVector();
    //auto proj  = problem0->getPolyhedralProjection();
    //auto bnd   = proj->getBoundConstraint();
    //auto lv    = ROL::makePtr<ROL::SingletonVector<RealT>>(); lv->zero();
    //auto lbnd  = ROL::makePtr<ROL::Bounds<RealT>>(*lv,false);
    //auto lcon  = ROL::makePtr<ROL::ScalarLinearConstraint<RealT>>(cost,budget);
    //auto lmul  = ROL::makePtr<ROL::SingletonVector<RealT>>();
    //auto dwpen = ROL::makePtr<ROL::OED::DoubleWellPenalty<RealT>>(1u);

    //std::vector<ROL::Ptr<ROL::Objective<RealT>>> ovec = {crit,dwpen};
    //std::vector<RealT> wvec = {1.0,0.0};

    //auto problem = ROL::makePtr<ROL::Problem<RealT>>(crit,x);
    //problem->addBoundConstraint(bnd);
    //problem->addLinearConstraint("Budget",lcon,lmul,lbnd);

    //// Setup ROL solver
    //std::clock_t timer = std::clock();
    //auto solver = ROL::makePtr<ROL::Solver<RealT>>(problem,*parlist);
    //solver->solve(*outStream);
    //*outStream << "  " << ocType << "-optimal design time:      "
    //           << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
    //           << " seconds" << std::endl;
    //factory->profile(*outStream);

    //RealT pen     = parlist->sublist("Problem").get("Initial Penalty Parameter",1.0);
    //RealT penscal = parlist->sublist("Problem").get("Penalty Parameter Update Scaling",10.0);
    //int contype   = parlist->sublist("Problem").get("Continuation Type",1);

    //RealT tau(1.0/penscal), A, B, C, c, p, v, v0, vtmp;
    //RealT tol(1e-2*std::sqrt(ROL::ROL_EPSILON<RealT>())), etol(tol);

    //dwpen->update(*x,ROL::UpdateType::Trial);
    //p = dwpen->value(*x,tol);
    //crit->update(*x,ROL::UpdateType::Trial);
    //c = crit->value(*x,tol);
    //v0 = c;

    //bool print = false;
    //if (print) {
    //  RealT logAlphaMax = 2.0;
    //  unsigned nalpha = 20;
    //  std::ofstream file;
    //  file.open("value_function_ex11.txt");
    //  file << std::scientific << std::setprecision(15);
    //  for (unsigned i=0u; i<nalpha; ++i) {
    //    RealT alpha = std::pow(static_cast<RealT>(10),static_cast<RealT>(i)*logAlphaMax/static_cast<RealT>(nalpha-1));
    //    // Generate optimization problem
    //    wvec[1] = alpha;
    //    auto obj = ROL::makePtr<ROL::LinearCombinationObjective<RealT>>(wvec,ovec);
    //    problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    //    problem->addBoundConstraint(bnd);
    //    problem->addLinearConstraint("Budget",lcon,lmul,lbnd);
    //    problem->setProjectionAlgorithm(*parlist);
    //    problem->finalize(false,true,*outStream);
    //    //test = factory->getDesign()->clone();
    //    //test->randomize(1,2);
    //    //problem->check(true,*outStream,test,0.1);

    //    // Setup ROL solver
    //    timer = std::clock();
    //    solver = ROL::makePtr<ROL::Solver<RealT>>(problem,*parlist);
    //    solver->solve(*outStream);
    //    *outStream << "  " << ocType << "-optimal design time:      "
    //               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
    //               << " seconds" << std::endl;
    //    factory->profile(*outStream);

    //    dwpen->update(*x,ROL::UpdateType::Trial);
    //    crit->update(*x,ROL::UpdateType::Trial);
    //    file << std::right << std::setw(25) << alpha;
    //    file << std::right << std::setw(25) << crit->value(*x,tol);
    //    file << std::right << std::setw(25) << dwpen->value(*x,tol);
    //    file << std::endl;
    //  }
    //  file.close();
    //  return 0;
    //}

    //if (contype==1) {
    //  auto xb  = x->clone();
    //  auto rnd = ROL::Elementwise::Round<RealT>();
    //  xb->set(*x); xb->applyUnary(rnd);
    //  crit->update(*xb,ROL::UpdateType::Trial);
    //  auto cb = crit->value(*xb,tol);
    //  pen = std::max(pen,(cb-c)/p);
    //  *outStream << "Initial Cost:  " << c << std::endl;
    //  *outStream << "Binary Cost:   " << cb << std::endl;
    //  *outStream << std::endl;
    //}

    //const RealT two(2), three(3);
    //for (unsigned i=0u; i<10u; ++i) {
    //  *outStream << "Binary Metric: " << p << std::endl;
    //  *outStream << "Penalty:       " << pen << std::endl;
    //  // Generate optimization problem
    //  wvec[1] = pen;
    //  auto obj = ROL::makePtr<ROL::LinearCombinationObjective<RealT>>(wvec,ovec);
    //  problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    //  problem->addBoundConstraint(bnd);
    //  problem->addLinearConstraint("Budget",lcon,lmul,lbnd);
    //  problem->setProjectionAlgorithm(*parlist);
    //  problem->finalize(false,true,*outStream);
    //  test = factory->getDesign()->clone();
    //  test->randomize(1,2);
    //  problem->check(true,*outStream,test,0.1);

    //  // Setup ROL solver
    //  timer = std::clock();
    //  solver = ROL::makePtr<ROL::Solver<RealT>>(problem,*parlist);
    //  solver->solve(*outStream);
    //  *outStream << "  " << ocType << "-optimal design time:      "
    //             << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
    //             << " seconds" << std::endl;
    //  factory->profile(*outStream);

    //  // Check for binary solution
    //  dwpen->update(*x,ROL::UpdateType::Trial);
    //  p = dwpen->value(*x,tol);
    //  
    //  if (p < etol) break;

    //  if (contype==1) {
    //    // Exact path following
    //    crit->update(*x,ROL::UpdateType::Trial);
    //    c = crit->value(*x,tol);
    //    v = c + pen*p;
    //    C = pen*pen*p/(c-v0);
    //    B = (C*(C+pen)*(v-v0))/pen;
    //    A = v0 + B/C;
    //    pen = B/(tau * std::abs(A-v))-C;
    //    tau /= penscal;
    //  }
    //  //else if (contype==2) {
    //  //  // Inexact path following
    //  //  vtmp = v;
    //  //  crit->update(*x,ROL::UpdateType::Trial);
    //  //  c = crit->value(*x,tol);
    //  //  v = c + pen*p;
    //  //  C = pen*pen*p/(c-v0);
    //  //  B = (C*(C+pen)*(v-v0))/pen;
    //  //  A = v0 + B/C;
    //  //  pen = std::max(pen * std::max(tau,two*p),std::pow(two*p,-(three/two)));
    //  //  t = v + p*(pen'-pen)-m(pen') = tau|vtmp-v|
    //  //}
    //  else {
    //    pen *= penscal;
    //  }
    //}
    //*outStream << "Binary Metric: " << p << std::endl;
    //*outStream << "Penalty:       " << pen << std::endl;
    //std::stringstream dname;
    //dname << ocType << "_optimal_design_ex11";
    //factory->printDesign(dname.str());

    std::stringstream dname;
    dname << ocType << "_optimal_design_ex11";
    bdalg->printDesign(dname.str());
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
