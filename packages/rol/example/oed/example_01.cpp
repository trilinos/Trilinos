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
#include "ROL_Solver.hpp"
#include "ROL_Stream.hpp"

#include "ROL_OED_Factory.hpp"
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
    for (unsigned i = 0; i <= deg_; ++i) {
      val  += theta[i] * xpow;
      xpow *= ROL::Objective<Real>::getParameter()[0];
    }
    return val;
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &theta, Real &tol) override {
    Real xpow(1);
    for (unsigned i = 0; i <= deg_; ++i) {
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
      //ptfile << std::scientific << std::setprecision(16);
      //wtfile << std::scientific << std::setprecision(16);
      ptfile << (ub - lb) * static_cast<RealT>(i) / static_cast<RealT>(nsamp - 1) + lb << std::endl;
      wtfile << static_cast<RealT>(1) / static_cast<RealT>(nsamp)                      << std::endl;
    }
    ptfile.close();
    wtfile.close();
    auto bman = ROL::makePtr<ROL::BatchManager<RealT>>();
    auto sampler = ROL::makePtr<ROL::UserInputGenerator<RealT>>("points.txt","weights.txt",nsamp,1,bman);
    auto isampler = ROL::makePtr<ROL::UserInputGenerator<RealT>>("pointsGL.txt","weightsGL.txt",11,1,bman);

    // Setup factory
    bool homNoise = true;
    std::string regType = "Least Squares";
    std::string ocType = parlist->sublist("OED").get("Optimality Type","A");
    auto type = ROL::OED::StringToRegressionType(regType);
    auto M = ROL::makePtr<ROL::OED::StdMomentOperator<RealT>>(type,homNoise,noise);
    bool addTik = parlist->sublist("Problem").get("Use Tikhonov",false);
    if (addTik) {
      RealT beta  = parlist->sublist("Problem").get("Tikhonov Parameter",1e-4);
      auto P = ROL::makePtr<RegularizationOperator<RealT>>(beta);
      M->setPerturbation(P);
    }
    auto factory = ROL::makePtr<ROL::OED::Factory<RealT>>(model,sampler,theta,M,*parlist);
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
    
    // Generate optimization problem
    auto problem = factory->get(*parlist,sampler);
    problem->setProjectionAlgorithm(*parlist);
    problem->finalize(false,true,*outStream);
    auto test = factory->getDesign()->clone();
    test->randomize(1,2);
    problem->check(true,*outStream,test,0.1);

    // Setup ROL solver
    std::clock_t timer = std::clock();
    auto solver = ROL::makePtr<ROL::Solver<RealT>>(problem,*parlist);
    solver->solve(*outStream);
    *outStream << "  " << ocType << "-optimal design time:      "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds" << std::endl;
    factory->profile(*outStream);
    std::stringstream dname;
    dname << ocType << "_optimal_design";
    factory->printDesign(dname.str());
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
