// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test ROL::OED::GreedyObjective interface.
*/

#include "ROL_StdObjective.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_UserInputGenerator.hpp"
#include "ROL_Solver.hpp"
#include "ROL_Stream.hpp"

#include "ROL_OED_GreedyAlgorithm.hpp"
#include "ROL_OED_Factory.hpp"
#include "ROL_OED_StdMomentOperator.hpp"
#include "ROL_OED_PrintDesign.hpp"

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

  RealT errtol = static_cast<RealT>(1e2)*std::sqrt(ROL::ROL_EPSILON<RealT>());

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
    auto factory = ROL::makePtr<ROL::OED::Factory<RealT>>(model,sampler,theta,M,*parlist);

    // Build greedy algorithm
    auto cost = factory->createDesignVector();
    cost->setScalar(static_cast<RealT>(1));
    RealT budget = parlist->sublist("Problem").get("Budget",5.0);
    if (parlist->sublist("Problem").get("Use Nonuniform Cost",false)) {
      const auto& cdata = ROL::staticPtrCast<ROL::ProbabilityVector<RealT>>(cost)->getVector();
      for (int i=sampler->start(); i<sampler->numMySamples(); ++i)
        (*cdata)[i] = static_cast<RealT>(i+1)/static_cast<RealT>(nsamp);
    }

    M->generateFactors(model,theta,sampler);

    auto w0 = factory->createDesignVector();
    auto w1 = factory->createDesignVector();
    w0->setScalar(static_cast<RealT>(1));
    w1->setScalar(static_cast<RealT>(1));
    auto& w1data = *ROL::staticPtrCast<ROL::ProbabilityVector<RealT>>(w1)->getVector();
    auto ts = ROL::makePtr<ROL::OED::TraceSampler<RealT>>(theta);
    std::vector<RealT> wts(deg+1,1);
    auto c  = theta->dual().clone(); c->setScalar(static_cast<RealT>(1));

    auto objGDA = ROL::makePtr<ROL::OED::GreedyObjectiveA<RealT>>(M,ts,wts,true);
    auto objGAA = ROL::makePtr<ROL::OED::GreedyObjectiveA<RealT>>(M,ts,wts,false);
    auto objA   = ROL::makePtr<ROL::OED::Hom::ObjectiveA<RealT>>(M,theta,ts,wts);
    
    const unsigned dim(w0->dimension());
    RealT r(static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX));
    unsigned ind(std::floor(r*static_cast<RealT>(dim)));
    w1data[ind] = static_cast<RealT>(0);
    objGDA->update(*w0);
    RealT valGDA0 = objGDA->value(*w0);
    RealT incGDA0 = objGDA->computeIncrement(*w0,ind);
    objGDA->update(*w1);
    RealT valGDA1 = objGDA->value(*w1);
    objGAA->update(*w1);
    RealT valGAA1 = objGAA->value(*w1);
    RealT incGAA1 = objGAA->computeIncrement(*w1,ind);
    objGAA->update(*w0);
    RealT valGAA0 = objGAA->value(*w0);
    objA->update(*w0,ROL::UpdateType::Trial);
    RealT valA0 = objA->value(*w0,errtol);
    objA->update(*w1,ROL::UpdateType::Trial);
    RealT valA1 = objA->value(*w1,errtol);

    RealT errA = std::abs(incGDA0-incGAA1);
    errA = std::max(std::abs(valA1-valA0-incGDA0),errA);
    errA = std::max(std::abs(valA1-valA0-incGAA1),errA);
    errA = std::max(std::abs(valA0-valGDA0),errA);
    errA = std::max(std::abs(valA0-valGAA0),errA);
    errA = std::max(std::abs(valA1-valGDA1),errA);
    errA = std::max(std::abs(valA1-valGAA1),errA);
    errA = std::max(std::abs(valGAA1-valGDA1),errA);
    errA = std::max(std::abs(valGAA0-valGDA0),errA);

    *outStream << std::endl;
    *outStream << "  A Greedy Deletion: " << valGDA0 << "  " << valGDA1 << "  " << incGDA0 << "  " << valGDA1-valGDA0 << std::endl;
    *outStream << "  A Greedy Addition: " << valGAA0 << "  " << valGAA1 << "  " << incGAA1 << "  " << valGAA1-valGAA0 << std::endl;
    *outStream << "  A Optimality:      " << valA0 << "  " << valA1 << "  " << valA1-valA0 << std::endl;
    *outStream << "  A Error:           " << errA << std::endl;

    auto objGDC = ROL::makePtr<ROL::OED::GreedyObjectiveC<RealT>>(M,c,true);
    auto objGAC = ROL::makePtr<ROL::OED::GreedyObjectiveC<RealT>>(M,c,false);
    auto objC   = ROL::makePtr<ROL::OED::Hom::ObjectiveC<RealT>>(M,c);
    
    objGDC->update(*w0);
    RealT valGDC0 = objGDC->value(*w0);
    RealT incGDC0 = objGDC->computeIncrement(*w0,ind);
    objGDC->update(*w1);
    RealT valGDC1 = objGDC->value(*w1);
    objGAC->update(*w1);
    RealT valGAC1 = objGAC->value(*w1);
    RealT incGAC1 = objGAC->computeIncrement(*w1,ind);
    objGAC->update(*w0);
    RealT valGAC0 = objGAC->value(*w0);
    objC->update(*w0,ROL::UpdateType::Trial);
    RealT valC0 = objC->value(*w0,errtol);
    objC->update(*w1,ROL::UpdateType::Trial);
    RealT valC1 = objC->value(*w1,errtol);

    RealT errC = std::abs(incGDC0-incGAC1);
    errC = std::max(std::abs(valC1-valC0-incGDC0),errC);
    errC = std::max(std::abs(valC1-valC0-incGAC1),errC);
    errC = std::max(std::abs(valC0-valGDC0),errC);
    errC = std::max(std::abs(valC0-valGAC0),errC);
    errC = std::max(std::abs(valC1-valGDC1),errC);
    errC = std::max(std::abs(valC1-valGAC1),errC);
    errC = std::max(std::abs(valGAC1-valGDC1),errC);
    errC = std::max(std::abs(valGAC0-valGDC0),errC);

    *outStream << std::endl;
    *outStream << "  C Greedy Deletion: " << valGDC0 << "  " << valGDC1 << "  " << incGDC0 << "  " << valGDC1-valGDC0 << std::endl;
    *outStream << "  C Greedy Addition: " << valGAC0 << "  " << valGAC1 << "  " << incGAC1 << "  " << valGAC1-valGAC0 << std::endl;
    *outStream << "  C Optimality:      " << valC0 << "  " << valC1 << "  " << valC1-valC0 << std::endl;
    *outStream << "  C Error:           " << errC << std::endl;

    auto objGDD = ROL::makePtr<ROL::OED::GreedyObjectiveD<RealT>>(M,true);
    auto objGAD = ROL::makePtr<ROL::OED::GreedyObjectiveD<RealT>>(M,false);
    auto objD   = ROL::makePtr<ROL::OED::Hom::ObjectiveD<RealT>>(M,theta);
    
    objGDD->update(*w0);
    RealT valGDD0 = objGDD->value(*w0);
    RealT incGDD0 = -std::log(1.0-objGDD->computeIncrement(*w0,ind));
    objGDD->update(*w1);
    RealT valGDD1 = objGDD->value(*w1);
    objGAD->update(*w1);
    RealT valGAD1 = objGAD->value(*w1);
    RealT incGAD1 = std::log(1.0+objGAD->computeIncrement(*w1,ind));
    objGAD->update(*w0);
    RealT valGAD0 = objGAD->value(*w0);
    objD->update(*w0,ROL::UpdateType::Trial);
    RealT valD0 = objD->value(*w0,errtol);
    objD->update(*w1,ROL::UpdateType::Trial);
    RealT valD1 = objD->value(*w1,errtol);

    RealT errD = std::abs(incGDD0-incGAD1);
    errD = std::max(std::abs(valD1-valD0-incGDD0),errD);
    errD = std::max(std::abs(valD1-valD0-incGAD1),errD);
    errD = std::max(std::abs(valD0-valGDD0),errD);
    errD = std::max(std::abs(valD0-valGAD0),errD);
    errD = std::max(std::abs(valD1-valGDD1),errD);
    errD = std::max(std::abs(valD1-valGAD1),errD);
    errD = std::max(std::abs(valGAD1-valGDD1),errD);
    errD = std::max(std::abs(valGAD0-valGDD0),errD);

    *outStream << std::endl;
    *outStream << "  D Greedy Deletion: " << valGDD0 << "  " << valGDD1 << "  " << incGDD0 << "  " << valGDD1-valGDD0 << std::endl;
    *outStream << "  D Greedy Addition: " << valGAD0 << "  " << valGAD1 << "  " << incGAD1 << "  " << valGAD1-valGAD0 << std::endl;
    *outStream << "  D Optimality:      " << valD0 << "  " << valD1 << "  " << valD1-valD0 << std::endl;
    *outStream << "  D Error:           " << errD << std::endl;

    auto objGDI = ROL::makePtr<ROL::OED::GreedyObjectiveI<RealT>>(M,M->getFactors(),sampler,ts,wts,true);
    auto objGAI = ROL::makePtr<ROL::OED::GreedyObjectiveI<RealT>>(M,M->getFactors(),sampler,ts,wts,false);
    auto objI   = ROL::makePtr<ROL::OED::Hom::ObjectiveI<RealT>>(M,M->getFactors(),sampler);
    
    objGDI->update(*w0);
    RealT valGDI0 = objGDI->value(*w0);
    RealT incGDI0 = objGDI->computeIncrement(*w0,ind);
    objGDI->update(*w1);
    RealT valGDI1 = objGDI->value(*w1);
    objGAI->update(*w1);
    RealT valGAI1 = objGAI->value(*w1);
    RealT incGAI1 = objGAI->computeIncrement(*w1,ind);
    objGAI->update(*w0);
    RealT valGAI0 = objGAI->value(*w0);
    objI->update(*w0,ROL::UpdateType::Trial);
    RealT valI0 = objI->value(*w0,errtol);
    objI->update(*w1,ROL::UpdateType::Trial);
    RealT valI1 = objI->value(*w1,errtol);

    RealT errI = std::abs(incGDI0-incGAI1);
    errI = std::max(std::abs(valI1-valI0-incGDI0),errI);
    errI = std::max(std::abs(valI1-valI0-incGAI1),errI);
    errI = std::max(std::abs(valI0-valGDI0),errI);
    errI = std::max(std::abs(valI0-valGAI0),errI);
    errI = std::max(std::abs(valI1-valGDI1),errI);
    errI = std::max(std::abs(valI1-valGAI1),errI);
    errI = std::max(std::abs(valGAI1-valGDI1),errI);
    errI = std::max(std::abs(valGAI0-valGDI0),errI);

    *outStream << std::endl;
    *outStream << "  I Greedy Deletion: " << valGDI0 << "  " << valGDI1 << "  " << incGDI0 << "  " << valGDI1-valGDI0 << std::endl;
    *outStream << "  I Greedy Addition: " << valGAI0 << "  " << valGAI1 << "  " << incGAI1 << "  " << valGAI1-valGAI0 << std::endl;
    *outStream << "  I Optimality:      " << valI0 << "  " << valI1 << "  " << valI1-valI0 << std::endl;
    *outStream << "  I Error:           " << errI << std::endl;
    *outStream << std::endl;

    if (errA>errtol || errC>errtol || errD>errtol || errI > errtol) errorFlag++;
  }

  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0) std::cout << "End Result: TEST FAILED\n";
  else                std::cout << "End Result: TEST PASSED\n";

  return 0;

}

