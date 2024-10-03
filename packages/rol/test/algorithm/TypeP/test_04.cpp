// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_03.cpp
    \brief Validate Trust Region algorithm.
*/

#include "ROL_TypeP_TrustRegionAlgorithm.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_l1Objective.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include <random>
#include <chrono>

template<typename Real>
class QuadraticTypeP_Test01 : public ROL::StdObjective<Real> {
private:
  int dim_;
  std::vector<Real> a_, b_;

public:
  QuadraticTypeP_Test01(int dim) : dim_(dim) {
    using seed_type = std::mt19937_64::result_type;
    seed_type const seed = 123;
    std::mt19937_64 eng{seed};
    std::uniform_real_distribution<Real> distA(0.0,5.0), distB(-10.0,10.0);
    a_.resize(dim);
    b_.resize(dim);
    for (int i = 0; i < dim; ++i) {
      a_[i] = distA(eng);
      b_[i] = distB(eng);
    }
  }

  Real value(const std::vector<Real> &x, Real &tol) {
    Real val(0);
    for (int i = 0; i < dim_; ++i)
      val += static_cast<Real>(0.5)*a_[i]*x[i]*x[i] + b_[i]*x[i];
    return val;
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    for (int i = 0; i < dim_; ++i)
      g[i] = a_[i]*x[i] + b_[i];
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    for (int i = 0; i < dim_; ++i)
      hv[i] = a_[i]*v[i];
  }

  void getSolution(std::vector<Real> &x, const std::vector<Real> &wts, const std::vector<Real> &y) const {
    for (int i = 0; i < dim_; ++i)
      x[i] = (std::min(wts[i], std::max(-wts[i], a_[i]*y[i] + b_[i])) - b_[i]) / a_[i];
  }
};

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a
  // (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try {
    RealT tol = 1e2*std::sqrt(ROL::ROL_EPSILON<RealT>());

    ROL::ParameterList list;
    list.sublist("General").set("Output Level",iprint);
    list.sublist("Step").set("Type","Trust Region");
    list.sublist("Status Test").set("Gradient Tolerance",1e-1*tol);
    list.sublist("Status Test").set("Constraint Tolerance",1e-1*tol);
    list.sublist("Status Test").set("Step Tolerance",1e-3*tol);
    list.sublist("Status Test").set("Iteration Limit", 50);
    int dim = 5;
    ROL::Ptr<ROL::StdVector<RealT>>        sol, wts, y;
    ROL::Ptr<QuadraticTypeP_Test01<RealT>> sobj;
    ROL::Ptr<ROL::l1Objective<RealT>>      nobj;
    ROL::Ptr<ROL::TypeP::TrustRegionAlgorithm<RealT>> algo;
    std::vector<RealT> data;
    RealT err(0);

    *outStream << std::endl << "Random Diagonal LASSO Test Problem" << std::endl << std::endl;
    ROL::Ptr<std::vector<RealT>> wtsP = ROL::makePtr<std::vector<RealT>>(dim);
    ROL::Ptr<std::vector<RealT>> yP   = ROL::makePtr<std::vector<RealT>>(dim);
    wts = ROL::makePtr<ROL::StdVector<RealT>>(wtsP);// wts->setSeed(234);
    y   = ROL::makePtr<ROL::StdVector<RealT>>(yP);  // y->setSeed(345);
    sol = ROL::makePtr<ROL::StdVector<RealT>>(dim);
    wts->randomize(static_cast<RealT>(0),static_cast<RealT>(1));
    y->randomize(static_cast<RealT>(-5),static_cast<RealT>(5));

    nobj = ROL::makePtr<ROL::l1Objective<RealT>>(wts,y);
    sobj = ROL::makePtr<QuadraticTypeP_Test01<RealT>>(dim);

    std::vector<RealT> xstar(dim);
    sobj->getSolution(xstar, *wtsP, *yP);
    RealT xmax(0);
    for (int i = 0; i < dim; ++i)
      xmax = std::max(xmax,std::abs(xstar[i]));

    // Check derivatives of smooth function
    ROL::Ptr<ROL::Vector<RealT>> xd = sol->clone();
    xd->randomize(-1.0,1.0);
    ROL::Ptr<ROL::Vector<RealT>> yd = sol->clone();
    yd->randomize(-1.0,1.0);
    ROL::Ptr<ROL::Vector<RealT>> zd = sol->clone();
    zd->randomize(-1.0,1.0);
    sobj->checkGradient(*xd,*yd,true,*outStream);
    sobj->checkHessVec(*xd,*yd,true,*outStream);
    sobj->checkHessSym(*xd,*yd,*zd,true,*outStream);

    list.sublist("Step").sublist("Trust Region").sublist("TRN").sublist("Solver").set("Subproblem Solver", "SPG");  
    sol->zero();
    algo = ROL::makePtr<ROL::TypeP::TrustRegionAlgorithm<RealT>>(list);
    auto begin = std::chrono::high_resolution_clock::now();
    algo->run(*sol,*sobj,*nobj,*outStream);
    auto end   = std::chrono::high_resolution_clock::now();
    *outStream << "  Optimization Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << " microseconds" << std::endl;

    err = static_cast<RealT>(0);
    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(sol)->getVector();
    *outStream << "  Result:   ";
    for (int i = 0; i < dim; ++i) {
      *outStream << "  x" << i+1 << " = " << data[i];
      err = std::max(err,std::abs(data[i]-xstar[i]));
    }
    *outStream << std::endl;
    *outStream << "  Truth:    ";
    for (int i = 0; i < dim; ++i) {
      *outStream << "  x" << i+1 << " = " << xstar[i];
    }
    *outStream << std::endl;
    *outStream << "  Max Relative Error = " << err/xmax << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    list.sublist("Step").sublist("Trust Region").sublist("TRN").sublist("Solver").set("Subproblem Solver", "Simplified SPG");  
    sol->zero();
    algo = ROL::makePtr<ROL::TypeP::TrustRegionAlgorithm<RealT>>(list);
    begin = std::chrono::high_resolution_clock::now();
    algo->run(*sol,*sobj,*nobj,*outStream);
    end   = std::chrono::high_resolution_clock::now();
    *outStream << "  Optimization Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << " microseconds" << std::endl;

    err = static_cast<RealT>(0);
    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(sol)->getVector();
    *outStream << "  Result:   ";
    for (int i = 0; i < dim; ++i) {
      *outStream << "  x" << i+1 << " = " << data[i];
      err = std::max(err,std::abs(data[i]-xstar[i]));
    }
    *outStream << std::endl;
    *outStream << "  Truth:    ";
    for (int i = 0; i < dim; ++i) {
      *outStream << "  x" << i+1 << " = " << xstar[i];
    }
    *outStream << std::endl;
    *outStream << "  Max Relative Error = " << err/xmax << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    list.sublist("Step").sublist("Trust Region").sublist("TRN").sublist("Solver").set("Subproblem Solver", "NCG");  
    sol->zero();
    algo = ROL::makePtr<ROL::TypeP::TrustRegionAlgorithm<RealT>>(list);
    begin = std::chrono::high_resolution_clock::now();
    algo->run(*sol,*sobj,*nobj,*outStream);
    end   = std::chrono::high_resolution_clock::now();
    *outStream << "  Optimization Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << " microseconds" << std::endl;

    err = static_cast<RealT>(0);
    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(sol)->getVector();
    *outStream << "  Result:   ";
    for (int i = 0; i < dim; ++i) {
      *outStream << "  x" << i+1 << " = " << data[i];
      err = std::max(err,std::abs(data[i]-xstar[i]));
    }
    *outStream << std::endl;
    *outStream << "  Truth:    ";
    for (int i = 0; i < dim; ++i) {
      *outStream << "  x" << i+1 << " = " << xstar[i];
    }
    *outStream << std::endl;
    *outStream << "  Max Relative Error = " << err/xmax << std::endl;
    errorFlag += (err > tol ? 1 : 0);
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
