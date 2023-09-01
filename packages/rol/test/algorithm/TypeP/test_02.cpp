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

/*! \file  test_02.cpp
    \brief Validate proximal gradient algorithm.
*/

#include "ROL_TypeP_SpectralGradientAlgorithm.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_l1Objective.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

template<typename Real>
class QuadraticTypeP_Test01 : public ROL::StdObjective<Real> {
private:
  int dim_;
  std::vector<Real> a_, b_;

public:
  QuadraticTypeP_Test01(int dim) : dim_(dim) {
    a_.resize(dim);
    b_.resize(dim);
    for (int i = 0; i < dim; ++i) {
      a_[i] = static_cast<Real>(5)*static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
      b_[i] = static_cast<Real>(20)*static_cast<Real>(rand())/static_cast<Real>(RAND_MAX) - static_cast<Real>(10);
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
    list.sublist("Step").set("Type","Line Search");
    list.sublist("Status Test").set("Gradient Tolerance",1e-7);
    list.sublist("Status Test").set("Constraint Tolerance",1e-8);
    list.sublist("Status Test").set("Step Tolerance",1e-12);
    list.sublist("Status Test").set("Iteration Limit", 10000);
 
    int dim = 5;
    ROL::Ptr<ROL::StdVector<RealT>>        sol, wts, y;
    ROL::Ptr<QuadraticTypeP_Test01<RealT>> sobj;
    ROL::Ptr<ROL::l1Objective<RealT>>      nobj;
    ROL::Ptr<ROL::TypeP::SpectralGradientAlgorithm<RealT>> algo;
    std::vector<RealT> data;
    RealT err(0);

    *outStream << std::endl << "Random Diagonal LASSO Test Problem" << std::endl << std::endl;
    ROL::Ptr<std::vector<RealT>> wtsP = ROL::makePtr<std::vector<RealT>>(dim);
    ROL::Ptr<std::vector<RealT>> yP   = ROL::makePtr<std::vector<RealT>>(dim);
    wts = ROL::makePtr<ROL::StdVector<RealT>>(wtsP);
    y   = ROL::makePtr<ROL::StdVector<RealT>>(yP);
    sol = ROL::makePtr<ROL::StdVector<RealT>>(dim);
    wts->randomize(static_cast<RealT>(0),static_cast<RealT>(1));
    y->randomize(static_cast<RealT>(-5),static_cast<RealT>(5));
    sol->zero();

    nobj = ROL::makePtr<ROL::l1Objective<RealT>>(wts,y);
    sobj = ROL::makePtr<QuadraticTypeP_Test01<RealT>>(dim);

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

    algo = ROL::makePtr<ROL::TypeP::SpectralGradientAlgorithm<RealT>>(list);
    algo->run(*sol,*sobj,*nobj,*outStream);

    std::vector<RealT> xstar(dim);
    sobj->getSolution(xstar, *wtsP, *yP);
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
    *outStream << "  Max-Error = " << err << std::endl;
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
