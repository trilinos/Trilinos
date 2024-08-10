// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

// ROL sample generators
#include "ROL_SparseGridGenerator.hpp"
#include "ROL_StdTeuchosBatchManager.hpp"

typedef double RealT;

template<class Real>
class Integrand {
public:
  virtual ~Integrand() {}
  virtual Real value(const std::vector<Real> &x) = 0;
};

template<class Real>
class TestIntegrand : public Integrand<Real> {
private:
  const Real coeff1_, coeff2_;
public:
  TestIntegrand(const Real &coeff1 = 10, const Real &coeff2 = 10)
    : coeff1_(coeff1), coeff2_(coeff2) {}
  Real value(const std::vector<Real> &x) {
    Real x1 = static_cast<Real>(0.5)*(x[0] + static_cast<Real>(1));
    Real x2 = static_cast<Real>(0.5)*(x[1] + static_cast<Real>(1));
    return coeff1_*std::exp(-x1*x1) + coeff2_*std::exp(-x2*x2);
  }
};

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Teuchos::DefaultComm<int>::getComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0 && Teuchos::rank<int>(*comm)==0)
    ROL::makePtrFromRef(std::cout);
  else
    ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(1);
    ROL::Ptr<ROL::Vector<RealT> > x = ROL::makePtr<ROL::StdVector<RealT>>(x_ptr);
    ROL::QuadratureInfo info;
    info.dim        = 2;
    info.maxLevel   = 7;
    info.rule1D.clear(); info.rule1D.resize(info.dim,ROL::QUAD_CLENSHAWCURTIS);    
    info.growth1D.clear(); info.growth1D.resize(info.dim,ROL::GROWTH_DEFAULT);
    info.normalized = true;
    info.adaptive   = true;
    info.print      = !Teuchos::rank<int>(*comm);
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<ROL::StdTeuchosBatchManager<RealT,int>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::SparseGridGenerator<RealT>>(bman,info);
    /**********************************************************************************************/
    /************************* CONSTRUCT INTEGRAND FUNCTION ***************************************/
    /**********************************************************************************************/
    RealT coeff1(10), coeff2(10);
    ROL::Ptr<Integrand<RealT> > func
      = ROL::makePtr<TestIntegrand<RealT>>(coeff1,coeff2);
    /**********************************************************************************************/
    /************************* ADAPTIVE QUADRATURE ************************************************/
    /**********************************************************************************************/
    RealT tol(1.e-7);
    for (int L = 0; L < 2; ++L) {
      sampler->update(*x);
      // Initial quadrature approximation
      RealT value(0), myvalue(0), ptvalue(0);
      for (int i = sampler->start(); i < sampler->numMySamples(); ++i) {
        ptvalue = func->value(sampler->getMyPoint(i));
        myvalue += sampler->getMyWeight(i) * ptvalue;
      }
      sampler->sumAll(&myvalue,&value,1);
      *outStream << "Initial integral value: " << value << std::endl << std::endl;
      // Adaptivity
      tol *= static_cast<RealT>(0.1);
      RealT error(tol+1), incvalue(0);
      myvalue = static_cast<RealT>(0);
      std::vector<RealT> ptvalues, errors;
      while (error > tol) {
        sampler->refine();
        for (int i = sampler->start(); i < sampler->numMySamples(); ++i) {
          ptvalue = func->value(sampler->getMyPoint(i));
          myvalue += sampler->getMyWeight(i) * ptvalue;
          ptvalues.push_back(ptvalue);
        }
        error = sampler->computeError(ptvalues);
        errors.push_back(error);
        ptvalues.clear();
      }
      sampler->sumAll(&myvalue,&incvalue,1);
      value += incvalue;
      sampler->setSamples(false);
      // Print result
      *outStream << "Absolute incremental errors" << std::endl;
      for (int i = 0; i < static_cast<int>(errors.size()); ++i) {
        *outStream << "  Step " << i << ": " << errors[i] << std::endl;
      }
      *outStream << std::endl;
      *outStream << "Relative incremental errors" << std::endl;
      for (int i = 0; i < static_cast<int>(errors.size()); ++i) {
        *outStream << "  Step " << i << ": " << errors[i]/errors[0] << std::endl;
      }
      *outStream << std::endl;
      RealT one(1), half(0.5), pi(M_PI);
      RealT value_true = (coeff1+coeff2)*std::sqrt(pi)*half*std::erf(one);
      *outStream << "Integrate f(x,y) = (10*exp(-x*x) + 10*exp(-y*y)) over [0,1]x[0,1]"
                 << std::endl;
      *outStream << "  True integral value is:        "
                 << value_true << std::endl;
      *outStream << "  Approximate integral value is: "
                 << value << std::endl;
      *outStream << "  Absolute Error:                "
                 << std::abs(value-value_true) << std::endl;
      *outStream << "  Relative Error:                "
                 << std::abs(value-value_true)/std::abs(value_true) << std::endl;
    }
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
