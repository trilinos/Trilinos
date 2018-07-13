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

#ifndef ROL_DISTRIBUTION_HPP
#define ROL_DISTRIBUTION_HPP

#include "ROL_Types.hpp"

#include <cmath>
#include <iostream>

namespace ROL {

template<class Real>
class Distribution {
public: 
  virtual ~Distribution(void) {}

  virtual Real evaluatePDF(const Real input) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Distribution): evaluatePDF not implemented!");
    return 0.;
  }

  virtual Real evaluateCDF(const Real input) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Distribution): evaluateCDF not implemented!");
    return 0.;
  }

  virtual Real integrateCDF(const Real input) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Distribution): integrateCDF not implemented!");
    return 0.;
  }

  virtual Real invertCDF(const Real input) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Distribution): invertCDF not implemented!");
    return 0.;
  }

  virtual Real moment(const size_t m) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Distribution): moment not implemented!");
    return 0.;
  }

  virtual Real lowerBound(void) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Distribution): lowerBound not implemented!");
    return 0.;
  }

  virtual Real upperBound(void) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Distribution): upperBound not implemented!");
    return 0.;
  }

  virtual void test(std::ostream &outStream = std::cout) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Distribution): test not implemented!");
  }
 
protected:
  void test(const std::vector<Real> &X, const std::vector<int> &T,
            std::ostream &outStream = std::cout ) const {
    size_t size = X.size();
    for ( size_t k = 0; k < size; k++ ) {
      if ( T[k] == 0 ) {
        test_onesided(X[k],outStream);
      }
      else {
        test_centered(X[k],outStream);
      }
    }
    size_t order = 2;
    test_moment(order,outStream);
  }

private:
  void test_onesided(const Real x, std::ostream &outStream = std::cout) const {
    Real X = x, vx = 0., vy = 0., dv = 0., t = 1., diff = 0., err = 0.;
    try {
      vx   = evaluateCDF(X);
      vy   = 0.0;
      dv   = evaluatePDF(X);
      t    = 1.0;
      diff = 0.0;
      err  = 0.0;
      outStream << std::scientific << std::setprecision(11);
      outStream << std::right << std::setw(20) << "CHECK DENSITY: f(x) = cdf(x) with x = "
                                               << X << " is correct?" << std::endl;
      outStream << std::right << std::setw(20) << "t"
                              << std::setw(20) << "f'(x)"
                              << std::setw(20) << "(f(x+t)-f(x))/t"
                              << std::setw(20) << "Error"
                              << std::endl;
      for (int i = 0; i < 13; i++) {
        vy = evaluateCDF(X+t);
        diff = (vy-vx)/t;
        err = std::abs(diff-dv);
        outStream << std::scientific << std::setprecision(11) << std::right
                  << std::setw(20) << t
                  << std::setw(20) << dv
                  << std::setw(20) << diff
                  << std::setw(20) << err
                  << std::endl;
        t *= 0.1;
      }
      outStream << std::endl;
    }
    catch(std::exception &e) {
      outStream << "Either evaluateCDF or evaluatePDF is not implemented!"
                << std::endl << std::endl;
    }
    // CHECK INTCDF
    try {
      vx   = integrateCDF(X);
      vy   = 0.0;
      dv   = evaluateCDF(X);
      t    = 1.0;
      diff = 0.0;
      err  = 0.0;
      outStream << std::scientific << std::setprecision(11);
      outStream << std::right << std::setw(20) << "CHECK DENSITY: f(x) = intcdf(x) with x = "
                                               << X << " is correct?" << std::endl;
      outStream << std::right << std::setw(20) << "t"
                              << std::setw(20) << "f'(x)"
                              << std::setw(20) << "(f(x+t)-f(x))/t"
                              << std::setw(20) << "Error"
                              << std::endl;
      for (int i = 0; i < 13; i++) {
        vy = integrateCDF(X+t);
        diff = (vy-vx)/t;
        err = std::abs(diff-dv);
        outStream << std::scientific << std::setprecision(11) << std::right
                  << std::setw(20) << t
                  << std::setw(20) << dv
                  << std::setw(20) << diff
                  << std::setw(20) << err
                  << std::endl;
        t *= 0.1;
      }
      outStream << std::endl;
    }
    catch(std::exception &e) {
      outStream << "Either evaluateCDF or integrateCDF is not implemented!"
                << std::endl << std::endl;
    }
    // CHECK INVCDF
    try {
      vx = evaluateCDF(X);
      vy = invertCDF(vx);
      err = std::abs(x-vy);
      outStream << std::scientific << std::setprecision(11);
      outStream << std::right << std::setw(20) << "CHECK DENSITY: f(x) = invcdf(x) with x = "
                                               << X << " is correct?" << std::endl;
      outStream << std::right << std::setw(20) << "cdf(x)"
                              << std::setw(20) << "invcdf(cdf(x))"
                              << std::setw(20) << "Error"
                              << std::endl;
      outStream << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << vx
                << std::setw(20) << vy
                << std::setw(20) << err
                << std::endl << std::endl;
    }
    catch(std::exception &e) {
      outStream << "Either evaluateCDF or invertCDF is not implemented!"
                << std::endl << std::endl;
    }
  }

  void test_centered(const Real x, std::ostream &outStream = std::cout) const {
    Real X  = x, vx = 0., vy = 0., dv = 0., t = 1., diff = 0., err = 0.;
    try {
      vx   = 0.0;
      vy   = 0.0;
      dv   = evaluatePDF(X);
      t    = 1.0;
      diff = 0.0;
      err  = 0.0;
      outStream << std::scientific << std::setprecision(11);
      outStream << std::right << std::setw(20) << "CHECK DENSITY: f(x) = cdf(x) with x = "
                                               << X << " is correct?" << std::endl;
      outStream << std::right << std::setw(20) << "t"
                              << std::setw(20) << "f'(x)"
                              << std::setw(20) << "(f(x+t)-f(x-t))/2t"
                              << std::setw(20) << "Error"
                              << std::endl;
      for (int i = 0; i < 13; i++) {
        vx = evaluateCDF(X+t);
        vy = evaluateCDF(X-t);
        diff = 0.5*(vx-vy)/t;
        err = std::abs(diff-dv);
        outStream << std::scientific << std::setprecision(11) << std::right
                  << std::setw(20) << t
                  << std::setw(20) << dv
                  << std::setw(20) << diff
                  << std::setw(20) << err
                  << std::endl;
        t *= 0.1;
      }
      outStream << "\n";
    }
    catch(std::exception &e) {
      outStream << "Either evaluateCDF or evaluatePDF is not implemented!"
                << std::endl << std::endl;
    }
    // CHECK INTCDF
    try {
      vx   = 0.0;
      vy   = 0.0;
      dv   = evaluateCDF(X);
      t    = 1.0;
      diff = 0.0;
      err  = 0.0;
      outStream << std::scientific << std::setprecision(11);
      outStream << std::right << std::setw(20) << "CHECK DENSITY: f(x) = intcdf(x) with x = "
                                               << X << " is correct?" << std::endl;
      outStream << std::right << std::setw(20) << "t"
                              << std::setw(20) << "f'(x)"
                              << std::setw(20) << "(f(x+t)-f(x-t))/2t"
                              << std::setw(20) << "Error"
                              << std::endl;
      for (int i = 0; i < 13; i++) {
        vx = integrateCDF(X+t);
        vy = integrateCDF(X-t);
        diff = 0.5*(vx-vy)/t;
        err = std::abs(diff-dv);
        outStream << std::scientific << std::setprecision(11) << std::right
                  << std::setw(20) << t
                  << std::setw(20) << dv
                  << std::setw(20) << diff
                  << std::setw(20) << err
                  << std::endl;
        t *= 0.1;
      }
      outStream << std::endl;
    }
    catch(std::exception &e) {
      outStream << "Either evaluateCDF or integrateCDF is not implemented!"
                << std::endl << std::endl;
    }
    // CHECK INVCDF
    try {
      vx = evaluateCDF(X);
      vy = invertCDF(vx);
      err = std::abs(X-vy);
      outStream << std::scientific << std::setprecision(11);
      outStream << std::right << std::setw(20) << "CHECK DENSITY: f(x) = invcdf(x) with x = "
                                               << X << " is correct?" << std::endl;
      outStream << std::right << std::setw(20) << "cdf(x)"
                              << std::setw(20) << "invcdf(cdf(x))"
                              << std::setw(20) << "Error"
                              << std::endl;
      outStream << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << vx
                << std::setw(20) << vy
                << std::setw(20) << err
                << std::endl << std::endl;
    }
    catch(std::exception &e) {
      outStream << "Either evaluateCDF or invertCDF is not implemented!"
                << std::endl << std::endl;
    }
  }

  void test_moment(const size_t order, std::ostream &outStream = std::cout) const {
    try {
      const size_t numPts = 10000;
      Real pt = 0., wt = 1./(Real)numPts;
      std::vector<Real> mVec(order,0.);
      for (size_t i = 0; i < numPts; i++) {
        pt = invertCDF((Real)rand()/(Real)RAND_MAX);
        mVec[0] += wt*pt;
        for (size_t q = 1; q < order; q++) {
          mVec[q] += wt*std::pow(pt,q+1);
        }
      }
      outStream << std::scientific << std::setprecision(0);
      outStream << std::right << std::setw(20) << "CHECK DENSITY: Check first " << order
                << " moments against Monte Carlo using " << numPts << " samples"
                << std::endl;
      outStream << std::setw(20) << "Error should be O(" << 1./std::sqrt(numPts) << ")" << std::endl;
      outStream << std::scientific << std::setprecision(11);
      for (size_t q = 0; q < order; q++) {
        outStream << std::setw(20) << "Error in " << q+1 << " moment: "
                  << std::abs(mVec[q]-moment(q+1)) << std::endl;
      }
      outStream << std::endl;
    }
    catch(std::exception &e) {
      outStream << "moment is not implemented!"
                << std::endl << std::endl;
    }
  }
};

}

#endif
