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

namespace ROL {

enum EDistribution {
  DISTRIBUTION_DIRAC = 0, 
  DISTRIBUTION_GAUSSIAN,
  DISTRIBUTION_UNIFORM, 
  DISTRIBUTION_LOGISTIC,
  DISTRIBUTION_TRIANGLE,
  DISTRIBUTION_PARABOLIC,
  DISTRIBUTION_RAISEDCOSINE,
  DISTRIBUTION_LAPLACE,
  DISTRIBUTION_CAUCHY,
  DISTRIBUTION_SMALE,
  DISTRIBUTION_LAST
};

template<class Real>
class Distribution {
private:
  EDistribution ed_;
  std::vector<Real> data_;

public: 

  Distribution(EDistribution ed) : ed_(ed) {
    switch(ed) {
      case DISTRIBUTION_DIRAC:
        data_.resize(1,0.0);
        data_[0] = 0.0;
        break;
      case DISTRIBUTION_GAUSSIAN:
        data_.resize(2,0.0);
        data_[0] = 0.0; // Mean 
        data_[1] = 1.0; // Standard Deviation
        break;
      case DISTRIBUTION_UNIFORM:
        data_.resize(2,0.0); 
        data_[0] = 0.0; // Lower Bound
        data_[1] = 1.0; // Upper Bound
        break;
      case DISTRIBUTION_LOGISTIC:
        data_[0] = 0.0; // Mean
        data_[1] = 1.0; // Variance Parameter
        break;
      case DISTRIBUTION_TRIANGLE:
        data_.resize(3,0.0); 
        data_[0] = 0.0; // Lower Bound
        data_[1] = 0.5; // Mid Point
        data_[2] = 1.0; // Upper Bound
        break;
      case DISTRIBUTION_PARABOLIC:
        data_.resize(2,0.0);
        data_[0] = 0.0; // Upper Bound
        data_[1] = 1.0; // Lower Bound
        break;
      case DISTRIBUTION_RAISEDCOSINE:
        data_.resize(2,0.0);
        data_[0] = 0.5; // Mean
        data_[1] = 0.5; // Variance Parameter
        break;
      case DISTRIBUTION_LAPLACE:
        data_.resize(2,0.0);
        data_[0] = 0.0;
        data_[1] = 1.0;
        break;
      case DISTRIBUTION_CAUCHY:
        data_.resize(2,0.0);
        data_[0] = 0.0;
        data_[1] = 1.0;
        break;
      case DISTRIBUTION_SMALE:
        data_.resize(2,0.0);
        data_[0] = 0.0;
        data_[1] = 1.0;
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                          ">>> ERROR (ROL::Distribution): Distribution not defined!");
    }
  }

  Distribution(EDistribution ed, std::vector<Real> &data) : ed_(ed), data_(data) {}

  Real pdf(Real input) {
    Real val = 0.0;
    switch(ed_) {
      case DISTRIBUTION_DIRAC:        val = dirac_pdf(input);        break;
      case DISTRIBUTION_GAUSSIAN:     val = gaussian_pdf(input);     break;
      case DISTRIBUTION_UNIFORM:      val = uniform_pdf(input);      break;
      case DISTRIBUTION_LOGISTIC:     val = logistic_pdf(input);     break;
      case DISTRIBUTION_TRIANGLE:     val = triangle_pdf(input);     break;
      case DISTRIBUTION_PARABOLIC:    val = parabolic_pdf(input);    break;
      case DISTRIBUTION_RAISEDCOSINE: val = raisedcosine_pdf(input); break;
      case DISTRIBUTION_LAPLACE:      val = laplace_pdf(input);      break;
      case DISTRIBUTION_CAUCHY:       val = cauchy_pdf(input);       break;
      case DISTRIBUTION_SMALE:        val = smale_pdf(input);        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                          ">>> ERROR (ROL::Distribution): Distribution not defined!");
    }
    return val;
  }

  Real cdf(Real input) {
    Real val = 0.0;
    switch(ed_) {
      case DISTRIBUTION_DIRAC:        val = dirac_cdf(input);        break;
      case DISTRIBUTION_GAUSSIAN:     val = gaussian_cdf(input);     break;
      case DISTRIBUTION_UNIFORM:      val = uniform_cdf(input);      break;
      case DISTRIBUTION_LOGISTIC:     val = logistic_cdf(input);     break;
      case DISTRIBUTION_TRIANGLE:     val = triangle_cdf(input);     break;
      case DISTRIBUTION_PARABOLIC:    val = parabolic_cdf(input);    break;
      case DISTRIBUTION_RAISEDCOSINE: val = raisedcosine_cdf(input); break;
      case DISTRIBUTION_LAPLACE:      val = laplace_cdf(input);      break;
      case DISTRIBUTION_CAUCHY:       val = cauchy_cdf(input);       break;
      case DISTRIBUTION_SMALE:        val = smale_cdf(input);        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                          ">>> ERROR (ROL::Distribution): Distribution not defined!");
    }
    return val;
  }

  Real intcdf(Real input) {
    Real val = 0.0;
    switch(ed_) {
      case DISTRIBUTION_DIRAC:        val = dirac_intcdf(input);        break;
      case DISTRIBUTION_GAUSSIAN:     val = gaussian_intcdf(input);     break;
      case DISTRIBUTION_UNIFORM:      val = uniform_intcdf(input);      break;
      case DISTRIBUTION_LOGISTIC:     val = logistic_intcdf(input);     break;
      case DISTRIBUTION_TRIANGLE:     val = triangle_intcdf(input);     break;
      case DISTRIBUTION_PARABOLIC:    val = parabolic_intcdf(input);    break;
      case DISTRIBUTION_RAISEDCOSINE: val = raisedcosine_intcdf(input); break;
      case DISTRIBUTION_LAPLACE:      val = laplace_intcdf(input);      break;
      case DISTRIBUTION_CAUCHY:       val = cauchy_intcdf(input);       break;
      case DISTRIBUTION_SMALE:        val = smale_intcdf(input);        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                          ">>> ERROR (ROL::Distribution): Distribution not defined!");
    }
    return val;
  }

  Real invcdf(Real input) {
    Real val = 0.0;
    switch(ed_) {
      case DISTRIBUTION_DIRAC:        val = dirac_invcdf(input);        break;
      case DISTRIBUTION_GAUSSIAN:     val = gaussian_invcdf(input);     break;
      case DISTRIBUTION_UNIFORM:      val = uniform_invcdf(input);      break;
      case DISTRIBUTION_LOGISTIC:     val = logistic_invcdf(input);     break;
      case DISTRIBUTION_TRIANGLE:     val = triangle_invcdf(input);     break;
      case DISTRIBUTION_PARABOLIC:    val = parabolic_invcdf(input);    break;
      case DISTRIBUTION_RAISEDCOSINE: val = raisedcosine_invcdf(input); break;
      case DISTRIBUTION_LAPLACE:      val = laplace_invcdf(input);      break;
      case DISTRIBUTION_CAUCHY:       val = cauchy_invcdf(input);       break;
      case DISTRIBUTION_SMALE:        val = smale_invcdf(input);        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                          ">>> ERROR (ROL::Distribution): Distribution not defined!");
    }
    return val;
  }
 
  void test(Real x) {
    int size;
    std::vector<Real> X;
    std::vector<int> T;
    switch (ed_) {
      case DISTRIBUTION_DIRAC:
      case DISTRIBUTION_GAUSSIAN:
      case DISTRIBUTION_LOGISTIC:
      case DISTRIBUTION_LAPLACE:
      case DISTRIBUTION_CAUCHY:
      case DISTRIBUTION_SMALE: 
        size = 1;
        X.resize(size,0.0);
        T.resize(size,0);
        X[0] = 4.0*(Real)rand()/(Real)RAND_MAX - 2.0;
        break;
      case DISTRIBUTION_UNIFORM:
      case DISTRIBUTION_PARABOLIC:
        size = 5;
        X.resize(size,0.0);
        T.resize(size,0);
        X[0] = data_[0]-4.0*(Real)rand()/(Real)RAND_MAX; 
        T[0] = 0;
        X[1] = data_[0]; 
        T[1] = 1;
        X[2] = (data_[1]-data_[0])*(Real)rand()/(Real)RAND_MAX + data_[0]; 
        T[2] = 0;
        X[3] = data_[1]; 
        T[3] = 1;
        X[4] = data_[1]+4.0*(Real)rand()/(Real)RAND_MAX; 
        T[4] = 0;
        break;
      case DISTRIBUTION_RAISEDCOSINE:
        size = 5;
        X.resize(size,0.0);
        T.resize(size,0);
        X[0] = data_[0]-data_[1]-4.0*(Real)rand()/(Real)RAND_MAX; 
        T[0] = 0;
        X[1] = data_[0]-data_[1]; 
        T[1] = 1;
        X[2] = (2.0*data_[1])*(Real)rand()/(Real)RAND_MAX + (data_[0]-data_[1]); 
        T[2] = 0;
        X[3] = data_[0]+data_[1]; 
        T[3] = 1;
        X[4] = data_[0]+data_[1]+4.0*(Real)rand()/(Real)RAND_MAX; 
        T[4] = 0;
        break;
      case DISTRIBUTION_TRIANGLE:
        size = 7;
        X.resize(size,0.0);
        T.resize(size,0);
        X[0] = data_[0]-4.0*(Real)rand()/(Real)RAND_MAX; 
        T[0] = 0;
        X[1] = data_[0]; 
        T[1] = 1;
        X[2] = (data_[1]-data_[0])*(Real)rand()/(Real)RAND_MAX + data_[0]; 
        T[2] = 0;
        X[3] = data_[1]; 
        T[3] = 1;
        X[4] = (data_[2]-data_[1])*(Real)rand()/(Real)RAND_MAX + data_[1]; 
        T[4] = 0;
        X[5] = data_[2]; 
        T[5] = 1;
        X[6] = data_[2]+4.0*(Real)rand()/(Real)RAND_MAX; 
        T[6] = 0;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                          ">>> ERROR (ROL::Distribution): Distribution not defined!");
    }
    for ( int k = 0; k < size; k++ ) {
      if ( T[k] == 0 ) {
        test_onesided(X[k]);
      }
      else {
        test_centered(X[k]);
      }
    }
  }

private:

  // Dirac Distribution 
  Real dirac_pdf(Real input) {
    return ((input==data_[0]) ? 1.0 : 0.0);
  }
  Real dirac_cdf(Real input) {
    return ((input >= data_[0]) ? 1.0 : 0.0);
  }
  Real dirac_intcdf(Real input) {
    return ((input < data_[0]) ? 0.0 : input);
  }
  Real dirac_invcdf(Real input) {
    return data_[0];
  }
  // Gaussian Distribution
  Real gaussian_pdf(Real input) {
    return std::exp(-std::pow(input-data_[0],2.0)/(2.0*data_[1]))/(std::sqrt(2.0*M_PI*data_[1]));
  }
  Real gaussian_cdf(Real input) {
    return 0.5*(1.0+erf((input-data_[0])/std::sqrt(2.0*data_[1])));
  }
  Real gaussian_intcdf(Real input) {
    // NOT IMPLEMENTED
    return ((input < data_[0]) ? 0.0 : input);
  }
  Real gaussian_invcdf(Real input) {
    std::vector<Real> coeff;
    Real x   = 2.0*input - 1.0;
    Real c   = 1.0;
    Real tmp = c * (std::sqrt(M_PI)/2.0 * x);
    Real val = tmp;
    coeff.push_back(c);
    int  cnt = 1;
    while (std::abs(tmp) > std::sqrt(ROL_EPSILON)*std::abs(val)) {
      c = 0.0;
      for ( unsigned i = 0; i < coeff.size(); i++ ) {
        c += coeff[i]*coeff[coeff.size()-1-i]/((i+1)*(2*i+1));
      }
      tmp  = c/(2.0*(Real)cnt+1.0) * std::pow(std::sqrt(M_PI)/2.0 * x,2.0*(Real)cnt+1.0);
      val += tmp;
      coeff.push_back(c);
      cnt++;
    }
    return std::sqrt(2*data_[1])*val + data_[0];
  }
  // Uniform Distribution
  Real uniform_pdf(Real input) {
    return ((input >= data_[0] && input <= data_[1]) ? 1.0/(data_[1]-data_[0]) : 0.0);
  }
  Real uniform_cdf(Real input) {
    return ((input < data_[0]) ? 0.0 : ((input > data_[1]) ? 1.0 :
            (input-data_[0])/(data_[1]-data_[0])));
  }
  Real uniform_intcdf(Real input) {
    return ((input < data_[0]) ? 0.0 : ((input > data_[1]) ? input - 0.5*(data_[0]+data_[1]) : 
              0.5*std::pow(input-data_[0],2.0)/(data_[1]-data_[0])));
  }
  Real uniform_invcdf(Real input) {
    return data_[0] + input*(data_[1]-data_[0]);
  }
  // Logistic Distribution
  Real logistic_pdf(Real input) {
    Real val = std::exp(-(input-data_[0])/data_[1]);
    return val/(data_[1]*std::pow(1.0+val,2.0));
  }
  Real logistic_cdf(Real input) {
    Real val = std::exp(-(input-data_[0])/data_[1]);
    return 1.0/(1.0+val);
  }
  Real logistic_intcdf(Real input) {
    Real val = std::exp(-(input-data_[0])/data_[1]);
    return (input-data_[0]) + data_[1]*std::log(1.0+val);
  }
  Real logistic_invcdf(Real input) {
    return data_[0] + data_[1]*std::log(input/(1.0-input));
  }
  // Triangle Distribution
  Real triangle_pdf(Real input) {
    Real a = data_[0], b = data_[1], c = data_[2];
    Real d1 = b-a, d2 = c-b, d = c-a;
    return ((input >= a && input < b) ? 2.0*(input-a)/(d*d1) :
           ((input >= b && input < c) ? 2.0*(c-input)/(d*d2) : 
             0.0));
  }
  Real triangle_cdf(Real input) {
    Real a = data_[0], b = data_[1], c = data_[2];
    Real d1 = b-a, d2 = c-b, d = c-a;
    return ((input < a) ? 0.0 : 
           ((input >= a && input < b) ? 
             std::pow(input-a,2.0)/(d*d1) :
           ((input >= b && input < c) ? 
             1.0-std::pow(c-input,2.0)/(d*d2) : 
             1.0)));
  }
  Real triangle_intcdf(Real input) {
    Real a = data_[0], b = data_[1], c = data_[2];
    Real d1 = b-a, d2 = c-b, d = c-a;
    return ((input < a) ? 0.0 : 
           ((input >= a && input < b) ? 
             std::pow(input-a,3.0)/(3.0*d*d1) : 
           ((input >= b && input < c) ?
             d1*d1/(3.0*d)+(input-b)+(std::pow(c-input,3.0)-d2*d2*d2)/(3.0*d*d2) :
             d1*d1/(3.0*d)+(input-b)-d2*d2/(3.0*d))));
  }
  Real triangle_invcdf(Real input) {
    Real a = data_[0], b = data_[1], c = data_[2];
    Real d1 = b-a, d2 = c-b, d = c-a;
    return ((input <= d1/d) ? a + std::sqrt(input*d1*d) : 
             c - std::sqrt((1.0-input)*d2*d));
  }
  // Parabolic Distribution
  Real parabolic_pdf(Real input) {
    Real scale = 6.0/std::pow(data_[1]-data_[0],3.0);
    return ((input >= data_[0] && input <= data_[1]) ? scale*(input-data_[0])*(data_[1]-input) : 0.0);
  }
  Real parabolic_cdf(Real input) {
    Real d1 = data_[1]-data_[0];
    Real d2 = d1*d1;
    Real d3 = d2*d1;
    Real v1 = input-data_[0];
    Real v2 = v1*v1;
    Real v3 = v1*v2;
    return ((input < data_[0]) ? 0.0 : ((input > data_[1]) ? 1.0 : 
            3.0*v2/d2 - 2.0*v3/d3));
  }
  Real parabolic_intcdf(Real input) {
    Real d1 = data_[1]-data_[0];
    Real d2 = d1*d1;
    Real d3 = d2*d1;
    Real v1 = input-data_[0];
    Real v2 = v1*v1;
    Real v3 = v1*v2;
    Real v4 = v1*v3;
    return ((input < data_[0]) ? 0.0 : 
           ((input > data_[1]) ? input - 0.5*d1 : 
             v3/d2 - 0.5*v4/d3));
  }
  Real parabolic_invcdf(Real input) {
    Real a  = data_[0] - data_[1];
    Real b  = data_[0] + data_[1];
    Real c  = 0.0;
    Real fa = parabolic_cdf(a) - input;
    Real fb = parabolic_cdf(b) - input;
    Real fc = 0.0;
    Real sa = ((fa < 0.0) ? -1.0 : ((fa > 0.0) ? 1.0 : 0.0));
    Real sb = ((fb < 0.0) ? -1.0 : ((fb > 0.0) ? 1.0 : 0.0));
    Real sc = 0.0;
    for (int i = 0; i < 100; i++) {
      c  = (a+b)*0.5;
      fc = parabolic_cdf(c) - input;
      sc = ((fc < 0.0) ? -1.0 : ((fc > 0.0) ? 1.0 : 0.0));
      if ( fc == 0.0 || (b-a)*0.5 < ROL_EPSILON ) {
        break;
      }
      if ( sc == sa ) { a = c; fa = fc; sa = sc; }
      else            { b = c; fb = fc; sb = sc; }
    } 
    return c;
  }
  // Raised Cosine Distribution
  Real raisedcosine_pdf(Real input) {
    Real a = data_[0]-data_[1];
    Real b = data_[0]+data_[1];
    return ((input >= a && input <= b) ? (1.0+std::cos(M_PI*(input-data_[0])/data_[1]))/(2.0*data_[1]) : 0.0);
  }
  Real raisedcosine_cdf(Real input) {
    Real a = data_[0]-data_[1];
    Real b = data_[0]+data_[1];
    return ((input < a) ? 0.0 : ((input > b) ? 1.0 : 
            0.5*(1.0+(input-data_[0])/data_[1]+std::sin(M_PI*(input-data_[0])/data_[1])/M_PI)));
  }
  Real raisedcosine_intcdf(Real input) {
    Real a = data_[0]-data_[1];
    Real b = data_[0]+data_[1];
    Real v = input-data_[0];
    return ((input < a) ? 0.0 : ((input > b) ? input - data_[1] : 
            0.5*(v+0.5*v*v/data_[1]-data_[1]*((std::cos(M_PI*v/data_[1])+1.0)/(M_PI*M_PI)-0.5))));
  }
  Real raisedcosine_invcdf(Real input) {
    Real a  = data_[0] - data_[1];
    Real b  = data_[0] + data_[1];
    Real c  = 0.0;
    Real fa = raisedcosine_cdf(a) - input;
    Real fb = raisedcosine_cdf(b) - input;
    Real fc = 0.0;
    Real sa = ((fa < 0.0) ? -1.0 : ((fa > 0.0) ? 1.0 : 0.0));
    Real sb = ((fb < 0.0) ? -1.0 : ((fb > 0.0) ? 1.0 : 0.0));
    Real sc = 0.0;
    for (int i = 0; i < 100; i++) {
      c  = (a+b)*0.5;
      fc = raisedcosine_cdf(c) - input;
      sc = ((fc < 0.0) ? -1.0 : ((fc > 0.0) ? 1.0 : 0.0));
      if ( fc == 0.0 || (b-a)*0.5 < ROL_EPSILON ) {
        break;
      }
      if ( sc == sa ) { a = c; fa = fc; sa = sc; }
      else            { b = c; fb = fc; sb = sc; }
    } 
    return c;
  }
  // Laplace Distribution
  Real laplace_pdf(Real input) {
    return 0.5*std::exp(-std::abs(input-data_[0])/data_[1])/data_[1];
  }
  Real laplace_cdf(Real input) {
    return ((input < data_[0]) ? 0.5*std::exp((input-data_[0])/data_[1]) : 
            1.0-0.5*std::exp(-(input-data_[0])/data_[1]));
  }
  Real laplace_intcdf(Real input) {
    return ((input < data_[0]) ? 0.5*data_[1]*std::exp((input-data_[0])/data_[1]) : 
            (input-data_[0])+0.5*data_[1]*std::exp(-(input-data_[0])/data_[1]));
  }
  Real laplace_invcdf(Real input) {
    Real sgn = ((input < 0.5) ? -1.0 : ((input > 0.5) ? 1.0 : 0.0));
    return data_[0] - data_[1]*sgn*std::log(1.0-2.0*std::abs(input-0.5));
  }
  // Cauchy Distribution
  Real cauchy_pdf(Real input) {
    return 1.0/(M_PI*data_[1]*(1.0+std::pow((input-data_[0])/data_[1],2.0)));
  }
  Real cauchy_cdf(Real input) {
    return 0.5+atan((input-data_[0])/data_[1])/M_PI;
  }
  Real cauchy_intcdf(Real input) {
    Real v = input-data_[0];
    return 0.5*input + (v*atan(v/data_[1]) - 0.5*data_[1]*std::log(v*v+data_[1]*data_[1]))/M_PI;
  }
  Real cauchy_invcdf(Real input) {
    return data_[0]+data_[1]*tan(M_PI*(input-0.5));
  }
  // Smale Distribution
  Real smale_pdf(Real input) {
    Real val  = std::pow(input-data_[0],2.0)+4.0*data_[1]*data_[1];
    Real root = std::sqrt(val);
    return 2.0*data_[1]*data_[1]/(val*root);
  }
  Real smale_cdf(Real input) {
    Real val  = std::pow(input-data_[0],2.0)+4.0*data_[1]*data_[1];
    Real root = std::sqrt(val);
    return 0.5*(1.0+input/root);
  }
  Real smale_intcdf(Real input) {
    Real val  = std::pow(input-data_[0],2.0)+4.0*data_[1]*data_[1];
    Real root = std::sqrt(val);
    return 0.5*(input+root);
  }
  Real smale_invcdf(Real input) {
    Real x   = data_[0];
    Real fx  = smale_cdf(x)-input;
    Real s   = 0.0;
    Real xs  = 0.0;
    Real a   = 1.0;
    Real tmp = 0.0;
    for (int i = 0; i < 100; i++) {
      if ( std::abs(fx) < ROL_EPSILON ) { break; }
      s   = -fx/smale_pdf(x);
      a   = 1.0;
      xs  = x + a*s;
      tmp = fx;
      fx  = smale_cdf(xs)-input;
      while ( std::abs(fx) > (1.0 - 1.e-4*a)*std::abs(tmp) ) {
        a *= 0.5;
        xs = x + a*s;
        fx = smale_cdf(xs)-input;
      }
      x = xs;
    }
    return x;
  }

  void test_onesided(Real x) {
    Real vx = cdf(x);
    Real vy = 0.0;
    Real dv = pdf(x);
    Real t = 1.0;
    Real diff = 0.0;
    Real err = 0.0;
    std::cout << std::scientific << std::setprecision(11);
    std::cout << std::right << std::setw(20) << "CHECK DENSITY: f(x) = cdf(x) with x = "
                                             << x << " is correct?\n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "f'(x)"
                            << std::setw(20) << "(f(x+t)-f(x))/t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = cdf(x+t);
      diff = (vy-vx)/t;
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= 0.1;
    }
    std::cout << "\n";
    // CHECK INTCDF
    vx = intcdf(x);
    vy = 0.0;
    dv = cdf(x);
    t = 1.0;
    diff = 0.0;
    err = 0.0;
    std::cout << std::scientific << std::setprecision(11);
    std::cout << std::right << std::setw(20) << "CHECK DENSITY: f(x) = intcdf(x) with x = "
                                             << x << " is correct?\n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "f'(x)"
                            << std::setw(20) << "(f(x+t)-f(x))/t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = intcdf(x+t);
      diff = (vy-vx)/t;
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= 0.1;
    }
    std::cout << "\n";
    // CHECK INVCDF
    vx = cdf(x);
    vy = invcdf(vx);
    err = std::abs(x-vy);
    std::cout << std::scientific << std::setprecision(11);
    std::cout << std::right << std::setw(20) << "CHECK DENSITY: f(x) = invcdf(x) with x = "
                                             << x << " is correct?\n";
    std::cout << std::right << std::setw(20) << "cdf(x)"
                            << std::setw(20) << "invcdf(cdf(x))"
                            << std::setw(20) << "Error"
                            << "\n";
    std::cout << std::scientific << std::setprecision(11) << std::right
              << std::setw(20) << vx
              << std::setw(20) << vy
              << std::setw(20) << err
              << "\n\n";
  }

  void test_centered(Real x) {
    Real vx = 0.0;
    Real vy = 0.0;
    Real dv = pdf(x);
    Real t = 1.0;
    Real diff = 0.0;
    Real err = 0.0;
    std::cout << std::scientific << std::setprecision(11);
    std::cout << std::right << std::setw(20) << "CHECK DENSITY: f(x) = cdf(x) with x = "
                                             << x << " is correct?\n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "f'(x)"
                            << std::setw(20) << "(f(x+t)-f(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vx = cdf(x+t);
      vy = cdf(x-t);
      diff = 0.5*(vx-vy)/t;
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= 0.1;
    }
    std::cout << "\n";
    // CHECK INTCDF
    vx = 0.0;
    vy = 0.0;
    dv = cdf(x);
    t = 1.0;
    diff = 0.0;
    err = 0.0;
    std::cout << std::scientific << std::setprecision(11);
    std::cout << std::right << std::setw(20) << "CHECK DENSITY: f(x) = intcdf(x) with x = "
                                             << x << " is correct?\n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "f'(x)"
                            << std::setw(20) << "(f(x+t)-f(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vx = intcdf(x+t);
      vy = intcdf(x-t);
      diff = 0.5*(vx-vy)/t;
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= 0.1;
    }
    std::cout << "\n";
    // CHECK INVCDF
    vx = cdf(x);
    vy = invcdf(vx);
    err = std::abs(x-vy);
    std::cout << std::scientific << std::setprecision(11);
    std::cout << std::right << std::setw(20) << "CHECK DENSITY: f(x) = invcdf(x) with x = "
                                             << x << " is correct?\n";
    std::cout << std::right << std::setw(20) << "cdf(x)"
                            << std::setw(20) << "invcdf(cdf(x))"
                            << std::setw(20) << "Error"
                            << "\n";
    std::cout << std::scientific << std::setprecision(11) << std::right
              << std::setw(20) << vx
              << std::setw(20) << vy
              << std::setw(20) << err
              << "\n\n";
  }
};

}

#endif
