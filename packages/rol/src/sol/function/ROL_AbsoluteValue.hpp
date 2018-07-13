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

#ifndef ROL_ABSOLUTEVALUE_HPP
#define ROL_ABSOLUTEVALUE_HPP

#include "ROL_Types.hpp"
#include "ROL_PositiveFunction.hpp"

namespace ROL {

enum EAbsoluteValue {
  ABSOLUTEVALUE_TRUE = 0,
  ABSOLUTEVALUE_SQUAREROOT,
  ABSOLUTEVALUE_SQRTDENOM,
  ABSOLUTEVALUE_C2,
  ABSOLUTEVALUE_LAST
};

template<class Real>
class AbsoluteValue : public PositiveFunction<Real> {
private:
  Real param_;
  EAbsoluteValue eav_;

public: 
  AbsoluteValue(Real param = 1., EAbsoluteValue eav = ABSOLUTEVALUE_TRUE)
    : param_(1./param), eav_(eav) {
    if ( eav != ABSOLUTEVALUE_TRUE && std::abs(param) < ROL_EPSILON<Real>() ) { param_ = 1.e2; }
  }

  AbsoluteValue(ROL::ParameterList &parlist) {
    Real param = parlist.get("Smoothing Parameter",1.);
    param_ = 1./((param > 0.) ? param : 1.);
    std::string type = parlist.get("Absolute Value Approximation","true");
    eav_ = ABSOLUTEVALUE_TRUE;
    if      ( type == "Square Root" )             { eav_ = ABSOLUTEVALUE_SQUAREROOT; }
    else if ( type == "Square Root Denominator" ) { eav_ = ABSOLUTEVALUE_SQRTDENOM;  }
    else if ( type == "C2")                       { eav_ = ABSOLUTEVALUE_C2;         }
    else if ( type == "true")                     { eav_ = ABSOLUTEVALUE_TRUE;       }
  }

  Real evaluate(Real input, int deriv) {
    Real val = 0.0;
    switch(eav_) {
      case ABSOLUTEVALUE_TRUE:       val = true_absolute_value(input,deriv);  break;
      case ABSOLUTEVALUE_SQUAREROOT: val = sqrt_absolute_value(input,deriv);  break;
      case ABSOLUTEVALUE_SQRTDENOM:  val = sqrtd_absolute_value(input,deriv); break;
      case ABSOLUTEVALUE_C2:         val = c2_absolute_value(input,deriv);    break;
      default:
        ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                          ">>> ERROR (ROL::AbsoluteValue): Absolute value approximation not defined!");
    }
    return val;
  }

private:
  Real true_absolute_value( Real input, int deriv ) {
    Real output = 0.0, e = 0.0;
    if ( std::abs(param_) > ROL_EPSILON<Real>() ) { e = 0.5/param_; }

    int region = 0;
    if ( input < -e )     { region = -1; }
    else if ( input > e ) { region = 1; }

    if ( deriv == 0 )      { output = std::abs(input); }
    else if ( deriv == 1 ) { output = (input < 0.0 ? -1.0   : 1.0); }
    else if ( deriv == 2 ) {
      if ( region == -1 )     { output = 0.0; }
      else if ( region == 0 ) { output = e*std::exp( -1.0/(1.0 - std::pow(std::abs(input/e),2.0)) ); }
      else if ( region == 1 ) { output = 0.0; }
    }
    return output;
  }

  Real sqrt_absolute_value( Real input, int deriv ) {
    Real output = 0.0;
    if ( deriv == 0 )      { output = std::sqrt(input*input + 1.0/param_); }
    else if ( deriv == 1 ) { output = input/std::sqrt(input*input+1.0/param_); }
    else if ( deriv == 2 ) { output = (1.0/param_)/std::pow(input*input+1.0/param_,1.5); }
    return output;
  }

  Real sqrtd_absolute_value( Real input, int deriv ) {
    Real output = 0.0;
    if ( deriv == 0 )      { output = input*input/std::sqrt(input*input + 1.0/param_); }
    else if ( deriv == 1 ) { output = (2.0/param_*input+std::pow(input,3.0)) /
                                    std::pow(input*input+1.0/param_,1.5); }
    else if ( deriv == 2 ) { output = ((2.0/param_-input*input)/param_) / 
                                    std::pow(input*input+1.0/param_,2.5); }
    return output;
  }

  Real c2_absolute_value( Real input, int deriv ) {
    Real output = 0.0, e = 1.0;
    if ( std::abs(param_) > ROL_EPSILON<Real>() ) { e = 0.5/param_; }

    int region = 0;
    if ( input < -e )     { region = -1; }
    else if ( input > e ) { region = 1; }

    if ( deriv == 0 )      { 
      if ( std::abs(region) == 1 ) { output = std::abs(input); }
      else if ( region == 0 ) { output = 1.875*std::pow(input*e,2.0) - 
                                         1.25 *std::pow(input*e,4.0) +
                                         0.375*std::pow(input*e,6.0); }
    }
    else if ( deriv == 1 ) { 
      if ( std::abs(region) == 1 ) { output = (input < 0.0 ? -1.0 : 1.0); }
      else if ( region == 0 ) { output = e*2.0*1.875*input*e - 
                                         e*4.0*1.25 *std::pow(input*e,3.0) +
                                         e*6.0*0.375*std::pow(input*e,5.0); }
    }
    else if ( deriv == 2 ) {
      if ( std::abs(region) == 1 ) { output = 0.0; }
      else if ( region == 0 ) { output =   e* 2.0*1.875 - 
                                         e*e*12.0*1.25 *std::pow(input*e,2.0) +
                                         e*e*30.0*0.375*std::pow(input*e,4.0); }
    }
    return output;
  }
};

}

#endif
