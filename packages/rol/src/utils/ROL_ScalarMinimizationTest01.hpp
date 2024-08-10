// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCALARMINIMIZATIONTEST01_H
#define ROL_SCALARMINIMIZATIONTEST01_H

/** \class ROL::ScalarMinimizationTest01
    \brief Implements test scalar function 01.
*/

#include "ROL_ScalarMinimizationTest.hpp"

namespace ROL { 

template<class Real>
class ScalarMinimizationTest01 : public ScalarMinimizationTest<Real> {
private:
  std::vector<Real> fvector_;
  std::vector<Real> xvector_;

  class ScalarFunctionTest01 : public ScalarFunction<Real> {
  public:
    Real value(const Real x) {
      Real zero(0), two(2), five(5);
      Real val = zero, I = zero;
      for (int i = 0; i < 20; i++) {
        I = (Real)(i+1);
        val += std::pow((two*I - five)/(x-(I*I)),two);
      }
      return val;
    }
    Real deriv(const Real x) {
      Real zero(0), two(2), three(3), five(5);
      Real val = zero, I = zero;
      for (int i = 0; i < 20; i++) {
        I = (Real)(i+1);
        val += std::pow((two*I - five),two)/std::pow((x-(I*I)),three);
      }
      return -two*val;
    }
  };

public:
  ScalarMinimizationTest01(ROL::ParameterList &parlist)
    : ScalarMinimizationTest<Real>(parlist) {
    fvector_.clear(); fvector_.resize(19,0);
    xvector_.clear(); xvector_.resize(19,0);

    fvector_[0]  = 3.6766990169; xvector_[0]  =   3.0229153;
    fvector_[1]  = 1.1118500100; xvector_[1]  =   6.6837536;
    fvector_[2]  = 1.2182217637; xvector_[2]  =  11.2387017;
    fvector_[3]  = 2.1621103109; xvector_[3]  =  19.6760001;
    fvector_[4]  = 3.0322905193; xvector_[4]  =  29.8282273;
    fvector_[5]  = 3.7583856477; xvector_[5]  =  41.9061162;
    fvector_[6]  = 4.3554103836; xvector_[6]  =  55.9535958;
    fvector_[7]  = 4.8482959563; xvector_[7]  =  71.9856656;
    fvector_[8]  = 5.2587585400; xvector_[8]  =  90.0088685;
    fvector_[9]  = 5.6036524295; xvector_[9]  = 110.0265327;
    fvector_[10] = 5.8956037976; xvector_[10] = 132.0405517;
    fvector_[11] = 6.1438861542; xvector_[11] = 156.0521144;
    fvector_[12] = 6.3550764593; xvector_[12] = 182.0620604;
    fvector_[13] = 6.5333662003; xvector_[13] = 210.0711010;
    fvector_[14] = 6.6803639849; xvector_[14] = 240.0800483;
    fvector_[15] = 6.7938538365; xvector_[15] = 272.0902669;
    fvector_[16] = 6.8634981053; xvector_[16] = 306.1051233;
    fvector_[17] = 6.8539024631; xvector_[17] = 342.1369454;
    fvector_[18] = 6.6008470481; xvector_[18] = 380.2687097;
  }

  bool test(std::ostream &stream = std::cout) {
    Real fx = 0, x = 0, A = 0, B = 0;
    int nfval = 0, ngrad = 0;
    bool flag = true;
    ROL::Ptr<ScalarFunction<Real> > f
      = ROL::makePtr<ScalarFunctionTest01>();

    stream << "\nScalar Minimization Test 01\n";
    stream << "  ";
    stream << std::setw(10) << std::left << "lower";
    stream << std::setw(10) << std::left << "upper";
    stream << std::setw(15) << std::left << "x-error";
    stream << std::setw(15) << std::left << "f(x)-error";
    stream << std::setw(10) << std::left << "nfval";
    stream << std::setw(10) << std::left << "ngrad";
    stream << "\n";
    for (int i = 0; i < 19; i++) {
      nfval = 0; ngrad = 0;
      A = (Real)((i+1)*(i+1));
      B = (Real)((i+2)*(i+2));
      ScalarMinimizationTest<Real>::run(fx,x,nfval,ngrad,*f,A,B);
      flag &= std::abs(x-xvector_[i]) < 10.0*std::sqrt(ROL_EPSILON<Real>())*xvector_[i];

      stream << "  ";
      stream << std::setw(10) << std::left << (int)A;
      stream << std::setw(10) << std::left << (int)B;
      stream << std::scientific << std::setprecision(6);
      stream << std::setw(15) << std::left << std::abs( x-xvector_[i])/xvector_[i];
      stream << std::setw(15) << std::left << std::abs(fx-fvector_[i])/fvector_[i];
      stream << std::setw(10) << std::left << nfval;
      stream << std::setw(10) << std::left << ngrad;
      stream << "\n";
    }
    stream << "\n";
    return flag;
  }
};

}

#endif
