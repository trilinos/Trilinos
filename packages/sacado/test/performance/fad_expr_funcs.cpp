// @HEADER
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software]; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation]; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY]; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library]; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// @HEADER

#include "fad_expr_funcs.hpp"

template <typename T>
void ExprFuncs<T>::mult1(const T x[], T& y) {
  y = 
    x[0]*x[1];
}

template <typename T>
void ExprFuncs<T>::mult2(const T x[], T& y) {
  y = 
    x[0]*x[1]*x[2];
}

template <typename T>
void ExprFuncs<T>::mult3(const T x[], T& y) {
  y = 
    x[0]*x[1]*x[2]*x[3];
}

template <typename T>
void ExprFuncs<T>::mult4(const T x[], T& y) {
  y = 
    x[0]*x[1]*x[2]*x[3]*x[4];
}

template <typename T>
void ExprFuncs<T>::mult5(const T x[], T& y) {
  y = 
    x[0]*x[1]*x[2]*x[3]*x[4]*x[5];
}

template <typename T>
void ExprFuncs<T>::mult10(const T x[], T& y) {
  y = 
    x[0]*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]*
    x[8]*x[9]*x[10];
}

template <typename T>
void ExprFuncs<T>::mult15(const T x[], T& y) {
  y = 
    x[0]*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]*
    x[8]*x[9]*x[10]*x[11]*x[12]*x[13]*x[14]*x[15];
}

template <typename T>
void ExprFuncs<T>::mult20(const T x[], T& y) {
  y = 
    x[0]*x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]*
    x[8]*x[9]*x[10]*x[11]*x[12]*x[13]*x[14]*x[15]*
    x[16]*x[17]*x[18]*x[19]* x[20];
}

template <typename T>
void ExprFuncs<T>::nest1(const T x[], T& y) {
  y = std::sin(x[0]);
}

template <typename T>
void ExprFuncs<T>::nest2(const T x[], T& y) {
  y = std::sin(std::sin(x[0]));
}

template <typename T>
void ExprFuncs<T>::nest3(const T x[], T& y) {
  y = std::sin(std::sin(std::sin(x[0])));
}

template <typename T>
void ExprFuncs<T>::nest4(const T x[], T& y) {
  y = std::sin(std::sin(std::sin(std::sin(x[0]))));
}

template <typename T>
void ExprFuncs<T>::nest5(const T x[], T& y) {
  y = std::sin(std::sin(std::sin(std::sin(std::sin(x[0])))));
}

template <typename T>
void ExprFuncs<T>::nest10(const T x[], T& y) {
  y = std::sin(std::sin(std::sin(std::sin(std::sin(
    std::sin(std::sin(std::sin(std::sin(std::sin(x[0]))))))))));
}

template <typename T>
void ExprFuncs<T>::nest15(const T x[], T& y) {
  y = std::sin(std::sin(std::sin(std::sin(std::sin(
    std::sin(std::sin(std::sin(std::sin(std::sin(
    std::sin(std::sin(std::sin(std::sin(std::sin(x[0])))))))))))))));
}

template <typename T>
void ExprFuncs<T>::nest20(const T x[], T& y) {
  y = std::sin(std::sin(std::sin(std::sin(std::sin(
    std::sin(std::sin(std::sin(std::sin(std::sin(
    std::sin(std::sin(std::sin(std::sin(std::sin(
    std::sin(std::sin(std::sin(std::sin(std::sin(x[0]))))))))))))))))))));
}

template <typename T>
void ExprFuncs<T>::add1(const T x[], T& y) {
  y = 
    x[0]+x[1];
}

template <typename T>
void ExprFuncs<T>::add2(const T x[], T& y) {
  y = 
    x[0]+x[1]+x[2];
}

template <typename T>
void ExprFuncs<T>::add3(const T x[], T& y) {
  y = 
    x[0]+x[1]+x[2]+x[3];
}

template <typename T>
void ExprFuncs<T>::add4(const T x[], T& y) {
  y = 
    x[0]+x[1]+x[2]+x[3]+x[4];
}

template <typename T>
void ExprFuncs<T>::add5(const T x[], T& y) {
  y = 
    x[0]+x[1]+x[2]+x[3]+x[4]+x[5];
}

template <typename T>
void ExprFuncs<T>::add10(const T x[], T& y) {
  y = 
    x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+
    x[8]+x[9]+x[10];
}

template <typename T>
void ExprFuncs<T>::add15(const T x[], T& y) {
  y = 
    x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+
    x[8]+x[9]+x[10]+x[11]+x[12]+x[13]+x[14]+x[15];
}

template <typename T>
void ExprFuncs<T>::add20(const T x[], T& y) {
  y = 
    x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+
    x[8]+x[9]+x[10]+x[11]+x[12]+x[13]+x[14]+x[15]+
    x[16]+x[17]+x[18]+x[19]+x[20];
}

template <typename T> 
const char* ExprFuncs<T>::mult_names[ExprFuncs<T>::nfunc] = { 
  "mult1",
  "mult2",
  "mult3",
  "mult4",
  "mult5",
  "mult10",
  "mult15",
  "mult20"
};

template <typename T> 
const char* ExprFuncs<T>::nest_names[ExprFuncs<T>::nfunc] = { 
  "nest1",
  "nest2",
  "nest3",
  "nest4",
  "nest5",
  "nest10",
  "nest15",
  "nest20"
};

template <typename T> 
const char* ExprFuncs<T>::add_names[ExprFuncs<T>::nfunc] = { 
  "add1",
  "add2",
  "add3",
  "add4",
  "add5",
  "add10",
  "add15",
  "add20"
};

template struct ExprFuncs< double >;
template struct ExprFuncs< Sacado::Fad::DFad<double> >;
template struct ExprFuncs< Sacado::ELRFad::DFad<double> >;
template struct ExprFuncs< Sacado::CacheFad::DFad<double> >;
template struct ExprFuncs< Sacado::ELRCacheFad::DFad<double> >;



