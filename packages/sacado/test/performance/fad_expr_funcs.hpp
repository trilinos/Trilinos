// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef FAD_EXPR_FUNCS_HPP
#define FAD_EXPR_FUNCS_HPP

#include "Sacado_MathFunctions.hpp"
#include "Sacado_Fad_DFad.hpp"
#include "Sacado_ELRFad_DFad.hpp"
#include "Sacado_CacheFad_DFad.hpp"
#include "Sacado_ELRCacheFad_DFad.hpp"

template <typename T> 
struct ExprFuncs {
  static const int nfunc = 8;
  static const char* mult_names[nfunc];
  static const char* nest_names[nfunc];
  static const char* add_names[nfunc];
  static const int nx_max = 21;

  static void mult0 (const T x[], T& y);
  static void mult1 (const T x[], T& y);
  static void mult2 (const T x[], T& y);
  static void mult3 (const T x[], T& y);
  static void mult4 (const T x[], T& y);
  static void mult5 (const T x[], T& y);
  static void mult10(const T x[], T& y);
  static void mult15(const T x[], T& y);
  static void mult20(const T x[], T& y);

  static void nest0 (const T x[], T& y);
  static void nest1 (const T x[], T& y);
  static void nest2 (const T x[], T& y);
  static void nest3 (const T x[], T& y);
  static void nest4 (const T x[], T& y);
  static void nest5 (const T x[], T& y);
  static void nest10(const T x[], T& y);
  static void nest15(const T x[], T& y);
  static void nest20(const T x[], T& y);

  static void add0 (const T x[], T& y);
  static void add1 (const T x[], T& y);
  static void add2 (const T x[], T& y);
  static void add3 (const T x[], T& y);
  static void add4 (const T x[], T& y);
  static void add5 (const T x[], T& y);
  static void add10(const T x[], T& y);
  static void add15(const T x[], T& y);
  static void add20(const T x[], T& y);
};

#endif
