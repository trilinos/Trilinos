// $Id$ 
// $Source$ 
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

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstdlib>

#include "Sacado_Random.hpp"

Sacado::Random::Random(double a_, double b_) :
  a(a_),
  b(b_),
  seed(static_cast<double>(rand()))
{
  // rand() can return 0 or 2147483647, so adjust seed if that happens
  if ((seed == 0.0) || (seed == 2147483647.0))
    seed = 1.0;
}

Sacado::Random::Random(double a_, double b_, int s) :
  a(a_),
  b(b_),
  seed(0.0)
{
  setSeed(s);
}

Sacado::Random::~Random()
{
}

void
Sacado::Random::setSeed(int s) {
  int ss = checkSeed("setSeed", s);
  srand(ss);
  seed = static_cast<double>(s);
}

double
Sacado::Random::number() {
  const double A = 16807.0;
  const double bigInt = 2147483647.0;
      
  seed = std::fmod(A*seed, bigInt);
  return (b-a)*(seed/bigInt) + a;
}

int
Sacado::Random::checkSeed(const std::string& func, int s) {
  if ((s < 1) || (s > 2147483646)) {
    std::cerr << "Error in Sacado::Random::" << s << "():  " 
	      << "supplied seed " 
	      << s << " is not an integer between 1 and 2147483646." 
	      << std::endl << "Using a seed of 1 instead." << std::endl;
    return 1;
  }
  else
    return s;
}
