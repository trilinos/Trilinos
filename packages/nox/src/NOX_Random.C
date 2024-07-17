// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Random.H"
#include <cstdlib>
#include <iostream>
#include <cmath>

// instantiate static data
double NOX::Random::seed = 1.0;

NOX::Random::Random()
{
  seed = static_cast<double>(std::rand());

  // rand() can return 0 or 2147483647, so adjust seed if that happens
  if ((seed == 0.0) || (seed == 2147483647.0))
    seed = 1.0;
}

NOX::Random::Random(int s)
{
  setSeed(s);
}

void NOX::Random::setSeed(int  s)
{
  int ss = checkSeed("setSeed", s);
  std::srand(ss);
  seed = static_cast<double>(s);
}

double NOX::Random::number()
{
  const double a = 16807.0;
  const double bigInt = 2147483647.0;

  seed = std::fmod(a*seed, bigInt);
  return 2.0*(seed/bigInt)-1.0;
}

int NOX::Random::checkSeed(const std::string& /* func */, int s)
{
  if ((s < 1) || (s > 2147483646)) {
    std::cerr << "Error in NOX::Random::" << s << "():  " << "supplied seed "
     << s << " is not an integer between 1 and 2147483646." << std::endl
     << "Using a seed of 1 instead." << std::endl;
    return 1;
  }
  else
    return s;
}
