// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "NOX_Random.H"

NOX::Random::Random() 
{
  seed = static_cast<double>(rand());

  // rand() can return 0 or 2147483647, so adjust seed if that happens
  if ((seed == 0.0) || (seed == 2147483647.0))
    seed = 1.0;
}

NOX::Random::Random(int s)
{
  seed = checkSeed("Random", s);
}

void NOX::Random::setSeed(int  s)
{
  seed = checkSeed("setSeed", s);
}

double NOX::Random::operator() ()
{
  const double a = 16807.0;
  const double bigInt = 2147483647.0;

  seed = fmod(a*seed, bigInt);
  return 2.0*(seed/bigInt)-1.0;
}

double NOX::Random::checkSeed(const string& func, int s)
{
  if ((s < 1) || (s > 2147483646)) {
    cerr << "Error in NOX::Random::" << s << "():  " << "supplied seed " 
	 << s << " is not an integer between 1 and 2147483646." << endl
	 << "Using a seed of 1 instead." << endl;
    return 1.0;
  }
  else
    return static_cast<double>(s);
}
