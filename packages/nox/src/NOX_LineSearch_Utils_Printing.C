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

#include "NOX_LineSearch_Utils_Printing.H"
#include "NOX_Utils.H"
#include "NOX_Common.H"

NOX::LineSearch::Utils::Printing::Printing(const NOX::Utils& u) :
  NOX::Utils(u)
{

}

NOX::LineSearch::Utils::Printing::~Printing()
{

}

void NOX::LineSearch::Utils::Printing::printStep(int n, double step, double oldf, double newf, const string s, bool unscaleF) const
{
  if (isPrintProcessAndType(NOX::Utils::InnerIteration)) 
  {
    cout << setw(3) << n << ":";
    cout << NOX::Utils::fill(1,' ') << "step = " << sciformat(step);
    if (unscaleF == true) {
      cout << NOX::Utils::fill(1,' ') << "oldf = " << sciformat(sqrt(2. * oldf));
      cout << NOX::Utils::fill(1,' ') << "newf = " << sciformat(sqrt(2. * newf));
    }
    else {
      cout << NOX::Utils::fill(1,' ') << "oldf = " << sciformat(oldf);
      cout << NOX::Utils::fill(1,' ') << "newf = " << sciformat(newf);
    }
    if (!s.empty()) 
    {
      cout << " " << s << "\n";
      cout << NOX::Utils::fill(72);
    }
    cout << endl;
  }
}
