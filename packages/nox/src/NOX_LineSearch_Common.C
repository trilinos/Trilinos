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

#include "NOX_LineSearch_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::LineSearch;

Common::Common(const NOX::Utils& u) :
  utils(u),
  vecPtr(NULL)
{
}

Common::~Common()
{
  delete vecPtr;
}

void Common::printStep(int n, double step, double oldf, double newf, const string s) const
{
  if (utils.isPrintProcessAndType(Utils::InnerIteration)) 
  {
    cout << setw(3) << n << ":";
    cout << Utils::fill(1,' ') << "step = " << utils.sciformat(step);
    cout << Utils::fill(1,' ') << "oldf = " << utils.sciformat(sqrt(2. * oldf));
    cout << Utils::fill(1,' ') << "newf = " << utils.sciformat(sqrt(2. * newf));
    if (!s.empty()) 
    {
      cout << " " << s << "\n";
      cout << Utils::fill(72);
    }
    cout << endl;
  }
}

double Common::computeSlope(const Abstract::Vector& dir, const Abstract::Group& grp) 
{
   if (grp.isGradient()) 
     return(dir.dot(grp.getGradient()));

  // Allocate space for vecPtr if necessary
  if (vecPtr == NULL) 
    vecPtr = dir.clone(ShapeCopy);

  // v = J * dir
  NOX::Abstract::Group::ReturnType status = grp.applyJacobian(dir,*vecPtr);
  
  if (status != NOX::Abstract::Group::Ok) 
  {
    cout << "NOX::LineSearch::Common::computeSlope -  Unable to apply Jacobian!" << endl;
    throw "NOX Error";
  }

  // Check that F exists
  if (!grp.isF()) 
  {
    cout << "NOX::LineSearch::Common::computeSlope - Invalid F" << endl;
    throw "NOX Error";
  }

  // Return <v, F> = F' * J * dir = <J'F, dir> = <g, dir>
  return(vecPtr->dot(grp.getF()));
}
