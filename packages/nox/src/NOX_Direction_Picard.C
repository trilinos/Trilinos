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

#ifdef WITH_PRERELEASE

#include "NOX_Direction_Picard.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Parameter_List.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

/* Some compilers (in particular the SGI and ASCI Red - TFLOP)
 * fail to find the max and min function.  Therfore we redefine them
 * here.
 */
#ifdef max
#undef max
#endif

#define max(a,b) ((a)>(b)) ? (a) : (b);

#ifdef min
#undef min
#endif

#define min(a,b) ((a)<(b)) ? (a) : (b);

using namespace NOX;
using namespace NOX::Direction;

Picard::Picard(const NOX::Utils& u, Parameter::List& params) :
  utils(u),
  paramsPtr(NULL)
{
  reset(params);
}


bool Picard::reset(Parameter::List& params) 
{
  return true;
}

Picard::~Picard() 
{
}


bool Picard::compute(Abstract::Vector& dir, Abstract::Group& soln,
                          const Solver::Generic& solver)
{
  Abstract::Group::ReturnType ok;

  niter = solver.getNumIterations();

  // Construct Residual and negate to give new search direction

  ok = soln.computeF();
  if (ok != Abstract::Group::Ok) {
    if (utils.isPrintProcessAndType(Utils::Warning))
      cout << "NOX::Direction::Picard::compute - Unable to compute F." 
           << endl;
    return false;
  }
  dir = soln.getF();  

  dir.scale(-1.0);

  return (ok == Abstract::Group::Ok);
}

#endif
