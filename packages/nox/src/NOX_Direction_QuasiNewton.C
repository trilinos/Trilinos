#ifdef WITH_PRERELEASE
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

#include "NOX_Direction_QuasiNewton.H" // class definition
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::Direction;

//------------------------------------------------------------

QuasiNewton::Update::Update(const Abstract::Vector& sin, const Abstract::Vector& yin) :
  sptr(sin.clone(DeepCopy)), 
  sref(*sptr),
  yptr(yin.clone(DeepCopy)), 
  yref(*yptr),
  rho(sref.dot(yref))
{
}

QuasiNewton::Update::~Update() 
{
  delete sptr;
  delete yptr;
}

void QuasiNewton::Update::reset(const Abstract::Vector& sin, const Abstract::Vector& yin) 
{
  sref = sin;
  yref = yin;
  rho = sref.dot(yref);
}

const Abstract::Vector& QuasiNewton::Update::s() const
{
  return sref;
}

const Abstract::Vector& QuasiNewton::Update::y() const
{
  return yref;
}

double QuasiNewton::Update::sdoty() const
{
  return rho;
}

//------------------------------------------------------------

QuasiNewton::Updates::Updates(int m) :
  maxSize(m)
{
}

QuasiNewton::Updates::~Updates()
{
  Update* u;
  deque<Update*>::iterator i;

  // Clean out update deque
  for (i = updateDeque.begin(); i != updateDeque.end(); ++i) {
    delete *i;
    *i = NULL;
  }

  // Clean out update deque
  for (i = recycleDeque.begin(); i != recycleDeque.end(); ++i) {
    delete *i;
    *i = NULL;
  }
  
}

void QuasiNewton::Updates::reset(int m)
{
  // Copy update pointer elements onto unused update deque
  copy(updateDeque.begin(), updateDeque.end(), recycleDeque.end());

  // Erase updateDeque
  updateDeque.clear();

  // Reset the max size
  maxSize = m;
}

void QuasiNewton::Updates::push_back(const Abstract::Vector& s, const Abstract::Vector& y)
{
  Update* updatePtr;

  // If the updateDeque is full, recycle and delete the first entry
  if (updateDeque.size() == maxSize) {
    recycleDeque.push_back(updateDeque.front());
    updateDeque.pop_front();
  }

  // Recycle or create a new entry
  if (!recycleDeque.empty()) {
    updatePtr = recycleDeque.front();
    recycleDeque.pop_front();
    updatePtr->reset(s,y);
  }
  else 
    updatePtr = new Update(s,y);

  // Push the new (or recycled) update onto the end of the deque.
  updateDeque.push_back(updatePtr);
}

QuasiNewton::UpdateConstIterator QuasiNewton::Updates::begin() const
{
  return updateDeque.begin();
}

QuasiNewton::UpdateConstIterator QuasiNewton::Updates::end() const
{
  return updateDeque.end();
}

QuasiNewton::UpdateConstReverseIterator QuasiNewton::Updates::rbegin() const
{
  return updateDeque.rbegin();
}

QuasiNewton::UpdateConstReverseIterator QuasiNewton::Updates::rend() const
{
  return updateDeque.rend();
}

bool QuasiNewton::Updates::empty() const
{
  return updateDeque.empty();
}

//------------------------------------------------------------

QuasiNewton::QuasiNewton(Parameter::List& p) :
  paramsPtr(NULL),
  stepDirPtr(NULL),
  nMemory(0),
  updates()
{
  reset(p);
}

QuasiNewton::~QuasiNewton()
{
  // Empty out the memory deque
  /*
  deque<tuplet*>::iterator i;
  for (i = memory.begin(); i != memory.end(); ++i) 
  {
    delete i->sPtr;
    delete i->yPtr;
  }
  */
}

bool QuasiNewton::reset(Parameter::List& p)
{
  paramsPtr = &p;
  return true;
}

bool QuasiNewton::compute(Abstract::Vector& dir, 
			   Abstract::Group& soln, 
			   const Solver::Generic& solver)
{
  // Compute F at current solution
  bool ok = soln.computeF();
  double normF = soln.getNormF();

  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::QuasiNewton::compute - Unable to compute F." << endl;
    return false;
  }

  // Compute Jacobian at current solution.
  ok = soln.computeJacobian();

  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::QuasiNewton::compute - Unable to compute Jacobian." << endl;
    return false;
  }
  
  // Compute the QuasiNewton direction
  
  return ok;
}


#endif
