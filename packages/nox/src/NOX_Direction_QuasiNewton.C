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

#include "NOX_Common.H"

#include "NOX_Direction_QuasiNewton.H" 
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::Direction;

//------------------------------------------------------------

QuasiNewton::Update::Update(const Abstract::Vector& newX, const Abstract::Vector& oldX, 
			    const Abstract::Vector& newG, const Abstract::Vector& oldG) :
  sptr(newX.clone(ShapeCopy)), 
  yptr(newX.clone(ShapeCopy))
{
  reset(newX, oldX, newG, oldG);
}

QuasiNewton::Update::~Update() 
{
  delete sptr;
  delete yptr;
}

void QuasiNewton::Update::reset(const Abstract::Vector& newX, const Abstract::Vector& oldX, 
				const Abstract::Vector& newG, const Abstract::Vector& oldG)
{
  sptr->update(1.0, newX, -1.0, oldX, 0.0);
  yptr->update(1.0, newG, -1.0, oldG, 0.0);
  rho = sptr->dot(*yptr);
}

const Abstract::Vector& QuasiNewton::Update::s() const
{
  return *sptr;
}

const Abstract::Vector& QuasiNewton::Update::y() const
{
  return *yptr;
}

double QuasiNewton::Update::sdoty() const
{
  return rho;
}

//------------------------------------------------------------

QuasiNewton::Updates::Updates(int m) :
  maxUpdates(m)
{
}

QuasiNewton::Updates::~Updates()
{
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
  maxUpdates = m;
}

void QuasiNewton::Updates::push_back(const Abstract::Vector& newX, const Abstract::Vector& oldX, 
				     const Abstract::Vector& newG, const Abstract::Vector& oldG)
{
  Update* updatePtr;

  // If the updateDeque is full, recycle and delete the first entry
  if (updateDeque.size() == maxUpdates) {
    recycleDeque.push_back(updateDeque.front());
    updateDeque.pop_front();
  }

  // Recycle or create a new entry
  if (!recycleDeque.empty()) {
    updatePtr = recycleDeque.front();
    recycleDeque.pop_front();
    updatePtr->reset(newX, oldX, newG, oldG);
  }
  else 
    updatePtr = new Update(newX, oldX, newG, oldG);

  // Push the new (or recycled) update onto the end of the deque.
  updateDeque.push_back(updatePtr);
}

const QuasiNewton::Update* QuasiNewton::Updates::back() const
{
  return updateDeque.back();
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
  paramsPtr(NULL)
{
  reset(p);
}

QuasiNewton::~QuasiNewton()
{
}

bool QuasiNewton::reset(Parameter::List& p)
{
  paramsPtr = &p;
  updates.reset(p.getParameter("Memory", 5));
  return true;
}

bool QuasiNewton::compute(Abstract::Vector& dir, 
			   Abstract::Group& soln, 
			   const Solver::Generic& solver)
{
  NOX::Abstract::Group::ReturnType status;
  
  // Compute F at current solution
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok) 
    throwError("compute", "Unable to compute F");

  // Compute Jacobian at current solution.
  status = soln.computeJacobian();
  if (status != NOX::Abstract::Group::Ok) 
    throwError("compute", "Unable to compute Jacobian");

  // Compute the gradient at the current solution
  status = soln.computeGradient();
  if (status != NOX::Abstract::Group::Ok) 
    throwError("compute", "Unable to compute gradient");

  // Push the old information onto the updates, but only after at least one previous iteration
  if (solver.getNumIterations() > 0) 
  {
    const Abstract::Group& oldSoln = solver.getPreviousSolutionGroup();
    if (oldSoln.isGradient())
      updates.push_back(soln.getX(), oldSoln.getX(), soln.getGradient(), oldSoln.getGradient());
  }

  // Calculate the QN direction
  dir = soln.getGradient();

  if (!updates.empty()) {

    deque<double> alpha;
    double a, b;
  
    for (UpdateConstReverseIterator i = updates.rbegin(); i != updates.rend(); i++) {

      const Update& u = *(*i);
      const double rho = 1.0 / u.sdoty();
      const Abstract::Vector& s = u.s();
      const Abstract::Vector& y = u.y();
      a = rho * dir.dot(s);
      dir.update(-1.0 * a, y, 1.0);
      alpha.push_front(a);
      
    }

    const Update& u = *(updates.back());
    const Abstract::Vector& y = u.y();
    double gamma = u.sdoty() / y.dot(y);
    dir.scale(gamma);

    for (UpdateConstIterator i = updates.begin(); i != updates.end(); i++) {
      
      const Update& u = *(*i);
      const double rho = 1.0 / u.sdoty();
      const Abstract::Vector& s = u.s();
      const Abstract::Vector& y = u.y();
      a = alpha.front();
      b = rho * dir.dot(y);
      dir.update(a - b, s, 1.0);
      alpha.pop_front();
      
    }
  }

  dir.scale(-1.0);

  return true;
}

void NOX::Direction::QuasiNewton::throwError(const string& functionName, const string& errorMsg)
{
    if (Utils::doPrint(Utils::Error))
      cerr << "NOX::Direction::QuasiNewton::" << functionName << " - " << errorMsg << endl;
    throw "NOX Error";
}
