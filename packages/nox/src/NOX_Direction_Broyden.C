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

#include "NOX_Common.H"

#include "NOX_Direction_Broyden.H" 
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_LineSearchBased.H"
#include "NOX_Utils.H"

//------------------------------------------------------------

NOX::Direction::Broyden::BroydenMemoryUnit::BroydenMemoryUnit() 
{
  sptr = NULL;
  lambda = 0;
  snormsqr = 0;
}

NOX::Direction::Broyden::BroydenMemoryUnit::~BroydenMemoryUnit() 
{
  delete sptr;
}

void NOX::Direction::Broyden::BroydenMemoryUnit::reset(const NOX::Abstract::Vector& d)
{
  // Copy d into s
  if (sptr == NULL)
    sptr = d.clone(DeepCopy);
  else
    *sptr = d;

  // snormsqr = |s|^2
  snormsqr = sptr->norm();
  snormsqr = snormsqr * snormsqr;

  lambda = 0.0;
}

void NOX::Direction::Broyden::BroydenMemoryUnit::setStep(double step)
{
  lambda = step;
  if (step != 1.0)
  {
    sptr->scale(step);
    snormsqr = step * step * snormsqr;
  }
}

const NOX::Abstract::Vector* NOX::Direction::Broyden::BroydenMemoryUnit::sPtr() const
{
  return sptr;
}

double NOX::Direction::Broyden::BroydenMemoryUnit::step() const
{
  return lambda;
}

double NOX::Direction::Broyden::BroydenMemoryUnit::sNormSqr() const
{
  return snormsqr;
}

//------------------------------------------------------------

NOX::Direction::Broyden::BroydenMemory::BroydenMemory() 
{
}

void NOX::Direction::Broyden::BroydenMemory::reset(NOX::Parameter::List& p)
{
  doRestart = p.getParameter("Memory Restart", false);

  int m = p.getParameter("Memory", 5);

  if (memory.size() < m)
    memory.resize(m);

  if (index.capacity() < m)
    index.reserve(m);

  index.resize(0);
}

NOX::Direction::Broyden::BroydenMemory::~BroydenMemory()
{
}

void NOX::Direction::Broyden::BroydenMemory::push(const NOX::Abstract::Vector& d)
{
  // Number of updates currently stored
  int m = index.size();
  
  // Adjust the index vector appropriately so that the last entry in
  // index points to the memory location to be used by the new update.

  if (m < memory.size()) // memory is not full - use memory location m
  {
    index.push_back(m);
  }
  else // memory is full
  {
    if (doRestart) // scrap the old updates and use the memory location 0
    {
      index.resize(0);
      index.push_back(0);
    }
    else // save the most recent m-1 updates and use the memory location that was freed
    {
      int k = index[0];		// save the index of the oldest update
      for (int i = 0; i < m - 1; i ++)
	index[i] = index[i+1];
      index[m-1] = k;		// reuse the oldest update
    }
  }

  // Save the vector d
  memory[index.back()].reset(d);
}

bool NOX::Direction::Broyden::BroydenMemory::empty() const
{
  return index.empty();
}

int NOX::Direction::Broyden::BroydenMemory::size() const
{
  return index.size();
}

NOX::Direction::Broyden::BroydenMemoryUnit& 
NOX::Direction::Broyden::BroydenMemory::operator[](int i)
{
  return memory[index[i]];
}

//------------------------------------------------------------

NOX::Direction::Broyden::Broyden(const NOX::Utils& u, Parameter::List& p) :
  utils(u),
  paramsPtr(NULL)
{
  reset(p);
}

NOX::Direction::Broyden::~Broyden()
{
}

bool NOX::Direction::Broyden::reset(Parameter::List& params)
{
  paramsPtr = &params;
  NOX::Parameter::List& p = params.sublist("Broyden");
  doComputeJacobian = p.getParameter("Compute Jacobian", false);
  doRightPreconditioning = p.getParameter("Apply Right Linear Preconditioning", false);
  if (doRightPreconditioning)
    precParamsPtr = &(p.sublist("Preconditioning"));    
  memory.reset(p);
  return true;
}

bool NOX::Direction::Broyden::compute(NOX::Abstract::Vector& dir, 
				      NOX::Abstract::Group& soln, 
				      const NOX::Solver::Generic& solver)
{
  throwError("compute", "This direction can only be used with a line search based solver.");
  return false;
}

bool NOX::Direction::Broyden::compute(NOX::Abstract::Vector& dir, 
				      NOX::Abstract::Group& soln, 
				      const NOX::Solver::LineSearchBased& solver)
{
  NOX::Abstract::Group::ReturnType status;
  
  // Scale the s-vector from the last iteration
  if (!memory.empty()) 
  {
    double step = solver.getStepSize();
    memory[memory.size() - 1].setStep(step);
  }

  // Compute F at current solution
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok) 
    throwError("compute", "Unable to compute F");

  // Optionally calculate the Jacobian
  if (doComputeJacobian)
  {
    status = soln.computeJacobian();
    if (status != NOX::Abstract::Group::Ok) 
      throwError("compute", "Unable to compute Jacobian");
  }

  // --- Calculate the Broyden direction ---
  
  if (doRightPreconditioning)	// d = - M * F
  {
    status = soln.applyRightPreconditioning(*precParamsPtr, soln.getF(), dir);
    if (status != NOX::Abstract::Group::Ok) 
      throwError("compute", "Unable to apply right preconditioning");
    dir.scale(-1.0);
   }
  else				// d = -F
    dir.update(-1.0, soln.getF(), 0.0);

  if (!memory.empty()) 
  {
    // Number of elements in the memory
    int m = memory.size();

    // Information corresponding to index i
    double step;
    const NOX::Abstract::Vector* sPtr;

    // Information corresponding to index i + 1 
    // (initialized for i = -1)
    double stepNext = memory[0].step();
    const NOX::Abstract::Vector* sPtrNext = memory[0].sPtr();

    // Intermediate storage
    double a, b, c, denom;

    for (int i = 0; i < m-1; i ++)
    {
      step = stepNext;
      sPtr = sPtrNext;
      stepNext = memory[i+1].step();
      sPtrNext = memory[i+1].sPtr();

      a = step / stepNext;
      b = step - 1;
      c = sPtr->dot(dir) / memory[i].sNormSqr();

      //cout << "a=" << a << " b=" << b << " c=" << c << endl;

      dir.update(a * c, *sPtrNext, b * c, *sPtr, 1.0);
    }

    step = stepNext;
    sPtr = sPtrNext;

    a = sPtr->dot(dir);		// <s,z>
    b = memory[m-1].sNormSqr();	// ||s||^2
    c = (step - 1) * a;		// (\lambda-1) <s,z>
    denom = b - step * a;	// ||s||^2 - \lambda <s,z>

    dir.update(c / denom, *sPtr, b / denom); 
  }

  //! Add this direction to the memory
  memory.push(dir);

  return true;
}

void NOX::Direction::Broyden::throwError(const string& functionName, const string& errorMsg)
{
    if (utils.isPrintProcessAndType(Utils::Error))
      cerr << "NOX::Direction::Broyden::" << functionName << " - " << errorMsg << endl;
    throw "NOX Error";
}


#endif

