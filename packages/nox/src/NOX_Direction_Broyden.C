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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
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

NOX::Direction::Broyden::BroydenMemory::~BroydenMemory()
{
}

void NOX::Direction::Broyden::BroydenMemory::reset(int m)
{
  mMax = m;

  if (memory.size() < (unsigned int) mMax)
    memory.resize(mMax);

  if (index.capacity() < (unsigned int) mMax)
    index.reserve(mMax);

  index.resize(0);
}

void NOX::Direction::Broyden::BroydenMemory::reset()
{
  index.resize(0);
}

void NOX::Direction::Broyden::BroydenMemory::push(const NOX::Abstract::Vector& d)
{
  // Do nothing if the Brodyen memory size is zero.
  if (mMax == 0)
    return;

  // Number of updates currently stored
  int m = index.size();
  
  // Adjust the index vector appropriately so that the last entry in
  // index points to the memory location to be used by the new update.
  if (m < mMax) // memory is not full - use memory location m
  {
    index.push_back(m);
  }
  else // memory is full, so recycle
  {
    int k = index[0];		// save the index of the oldest update
    for (int i = 0; i < m - 1; i ++)
      index[i] = index[i+1];
    index[m-1] = k;		// reuse the oldest update
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
  lsParamsPtr(NULL),
  oldJacobianGrpPtr(NULL),
  inexactNewtonUtils(u, p)
{
  reset(p);
}

NOX::Direction::Broyden::~Broyden()
{
  delete oldJacobianGrpPtr;
}

bool NOX::Direction::Broyden::reset(Parameter::List& params)
{
  NOX::Parameter::List&  p = params.sublist("Broyden");

  // Save a pointer to the Linear Solver sublist
  lsParamsPtr = &p.sublist("Linear Solver"); 

  // Set the default linear solver tolerance
  //if (!lsParamsPtr->isParameter("Tolerance"))
  //lsParamsPtr->getParameter("Tolerance", 1.0e-4);

  // Reset the inexact Newton Utilities (including linear solve tolerance)
  inexactNewtonUtils.reset(utils, params);

  // Get the restart frequency
  cntMax = p.getParameter("Restart Frequency", 10);

  // Get the maximum convergence rate
  maxConvRate = p.getParameter("Max Convergence Rate", 1.0);

  // Get the memory size
  memorySizeMax = p.getParameter("Memory", cntMax);
  
  // Reset the memory
  memory.reset(memorySizeMax);

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
  // Return value for group operations (temp variable)
  NOX::Abstract::Group::ReturnType status;
  
  // Compute F at current solution
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok) 
    throwError("compute", "Unable to compute F");

  // Check for restart
  if (doRestart(soln, solver))
  {
    // Reset memory
    memory.reset();

    // Update group
    if (oldJacobianGrpPtr == NULL)
      oldJacobianGrpPtr = soln.clone(NOX::DeepCopy);
    else
      // RPP - update the entire group (this grabs state vectors in xyce).
      // Otherwise, xyce is forced to recalculate F at each iteration.
      //oldJacobianGrpPtr->setX(soln.getX());
      *oldJacobianGrpPtr = soln;

    // Calcuate new Jacobian
    if (utils.isPrintProcessAndType(NOX::Utils::Details))
      cout << "       Recomputing Jacobian" << endl;
 
    status = oldJacobianGrpPtr->computeJacobian();
    if (status != NOX::Abstract::Group::Ok) 
      throwError("compute", "Unable to compute Jacobian");

    // Reset counter
    cnt = 0;
  }

  // If necesary, scale the s-vector from the last iteration
  if (!memory.empty()) 
  {
    double step = solver.getStepSize();
    memory[memory.size() - 1].setStep(step);
  }

  // --- Calculate the Broyden direction ---

  // Compute inexact forcing term if requested.
  inexactNewtonUtils.computeForcingTerm(soln, 
					solver.getPreviousSolutionGroup(),
					solver.getNumIterations(),
					solver);

  // dir = - J_old^{-1} * F
  cnt ++;
  status = oldJacobianGrpPtr->applyJacobianInverse(*lsParamsPtr, 
						   soln.getF(), 
						   dir);
  if (status != NOX::Abstract::Group::Ok) 
    throwError("compute", "Unable to apply Jacobian inverse");
  dir.scale(-1.0);

  // Apply the Broyden modifications to the old Jacobian (implicitly)
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

bool NOX::Direction::Broyden::doRestart(NOX::Abstract::Group& soln, 
					const NOX::Solver::LineSearchBased& solver)
{
  // Test 1 - First iteration!
  if (solver.getNumIterations() == 0)
    return true;

  // Test 2 - Frequency
  if (cnt >= cntMax)
    return true;

  // Test 3 - Last step was zero!
  if (solver.getStepSize() == 0.0)
    return true;

  // Test 4 - Check for convergence rate
  convRate = soln.getNormF() / solver.getPreviousSolutionGroup().getNormF();
  if (convRate > maxConvRate)
    return true;

  return false;
}


void NOX::Direction::Broyden::throwError(const string& functionName, const string& errorMsg)
{
    if (utils.isPrintProcessAndType(NOX::Utils::Error))
      cerr << "NOX::Direction::Broyden::" << functionName << " - " << errorMsg << endl;
    throw "NOX Error";
}




#endif

