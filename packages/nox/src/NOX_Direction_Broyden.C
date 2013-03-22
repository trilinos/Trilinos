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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_Common.H"

#include "NOX_Direction_Broyden.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_LineSearchBased.H"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"

//------------------------------------------------------------

NOX::Direction::Broyden::BroydenMemoryUnit::BroydenMemoryUnit() 
{
  lambda = 0;
  snormsqr = 0;
}

NOX::Direction::Broyden::BroydenMemoryUnit::~BroydenMemoryUnit() 
{

}

void NOX::Direction::Broyden::BroydenMemoryUnit::
reset(const NOX::Abstract::Vector& d)
{
  // Copy d into s
  if (Teuchos::is_null(sptr))
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

Teuchos::RCP<const NOX::Abstract::Vector> NOX::Direction::Broyden::
BroydenMemoryUnit::sPtr() const
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

void NOX::Direction::Broyden::BroydenMemory::
push(const NOX::Abstract::Vector& d)
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

NOX::Direction::Broyden::
Broyden(const Teuchos::RCP<NOX::GlobalData>& gd, 
	Teuchos::ParameterList& p) :
  lsParamsPtr(NULL),
  inexactNewtonUtils(gd, p)
{
  reset(gd, p);
}

NOX::Direction::Broyden::~Broyden()
{

}

bool NOX::Direction::Broyden::
reset(const Teuchos::RCP<NOX::GlobalData>& gd, 
      Teuchos::ParameterList& params)
{
  globalDataPtr = gd;
  utils = gd->getUtils();

  Teuchos::ParameterList&  p = params.sublist("Broyden");

  // Save a pointer to the Linear Solver sublist
  lsParamsPtr = &p.sublist("Linear Solver"); 

  // Set the default linear solver tolerance
  //if (!lsParamsPtr->isParameter("Tolerance"))
  //lsParamsPtr->get("Tolerance", 1.0e-4);

  // Reset the inexact Newton Utilities (including linear solve tolerance)
  inexactNewtonUtils.reset(gd, params);

  // Get the restart frequency
  cntMax = p.get("Restart Frequency", 10);

  // Get the maximum convergence rate
  maxConvRate = p.get("Max Convergence Rate", 1.0);

  // Get the memory size
  memorySizeMax = p.get("Memory", cntMax);
  
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
    if (Teuchos::is_null(oldJacobianGrpPtr))
      oldJacobianGrpPtr = soln.clone(NOX::DeepCopy);
    else
      // RPP - update the entire group (this grabs state vectors in xyce).
      // Otherwise, xyce is forced to recalculate F at each iteration.
      //oldJacobianGrpPtr->setX(soln.getX());
      *oldJacobianGrpPtr = soln;

    // Calcuate new Jacobian
    if (utils->isPrintType(NOX::Utils::Details))
      utils->out() << "       Recomputing Jacobian" << std::endl;
 
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
    Teuchos::RCP<const NOX::Abstract::Vector> sPtr;

    // Information corresponding to index i + 1 
    // (initialized for i = -1)
    double stepNext = memory[0].step();
    Teuchos::RCP<const NOX::Abstract::Vector> sPtrNext = 
      memory[0].sPtr();

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
      c = sPtr->innerProduct(dir) / memory[i].sNormSqr();

      dir.update(a * c, *sPtrNext, b * c, *sPtr, 1.0);
    }

    step = stepNext;
    sPtr = sPtrNext;

    a = sPtr->innerProduct(dir);		// <s,z>
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


void NOX::Direction::Broyden::throwError(const std::string& functionName, const std::string& errorMsg)
{
    if (utils->isPrintType(NOX::Utils::Error))
      utils->err() << "NOX::Direction::Broyden::" << functionName << " - " << errorMsg << std::endl;
    throw "NOX Error";
}
