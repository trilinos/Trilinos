// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Common.H"

#ifdef WITH_PRERELEASE

#include "NOX_Direction_QuasiNewton.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"

//------------------------------------------------------------

NOX::Direction::QuasiNewton::MemoryUnit::MemoryUnit()
{
  sdotyValue = 0;
  ydotyValue = 0;
  rhoValue = 0;
}

NOX::Direction::QuasiNewton::MemoryUnit::~MemoryUnit()
{

}

void NOX::Direction::QuasiNewton::MemoryUnit::reset(const Abstract::Vector& newX,
                             const Abstract::Vector& oldX,
                             const Abstract::Vector& newG,
                             const Abstract::Vector& oldG)
{
  if (Teuchos::is_null(sPtr))
  {
    sPtr = newX.clone(ShapeCopy);
    yPtr = newX.clone(ShapeCopy);
  }

  sPtr->update(1.0, newX, -1.0, oldX, 0.0);
  yPtr->update(1.0, newG, -1.0, oldG, 0.0);
  sdotyValue = sPtr->innerProduct(*yPtr);
  ydotyValue = yPtr->innerProduct(*yPtr);
  rhoValue = 1.0 / sdotyValue;
}

const NOX::Abstract::Vector& NOX::Direction::QuasiNewton::MemoryUnit::s() const
{
  return *sPtr;
}

const NOX::Abstract::Vector& NOX::Direction::QuasiNewton::MemoryUnit::y() const
{
  return *yPtr;
}

double NOX::Direction::QuasiNewton::MemoryUnit::sdoty() const
{
  return sdotyValue;
}

double NOX::Direction::QuasiNewton::MemoryUnit::ydoty() const
{
  return ydotyValue;
}

double NOX::Direction::QuasiNewton::MemoryUnit::rho() const
{
  return rhoValue;
}

//------------------------------------------------------------

NOX::Direction::QuasiNewton::Memory::Memory(int m)
{
  reset(m);
}

void NOX::Direction::QuasiNewton::Memory::reset(int m)
{
  memory.resize(m);
  index.resize(0);
  index.reserve(m);
}

NOX::Direction::QuasiNewton::Memory::~Memory()
{
}

void NOX::Direction::QuasiNewton::Memory::add(const NOX::Abstract::Vector& newX,
                           const NOX::Abstract::Vector& oldX,
                           const NOX::Abstract::Vector& newG,
                           const NOX::Abstract::Vector& oldG)
{
  int m = index.size();

  if (m < (int) memory.size())
  {
    index.push_back(m);
  }
  else
  {
    // rotate(index, index + (m-1), index + m);
    int k = index[0];
    for (int i = 0; i < m - 1; i ++)
      index[i] = index[i+1];
    index[m-1] = k;
  }

  memory[index.back()].reset(newX,oldX,newG,oldG);
}

bool NOX::Direction::QuasiNewton::Memory::empty() const
{
  return index.empty();
}

int NOX::Direction::QuasiNewton::Memory::size() const
{
  return index.size();
}

const NOX::Direction::QuasiNewton::MemoryUnit&
NOX::Direction::QuasiNewton::Memory::operator[](int i) const
{
  return memory[index[i]];
}

//------------------------------------------------------------

NOX::Direction::QuasiNewton::
QuasiNewton(const Teuchos::RCP<NOX::GlobalData>& gd,
        Teuchos::ParameterList& p) :
  paramsPtr(NULL)
{
  reset(gd, p);
}

NOX::Direction::QuasiNewton::~QuasiNewton()
{
}

bool NOX::Direction::QuasiNewton::
reset(const Teuchos::RCP<NOX::GlobalData>& gd,
      Teuchos::ParameterList& params)
{
  globalDataPtr = gd;
  utils = gd->getUtils();
  paramsPtr = &params;
  Teuchos::ParameterList& p = params.sublist("Quasi-Newton");
  memory.reset(p.get("Memory", 5));
  return true;
}

bool NOX::Direction::QuasiNewton::compute(NOX::Abstract::Vector& dir,
                      NOX::Abstract::Group& soln,
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

  // Push the old information onto the memory, but only after at least one previous iteration
  if (solver.getNumIterations() > 0)
  {
    const NOX::Abstract::Group& oldSoln = solver.getPreviousSolutionGroup();
    if (oldSoln.isGradient())
      memory.add(soln.getX(), oldSoln.getX(), soln.getGradient(), oldSoln.getGradient());
  }

  // *** Calculate the QN direction ***

  // d = -g
  dir = soln.getGradient();
  dir.scale(-1.0);

  if (!memory.empty())
  {

    int m = memory.size();
    std::vector<double> alpha(m);
    double beta;

    for (int i = m-1; i >= 0; i --)
    {
      alpha[i] = memory[i].rho() * dir.innerProduct( memory[i].s() );
      dir.update(-1.0 * alpha[i], memory[i].y(), 1.0);
    }

    dir.scale( memory[m-1].sdoty() / memory[m-1].ydoty() );

    for (int i = 0; i < m; i ++)
    {
      beta = memory[i].rho() * dir.innerProduct( memory[i].y() );
      dir.update(alpha[i] - beta, memory[i].s(), 1.0);
    }
  }

  return true;
}

bool NOX::Direction::QuasiNewton::
compute(NOX::Abstract::Vector& dir,
    NOX::Abstract::Group& soln,
    const Solver::LineSearchBased& solver)
{
  return NOX::Direction::Generic::compute( dir, soln, solver );
}

void NOX::Direction::QuasiNewton::throwError(const std::string& functionName,
                         const std::string& errorMsg)
{
  if (utils->isPrintType(Utils::Error))
    utils->err() << "NOX::Direction::QuasiNewton::" << functionName
         << " - " << errorMsg << std::endl;
  throw std::runtime_error("NOX Error");
}


#endif

