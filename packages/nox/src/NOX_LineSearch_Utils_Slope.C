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

#include "NOX_LineSearch_Utils_Slope.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_GlobalData.H"

NOX::LineSearch::Utils::Slope::
Slope(const Teuchos::RCP<NOX::GlobalData>& gd) :
  utils(*(gd->getUtils()))
{

}

NOX::LineSearch::Utils::Slope::~Slope()
{

}

void NOX::LineSearch::Utils::Slope::
reset(const Teuchos::RCP<NOX::GlobalData>& gd)
{
  utils = *(gd->getUtils());
}

double NOX::LineSearch::Utils::Slope::
computeSlope(const Abstract::Vector& dir, const Abstract::Group& grp) 
{
   if (grp.isGradient()) 
     return(dir.innerProduct(grp.getGradient()));

  // Allocate space for vecPtr if necessary
   if (Teuchos::is_null(vecPtr)) 
     vecPtr = dir.clone(ShapeCopy);

  // v = J * dir
  NOX::Abstract::Group::ReturnType status = grp.applyJacobian(dir,*vecPtr);
  
  if (status != NOX::Abstract::Group::Ok) 
  {
    utils.out() << "NOX::LineSearch::Utils::Slope::computeSlope -  Unable to apply Jacobian!" << std::endl;
    throw "NOX Error";
  }

  // Check that F exists
  if (!grp.isF()) 
  {
    utils.out() << "NOX::LineSearch::Utils::Slope::computeSlope - Invalid F" << std::endl;
    throw "NOX Error";
  }

  // Return <v, F> = F' * J * dir = <J'F, dir> = <g, dir>
  return(vecPtr->innerProduct(grp.getF()));
}

double NOX::LineSearch::Utils::Slope::
computeSlopeWithOutJac(const Abstract::Vector& dir, 
		       const Abstract::Group& grp) 
{
  // Allocate space for vecPtr and grpPtr if necessary
  if (Teuchos::is_null(vecPtr)) 
    vecPtr = dir.clone(ShapeCopy);
  if (Teuchos::is_null(grpPtr))
    grpPtr = grp.clone(ShapeCopy);

  // Check that F exists
  if (!grp.isF()) 
  {
    utils.out() << "NOX::LineSearch::Utils::Slope::computeSlope - Invalid F" << std::endl;
    throw "NOX Error";
  }

  // Compute the perturbation parameter
  double lambda = 1.0e-6;
  double denominator = dir.norm();

  // Don't divide by zero
  if (denominator == 0.0)
    denominator = 1.0;

  double eta = lambda * (lambda + grp.getX().norm() / denominator);

  // Don't divide by zero
  if (eta == 0.0)
    eta = 1.0e-6;

  // Perturb the solution vector
  vecPtr->update(eta, dir, 1.0, grp.getX(), 0.0);

  // Compute the new F --> F(x + eta * dir)
  grpPtr->setX(*vecPtr);  
  grpPtr->computeF();

  // Compute Js = (F(x + eta * dir) - F(x))/eta
  vecPtr->update(-1.0/eta, grp.getF(), 1.0/eta, grpPtr->getF(), 0.0);
  
  return(vecPtr->innerProduct(grp.getF()));
}
