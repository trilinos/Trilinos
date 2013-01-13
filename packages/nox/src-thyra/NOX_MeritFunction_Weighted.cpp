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

#include "NOX_MeritFunction_Weighted.hpp"
#include "NOX_Thyra_Vector.H"
#include "NOX_Thyra_Group.H"
#include "Teuchos_Assert.hpp"

NOX::Thyra::WeightedMeritFunction::
WeightedMeritFunction(const Teuchos::RCP<const ::Thyra::VectorBase<double> > weights, 
		      bool optimizeSlope) :
  name_("WeightedMeritFunction"),
  optimizeSlopeCalc_(optimizeSlope)
{
  using Teuchos::rcp;

  // const_cast is nasty.  Need to add support for const constructor
  // in NOX::Thyra::Vector class.
  weights_ = rcp(new NOX::Thyra::Vector(Teuchos::rcp_const_cast< ::Thyra::VectorBase<double> >(weights))); // no allocation
  tmp1_ = rcp(new NOX::Thyra::Vector(*weights_)); // performs allocation
  tmp2_ = rcp(new NOX::Thyra::Vector(*weights_)); // performs allocation
}

NOX::Thyra::WeightedMeritFunction::
WeightedMeritFunction(const WeightedMeritFunction& source) :
  weights_(source.weights_),
  name_(source.name_),
  optimizeSlopeCalc_(source.optimizeSlopeCalc_)
{
  tmp1_ = Teuchos::rcp(new NOX::Thyra::Vector(*weights_));
  tmp2_ = Teuchos::rcp(new NOX::Thyra::Vector(*weights_));
}

NOX::Thyra::WeightedMeritFunction::~WeightedMeritFunction()
{ }

const std::string& NOX::Thyra::WeightedMeritFunction::name() const
{
  return name_;
}

std::ostream& NOX::Thyra::WeightedMeritFunction::print(std::ostream& os, int indent) const
{
  return os;
}

double NOX::Thyra::WeightedMeritFunction::
computef(const NOX::Abstract::Group& group) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  TEUCHOS_ASSERT(group.isF());
  
  RCP<const NOX::Thyra::Vector> f =
    Teuchos::rcp_dynamic_cast<const NOX::Thyra::Vector>(group.getFPtr(),true);

  ::Thyra::copy(f->getThyraVector(), outArg(tmp1_->getThyraVector()));
  ::Thyra::ele_wise_scale(weights_->getThyraVector(), outArg(tmp1_->getThyraVector()));

  double tmp = tmp1_->norm();

  double value = 0.5 * tmp * tmp;

  return value;
}

double NOX::Thyra::WeightedMeritFunction::
computeSlope(const NOX::Abstract::Vector& dir,
	     const NOX::Abstract::Group& group) const
{
  double value = 0.0;

  //const NOX::Thyra::Vector& epDir = dynamic_cast<const NOX::Thyra::Vector&>(dir);

  if (optimizeSlopeCalc_) {  // use matrix-free to evaluate F^T Js

    if ( is_null(tmpGrpPtr) ) {
      tmpGrpPtr = Teuchos::rcp_dynamic_cast<NOX::Thyra::Group>
	(group.clone(NOX::ShapeCopy),true);
    }
    
    // these norms are implicitly weighted
    double numerator = group.getX().norm();
    double denominator = dir.norm();
 
    if (denominator < 1.0e-12) {
      std::cout << "WARNING: NOX::Thyra::WeightedMeritFunction::computeSlope()\n"
		<< "denominator in perturbation < 1.0e-12! Setting to 1.0!" << std::endl;
      denominator = 1.0;
    }

    double eta = 1.0e-6 * (1.0e-6 + numerator/denominator);

    *tmp1_ = dir;

    tmpGrpPtr->computeX(group, dir, eta);
    
    tmpGrpPtr->computeF();

    *tmp1_ = tmpGrpPtr->getF();

    tmp1_->update(-1.0, group.getF(), 1.0);
    
    tmp1_->scale(1.0/eta);

//     scalingPtr->applyLeftScaling(tmp1_->getEpetraVector(), 
// 			   tmp1_->getEpetraVector());

    ::Thyra::ele_wise_scale(weights_->getThyraVector(), outArg(tmp1_->getThyraVector()));

    *tmp2_ = group.getF();

//     scalingPtr->applyLeftScaling(tmp2_->getEpetraVector(), 
// 			   tmp2_->getEpetraVector());

    ::Thyra::ele_wise_scale(weights_->getThyraVector(), outArg(tmp2_->getThyraVector()));

    value = tmp2_->innerProduct(*tmp1_);

  }

  else {   // explicitly compute F^T Js

    group.applyJacobian(dir, *tmp1_);
    
//     scalingPtr->applyLeftScaling(tmp1_->getEpetraVector(), 
// 			   tmp1_->getEpetraVector());

    ::Thyra::ele_wise_scale(weights_->getThyraVector(), outArg(tmp1_->getThyraVector()));

    *tmp2_ = group.getF();

//     scalingPtr->applyLeftScaling(tmp2_->getEpetraVector(), 
// 			   tmp2_->getEpetraVector());

    ::Thyra::ele_wise_scale(weights_->getThyraVector(), outArg(tmp2_->getThyraVector()));

    value = tmp2_->innerProduct(*tmp1_);

  }

  return value;
}

void NOX::Thyra::WeightedMeritFunction::
computeGradient(const NOX::Abstract::Group& group,
		NOX::Abstract::Vector& result) const
{
  TEUCHOS_ASSERT(group.isF());
  TEUCHOS_ASSERT(group.isJacobian());

  const NOX::Thyra::Vector& f = 
    dynamic_cast<const NOX::Thyra::Vector&>(group.getF());
  
  NOX::Thyra::Vector& tResult = dynamic_cast<NOX::Thyra::Vector&>(result);

  // Apply D squared
//   scalingPtr->applyLeftScaling(fVec.getEpetraVector(), 
// 			   tmp1_->getEpetraVector());
//   scalingPtr->applyLeftScaling(tmp1_->getEpetraVector(), 
// 			   tmp1_->getEpetraVector());

  ::Thyra::copy(f.getThyraVector(), outArg(tmp1_->getThyraVector()));
  ::Thyra::ele_wise_scale(weights_->getThyraVector(), outArg(tmp1_->getThyraVector()));
  ::Thyra::ele_wise_scale(weights_->getThyraVector(), outArg(tmp1_->getThyraVector()));

  NOX::Abstract::Group::ReturnType status = 
    group.applyJacobianTranspose(*tmp1_, tResult);

  TEUCHOS_TEST_FOR_EXCEPTION(status != NOX::Abstract::Group::Ok, std::logic_error,
			     "applyJacobianTranspose() failed!");
}

double NOX::Thyra::WeightedMeritFunction::
computeQuadraticModel(const NOX::Abstract::Vector& dir,
		      const NOX::Abstract::Group& group) const
{
  double f = this->computef(group);

  this->computeGradient(group, *tmp2_);

  double gradfTd = tmp2_->innerProduct(dir);

  group.applyJacobian(dir, *tmp2_);

//   scalingPtr->applyLeftScaling(tmp2_->getEpetraVector(), 
// 			   tmp2_->getEpetraVector());

  ::Thyra::ele_wise_scale(weights_->getThyraVector(), outArg(tmp2_->getThyraVector()));
  
  double dTgrad2fd = 0.5 * tmp2_->innerProduct(*tmp2_);

  double m = f + gradfTd + dTgrad2fd; 

  return m;
}

void NOX::Thyra::WeightedMeritFunction::
computeQuadraticMinimizer(const NOX::Abstract::Group &grp, 
			  NOX::Abstract::Vector &result) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
			     "NOX::Thyra::WeightedMeritFunction::computeQuadraticMinimizer() method has not been implemented yet!");
}

bool NOX::Thyra::WeightedMeritFunction::
computeSteepestDescentDir(const NOX::Abstract::Group& group,
			  NOX::Abstract::Vector& result) const
{
  TEUCHOS_ASSERT(group.isF());
  TEUCHOS_ASSERT(group.isJacobian());

  NOX::Thyra::Vector& dir = dynamic_cast<NOX::Thyra::Vector&>(result);

  NOX::Abstract::Group::ReturnType status;

  // Compute Steepest Descent direction
  this->computeGradient(group, dir);

  // Scale by the quadratic minimization model
  status = group.applyJacobian(dir, *tmp2_);
  TEUCHOS_TEST_FOR_EXCEPTION(status != NOX::Abstract::Group::Ok, std::logic_error,
			     "applyJacobian() failed!");

//   scalingPtr->applyLeftScaling(tmp2_->getEpetraVector(), 
//   		   tmp2_->getEpetraVector());

  ::Thyra::ele_wise_scale(weights_->getThyraVector(), outArg(tmp2_->getThyraVector()));
  
  dir.scale( -1.0 * dir.innerProduct(dir) / tmp2_->innerProduct(*tmp2_) );
  
  return true;
}
