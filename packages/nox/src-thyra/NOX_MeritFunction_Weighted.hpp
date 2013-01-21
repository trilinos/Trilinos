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

#ifndef NOX_MERIT_FUNCTION_WEIGHTED_HPP
#define NOX_MERIT_FUNCTION_WEIGHTED_HPP

#include "NOX_MeritFunction_Generic.H" // base class
#include "Teuchos_RCP.hpp"
#include <string>

namespace Thyra {
  template<typename Scalar> class VectorBase;
}

namespace NOX {
namespace Thyra {

class Vector;
class Group;

/** \brief Implementation of merit function for implicitly weighted norm.
    
    NOTE: The nox vectors in this object are always unweighted (we
    always apply the weights explicitly)!  Be careful about using
    norms and innerProducts from incoming nox objects as these will
    have the implicit weighting attached.
*/
class WeightedMeritFunction : public NOX::MeritFunction::Generic {
  
public:			

  //! Constructor.
  explicit WeightedMeritFunction(const Teuchos::RCP<const ::Thyra::VectorBase<double> > weights, 
				 bool optimizeSlopeCalc = true);

  //! Copy constructor.
  explicit WeightedMeritFunction(const WeightedMeritFunction& source);

  //! Destructor.
  ~WeightedMeritFunction();

  virtual const std::string& name() const;
  
  virtual std::ostream& print(std::ostream& os, int indent=0) const;
  
  virtual double computef(const NOX::Abstract::Group& group) const;

  virtual void computeGradient(const NOX::Abstract::Group& group,
			       NOX::Abstract::Vector& result) const;
  
  virtual double computeSlope(const NOX::Abstract::Vector& dir,
			      const NOX::Abstract::Group& group) const;
  
  virtual double computeQuadraticModel(const NOX::Abstract::Vector& dir,
				       const NOX::Abstract::Group& group) const;
  
  virtual void 	computeQuadraticMinimizer(const NOX::Abstract::Group &grp, 
					  NOX::Abstract::Vector &result) const;
  
  virtual bool computeSteepestDescentDir(const NOX::Abstract::Group& group,
					 NOX::Abstract::Vector& result) const;

private:
  
  Teuchos::RCP<const NOX::Thyra::Vector> weights_;
  const std::string name_;
  bool optimizeSlopeCalc_;

  mutable Teuchos::RCP<NOX::Thyra::Vector> tmp1_;
  mutable Teuchos::RCP<NOX::Thyra::Vector> tmp2_;
  mutable Teuchos::RCP<NOX::Thyra::Group> tmpGrpPtr;
  
};
} // namespace Thyra
} // namespace NOX

#endif
