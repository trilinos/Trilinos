// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

  virtual void     computeQuadraticMinimizer(const NOX::Abstract::Group &grp,
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
