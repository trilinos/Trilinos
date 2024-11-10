// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SCALED_MODEL_EVALUATOR_HPP
#define THYRA_SCALED_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"

namespace Thyra {


/** \brief This class decorates a ModelEvaluator and returns scaled
 * residual and Jacobian values.
 *
 * Given a scaling vector <tt>s</tt>, this object is treated as a diagonal
 * scaling matrix and applied to <tt>x -> Sf(x)</tt> and <tt>x -> sW</tt>.
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class ScaledModelEvaluator : 
    virtual public ModelEvaluatorDelegatorBase<Scalar>
{
public:
  
  /** \brief Constructs to uninitialized */
  ScaledModelEvaluator();
  
  /** \brief . */
  std::string description() const;

  /** \brief . */
  void set_f_scaling(const RCP<const Thyra::VectorBase<Scalar> >& f_scaling);

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}
  
private:

  //* Diagonal scaling vector */
  RCP<const Thyra::VectorBase<Scalar> > f_scaling_;

};


/** \brief Nonmember constructor. */
template<class Scalar>
RCP<ScaledModelEvaluator<Scalar> >
createNonconstScaledModelEvaluator(const RCP<ModelEvaluator<Scalar > > &model)
{
  RCP<ScaledModelEvaluator<Scalar> > srme(new ScaledModelEvaluator<Scalar>);
  srme->initialize(model);
  return srme;
}


} // namespace Thyra


#endif // THYRA_SCALED_MODEL_EVALUATOR_HPP
