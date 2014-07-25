//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef RYTHMOS_SINGLE_RESIDUAL_MODEL_EVALUATOR_BASE_HPP
#define RYTHMOS_SINGLE_RESIDUAL_MODEL_EVALUATOR_BASE_HPP


#include "Thyra_ModelEvaluator.hpp"


namespace Rythmos {


/** \brief Base class mix-in interface for single-residual model evaluators.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class SingleResidualModelEvaluatorBase
  : virtual public Thyra::ModelEvaluator<Scalar>
{
public:

  /** \name Intialization */
  //@{
  
  /** \brief . */
  virtual void initializeSingleResidualModel(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &daeModel,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
    const Scalar &coeff_x_dot,
    const RCP<const Thyra::VectorBase<Scalar> > &x_dot_base,
    const Scalar &coeff_x,
    const RCP<const Thyra::VectorBase<Scalar> > &x_base,
    const Scalar &t_base,
    const RCP<const Thyra::VectorBase<Scalar> > &x_bar_init
    ) = 0;

  //@}

  /** \name Query */
  //@{

  /** \brief . */
  virtual Scalar get_coeff_x_dot() const = 0;

  /** \brief . */
  virtual RCP<const Thyra::VectorBase<Scalar> >
  get_x_dot_base() const = 0;

  /** \brief . */
  virtual Scalar get_coeff_x() const = 0;

  /** \brief . */
  virtual RCP<const Thyra::VectorBase<Scalar> >
  get_x_base() const = 0;

  /** \brief . */
  virtual Scalar get_t_base() const = 0;

  //@}

};


} // namespace Rythmos


#endif // RYTHMOS_SINGLE_RESIDUAL_MODEL_EVALUATOR_BASE_HPP
