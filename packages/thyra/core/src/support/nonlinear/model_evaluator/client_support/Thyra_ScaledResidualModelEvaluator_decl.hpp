// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_SCALED_RESIDUAL_MODEL_EVALUATOR_HPP
#define THYRA_SCALED_RESIDUAL_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"

namespace Thyra {


/** \brief This class decorates a ModelEvaluator and returns scaled
 * residual and Jacobian values.  Given a scaling vector <tt>s</tt>, this
 * object is treated as a diagonal scaling matrix and applied to <tt>x
 * -> Sf(x)</tt> and <tt>x -> sW</tt>.
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class ScaledResidualModelEvaluator : 
    virtual public ModelEvaluatorDelegatorBase<Scalar> {
public:
  
  //* Constructs to uninitialized */
  ScaledResidualModelEvaluator();

  //* Calls initialize() from ModelEvaluatorDelegatorBase. */
  ScaledResidualModelEvaluator(const RCP< ModelEvaluator< Scalar > > &model);

  //* Calls initialize() from ModelEvaluatorDelegatorBase. */
  ScaledResidualModelEvaluator(const RCP< const ModelEvaluator< Scalar > > &model);
  
  std::string description() const;

  void set_f_scaling(const RCP<Thyra::VectorBase<Scalar> >& f_scaling);

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
  RCP<Thyra::VectorBase<Scalar> > f_scaling_;

};

} // namespace Thyra


#endif // THYRA_SCALED_RESIDUAL_MODEL_EVALUATOR_HPP
