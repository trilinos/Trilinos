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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef RYTHMOS_FORWARD_SENSITIVITY_MODEL_EVALUATOR_BASE_HPP
#define RYTHMOS_FORWARD_SENSITIVITY_MODEL_EVALUATOR_BASE_HPP


#include "Rythmos_IntegratorBase.hpp"
#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_PhysicallyBlockedLinearOpWithSolveBase.hpp" // Interface
#include "Thyra_DefaultBlockedTriangularLinearOpWithSolve.hpp" // Implementation
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_Assert.hpp"


namespace Rythmos {


/** \brief Forward sensitivity transient <tt>ModelEvaluator</tt> node
 * interface class.
 *
 * This base interface class provides the generalities of a linear forward
 * sensitivity model evaluator for a differential equation.
 *  
 * There are two derived classes which implement an implicit DAE and explicit
 * ODE formulation of the sensitivity equations.
 */
template<class Scalar>
class ForwardSensitivityModelEvaluatorBase
  : virtual public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  /** \brief Initialize the structure of the model.
   *
   * \param stateModel [in,persisting] The ModelEvaluator that defines the
   * parameterized state model <tt>f(x_dot,x,p)</tt>.
   *
   * \param p_index [in] The index of the parameter subvector in
   * <tt>stateModel</tt> for which sensitivities will be computed for.
   *
   * This function only intializes the spaces etc. needed to define structure
   * of the problem.  <tt>*this</tt> model object is not fully initialized at
   * this point in that <tt>evalModel()</tt> will not work yet and will thrown
   * exceptions if called.  The function <tt>initalizeState()</tt> must be
   * called later in order to fully initalize the model.
   */
  virtual void initializeStructure(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &stateModel,
    const int p_index
    ) = 0;

  /** \brief Initialize the structure of the model for an initial condition
   * only sensitivity problem.
   *
   * \param stateModel [in,persisting] The ModelEvaluator that defines the
   * parameterized state model <tt>f(x_dot,x,p)</tt>.
   *
   * \param p_space [in] The vector space for the parameterized initial
   * condition parameters.
   *
   * This function only intializes the spaces etc. needed to define structure
   * of the problem.  <tt>*this</tt> model object is not fully initialized at
   * this point in that <tt>evalModel()</tt> will not work yet and will thrown
   * exceptions if called.  The function <tt>initalizeState()</tt> must be
   * called later in order to fully initalize the model.
   */
  virtual void initializeStructureInitCondOnly(
    const RCP<const Thyra::ModelEvaluator<Scalar> >& stateModel,
    const RCP<const Thyra::VectorSpaceBase<Scalar> >& p_space
    ) = 0;
  
  /** \brief . */
  virtual RCP<const Thyra::ModelEvaluator<Scalar> >
  getStateModel() const =0;

  /** \brief . */
  virtual RCP<Thyra::ModelEvaluator<Scalar> >
  getNonconstStateModel() const =0;
  
  /** \brief . */
  virtual int get_p_index() const = 0;
  
  /** \brief . */
  virtual RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> >
  get_s_bar_space() const = 0;
  
  /** \brief . */
  virtual RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space() const = 0;

  /** \brief . */
  virtual void initializePointState(
      Ptr<StepperBase<Scalar> > stateStepper,
      bool forceUpToDateW
      ) =0;

};


/** \brief Create a wrapped non-const s_bar vector object given a non-const S
 * multi-vector object.
 */
template<class Scalar>
RCP<const Thyra::VectorBase<Scalar> > create_s_bar_given_S(
  const ForwardSensitivityModelEvaluatorBase<Scalar> &fwdSensModel,
  const RCP<Thyra::MultiVectorBase<Scalar> > &S
  )
{
  return Thyra::multiVectorProductVector(fwdSensModel.get_s_bar_space(), S); 
}


/** \brief Create a wrapped const s_bar vector object given a const S
 * multi-vector object.
 */
template<class Scalar>
RCP<const Thyra::VectorBase<Scalar> > create_s_bar_given_S(
  const ForwardSensitivityModelEvaluatorBase<Scalar> &fwdSensModel,
  const RCP<const Thyra::MultiVectorBase<Scalar> > &S
  )
{
  return Thyra::multiVectorProductVector(fwdSensModel.get_s_bar_space(), S); 
}


} // namespace Rythmos


#endif // RYTHMOS_FORWARD_SENSITIVITY_MODEL_EVALUATOR_BASE_HPP
