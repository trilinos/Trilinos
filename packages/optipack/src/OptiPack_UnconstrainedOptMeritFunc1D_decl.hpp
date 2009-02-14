/*
// @HEADER
// ***********************************************************************
// 
//    OptiPack: Collection of simple Thyra-based Optimization ANAs
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef OPTIPACK_UNCONSTRAINED_OPT_MERIT_FUNC_1D_DECL_HPP
#define OPTIPACK_UNCONSTRAINED_OPT_MERIT_FUNC_1D_DECL_HPP


#include "OptiPack_Types.hpp"
#include "GlobiPack_MeritFunc1DBase.hpp"
#include "OptiPack_LineSearchPointEvaluatorBase.hpp"
#include "Thyra_OperatorVectorTypes.hpp"


namespace OptiPack {


/** \brief Concreate subclass for unconstrained optimization objective
 * function.
 *
 * This subclass turns a response-only ModelEvaluator for an unconstrained
 * optimization problem <tt>g(p)</tt> into a 1D merit function.
 *
 * ToDo: Finish Documentation!
 */
template<typename Scalar>
class UnconstrainedOptMeritFunc1D
  : public GlobiPack::MeritFunc1DBase<typename ScalarTraits<Scalar>::magnitudeType>
{
public:

  /** \brief . */
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructor/Initializers/Accessors */
  //@{

  /** \brief Construct with default parameters. */
  UnconstrainedOptMeritFunc1D();

  /** \brief Set the model. */
  void setModel(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
    const int paramIndex,
    const int responseIndex
    );

  /** \brief Set the evaluation qualities.
   *
   * \param pointEvaluator [persisting non-modified] Evaluates p = p(alpha).
   *
   * \param p [persisting modified] Point that is updated with every call to
   * <tt>eval(...)</tt>.
   *
   * \param g_vec [persisting modified] The 1D vector used to store the
   * objective value computed by <tt>model</tt>.
   *
   * \param g_grad_vec [persisting modified] If <tt>!is_null(g_grad_vec)</tt>,
   * then this will be updated when <tt>eval(...)</tt> is called.
   *
   * \param baseDeriv [in] If <tt>!is_null(baseDeriv)</tt>, then gives the
   * value to be returned in <tt>this->baseDeriv()</tt>.
   *
   * <b>Postconditions:</b><ul>
   *
   * <li> [<tt>!is_null(g_grad_vec)</tt>]
   * <tt>this->supportsDerivEvals()==true</tt>
   *
   * <li> [<tt>!is_null(baseDeriv)</tt>]
   * <tt>this->supportsBaseDeriv()==true</tt>
   *
   * </ul>
   */
  void setEvaluationQuantities(
    const RCP<const LineSearchPointEvaluatorBase<Scalar> > &pointEvaluator,
    const RCP<Thyra::VectorBase<Scalar> > &p,
    const RCP<Thyra::VectorBase<Scalar> > &g_vec,
    const RCP<Thyra::VectorBase<Scalar> > &g_grad_vec
    );

  //@}

  /** \name Overridden from MeritFunc1DBase. */
  //@{
  
  /** \brief . */
  virtual bool supportsDerivEvals() const;
  /** \brief . */
  virtual void eval( const ScalarMag &alpha, const Ptr<ScalarMag> &phi,
    const Ptr<ScalarMag> &Dphi ) const;

  //@}

private:

  // //////////////////////
  // Private data members

  RCP<const Thyra::ModelEvaluator<Scalar> > model_;
  int paramIndex_;
  int responseIndex_;

  RCP<const LineSearchPointEvaluatorBase<Scalar> > pointEvaluator_;
  RCP<Thyra::VectorBase<Scalar> > p_;
  RCP<Thyra::VectorBase<Scalar> > g_vec_;
  RCP<Thyra::VectorBase<Scalar> > g_grad_vec_;

};


/** \brief Nonmember constructor.
 *
 * \relates UnconstrainedOptMeritFunc1D
 */
template<typename Scalar>
const RCP<UnconstrainedOptMeritFunc1D<Scalar> >
unconstrainedOptMeritFunc1D(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
  const int paramIndex,
  const int responseIndex
  )
{
  const RCP<UnconstrainedOptMeritFunc1D<Scalar> > meritFunc = 
    Teuchos::rcp(new UnconstrainedOptMeritFunc1D<Scalar>);
  meritFunc->setModel(model, paramIndex, responseIndex);
  return meritFunc;
}


} // namespace OptiPack


#endif // OPTIPACK_UNCONSTRAINED_OPT_MERIT_FUNC_1D_DECL_HPP
