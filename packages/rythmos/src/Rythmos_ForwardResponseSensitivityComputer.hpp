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

#ifndef RYTHMOS_FORWARD_RESPONSE_SENSITIVITY_COMPUTER_HPP
#define RYTHMOS_FORWARD_RESPONSE_SENSITIVITY_COMPUTER_HPP


#include "Rythmos_Types.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"


namespace Rythmos {


/** \brief Concrete utility class for computing (assembling) forward transient
 * response sensitivities.
 *
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class ForwardResponseSensitivityComputer
  : public Teuchos::VerboseObject<ForwardResponseSensitivityComputer<Scalar> >
{
public:
  
  /** \brief . */
  ForwardResponseSensitivityComputer();

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, dumpSensitivities );

  /** \brief Set the response function for the first time.
   *
   * \param responseFunc [in,persisting] The response function that gives the
   * structure of response.
   *
   * \param basePoint [in] The base point for the calculation of the response
   * function.  Note that this must also include the current values of the
   * parameters!  This can be empty as long as it will be given later.
   *
   * \param p_index [in] The index of the parameter subvector in the response
   * function.
   *
   * \param g_index [in] The index of the response function(s).
   *
   * This sets the structure of the response function.
   */
  void setResponseFunction(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &responseFunc,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
    const int p_index,
    const int g_index
    );

  /** \brief Reset the point-specific response function along with its base
   * point.
   *
   */
  void resetResponseFunction(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &responseFunc,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint
    );

  /** \brief . */
  const RCP<Thyra::VectorBase<Scalar> > create_g_hat() const;

  /** \brief . */
  const RCP<Thyra::MultiVectorBase<Scalar> > create_D_g_hat_D_p() const;

  /** \brief Compute the reduced response at a point
   * (xdot,x,t).
   *
   * \param xdot [in,optional]
   *
   * \param x [in] The x vector.
   *
   * \param t [in] The time point
   *
   * \param g_hat [out,optional] The output response function, if set.  This
   * can be created by calling this->create_g_hat().
   */
  void computeResponse(
    const Thyra::VectorBase<Scalar> *x_dot,
    const Thyra::VectorBase<Scalar> &x,
    const Scalar t,
    Thyra::VectorBase<Scalar> *g_hat
    ) const;

  /** \brief Compute the reduced sensitivity and perhaps the response itself
   * at a point (xdot,x,t).
   *
   * \param x_dot [in,optional]
   *
   * \param S_dot [in,optional] 
   *
   * \param x [in]
   *
   * \param S [in]
   *
   * \param t [in]
   *
   * \param g_hat [out,optional] The output response function, if set.  This
   * can be created by calling this->create_g_hat().
   *
   * \param D_g_hat_D_p [out] The output response function reduced parameter
   * derivative.  This can be created by calling this->create_D_g_hat_D_p().
   */
  void computeResponseAndSensitivity(
    const Thyra::VectorBase<Scalar> *x_dot,
    const Thyra::MultiVectorBase<Scalar> *S_dot,
    const Thyra::VectorBase<Scalar> &x,
    const Thyra::MultiVectorBase<Scalar> &S,
    const Scalar t,
    Thyra::VectorBase<Scalar> *g_hat,
    Thyra::MultiVectorBase<Scalar> *D_g_hat_D_p
    ) const;

private: // Data members

  RCP<const Thyra::ModelEvaluator<Scalar> > responseFunc_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;
  int p_index_;
  int g_index_;

  RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  RCP<const Thyra::VectorSpaceBase<Scalar> > g_space_;

  bool response_func_supports_x_dot_;
  bool response_func_supports_D_x_dot_;
  bool response_func_supports_D_p_;

  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> responseInArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> responseOutArgs_;

  mutable RCP<Thyra::LinearOpBase<Scalar> > D_g_D_x_dot_;
  mutable RCP<Thyra::LinearOpBase<Scalar> > D_g_D_x_;
  mutable RCP<Thyra::MultiVectorBase<Scalar> > D_g_D_p_;

private: // Functions

  void clearCache();

  void createCache(const bool computeSens) const;

  void computeResponseAndSensitivityImpl(
    const Thyra::VectorBase<Scalar> *x_dot,
    const Thyra::MultiVectorBase<Scalar> *S_dot,
    const Thyra::VectorBase<Scalar> &x,
    const Thyra::MultiVectorBase<Scalar> *S,
    const Scalar t,
    Thyra::VectorBase<Scalar> *g_hat,
    Thyra::MultiVectorBase<Scalar> *D_g_hat_D_p
    ) const;

};


//
// Implementations
//


template<class Scalar>
ForwardResponseSensitivityComputer<Scalar>::ForwardResponseSensitivityComputer()
  :dumpSensitivities_(false),
   p_index_(-1),
   g_index_(-1),
   response_func_supports_x_dot_(false),
   response_func_supports_D_x_dot_(false),
   response_func_supports_D_p_(false)
{}


template<class Scalar>
void ForwardResponseSensitivityComputer<Scalar>::setResponseFunction(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &responseFunc,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
  const int p_index,
  const int g_index
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

  // ToDo: Validate input!

  responseFunc_ = responseFunc;
  basePoint_ = basePoint;
  p_index_ = p_index;
  g_index_ = g_index;

  p_space_ = responseFunc_->get_p_space(p_index_);
  g_space_ = responseFunc_->get_g_space(g_index_);


  MEB::InArgs<Scalar>
    responseInArgs = responseFunc_->createInArgs();
  response_func_supports_x_dot_ =
    responseInArgs.supports(MEB::IN_ARG_x_dot);
  MEB::OutArgs<Scalar>
    responseOutArgs = responseFunc_->createOutArgs();
  response_func_supports_D_x_dot_ =
    !responseOutArgs.supports(MEB::OUT_ARG_DgDx_dot,g_index_).none();
  response_func_supports_D_p_ =
    !responseOutArgs.supports(MEB::OUT_ARG_DgDp,g_index_,p_index_).none();

  clearCache();

}


template<class Scalar>
void ForwardResponseSensitivityComputer<Scalar>::resetResponseFunction(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &responseFunc,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint
  )
{
  // ToDo: Validate that responseFunc is the same structure as the one already
  // set!
  responseFunc_ = responseFunc;
  basePoint_ = basePoint;
}


template<class Scalar>
const RCP<Thyra::VectorBase<Scalar> >
ForwardResponseSensitivityComputer<Scalar>::create_g_hat() const
{
  return Thyra::createMember(g_space_);
}


template<class Scalar>
const RCP<Thyra::MultiVectorBase<Scalar> >
ForwardResponseSensitivityComputer<Scalar>::create_D_g_hat_D_p() const
{
  return Thyra::createMembers(g_space_,p_space_->dim());
}


template<class Scalar>
void ForwardResponseSensitivityComputer<Scalar>::computeResponse(
  const Thyra::VectorBase<Scalar> *x_dot,
  const Thyra::VectorBase<Scalar> &x,
  const Scalar t,
  Thyra::VectorBase<Scalar> *g_hat
  ) const
{
  computeResponseAndSensitivityImpl(x_dot,0,x,0,t,g_hat,0);
}


template<class Scalar>
void ForwardResponseSensitivityComputer<Scalar>::computeResponseAndSensitivity(
  const Thyra::VectorBase<Scalar> *x_dot,
  const Thyra::MultiVectorBase<Scalar> *S_dot,
  const Thyra::VectorBase<Scalar> &x,
  const Thyra::MultiVectorBase<Scalar> &S,
  const Scalar t,
  Thyra::VectorBase<Scalar> *g_hat,
  Thyra::MultiVectorBase<Scalar> *D_g_hat_D_p
  ) const
{
  computeResponseAndSensitivityImpl(x_dot,S_dot,x,&S,t,g_hat,D_g_hat_D_p);
}


// private


template<class Scalar>
void ForwardResponseSensitivityComputer<Scalar>::clearCache()
{
  D_g_D_x_dot_ = Teuchos::null;
  D_g_D_x_ = Teuchos::null;
  D_g_D_p_ = Teuchos::null;
}


template<class Scalar>
void ForwardResponseSensitivityComputer<Scalar>::createCache(
  const bool computeSens
  ) const
{
  if (computeSens) {
    if (response_func_supports_D_x_dot_ && is_null(D_g_D_x_dot_))
      D_g_D_x_dot_ = responseFunc_->create_DgDx_dot_op(g_index_);
    D_g_D_x_ = responseFunc_->create_DgDx_op(g_index_);
    if (response_func_supports_D_p_ && is_null(D_g_D_p_))
      D_g_D_p_ = createMembers(g_space_,p_space_->dim());
  }
}


template<class Scalar>
void ForwardResponseSensitivityComputer<Scalar>::computeResponseAndSensitivityImpl(
  const Thyra::VectorBase<Scalar> *x_dot,
  const Thyra::MultiVectorBase<Scalar> *S_dot,
  const Thyra::VectorBase<Scalar> &x,
  const Thyra::MultiVectorBase<Scalar> *S,
  const Scalar t,
  Thyra::VectorBase<Scalar> *g_hat,
  Thyra::MultiVectorBase<Scalar> *D_g_hat_D_p
  ) const
{

  using Teuchos::rcp;
  using Teuchos::ptr;
  using Teuchos::includesVerbLevel;
  typedef ScalarTraits<Scalar> ST;
  using Thyra::apply;
  using Thyra::Vp_V;
  typedef Thyra::ModelEvaluatorBase MEB;

  //
  // A) Setup for output
  //

  const RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  const bool trace =
    out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW);
  const bool print_norms =
    out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM);
  const bool print_x =
    out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_EXTREME);

  //
  // B) Initialize storage
  //

  const bool computeSens = ( D_g_hat_D_p != 0 );
  createCache(computeSens);

  //
  // C) Evaluate the response function
  //

  //
  // C.1) Setup input/output and evaluate the response function
  //

  if (trace)
    *out << "\nEvaluating response function at time t = " << t << " ...\n";
   
  // C.1.a) Setup the input arguments
   
  MEB::InArgs<Scalar> responseInArgs = responseFunc_->createInArgs();
  responseInArgs.setArgs(basePoint_);
  responseInArgs.set_x(rcp(&x,false));
  if (response_func_supports_x_dot_)
    responseInArgs.set_x_dot(rcp(x_dot,false));
      
  // C.1.b) Setup output arguments
  
  MEB::OutArgs<Scalar> responseOutArgs = responseFunc_->createOutArgs();
  
  if (g_hat)
    responseOutArgs.set_g(g_index_,rcp(g_hat,false));
  
  if (computeSens) {
        
    // D_g_D_x_dot
    if (response_func_supports_D_x_dot_) {
      responseOutArgs.set_DgDx_dot(
        g_index_,
        MEB::Derivative<Scalar>(D_g_D_x_dot_)
        );
    }
    
    // D_g_D_x
    responseOutArgs.set_DgDx(
      g_index_,
      MEB::Derivative<Scalar>(D_g_D_x_)
      );
    
    // D_g_D_p
    if (response_func_supports_D_p_) {
      responseOutArgs.set_DgDp(
        g_index_, p_index_,
        MEB::Derivative<Scalar>(D_g_D_p_,MEB::DERIV_MV_BY_COL)
        );
    }
    
  }
      
  // C.1.c) Evaluate the response function k

  {
#ifdef ENABLE_RYTHMOS_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("Rythmos:ForwardResponseSensitivityComputer::evalModel: evalResponse");
#endif
    responseFunc_->evalModel( responseInArgs, responseOutArgs );
  }
  
  // C.1.d) Print the outputs just coputed
  
  if (print_norms) {
    if (g_hat)
      *out << "\n||g_hat||inf = " << norm_inf(*g_hat) << endl;
    if (computeSens && !is_null(D_g_D_p_))
      *out << "\n||D_g_D_p_||inf = " << norms_inf(*D_g_D_p_) << endl;
  }
  
  if ( g_hat && (dumpSensitivities_ || print_x) )
    *out << "\ng_hat = " << *g_hat;
      
  if (computeSens && print_x) {
    if (!is_null(D_g_D_x_dot_))
      *out << "\nD_g_D_x_dot = " << *D_g_D_x_ << endl;
    if (!is_null(D_g_D_x_))
      *out << "\nD_g_D_x = " << *D_g_D_x_ << endl;
    if (!is_null(D_g_D_p_))
      *out << "\nD_g_D_p = " << *D_g_D_p_ << endl;
  }
  
  //
  // C.2) Assemble the output response function sensitivity D_d_hat_D_p
  //

  // D_g_hat_D_p = DgDx_dot * S_dot + DgDx * S + DgDp 
  
  if (computeSens) {

#ifdef ENABLE_RYTHMOS_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("Rythmos:ForwardResponseSensitivityComputer::evalModel: computeSens");
#endif
    
    if (trace)
      *out << "\nD_g_hat_D_p = DgDx_dot * S_dot + DgDx * S + DgDp ...\n";
    
    assign( ptr(D_g_hat_D_p), ST::zero() );
      
    // D_g_hat_D_p += DgDx_dot * S_dot
    if (response_func_supports_D_x_dot_) {
      apply( *D_g_D_x_dot_, Thyra::NOTRANS, *S_dot,
        ptr(D_g_hat_D_p), ST::one(), ST::one() );
    }
    
    // D_g_hat_D_p += DgDx * S
    apply( *D_g_D_x_, Thyra::NOTRANS, *S,
      ptr(D_g_hat_D_p), ST::one(), ST::one() );
    
    // D_g_hat_D_p += DgDp
    if (response_func_supports_D_p_) {
      Vp_V( ptr(D_g_hat_D_p), *D_g_D_p_ );
    }
      
    if (dumpSensitivities_ || print_x)
      *out << "\nD_g_hat_D_p = "
           << Teuchos::describe(*D_g_hat_D_p,Teuchos::VERB_EXTREME);
    
  }

}


} // namespace Rythmos


#endif // RYTHMOS_FORWARD_RESPONSE_SENSITIVITY_COMPUTER_HPP
