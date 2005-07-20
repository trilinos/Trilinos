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

#ifndef THYRA_MODEL_EVALUATOR_HPP
#define THYRA_MODEL_EVALUATOR_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"

namespace Thyra {

/** \brief . */
class ModelEvaluatorBase {
public:

  /** \name Public types */
  //@{

  /** \brief.  */
  enum EInArgsMembers {
    IN_ARG_x_dot
    ,IN_ARG_x
    ,IN_ARG_t
    ,IN_ARG_alpha
    ,IN_ARG_beta
  };
  static const int NUM_E_IN_ARGS_MEMBERS=5;

  /** \brief . */
  template<class Scalar>
  class InArgs {
  public:
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    InArgs();
    void set_x_dot( const Teuchos::RefCountPtr<const VectorBase<Scalar> > &x_dot );
    Teuchos::RefCountPtr<const VectorBase<Scalar> > get_x_dot() const;
    void set_x( const Teuchos::RefCountPtr<const VectorBase<Scalar> > &x );
    Teuchos::RefCountPtr<const VectorBase<Scalar> > get_x() const;
    void set_t( ScalarMag t );
    ScalarMag get_t() const;
    void set_alpha( Scalar alpha );
    Scalar get_alpha() const;
    void set_beta( Scalar beta );
    Scalar get_beta() const;
    bool supports(EInArgsMembers arg) const;
  protected:
    /** \brief . */
    void _setSupports( EInArgsMembers arg, bool supports );
  private:
    Teuchos::RefCountPtr<const VectorBase<Scalar> >  x_dot_;
    Teuchos::RefCountPtr<const VectorBase<Scalar> >  x_;
    ScalarMag                                        t_;
    Scalar                                           alpha_;
    Scalar                                           beta_;
    bool supports_[NUM_E_IN_ARGS_MEMBERS];
    void assert_supports(EInArgsMembers arg) const;
  };

  /** \brief.  */
  enum EOutArgsMembers {
    OUT_ARG_f
    ,OUT_ARG_W
  };
  static const int NUM_E_OUT_ARGS_MEMBERS=2;

  /** \brief . */
  template<class Scalar>
  class OutArgs {
  public:
    OutArgs();
    void set_f( const Teuchos::RefCountPtr<VectorBase<Scalar> > &f );
    Teuchos::RefCountPtr<VectorBase<Scalar> > get_f() const;
    void set_W( const Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > &W );
    Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > get_W() const;
    bool supports(EOutArgsMembers arg) const;
  protected:
    /** \brief . */
    void _setSupports( EOutArgsMembers arg, bool supports );
  private:
    Teuchos::RefCountPtr<VectorBase<Scalar> >             f_;
    Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >  W_;
    bool supports_[NUM_E_OUT_ARGS_MEMBERS];
    void assert_supports(EOutArgsMembers arg) const;
  };

  //@}

#ifdef HAVE_PROTECTED_NESTED_TEMPLATE_CLASS_ACCESS // Added since at least gcc 3.3.4 does not do the right thing here!
protected:
#endif

  /** \name Protected types */
  //@{


  /** \brief . */
  template<class Scalar>
  class InArgsSetup : public InArgs<Scalar> {
  public:
    /** \brief . */
    void setSupports( EInArgsMembers arg, bool supports = true ) { this->_setSupports(arg,supports); }
  };

  /** \brief . */
  template<class Scalar>
  class OutArgsSetup : public OutArgs<Scalar> {
  public:
    /** \brief . */
    void setSupports( EOutArgsMembers arg, bool supports = true ) { this->_setSupports(arg,supports); }
  };

  //@}

};

/** \brief Base interface for evaluating a stateless "model".
 *
 *
 * ToDo: Finish Documentation!
 */
template<class Scalar>
class ModelEvaluator : public ModelEvaluatorBase {
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \brief . */
  virtual ~ModelEvaluator() {}

  /** \name Pure virtual functions that must be overridden by subclasses. */
  //@{

  /** \breif . */
  virtual Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_x_space() const = 0;

  /** \breif . */
  virtual Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_f_space() const = 0;

  /** \brief . */
  virtual InArgs<Scalar> createInArgs() const = 0;

  /** \brief . */
  virtual OutArgs<Scalar> createOutArgs() const = 0;

  /** \brief . */
  virtual void evalModel( const InArgs<Scalar>& inArgs, const OutArgs<Scalar>& outArgs ) const = 0;

  //@}

  /** \name Virtual functions with default implementations. */
  //@{

  /** \brief Return an optional initial guess for x.
   *
   * If an initial guess is not supported then <tt>return.get()==NULL</tt>.
   */
  virtual Teuchos::RefCountPtr<const VectorBase<Scalar> > get_x_init() const;

  /** \brief Return an optional initial guess for t.
   *
   * If an initial guess is not supported then <tt>return==0.0</tt>.
   */
  virtual ScalarMag get_t_init() const;

  /** \brief If supported, create a <tt>LinearOpWithSolveBase</tt> object for
   * <tt>W</tt> to be evaluated.
   *
   * The default implementation returns <tt>return.get()==NULL</tt>
   * (i.e. implicit solvers are not supported by default).
   */
  virtual Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > create_W() const;

  //@}

};

// //////////////////////////////////
// Helper functions

/** \brief . */
template<class Scalar>
inline 
void eval_f(
  const ModelEvaluator<Scalar>                                    &model
  ,const VectorBase<Scalar>                                       &x
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag   &t
  ,VectorBase<Scalar>                                             *f
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();

  inArgs.set_x(Teuchos::rcp(&x,false));
  if(inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(t);

  outArgs.set_f(Teuchos::rcp(f,false));

  model.evalModel(inArgs,outArgs);

}

/** \brief . */
template<class Scalar>
inline 
void eval_f(
  const ModelEvaluator<Scalar>                                    &model
  ,const VectorBase<Scalar>                                       &x_dot
  ,const VectorBase<Scalar>                                       &x
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag   &t
  ,VectorBase<Scalar>                                             *f
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();

  inArgs.set_x_dot(Teuchos::rcp(&x_dot,false));
  inArgs.set_x(Teuchos::rcp(&x,false));
  if(inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(t);

  outArgs.set_f(Teuchos::rcp(f,false));

  model.evalModel(inArgs,outArgs);

}

/** \brief . */
template<class Scalar>
inline 
void eval_f_W(
  const ModelEvaluator<Scalar>                                    &model
  ,const VectorBase<Scalar>                                       &x_dot
  ,const VectorBase<Scalar>                                       &x
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag   &t
  ,const Scalar                                                   &alpha
  ,const Scalar                                                   &beta
  ,VectorBase<Scalar>                                             *f
  ,LinearOpWithSolveBase<Scalar>                                  *W
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();

  inArgs.set_x_dot(Teuchos::rcp(&x_dot,false));
  inArgs.set_x(Teuchos::rcp(&x,false));
  if(inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(t);
  inArgs.set_alpha(alpha);
  inArgs.set_beta(beta);

  if(f) outArgs.set_f(Teuchos::rcp(f,false));
  if(W) outArgs.set_W(Teuchos::rcp(W,false));

  model.evalModel(inArgs,outArgs);

}

/** \brief . */
template<class Scalar>
inline 
void eval_f_W(
  const ModelEvaluator<Scalar>                                    &model
  ,const VectorBase<Scalar>                                       &x
  ,const Scalar                                                   &beta
  ,VectorBase<Scalar>                                             *f
  ,LinearOpWithSolveBase<Scalar>                                  *W
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();

  inArgs.set_x(Teuchos::rcp(&x,false));
  inArgs.set_beta(beta);

  outArgs.set_f(Teuchos::rcp(f,false));
  if(W) outArgs.set_W(Teuchos::rcp(W,false));

  model.evalModel(inArgs,outArgs);

}


// //////////////////////////////////
// Definitions

// ModelEvaluatorBase::InArgs

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>::InArgs()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Teuchos::ScalarTraits<typename ST::magnitudeType> SMT;
  std::fill_n(&supports_[0],NUM_E_IN_ARGS_MEMBERS,false);
  t_     = SMT::zero();
  alpha_ = ST::zero();
  beta_  = ST::zero();
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x_dot( const Teuchos::RefCountPtr<const VectorBase<Scalar> > &x_dot )
{ assert_supports(IN_ARG_x_dot); x_dot_ = x_dot; }

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_x_dot() const
{ assert_supports(IN_ARG_x_dot); return x_dot_; }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x( const Teuchos::RefCountPtr<const VectorBase<Scalar> > &x )
{ assert_supports(IN_ARG_x); x_ = x; }

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_x() const
{ assert_supports(IN_ARG_x); return x_; }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_t( ScalarMag t )
{ assert_supports(IN_ARG_t); t_ = t; }

template<class Scalar>
typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag
ModelEvaluatorBase::InArgs<Scalar>::get_t() const
{ assert_supports(IN_ARG_t); return t_; }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_alpha( Scalar alpha )
{ assert_supports(IN_ARG_alpha); alpha_ = alpha; }

template<class Scalar>
Scalar ModelEvaluatorBase::InArgs<Scalar>::get_alpha() const
{ assert_supports(IN_ARG_alpha); return alpha_; }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_beta( Scalar beta )
{ assert_supports(IN_ARG_beta); beta_ = beta; }

template<class Scalar>
Scalar ModelEvaluatorBase::InArgs<Scalar>::get_beta() const
{ assert_supports(IN_ARG_beta); return beta_; }

template<class Scalar>
bool ModelEvaluatorBase::InArgs<Scalar>::supports(EInArgsMembers arg) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error,"Error, arg="<<arg<<" is invalid!");
#endif
  return supports_[arg];
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setSupports( EInArgsMembers arg, bool supports )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error,"Error, arg="<<arg<<" is invalid!");
#endif
  supports_[arg] = supports;
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::assert_supports(EInArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"Thyra::ModelEvaluatorBase::InArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(arg), Error, "
    "The argument arg = " << arg << " is not supported!"
    );
}

// ModelEvaluatorBase::OutArgs

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::OutArgs()
{ std::fill_n(&supports_[0],NUM_E_OUT_ARGS_MEMBERS,false); }

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_f( const Teuchos::RefCountPtr<VectorBase<Scalar> > &f )
{ f_ = f; }

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_f() const { return f_; }

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_W( const Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > &W )
{ W_ = W; }

template<class Scalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_W() const { return W_; }

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(EOutArgsMembers arg) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error,"Error, arg="<<arg<<" is invalid!");
#endif
  return supports_[arg];
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports( EOutArgsMembers arg, bool supports )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error,"Error, arg="<<arg<<" is invalid!");
#endif
  supports_[arg] = supports;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(EOutArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(arg), Error, "
    "The argument arg = " << arg << " is not supported!"
    );
}

// ModelEvaluator

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
ModelEvaluator<Scalar>::get_x_init() const
{ return Teuchos::null; }

template<class Scalar>
typename ModelEvaluator<Scalar>::ScalarMag
ModelEvaluator<Scalar>::get_t_init() const
{ return 0.0; }

template<class Scalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
ModelEvaluator<Scalar>::create_W() const
{ return Teuchos::null; }

} // namespace Thyra

#endif // THYRA_MODEL_EVALUATOR_HPP
