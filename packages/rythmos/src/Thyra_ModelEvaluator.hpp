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

#include "Thyra_OperatorVectorTypes.hpp"

namespace Thyra {

/** \brief . */
class ModelEvaluatorBase {
public:

  /** \name Public types */
  //@{

  /** \brief.  */
  enum EInArgsMembers {
    IN_ARG_x
    ,IN_ARG_t
  };
  static const int NUM_E_IN_ARGS_MEMBERS=2;

  /** \brief . */
  template<class Scalar>
  class InArgs {
  public:
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    InArgs();
    void set_x( const Teuchos::RefCountPtr<const VectorBase<Scalar> > &x );
    Teuchos::RefCountPtr<const VectorBase<Scalar> > get_x() const;
    void set_t( ScalarMag t );
    ScalarMag get_t() const;
    bool supports(EInArgsMembers arg);
  protected:
    /** \brief . */
    void _setSupports( EInArgsMembers arg, bool supports );
  private:
    Teuchos::RefCountPtr<const VectorBase<Scalar> >  x_;
    Scalar                                           t_;
    bool supports_[NUM_E_IN_ARGS_MEMBERS];
  };

  /** \brief.  */
  enum EOutArgsMembers {
    OUT_ARG_f
  };
  static const int NUM_E_OUT_ARGS_MEMBERS=1;

  /** \brief . */
  template<class Scalar>
  class OutArgs {
  public:
    OutArgs();
    void set_f( const Teuchos::RefCountPtr<VectorBase<Scalar> > &f );
    Teuchos::RefCountPtr<VectorBase<Scalar> > get_f() const;
    bool supports(EOutArgsMembers arg);
  protected:
    /** \brief . */
    void _setSupports( EOutArgsMembers arg, bool supports );
  private:
    Teuchos::RefCountPtr<VectorBase<Scalar> >  f_;
    bool supports_[NUM_E_OUT_ARGS_MEMBERS];
  };

  //@}

protected:

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

  //@}

};

// //////////////////////////////////
// Definitions

// ModelEvaluatorBase::InArgs

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>::InArgs()
{ std::fill_n(&supports_[0],NUM_E_IN_ARGS_MEMBERS,false); }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x( const Teuchos::RefCountPtr<const VectorBase<Scalar> > &x )
{ x_ = x; }

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_x() const
{ return x_; }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_t( ScalarMag t )
{ t_ = t; }

template<class Scalar>
typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag
ModelEvaluatorBase::InArgs<Scalar>::get_t() const
{ return t_; }

template<class Scalar>
bool ModelEvaluatorBase::InArgs<Scalar>::supports(EInArgsMembers arg)
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
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(EOutArgsMembers arg)
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

// ModelEvaluator

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
ModelEvaluator<Scalar>::get_x_init() const
{ return Teuchos::null; }

template<class Scalar>
typename ModelEvaluator<Scalar>::ScalarMag
ModelEvaluator<Scalar>::get_t_init() const
{ return 0.0; }

} // namespace Thyra

#endif // THYRA_MODEL_EVALUATOR_HPP
