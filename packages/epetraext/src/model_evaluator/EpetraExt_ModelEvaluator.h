// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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

#ifndef EPETRA_EXT_MODEL_EVALUATOR_HPP
#define EPETRA_EXT_MODEL_EVALUATOR_HPP

#include "Teuchos_RefCountPtr.hpp"

class Epetra_Map;
class Epetra_Vector;
class Epetra_Operator;

namespace EpetraExt {

/** \brief Base interface for evaluating a stateless "model".
 *
 * ToDo: Finish Documentation!
 */
class ModelEvaluator {
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
  class InArgs {
  public:
    /** \brief. */
    InArgs();
    void set_x_dot( const Teuchos::RefCountPtr<const Epetra_Vector> &x_dot );
    Teuchos::RefCountPtr<const Epetra_Vector> get_x_dot() const;
    void set_x( const Teuchos::RefCountPtr<const Epetra_Vector> &x );
    Teuchos::RefCountPtr<const Epetra_Vector> get_x() const;
    void set_t( double t );
    double get_alpha() const;
    void set_alpha( double alpha );
    double get_beta() const;
    void set_beta( double beta );
    double get_t() const;
    bool supports(EInArgsMembers arg) const;
  protected:
    /** \brief . */
    void _setSupports( EInArgsMembers arg, bool supports );
  private:
    Teuchos::RefCountPtr<const Epetra_Vector>  x_dot_;
    Teuchos::RefCountPtr<const Epetra_Vector>  x_;
    double                                     t_;
    double                                     alpha_;
    double                                     beta_;
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
  class OutArgs {
  public:
    /** \brief. */
    OutArgs();
    void set_f( const Teuchos::RefCountPtr<Epetra_Vector> &f );
    Teuchos::RefCountPtr<Epetra_Vector> get_f() const;
    void set_W( const Teuchos::RefCountPtr<Epetra_Operator> &W );
    Teuchos::RefCountPtr<Epetra_Operator> get_W() const;
    bool supports(EOutArgsMembers arg) const;
  protected:
    /** \brief . */
    void _setSupports( EOutArgsMembers arg, bool supports );
  private:
    Teuchos::RefCountPtr<Epetra_Vector>    f_;
    Teuchos::RefCountPtr<Epetra_Operator>  W_;
    bool supports_[NUM_E_OUT_ARGS_MEMBERS];
    void assert_supports(EOutArgsMembers arg) const;
  };

  //@}

  /** \brief . */
  virtual ~ModelEvaluator() {}

  /** \name Pure virtual functions that must be overridden by subclasses. */
  //@{

  /** \breif . */
  virtual Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const = 0;

  /** \breif . */
  virtual Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const = 0;

  /** \brief . */
  virtual InArgs createInArgs() const = 0;

  /** \brief . */
  virtual OutArgs createOutArgs() const = 0;

  /** \brief . */
  virtual void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const = 0;

  //@}

  /** \name Virtual functions with default implementations. */
  //@{

  /** \brief Return an optional initial guess for x.
   *
   * If an initial guess is not supported then <tt>return.get()==NULL</tt>.
   */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_x_init() const;

  /** \brief Return an optional initial guess for t.
   *
   * If an initial guess is not supported then <tt>return==0.0</tt>.
   */
  virtual double get_t_init() const;

  /** \brief If supported, create a <tt>Epetra_Operator</tt> object for
   * <tt>W</tt> to be evaluated.
   *
   * The default implementation returns <tt>return.get()==NULL</tt>
   * (i.e. implicit solvers are not supported by default).
   */
  virtual Teuchos::RefCountPtr<Epetra_Operator> create_W() const;

  //@}

protected:

  /** \name Protected types */
  //@{

  /** \brief . */
  class InArgsSetup : public InArgs {
  public:
    /** \brief . */
    void setSupports( EInArgsMembers arg, bool supports = true ) { _setSupports(arg,supports); }
  };

  /** \brief . */
  class OutArgsSetup : public OutArgs {
  public:
    /** \brief . */
    void setSupports( EOutArgsMembers arg, bool supports = true ) { _setSupports(arg,supports); }
  };

  //@}

};

// ///////////////////////////
// Definitions

// ModelEvaluator::InArgs

inline
ModelEvaluator::InArgs::InArgs()
{
  std::fill_n(&supports_[0],NUM_E_IN_ARGS_MEMBERS,false);
  t_     = 0.0;
  alpha_ = 0.0;
  beta_  = 0.0;
}

inline
void ModelEvaluator::InArgs::set_x_dot( const Teuchos::RefCountPtr<const Epetra_Vector> &x_dot )
{ assert_supports(IN_ARG_x_dot); x_dot_ = x_dot; }

inline
Teuchos::RefCountPtr<const Epetra_Vector> ModelEvaluator::InArgs::get_x_dot() const
{ assert_supports(IN_ARG_x_dot); return x_dot_; }

inline
void ModelEvaluator::InArgs::set_x( const Teuchos::RefCountPtr<const Epetra_Vector> &x )
{ assert_supports(IN_ARG_x); x_ = x; }

inline
Teuchos::RefCountPtr<const Epetra_Vector> ModelEvaluator::InArgs::get_x() const
{ assert_supports(IN_ARG_x); return x_; }

inline
void ModelEvaluator::InArgs::set_t( double t )
{ assert_supports(IN_ARG_t); t_ = t; }

inline
double ModelEvaluator::InArgs::get_t() const
{ assert_supports(IN_ARG_t); return t_; }

inline
void ModelEvaluator::InArgs::set_alpha( double alpha )
{ assert_supports(IN_ARG_alpha); alpha_ = alpha; }

inline
double ModelEvaluator::InArgs::get_alpha() const
{ assert_supports(IN_ARG_alpha); return alpha_; }

inline
void ModelEvaluator::InArgs::set_beta( double beta )
{ assert_supports(IN_ARG_beta); beta_ = beta; }

inline
double ModelEvaluator::InArgs::get_beta() const
{ assert_supports(IN_ARG_beta); return beta_; }

inline
bool ModelEvaluator::InArgs::supports(EInArgsMembers arg) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error,"Error, arg="<<arg<<" is invalid!");
#endif
  return supports_[arg];
}

inline
void ModelEvaluator::InArgs::_setSupports( EInArgsMembers arg, bool supports )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error,"Error, arg="<<arg<<" is invalid!");
#endif
  supports_[arg] = supports;
}

inline
void ModelEvaluator::InArgs::assert_supports(EInArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"EpetraExt::ModelEvaluator::InArgs::assert_supports(arg), Error, "
    "The argument arg = " << arg << " is not supported!"
    );
}

// ModelEvaluator::OutArgs

inline
ModelEvaluator::OutArgs::OutArgs()
{
  std::fill_n(&supports_[0],NUM_E_OUT_ARGS_MEMBERS,false);
}

inline
void ModelEvaluator::OutArgs::set_f( const Teuchos::RefCountPtr<Epetra_Vector> &f ) { f_ = f; }

inline
Teuchos::RefCountPtr<Epetra_Vector> ModelEvaluator::OutArgs::get_f() const { return f_; }

inline
void ModelEvaluator::OutArgs::set_W( const Teuchos::RefCountPtr<Epetra_Operator> &W ) { W_ = W; }

inline
Teuchos::RefCountPtr<Epetra_Operator> ModelEvaluator::OutArgs::get_W() const { return W_; }

inline
bool ModelEvaluator::OutArgs::supports(EOutArgsMembers arg) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error,"Error, arg="<<arg<<" is invalid!");
#endif
  return supports_[arg];
}

inline
void ModelEvaluator::OutArgs::_setSupports( EOutArgsMembers arg, bool supports )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error,"Error, arg="<<arg<<" is invalid!");
#endif
  supports_[arg] = supports;
}

inline
void ModelEvaluator::OutArgs::assert_supports(EOutArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(arg), Error, "
    "The argument arg = " << arg << " is not supported!"
    );
}

// ModelEvaluator

inline
Teuchos::RefCountPtr<const Epetra_Vector>
ModelEvaluator::get_x_init() const
{ return Teuchos::null; }

inline
double ModelEvaluator::get_t_init() const
{ return 0.0; }

inline
Teuchos::RefCountPtr<Epetra_Operator>
ModelEvaluator::create_W() const
{ return Teuchos::null; }

} // namespace EpetraExt

#endif // EPETRA_EXT_MODEL_EVALUATOR_HPP
