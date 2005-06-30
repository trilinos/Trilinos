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
    IN_ARG_x
    ,IN_ARG_t
  };
  static const int NUM_E_IN_ARGS_MEMBERS=2;

  /** \brief . */
  class InArgs {
  public:
    /** \brief. */
    InArgs();
    void set_x( const Teuchos::RefCountPtr<const Epetra_Vector> &x );
    Teuchos::RefCountPtr<const Epetra_Vector> get_x() const;
    void set_t( double t );
    double get_t() const;
    bool supports(EInArgsMembers arg);
  protected:
    /** \brief . */
    void _setSupports( EInArgsMembers arg, bool supports );
  private:
    Teuchos::RefCountPtr<const Epetra_Vector>  x_;
    double                                     t_;
    bool supports_[NUM_E_IN_ARGS_MEMBERS];
  };

  /** \brief.  */
  enum EOutArgsMembers {
    OUT_ARG_f
  };
  static const int NUM_E_OUT_ARGS_MEMBERS=1;

  /** \brief . */
  class OutArgs {
  public:
    /** \brief. */
    OutArgs();
    void set_f( const Teuchos::RefCountPtr<Epetra_Vector> &f );
    Teuchos::RefCountPtr<Epetra_Vector> get_f() const;
    bool supports(EOutArgsMembers arg);
  protected:
    /** \brief . */
    void _setSupports( EOutArgsMembers arg, bool supports );
  private:
    Teuchos::RefCountPtr<Epetra_Vector>  f_;
    bool supports_[NUM_E_OUT_ARGS_MEMBERS];
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

ModelEvaluator::InArgs::InArgs() { std::fill_n(&supports_[0],NUM_E_IN_ARGS_MEMBERS,false); }

void ModelEvaluator::InArgs::set_x( const Teuchos::RefCountPtr<const Epetra_Vector> &x ) { x_ = x; }

Teuchos::RefCountPtr<const Epetra_Vector> ModelEvaluator::InArgs::get_x() const { return x_; }

void ModelEvaluator::InArgs::set_t( double t ) { t_ = t; }

double ModelEvaluator::InArgs::get_t() const { return t_; }

bool ModelEvaluator::InArgs::supports(EInArgsMembers arg)
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error,"Error, arg="<<arg<<" is invalid!");
#endif
  return supports_[arg];
}

void ModelEvaluator::InArgs::_setSupports( EInArgsMembers arg, bool supports )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error,"Error, arg="<<arg<<" is invalid!");
#endif
  supports_[arg] = supports;
}

// ModelEvaluator::OutArgs

ModelEvaluator::OutArgs::OutArgs() { std::fill_n(&supports_[0],NUM_E_OUT_ARGS_MEMBERS,false); }

void ModelEvaluator::OutArgs::set_f( const Teuchos::RefCountPtr<Epetra_Vector> &f ) { f_ = f; }

Teuchos::RefCountPtr<Epetra_Vector> ModelEvaluator::OutArgs::get_f() const { return f_; }

bool ModelEvaluator::OutArgs::supports(EOutArgsMembers arg)
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error,"Error, arg="<<arg<<" is invalid!");
#endif
  return supports_[arg];
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsMembers arg, bool supports )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error,"Error, arg="<<arg<<" is invalid!");
#endif
  supports_[arg] = supports;
}

// ModelEvaluator

Teuchos::RefCountPtr<const Epetra_Vector>
ModelEvaluator::get_x_init() const
{ return Teuchos::null; }

double ModelEvaluator::get_t_init() const
{ return 0.0; }

} // namespace EpetraExt

#endif // EPETRA_EXT_MODEL_EVALUATOR_HPP
