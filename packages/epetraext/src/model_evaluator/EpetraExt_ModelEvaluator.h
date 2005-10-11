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
#include "Teuchos_Describable.hpp"

class Epetra_Map;
class Epetra_Vector;
class Epetra_Operator;

namespace EpetraExt {

/** \brief Base interface for evaluating a stateless "model".
 *
 * ToDo: Finish Documentation!
 */
class ModelEvaluator : virtual public Teuchos::Describable {
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
    /** \brief .  */
    int Np() const;
    /** \brief. */
    void set_x_dot( const Teuchos::RefCountPtr<const Epetra_Vector> &x_dot );
    /** \brief. */
    Teuchos::RefCountPtr<const Epetra_Vector> get_x_dot() const;
    /** \brief. */
    void set_x( const Teuchos::RefCountPtr<const Epetra_Vector> &x );
    /** \brief. */
    Teuchos::RefCountPtr<const Epetra_Vector> get_x() const;
    /** \brief. */
    void set_p( int l, const Teuchos::RefCountPtr<const Epetra_Vector> &p_l );
    /** \brief. */
    Teuchos::RefCountPtr<const Epetra_Vector> get_p(int l) const;
    /** \brief. */
    void set_t( double t );
    /** \brief. */
    double get_alpha() const;
    /** \brief. */
    void set_alpha( double alpha );
    /** \brief. */
    double get_beta() const;
    /** \brief. */
    void set_beta( double beta );
    /** \brief. */
    double get_t() const;
    /** \brief. */
    bool supports(EInArgsMembers arg) const;
  protected:
    /** \brief . */
    void _setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void _set_Np(int Np);
    /** \brief . */
    void _setSupports( EInArgsMembers arg, bool supports );
  private:
    // types
    typedef std::vector<Teuchos::RefCountPtr<const Epetra_Vector> > p_t;
    // data
    std::string                                modelEvalDescription_;
    Teuchos::RefCountPtr<const Epetra_Vector>  x_dot_;
    Teuchos::RefCountPtr<const Epetra_Vector>  x_;
    p_t                                        p_;
    double                                     t_;
    double                                     alpha_;
    double                                     beta_;
    bool supports_[NUM_E_IN_ARGS_MEMBERS];
    // functions
    void assert_supports(EInArgsMembers arg) const;
    void assert_l(int l) const;
  };

  /** \brief.  */
  enum EOutArgsMembers {
    OUT_ARG_f
    ,OUT_ARG_W
  };
  static const int NUM_E_OUT_ARGS_MEMBERS=2;

  /** \brief . */
  enum EDerivativeLinearity {
    DERIV_LINEARITY_UNKNOWN      ///< .
    ,DERIV_LINEARITY_CONST       ///< .
    ,DERIV_LINEARITY_NONCONST    ///< .
  };
  /** \brief . */
  enum ERankStatus {
    DERIV_RANK_UNKNOWN       ///< .
    ,DERIV_RANK_FULL         ///< .
    ,DERIV_RANK_DEFICIENT    ///< .
  };

  /** \breif . */
  struct DerivativeProperties {
    /** \breif . */
    EDerivativeLinearity     linearity;
    /** \breif . */
    ERankStatus              rank;
    /** \breif . */
    bool                     supportsAdjoint;
    /** \brief . */
    DerivativeProperties()
      :linearity(DERIV_LINEARITY_UNKNOWN),rank(DERIV_RANK_UNKNOWN),supportsAdjoint(false) {}
    /** \brief . */
    DerivativeProperties(
      EDerivativeLinearity in_linearity, ERankStatus in_rank, bool in_supportsAdjoint
      ):linearity(in_linearity),rank(in_rank),supportsAdjoint(in_supportsAdjoint) {}
  };

  /** \brief . */
  class OutArgs {
  public:
    /** \brief. */
    OutArgs();
    /** \brief .  */
    int Ng() const;
    /** \brief. */
    void set_f( const Teuchos::RefCountPtr<Epetra_Vector> &f );
    /** \brief. */
    Teuchos::RefCountPtr<Epetra_Vector> get_f() const;
    /** \brief Set <tt>g(j)</tt> where <tt>1 <= j && j <= this->Ng()</tt>.  */
    void set_g( int j, const Teuchos::RefCountPtr<Epetra_Vector> &g_j );
    /** \brief Get <tt>g(j)</tt> where <tt>1 <= j && j <= this->Ng()</tt>.  */
    Teuchos::RefCountPtr<Epetra_Vector> get_g(int j) const;
    /** \brief. */
    void set_W( const Teuchos::RefCountPtr<Epetra_Operator> &W );
    /** \brief. */
    Teuchos::RefCountPtr<Epetra_Operator> get_W() const;
    /** \brief . */
    DerivativeProperties get_W_properties() const;
    /** \brief. */
    bool supports(EOutArgsMembers arg) const;
  protected:
    /** \brief . */
    void _setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void _set_Ng(int Ng);
    /** \brief . */
    void _setSupports( EOutArgsMembers arg, bool supports );
    /** \brief . */
    void _set_W_properties( const DerivativeProperties &W_properties );
  private:
    // types
    typedef std::vector<Teuchos::RefCountPtr<Epetra_Vector> > g_t;
    // data
    std::string                            modelEvalDescription_;
    Teuchos::RefCountPtr<Epetra_Vector>    f_;
    g_t                                    g_;
    Teuchos::RefCountPtr<Epetra_Operator>  W_;
    DerivativeProperties                   W_properties_;
    bool supports_[NUM_E_OUT_ARGS_MEMBERS];
    // functions
    void assert_supports(EOutArgsMembers arg) const;
    void assert_j(int j) const;
  };

  //@}

  /** \name Destructor */
  //@{

  /** \brief . */
  virtual ~ModelEvaluator();

  //@}

  /** \name Vector maps */
  //@{

  /** \breif . */
  virtual Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const = 0;

  /** \breif . */
  virtual Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const = 0;

  /** \breif . */
  virtual Teuchos::RefCountPtr<const Epetra_Map> get_p_map(int l) const;

  /** \breif . */
  virtual Teuchos::RefCountPtr<const Epetra_Map> get_g_map(int j) const;

  //@}

  /** \name Initial guesses for variables/parameters */
  //@{

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_x_init() const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_p_init(int l) const;

  /** \brief . */
  virtual double get_t_init() const;

  //@}

  /** \name Bounds for variables/parameters */
  //@{

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_x_lower_bounds() const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_x_upper_bounds() const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_p_lower_bounds(int l) const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_p_upper_bounds(int l) const;

  /** \brief . */
  virtual double get_t_lower_bound() const;

  /** \brief . */
  virtual double get_t_upper_bound() const;

  //@}

  /** \name Factory functions for creating derivative objects */
  //@{

  /** \brief If supported, create a <tt>Epetra_Operator</tt> object for
   * <tt>W</tt> to be evaluated.
   *
   * The default implementation returns <tt>return.get()==NULL</tt>
   * (i.e. implicit solvers are not supported by default).
   */
  virtual Teuchos::RefCountPtr<Epetra_Operator> create_W() const;

  //@}

  /** \name Computational functions */
  //@{

  /** \brief . */
  virtual InArgs createInArgs() const = 0;

  /** \brief . */
  virtual OutArgs createOutArgs() const = 0;

  /** \brief . */
  virtual void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const = 0;

  //@}

protected:

  /** \name Protected types */
  //@{

  /** \brief . */
  class InArgsSetup : public InArgs {
  public:
    /** \brief . */
    void setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void set_Np(int Np);
    /** \brief . */
    void setSupports( EInArgsMembers arg, bool supports = true );
  };

  /** \brief . */
  class OutArgsSetup : public OutArgs {
  public:
    /** \brief . */
    void setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void set_Ng(int Ng);
    /** \brief . */
    void setSupports( EOutArgsMembers arg, bool supports = true );
    /** \brief . */
    void set_W_properties( const DerivativeProperties &W_properties );
  };

  //@}

};

// ///////////////////////////
// Inline Functions

//
// ModelEvaluator::InArgs
//

inline
int ModelEvaluator::InArgs::Np() const
{ return p_.size(); }

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
void ModelEvaluator::InArgs::set_p( int l, const Teuchos::RefCountPtr<const Epetra_Vector> &p_l )
{ assert_l(l); p_[l-1] = p_l; }

inline
Teuchos::RefCountPtr<const Epetra_Vector> ModelEvaluator::InArgs::get_p(int l) const
{ assert_l(l); return p_[l-1]; }

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
void ModelEvaluator::InArgs::_setModelEvalDescription( const std::string &modelEvalDescription )
{
  modelEvalDescription_ = modelEvalDescription;
}

inline
void ModelEvaluator::InArgs::_set_Np(int Np)
{
  p_.resize(Np);
}

//
// ModelEvaluator::OutArgs
//

inline
int ModelEvaluator::OutArgs::Ng() const
{ 
  return g_.size();
}

inline
void ModelEvaluator::OutArgs::set_f( const Teuchos::RefCountPtr<Epetra_Vector> &f ) { f_ = f; }

inline
Teuchos::RefCountPtr<Epetra_Vector> ModelEvaluator::OutArgs::get_f() const { return f_; }

inline
void ModelEvaluator::OutArgs::set_g( int j, const Teuchos::RefCountPtr<Epetra_Vector> &g_j )
{
  assert_j(j);
  g_[j-1] = g_j;
}

inline
Teuchos::RefCountPtr<Epetra_Vector> ModelEvaluator::OutArgs::get_g(int j) const
{
  assert_j(j);
  return g_[j-1];
}

inline
void ModelEvaluator::OutArgs::set_W( const Teuchos::RefCountPtr<Epetra_Operator> &W ) { W_ = W; }

inline
Teuchos::RefCountPtr<Epetra_Operator> ModelEvaluator::OutArgs::get_W() const { return W_; }

inline
ModelEvaluator::DerivativeProperties ModelEvaluator::OutArgs::get_W_properties() const
{
  return W_properties_;
}

inline
void ModelEvaluator::OutArgs::_setModelEvalDescription( const std::string &modelEvalDescription )
{
  modelEvalDescription_ = modelEvalDescription;
}

inline
void ModelEvaluator::OutArgs::_set_Ng(int Ng)
{
  g_.resize(Ng);
}

inline
void ModelEvaluator::OutArgs::_set_W_properties( const DerivativeProperties &W_properties )
{
  W_properties_ = W_properties;
}

//
// ModelEvaluatorBase::InArgsSetup
//

inline
void ModelEvaluator::InArgsSetup::setModelEvalDescription( const std::string &modelEvalDescription )
{
  this->_setModelEvalDescription(modelEvalDescription);
}

inline
void ModelEvaluator::InArgsSetup::set_Np(int Np)
{ this->_set_Np(Np); }

inline
void ModelEvaluator::InArgsSetup::setSupports( EInArgsMembers arg, bool supports )
{ this->_setSupports(arg,supports); }

//
// ModelEvaluatorBase::OutArgsSetup
//

inline
void ModelEvaluator::OutArgsSetup::setModelEvalDescription( const std::string &modelEvalDescription )
{
  this->_setModelEvalDescription(modelEvalDescription);
}

inline
void ModelEvaluator::OutArgsSetup::set_Ng(int Ng)
{ this->_set_Ng(Ng); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsMembers arg, bool supports )
{ this->_setSupports(arg,supports); }

inline
void ModelEvaluator::OutArgsSetup::set_W_properties( const DerivativeProperties &W_properties )
{ this->_set_W_properties(W_properties); }

} // namespace EpetraExt

#endif // EPETRA_EXT_MODEL_EVALUATOR_HPP
