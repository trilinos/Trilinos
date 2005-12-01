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

#include "Teuchos_Describable.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"

#include "Teuchos_Polynomial.hpp"
#include "Thyra_PolynomialVectorTraits.hpp"

namespace Thyra {

/** \brief Base subclass for <tt>ModelEvaluator</tt> that defines some basic
 * types.
 *
 * ToDo: Finish documentation!
 */
class ModelEvaluatorBase : virtual public Teuchos::Describable {
public:

  /** \name Public types */
  //@{

  /** \brief .  */
  enum EInArgsMembers {
    IN_ARG_x_dot ///< .
    ,IN_ARG_x ///< .
    ,IN_ARG_x_dot_poly ///< .
    ,IN_ARG_x_poly ///< .
    ,IN_ARG_t ///< .
    ,IN_ARG_alpha ///< .
    ,IN_ARG_beta ///< .
  };
  /** \brief .  */
  static const int NUM_E_IN_ARGS_MEMBERS=7;

  /** \brief . */
  template<class Scalar>
  class InArgs {
  public:
    /** \brief .  */
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    /** \brief .  */
    InArgs();
    /** \brief .  */
    int Np() const;
    /** \brief .  */
    void set_x_dot( const Teuchos::RefCountPtr<const VectorBase<Scalar> > &x_dot );
    /** \brief .  */
    Teuchos::RefCountPtr<const VectorBase<Scalar> > get_x_dot() const;
    /** \brief .  */
    void set_x( const Teuchos::RefCountPtr<const VectorBase<Scalar> > &x );
    /** \brief .  */
    Teuchos::RefCountPtr<const VectorBase<Scalar> > get_x() const;
    /** \brief .  */
    void set_x_poly( const Teuchos::RefCountPtr<const Teuchos::Polynomial< VectorBase<Scalar> > > &x_poly );
    /** \brief .  */
    Teuchos::RefCountPtr<const Teuchos::Polynomial< VectorBase<Scalar> > > get_x_poly() const;
    /** \brief .  */
    void set_x_dot_poly( const Teuchos::RefCountPtr<const Teuchos::Polynomial< VectorBase<Scalar> > > &x_dot_poly );
    /** \brief .  */
    Teuchos::RefCountPtr<const Teuchos::Polynomial< VectorBase<Scalar> > > get_x_dot_poly() const;
    /** \brief Set <tt>p(l)</tt> where <tt>1 <= l && l <= this->Np()</tt>.  */
    void set_p( int l, const Teuchos::RefCountPtr<const VectorBase<Scalar> > &p_l );
    /** \brief Get <tt>p(l)</tt> where <tt>1 <= l && l <= this->Np()</tt>.  */
    Teuchos::RefCountPtr<const VectorBase<Scalar> > get_p(int l) const;
    /** \brief .  */
    void set_t( ScalarMag t );
    /** \brief .  */
    ScalarMag get_t() const;
    /** \brief .  */
    void set_alpha( Scalar alpha );
    /** \brief .  */
    Scalar get_alpha() const;
    /** \brief .  */
    void set_beta( Scalar beta );
    /** \brief .  */
    Scalar get_beta() const;
    /** \brief .  */
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
    typedef std::vector<Teuchos::RefCountPtr<const VectorBase<Scalar> > > p_t;
    // data
    std::string                                      modelEvalDescription_;
    Teuchos::RefCountPtr<const VectorBase<Scalar> >  x_dot_;
    Teuchos::RefCountPtr<const VectorBase<Scalar> >  x_;
    Teuchos::RefCountPtr<const Teuchos::Polynomial< VectorBase<Scalar> > > x_dot_poly_;
    Teuchos::RefCountPtr<const Teuchos::Polynomial< VectorBase<Scalar> > > x_poly_;
    p_t                                              p_;
    ScalarMag                                        t_;
    Scalar                                           alpha_;
    Scalar                                           beta_;
    bool supports_[NUM_E_IN_ARGS_MEMBERS];
    // functions
    void assert_supports(EInArgsMembers arg) const;
    void assert_l(int l) const;
  };

  /** \brief . */
  enum EDerivativeMultiVectorOrientation {
    DERIV_MV_BY_COL           ///< .
    ,DERIV_TRANS_MV_BY_ROW    ///< .
  };

  /** \brief . */
  enum EDerivativeLinearOp { DERIV_LINEAR_OP};

  /** \brief . */
  class DerivativeSupport {
  public:
    /** \brief . */
    DerivativeSupport()
      :supportsLinearOp_(false), supportsMVByCol_(false), supportsTransMVByRow_(false)
      {}
    /** \brief . */
    DerivativeSupport( EDerivativeLinearOp )
      :supportsLinearOp_(true), supportsMVByCol_(false), supportsTransMVByRow_(false)
      {}
    /** \brief . */
    DerivativeSupport( EDerivativeMultiVectorOrientation mvOrientation )
      :supportsLinearOp_(false), supportsMVByCol_(mvOrientation==DERIV_MV_BY_COL)
      ,supportsTransMVByRow_(mvOrientation==DERIV_TRANS_MV_BY_ROW)
      {}
    /** \brief . */
    DerivativeSupport& plus(EDerivativeLinearOp)
      { supportsLinearOp_ = true; return *this; }
    /** \brief . */
    DerivativeSupport& plus(EDerivativeMultiVectorOrientation mvOrientation)
      {
        switch(mvOrientation) {
          case DERIV_MV_BY_COL: supportsMVByCol_ = true; break;
          case DERIV_TRANS_MV_BY_ROW: supportsTransMVByRow_ = true; break;
          default: TEST_FOR_EXCEPT(true);
        }
        return *this;
      }
    /** \brief . */
    bool none() const
      { return ( !supportsLinearOp_ && !supportsMVByCol_ && !supportsTransMVByRow_ ); }
    /** \brief . */
    bool supports(EDerivativeLinearOp) const
      { return supportsLinearOp_; }
    /** \brief . */
    bool supports(EDerivativeMultiVectorOrientation mvOrientation) const
      {
        switch(mvOrientation) {
          case DERIV_MV_BY_COL: return supportsMVByCol_;
          case DERIV_TRANS_MV_BY_ROW: return supportsTransMVByRow_;
          default: TEST_FOR_EXCEPT(true);
        }
        return false; // Will never be called!
      }
  private:
    bool supportsLinearOp_;
    bool supportsMVByCol_;
    bool supportsTransMVByRow_;
  public:
  };
  
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

  /** \brief . */
  struct DerivativeProperties {
    /** \brief . */
    EDerivativeLinearity     linearity;
    /** \brief . */
    ERankStatus              rank;
    /** \brief . */
    bool                     supportsAdjoint;
    /** \brief . */
    DerivativeProperties()
      :linearity(DERIV_LINEARITY_UNKNOWN),rank(DERIV_RANK_UNKNOWN),supportsAdjoint(false) {}
    /** \brief . */
    DerivativeProperties(
      EDerivativeLinearity in_linearity, ERankStatus in_rank, bool in_supportsAdjoint
      ):linearity(in_linearity),rank(in_rank),supportsAdjoint(in_supportsAdjoint) {}
  };

  /** \brief Simple aggregate class for a derivative object represented as a
   * column-wise multi-vector or its transpose as a row-wise multi-vector.
   */
  template<class Scalar>
  class DerivativeMultiVector {
  public:
    /** \brief . */
    DerivativeMultiVector() {}
    /** \brief . */
    DerivativeMultiVector(
      const Teuchos::RefCountPtr<MultiVectorBase<Scalar> >  &mv
      ,const EDerivativeMultiVectorOrientation              orientation = DERIV_MV_BY_COL
      ) : mv_(mv), orientation_(orientation) {}
    /** \brief . */
    void changeOrientation( const EDerivativeMultiVectorOrientation orientation )
      { orientation_ = orientation; };
    /** \brief . */
    Teuchos::RefCountPtr<MultiVectorBase<Scalar> > getMultiVector() const
      { return mv_; }
    /** \brief . */
    EDerivativeMultiVectorOrientation getOrientation() const
      { return orientation_; }
  private:
    Teuchos::RefCountPtr<MultiVectorBase<Scalar> >   mv_;
    EDerivativeMultiVectorOrientation                orientation_;
  };

  /** \brief Simple aggregate class that stores a derivative object
   * as a general linear operator or as a multi-vector.
   */
  template<class Scalar>
  class Derivative {
  public:
    /** \brief . */
    Derivative() {}
    /** \brief . */
    Derivative( const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &lo )
      : lo_(lo) {}
    /** \brief . */
    Derivative( const DerivativeMultiVector<Scalar> &dmv )
      : dmv_(dmv) {}
    /** \brief . */
    bool isEmpty() const
      { return ( lo_.get()==NULL && dmv_.getMultiVector().get()==NULL ); }
    /** \brief . */
    Teuchos::RefCountPtr<LinearOpBase<Scalar> > getLinearOp() const
      { return lo_; }
    /** \brief . */
    DerivativeMultiVector<Scalar> getDerivativeMultiVector() const
      { return dmv_; }
  private:
    Teuchos::RefCountPtr<LinearOpBase<Scalar> >   lo_;
    DerivativeMultiVector<Scalar>                 dmv_;
  };

  /** \brief .  */
  enum EOutArgsMembers {
    OUT_ARG_f       ///< .
    ,OUT_ARG_W      ///< .
    ,OUT_ARG_f_poly ///< .
  };
  /** \brief .  */
  static const int NUM_E_OUT_ARGS_MEMBERS=3;

  /** \brief . */
  enum EOutArgsDfDp {
    OUT_ARG_DfDp   ///< .
  };

  /** \brief . */
  enum EOutArgsDgDx {
    OUT_ARG_DgDx   ///< .
  };

  /** \brief . */
  enum EOutArgsDgDp {
    OUT_ARG_DgDp   ///< .
  };
  
  /** \brief . */
  template<class Scalar>
  class OutArgs {
  public:
    /** \brief .  */
    OutArgs();
    /** \brief .  */
    int Np() const;
    /** \brief .  */
    int Ng() const;
    /** \brief .  */
    bool supports(EOutArgsMembers arg) const;
    /** \brief <tt>1 <= l && l <= Np()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDfDp arg, int l) const;
    /** \brief <tt>1 <= j && j <= Ng()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDx arg, int j) const;
    /** \brief <tt>1 <= j && j <= Ng()</tt> and <tt>1 <= l && l <= Np()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDp arg, int j, int l) const;
    /** \brief .  */
    void set_f( const Teuchos::RefCountPtr<VectorBase<Scalar> > &f );
    /** \brief .  */
    Teuchos::RefCountPtr<VectorBase<Scalar> > get_f() const;
    /** \brief .  */
    void set_g( int j, const Teuchos::RefCountPtr<VectorBase<Scalar> > &g_j );
    /** \brief .  */
    Teuchos::RefCountPtr<VectorBase<Scalar> > get_g(int j) const;
    /** \brief .  */
    void set_W( const Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > &W );
    /** \brief .  */
    Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > get_W() const;
    /** \brief . */
    DerivativeProperties get_W_properties() const;
    /** \brief .  */
    void set_DfDp(int l,  const Derivative<Scalar> &DfDp_l);
    /** \brief .  */
    Derivative<Scalar> get_DfDp(int l) const;
    /** \brief . */
    DerivativeProperties get_DfDp_properties(int l) const;
    /** \brief .  */
    void set_DgDx(int j, const Derivative<Scalar> &DgDx_j);
    /** \brief .  */
    Derivative<Scalar> get_DgDx(int j) const;
    /** \brief . */
    DerivativeProperties get_DgDx_properties(int j) const;
    /** \brief .  */
    void set_DgDp( int j, int l, const Derivative<Scalar> &DgDp_j_l );
    /** \brief .  */
    Derivative<Scalar> get_DgDp(int j, int l) const;
    /** \brief . */
    DerivativeProperties get_DgDp_properties(int j, int l) const;
    /** \brief .  */
    void set_f_poly( const Teuchos::RefCountPtr<Teuchos::Polynomial< VectorBase<Scalar> > > &f_poly );
    /** \brief .  */
    Teuchos::RefCountPtr<Teuchos::Polynomial< VectorBase<Scalar> > > get_f_poly() const;
  protected:
    /** \brief . */
    void _setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void _set_Np_Ng(int Np, int Ng);
    /** \brief . */
    void _setSupports( EOutArgsMembers arg, bool supports );
    /** \brief . */
    void _setSupports( EOutArgsDfDp arg, int l, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDx arg, int j, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDp arg, int j, int l, const DerivativeSupport& );
    /** \brief . */
    void _set_W_properties( const DerivativeProperties &properties );
    /** \brief . */
    void _set_DfDp_properties( int l, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDx_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDp_properties( int j, int l, const DerivativeProperties &properties );
  private:
    // types
    typedef std::vector<Teuchos::RefCountPtr<VectorBase<Scalar> > >     g_t;
    typedef std::vector<Derivative<Scalar> >                            deriv_t;
    typedef std::vector<DerivativeProperties>                           deriv_properties_t;
    typedef std::vector<DerivativeSupport>                              supports_t;
    // data
    std::string                                           modelEvalDescription_;
    bool                                                  supports_[NUM_E_OUT_ARGS_MEMBERS];
    supports_t                                            supports_DfDp_;   // Np
    supports_t                                            supports_DgDx_;   // Ng
    supports_t                                            supports_DgDp_;   // Ng x Np
    Teuchos::RefCountPtr<VectorBase<Scalar> >             f_;
    g_t                                                   g_;               // Ng
    Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >  W_;
    DerivativeProperties                                  W_properties_;
    deriv_t                                               DfDp_;            // Np
    deriv_properties_t                                    DfDp_properties_; // Np
    deriv_t                                               DgDx_;            // Ng
    deriv_properties_t                                    DgDx_properties_; // Ng
    deriv_t                                               DgDp_;            // Ng x Np
    deriv_properties_t                                    DgDp_properties_; // Ng x Np
    Teuchos::RefCountPtr<Teuchos::Polynomial< VectorBase<Scalar> > > f_poly_;
    // functions
    void assert_supports(EOutArgsMembers arg) const;
    void assert_supports(EOutArgsDfDp arg, int l) const;
    void assert_supports(EOutArgsDgDx arg, int j) const;
    void assert_supports(EOutArgsDgDp arg, int j, int l) const;
    void assert_l(int l) const;
    void assert_j(int j) const;
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
    void setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void set_Np(int Np);
    /** \brief . */
    void setSupports( EInArgsMembers arg, bool supports = true );
  };

  /** \brief . */
  template<class Scalar>
  class OutArgsSetup : public OutArgs<Scalar> {
  public:
    /** \brief . */
    void setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void set_Np_Ng(int Np, int Ng);
    /** \brief . */
    void setSupports( EOutArgsMembers arg, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgsDfDp arg, int l, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDx arg, int j, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDp arg, int j, int l, const DerivativeSupport& );
    /** \brief . */
    void set_W_properties( const DerivativeProperties &properties );
    /** \brief . */
    void set_DfDp_properties( int l, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDx_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDp_properties( int j, int l, const DerivativeProperties &properties );
  };

  //@}

};

/** \defgroup Thyra_MEB_helper_functions_grp Helper functions for Thyra::ModelEvaluatorBase.
 *
 */
//@{

/** \brief . */
std::string toString(ModelEvaluatorBase::EInArgsMembers);

/** \brief . */
std::string toString(ModelEvaluatorBase::EOutArgsMembers);

//@}

/** \brief Base interface for evaluating a stateless "model" that can be
 * mapped into a number of different types of problems.
 *
 * \section Thyra_ME_outline_sec Outline
 *
 * <ul>
 * <li>\ref Thyra_ME_intro_sec
 * <li>\ref Thyra_ME_problem_types_sec
 *     <ul>
 *     <li>\ref Thyra_ME_nonlinear_equations_sec
 *     <li>\ref Thyra_ME_explicit_ode_sec
 *     <li>\ref Thyra_ME_implicit_dae_sec
 *     <li>\ref Thyra_ME_unconstrained_optimization_sec
 *     <li>\ref Thyra_ME_equality_constrained_optimization_sec
 *     </ul>
 * <li>\ref Thyra_ME_derivatives_sec
 * <li>\ref Thyra_ME_dev_notes_sec
 * </ul>
 *
 * \section Thyra_ME_intro_sec Introduction
 *
 * The model represented by this interface is composed of the following
 * general functions:
 *
 * <ul>
 *
 * <li><b>State vector function:</b>
 *
 * <tt>(x_dot,x,{p(l)},t}) -> f</tt>
 *
 * <li><b>Auxiliary response vector functions:</b>
 *
 * <tt>(x_dot,x,{p(l)},t}) -> g(j)</tt>,
 *
 * for <tt>j=1...Ng</tt>
 *
 * </ul>
 *
 * given the general input variables/parameters:
 *
 * <ul>
 *
 * <li><b>State variables vector:</b>
 *
 * <tt>x</tt>
 *
 * <li><b>State variables derivative w.r.t. <tt>t</tt> vector:</b>
 *
 * <tt>x_dot</tt>
 *
 * <li><b>Auxiliary parameter vectors:</b>
 *
 * <tt>p(l)</tt>,
 *
 * for <tt>l=1...Np</tt>
 *
 * <li>Time point (or some other independent variable):
 *
 * <tt>t</tt>
 *
 * </ul>
 *
 * Above, the notation <tt>{p(l)}</tt> is shorthand for the set of parameter
 * vectors <tt>{ p(1), p(2), ..., p(Np) }</tt>.
 *
 * All of the above variables/parameters and functions are represented as
 * abstract <tt>Thyra::VectorBase</tt> objects.  The vector spaces associated
 * with these vector quantities are returned by <tt>get_space_x()</tt>,
 * <tt>get_p_space(int)</tt>, <tt>get_space_f()</tt>, and
 * <tt>get_g_space(int)</tt>.
 *
 * These functions all get computed at the same time in the function
 * <tt>evalModel()</tt>
 *
 * \section Thyra_ME_problem_types_sec Problem Types
 *
 * There are a number of different types of mathematical problems that can be
 * formulated using this interface.
 *
 * \subsection Thyra_ME_nonlinear_equations_sec Nonlinear Equations
 *
 
 \verbatim

  f(x) = 0.

 \endverbatim 
 
 * \subsection Thyra_ME_explicit_ode_sec Explicit ODEs
 *
 
 \verbatim

  x_dot = f(x,t).

 \endverbatim 

 * Here the argument <tt>t</tt> may or may not be accepted by <tt>*this</tt>.
 *
 * \subsection Thyra_ME_implicit_dae_sec Implicit ODEs or DAEs
 
 \verbatim

  f(x_dot,x,t) = 0.

 \endverbatim
 
 * Whether the problem is an implicit ODE or DAE is determined by the nature
 * of the derivative matrix <tt>d(f)/d(x_dot)</tt>:
 *
 * <ul>
 * <li> ODE: <tt>d(f)/d(x_dot)</tt> is full rank
 * <li> DAE: <tt>d(f)/d(x_dot)</tt> is not full rank
 * </ul>
 *
 * This information is given by the function ??? (ToDo: Add this function!)
 *
 * Here the argument <tt>t</tt> may or may not be accepted by <tt>*this</tt>.
 *
 * \subsection Thyra_ME_unconstrained_optimization_sec Unconstrained optimization
 
 \verbatim

  min g(x,{p(l)})

 \endverbatim

 * where the objective function <tt>g(x,{p(l)})</tt> is some aggregated
 * function built from some subset of the the auxiliary response functions
 * <tt>g(j)(x,{p(l)})</tt>, for <tt>j=1...Ng</tt>.
 *
 * \subsection Thyra_ME_equality_constrained_optimization_sec Equality constrained optimization
 
 \verbatim

  min g(x,{p(l)})

  s.t. f(x,{p(l)}) = 0

 \endverbatim 

 * where the objective function <tt>g(x,{p(l)})</tt> is some aggregated
 * function built from some subset of the the auxiliary response functions
 * <tt>g(j)(x,{p(l)})</tt>, for <tt>j=1...Ng</tt>.
 *
 * \subsection Thyra_ME_general_constrained_optimization_sec Equality constrained optimization
 
 \verbatim

  min g(x,{p(l)})

  s.t. f(x,{p(l)}) = 0
       r(x,{p(l)}) = 0
       hL <= h(x,{p(l)}) <= hU
       xL <= x <= xU
       pL(l) <= p(l) <= pU(l)

 \endverbatim 

 * where the objective function <tt>g(x,{p(l)})</tt> and the auxiliary
 * equality <tt>r(x,{p(l)})</tt> and inequality <tt>h(x,{p(l)})</tt>
 * constraint functions are aggregated functions built from some subset of the
 * the auxiliary response functions <tt>g(j)(x,{p(l)})</tt>, for
 * <tt>j=1...Ng</tt>.  The auxiliary response functions for a particualar
 * model can be interpreted in a wide variety of ways and can be mapped into a
 * number of different optimization problems.
 *
 * \section Thyra_ME_derivatives_sec Function derivatives and sensitivities
 *
 * A model can also optionally support various derivatives of the underlying
 * model functions.  The primary use for these derivatives is in the
 * computation of varyious types of sensitivities.  Specifically, direct and
 * adjoint sensitivities will be considered.
 *
 * To illustrate the issues involved, consider a single auxiliary parameter
 * <tt>p</tt> and a single auxiliary response function <tt>g</tt> of the form
 * <tt>(x,p) => g</tt>.  Assuming that <tt>(x,p) ==> f</tt> defines the state
 * equation <tt>f(x,p)=0</tt> and that <tt>D(f)/D(x)</tt> is full rank, then
 * <tt>f(x,p)=0</tt> defines the implicit function <tt>p ==> x(p)</tt>.  Given
 * this implicit function, the reduced auxiliary function is
 *
 * <tt>g_hat(p) = g(x(p),p)</tt>
 *
 * The reduced derivative <tt>D(g_hat)/D(p)</tt> is given as:
 *
 * <tt>D(g_hat)/D(p) = D(g)/D(x) * D(x)/D(p) + D(g)/D(p)</tt>
 *
 * where <tt>D(x)/D(p) = - [D(f)/D(x)]^{-1} * [D(f)/D(p)]</tt>
 *
 * Restated, the reduced derivative <tt>D(g_hat)/D(p)</tt> is given as:
 *
 * <tt>D(g_hat)/D(p) = - [D(g)/D(x)] * [D(f)/D(x)]^{-1} * [D(f)/D(p)] + D(g)/D(p)</tt>
 *
 * The reduced derivative <tt>D(g_hat)/D(p)</tt> can be computed using the
 * direct or the adjoint approahces.
 *
 * <ul>
 *
 * <li><b>State function derivatives</b>
 *
 *     <ul>
 *     
 *     <li><b>State variable derivatives</b>
 *
 *     <tt>W = alpha*D(f)/D(x_dot) + beta*D(f)/D(x)</tt>
 *
 *     This is derivative operator is a special object that is derived from
 *     the <tt>LinearOpWithSolveBase</tt> interface and therefore supports
 *     linear solves.  Objects of this type are created with the function
 *     <tt>create_W()</tt> before they are computed in <tt>evalModel()</tt>.
 *
 *     <li><b>State variable Taylor coefficients</b>
 *
 *     <tt>x_poly =</tt> \f$\sum_{i=0}^d x_i t^i\f$,
 *
 *     <tt>x_dot_poly =</tt> \f$\sum_{i=0}^d-1 (i+1)x_{i+1} t^i\f$,
 *
 *     <tt>f_poly =</tt> \f$\sum_{i=0}^{d-1} f_i t^i\f$.
 *
 *     If <tt>x_poly</tt> is a given polynomial of degree \f$d\f$ and
 *     <tt>x_dot_poly =</tt> \f$d(\f$<tt>x_poly</tt>\f$)/dt\f$, then <tt>f_poly
 *     = f(x_poly, x_dot_poly, t) +</tt> \f$O(t^{d})\f$ where \f$f_i =
 *     \frac{1}{k!} \frac{d^k}{dt^k} f(dx/dt ,x(t),t), i=0,\ldots,d-1\f$ are
 *     the Taylor series coefficient of \f$f\f$.  The polynomials
 *     <tt>x_poly</tt>, <tt>x_dot_poly</tt>, and <tt>f_poly</tt> are
 *     represented by <tt>Teuchos::Polynomial</tt> objects where each
 *     coefficient is a <tt>Thyra::VectorBase</tt> object.  The Taylor series
 *     coefficients of \f$f\f$ can easily be computed using automatic
 *     differentiation, but this is generally the only means of doing so.
 *     
 *     <li><b>Auxiliary parameter derivatives</b>
 *
 *     <tt>DfDp(l) = D(f)/D(p(l))</tt>, for <tt>l=1...Np</tt>
 *
 *     These are derivative objects that represent the derivative of the state
 *     function <tt>f</tt> with respect to the auxiliary parameters
 *     <tt>p(l)</tt>.  This object can either be represented and manipulated
 *     in a very abstract fashion through the <tt>LinearOpBase</tt> interface
 *     interface or as a multi-vector object.
 *
 *     </ul>
 *
 * <li><b>Auxiliary response function derivatives</b>
 *
 * </ul>
 *
 * \section Thyra_ME_dev_notes_sec Notes to subclass devleopers
 *
 * This interface is setup so that the default problem type a set of nonlinear
 * equaitons <tt>f(x)=0</tt>.  Therefore, the only pure virtual response
 * functions that need to be overrriden to create a concrete subclass are
 * <tt>get_x_space()</tt>, <tt>get_f_space()</tt>, <tt>createInArgs()</tt>,
 * <tt>createOutArgs()</tt>, and <tt>evalModel()</tt>.  All of the other
 * virtual functions have default implementations that are appropriate for
 * this type of problem.  Supporting other problem types involves overridding
 * various sets of virtual functions (see able).
 *
 * ToDo: Finish Documentation!
 */
template<class Scalar>
class ModelEvaluator : public ModelEvaluatorBase {
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Basic inforamtion */
  //@{

	/** \brief Return the number of sets of auxiliary parameters.
   *
	 * If this function returns 0, then there are no auxiliary parameters.
	 *
	 * The default implementation returns 0.
	 */
	virtual int Np() const;

	/** \brief Return the number of sets of auxiliary response functions.
   *
	 * If this function returns 0, then there are no auxiliary response
	 * functions.
	 *
	 * The default implementation returns 0.
	 */
	virtual int Ng() const;

  //@}

  /** \name Vector spaces */
  //@{

  /** \brief . */
  virtual Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_x_space() const = 0;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_f_space() const = 0;

	/** \brief Return the vector space for the auxiliary parameters
	 * <tt>p(l)</tt>.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->Np() > 0</tt>
	 * <li><tt>1 <= l <= this->Np()</tt>
	 * </ul>
	 *
	 * <b>Postconditions:</b><ul>
	 * <li> <tt>return.get()!=NULL</tt>
	 * </ul>
	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Np()==0</tt>.
	 */
	virtual Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_p_space(int l) const;

	/** \brief Return the vector space for the auxiliary response functions
	 * <tt>g(j)</tt>.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->Ng() > 0</tt>
	 * <li><tt>1 <= j <= this->Ng()</tt>
	 * </ul>
	 *
	 * <b>Postconditions:</b><ul>
	 * <li> <tt>return.get()!=NULL</tt>
	 * </ul>
	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Ng()==0</tt>.
	 */
	virtual Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_g_space(int j) const;

  //@}

  /** \name Initial guesses for variables/parameters */
  //@{

  /** \brief Return an optional initial guess for x.
   *
   * If an initial guess is not supported then <tt>return.get()==NULL</tt>.
   *
	 * The default implementation returns <tt>return.get()==NULL</tt>.
   */
  virtual Teuchos::RefCountPtr<const VectorBase<Scalar> > get_x_init() const;

	/** \brief Return an initial guess for <tt>p(l)</tt>
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->Np() > 0</tt>
	 * <li><tt>1 <= l && l <= this->Np()</tt>
	 * </ul>
	 *
	 * <b>Postconditions:</b><ul>
	 * <li> If <tt>return.get()!=NULL</tt> then <tt>*return</tt> gives the
   *      initial guess for <tt>p(l)</tt>.
	 * </ul>
	 *
   * If an initial guess is not supported then <tt>return.get()==NULL</tt>.
   *
	 * The default implementation returns <tt>return.get()==NULL</tt>.
	 */
	virtual Teuchos::RefCountPtr<const VectorBase<Scalar> > get_p_init(int l) const;

  /** \brief Return an optional initial guess for t.
   *
   * If an initial guess is not supported then <tt>return==0.0</tt>.
   *
	 * The default implementation returns <tt>return==0.0</tt>.
   */
  virtual ScalarMag get_t_init() const;

  //@}

  /** \name Bounds for variables/parameters */
  //@{

  /** \brief Return lower bounds for <tt>x</tt> if supported.
   *
   * If these bounds are not supported then <tt>return.get()==NULL</tt>.
   */
  virtual Teuchos::RefCountPtr<const VectorBase<Scalar> > get_x_lower_bounds() const;

  /** \brief Return upper bounds for <tt>x</tt> if supported.
   *
   * If these bounds are not supported then <tt>return.get()==NULL</tt>.
   */
  virtual Teuchos::RefCountPtr<const VectorBase<Scalar> > get_x_upper_bounds() const;

	/** \brief Return lower bounds for <tt>p(l)</tt> if supported.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->Np() > 0</tt>
	 * <li><tt>1 <= l <= this->Np()</tt>
	 * </ul>
	 *
	 * <b>Postconditions:</b><ul>
	 * <li> If <tt>return.get()!=NULL</tt> then <tt>*return</tt> gives the bounds.
	 * </ul>
	 */
	virtual Teuchos::RefCountPtr<const VectorBase<Scalar> > get_p_lower_bounds(int l) const;

	/** \brief Return upper bounds for <tt>p(l)</tt> if supported.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->Np() > 0</tt>
	 * <li><tt>1 <= l <= this->Np()</tt>
	 * </ul>
	 *
	 * <b>Postconditions:</b><ul>
	 * <li> If <tt>return.get()!=NULL</tt> then <tt>*return</tt> gives the bounds.
	 * </ul>
	 */
	virtual Teuchos::RefCountPtr<const VectorBase<Scalar> > get_p_upper_bounds(int l) const;

  /** \brief Return lower bound for <tt>t</tt> if supported.
   */
  virtual ScalarMag get_t_lower_bound() const;

  /** \brief Return upper bound for <tt>t</tt> if supported.
   */
  virtual ScalarMag get_t_upper_bound() const;

  //@}

  /** \name Factory functions for creating derivative objects */
  //@{

  /** \brief If supported, create a <tt>LinearOpWithSolveBase</tt> object for
   * <tt>W</tt> to be evaluated.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->createOutArgs().supports(OUT_ARG_W)==true</tt>
	 * </ul>
   *
   * The default implementation returns <tt>return.get()==NULL</tt>
   * (i.e. implicit solvers are not supported by default).
   */
  virtual Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > create_W() const;

  /** \brief If supported, create a linear operator derivative object for
   * <tt>D(f)/D(p(l))</tt>.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->Np() > 0</tt>
	 * <li><tt>1 <= l <= this->Np()</tt>
   * <li><tt>outArgs.supports_DfDp(l).supports(DERIV_LINEAR_OP)==true</tt>,
   *     where <tt>outArgs = this->createOutArgs()</tt>
	 * </ul>
   *
   * The default implementation returns <tt>return.get()==NULL</tt>
   * (i.e. no sensitivities are supported by default).
   */
  virtual Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DfDp_op(int l) const;

  /** \brief If supported, create a multi-vector derivative object for
   * <tt>D(f)/D(p(l))</tt> or <tt>D(f)/D(p(l))^T</tt>.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->Np() > 0</tt>
	 * <li><tt>1 <= l <= this->Np()</tt>
   * <li><tt>outArgs.supports_DfDp(l).supports(orientation)==true</tt>,
   *     where <tt>outArgs = this->createOutArgs()</tt>
	 * </ul>
   *
   * The default implementation returns <tt>return.getMultiVector().get()==NULL</tt>
   * (i.e. no sensitivities are supported by default).
   */
  virtual DerivativeMultiVector<Scalar> create_DfDp_mv(int l, EDerivativeMultiVectorOrientation orientation) const;

  // ToDo: Add functions for creating D(g(j))/D(x_dot) if needed!

  /** \brief If supported, create a linear operator derivative object for
   * <tt>D(g(j))/D(x)</tt>.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->Ng() > 0</tt>
	 * <li><tt>1 <= j <= this->Ng()</tt>
   * <li><tt>outArgs.supports_DgDx(j).supports(DERIV_LINEAR_OP)==true</tt>,
   *     where <tt>outArgs = this->createOutArgs()</tt>
	 * </ul>
   *
   * The default implementation returns <tt>return.get()==NULL</tt>
   * (i.e. no sensitivities are supported by default).
   */
  virtual Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DgDx_op(int j) const;

  /** \brief If supported, create a multi-vector derivative object for
   * <tt>D(g(j))/D(x)</tt> or <tt>D(g(j))/D(x)^T</tt>.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->Ng() > 0</tt>
	 * <li><tt>1 <= j <= this->Ng()</tt>
   * <li><tt>outArgs.supports_DgDx(j).supports(orientation)==true</tt>,
   *     where <tt>outArgs = this->createOutArgs()</tt>
	 * </ul>
   *
   * The default implementation returns <tt>return.getMultiVector().get()==NULL</tt>
   * (i.e. no sensitivities are supported by default).
   */
  virtual DerivativeMultiVector<Scalar> create_DgDx_mv(int j, EDerivativeMultiVectorOrientation orientation) const;

  /** \brief If supported, create a linear operator derivative object for
   * <tt>D(g(j))/D(p(l))</tt>.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->Ng() > 0</tt>
	 * <li><tt>this->Np() > 0</tt>
	 * <li><tt>1 <= j <= this->Ng()</tt>
	 * <li><tt>1 <= l <= this->Np()</tt>
   * <li><tt>outArgs.supports_DgDp(j,l).supports(DERIV_LINEAR_OP)==true</tt>,
   *     where <tt>outArgs = this->createOutArgs()</tt>
	 * </ul>
   *
   * The default implementation returns <tt>return.get()==NULL</tt>
   * (i.e. no sensitivities are supported by default).
   */
  virtual Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DgDp_op( int j, int l ) const;

  /** \brief If supported, create a multi-vector derivative object for
   * <tt>D(g(j))/D(p(l))</tt> or <tt>D(g(j))/D(p(l))^T</tt>.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->Ng() > 0</tt>
	 * <li><tt>this->Np() > 0</tt>
	 * <li><tt>1 <= j <= this->Ng()</tt>
	 * <li><tt>1 <= l <= this->Np()</tt>
   * <li><tt>outArgs.supports_DgDp(j,l).supports(orientation)==true</tt>,
   *     where <tt>outArgs = this->createOutArgs()</tt>
	 * </ul>
   *
   * The default implementation returns <tt>return.getMultiVector().get()==NULL</tt>
   * (i.e. no sensitivities are supported by default).
   */
  virtual DerivativeMultiVector<Scalar> create_DgDp_mv( int j, int l, EDerivativeMultiVectorOrientation orientation ) const;
  
  //@}

  /** \name Computational functions */
  //@{

  /** \brief . */
  virtual ModelEvaluatorBase::InArgs<Scalar> createInArgs() const = 0;

  /** \brief . */
  virtual ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const = 0;

  /** \brief . */
  virtual void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar>       &inArgs
    ,const ModelEvaluatorBase::OutArgs<Scalar>     &outArgs
    ) const = 0;

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

/** \brief . */
template<class Scalar>
inline 
void eval_f_poly(
  const ModelEvaluator<Scalar>                                    &model
  ,const Teuchos::Polynomial< VectorBase<Scalar> >                &x_poly
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag   &t
  ,Teuchos::Polynomial< VectorBase<Scalar> >                      *f_poly
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();

  inArgs.set_x_poly(Teuchos::rcp(&x_poly,false));
  if(inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(t);

  outArgs.set_f_poly(Teuchos::rcp(f_poly,false));

  model.evalModel(inArgs,outArgs);

}

/** \brief . */
template<class Scalar>
inline 
void eval_f_poly(
  const ModelEvaluator<Scalar>                                    &model
  ,const Teuchos::Polynomial< VectorBase<Scalar> >                &x_dot_poly
  ,const VectorBase<Scalar>                                       &x_poly
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag   &t
  ,Teuchos::Polynomial< VectorBase<Scalar> >                      *f_poly
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();

  inArgs.set_x_dot_poly(Teuchos::rcp(&x_dot_poly,false));
  inArgs.set_x_poly(Teuchos::rcp(&x_poly,false));
  if(inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(t);

  outArgs.set_f_poly(Teuchos::rcp(f_poly,false));

  model.evalModel(inArgs,outArgs);

}

} // namespace Thyra

// //////////////////////////////////
// Inline Defintions

//
// Thyra_MEB_helper_functions_grp
//

inline
std::string Thyra::toString(ModelEvaluatorBase::EInArgsMembers arg)
{
  switch(arg) {
    case ModelEvaluatorBase::IN_ARG_x_dot:
      return "OUT_ARG_x_dot";
    case ModelEvaluatorBase::IN_ARG_x:
      return "OUT_ARG_x";
    case ModelEvaluatorBase::IN_ARG_t:
      return "OUT_ARG_t";
    case ModelEvaluatorBase::IN_ARG_alpha:
      return "OUT_ARG_alpha";
    case ModelEvaluatorBase::IN_ARG_beta:
      return "OUT_ARG_beta";
    default:
      TEST_FOR_EXCEPT(true);
  }
  return ""; // Will never be executed!
}

inline
std::string Thyra::toString(ModelEvaluatorBase::EOutArgsMembers arg)
{
  switch(arg) {
    case ModelEvaluatorBase::OUT_ARG_f:
      return "OUT_ARG_f";
    case ModelEvaluatorBase::OUT_ARG_W:
      return "OUT_ARG_W";
    default:
      TEST_FOR_EXCEPT(true);
  }
  return ""; // Will never be executed!
}

// //////////////////////////////////
// Definitions

namespace Thyra {

//
// ModelEvaluatorBase::InArgs
//

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
int ModelEvaluatorBase::InArgs<Scalar>::Np() const
{ return p_.size(); }

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
void ModelEvaluatorBase::InArgs<Scalar>::set_x_dot_poly( const Teuchos::RefCountPtr<const Teuchos::Polynomial< VectorBase<Scalar> > > &x_dot_poly )
{ assert_supports(IN_ARG_x_dot_poly); x_dot_poly_ = x_dot_poly; }

template<class Scalar>
Teuchos::RefCountPtr<const Teuchos::Polynomial< VectorBase<Scalar> > >
ModelEvaluatorBase::InArgs<Scalar>::get_x_dot_poly() const
{ assert_supports(IN_ARG_x_dot_poly); return x_dot_poly_; }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x_poly( const Teuchos::RefCountPtr<const Teuchos::Polynomial< VectorBase<Scalar> > > &x_poly )
{ assert_supports(IN_ARG_x_poly); x_poly_ = x_poly; }

template<class Scalar>
Teuchos::RefCountPtr<const Teuchos::Polynomial< VectorBase<Scalar> > >
ModelEvaluatorBase::InArgs<Scalar>::get_x_poly() const
{ assert_supports(IN_ARG_x_poly); return x_poly_; }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_p( int l, const Teuchos::RefCountPtr<const VectorBase<Scalar> > &p_l )
{ assert_l(l); p_[l-1] = p_l; }

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_p(int l) const
{ assert_l(l); return p_[l-1]; }

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
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"*this = \'"<<modelEvalDescription_<<"\': Error, arg="<<arg<<" is invalid!"
    );
  return supports_[arg];
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setModelEvalDescription( const std::string &modelEvalDescription )
{ modelEvalDescription_ = modelEvalDescription; }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_set_Np(int Np)
{
  p_.resize(Np);
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setSupports( EInArgsMembers arg, bool supports )
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"*this = \'"<<modelEvalDescription_<<"\': Error, arg="<<arg<<" is invalid!");
  supports_[arg] = supports;
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::assert_supports(EInArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"Thyra::ModelEvaluatorBase::InArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(arg): "
    "*this = \'"<<modelEvalDescription_<<"\': Error, "
    "The argument arg = " << toString(arg) << " is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::assert_l(int l) const
{
  TEST_FOR_EXCEPTION(
    !( 1 <= l && l <= Np() ), std::logic_error
    ,"Thyra::ModelEvaluatorBase::InArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_l(l): "
    " *this = \'"<<modelEvalDescription_<<"\': Error, "
    "The parameter l = " << l << " is not in the range [1,"<<Np()<<"]!"
    );
}

//
// ModelEvaluatorBase::OutArgs
//

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::OutArgs()
{ std::fill_n(&supports_[0],NUM_E_OUT_ARGS_MEMBERS,false); }

template<class Scalar>
int ModelEvaluatorBase::OutArgs<Scalar>::Np() const
{ return DfDp_.size(); }

template<class Scalar>
int ModelEvaluatorBase::OutArgs<Scalar>::Ng() const
{ return g_.size(); }

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(EOutArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"*this = \'"<<modelEvalDescription_<<"\': Error, arg="<<arg<<" is invalid!"
    );
  return supports_[arg];
}

template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(EOutArgsDfDp arg, int l) const
{
  assert_l(l);
  return supports_DfDp_[l-1];
}

template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(EOutArgsDgDx arg, int j) const
{
  assert_j(j);
  return supports_DgDx_[j-1];
}

template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(EOutArgsDgDp arg, int j, int l) const
{
  assert_j(j);
  assert_l(l);
  return supports_DgDp_[ (j-1)*Np() + (l-1) ];
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_f( const Teuchos::RefCountPtr<VectorBase<Scalar> > &f )
{
  assert_supports(OUT_ARG_f);
  f_ = f;
}

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_f() const
{
  assert_supports(OUT_ARG_f);
  return f_;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_g( int j, const Teuchos::RefCountPtr<VectorBase<Scalar> > &g_j )
{
  assert_j(j);
  g_[j-1] = g_j;
}

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_g(int j) const
{ 
  assert_j(j);
  return g_[j-1];
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_W( const Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > &W )
{
  assert_supports(OUT_ARG_W);
  W_ = W;
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_W() const
{
  assert_supports(OUT_ARG_W);
  return W_;
}

template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_W_properties() const
{
  return W_properties_;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DfDp( int l, const Derivative<Scalar> &DfDp_l )
{
  assert_supports(OUT_ARG_DfDp,l);
  DfDp_[l-1] = DfDp_l;
}

template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::get_DfDp(int l) const
{
  assert_supports(OUT_ARG_DfDp,l);
  return DfDp_[l-1];
}

template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DfDp_properties(int l) const
{
  assert_supports(OUT_ARG_DfDp,l);
  return DfDp_properties_[l-1];
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DgDx( int j, const Derivative<Scalar> &DgDx_j )
{
  assert_supports(OUT_ARG_DgDx,j);
  DgDx_[j-1] = DgDx_j;
}

template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx(int j) const
{
  assert_supports(OUT_ARG_DgDx,j);
  return DgDx_[j-1];
}

template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx,j);
  return DgDx_properties_[j-1];
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DgDp( int j, int l, const Derivative<Scalar> &DgDp_j_l )
{
  assert_supports(OUT_ARG_DgDp,j,l);
  DgDp_[ (j-1)*Np() + (l-1) ] = DgDp_j_l;
}

template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDp(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp,j,l);
  return DgDp_[ (j-1)*Np() + (l-1) ];
}

template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDp_properties(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp,j,l);
  return DgDp_properties_[ (j-1)*Np() + (l-1) ];
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_f_poly( const Teuchos::RefCountPtr<Teuchos::Polynomial< VectorBase<Scalar> > > &f_poly )
{ f_poly_ = f_poly; }

template<class Scalar>
Teuchos::RefCountPtr<Teuchos::Polynomial< VectorBase<Scalar> > >
ModelEvaluatorBase::OutArgs<Scalar>::get_f_poly() const
{ return f_poly_; }

// protected

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setModelEvalDescription( const std::string &modelEvalDescription )
{ modelEvalDescription_ = modelEvalDescription; }

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_Np_Ng(int Np, int Ng)
{
  if(Np) {
    supports_DfDp_.resize(Np);
    DfDp_.resize(Np);                 std::fill_n(DfDp_.begin(),Np,Derivative<Scalar>());
    DfDp_properties_.resize(Np);      std::fill_n(DfDp_properties_.begin(),Np,DerivativeProperties());
  }
  if(Ng) {
    g_.resize(Ng);                    std::fill_n(g_.begin(),Ng,Teuchos::null);
    supports_DgDx_.resize(Ng);
    DgDx_.resize(Ng);                 std::fill_n(DgDx_.begin(),Ng,Derivative<Scalar>());
    DgDx_properties_.resize(Ng);      std::fill_n(DgDx_properties_.begin(),Ng,DerivativeProperties());
  }
  if(Np && Ng) {
    const int NpNg = Np*Ng;
    supports_DgDp_.resize(NpNg);
    DgDp_.resize(NpNg);                 std::fill_n(DgDp_.begin(),NpNg,Derivative<Scalar>());
    DgDp_properties_.resize(NpNg);      std::fill_n(DgDp_properties_.begin(),NpNg,DerivativeProperties());
  }
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports( EOutArgsMembers arg, bool supports )
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"*this = \'"<<modelEvalDescription_<<"\': Error, arg="<<arg<<" is invalid!"
    );
  supports_[arg] = supports;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports( EOutArgsDfDp arg, int l, const DerivativeSupport& supports )
{
  assert_l(l);
  supports_DfDp_[l-1] = supports;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports( EOutArgsDgDx arg, int j, const DerivativeSupport& supports )
{
  assert_j(j);
  supports_DgDx_[j-1] = supports;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports( EOutArgsDgDp arg, int j, int l, const DerivativeSupport& supports )
{
  assert_j(j);
  assert_l(l);
  supports_DgDp_[ (j-1)*Np() + (l-1) ] = supports;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_W_properties( const DerivativeProperties &properties )
{
  W_properties_ = properties;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DfDp_properties( int l, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DfDp,l);
  DfDp_properties_[l-1] = properties;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DgDx_properties( int j, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDx,j);
  DgDx_properties_[j-1] = properties;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DgDp_properties( int j, int l, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDp,j,l);
  DgDp_properties_[ (j-1)*Np() + (l-1) ] = properties;
}

// private

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(EOutArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(arg): "
    "*this = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument arg = " << toString(arg) << " is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(EOutArgsDfDp arg, int l) const
{
  assert_l(l);
  TEST_FOR_EXCEPTION(
    supports_DfDp_[l-1].none(), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(OUT_ARG_DfDp,l): "
    "*this = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DfDp(l) with index l = " << l << " is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(EOutArgsDgDx arg, int j) const
{
  assert_j(j);
  TEST_FOR_EXCEPTION(
    supports_DgDx_[j-1].none(), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(OUT_ARG_DgDx,j): "
    "*this = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx(j) with index j = " << j << " is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(EOutArgsDgDp arg, int j, int l) const
{
  assert_j(j);
  TEST_FOR_EXCEPTION(
    supports_DgDx_[ (j-1)*Np() + (l-1) ].none(), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(OUT_ARG_DgDp,j,l): "
    "*this = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDp(j,l) with indexes j = " << j << " and l = " << l << " is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_l(int l) const
{
  TEST_FOR_EXCEPTION(
    !( 1 <= l && l <= Np() ), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_l(l): "
    "*this = \'"<<modelEvalDescription_<<"\': Error, "
    "The parameter subvector p(l) index l = " << l << " is not in the range [1,"<<Np()<<"]!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_j(int j) const
{
  TEST_FOR_EXCEPTION(
    !( 1 <= j && j <= Ng() ), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_j(j): "
    "*this = \'"<<modelEvalDescription_<<"\': Error, "
    "The auxiliary function g(j) index j = " << j << " is not in the range [1,"<<Ng()<<"]!"
    );
}

//
// ModelEvaluatorBase::InArgsSetup
//

template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::setModelEvalDescription( const std::string &modelEvalDescription )
{ this->_setModelEvalDescription(modelEvalDescription); }

template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::set_Np(int Np)
{ this->_set_Np(Np); }

template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::setSupports( EInArgsMembers arg, bool supports )
{ this->_setSupports(arg,supports); }

//
// ModelEvaluatorBase::OutArgsSetup
//

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setModelEvalDescription( const std::string &modelEvalDescription )
{ this->_setModelEvalDescription(modelEvalDescription); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_Np_Ng(int Np, int Ng)
{ this->_set_Np_Ng(Np,Ng); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports( EOutArgsMembers arg, bool supports )
{ this->_setSupports(arg,supports); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports( EOutArgsDfDp arg, int l, const DerivativeSupport& supports )
{ this->_setSupports(arg,l,supports); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports( EOutArgsDgDx arg, int j, const DerivativeSupport& supports )
{ this->_setSupports(arg,j,supports); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports( EOutArgsDgDp arg, int j, int l, const DerivativeSupport& supports )
{ this->_setSupports(arg,j,l,supports); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_W_properties( const DerivativeProperties &properties )
{ this->_set_W_properties(properties); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DfDp_properties( int l, const DerivativeProperties &properties )
{
  this->_set_DfDp_properties(l,properties);
}

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DgDx_properties( int j, const DerivativeProperties &properties )
{
  this->_set_DgDx_properties(j,properties);
}

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DgDp_properties( int j, int l, const DerivativeProperties &properties )
{
  this->_set_DgDp_properties(j,l,properties);
}

//
// ModelEvaluator
//

// Basic inforamtion

template<class Scalar>
int ModelEvaluator<Scalar>::Np() const
{ return 0; }

template<class Scalar>
int ModelEvaluator<Scalar>::Ng() const
{ return 0; }

// Vector spaces

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
ModelEvaluator<Scalar>::get_p_space(int l) const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"ModelEvaluator<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::get_p_space(l): "
    "Error, this function was not overridden in *this = \'"<<this->description()<<"\'!"
		);
	return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
ModelEvaluator<Scalar>::get_g_space(int j) const
{
	TEST_FOR_EXCEPTION(
		true,std::logic_error
		,"ModelEvaluator<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::get_g_space(j): "
    " Error, this function was not overridden in \'"
    <<this->description()<<"\'!"
		);
	return Teuchos::null;
}

// Initial guesses for variables/parameters

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
ModelEvaluator<Scalar>::get_x_init() const
{ return Teuchos::null; }

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
ModelEvaluator<Scalar>::get_p_init(int l) const
{ return Teuchos::null; }

template<class Scalar>
typename ModelEvaluator<Scalar>::ScalarMag
ModelEvaluator<Scalar>::get_t_init() const
{ return 0.0; }

// Bounds for variables/parameters

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
ModelEvaluator<Scalar>::get_x_lower_bounds() const
{ return Teuchos::null; }

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
ModelEvaluator<Scalar>::get_x_upper_bounds() const
{ return Teuchos::null; }

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
ModelEvaluator<Scalar>::get_p_lower_bounds(int l) const
{ return Teuchos::null; }

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
ModelEvaluator<Scalar>::get_p_upper_bounds(int l) const
{ return Teuchos::null; }

template<class Scalar>
typename ModelEvaluator<Scalar>::ScalarMag
ModelEvaluator<Scalar>::get_t_lower_bound() const
{ return 0.0; }

template<class Scalar>
typename ModelEvaluator<Scalar>::ScalarMag
ModelEvaluator<Scalar>::get_t_upper_bound() const
{ return 0.0; }

// Factory functions for creating derivative objects

template<class Scalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
ModelEvaluator<Scalar>::create_W() const
{ return Teuchos::null; }

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
ModelEvaluator<Scalar>::create_DfDp_op(int l) const
{ return Teuchos::null; }

template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
ModelEvaluator<Scalar>::create_DfDp_mv(int l, EDerivativeMultiVectorOrientation orientation) const
{ return DerivativeMultiVector<Scalar>(); }

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
ModelEvaluator<Scalar>::create_DgDx_op(int j) const
{ return Teuchos::null; }

template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
ModelEvaluator<Scalar>::create_DgDx_mv(int j, EDerivativeMultiVectorOrientation orientation) const
{ return DerivativeMultiVector<Scalar>(); }

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
ModelEvaluator<Scalar>::create_DgDp_op( int j, int l ) const
{ return Teuchos::null; }

template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
ModelEvaluator<Scalar>::create_DgDp_mv( int j, int l, EDerivativeMultiVectorOrientation orientation ) const
{ return DerivativeMultiVector<Scalar>(); }

} // namespace Thyra


#endif // THYRA_MODEL_EVALUATOR_HPP
