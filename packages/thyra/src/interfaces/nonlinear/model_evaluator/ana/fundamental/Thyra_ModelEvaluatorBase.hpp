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

#ifndef THYRA_MODEL_EVALUATOR_BASE_HPP
#define THYRA_MODEL_EVALUATOR_BASE_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_Polynomial.hpp"
#include "Teuchos_Array.hpp"
#include "Thyra_PolynomialVectorTraits.hpp"

namespace Thyra {

/** \brief Base subclass for <tt>ModelEvaluator</tt> that defines some basic
 * types.
 *
 * ToDo: Finish documentation!
 */
class ModelEvaluatorBase
  : virtual public Teuchos::Describable
  , virtual public Teuchos::VerboseObject<ModelEvaluatorBase>
{
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
  class InArgs : public Teuchos::Describable {
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
    /** \brief Set <tt>p(l)</tt> where <tt>0 <= l && l < this->Np()</tt>.  */
    void set_p( int l, const Teuchos::RefCountPtr<const VectorBase<Scalar> > &p_l );
    /** \brief Get <tt>p(l)</tt> where <tt>0 <= l && l < this->Np()</tt>.  */
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
    /** \brief Set non-null arguments (does not overwrite non-NULLs with NULLs) .  */
    void setArgs( const InArgs<Scalar>& inArgs, bool ignoreUnsupported = false );
    /** \brief . */
    std::string description() const;
    /** \brief . */
    void describe( Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel ) const;
  protected:
    /** \brief . */
    void _setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void _set_Np(int Np);
    /** \brief . */
    void _setSupports( EInArgsMembers arg, bool supports );
    /** \brief . */
    void _setSupports( const InArgs<Scalar>& inArgs );
    /** \brief . */
    void _setUnsupportsAndRelated( EInArgsMembers arg );
  private:
    // types
    typedef Teuchos::Array<Teuchos::RefCountPtr<const VectorBase<Scalar> > > p_t;
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
  enum EDerivativeLinearOp {
    DERIV_LINEAR_OP ///< .
  };

  /** \brief Determines the forms of a general derivative that are
   * supported. */
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
      ) : mv_(mv.assert_not_null()), orientation_(orientation) {}
    /** \brief . */
    void changeOrientation( const EDerivativeMultiVectorOrientation orientation )
      { orientation_ = orientation; };
    /** \brief . */
    const DerivativeMultiVector<Scalar>& assert_not_null() const
      { mv_.assert_not_null(); return *this; }
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
      : lo_(lo.assert_not_null()) {}
    /** \brief . */
    Derivative( const DerivativeMultiVector<Scalar> &dmv )
      : dmv_(dmv) {}
    /** \brief . */
    bool isEmpty() const
      { return ( lo_.get()==NULL && dmv_.getMultiVector().get()==NULL ); }
    /** \brief . */
    const Derivative<Scalar>& assert_not_null() const
      { dmv_.assert_not_null(); lo_.assert_not_null(); return *this; }
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
    ,OUT_ARG_W_op   ///< .
    ,OUT_ARG_f_poly ///< .
  };
  /** \brief .  */
  static const int NUM_E_OUT_ARGS_MEMBERS=4;

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
  class OutArgs : public Teuchos::Describable {
  public:
    /** \brief .  */
    OutArgs();
    /** \brief .  */
    int Np() const;
    /** \brief .  */
    int Ng() const;
    /** \brief .  */
    bool supports(EOutArgsMembers arg) const;
    /** \brief <tt>0 <= l && l < Np()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDfDp arg, int l) const;
    /** \brief <tt>0 <= j && j < Ng()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDx arg, int j) const;
    /** \brief <tt>0 <= j && j < Ng()</tt> and <tt>0 <= l && l < Np()</tt>.  */
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
    /** \brief .  */
    void set_W_op( const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &W_op );
    /** \brief .  */
    Teuchos::RefCountPtr<LinearOpBase<Scalar> > get_W_op() const;
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
    /** \brief .  */
    void setArgs( const OutArgs<Scalar>& outArgs, bool ignoreUnsupported = false );
    /** \brief . */
    void setFailed() const;
    /** \brief . */
    bool isFailed() const;
    /** \brief . */
    std::string description() const;
    /** \brief . */
    void describe( Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel ) const;
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
    /** \brief . */
    void _setSupports( const OutArgs<Scalar>& outArgs );
    /** \brief . */
    void _setUnsupportsAndRelated( EInArgsMembers arg );
    /** \brief . */
    void _setUnsupportsAndRelated( EOutArgsMembers arg );
  private:
    // types
    typedef Teuchos::Array<Teuchos::RefCountPtr<VectorBase<Scalar> > >     g_t;
    typedef Teuchos::Array<Derivative<Scalar> >                            deriv_t;
    typedef Teuchos::Array<DerivativeProperties>                           deriv_properties_t;
    typedef Teuchos::Array<DerivativeSupport>                              supports_t;
    // data
    std::string                                           modelEvalDescription_;
    bool                                                  supports_[NUM_E_OUT_ARGS_MEMBERS];
    supports_t                                            supports_DfDp_;   // Np
    supports_t                                            supports_DgDx_;   // Ng
    supports_t                                            supports_DgDp_;   // Ng x Np
    Teuchos::RefCountPtr<VectorBase<Scalar> >             f_;
    g_t                                                   g_;               // Ng
    Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >  W_;
    Teuchos::RefCountPtr<LinearOpBase<Scalar> >           W_op_;
    DerivativeProperties                                  W_properties_;
    deriv_t                                               DfDp_;            // Np
    deriv_properties_t                                    DfDp_properties_; // Np
    deriv_t                                               DgDx_;            // Ng
    deriv_properties_t                                    DgDx_properties_; // Ng
    deriv_t                                               DgDp_;            // Ng x Np
    deriv_properties_t                                    DgDp_properties_; // Ng x Np
    Teuchos::RefCountPtr<Teuchos::Polynomial< VectorBase<Scalar> > > f_poly_;
    mutable bool                                          isFailed_;
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
    InArgsSetup();
    /** \brief . */
    InArgsSetup( const InArgs<Scalar>& );
    /** \brief . */
    void setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void set_Np(int Np);
    /** \brief . */
    void setSupports( EInArgsMembers arg, bool supports = true );
    /** \brief . */
    void setSupports( const InArgs<Scalar>& inArgs );
    /** \brief . */
    void setUnsupportsAndRelated( EInArgsMembers arg );
  };

  /** \brief . */
  template<class Scalar>
  class OutArgsSetup : public OutArgs<Scalar> {
  public:
    /** \brief . */
    OutArgsSetup();
    /** \brief . */
    OutArgsSetup( const OutArgs<Scalar>& );
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
    /** \brief . */
    void setSupports( const OutArgs<Scalar>& outArgs );
   /** \brief . */
    void setUnsupportsAndRelated( EInArgsMembers arg );
    /** \brief . */
    void setUnsupportsAndRelated( EOutArgsMembers arg );
   };

  //@}

};

/** \defgroup Thyra_MEB_helper_functions_grp Helper functions for Thyra::ModelEvaluatorBase.
 *
 */
//@{

/** \relates ModelEvaluatorBase */
std::string toString(ModelEvaluatorBase::EInArgsMembers);

/** \relates ModelEvaluatorBase */
std::string toString(ModelEvaluatorBase::EOutArgsMembers);

/** \relates ModelEvaluatorBase */
std::string toString(ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation);

//@}

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
      return "IN_ARG_x_dot";
    case ModelEvaluatorBase::IN_ARG_x:
      return "IN_ARG_x";
    case ModelEvaluatorBase::IN_ARG_t:
      return "IN_ARG_t";
    case ModelEvaluatorBase::IN_ARG_alpha:
      return "IN_ARG_alpha";
    case ModelEvaluatorBase::IN_ARG_beta:
      return "IN_ARG_beta";
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

inline
std::string Thyra::toString(ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation)
{
  switch(orientation) {
    case ModelEvaluatorBase::DERIV_MV_BY_COL:
      return "DERIV_MV_BY_COL";
    case ModelEvaluatorBase::DERIV_TRANS_MV_BY_ROW:
      return "DERIV_TRANS_MV_BY_ROW";
    default:
      TEST_FOR_EXCEPT(true);
  }
  return ""; // Should never execute this!
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
{ assert_l(l); p_[l] = p_l; }

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_p(int l) const
{ assert_l(l); return p_[l]; }

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
    ,"model = \'"<<modelEvalDescription_<<"\': Error, arg="<<toString(arg)<<" is invalid!"
    );
  return supports_[arg];
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::setArgs( const InArgs<Scalar>& inArgs, bool ignoreUnsupported )
{
  if( inArgs.supports(IN_ARG_x_dot) && inArgs.get_x_dot().get() ) {
    if(supports(IN_ARG_x_dot) || !ignoreUnsupported)
      set_x_dot(inArgs.get_x_dot());
  }
  if( inArgs.supports(IN_ARG_x) && inArgs.get_x().get() ) {
    if(supports(IN_ARG_x) || !ignoreUnsupported)
      set_x(inArgs.get_x());
  }
  if( inArgs.supports(IN_ARG_x_dot_poly) && inArgs.get_x_dot_poly().get() ) {
    if(supports(IN_ARG_x_dot_poly) || !ignoreUnsupported)
      set_x_dot_poly(inArgs.get_x_dot_poly());
  }
  if( inArgs.supports(IN_ARG_x_poly) && inArgs.get_x_poly().get() ) {
    if(supports(IN_ARG_x_poly) || !ignoreUnsupported)
      set_x_poly(inArgs.get_x_poly());
  }
  const int min_Np = TEUCHOS_MIN(this->Np(),inArgs.Np());
  for( int l = 0; l < min_Np; ++l ) {
    if(inArgs.get_p(l).get())
      set_p(l,inArgs.get_p(l));
  }
  if( inArgs.supports(IN_ARG_t) ) {
    if(supports(IN_ARG_t) || !ignoreUnsupported)
      set_t(inArgs.get_t());
  }
  if( inArgs.supports(IN_ARG_alpha) ) {
    if(supports(IN_ARG_alpha) || !ignoreUnsupported)
      set_alpha(inArgs.get_alpha());
  }
  if( inArgs.supports(IN_ARG_beta) ) {
    if(supports(IN_ARG_beta) || !ignoreUnsupported)
      set_beta(inArgs.get_beta());
  }
}

template<class Scalar>
std::string ModelEvaluatorBase::InArgs<Scalar>::description() const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  std::ostringstream oss;
  oss
    << "Thyra::ModelEvaluatorBase::InArgs<"<<ST::name()<<">"
    << "{"
    << "model="<<modelEvalDescription_
    << ",Np="<<Np()
    << "}";
  return oss.str();
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::describe(
  Teuchos::FancyOStream &out_arg, const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::OSTab;
  typedef Teuchos::RefCountPtr<const VectorBase<Scalar> > CV_ptr;
  if(verbLevel == Teuchos::VERB_NONE)
    return;
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::rcp(&out_arg,false);
  OSTab tab(out);
  *out <<"Thyra::ModelEvaluatorBase::InArgs<"<<ST::name()<<">:\n";
  tab.incrTab();
  *out <<"model = " << modelEvalDescription_ << "\n";
  *out <<"Np = " << Np() << "\n";
  switch(verbLevel) {
    case Teuchos::VERB_LOW:
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      if(this->supports(IN_ARG_x_dot) ) {
        *out << "x_dot =";
        CV_ptr x_dot = this->get_x_dot();
        if(x_dot.get())
          *out << "\n" << Teuchos::describe(*x_dot,verbLevel);
        else
          *out << " NULL\n";
      }
      if(this->supports(IN_ARG_x) ) {
        *out << "x =";
        CV_ptr x = this->get_x();
        if(x.get())
          *out << "\n" << Teuchos::describe(*x,verbLevel);
        else
          *out << " NULL\n";
      }
      for( int l = 0; l < Np(); ++l ) {
        *out << "p("<<l<<") =";
        CV_ptr p_l = this->get_p(l);
        if(p_l.get())
          *out << "\n" << Teuchos::describe(*p_l,verbLevel);
        else
          *out << " NULL\n";
      }
      // ToDo: Add output for more objects!
      break;
    }
    default:
      TEST_FOR_EXCEPT(true);
  }
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
    ,"model = \'"<<modelEvalDescription_<<"\': Error, arg="<<toString(arg)<<" is invalid!");
  supports_[arg] = supports;
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setSupports( const InArgs<Scalar>& inArgs )
{
  std::copy( &inArgs.supports_[0], &inArgs.supports_[0] + NUM_E_IN_ARGS_MEMBERS, &supports_[0] );
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setUnsupportsAndRelated( EInArgsMembers arg )
{
  this->_setSupports(arg,false);
  switch(arg) {
    case IN_ARG_x: {
      this->_setSupports(IN_ARG_x_dot,false);
      this->_setSupports(IN_ARG_x_dot_poly,false);
      this->_setSupports(IN_ARG_alpha,false);
      this->_setSupports(IN_ARG_beta,false);
      break;
    }
    default:
      TEST_FOR_EXCEPTION(
        true ,std::logic_error
        ,"Error, can handle args other than IN_ARG_x yet!"
        );
  }
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::assert_supports(EInArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"Thyra::ModelEvaluatorBase::InArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(arg): "
    "model = \'"<<modelEvalDescription_<<"\': Error, "
    "The argument arg = " << toString(arg) << " is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::assert_l(int l) const
{
  TEST_FOR_EXCEPTION(
    !( 0 <= l && l < Np() ), std::logic_error
    ,"Thyra::ModelEvaluatorBase::InArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_l(l): "
    " model = \'"<<modelEvalDescription_<<"\': Error, "
    "The parameter l = " << l << " is not in the range [0,"<<Np()-1<<"]!"
    );
}

//
// ModelEvaluatorBase::OutArgs
//

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::OutArgs()
  :isFailed_(false)
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
    ,"model = \'"<<modelEvalDescription_<<"\': Error, arg="<<toString(arg)<<" is invalid!"
    );
  return supports_[arg];
}

template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(EOutArgsDfDp arg, int l) const
{
  assert_l(l);
  return supports_DfDp_[l];
}

template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(EOutArgsDgDx arg, int j) const
{
  assert_j(j);
  return supports_DgDx_[j];
}

template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(EOutArgsDgDp arg, int j, int l) const
{
  assert_j(j);
  assert_l(l);
  return supports_DgDp_[ j*Np() + l ];
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
  g_[j] = g_j;
}

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_g(int j) const
{ 
  assert_j(j);
  return g_[j];
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
void ModelEvaluatorBase::OutArgs<Scalar>::set_W_op( const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &W_op )
{
  assert_supports(OUT_ARG_W_op);
  W_op_ = W_op;
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_W_op() const
{
  assert_supports(OUT_ARG_W_op);
  return W_op_;
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
  DfDp_[l] = DfDp_l;
}

template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::get_DfDp(int l) const
{
  assert_supports(OUT_ARG_DfDp,l);
  return DfDp_[l];
}

template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DfDp_properties(int l) const
{
  assert_supports(OUT_ARG_DfDp,l);
  return DfDp_properties_[l];
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DgDx( int j, const Derivative<Scalar> &DgDx_j )
{
  assert_supports(OUT_ARG_DgDx,j);
  DgDx_[j] = DgDx_j;
}

template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx(int j) const
{
  assert_supports(OUT_ARG_DgDx,j);
  return DgDx_[j];
}

template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx,j);
  return DgDx_properties_[j];
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DgDp( int j, int l, const Derivative<Scalar> &DgDp_j_l )
{
  assert_supports(OUT_ARG_DgDp,j,l);
  DgDp_[ j*Np() + l ] = DgDp_j_l;
}

template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDp(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp,j,l);
  return DgDp_[ j*Np() + l ];
}

template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDp_properties(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp,j,l);
  return DgDp_properties_[ j*Np() + l ];
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_f_poly( const Teuchos::RefCountPtr<Teuchos::Polynomial< VectorBase<Scalar> > > &f_poly )
{ f_poly_ = f_poly; }

template<class Scalar>
Teuchos::RefCountPtr<Teuchos::Polynomial< VectorBase<Scalar> > >
ModelEvaluatorBase::OutArgs<Scalar>::get_f_poly() const
{ return f_poly_; }


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::setArgs( const OutArgs<Scalar>& outArgs, bool ignoreUnsupported )
{
  const int min_Np = TEUCHOS_MIN(this->Np(),outArgs.Np());
  const int min_Ng = TEUCHOS_MIN(this->Ng(),outArgs.Ng());
  if( outArgs.supports(OUT_ARG_f) && outArgs.get_f().get() ) {
    if(supports(OUT_ARG_f) || !ignoreUnsupported)
      set_f(outArgs.get_f());
  }
  for( int j = 0; j < min_Ng; ++j ) {
    if(outArgs.get_g(j).get())
      set_g(j,outArgs.get_g(j));
  }
  if( outArgs.supports(OUT_ARG_W) && outArgs.get_W().get() ) {
    if(supports(OUT_ARG_W) || !ignoreUnsupported)
      set_W(outArgs.get_W());
  }
  for( int l = 0; l < min_Np; ++l ) {
    if( !outArgs.supports(OUT_ARG_DfDp,l).none() && !outArgs.get_DfDp(l).isEmpty() ) {
      if(!supports(OUT_ARG_DfDp,l).none() || !ignoreUnsupported)
        set_DfDp(l,outArgs.get_DfDp(l));
    }
  }
  for( int j = 0; j < min_Ng; ++j ) {
    if( !outArgs.supports(OUT_ARG_DgDx,j).none() && !outArgs.get_DgDx(j).isEmpty() ) {
      if(!supports(OUT_ARG_DgDx,j).none() || !ignoreUnsupported)
        set_DgDx(j,outArgs.get_DgDx(j));
    }
  }
  for( int l = 0; l < min_Np; ++l ) {
    for( int j = 0; j < min_Ng; ++j ) {
      if( !outArgs.supports(OUT_ARG_DgDp,j,l).none() && !outArgs.get_DgDp(j,l).isEmpty() ) {
        if(!supports(OUT_ARG_DgDp,j,l).none() || !ignoreUnsupported)
          set_DgDp(j,l,outArgs.get_DgDp(j,l));
      }
    }
  }
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::setFailed() const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  isFailed_ = true;
  if( this->supports(OUT_ARG_f) && this->get_f().get() ) {
    assign(&*this->get_f(),ST::nan());
  }
  for( int j = 0; j < this->Ng(); ++j ) {
    if(this->get_g(j).get())
      assign(&*this->get_g(j),ST::nan());
  }
  // ToDo: Set other objects to NaN as well!
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::isFailed() const
{
  return isFailed_;
}

template<class Scalar>
std::string ModelEvaluatorBase::OutArgs<Scalar>::description() const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  std::ostringstream oss;
  oss
    << "Thyra::ModelEvaluatorBase::OutArgs<"<<ST::name()<<">"
    << "{"
    << "model="<<modelEvalDescription_
    << ",Np="<<Np()
    << ",Ng="<<Ng()
    << "}";
  return oss.str();
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::describe(
  Teuchos::FancyOStream &out_arg, const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Teuchos::RefCountPtr<const VectorBase<Scalar> > CV_ptr;
  typedef Teuchos::RefCountPtr<const LinearOpWithSolveBase<Scalar> > CLOWS_ptr;
  typedef ModelEvaluatorBase MEB;
  if(verbLevel == Teuchos::VERB_NONE)
    return;
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::rcp(&out_arg,false);
  OSTab tab(out);
  *out <<"Thyra::ModelEvaluatorBase::OutArgs<"<<ST::name()<<">:\n";
  tab.incrTab();
  *out <<"model = " << modelEvalDescription_ << "\n";
  *out <<"Np = " << Np() << "\n";
  *out <<"Ng = " << Ng() << "\n";
  switch(verbLevel) {
    case Teuchos::VERB_LOW:
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      if(this->supports(OUT_ARG_f) ) {
        *out << "f =";
        CV_ptr f = this->get_f();
        if(f.get())
          *out << "\n" << Teuchos::describe(*f,verbLevel);
        else
          *out << " NULL\n";
      }
      for( int j = 0; j < Ng(); ++j ) {
        *out << "g("<<j<<") =";
        CV_ptr g_j = this->get_g(j);
        if(g_j.get())
          *out << "\n" << Teuchos::describe(*g_j,verbLevel);
        else
          *out << " NULL\n";
      }
      if(this->supports(OUT_ARG_W) ) {
        *out << "W =";
        CLOWS_ptr W = this->get_W();
        if(W.get())
          *out << "\n" << Teuchos::describe(*W,verbLevel);
        else
          *out << " NULL\n";
      }
      for( int l = 0; l < Np(); ++l ) {
        if(!this->supports(OUT_ARG_DfDp,l).none()) {
          *out << "DfDp("<<l<<") =";
          const MEB::Derivative<Scalar> DfDp_l = this->get_DfDp(l);
          if(!DfDp_l.isEmpty()) {
            *out << "\n";
            if(DfDp_l.getLinearOp().get()) {
              *out << Teuchos::describe(*DfDp_l.getLinearOp(),verbLevel);
            }
            else {
              OSTab(out).o()
                << "orientation="
                << toString(DfDp_l.getDerivativeMultiVector().getOrientation()) << "\n";
              *out << Teuchos::describe(*DfDp_l.getDerivativeMultiVector().getMultiVector(),verbLevel);
            }
          }
          else {
            *out << " NULL\n";
          }
        }
      }
      // ToDo: Add output for more objects?
      break;
    }
    default:
      TEST_FOR_EXCEPT(true);
  }
}

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
    ,"model = \'"<<modelEvalDescription_<<"\': Error, arg="<<toString(arg)<<" is invalid!"
    );
  supports_[arg] = supports;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports( EOutArgsDfDp arg, int l, const DerivativeSupport& supports )
{
  assert_l(l);
  supports_DfDp_[l] = supports;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports( EOutArgsDgDx arg, int j, const DerivativeSupport& supports )
{
  assert_j(j);
  supports_DgDx_[j] = supports;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports( EOutArgsDgDp arg, int j, int l, const DerivativeSupport& supports )
{
  assert_j(j);
  assert_l(l);
  supports_DgDp_[ j*Np()+ l ] = supports;
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
  DfDp_properties_[l] = properties;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DgDx_properties( int j, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDx,j);
  DgDx_properties_[j] = properties;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DgDp_properties( int j, int l, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDp,j,l);
  DgDp_properties_[ j*Np()+ l ] = properties;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports( const OutArgs<Scalar>& outArgs )
{
  typedef ModelEvaluatorBase MEB;
  const int Np = TEUCHOS_MIN(this->Np(),outArgs.Np()); 
  const int Ng = TEUCHOS_MIN(this->Ng(),outArgs.Ng()); 
  std::copy( &outArgs.supports_[0], &outArgs.supports_[0] + NUM_E_OUT_ARGS_MEMBERS, &supports_[0] );
  for( int l = 0; l < Np; ++l ) {
    DerivativeSupport ds = outArgs.supports(MEB::OUT_ARG_DfDp,l);
    this->_setSupports(MEB::OUT_ARG_DfDp,l,ds);
    if(!ds.none()) this->_set_DfDp_properties(l,outArgs.get_DfDp_properties(l));
  }
  for( int j = 0; j < Ng; ++j ) {
    DerivativeSupport ds = outArgs.supports(MEB::OUT_ARG_DgDx,j);
    this->_setSupports(MEB::OUT_ARG_DgDx,j,ds);
    if(!ds.none()) this->_set_DgDx_properties(j,outArgs.get_DgDx_properties(j));
  }
  for( int j = 0; j < Ng; ++j ) for( int l = 0; l < Np; ++l ) {
    DerivativeSupport ds = outArgs.supports(MEB::OUT_ARG_DgDp,j,l);
    this->_setSupports(MEB::OUT_ARG_DgDp,j,l,ds);
    if(!ds.none()) this->_set_DgDp_properties(j,l,outArgs.get_DgDp_properties(j,l));
  }
  if(this->supports(OUT_ARG_W))
    this->_set_W_properties(outArgs.get_W_properties());
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setUnsupportsAndRelated( EInArgsMembers arg )
{
  switch(arg) {
    case IN_ARG_x: {
      const int Ng = this->Ng();
      for( int j = 0; j < Ng; ++j )
        this->_setSupports(OUT_ARG_DgDx,j,DerivativeSupport());
      break;
    }
    default:
      TEST_FOR_EXCEPTION(
        true ,std::logic_error
        ,"Error, can handle args other than IN_ARG_x yet!"
        );
  }
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setUnsupportsAndRelated( EOutArgsMembers arg )
{
  this->_setSupports(arg,false);
  switch(arg) {
    case OUT_ARG_f: {
      this->_setSupports(OUT_ARG_W,false);
      this->_setSupports(OUT_ARG_W_op,false);
      this->_setSupports(OUT_ARG_f_poly,false);
      const int Np = this->Np();
      for( int l = 0; l < Np; ++l )
        this->_setSupports(OUT_ARG_DfDp,l,DerivativeSupport());
      break;
    }
    default:
      TEST_FOR_EXCEPTION(
        true ,std::logic_error
        ,"Error, can handle args other than OUT_ARG_f yet!"
        );
  }
}

// private

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(EOutArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !this->supports(arg), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(arg): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument arg = " << toString(arg) << " is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(EOutArgsDfDp arg, int l) const
{
  TEST_FOR_EXCEPTION(
    this->supports(arg,l).none(), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(OUT_ARG_DfDp,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DfDp(l) with index l = " << l << " is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(EOutArgsDgDx arg, int j) const
{
  TEST_FOR_EXCEPTION(
    this->supports(arg,j).none(), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(OUT_ARG_DgDx,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx(j) with index j = " << j << " is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(EOutArgsDgDp arg, int j, int l) const
{
  TEST_FOR_EXCEPTION(
    this->supports(arg,j,l).none(), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(OUT_ARG_DgDp,j,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDp(j,l) with indexes j = " << j << " and l = " << l << " is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_l(int l) const
{
  TEST_FOR_EXCEPTION(
    !( 0 <= l && l < Np() ), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_l(l): "
    "model = \'"<<modelEvalDescription_<<"\': Error, "
    "The parameter subvector p(l) index l = " << l << " is not in the range [0,"<<Np()-1<<"]!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_j(int j) const
{
  TEST_FOR_EXCEPTION(
    !( 0 <= j && j < Ng() ), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<" << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_j(j): "
    "model = \'"<<modelEvalDescription_<<"\': Error, "
    "The auxiliary function g(j) index j = " << j << " is not in the range [0,"<<Ng()-1<<"]!"
    );
}

//
// ModelEvaluatorBase::InArgsSetup
//

template<class Scalar>
ModelEvaluatorBase::InArgsSetup<Scalar>::InArgsSetup()
{}

template<class Scalar>
ModelEvaluatorBase::InArgsSetup<Scalar>::InArgsSetup( const InArgs<Scalar>& inArgs )
  :InArgs<Scalar>(inArgs)
{}

template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::setModelEvalDescription( const std::string &modelEvalDescription )
{ this->_setModelEvalDescription(modelEvalDescription); }

template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::set_Np(int Np)
{ this->_set_Np(Np); }

template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::setSupports( EInArgsMembers arg, bool supports )
{ this->_setSupports(arg,supports); }

template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::setSupports( const InArgs<Scalar>& inArgs )
{ this->_setSupports(inArgs); }

template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::setUnsupportsAndRelated( EInArgsMembers arg )
{ this->_setUnsupportsAndRelated(arg); }

//
// ModelEvaluatorBase::OutArgsSetup
//

template<class Scalar>
ModelEvaluatorBase::OutArgsSetup<Scalar>::OutArgsSetup()
{}

template<class Scalar>
ModelEvaluatorBase::OutArgsSetup<Scalar>::OutArgsSetup( const OutArgs<Scalar>& outArgs )
  :OutArgs<Scalar>(outArgs)
{}

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
{ this->_set_DfDp_properties(l,properties); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DgDx_properties( int j, const DerivativeProperties &properties )
{ this->_set_DgDx_properties(j,properties); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DgDp_properties( int j, int l, const DerivativeProperties &properties )
{ this->_set_DgDp_properties(j,l,properties); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports( const OutArgs<Scalar>& outArgs )
{ this->_setSupports(outArgs); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setUnsupportsAndRelated( EInArgsMembers arg )
{ this->_setUnsupportsAndRelated(arg); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setUnsupportsAndRelated( EOutArgsMembers arg )
{ this->_setUnsupportsAndRelated(arg); }

} // namespace Thyra

#endif // THYRA_MODEL_EVALUATOR_BASE_HPP
