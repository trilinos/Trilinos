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

#ifndef THYRA_MODEL_EVALUATOR_BASE_DECL_HPP
#define THYRA_MODEL_EVALUATOR_BASE_DECL_HPP


#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_Assert.hpp"

#ifdef HAVE_THYRA_ME_POLYNOMIAL
#  include "Teuchos_Polynomial.hpp"
#endif


namespace Thyra {


/** \brief Base subclass for <tt>ModelEvaluator</tt> that defines some basic
 * types.
 *
 * This non-templated base class is used for two very important reasons.
 *
 * First, a non-templated base class holding templated nested classes makes it
 * easier for client to form the names of the nested classes.  This also makes
 * it easier to access non-tempated enum types and values as well.  While most
 * of these nested types could have been defined outside of a base class, by
 * putting them in a base class, we get better namespace scoping and can
 * therefore use shorter names.
 *
 * Second, there are some protected nested classes (i.e. <tt>InArgsSetup</tt>
 * and <tt>OutArgsSetup</tt>0 that only subclasses should be able to access.
 * This makes the design very secure to help avoid bad usage of the nested
 * classes.
 *
 * \ingroup Thyra_Nonlinear_model_evaluator_interfaces_code_grp
 */
class ModelEvaluatorBase
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<ModelEvaluatorBase>
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

  /** \brief Concrete aggregate class for all output arguments computable by a
   * <tt>ModelEvaluator</tt> subclass object.
   *
   * The set of supported objects is returned from the <tt>supports()</tt>
   * function.
   *
   * A client can not directly set what input arguments are supported or not
   * supported.  Only a subclass of <tt>ModelEvaluator</tt> can do that
   * (through the <tt>InArgsSetup</tt> subclass).
   */
  template<class Scalar>
  class InArgs : public Teuchos::Describable {
  public:
    /** \brief .  */
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    /** \brief .  */
    InArgs();
    /** \brief Return the number of parameter subvectors <tt>p(l)</tt>
     * supported (<tt>Np >= 0</tt>).  */
    int Np() const;
    /** \brief Determines if an input argument is supported or not.  */
    bool supports(EInArgsMembers arg) const;
    /** \brief Precondition: <tt>supports(IN_ARG_x_dot)==true</tt>.  */
    void set_x_dot( const RCP<const VectorBase<Scalar> > &x_dot );
    /** \brief Precondition: <tt>supports(IN_ARG_x_dot)==true</tt>.  */
    RCP<const VectorBase<Scalar> > get_x_dot() const;
    /** \brief Precondition: <tt>supports(IN_ARG_x)==true</tt>.  */
    void set_x( const RCP<const VectorBase<Scalar> > &x );
    /** \brief Precondition: <tt>supports(IN_ARG_x)==true</tt>.  */
    RCP<const VectorBase<Scalar> > get_x() const;
#ifdef HAVE_THYRA_ME_POLYNOMIAL
    /** \brief Precondition: <tt>supports(IN_ARG_x_poly)==true</tt>.  */
    void set_x_poly( 
      const RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > &x_poly );
    /** \brief Precondition: <tt>supports(IN_ARG_x)==true</tt>.  */
    RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > get_x_poly() const;
    /** \brief Precondition: <tt>supports(IN_ARG_x_dot_poly)==true</tt>.  */
    void set_x_dot_poly(
      const RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > &x_dot_poly );
    /** \brief Precondition: <tt>supports(IN_ARG_x_dot_poly)==true</tt>.  */
    RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > get_x_dot_poly() const;
#endif // HAVE_THYRA_ME_POLYNOMIAL
    /** \brief Set <tt>p(l)</tt> where <tt>0 <= l && l < this->Np()</tt>.  */
    void set_p( int l, const RCP<const VectorBase<Scalar> > &p_l );
    /** \brief Get <tt>p(l)</tt> where <tt>0 <= l && l < this->Np()</tt>.  */
    RCP<const VectorBase<Scalar> > get_p(int l) const;
    /** \brief Precondition: <tt>supports(IN_ARG_t)==true</tt>.  */
    void set_t( ScalarMag t );
    /** \brief .Precondition: <tt>supports(IN_ARG_t)==true</tt>  */
    ScalarMag get_t() const;
    /** \brief Precondition: <tt>supports(IN_ARG_alpha)==true</tt>.  */
    void set_alpha( Scalar alpha );
    /** \brief Precondition: <tt>supports(IN_ARG_alph)==true</tt>.  */
    Scalar get_alpha() const;
    /** \brief Precondition: <tt>supports(IN_ARG_beta)==true</tt>.  */
    void set_beta( Scalar beta );
    /** \brief Precondition: <tt>supports(IN_ARG_beta)==true</tt>.  */
    Scalar get_beta() const;
    /** \brief Set non-null arguments (does not overwrite non-NULLs with
     * NULLs) .  */
    void setArgs(
      const InArgs<Scalar>& inArgs, bool ignoreUnsupported = false,
      bool cloneObjects = false
      );
    /** \brief Assert that two InArgs objects have the same support.*/
    void assertSameSupport( const InArgs<Scalar> &inArgs ) const;
    /** \brief . */
    std::string modelEvalDescription() const;
    /** \brief . */
    std::string description() const;
    /** \brief Create a more detailed description along about this object and
     * the ModelEvaluator that created it. */
    void describe(
      Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel
      ) const;
  protected:
    /** \brief . */
    void _setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void _set_Np(int Np);
    /** \brief . */
    void _setSupports( EInArgsMembers arg, bool supports );
    /** \brief . */
    void _setSupports( const InArgs<Scalar>& inputInArgs, const int Np );
    /** \brief . */
    void _setUnsupportsAndRelated( EInArgsMembers arg );
  private:
    // types
    typedef Teuchos::Array<RCP<const VectorBase<Scalar> > > p_t;
    // data
    std::string modelEvalDescription_;
    RCP<const VectorBase<Scalar> > x_dot_;
    RCP<const VectorBase<Scalar> > x_;
#ifdef HAVE_THYRA_ME_POLYNOMIAL
    RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > x_dot_poly_;
    RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > x_poly_;
#endif // HAVE_THYRA_ME_POLYNOMIAL
    p_t p_;
    ScalarMag t_;
    Scalar alpha_;
    Scalar beta_;
    bool supports_[NUM_E_IN_ARGS_MEMBERS];
    // functions
    void assert_supports(EInArgsMembers arg) const;
    void assert_l(int l) const;
  };

  /** \brief . */
  enum EDerivativeMultiVectorOrientation {
    DERIV_MV_JACOBIAN_FORM, ///< Jacobian form DhDz (nz columns of h_space vectors)
    DERIV_MV_GRADIENT_FORM, ///< Gradient form DhDz^T (nh columns of z_space vectors)
    DERIV_MV_BY_COL = DERIV_MV_JACOBIAN_FORM, ///< Deprecated!
    DERIV_TRANS_MV_BY_ROW = DERIV_MV_GRADIENT_FORM ///< Deprecated!
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
    /** \brief . */
    bool isSameSupport(const DerivativeSupport &derivSupport) const
      {
        return (
          supportsLinearOp_ == derivSupport.supportsLinearOp_
          && supportsMVByCol_ == derivSupport.supportsMVByCol_
          && supportsTransMVByRow_ == derivSupport.supportsTransMVByRow_
          );
      } 
    /** \brief . */
    std::string description() const;
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

  /** \brief Simple public strict containing properties of a derivative
   * object. */
  struct DerivativeProperties {
    /** \brief . */
    EDerivativeLinearity     linearity;
    /** \brief . */
    ERankStatus              rank;
    /** \brief . */
    bool                     supportsAdjoint;
    /** \brief . */
    DerivativeProperties()
      :linearity(DERIV_LINEARITY_UNKNOWN),
       rank(DERIV_RANK_UNKNOWN),supportsAdjoint(false)
      {}
    /** \brief . */
    DerivativeProperties(
      EDerivativeLinearity in_linearity, ERankStatus in_rank,
      bool in_supportsAdjoint
      )
      :linearity(in_linearity),rank(in_rank),
       supportsAdjoint(in_supportsAdjoint)
      {}
  };

  /** \brief Simple aggregate class for a derivative object represented as a
   * column-wise multi-vector or its transpose as a row-wise multi-vector.
   */
  template<class Scalar>
  class DerivativeMultiVector {
  public:
    /** \brief . */
    DerivativeMultiVector()
      :orientation_(DERIV_MV_BY_COL)
      {}
    /** \brief . */
    DerivativeMultiVector(
      const RCP<MultiVectorBase<Scalar> > &mv
      ,const EDerivativeMultiVectorOrientation orientation = DERIV_MV_BY_COL
      ) : mv_(mv.assert_not_null()), orientation_(orientation) {}
    /** \brief . */
    void changeOrientation( const EDerivativeMultiVectorOrientation orientation )
      { orientation_ = orientation; };
    /** \brief . */
    const DerivativeMultiVector<Scalar>& assert_not_null() const
      { mv_.assert_not_null(); return *this; }
    /** \brief . */
    RCP<MultiVectorBase<Scalar> > getMultiVector() const
      { return mv_; }
    /** \brief . */
    EDerivativeMultiVectorOrientation getOrientation() const
      { return orientation_; }
    /** \brief . */
    std::string description() const;
    /** \brief . */
    void describe( 
      Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel
      ) const;
  private:
    RCP<MultiVectorBase<Scalar> > mv_;
    EDerivativeMultiVectorOrientation orientation_;
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
    Derivative( const RCP<LinearOpBase<Scalar> > &lo )
      : lo_(lo.assert_not_null()) {}
    /** \brief . */
    Derivative(
      const RCP<MultiVectorBase<Scalar> > &mv,
      const EDerivativeMultiVectorOrientation orientation = DERIV_MV_BY_COL
      ) : dmv_(mv,orientation) {}
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
    RCP<LinearOpBase<Scalar> > getLinearOp() const
      { return lo_; }
    /** \brief . */
    RCP<MultiVectorBase<Scalar> > getMultiVector() const
      { return dmv_.getMultiVector(); }
    /** \brief . */
    EDerivativeMultiVectorOrientation getMultiVectorOrientation() const
      { return dmv_.getOrientation(); }
    /** \brief . */
    DerivativeMultiVector<Scalar> getDerivativeMultiVector() const
      { return dmv_; }
    /** \brief Returns true if the form of the derivative contained here is
     * supported by deriveSupport.
     */
    bool isSupportedBy( const DerivativeSupport &derivSupport ) const
      {
        // If there is not derivative support then we will return false!
        if (derivSupport.none())
          return false;
        if (!is_null(getMultiVector())) {
          return derivSupport.supports(getMultiVectorOrientation());
        }
        else if(!is_null(getLinearOp())) {
          return derivSupport.supports(DERIV_LINEAR_OP);
        }
        // If nothing is set then of course we support that!
        return true;
      }
    /** \brief . */
    std::string description() const;
    /** \brief . */
    void describe( 
      Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel
      ) const;
  private:
    RCP<LinearOpBase<Scalar> > lo_;
    DerivativeMultiVector<Scalar> dmv_;
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
  enum EOutArgsDgDx_dot {
    OUT_ARG_DgDx_dot   ///< .
  };

  /** \brief . */
  enum EOutArgsDgDx {
    OUT_ARG_DgDx   ///< .
  };

  /** \brief . */
  enum EOutArgsDgDp {
    OUT_ARG_DgDp   ///< .
  };
  
  /** \brief Concrete aggregate class for all output arguments computable by a
   * <tt>ModelEvaluator</tt> subclass object.
   *
   * Note that <tt>const OutArgs</tt> object means that a client can not
   * change what output objects are being pointed to but they can still change
   * the states of the contained objects.  This is slight variation on the
   * concept of logical const-ness in C++ but it is totally consistent with
   * the vary nature of this class.
   *
   * In addition to storing the output objects themselves, this class also
   * allows the storage if the properties of some of the objects as well.
   * Therefore, objects of this type are used to communicate a lot of
   * different information about the output functions and derivatives
   * supported by a model.  It tells clients what functions and derivatives
   * are supported and what the know properties are.
   *
   * A client can not directly set what input arguments are supported or not
   * supported.  Only a subclass of <tt>ModelEvaluator</tt> can do that
   * (through the <tt>OutArgsSetup</tt> subclass).
   */
  template<class Scalar>
  class OutArgs : public Teuchos::Describable {
  public:
    /** \brief .  */
    OutArgs();
    /** \brief Return the number of parameter subvectors <tt>p(l)</tt>
     * supported (<tt>Np >= 0</tt>).  */
    int Np() const;
    /** \brief Return the number of axillary response functions
     * <tt>g(j)(...)</tt> supported (<tt>Ng >= 0</tt>).  */
    int Ng() const;
    /** \brief Determine if an input argument is supported or not.  */
    bool supports(EOutArgsMembers arg) const;
    /** \brief Determine if <tt>DfDp(l)</tt> is supported or not, where <tt>0
     * <= l && l < Np()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDfDp arg, int l) const;
    /** \brief Determine if <tt>DgDx_dot(j)</tt> is supported or not, <tt>0 <= j
     * && j < Ng()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDx_dot arg, int j) const;
    /** \brief Determine if <tt>DgDx(j)</tt> is supported or not, <tt>0 <= j
     * && j < Ng()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDx arg, int j) const;
    /** \brief Determine if <tt>DgDp(j,l)</tt> is supported or not, <tt>0 <= j
     * && j < Ng()</tt> and <tt>0 <= l && l < Np()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDp arg, int j, int l) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_f)==true</tt>.  */
    void set_f( const RCP<VectorBase<Scalar> > &f );
    /** \brief Precondition: <tt>supports(OUT_ARG_f)==true</tt>.  */
    RCP<VectorBase<Scalar> > get_f() const;
    /** \brief Precondition: <tt>supports(OUT_ARG_g)==true</tt>.  */
    void set_g( int j, const RCP<VectorBase<Scalar> > &g_j );
    /** \brief Precondition: <tt>supports(OUT_ARG_g)==true</tt>..  */
    RCP<VectorBase<Scalar> > get_g(int j) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_W)==true</tt>.  */
    void set_W( const RCP<LinearOpWithSolveBase<Scalar> > &W );
    /** \brief Precondition: <tt>supports(OUT_ARG_W)==true</tt>.  */
    RCP<LinearOpWithSolveBase<Scalar> > get_W() const;
    /** \brief Precondition: <tt>supports(OUT_ARG_W_op)==true</tt>.  */
    void set_W_op( const RCP<LinearOpBase<Scalar> > &W_op );
    /** \brief Precondition: <tt>supports(OUT_ARG_W_op)==true</tt>.  */
    RCP<LinearOpBase<Scalar> > get_W_op() const;
    /** \brief Return the known properties of <tt>W</tt> (precondition:
     * <tt>supports(OUT_ARG_f)==true</tt>). */
    DerivativeProperties get_W_properties() const;
    /** \brief Precondition: <tt>supports(OUT_ARG_DfDp,l)==true</tt>.  */
    void set_DfDp(int l,  const Derivative<Scalar> &DfDp_l);
    /** \brief Precondition: <tt>supports(OUT_ARG_DfDp,l)==true</tt>.  */
    Derivative<Scalar> get_DfDp(int l) const;
    /** \brief Return the know properties of <tt>DfDp(l)</tt> (precondition:
     * <tt>supports(OUT_ARG_DfDp,l)==true</tt>). */
    DerivativeProperties get_DfDp_properties(int l) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_DgDx_dot,j)==true</tt>.  */
    void set_DgDx_dot(int j, const Derivative<Scalar> &DgDx_dot_j);
    /** \brief Precondition: <tt>supports(OUT_ARG_DgDx_dot,j)==true</tt>.  */
    Derivative<Scalar> get_DgDx_dot(int j) const;
    /** \brief Return the know properties of <tt>DgDx_dot(j)</tt> (precondition:
     * <tt>supports(OUT_ARG_DgDx_dot,j)==true</tt>). */
    DerivativeProperties get_DgDx_dot_properties(int j) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_DgDx,j)==true</tt>.  */
    void set_DgDx(int j, const Derivative<Scalar> &DgDx_j);
    /** \brief Precondition: <tt>supports(OUT_ARG_DgDx,j)==true</tt>.  */
    Derivative<Scalar> get_DgDx(int j) const;
    /** \brief Return the know properties of <tt>DgDx(j)</tt> (precondition:
     * <tt>supports(OUT_ARG_DgDx,j)==true</tt>). */
    DerivativeProperties get_DgDx_properties(int j) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_DgDp,j,l)==true</tt>.  */
    void set_DgDp( int j, int l, const Derivative<Scalar> &DgDp_j_l );
    /** \brief Precondition: <tt>supports(OUT_ARG_DgDp,j,l)==true</tt>.  */
    Derivative<Scalar> get_DgDp(int j, int l) const;
    /** \brief Return the know properties of <tt>DgDp(j,l)</tt> (precondition:
     * <tt>supports(OUT_ARG_DgDp,j,l)==true</tt>). */
    DerivativeProperties get_DgDp_properties(int j, int l) const;
#ifdef HAVE_THYRA_ME_POLYNOMIAL
    /** \brief Precondition: <tt>supports(OUT_ARG_f_poly)==true</tt>.  */
    void set_f_poly( const RCP<Teuchos::Polynomial< VectorBase<Scalar> > > &f_poly );
    /** \brief Precondition: <tt>supports(OUT_ARG_f_poly)==true</tt>.  */
    RCP<Teuchos::Polynomial< VectorBase<Scalar> > > get_f_poly() const;
#endif // HAVE_THYRA_ME_POLYNOMIAL
    /** \brief Set all arguments fron <tt>outArgs</tt> into <tt>*this</tt>.
     *
     * If <tt>ignoreUnsupported==true</tt>, then arguments in <tt>outArgs</tt>
     * that are not supported in <tt>*this</tt> will be ignored.  Othereise,
     * if an unsupported argument is set, then an exception will be thrown.*/
    void setArgs( const OutArgs<Scalar>& outArgs, bool ignoreUnsupported = false );
    /** \brief Set that the evaluation as a whole failed.
     *
     * Note that this function is declared as <tt>const</tt> even through it
     * technically changes the state of <tt>*this</tt> object.  This was done
     * so that this property could be set by a <tt>ModelEvaluator</tt>
     * subclass in <tt>evalModel()</tt> which takes a <tt>const OutArgs</tt>
     * object.  This is consistent with the behavior of the rest of a
     * <tt>const OutArgs</tt> object in that a client is allowed to change the
     * state of objects through a <tt>const OutArgs</tt> object, they just
     * can't change what objects are pointed to.
     */
    void setFailed() const;
    /** \brief Return if the evaluation failed or not.
     *
     * If the evaluation failed, no assumptions should be made at all about
     * the state of the output objects.
     */
    bool isFailed() const;
    /** \brief . */
    bool isEmpty() const;
    /** \brief Assert that two OutArgs objects have the same support.*/
    void assertSameSupport( const OutArgs<Scalar> &outArgs ) const;
    /** \brief . */
    std::string modelEvalDescription() const;
    /** \brief . */
    std::string description() const;
    /** \brief Create a more detailed description along about this object and
     * the ModelEvaluator that created it. */
    void describe(
      Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel
      ) const;
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
    void _setSupports( EOutArgsDgDx_dot arg, int j, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDx arg, int j, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDp arg, int j, int l, const DerivativeSupport& );
    /** \brief . */
    void _set_W_properties( const DerivativeProperties &properties );
    /** \brief . */
    void _set_DfDp_properties( int l, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDx_dot_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDx_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDp_properties( int j, int l, const DerivativeProperties &properties );
    /** \brief . */
    void _setSupports( const OutArgs<Scalar>& inputOutArgs );
    /** \brief . */
    void _setUnsupportsAndRelated( EInArgsMembers arg );
    /** \brief . */
    void _setUnsupportsAndRelated( EOutArgsMembers arg );
  private:
    // types
    typedef Teuchos::Array<RCP<VectorBase<Scalar> > > g_t;
    typedef Teuchos::Array<Derivative<Scalar> > deriv_t;
    typedef Teuchos::Array<DerivativeProperties> deriv_properties_t;
    typedef Teuchos::Array<DerivativeSupport> supports_t;
    // data
    std::string modelEvalDescription_;
    bool supports_[NUM_E_OUT_ARGS_MEMBERS];
    supports_t supports_DfDp_; // Np
    supports_t supports_DgDx_dot_; // Ng
    supports_t supports_DgDx_; // Ng
    supports_t supports_DgDp_; // Ng x Np
    RCP<VectorBase<Scalar> > f_;
    g_t g_; // Ng
    RCP<LinearOpWithSolveBase<Scalar> > W_;
    RCP<LinearOpBase<Scalar> > W_op_;
    DerivativeProperties W_properties_;
    deriv_t DfDp_; // Np
    deriv_properties_t DfDp_properties_; // Np
    deriv_t DgDx_dot_; // Ng
    deriv_t DgDx_; // Ng
    deriv_properties_t DgDx_dot_properties_; // Ng
    deriv_properties_t DgDx_properties_; // Ng
    deriv_t DgDp_; // Ng x Np
    deriv_properties_t DgDp_properties_; // Ng x Np
#ifdef HAVE_THYRA_ME_POLYNOMIAL
   RCP<Teuchos::Polynomial< VectorBase<Scalar> > > f_poly_;
#endif // HAVE_THYRA_ME_POLYNOMIAL
    mutable bool isFailed_;
    // functions
    void assert_supports(EOutArgsMembers arg) const;
    void assert_supports(
      EOutArgsDfDp arg, int l,
      const Derivative<Scalar> &deriv = Derivative<Scalar>()
      ) const;
    void assert_supports(
      EOutArgsDgDx_dot arg, int j,
      const Derivative<Scalar> &deriv = Derivative<Scalar>()
      ) const;
    void assert_supports(
      EOutArgsDgDx arg, int j,
      const Derivative<Scalar> &deriv = Derivative<Scalar>()
      ) const;
    void assert_supports(
      EOutArgsDgDp arg, int j, int l,
      const Derivative<Scalar> &deriv = Derivative<Scalar>()
      ) const;
    void assert_l(int l) const;
    void assert_j(int j) const;
  };

  //@}

// Added since at least gcc 3.3.4 does not do the right thing here!
#ifdef HAVE_PROTECTED_NESTED_TEMPLATE_CLASS_ACCESS
protected:
#endif

  /** \name Protected types */
  //@{

  /** \brief Protected subclass of <tt>InArgs</tt> that only
   * <tt>ModelEvaluator</tt> subclasses can access to set up the selection of
   * supported input arguments.
   *
   * Objects of this type must be created in overrides of
   * <tt>ModelEvaluator::createInArgs()</tt>.
   */
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
    void setSupports( const InArgs<Scalar>& inputInArgs, const int Np = -1 );
    /** \brief . */
    void setUnsupportsAndRelated( EInArgsMembers arg );
  };

  /** \brief Protected subclass of <tt>OutArgs</tt> that only
   * <tt>ModelEvaluator</tt> subclasses can access to set up the selection of
   * supported input arguments.
   *
   * Objects of this type must be created in overrides of
   * <tt>ModelEvaluator::createOutArgs()</tt>.
   */
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
    void setSupports(EOutArgsDgDx_dot arg, int j, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDx arg, int j, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDp arg, int j, int l, const DerivativeSupport& );
    /** \brief . */
    void set_W_properties( const DerivativeProperties &properties );
    /** \brief . */
    void set_DfDp_properties( int l, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDx_dot_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDx_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDp_properties( int j, int l, const DerivativeProperties &properties );
    /** \brief . */
    void setSupports( const OutArgs<Scalar>& inputOutArgs );
   /** \brief . */
    void setUnsupportsAndRelated( EInArgsMembers arg );
    /** \brief Must be called after the above function. */
    void setUnsupportsAndRelated( EOutArgsMembers arg );
   };

  //@}

};


/** \relates ModelEvaluatorBase */
std::string toString(ModelEvaluatorBase::EInArgsMembers);


/** \relates ModelEvaluatorBase */
std::string toString(ModelEvaluatorBase::EOutArgsMembers);


/** \relates ModelEvaluatorBase */
std::string toString(
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  );


/** \relates ModelEvaluatorBase */
ModelEvaluatorBase::EDerivativeMultiVectorOrientation
getOtherDerivativeMultiVectorOrientation(
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  );


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
    case ModelEvaluatorBase::IN_ARG_x_dot_poly:
      return "IN_ARG_x_dot_poly";
    case ModelEvaluatorBase::IN_ARG_x_poly:
      return "IN_ARG_x_poly";
    case ModelEvaluatorBase::IN_ARG_t:
      return "IN_ARG_t";
    case ModelEvaluatorBase::IN_ARG_alpha:
      return "IN_ARG_alpha";
    case ModelEvaluatorBase::IN_ARG_beta:
      return "IN_ARG_beta";
#ifdef TEUCHOS_DEBUG
    default:
      TEST_FOR_EXCEPT(true);
#endif
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
    case ModelEvaluatorBase::OUT_ARG_W_op:
      return "OUT_ARG_W_op";
    case ModelEvaluatorBase::OUT_ARG_f_poly:
      return "OUT_ARG_f_poly";
#ifdef TEUCHOS_DEBUG
    default:
      TEST_FOR_EXCEPT(true);
#endif
  }
  return ""; // Will never be executed!
}


inline
std::string Thyra::toString(
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  )
{
  switch(orientation) {
    case ModelEvaluatorBase::DERIV_MV_BY_COL:
      return "DERIV_MV_BY_COL";
    case ModelEvaluatorBase::DERIV_TRANS_MV_BY_ROW:
      return "DERIV_TRANS_MV_BY_ROW";
#ifdef TEUCHOS_DEBUG
    default:
      TEST_FOR_EXCEPT(true);
#endif
  }
  return ""; // Should never execute this!
}


inline
Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation
Thyra::getOtherDerivativeMultiVectorOrientation(
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation
  )
{
  switch(orientation) {
    case ModelEvaluatorBase::DERIV_MV_BY_COL:
      return ModelEvaluatorBase::DERIV_TRANS_MV_BY_ROW;
    case ModelEvaluatorBase::DERIV_TRANS_MV_BY_ROW:
      return ModelEvaluatorBase::DERIV_MV_BY_COL;
#ifdef TEUCHOS_DEBUG
    default:
      TEST_FOR_EXCEPT(true);
#endif
  }
  return ModelEvaluatorBase::DERIV_MV_BY_COL; // Should never execute this!
}


#endif // THYRA_MODEL_EVALUATOR_BASE_DECL_HPP
