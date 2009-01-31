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
#include "Thyra_PolynomialVectorTraits.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_Polynomial.hpp"
#include "Teuchos_Assert.hpp"


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
    RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > x_dot_poly_;
    RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > x_poly_;
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
    /** \brief Precondition: <tt>supports(OUT_ARG_f_poly)==true</tt>.  */
    void set_f_poly( const RCP<Teuchos::Polynomial< VectorBase<Scalar> > > &f_poly );
    /** \brief Precondition: <tt>supports(OUT_ARG_f_poly)==true</tt>.  */
    RCP<Teuchos::Polynomial< VectorBase<Scalar> > > get_f_poly() const;
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
    RCP<Teuchos::Polynomial< VectorBase<Scalar> > > f_poly_;
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


/** \defgroup Thyra_MEB_helper_functions_grp Helper functions for Thyra::ModelEvaluatorBase.
 *
 */
//@{


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


// //////////////////////////////////
// Definitions


namespace Thyra {


namespace ModelEvaluatorHelperPack {


template<class Scalar>
inline
RCP<const Thyra::VectorBase<Scalar> >
condCloneVec(
  const RCP<const Thyra::VectorBase<Scalar> > &vec,
  bool cloneObject
  )
{
  if(cloneObject)
    return vec->clone_v();
  return vec;
}


} // namespace ModelEvaluatorHelperPack


//
// ModelEvaluatorBase::DerivativeSupport
//


inline
std::string
ModelEvaluatorBase::DerivativeSupport::description() const
{
  std::ostringstream oss;
  oss << "DerivativeSupport{";
  if (none()) {
    oss << "none";
  }
  else {
    bool wroteOutput = false;
    if (supportsLinearOp_) {
      oss << "DERIV_LINEAR_OP";
      wroteOutput = true;
    }
    if (supportsMVByCol_) {
      oss << (wroteOutput?",":"") << toString(DERIV_MV_BY_COL);
      wroteOutput = true;
    }
    if (supportsTransMVByRow_) {
      oss << (wroteOutput?",":"") << toString(DERIV_TRANS_MV_BY_ROW);
      wroteOutput = true;
    }
  }
  oss << "}";
  return oss.str();
}
// 2007/09/08: rabartl: ToDo: Above: I really should move this function
// definition into a *.cpp file since it is not templated and it is too long
// to be inlined.  I am just making it inline so I can avoid this for now.


//
// ModelEvaluatorBase::InArgs
//


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>::InArgs()
  :modelEvalDescription_("WARNING!  THIS INARGS OBJECT IS UNINITALIZED!")
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
bool ModelEvaluatorBase::InArgs<Scalar>::supports(EInArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_
    <<"\': Error, arg="<<toString(arg)<<" is invalid!"
    );
  return supports_[arg];
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x_dot(
  const RCP<const VectorBase<Scalar> > &x_dot
  )
{ assert_supports(IN_ARG_x_dot); x_dot_ = x_dot; }


template<class Scalar>
RCP<const VectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_x_dot() const
{ assert_supports(IN_ARG_x_dot); return x_dot_; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x(
  const RCP<const VectorBase<Scalar> > &x
  )
{ assert_supports(IN_ARG_x); x_ = x; }


template<class Scalar>
RCP<const VectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_x() const
{ assert_supports(IN_ARG_x); return x_; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x_dot_poly(
  const RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > &x_dot_poly
  )
{ assert_supports(IN_ARG_x_dot_poly); x_dot_poly_ = x_dot_poly; }


template<class Scalar>
RCP<const Teuchos::Polynomial< VectorBase<Scalar> > >
ModelEvaluatorBase::InArgs<Scalar>::get_x_dot_poly() const
{ assert_supports(IN_ARG_x_dot_poly); return x_dot_poly_; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x_poly(
  const RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > &x_poly
  )
{ assert_supports(IN_ARG_x_poly); x_poly_ = x_poly; }


template<class Scalar>
RCP<const Teuchos::Polynomial< VectorBase<Scalar> > >
ModelEvaluatorBase::InArgs<Scalar>::get_x_poly() const
{ assert_supports(IN_ARG_x_poly); return x_poly_; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_p(
  int l, const RCP<const VectorBase<Scalar> > &p_l
  )
{ assert_l(l); p_[l] = p_l; }


template<class Scalar>
RCP<const VectorBase<Scalar> >
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
void ModelEvaluatorBase::InArgs<Scalar>::setArgs(
  const InArgs<Scalar>& inArgs, bool ignoreUnsupported, bool cloneObjects
  )
{
  using ModelEvaluatorHelperPack::condCloneVec;
  if( inArgs.supports(IN_ARG_x_dot) && inArgs.get_x_dot().get() ) {
    if(supports(IN_ARG_x_dot) || !ignoreUnsupported)
      set_x_dot(condCloneVec(inArgs.get_x_dot(),cloneObjects));
  }
  if( inArgs.supports(IN_ARG_x) && inArgs.get_x().get() ) {
    if(supports(IN_ARG_x) || !ignoreUnsupported)
      set_x(condCloneVec(inArgs.get_x(),cloneObjects));
  }
  if( inArgs.supports(IN_ARG_x_dot_poly) && inArgs.get_x_dot_poly().get() ) {
    if(supports(IN_ARG_x_dot_poly) || !ignoreUnsupported) {
      TEST_FOR_EXCEPT(
        cloneObjects && "Have not implemented cloning for x_dot_poly yet!" );
      set_x_dot_poly(inArgs.get_x_dot_poly());
    }
  }
  if( inArgs.supports(IN_ARG_x_poly) && inArgs.get_x_poly().get() ) {
    if(supports(IN_ARG_x_poly) || !ignoreUnsupported) {
      TEST_FOR_EXCEPT(
        cloneObjects && "Have not implemented cloning for x_poly yet!" );
      set_x_poly(inArgs.get_x_poly());
    }
  }
  const int min_Np = TEUCHOS_MIN(this->Np(),inArgs.Np());
  for( int l = 0; l < min_Np; ++l ) {
    if(inArgs.get_p(l).get())
      set_p(l,condCloneVec(inArgs.get_p(l),cloneObjects));
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
void ModelEvaluatorBase::InArgs<Scalar>::assertSameSupport(
  const InArgs<Scalar> &inArgs
  ) const
{
  for ( int inArg_i = 0; inArg_i < NUM_E_IN_ARGS_MEMBERS; ++inArg_i ) {
    const EInArgsMembers inArg_arg = static_cast<EInArgsMembers>(inArg_i);
    const std::string inArg_name = toString(inArg_arg);
    TEST_FOR_EXCEPTION(
      supports(inArg_arg) != inArgs.supports(inArg_arg), std::logic_error,
      "Error, the input argument "<<inArg_name<<" with support "<<inArgs.supports(inArg_arg)<<"\n"
      "in the InArgs object for the model:\n\n"
      "  "<<inArgs.modelEvalDescription()<<"\n\n"
      "is not the same the argument "<<inArg_name<<" with support "<<supports(inArg_arg)<<"\n"
      "in the InArgs object for the model:\n\n"
      "  "<<modelEvalDescription()<<"\n\n"
      "and these two InArgs objects are not compatible!"
      );
  }
  TEUCHOS_ASSERT_EQUALITY( this->Np(), inArgs.Np() );
}


template<class Scalar>
std::string ModelEvaluatorBase::InArgs<Scalar>::modelEvalDescription() const
{
  return modelEvalDescription_;
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
  using std::endl;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::OSTab;
  using Teuchos::describe;
  using Teuchos::includesVerbLevel;
  typedef RCP<const VectorBase<Scalar> > CV_ptr;

  if(verbLevel == Teuchos::VERB_NONE)
    return;

  RCP<Teuchos::FancyOStream>
    out = Teuchos::rcp(&out_arg,false);
  const bool dump_x = includesVerbLevel(verbLevel,Teuchos::VERB_HIGH);
  const Teuchos::EVerbosityLevel x_verbLevel =
    dump_x?Teuchos::VERB_EXTREME:verbLevel;
  const bool print_x_nrm = includesVerbLevel(verbLevel,Teuchos::VERB_LOW);
  const bool dump_p = includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM);
  const Teuchos::EVerbosityLevel p_verbLevel =
    dump_p?Teuchos::VERB_EXTREME:verbLevel;
  const bool print_p_nrm = includesVerbLevel(verbLevel,Teuchos::VERB_LOW);
  OSTab tab(out);

  *out <<"Thyra::ModelEvaluatorBase::InArgs<"<<ST::name()<<">:\n";
  tab.incrTab();

  *out <<"model = " << modelEvalDescription_ << "\n";
  *out <<"Np = " << Np() << "\n";

  CV_ptr x_dot;
  if ( this->supports(IN_ARG_x_dot) && !is_null(x_dot=get_x_dot()) ) {
    *out << "x_dot = " << Teuchos::describe(*x_dot,x_verbLevel);
    if (print_x_nrm)
      *out << "||x_dot|| = " << norm(*x_dot) << endl;
  }

  CV_ptr x;
  if ( this->supports(IN_ARG_x) && !is_null(x=get_x()) ) {
    *out << "x = " << Teuchos::describe(*x,x_verbLevel);
    if (print_x_nrm)
      *out << "||x|| = " << norm(*x) << endl;
  }
  
  if (print_x_nrm) {
    for( int l = 0; l < Np(); ++l ) {
      CV_ptr p_l;
      if ( !is_null(p_l = this->get_p(l)) ) {
        *out << "p("<<l<<") = " << Teuchos::describe(*p_l,p_verbLevel);
        if (print_p_nrm)
          *out << "||p("<<l<<")|| = " << norm(*p_l) << endl;
      }
    }
  }
  
  if (includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM)) {
    if (this->supports(IN_ARG_t)) {
      *out << "t = " << t_ << endl;
    }
    if (this->supports(IN_ARG_alpha)) {
      *out << "alpha = " << alpha_ << endl;
    }
    if (this->supports(IN_ARG_beta)) {
      *out << "beta = " << beta_ << endl;
    }
  }
  
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setModelEvalDescription(
  const std::string &modelEvalDescription
  )
{
  modelEvalDescription_ = modelEvalDescription;
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_set_Np(int Np)
{
  p_.resize(Np);
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setSupports(
  EInArgsMembers arg, bool supports
  )
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_
    <<"\': Error, arg="<<toString(arg)<<" is invalid!");
  supports_[arg] = supports;
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setSupports(
  const InArgs<Scalar>& inArgs, const int Np
  )
{
  std::copy(
    &inArgs.supports_[0],
    &inArgs.supports_[0] + NUM_E_IN_ARGS_MEMBERS, &supports_[0] );
  this->_set_Np( Np >= 0 ? Np : inArgs.Np() );
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setUnsupportsAndRelated(
  EInArgsMembers arg
  )
{
  switch(arg) {
    case IN_ARG_x: {
      this->_setSupports(IN_ARG_x_dot,false);
      this->_setSupports(IN_ARG_x_dot_poly,false);
      this->_setSupports(IN_ARG_alpha,false);
      this->_setSupports(IN_ARG_beta,false);
      break;
    }
#ifdef TEUCHOS_DEBUG
    default:
      TEST_FOR_EXCEPTION(
        true ,std::logic_error
        ,"Error, can handle args other than IN_ARG_x yet!"
        );
#endif
  }
  this->_setSupports(arg,false);
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::assert_supports(
  EInArgsMembers arg
  ) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"Thyra::ModelEvaluatorBase::InArgs<"
    << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(arg): "
    "model = \'"<<modelEvalDescription_<<"\': Error, "
    "The argument arg = " << toString(arg) << " is not supported!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::assert_l(int l) const
{
  TEST_FOR_EXCEPTION(
    !( 0 <= l && l < Np() ), std::logic_error
    ,"Thyra::ModelEvaluatorBase::InArgs<Scalar>::assert_l(l):\n\n"
    " model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The parameter l = " << l << " is not in the range [0,"<<Np()<<")!"
    );
}


//
// ModelEvaluatorBase::DerivativeMultiVector
//


template<class Scalar>
std::string ModelEvaluatorBase::DerivativeMultiVector<Scalar>::description() const
{
  using std::endl;
  std::ostringstream oss;
  oss << "DerivativeMultiVector{";
  if (is_null(getMultiVector())) {
    oss << "NULL";
  }
  else {
    oss
      << "multiVec=" << getMultiVector()->description()
      << ",orientation=" << toString(getOrientation());
  }
  oss << "}";
  return oss.str();
}


template<class Scalar>
void ModelEvaluatorBase::DerivativeMultiVector<Scalar>::describe(
  Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using std::endl;
  using Teuchos::describe;
  Teuchos::OSTab tab1(out);
  out << "DerivativeMultiVector\n";
  Teuchos::OSTab tab2(out);
  out
    << "multiVec = "
    << describe(*getMultiVector(),verbLevel)
    << "orientation = "
    << toString(getOrientation()) << endl;
}


// 2007/06/12: rabartl: The above description() and describe(...) functions
// have to be defined here and not in the class DerivativeMultiVector since it
// relies on the non-member function
// toString(ModelEvaluatorBase::EDerivativeMultiVectorOrientation) which is
// defined after the class definition for ModelEvaluatorBase.  This was caught
// by the intel compiler.  I am not sure why this worked with gcc.


//
// ModelEvaluatorBase::Derivative
//


template<class Scalar>
std::string
ModelEvaluatorBase::Derivative<Scalar>::description() const
{
  using std::endl;
  std::ostringstream oss;
  oss << "Derivative{";
  if (isEmpty()) {
    oss << "NULL";
  }
  else if (!is_null(getLinearOp())) {
    oss << "linearOp=" << getLinearOp()->description();
  }
  else {
    oss << "derivMultiVec=" << getDerivativeMultiVector().description();
  }
  oss << "}";
  return oss.str();
}


template<class Scalar>
void ModelEvaluatorBase::Derivative<Scalar>::describe( 
  Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using std::endl;
  using Teuchos::describe;
  Teuchos::OSTab tab1(out);
  out << "Derivative:";
  if (isEmpty()) {
    out << " NULL\n";
  }
  else if (!is_null(getLinearOp())) {
    out
      << endl
      << "linearOp = " << describe(*getLinearOp(),verbLevel);
  }
  else {
    out
      << endl
      << "derivMultiVec = ";
    getDerivativeMultiVector().describe(out,verbLevel);
  }
}


//
// ModelEvaluatorBase::OutArgs
//


template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::OutArgs()
  :modelEvalDescription_("WARNING!  THIS OUTARGS OBJECT IS UNINITALIZED!"),
   isFailed_(false)
{ std::fill_n(&supports_[0],NUM_E_OUT_ARGS_MEMBERS,false); }


template<class Scalar>
int ModelEvaluatorBase::OutArgs<Scalar>::Np() const
{ return DfDp_.size(); }


template<class Scalar>
int ModelEvaluatorBase::OutArgs<Scalar>::Ng() const
{ return g_.size(); }


template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsMembers arg
  ) const
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_
    <<"\': Error, arg="<<toString(arg)<<" is invalid!"
    );
  return supports_[arg];
}


template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsDfDp arg, int l
  ) const
{
  assert_l(l);
  return supports_DfDp_[l];
}


template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsDgDx_dot arg, int j
  ) const
{
  assert_j(j);
  return supports_DgDx_dot_[j];
}


template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsDgDx arg, int j
  ) const
{
  assert_j(j);
  return supports_DgDx_[j];
}


template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsDgDp arg, int j, int l
  ) const
{
  assert_j(j);
  assert_l(l);
  return supports_DgDp_[ j*Np() + l ];
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_f( 
  const RCP<VectorBase<Scalar> > &f
  )
{
  assert_supports(OUT_ARG_f);
  f_ = f;
}


template<class Scalar>
RCP<VectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_f() const
{
  assert_supports(OUT_ARG_f);
  return f_;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_g(
  int j, const RCP<VectorBase<Scalar> > &g_j
  )
{
  assert_j(j);
  g_[j] = g_j;
}


template<class Scalar>
RCP<VectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_g(int j) const
{ 
  assert_j(j);
  return g_[j];
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_W(
  const RCP<LinearOpWithSolveBase<Scalar> > &W
  )
{
  assert_supports(OUT_ARG_W);
  W_ = W;
}


template<class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_W() const
{
  assert_supports(OUT_ARG_W);
  return W_;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_W_op(
  const RCP<LinearOpBase<Scalar> > &W_op
  )
{
  assert_supports(OUT_ARG_W_op);
  W_op_ = W_op;
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_W_op() const
{
  assert_supports(OUT_ARG_W_op);
  return W_op_;
}


template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_W_properties() const
{
  assert_supports(OUT_ARG_f);
  return W_properties_;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DfDp(
  int l, const Derivative<Scalar> &DfDp_l
  )
{
  assert_supports(OUT_ARG_DfDp,l,DfDp_l);
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
void ModelEvaluatorBase::OutArgs<Scalar>::set_DgDx_dot(
  int j, const Derivative<Scalar> &DgDx_dot_j
  )
{
  assert_supports(OUT_ARG_DgDx_dot,j,DgDx_dot_j);
  DgDx_dot_[j] = DgDx_dot_j;
}


template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx_dot(int j) const
{
  assert_supports(OUT_ARG_DgDx_dot,j);
  return DgDx_dot_[j];
}


template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx_dot_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx_dot,j);
  return DgDx_dot_properties_[j];
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DgDx(
  int j, const Derivative<Scalar> &DgDx_j
  )
{
  assert_supports(OUT_ARG_DgDx,j,DgDx_j);
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
void ModelEvaluatorBase::OutArgs<Scalar>::set_DgDp(
  int j, int l, const Derivative<Scalar> &DgDp_j_l
  )
{
  assert_supports(OUT_ARG_DgDp,j,l,DgDp_j_l);
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
void ModelEvaluatorBase::OutArgs<Scalar>::set_f_poly(
  const RCP<Teuchos::Polynomial< VectorBase<Scalar> > > &f_poly
  )
{
  f_poly_ = f_poly;
}


template<class Scalar>
RCP<Teuchos::Polynomial< VectorBase<Scalar> > >
ModelEvaluatorBase::OutArgs<Scalar>::get_f_poly() const
{
  return f_poly_;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::setArgs(
  const OutArgs<Scalar>& inputOutArgs, bool ignoreUnsupported
  )
{
  typedef ModelEvaluatorBase MEB;
  const int min_Np = TEUCHOS_MIN(this->Np(),inputOutArgs.Np());
  const int min_Ng = TEUCHOS_MIN(this->Ng(),inputOutArgs.Ng());
  // f
  if ( inputOutArgs.supports(OUT_ARG_f) && inputOutArgs.get_f().get() ) {
    if ( supports(OUT_ARG_f) || !ignoreUnsupported )
      set_f(inputOutArgs.get_f());
  }
  // f_poly
  if ( inputOutArgs.supports(OUT_ARG_f_poly) && inputOutArgs.get_f_poly().get() ) {
    if ( supports(OUT_ARG_f_poly) || !ignoreUnsupported )
      set_f_poly(inputOutArgs.get_f_poly());
  }
  // g(j)
  for ( int j = 0; j < min_Ng; ++j ) {
    if ( inputOutArgs.get_g(j).get() )
      set_g(j,inputOutArgs.get_g(j));
  }
  // W
  if( inputOutArgs.supports(OUT_ARG_W) && inputOutArgs.get_W().get() ) {
    if ( supports(OUT_ARG_W) || !ignoreUnsupported )
      set_W(inputOutArgs.get_W());
  }
  // W_op
  if( inputOutArgs.supports(OUT_ARG_W_op) && inputOutArgs.get_W_op().get() ) {
    if ( supports(OUT_ARG_W_op) || !ignoreUnsupported )
      set_W_op(inputOutArgs.get_W_op());
  }
  // DfDp(l)
  for ( int l = 0; l < min_Np; ++l ) {
    MEB::Derivative<Scalar> DfDp_l;
    if ( !inputOutArgs.supports(OUT_ARG_DfDp,l).none()
      && !(DfDp_l=inputOutArgs.get_DfDp(l)).isEmpty() )
    {
      if ( DfDp_l.isSupportedBy(supports(OUT_ARG_DfDp,l)) || !ignoreUnsupported )
        set_DfDp(l,DfDp_l);
    }
  }
  // DgDx_dot(j) and DgDx(j)
  for ( int j = 0; j < min_Ng; ++j ) {
    // DgDx_dot(j)
    MEB::Derivative<Scalar> DgDx_dot_j;
    if ( !inputOutArgs.supports(OUT_ARG_DgDx_dot,j).none()
      && !(DgDx_dot_j=inputOutArgs.get_DgDx_dot(j)).isEmpty() )
    {
      if( DgDx_dot_j.isSupportedBy(supports(OUT_ARG_DgDx_dot,j)) || !ignoreUnsupported )
        set_DgDx_dot(j,DgDx_dot_j);
    }
    // DgDx(j)
    MEB::Derivative<Scalar> DgDx_j;
    if ( !inputOutArgs.supports(OUT_ARG_DgDx,j).none()
      && !(DgDx_j=inputOutArgs.get_DgDx(j)).isEmpty() ) {
      if ( DgDx_j.isSupportedBy(supports(OUT_ARG_DgDx,j)) || !ignoreUnsupported )
        set_DgDx(j,DgDx_j);
    }
  }
  // DgDp(j,l)
  for ( int l = 0; l < min_Np; ++l ) {
    for ( int j = 0; j < min_Ng; ++j ) {
      MEB::Derivative<Scalar> DgDp_j_l;
      if ( !inputOutArgs.supports(OUT_ARG_DgDp,j,l).none() 
        && !(DgDp_j_l=inputOutArgs.get_DgDp(j,l)).isEmpty() )
      {
        if ( DgDp_j_l.isSupportedBy(supports(OUT_ARG_DgDp,j,l)) || !ignoreUnsupported )
          set_DgDp(j,l,DgDp_j_l);
      }
    }
  }
  // ToDo: Add more args as needed!
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
bool ModelEvaluatorBase::OutArgs<Scalar>::isEmpty() const
{
  if (!is_null(f_))
    return false;
  if (!is_null(W_))
    return false;
  if (!is_null(W_op_))
    return false;
  for ( int l = 0; l < Np(); ++l ) {
    if (!DfDp_[l].isEmpty())
      return false;
  }
  if (!is_null(f_poly_))
    return false;
  for ( int j = 0; j < Ng(); ++j ) {
    if (!is_null(g_[j]))
      return false;
    if (!DgDx_dot_[j].isEmpty())
      return false;
    if (!DgDx_[j].isEmpty())
      return false;
    for ( int l = 0; l < Np(); ++l ) {
      if (!DgDp_[j*Np()+l].isEmpty())
        return false;
    }
  }
  return true;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assertSameSupport(
  const OutArgs<Scalar> &outArgs
  ) const
{

  for ( int outArg_i = 0; outArg_i < NUM_E_OUT_ARGS_MEMBERS; ++outArg_i ) {
    const EOutArgsMembers outArg_arg = static_cast<EOutArgsMembers>(outArg_i);
    const std::string outArg_name = toString(outArg_arg);
    TEST_FOR_EXCEPTION(
      supports(outArg_arg) != outArgs.supports(outArg_arg), std::logic_error,
      "Error, the output argument "<<outArg_name<<" with support "<<outArgs.supports(outArg_arg)<<"\n"
      "in the OutArgs object for the model:\n\n"
      "  "<<outArgs.modelEvalDescription()<<"\n\n"
      "is not the same the argument "<<outArg_name<<" with support "<<supports(outArg_arg)<<"\n"
      "in the OutArgs object for the model:\n\n"
      "  "<<modelEvalDescription()<<"\n\n"
      "and these two OutArgs objects are not compatible!"
      );
  }

  const int Np = this->Np();
  const int Ng = this->Ng();
  TEUCHOS_ASSERT_EQUALITY( Np, outArgs.Np() );
  TEUCHOS_ASSERT_EQUALITY( Ng, outArgs.Ng() );

  if (supports(OUT_ARG_f)) {
    for ( int l = 0; l < Np; ++l ) {
      TEST_FOR_EXCEPTION(
        !supports(OUT_ARG_DfDp,l).isSameSupport(outArgs.supports(OUT_ARG_DfDp,l)),
        std::logic_error,
        "Error, the support for DfDp("<<l<<") is not the same for the models\n\n"
        "  "<<outArgs.modelEvalDescription()<<"\n\n"
        "and:\n\n"
        "  "<<modelEvalDescription()<<"\n\n"
        "and these two OutArgs objects are not compatible!"
        );
    }
  }

  for ( int j = 0; j < Ng; ++j ) {
    TEST_FOR_EXCEPTION(
      !supports(OUT_ARG_DgDx_dot,j).isSameSupport(outArgs.supports(OUT_ARG_DgDx_dot,j)),
      std::logic_error,
      "Error, the support for DgDx_dot("<<j<<") is not the same for the models\n\n"
      "  "<<outArgs.modelEvalDescription()<<"\n\n"
      "and:\n\n"
      "  "<<modelEvalDescription()<<"\n\n"
      "and these two OutArgs objects are not compatible!"
      );
    TEST_FOR_EXCEPTION(
      !supports(OUT_ARG_DgDx,j).isSameSupport(outArgs.supports(OUT_ARG_DgDx,j)),
      std::logic_error,
      "Error, the support for DgDx("<<j<<") is not the same for the models\n\n"
      "  "<<outArgs.modelEvalDescription()<<"\n\n"
      "and:\n\n"
      "  "<<modelEvalDescription()<<"\n\n"
      "and these two OutArgs objects are not compatible!"
      );
    for ( int l = 0; l < Np; ++l ) {
      TEST_FOR_EXCEPTION(
        !supports(OUT_ARG_DgDp,j,l).isSameSupport(outArgs.supports(OUT_ARG_DgDp,j,l)),
        std::logic_error,
        "Error, the support for DgDp("<<j<<","<<l<<") is not the same for the models\n\n"
        "  "<<outArgs.modelEvalDescription()<<"\n\n"
        "and:\n\n"
        "  "<<modelEvalDescription()<<"\n\n"
        "and these two OutArgs objects are not compatible!"
        );
    }
  }
}


template<class Scalar>
std::string ModelEvaluatorBase::OutArgs<Scalar>::modelEvalDescription() const
{
  return modelEvalDescription_;
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
  using Teuchos::describe;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef RCP<const VectorBase<Scalar> > CV_ptr;
  typedef RCP<const LinearOpBase<Scalar> > CLO_ptr;
  typedef RCP<const LinearOpWithSolveBase<Scalar> > CLOWS_ptr;
  typedef ModelEvaluatorBase MEB;
  typedef MEB::Derivative<Scalar> Deriv;

  if( verbLevel == Teuchos::VERB_NONE && verbLevel == Teuchos::VERB_DEFAULT )
    return;

  RCP<Teuchos::FancyOStream>
    out = Teuchos::rcp(&out_arg,false);
  OSTab tab(out);

  *out <<"Thyra::ModelEvaluatorBase::OutArgs<"<<ST::name()<<">:\n";
  tab.incrTab();

  *out <<"model = " << modelEvalDescription_ << "\n";
  *out <<"Np = " << Np() << "\n";
  *out <<"Ng = " << Ng() << "\n";

  CV_ptr f;
  if (this->supports(OUT_ARG_f) && !is_null(f=get_f()) ) {
    *out << "f = " << Teuchos::describe(*f,verbLevel);
  }
  
  for( int j = 0; j < Ng(); ++j ) {
    CV_ptr g_j;
    if (!is_null(g_j=this->get_g(j)))
      *out << "g("<<j<<") = " << Teuchos::describe(*g_j,verbLevel);
  }
  
  CLOWS_ptr W;
  if ( this->supports(OUT_ARG_W) && !is_null(W=get_W()) ) {
    *out << "W = " << Teuchos::describe(*W,verbLevel);
  }
  
  CLO_ptr W_op;
  if ( this->supports(OUT_ARG_W_op) && !is_null(W_op=get_W_op()) ) {
    *out << "W_op = " << Teuchos::describe(*W_op,verbLevel);
  }
  
  for( int l = 0; l < Np(); ++l ) {
    Deriv DfDp_l;
    if (
      !this->supports(OUT_ARG_DfDp,l).none()
      && !(DfDp_l=get_DfDp(l)).isEmpty()
      )
    {
      *out << "DfDp("<<l<<") = ";
      DfDp_l.describe(*out,verbLevel);
    }
  }
  
  for( int j = 0; j < Ng(); ++j ) {
    
    Deriv DgDx_dot_j;
    if (
      !this->supports(OUT_ARG_DgDx_dot,j).none()
      && !(DgDx_dot_j=get_DgDx_dot(j)).isEmpty()
      )
    {
      *out << "DgDx_dot("<<j<<") = ";
      DgDx_dot_j.describe(*out,verbLevel);
    }
    
    Deriv DgDx_j;
    if (
      !this->supports(OUT_ARG_DgDx,j).none()
      && !(DgDx_j=get_DgDx(j)).isEmpty()
      )
    {
      *out << "DgDx("<<j<<") = ";
      DgDx_j.describe(*out,verbLevel);
    }
    
    for( int l = 0; l < Np(); ++l ) {
      
      Deriv DgDp_j_l;
      if (
        !this->supports(OUT_ARG_DgDp,j,l).none()
        && !(DgDp_j_l=get_DgDp(j,l)).isEmpty()
        )
      {
        *out << "DgDp("<<j<<","<<l<<") = ";
        DgDp_j_l.describe(*out,verbLevel);
      }
    }
    
  }
  
  // ToDo: Add output for more objects?

}


// protected


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setModelEvalDescription(
  const std::string &modelEvalDescription
  )
{ modelEvalDescription_ = modelEvalDescription; }

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_Np_Ng(int Np, int Ng)
{
  if(Np) {
    supports_DfDp_.resize(Np);
    DfDp_.resize(Np); std::fill_n(DfDp_.begin(),Np,Derivative<Scalar>());
    DfDp_properties_.resize(Np); std::fill_n(DfDp_properties_.begin(),Np,DerivativeProperties());
  }
  if(Ng) {
    g_.resize(Ng); std::fill_n(g_.begin(),Ng,Teuchos::null);
    supports_DgDx_dot_.resize(Ng);
    DgDx_dot_.resize(Ng); std::fill_n(DgDx_dot_.begin(),Ng,Derivative<Scalar>());
    DgDx_dot_properties_.resize(Ng); std::fill_n(DgDx_dot_properties_.begin(),Ng,DerivativeProperties());
    supports_DgDx_.resize(Ng);
    DgDx_.resize(Ng); std::fill_n(DgDx_.begin(),Ng,Derivative<Scalar>());
    DgDx_properties_.resize(Ng); std::fill_n(DgDx_properties_.begin(),Ng,DerivativeProperties());
  }
  if(Np && Ng) {
    const int NpNg = Np*Ng;
    supports_DgDp_.resize(NpNg);
    DgDp_.resize(NpNg); std::fill_n(DgDp_.begin(),NpNg,Derivative<Scalar>());
    DgDp_properties_.resize(NpNg); std::fill_n(DgDp_properties_.begin(),NpNg,DerivativeProperties());
  }
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsMembers arg, bool supports )
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_
    <<"\': Error, arg="<<toString(arg)<<" is invalid!"
    );
  supports_[arg] = supports;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsDfDp arg, int l, const DerivativeSupport& supports
  )
{
  assert_supports(OUT_ARG_f);
  assert_l(l);
  supports_DfDp_[l] = supports;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsDgDx_dot arg, int j, const DerivativeSupport& supports
  )
{
  assert_j(j);
  supports_DgDx_dot_[j] = supports;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsDgDx arg, int j, const DerivativeSupport& supports
  )
{
  assert_j(j);
  supports_DgDx_[j] = supports;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsDgDp arg, int j, int l, const DerivativeSupport& supports
  )
{
  assert_j(j);
  assert_l(l);
  supports_DgDp_[ j*Np()+ l ] = supports;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_W_properties( 
  const DerivativeProperties &properties
  )
{
  W_properties_ = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DfDp_properties(
  int l, const DerivativeProperties &properties
  )
{
  assert_supports(OUT_ARG_DfDp,l);
  DfDp_properties_[l] = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DgDx_dot_properties(
  int j, const DerivativeProperties &properties
  )
{
  assert_supports(OUT_ARG_DgDx_dot,j);
  DgDx_dot_properties_[j] = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DgDx_properties(
  int j, const DerivativeProperties &properties
  )
{
  assert_supports(OUT_ARG_DgDx,j);
  DgDx_properties_[j] = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DgDp_properties(
  int j, int l, const DerivativeProperties &properties
  )
{
  assert_supports(OUT_ARG_DgDp,j,l);
  DgDp_properties_[ j*Np()+ l ] = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  const OutArgs<Scalar>& inputOutArgs
  )
{
  typedef ModelEvaluatorBase MEB;
  const int Np = TEUCHOS_MIN(this->Np(),inputOutArgs.Np()); 
  const int Ng = TEUCHOS_MIN(this->Ng(),inputOutArgs.Ng()); 
  std::copy(
    &inputOutArgs.supports_[0],
    &inputOutArgs.supports_[0] + NUM_E_OUT_ARGS_MEMBERS, &supports_[0] );
  for( int l = 0; l < Np; ++l ) {
    DerivativeSupport ds = inputOutArgs.supports(MEB::OUT_ARG_DfDp,l);
    if (!ds.none()) {
      this->_setSupports(MEB::OUT_ARG_DfDp,l,ds);
      this->_set_DfDp_properties(l,inputOutArgs.get_DfDp_properties(l));
    }
  }
  for( int j = 0; j < Ng; ++j ) {
    DerivativeSupport ds = inputOutArgs.supports(MEB::OUT_ARG_DgDx_dot,j);
    this->_setSupports(MEB::OUT_ARG_DgDx_dot,j,ds);
    if(!ds.none()) this->_set_DgDx_dot_properties(j,inputOutArgs.get_DgDx_dot_properties(j));
  }
  for( int j = 0; j < Ng; ++j ) {
    DerivativeSupport ds = inputOutArgs.supports(MEB::OUT_ARG_DgDx,j);
    this->_setSupports(MEB::OUT_ARG_DgDx,j,ds);
    if(!ds.none()) this->_set_DgDx_properties(j,inputOutArgs.get_DgDx_properties(j));
  }
  for( int j = 0; j < Ng; ++j ) for( int l = 0; l < Np; ++l ) {
    DerivativeSupport ds = inputOutArgs.supports(MEB::OUT_ARG_DgDp,j,l);
    this->_setSupports(MEB::OUT_ARG_DgDp,j,l,ds);
    if(!ds.none()) this->_set_DgDp_properties(j,l,inputOutArgs.get_DgDp_properties(j,l));
  }
  if(this->supports(OUT_ARG_W) || this->supports(OUT_ARG_W_op))
    this->_set_W_properties(inputOutArgs.get_W_properties());
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setUnsupportsAndRelated(
  EInArgsMembers arg
  )
{
  switch(arg) {
    case IN_ARG_x: {
      const int Ng = this->Ng();
      for( int j = 0; j < Ng; ++j ) {
        this->_setSupports(OUT_ARG_DgDx_dot,j,DerivativeSupport());
        this->_setSupports(OUT_ARG_DgDx,j,DerivativeSupport());
      }
      break;
    }
#ifdef TEUCHOS_DEBUG
    default:
      TEST_FOR_EXCEPTION(
        true ,std::logic_error
        ,"Error, can handle args other than IN_ARG_x yet!"
        );
#endif
  }
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setUnsupportsAndRelated(
  EOutArgsMembers arg
  )
{
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
#ifdef TEUCHOS_DEBUG
    default:
      TEST_FOR_EXCEPTION(
        true ,std::logic_error
        ,"Error, can handle args other than OUT_ARG_f yet!"
        );
#endif
  }
  this->_setSupports(arg,false);
}


// private


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(EOutArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !this->supports(arg), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(arg):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument arg = " << toString(arg) << " is not supported!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgsDfDp arg, int l, const Derivative<Scalar> &deriv
  ) const
{
  const DerivativeSupport derivSupport = this->supports(arg,l);
  TEST_FOR_EXCEPTION(
    !deriv.isSupportedBy(derivSupport), std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_DfDp,l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument DfDp("<<l<<") = " << deriv.description() << "\n"
    "is not supported!\n\n"
    "The supported types include " << derivSupport.description() << "!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgsDgDx_dot arg, int j, const Derivative<Scalar> &deriv
  ) const
{
  const DerivativeSupport derivSupport = this->supports(arg,j);
  TEST_FOR_EXCEPTION(
    !deriv.isSupportedBy(derivSupport), std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_DgDx_dot,j):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument DgDx_dot("<<j<<") = " << deriv.description() << "\n"
    "is not supported!\n\n"
    "The supported types include " << derivSupport.description() << "!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgsDgDx arg, int j, const Derivative<Scalar> &deriv
  ) const
{
  const DerivativeSupport derivSupport = this->supports(arg,j);
  TEST_FOR_EXCEPTION(
    !deriv.isSupportedBy(derivSupport), std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_DgDx,j):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument DgDx("<<j<<") = " << deriv.description() << "\n"
    "is not supported!\n\n"
    "The supported types include " << derivSupport.description() << "!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgsDgDp arg, int j, int l, const Derivative<Scalar> &deriv
  ) const
{
  const DerivativeSupport derivSupport = this->supports(arg,j,l);
  TEST_FOR_EXCEPTION(
    !deriv.isSupportedBy(derivSupport), std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_DgDp,j,l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument DgDp("<<j<<","<<l<<") = " << deriv.description() << "\n"
    "is not supported!\n\n"
    "The supported types include " << derivSupport.description() << "!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_l(int l) const
{
  TEST_FOR_EXCEPTION(
    !( 0 <= l && l < Np() ), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_l(l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error,  The parameter subvector p("<<l<<")"
    " is not in the range [0,"<<Np()<<")!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_j(int j) const
{
  TEST_FOR_EXCEPTION(
    !( 0 <= j && j < Ng() ), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_j(j):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The auxiliary function g("<<j<<")"
    " is not in the range [0,"<<Ng()<<")!"
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
void ModelEvaluatorBase::InArgsSetup<Scalar>::setSupports(
  const InArgs<Scalar>& inArgs, const int Np
  )
{ this->_setSupports(inArgs,Np); }


template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::setUnsupportsAndRelated(
  EInArgsMembers arg
  )
{ this->_setUnsupportsAndRelated(arg); }


//
// ModelEvaluatorBase::OutArgsSetup
//


template<class Scalar>
ModelEvaluatorBase::OutArgsSetup<Scalar>::OutArgsSetup()
{}


template<class Scalar>
ModelEvaluatorBase::OutArgsSetup<Scalar>::OutArgsSetup(
  const OutArgs<Scalar>& inputOutArgs
  )
  :OutArgs<Scalar>(inputOutArgs)
{}


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setModelEvalDescription(
  const std::string &modelEvalDescription
  )
{ this->_setModelEvalDescription(modelEvalDescription); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_Np_Ng(int Np, int Ng)
{ this->_set_Np_Ng(Np,Ng); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgsMembers arg, bool supports
  )
{ this->_setSupports(arg,supports); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgsDfDp arg, int l, const DerivativeSupport& supports
  )
{ this->_setSupports(arg,l,supports); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports( 
  EOutArgsDgDx_dot arg, int j, const DerivativeSupport& supports
  )
{ this->_setSupports(arg,j,supports); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgsDgDx arg, int j, const DerivativeSupport& supports
  )
{ this->_setSupports(arg,j,supports); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgsDgDp arg, int j, int l, const DerivativeSupport& supports
  )
{ this->_setSupports(arg,j,l,supports); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_W_properties(
  const DerivativeProperties &properties
  )
{ this->_set_W_properties(properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DfDp_properties(
  int l, const DerivativeProperties &properties
  )
{ this->_set_DfDp_properties(l,properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DgDx_dot_properties(
  int j, const DerivativeProperties &properties
  )
{ this->_set_DgDx_dot_properties(j,properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DgDx_properties(
  int j, const DerivativeProperties &properties
  )
{ this->_set_DgDx_properties(j,properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DgDp_properties(
  int j, int l, const DerivativeProperties &properties
  )
{ this->_set_DgDp_properties(j,l,properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  const OutArgs<Scalar>& inputOutArgs
  )
{ this->_setSupports(inputOutArgs); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setUnsupportsAndRelated(
  EInArgsMembers arg
  )
{ this->_setUnsupportsAndRelated(arg); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setUnsupportsAndRelated(
  EOutArgsMembers arg
  )
{ this->_setUnsupportsAndRelated(arg); }


} // namespace Thyra


#endif // THYRA_MODEL_EVALUATOR_BASE_HPP
