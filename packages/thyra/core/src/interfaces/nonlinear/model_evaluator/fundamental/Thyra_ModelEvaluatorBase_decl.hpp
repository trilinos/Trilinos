// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_MODEL_EVALUATOR_BASE_DECL_HPP
#define THYRA_MODEL_EVALUATOR_BASE_DECL_HPP


#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_any.hpp"
#include "Teuchos_Assert.hpp"

#ifdef HAVE_THYRA_ME_POLYNOMIAL
#  include "Teuchos_Polynomial.hpp"
#endif

namespace Stokhos {
  class ProductEpetraVector;
  class ProductEpetraMultiVector;
  class ProductEpetraOperator;
}

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
    IN_ARG_x_dot_dot ///< .
    ,IN_ARG_x_dot ///< .
    ,IN_ARG_x ///< .
    ,IN_ARG_x_dot_poly ///< .
    ,IN_ARG_x_poly ///< .
    ,IN_ARG_x_dot_mp ///< .
    ,IN_ARG_x_mp ///< .
    ,IN_ARG_t ///< .
    ,IN_ARG_alpha ///< .
    ,IN_ARG_beta ///< .
    ,IN_ARG_W_x_dot_dot_coeff ///< .
    ,IN_ARG_step_size///< .
    ,IN_ARG_stage_number///< .
  };
  /** \brief .  */
  static const int NUM_E_IN_ARGS_MEMBERS=13;

  /** \brief .  */
  enum EInArgs_p_mp {
    IN_ARG_p_mp ///< .
  };

  /** \brief Concrete aggregate class for all input arguments computable by a
   * <tt>ModelEvaluator</tt> subclass object.
   *
   * The set of supported objects is returned from the <tt>supports()</tt>
   * function.
   *
   * A client can not directly set what input arguments are supported or not
   * supported.  Only a subclass of <tt>ModelEvaluator</tt> can do that
   * (through the <tt>InArgsSetup</tt> subclass).
   *
   * Extended types: The methods for %supports(), set() and get() that
   * are templated on the <tt>ObjectType</tt> are used to support what
   * are called extended types. This functionality allows developers
   * to inject new objects into the ModelEvaluator evaluations. This
   * can be used for expermenting with new capabilities that require
   * adding new/special arguments to the inArgs. This also can be used
   * to reduce the clutter and complexity of the model evaluator
   * inArgs object. For example, users could create their own inArgs
   * for a stochastic computation:
   *
   * \code{.cpp}
   *   struct StochasticInArgs {
   *     StochLib::MultiVector<Scalar> stochastic_x;
   *     ...
   *   };
   * \endcode
   *
   * This object could then be used in a model evaluator without
   * changing the base classes in Thyra:
   *
   * \code{.cpp}
   * auto inArgs = model->createInArgs();
   * assert(inArgs.template supports<StochasticInArgs<Scalar>>());
   *
   * RCP<StochasticInArgs> stochInArgs = rcp(new StochasticInArgs<Scalar>);
   * 
   * stochInArgs->stochastic_x = ... // create and set objects
   * ...
   *
   * inArgs.set(stochInArgs);
   * \endcode
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
    /** \brief Return the number of axillary response functions
     * <tt>g(j)(...)</tt> supported (<tt>Ng >= 0</tt>).  */
    int Ng() const;
    /** \brief Determines if an input argument is supported or not.  */
    bool supports(EInArgsMembers arg) const;
    /** \brief Precondition: <tt>supports(IN_ARG_x_dot_dot)==true</tt>.  */
    void set_x_dot_dot( const RCP<const VectorBase<Scalar> > &x_dot_dot );
    /** \brief Precondition: <tt>supports(IN_ARG_x_dot_dot)==true</tt>.  */
    RCP<const VectorBase<Scalar> > get_x_dot_dot() const;
    /** \brief Precondition: <tt>supports(IN_ARG_x_dot)==true</tt>.  */
    void set_x_dot( const RCP<const VectorBase<Scalar> > &x_dot );
    /** \brief Precondition: <tt>supports(IN_ARG_x_dot)==true</tt>.  */
    RCP<const VectorBase<Scalar> > get_x_dot() const;
    /** \brief Precondition: <tt>supports(IN_ARG_x)==true</tt>.  */
    void set_x( const RCP<const VectorBase<Scalar> > &x );
    /** \brief Precondition: <tt>supports(IN_ARG_x)==true</tt>.  */
    RCP<const VectorBase<Scalar> > get_x() const;

    /** \brief Precondition: <tt>supports(IN_ARG_x)==true</tt>.  */
    void set_x_direction( const RCP<const MultiVectorBase<Scalar> > &x_direction );
    /** \brief Precondition: <tt>supports(IN_ARG_x)==true</tt>.  */
    void set_p_direction( int l, const RCP<const MultiVectorBase<Scalar> > &p_direction_l );
    /** \brief Precondition: <tt>supports(IN_ARG_x)==true</tt>.  */
    RCP<const MultiVectorBase<Scalar> > get_x_direction() const;
    /** \brief Get <tt>p(l)</tt> where <tt>0 <= l && l < this->Np()</tt>.  */
    RCP<const MultiVectorBase<Scalar> > get_p_direction(int l) const;

    /** \brief Precondition: <tt>supports(IN_ARG_x)==true</tt>.  */
    void set_f_multiplier( const RCP<const VectorBase<Scalar> > &f_multiplier );
    /** \brief Precondition: <tt>supports(IN_ARG_x)==true</tt>.  */
    RCP<const VectorBase<Scalar> > get_f_multiplier() const;
    /** \brief Precondition: <tt>supports(IN_ARG_x)==true</tt>.  */
    void set_g_multiplier( int j, const RCP<const VectorBase<Scalar> > &g_multiplier );
    /** \brief Precondition: <tt>supports(IN_ARG_x)==true</tt>.  */
    RCP<const VectorBase<Scalar> > get_g_multiplier(int j) const;

    /** \brief Determines if an extended input argument of type <tt>ObjectType</tt> is supported. */
    template<typename ObjectType>
    bool supports() const;
    /** \brief Set an extended input object of type <tt>ObjectType</tt> in the InArgs. Precondition: <tt>supports()==true</tt>. */
    template<typename ObjectType>
    void set(const RCP<const ObjectType>& uo);
    /** \brief Get an extended input object of type <tt>ObjectType>/tt> from the InArgs. Precondition: <tt>supports()==true</tt>. */
    template<typename ObjectType>
    RCP<const ObjectType> get() const;
    
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


    /** \brief Precondition: <tt>supports(IN_ARG_x_dot_mp)==true</tt>.  */
    void set_x_dot_mp( const RCP<const Stokhos::ProductEpetraVector > &x_dot_mp );
    /** \brief Precondition: <tt>supports(IN_ARG_x_dotmp)==true</tt>.  */
    RCP<const Stokhos::ProductEpetraVector > get_x_dot_mp() const;

    /** \brief Precondition: <tt>supports(IN_ARG_x_mp)==true</tt>.  */
    void set_x_mp( const RCP<const Stokhos::ProductEpetraVector > &x_mp );
    /** \brief Precondition: <tt>supports(IN_ARG_x_mp)==true</tt>.  */
    RCP<const Stokhos::ProductEpetraVector > get_x_mp() const;

    void set_p_mp( int l, const RCP<const Stokhos::ProductEpetraVector > &p_mp_l );
    RCP<const Stokhos::ProductEpetraVector > get_p_mp(int l) const;
    /** Whether p_mp is supported for parameter vector l */
    bool supports(EInArgs_p_mp arg, int l) const;


    /** \brief Precondition: <tt>supports(IN_ARG_t)==true</tt>.  */
    void set_t( ScalarMag t );
    /** \brief .Precondition: <tt>supports(IN_ARG_t)==true</tt>  */
    ScalarMag get_t() const;
    /** \brief Precondition: <tt>supports(IN_ARG_alpha)==true</tt>.  */
    void set_alpha( Scalar alpha );
    /** \brief Precondition: <tt>supports(IN_ARG_alpha)==true</tt>.  */
    Scalar get_alpha() const;
    /** \brief Precondition: <tt>supports(IN_ARG_beta)==true</tt>.  */
    void set_beta( Scalar beta );
    /** \brief Precondition: <tt>supports(IN_ARG_beta)==true</tt>.  */
    Scalar get_beta() const;
    /** \brief Precondition: <tt>supports(IN_ARG_W_x_dot_dot_coeff)==true</tt>.  */
    void set_W_x_dot_dot_coeff( Scalar W_x_dot_dot_coeff );
    /** \brief Precondition: <tt>supports(IN_ARG_W_x_dot_dot_coeff)==true</tt>.  */
    Scalar get_W_x_dot_dot_coeff() const;
    /** \brief Precondition: <tt>supports(IN_ARG_step_size)==true</tt>.  */
    void set_step_size( Scalar step_size);
    /** \brief Precondition: <tt>supports(IN_ARG_step_size)==true</tt>.  */
    Scalar get_step_size() const;
    /** \brief Precondition: <tt>supports(IN_ARG_stage_number)==true</tt>.  */
    void set_stage_number( Scalar stage_number);
    /** \brief Precondition: <tt>supports(IN_ARG_stage_number)==true</tt>.  */
    Scalar get_stage_number() const;
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
    void _set_Np_Ng(int Np, int Ng);
    /** \brief . */
    void _setSupports( EInArgsMembers arg, bool supports );
    /** \brief . */
    void _setSupports( EInArgs_p_mp arg, int l, bool supports );
    /** \brief . */
    void _setSupports( const InArgs<Scalar>& inputInArgs, const int Np );
    /** \brief . */
    template<typename ObjectType>
    void _setSupports( const bool supports );
    /** \brief . */
    void _setUnsupportsAndRelated( EInArgsMembers arg );
  private:
    // types
    typedef Teuchos::Array<RCP<const VectorBase<Scalar> > > p_t;
    typedef Teuchos::Array<RCP<const MultiVectorBase<Scalar> > > p_direction_t;
    // data
    std::string modelEvalDescription_;
    RCP<const VectorBase<Scalar> > x_dot_dot_;
    RCP<const VectorBase<Scalar> > x_dot_;
    RCP<const VectorBase<Scalar> > x_;
    RCP<const MultiVectorBase<Scalar> > x_direction_;
    RCP<const Stokhos::ProductEpetraVector > x_dot_mp_;
    RCP<const Stokhos::ProductEpetraVector > x_mp_;

    RCP<const VectorBase<Scalar> > f_multiplier_;
    p_t g_multiplier_;
    Teuchos::Array< RCP< const Stokhos::ProductEpetraVector > > p_mp_;
#ifdef HAVE_THYRA_ME_POLYNOMIAL
    RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > x_dot_poly_;
    RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > x_poly_;
#endif // HAVE_THYRA_ME_POLYNOMIAL
    p_t p_;
    p_direction_t p_direction_;
    ScalarMag t_;
    Scalar alpha_;
    Scalar beta_;
    Scalar W_x_dot_dot_coeff_;
    Scalar step_size_;
    Scalar stage_number_;
    bool supports_[NUM_E_IN_ARGS_MEMBERS];
    Teuchos::Array<bool> supports_p_mp_; //Np
    // functions
    void assert_supports(EInArgsMembers arg) const;
    void assert_supports(EInArgs_p_mp arg, int l) const;
    void assert_l(int l) const;
    void assert_j(int j) const;

    std::map<std::string,Teuchos::any> extended_inargs_;
  };

  /** \brief The type of an evaluation. */
  enum EEvalType {
    EVAL_TYPE_EXACT = 0, /// Do an exact evaluation (default)
    EVAL_TYPE_APPROX_DERIV, /// An approx. eval. for a F.D. deriv.
    EVAL_TYPE_VERY_APPROX_DERIV /// An approx. eval. for a F.D. prec.
  };

  /** \brief Type to embed evaluation accuracy with an RCP-managed object.
   *
   * This type derives from Teuchos::RCP and therefore is a drop in
   * replacement for an RCP object as it implicitly converts to an from such a
   * type.
   */
  template<class ObjType>
  class Evaluation : public RCP<ObjType> {
  public:
    /** \brief . */
    Evaluation(Teuchos::ENull)
      : evalType_(EVAL_TYPE_EXACT)  {}
    /** \brief . */
    Evaluation() : evalType_(EVAL_TYPE_EXACT) {}
    /** \brief Implicit conversion from RCP<ObjType>. */
    Evaluation( const RCP<ObjType> &obj )
      : RCP<ObjType>(obj), evalType_(EVAL_TYPE_EXACT) {}
    /** \brief . */
    Evaluation( const RCP<ObjType> &obj, EEvalType evalType )
      : RCP<ObjType>(obj), evalType_(evalType) {}
    /** \brief . */
    EEvalType getType() const { return evalType_; }
    /** \brief . */
    void reset( const RCP<ObjType> &obj, EEvalType evalType ) 
    { this->operator=(obj); evalType_ = evalType; }
   private:
    EEvalType evalType_;
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
          default: TEUCHOS_TEST_FOR_EXCEPT(true);
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
          default: TEUCHOS_TEST_FOR_EXCEPT(true);
        }
        TEUCHOS_UNREACHABLE_RETURN(false);
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

  /** \brief Simple aggregate class for a derivative object represented as a
   * column-wise multi-vector or its transpose as a row-wise multi-vector.
   */
  class MPDerivativeMultiVector {
  public:
    /** \brief . */
    MPDerivativeMultiVector()
      :orientation_(DERIV_MV_BY_COL)
      {}
    /** \brief . */
    MPDerivativeMultiVector(
      const RCP<Stokhos::ProductEpetraMultiVector > &mv
      ,const EDerivativeMultiVectorOrientation orientation = DERIV_MV_BY_COL
      ,const Teuchos::Array<int> &paramIndexes = Teuchos::Array<int>()
      ) : mv_(mv.assert_not_null()), orientation_(orientation), paramIndexes_(paramIndexes) {}
    /** \brief . */
    void changeOrientation( const EDerivativeMultiVectorOrientation orientation )
      { orientation_ = orientation; };
    /** \brief . */
    const MPDerivativeMultiVector& assert_not_null() const
      { mv_.assert_not_null(); return *this; }
    /** \brief . */
    RCP<Stokhos::ProductEpetraMultiVector > getMultiVector() const
      { return mv_; }
    /** \brief . */
    EDerivativeMultiVectorOrientation getOrientation() const
      { return orientation_; }
    /** \brief . */
    const Teuchos::Array<int>& getParamIndexes() const
      { return paramIndexes_; }
    /** \brief . */
    std::string description() const {return "\n";}
    /** \brief . */
    void describe( 
      Teuchos::FancyOStream &/* out */, const Teuchos::EVerbosityLevel /* verbLevel */
      ) const {}
  private:
    RCP<Stokhos::ProductEpetraMultiVector > mv_;
    EDerivativeMultiVectorOrientation orientation_;
    Teuchos::Array<int> paramIndexes_;
  };

  /** \brief Simple aggregate class that stores a derivative object
   * as a general linear operator or as a multi-vector.
   */
  class MPDerivative {
  public:
    /** \brief . */
    MPDerivative() {}
    /** \brief . */
    MPDerivative( const RCP<Stokhos::ProductEpetraOperator > &lo )
      : lo_(lo.assert_not_null()) {}
    /** \brief . */
    MPDerivative(
      const RCP<Stokhos::ProductEpetraMultiVector > &mv,
      const EDerivativeMultiVectorOrientation orientation = DERIV_MV_BY_COL
      ) : dmv_(mv,orientation) {}
    /** \brief . */
    MPDerivative( const MPDerivativeMultiVector &dmv )
      : dmv_(dmv) {}
    /** \brief . */
    bool isEmpty() const
      { return ( lo_.get()==NULL && dmv_.getMultiVector().get()==NULL ); }
    /** \brief . */
    const MPDerivative& assert_not_null() const
      { dmv_.assert_not_null(); lo_.assert_not_null(); return *this; }
    /** \brief . */
    RCP<Stokhos::ProductEpetraOperator > getLinearOp() const
      { return lo_; }
    /** \brief . */
    RCP<Stokhos::ProductEpetraMultiVector > getMultiVector() const
      { return dmv_.getMultiVector(); }
    /** \brief . */
    EDerivativeMultiVectorOrientation getMultiVectorOrientation() const
      { return dmv_.getOrientation(); }
    /** \brief . */
    MPDerivativeMultiVector getDerivativeMultiVector() const
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
    std::string description() const {return "\n";}
    /** \brief . */
    void describe( 
      Teuchos::FancyOStream &/* out */, const Teuchos::EVerbosityLevel /* verbLevel */
      ) const {}
  private:
    RCP<Stokhos::ProductEpetraOperator > lo_;
    MPDerivativeMultiVector dmv_;
  };

  /** \brief .  */
  enum EOutArgsMembers {
    OUT_ARG_f,  ///< .
    OUT_ARG_W,  ///< .
    OUT_ARG_f_mp,  ///< .
    OUT_ARG_W_mp,  ///< .
    OUT_ARG_W_op,  ///< .
    OUT_ARG_W_prec, ///< .
    OUT_ARG_f_poly  ///< .
  };
  /** \brief .  */
  static const int NUM_E_OUT_ARGS_MEMBERS=7;

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

#ifdef Thyra_BUILD_HESSIAN_SUPPORT
  /** \brief . */
  enum EOutArgs_hess_vec_prod_f_xx {
    OUT_ARG_hess_vec_prod_f_xx   ///< .
  };

  /** \brief . */
  enum EOutArgs_hess_vec_prod_f_xp {
    OUT_ARG_hess_vec_prod_f_xp   ///< .
  };

  /** \brief . */
  enum EOutArgs_hess_vec_prod_f_px {
    OUT_ARG_hess_vec_prod_f_px   ///< .
  };

  /** \brief . */
  enum EOutArgs_hess_vec_prod_f_pp {
    OUT_ARG_hess_vec_prod_f_pp   ///< .
  };

  /** \brief . */
  enum EOutArgs_hess_vec_prod_g_xx {
    OUT_ARG_hess_vec_prod_g_xx   ///< .
  };

  /** \brief . */
  enum EOutArgs_hess_vec_prod_g_xp {
    OUT_ARG_hess_vec_prod_g_xp   ///< .
  };

  /** \brief . */
  enum EOutArgs_hess_vec_prod_g_px {
    OUT_ARG_hess_vec_prod_g_px   ///< .
  };

  /** \brief . */
  enum EOutArgs_hess_vec_prod_g_pp {
    OUT_ARG_hess_vec_prod_g_pp   ///< .
  };

  /** \brief . */
  enum EOutArgs_hess_f_xx {
    OUT_ARG_hess_f_xx   ///< .
  };

  /** \brief . */
  enum EOutArgs_hess_f_xp {
    OUT_ARG_hess_f_xp   ///< .
  };

  /** \brief . */
  enum EOutArgs_hess_f_pp {
    OUT_ARG_hess_f_pp   ///< .
  };

  /** \brief . */
  enum EOutArgs_hess_g_xx {
    OUT_ARG_hess_g_xx   ///< .
  };

  /** \brief . */
  enum EOutArgs_hess_g_xp {
    OUT_ARG_hess_g_xp   ///< .
  };

  /** \brief . */
  enum EOutArgs_hess_g_pp {
    OUT_ARG_hess_g_pp   ///< .
  };

  /** \brief . */
  enum EOutArgs_H_xx {
    OUT_ARG_H_xx   ///< .
  };

  /** \brief . */
  enum EOutArgs_H_xp {
    OUT_ARG_H_xp   ///< .
  };

  /** \brief . */
  enum EOutArgs_H_pp {
    OUT_ARG_H_pp   ///< .
  };
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

  /** \brief . */
  enum EOutArgsDfDp_mp {
    OUT_ARG_DfDp_mp   ///< .
  };

  /** \brief . */
  enum EOutArgs_g_mp {
    OUT_ARG_g_mp   ///< .
  };

  /** \brief . */
  enum EOutArgsDgDx_dot_mp {
    OUT_ARG_DgDx_dot_mp   ///< .
  };

  /** \brief . */
  enum EOutArgsDgDx_mp {
    OUT_ARG_DgDx_mp   ///< .
  };

  /** \brief . */
  enum EOutArgsDgDp_mp {
    OUT_ARG_DgDp_mp   ///< .
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
   *
   * Extended types: The methods for %supports(), set() and get() that
   * are templated on the <tt>ObjectType</tt> are used to support what
   * are called extended types. This functionality allows developers
   * to inject new objects into the ModelEvaluator evaluations. This
   * can be used for expermenting with new capabilities that require
   * adding new/special arguments to the outArgs. This also can be
   * used to reduce the clutter and complexity of the model evaluator
   * outArgs object. For example, users could create their own outArgs
   * for a stochastic computation:
   *
   * \code{.cpp}
   *   struct StochasticOutArgs {
   *     StochLib::MultiVector<Scalar> stochastic_f;
   *     ...
   *   };
   * \endcode
   *
   * This object could then be used in a model evaluator without
   * changing the base classes in Thyra:
   *
   * \code{.cpp}
   * auto outArgs = model->createOutArgs();
   * assert(outArgs.template supports<StochasticOutArgs<Scalar>>());
   *
   * RCP<StochasticOutArgs> stochOutArgs = rcp(new StochasticOutArgs<Scalar>);
   * 
   * stochOutArgs->stochastic_f = ... // create and set objects
   * ...
   *
   * outArgs.set(stochOutArgs);
   * \endcode
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

#ifdef Thyra_BUILD_HESSIAN_SUPPORT
    /** \brief Determine if <tt>hess_vec_prod_f_xx</tt> is supported or not.  */
    bool supports(EOutArgs_hess_vec_prod_f_xx arg) const;
    /** \brief Determine if <tt>hess_vec_prod_f_xp(l)</tt> is supported or not, <tt>0 <= l 
     * && l < Np()</tt>.  */
    bool supports(EOutArgs_hess_vec_prod_f_xp arg, int l) const;
    /** \brief Determine if <tt>hess_vec_prod_f_px(l)</tt> is supported or not, <tt>0 <= l
     * && l < Np()</tt>.  */
    bool supports(EOutArgs_hess_vec_prod_f_px arg, int l) const;
    /** \brief Determine if <tt>hess_vec_prod_f_pp(l1,l2)</tt> is supported or not, <tt>0 <= l1
     * && l1 < Np()</tt> and <tt>0 <= l2 && l2 < Np()</tt>.  */
    bool supports(EOutArgs_hess_vec_prod_f_pp arg, int l1, int l2) const;
    /** \brief Determine if <tt>hess_vec_prod_g_xx(j,l)</tt> is supported or not, <tt>0 <= j
     * && j < Ng()</tt>.  */
    bool supports(EOutArgs_hess_vec_prod_g_xx arg, int j) const;
    /** \brief Determine if <tt>hess_vec_prod_g_xp(j,l)</tt> is supported or not, <tt>0 <= j
     * && j < Ng()</tt> and <tt>0 <= l && l < Np()</tt>.  */
    bool supports(EOutArgs_hess_vec_prod_g_xp arg, int j, int l) const;
    /** \brief Determine if <tt>hess_vec_prod_g_px(j,l)</tt> is supported or not, <tt>0 <= j
     * && j < Ng()</tt> and <tt>0 <= l && l < Np()</tt>.  */
    bool supports(EOutArgs_hess_vec_prod_g_px arg, int j, int l) const;
    /** \brief Determine if <tt>hess_vec_prod_g_pp(j,l1,l2)</tt> is supported or not, <tt>0 <= j
     * && j < Ng()</tt>, <tt>0 <= l1 && l1 < Np()</tt>, and <tt>0 <= l2 && l2 < Np()</tt>.  */
    bool supports(EOutArgs_hess_vec_prod_g_pp arg, int j, int l1, int l2) const;
    /** \brief Determine if <tt>hess_f_xx</tt> is supported or not.  */
    bool supports(EOutArgs_hess_f_xx arg) const;
    /** \brief Determine if <tt>hess_f_xp(l)</tt> is supported or not, <tt>0 <= l
     * && l < Np()</tt>.  */
    bool supports(EOutArgs_hess_f_xp arg, int l) const;
    /** \brief Determine if <tt>hess_f_pp(l1,l2)</tt> is supported or not, <tt>0 <= l1
     * && l1 < Np()</tt> and <tt>0 <= l2 && l2 < Np()</tt>.  */
    bool supports(EOutArgs_hess_f_pp arg, int l1, int l2) const;
    /** \brief Determine if <tt>hess_g_xx(j,l)</tt> is supported or not, <tt>0 <= j
     * && j < Ng()</tt>.  */
    bool supports(EOutArgs_hess_g_xx arg, int j) const;
    /** \brief Determine if <tt>hess_g_xp(j,l)</tt> is supported or not, <tt>0 <= j
     * && j < Ng()</tt> and <tt>0 <= l && l < Np()</tt>.  */
    bool supports(EOutArgs_hess_g_xp arg, int j, int l) const;
    /** \brief Determine if <tt>hess_g_pp(j,l1,l2)</tt> is supported or not, <tt>0 <= j
     * && j < Ng()</tt>, <tt>0 <= l1 && l1 < Np()</tt>, and <tt>0 <= l2 && l2 < Np()</tt>.  */
    bool supports(EOutArgs_hess_g_pp arg, int j, int l1, int l2) const;
    /** \brief Determine if <tt>H_xx</tt> is supported or not.  */
    bool supports(EOutArgs_H_xx arg) const;
    /** \brief Determine if <tt>H_xp(l)</tt> is supported or not, <tt>0 <= l
     * && l < Np()</tt>.  */
    bool supports(EOutArgs_H_xp arg, int l) const;
    /** \brief Determine if <tt>H_pp(l1,l2)</tt> is supported or not, <tt>0 <= l1
     * && l1 < Np()</tt> and <tt>0 <= l2 && l2 < Np()</tt>.  */
    bool supports(EOutArgs_H_pp arg, int l1, int l2) const;
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

    /** \brief Precondition: <tt>supports(OUT_ARG_f)==true</tt>.  */
    void set_f( const Evaluation<VectorBase<Scalar> > &f );
    /** \brief Precondition: <tt>supports(OUT_ARG_f)==true</tt>.  */
    Evaluation<VectorBase<Scalar> > get_f() const;
    /** \brief Precondition: <tt>supports(OUT_ARG_g)==true</tt>.  */
    void set_g( int j, const Evaluation<VectorBase<Scalar> > &g_j );
    /** \brief Precondition: <tt>supports(OUT_ARG_g)==true</tt>..  */
    Evaluation<VectorBase<Scalar> > get_g(int j) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_W)==true</tt>.  */
    void set_W( const RCP<LinearOpWithSolveBase<Scalar> > &W );
    /** \brief Precondition: <tt>supports(OUT_ARG_W)==true</tt>.  */
    RCP<LinearOpWithSolveBase<Scalar> > get_W() const;

    /** \brief Determines if an extended output argument of type <tt>ObjectType</tt> is supported. */
    template<typename ObjectType>
    bool supports() const;
    /** \brief Set an extended output argument of type <tt>ObjectType</tt> in OutArgs. Precondition: <tt>supports()==true</tt>. */
    template<typename ObjectType>
    void set(const RCP<const ObjectType>& uo);
    /** \brief Get an extended output argument of type <tt>ObjectType</tt> from OutArgs. Precondition: <tt>supports()==true</tt>. */
    template<typename ObjectType>
    RCP<const ObjectType> get() const;

    const DerivativeSupport& supports(EOutArgsDfDp_mp arg, int l) const;
    bool supports(EOutArgs_g_mp arg, int j) const;
    const DerivativeSupport& supports(EOutArgsDgDx_dot_mp arg, int j) const;
    const DerivativeSupport& supports(EOutArgsDgDx_mp arg, int j) const;
    const DerivativeSupport& supports(EOutArgsDgDp_mp arg, int j, int l) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_f_mp)==true</tt>.  */
    void set_f_mp( const RCP<Stokhos::ProductEpetraVector> &f_mp );
    /** \brief Precondition: <tt>supports(OUT_ARG_f_mp)==true</tt>.  */
    RCP<Stokhos::ProductEpetraVector> get_f_mp() const;
    /** \brief Precondition: <tt>supports(OUT_ARG_g_mp)==true</tt>.  */
    void set_g_mp( int j, const RCP<Stokhos::ProductEpetraVector> &g_mp_j );
    /** \brief Precondition: <tt>supports(OUT_ARG_g_mp)==true</tt>..  */
    RCP<Stokhos::ProductEpetraVector> get_g_mp(int j) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_W_mp)==true</tt>.  */
    void set_W_mp( const RCP<Stokhos::ProductEpetraOperator> &W_mp );
    /** \brief Precondition: <tt>supports(OUT_ARG_W_mp)==true</tt>.  */
    RCP<Stokhos::ProductEpetraOperator> get_W_mp() const;

    /** \brief Precondition: <tt>supports(OUT_ARG_W_op)==true</tt>.  */
    void set_W_op( const RCP<LinearOpBase<Scalar> > &W_op );
    /** \brief Precondition: <tt>supports(OUT_ARG_W_op)==true</tt>.  */
    RCP<LinearOpBase<Scalar> > get_W_op() const;
    /** \brief Precondition: <tt>supports(OUT_ARG_W_op)==true</tt>.  */
    void set_W_prec( const RCP<PreconditionerBase<Scalar> > &W_prec );
    /** \brief Precondition: <tt>supports(OUT_ARG_W_op)==true</tt>.  */
    RCP<PreconditionerBase<Scalar> > get_W_prec() const;
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

#ifdef Thyra_BUILD_HESSIAN_SUPPORT

    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_f_xx)==true</tt>.  */
    void set_hess_vec_prod_f_xx(const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_f_xx);
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_f_xp,l)==true</tt>.  */
    void set_hess_vec_prod_f_xp(int l, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_f_xp_l);
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_f_px,l)==true</tt>.  */
    void set_hess_vec_prod_f_px(int l, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_f_px_l);
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_f_pp,l1,l2)==true</tt>.  */
    void set_hess_vec_prod_f_pp(int l1, int l2, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_f_pp_l1_l2);

    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_f_xx)==true</tt>.  */
    RCP<MultiVectorBase<Scalar> > get_hess_vec_prod_f_xx() const;
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_f_xp,l)==true</tt>.  */
    RCP<MultiVectorBase<Scalar> > get_hess_vec_prod_f_xp(int l) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_f_px,l)==true</tt>.  */
    RCP<MultiVectorBase<Scalar> > get_hess_vec_prod_f_px(int l) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_f_pp,l1,l2)==true</tt>.  */
    RCP<MultiVectorBase<Scalar> > get_hess_vec_prod_f_pp(int l1, int l2) const;

    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_g_xx,j)==true</tt>.  */
    void set_hess_vec_prod_g_xx(int j, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_g_xx_j);
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_g_xp,j,l)==true</tt>.  */
    void set_hess_vec_prod_g_xp(int j, int l, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_g_xp_j_l);
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_g_px,j,l)==true</tt>.  */
    void set_hess_vec_prod_g_px(int j, int l, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_g_px_j_l);
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_g_pp,j,l1,l2)==true</tt>.  */
    void set_hess_vec_prod_g_pp(int j, int l1, int l2, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_g_pp_j_l1_l2);

    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_g_xx,j)==true</tt>.  */
    RCP<MultiVectorBase<Scalar> > get_hess_vec_prod_g_xx(int j) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_g_xp,j,l)==true</tt>.  */
    RCP<MultiVectorBase<Scalar> > get_hess_vec_prod_g_xp(int j, int l) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_g_px,j,l)==true</tt>.  */
    RCP<MultiVectorBase<Scalar> > get_hess_vec_prod_g_px(int j, int l) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_vec_prod_g_pp,j,l1,l2)==true</tt>.  */
    RCP<MultiVectorBase<Scalar> > get_hess_vec_prod_g_pp(int j, int l1, int l2) const;

    /** \brief Precondition: <tt>supports(OUT_ARG_hess_f_xx)==true</tt>.  */
    void set_hess_f_xx(const RCP<LinearOpBase<Scalar> > &hess_f_xx);
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_f_xp,l)==true</tt>.  */
    void set_hess_f_xp(int l, const RCP<LinearOpBase<Scalar> > &hess_f_xp_l);
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_f_pp,l1,l2)==true</tt>.  */
    void set_hess_f_pp(int l1, int l2, const RCP<LinearOpBase<Scalar> > &hess_f_pp_l1_l2);
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_g_xx,j)==true</tt>.  */
    void set_hess_g_xx(int j, const RCP<LinearOpBase<Scalar> > &hess_g_xx_j);
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_g_xp,j,l)==true</tt>.  */
    void set_hess_g_xp(int j, int l, const RCP<LinearOpBase<Scalar> > &hess_g_xp_j_l);
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_g_pp,j,l1,l2)==true</tt>.  */
    void set_hess_g_pp(int j, int l1, int l2, const RCP<LinearOpBase<Scalar> > &hess_g_pp_j_l1_l2);
    /** \brief Precondition: <tt>supports(OUT_ARG_H_xx)==true</tt>.  */
    void set_H_xx(const RCP<LinearOpBase<Scalar> > &H_xx);
    /** \brief Precondition: <tt>supports(OUT_ARG_H_xp,l)==true</tt>.  */
    void set_H_xp(int l, const RCP<LinearOpBase<Scalar> > &H_xp_l);
    /** \brief Precondition: <tt>supports(OUT_ARG_H_pp,l1,l2)==true</tt>.  */
    void set_H_pp(int l1, int l2, const RCP<LinearOpBase<Scalar> > &H_pp_l1_l2);

    /** \brief Precondition: <tt>supports(OUT_ARG_hess_f_xx)==true</tt>.  */
    RCP<LinearOpBase<Scalar> > get_hess_f_xx() const;
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_f_xp,l)==true</tt>.  */
    RCP<LinearOpBase<Scalar> > get_hess_f_xp(int l) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_f_pp,l1,l2)==true</tt>.  */
    RCP<LinearOpBase<Scalar> > get_hess_f_pp(int l1, int l2) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_g_xx,j)==true</tt>.  */
    RCP<LinearOpBase<Scalar> > get_hess_g_xx(int j) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_g_xp,j,l)==true</tt>.  */
    RCP<LinearOpBase<Scalar> > get_hess_g_xp(int j, int l) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_hess_g_pp,j,l1,l2)==true</tt>.  */
    RCP<LinearOpBase<Scalar> > get_hess_g_pp(int j, int l1, int l2) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_H_xx)==true</tt>.  */
    RCP<LinearOpBase<Scalar> > get_H_xx() const;
    /** \brief Precondition: <tt>supports(OUT_ARG_H_xp,l)==true</tt>.  */
    RCP<LinearOpBase<Scalar> > get_H_xp(int l) const;
    /** \brief Precondition: <tt>supports(OUT_ARG_H_pp,l1,l2)==true</tt>.  */
    RCP<LinearOpBase<Scalar> > get_H_pp(int l1, int l2) const;

#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

    void set_DfDp_mp(int l,  const MPDerivative &DfDp_mp_l);
    MPDerivative get_DfDp_mp(int l) const;
    DerivativeProperties get_DfDp_mp_properties(int l) const;
    void set_DgDx_dot_mp(int j, const MPDerivative &DgDx_dot_mp_j);
    MPDerivative get_DgDx_dot_mp(int j) const;
    DerivativeProperties get_DgDx_dot_mp_properties(int j) const;
    void set_DgDx_mp(int j, const MPDerivative &DgDx_mp_j);
    MPDerivative get_DgDx_mp(int j) const;
    DerivativeProperties get_DgDx_mp_properties(int j) const;
    void set_DgDp_mp( int j, int l, const MPDerivative &DgDp_mp_j_l );
    MPDerivative get_DgDp_mp(int j, int l) const;
    DerivativeProperties get_DgDp_mp_properties(int j, int l) const;

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
    template<typename ObjectType>
    void _setSupports( const bool supports );
    /** \brief . */
    void _setSupports( EOutArgsDfDp arg, int l, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDx_dot arg, int j, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDx arg, int j, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDp arg, int j, int l, const DerivativeSupport& );

#ifdef Thyra_BUILD_HESSIAN_SUPPORT

    /** \brief . */
    void _setSupports( EOutArgs_hess_vec_prod_f_xx arg, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_hess_vec_prod_f_xp arg, int l, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_hess_vec_prod_f_px arg, int l, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_hess_vec_prod_f_pp arg, int l1, int l2, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_hess_vec_prod_g_xx arg, int j, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_hess_vec_prod_g_xp arg, int j, int l, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_hess_vec_prod_g_px arg, int j, int l, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_hess_vec_prod_g_pp arg, int j, int l1, int l2, bool supports );

    /** \brief . */
    void _setSupports( EOutArgs_hess_f_xx arg, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_hess_f_xp arg, int l, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_hess_f_pp arg, int l1, int l2, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_hess_g_xx arg, int j, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_hess_g_xp arg, int j, int l, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_hess_g_pp arg, int j, int l1, int l2, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_H_xx arg, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_H_xp arg, int l, bool supports );
    /** \brief . */
    void _setSupports( EOutArgs_H_pp arg, int l1, int l2, bool supports );

#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

    void _setSupports( EOutArgs_g_mp arg, int j, bool supports );
    void _setSupports( EOutArgsDfDp_mp arg, int l, const DerivativeSupport& );
    void _setSupports( EOutArgsDgDx_dot_mp arg, int j, const DerivativeSupport& );
    void _setSupports( EOutArgsDgDx_mp arg, int j, const DerivativeSupport& );
    void _setSupports( EOutArgsDgDp_mp arg, int j, int l, const DerivativeSupport& );

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

    void _set_DfDp_mp_properties( int l, const DerivativeProperties &properties );
    void _set_DgDx_dot_mp_properties( int j, const DerivativeProperties &properties );
    void _set_DgDx_mp_properties( int j, const DerivativeProperties &properties );
    void _set_DgDp_mp_properties( int j, int l, const DerivativeProperties &properties );

    /** \brief . */
    void _setSupports( const OutArgs<Scalar>& inputOutArgs );
    /** \brief . */
    void _setUnsupportsAndRelated( EInArgsMembers arg );
    /** \brief . */
    void _setUnsupportsAndRelated( EOutArgsMembers arg );
#ifdef Thyra_BUILD_HESSIAN_SUPPORT
    /** \brief . */
    void _setHessianSupports( const bool supports );
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

  private:
    // types
    typedef Teuchos::Array<Evaluation<VectorBase<Scalar> > > g_t;
    typedef Teuchos::Array<Derivative<Scalar> > deriv_t;
    typedef Teuchos::Array<DerivativeProperties> deriv_properties_t;
    typedef Teuchos::Array<DerivativeSupport> supports_t;
#ifdef Thyra_BUILD_HESSIAN_SUPPORT
    typedef Teuchos::Array<RCP<LinearOpBase<Scalar> > > hess_t;
    typedef Teuchos::Array<RCP<MultiVectorBase<Scalar> > > hess_vec_t;
    typedef Teuchos::Array<bool> hess_supports_t;
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

    // data
    std::string modelEvalDescription_;
    bool supports_[NUM_E_OUT_ARGS_MEMBERS];
    supports_t supports_DfDp_; // Np
    supports_t supports_DgDx_dot_; // Ng
    supports_t supports_DgDx_; // Ng
    supports_t supports_DgDp_; // Ng x Np

#ifdef Thyra_BUILD_HESSIAN_SUPPORT

    bool supports_hess_f_xx_;
    hess_supports_t supports_hess_f_xp_;
    hess_supports_t supports_hess_f_pp_;
    hess_supports_t supports_hess_g_xx_;
    hess_supports_t supports_hess_g_xp_;
    hess_supports_t supports_hess_g_pp_;
    bool supports_H_xx_;
    hess_supports_t supports_H_xp_;
    hess_supports_t supports_H_pp_;

    bool supports_hess_vec_prod_f_xx_;
    hess_supports_t supports_hess_vec_prod_f_xp_;
    hess_supports_t supports_hess_vec_prod_f_px_;
    hess_supports_t supports_hess_vec_prod_f_pp_;
    hess_supports_t supports_hess_vec_prod_g_xx_;
    hess_supports_t supports_hess_vec_prod_g_xp_;
    hess_supports_t supports_hess_vec_prod_g_px_;
    hess_supports_t supports_hess_vec_prod_g_pp_;

#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

    Evaluation<VectorBase<Scalar> > f_;
    g_t g_; // Ng
    RCP<LinearOpWithSolveBase<Scalar> > W_;
    RCP<LinearOpBase<Scalar> > W_op_;
    RCP<PreconditionerBase<Scalar> > W_prec_;
    DerivativeProperties W_properties_;
    deriv_t DfDp_; // Np
    deriv_properties_t DfDp_properties_; // Np
    deriv_t DgDx_dot_; // Ng
    deriv_t DgDx_; // Ng
    deriv_properties_t DgDx_dot_properties_; // Ng
    deriv_properties_t DgDx_properties_; // Ng
    deriv_t DgDp_; // Ng x Np
    deriv_properties_t DgDp_properties_; // Ng x Np

#ifdef Thyra_BUILD_HESSIAN_SUPPORT

    RCP<LinearOpBase<Scalar> > hess_f_xx_;
    hess_t hess_f_xp_;
    hess_t hess_f_pp_;
    hess_t hess_g_xx_;
    hess_t hess_g_xp_;
    hess_t hess_g_pp_;
    RCP<LinearOpBase<Scalar> > H_xx_;
    hess_t H_xp_;
    hess_t H_pp_;

    RCP<MultiVectorBase<Scalar> > hess_vec_prod_f_xx_;
    hess_vec_t hess_vec_prod_f_xp_; // Np
    hess_vec_t hess_vec_prod_f_px_; // Np
    hess_vec_t hess_vec_prod_f_pp_; // Np x Np

    hess_vec_t hess_vec_prod_g_xx_; // Ng
    hess_vec_t hess_vec_prod_g_xp_; // Ng x Np
    hess_vec_t hess_vec_prod_g_px_; // Ng x Np
    hess_vec_t hess_vec_prod_g_pp_; // Ng x Np x Np

#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

    Teuchos::Array<bool> supports_g_mp_; //Ng
    supports_t supports_DfDp_mp_; // Np_mp
    supports_t supports_DgDx_dot_mp_; // Ng_mp
    supports_t supports_DgDx_mp_; // Ng_mp
    supports_t supports_DgDp_mp_; // Ng_mp x Np_mp
    Teuchos::Array< RCP< Stokhos::ProductEpetraVector > > g_mp_;
    RCP<Stokhos::ProductEpetraVector> f_mp_;
    RCP<Stokhos::ProductEpetraOperator> W_mp_;
    Teuchos::Array<MPDerivative> DfDp_mp_;
    Teuchos::Array<MPDerivative> DgDx_dot_mp_;
    Teuchos::Array<MPDerivative> DgDx_mp_;
    Teuchos::Array<MPDerivative> DgDp_mp_;
    deriv_properties_t DfDp_mp_properties_;
    deriv_properties_t DgDx_dot_mp_properties_;
    deriv_properties_t DgDx_mp_properties_;
    deriv_properties_t DgDp_mp_properties_;

#ifdef HAVE_THYRA_ME_POLYNOMIAL
   RCP<Teuchos::Polynomial< VectorBase<Scalar> > > f_poly_;
#endif // HAVE_THYRA_ME_POLYNOMIAL
    mutable bool isFailed_;

    std::map<std::string,Teuchos::any> extended_outargs_;

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

#ifdef Thyra_BUILD_HESSIAN_SUPPORT

    void assert_supports(
      EOutArgs_hess_vec_prod_f_xx arg
      ) const;
    void assert_supports(
      EOutArgs_hess_vec_prod_f_xp arg, int l
      ) const;
    void assert_supports(
      EOutArgs_hess_vec_prod_f_px arg, int l
      ) const;
    void assert_supports(
      EOutArgs_hess_vec_prod_f_pp arg, int l1, int l2
      ) const;
    void assert_supports(
      EOutArgs_hess_vec_prod_g_xx arg, int j
      ) const;
    void assert_supports(
      EOutArgs_hess_vec_prod_g_xp arg, int j, int l
      ) const;
    void assert_supports(
      EOutArgs_hess_vec_prod_g_px arg, int j, int l
      ) const;
    void assert_supports(
      EOutArgs_hess_vec_prod_g_pp arg, int j, int l1, int l2
      ) const;

    void assert_supports(
      EOutArgs_hess_f_xx arg
      ) const;
    void assert_supports(
      EOutArgs_hess_f_xp arg, int l
      ) const;
    void assert_supports(
      EOutArgs_hess_f_pp arg, int l1, int l2
      ) const;
    void assert_supports(
      EOutArgs_hess_g_xx arg, int j
      ) const;
    void assert_supports(
      EOutArgs_hess_g_xp arg, int j, int l
      ) const;
    void assert_supports(
      EOutArgs_hess_g_pp arg, int j, int l1, int l2
      ) const;
    void assert_supports(
      EOutArgs_H_xx arg
      ) const;
    void assert_supports(
      EOutArgs_H_xp arg, int l
      ) const;
    void assert_supports(
      EOutArgs_H_pp arg, int l1, int l2
      ) const;

#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

    void assert_supports(EOutArgs_g_mp arg, int j) const;
    void assert_supports(
      EOutArgsDfDp_mp arg, int l,
      const MPDerivative &deriv = MPDerivative()
      ) const;
    void assert_supports(
      EOutArgsDgDx_dot_mp arg, int j,
      const MPDerivative &deriv = MPDerivative()
      ) const;
    void assert_supports(
      EOutArgsDgDx_mp arg, int j,
      const MPDerivative &deriv = MPDerivative()
      ) const;
    void assert_supports(
      EOutArgsDgDp_mp arg, int j, int l,
      const MPDerivative &deriv = MPDerivative()
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
    void set_Np_Ng(int Np, int Ng);
    /** \brief . */
    void setSupports( EInArgsMembers arg, bool supports = true );
    /** \brief . */
    void setSupports( const InArgs<Scalar>& inputInArgs, const int Np = -1 );
    /** \brief Set support for specific extended data types. */
    template<typename ObjectType>
    void setSupports(const bool supports = true);
    /** \brief . */
    void setUnsupportsAndRelated( EInArgsMembers arg );

    void setSupports( EInArgs_p_mp arg, int l, bool supports);
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

#ifdef Thyra_BUILD_HESSIAN_SUPPORT
    /** \brief . */
    void setSupports(EOutArgs_hess_vec_prod_f_xx arg, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_hess_vec_prod_f_xp arg, int l, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_hess_vec_prod_f_px arg, int l, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_hess_vec_prod_f_pp arg, int l1, int l2, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_hess_vec_prod_g_xx arg, int j, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_hess_vec_prod_g_xp arg, int j, int l, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_hess_vec_prod_g_px arg, int j, int l, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_hess_vec_prod_g_pp arg, int j, int l1, int l2, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_hess_f_xx arg, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_hess_f_xp arg, int l, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_hess_f_pp arg, int l1, int l2, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_hess_g_xx arg, int j, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_hess_g_xp arg, int j, int l, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_hess_g_pp arg, int j, int l1, int l2, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_H_xx arg, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_H_xp arg, int l, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgs_H_pp arg, int l1, int l2, bool supports = true );
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

    /** \brief Set support for specific extended data types. */
    template<typename ObjectType>
    void setSupports(const bool supports = true);

    void setSupports( EOutArgs_g_mp arg, int j, bool supports);
    void setSupports(EOutArgsDfDp_mp arg, int l, const DerivativeSupport& );
    void setSupports(EOutArgsDgDx_dot_mp arg, int j, const DerivativeSupport& );
    void setSupports(EOutArgsDgDx_mp arg, int j, const DerivativeSupport& );
    void setSupports(EOutArgsDgDp_mp arg, int j, int l, const DerivativeSupport& );

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

    void set_DfDp_mp_properties( int l, const DerivativeProperties &properties );
    void set_DgDx_dot_mp_properties( int j, const DerivativeProperties &properties );
    void set_DgDx_mp_properties( int j, const DerivativeProperties &properties );
    void set_DgDp_mp_properties( int j, int l, const DerivativeProperties &properties );

    /** \brief . */
    void setSupports( const OutArgs<Scalar>& inputOutArgs );
   /** \brief . */
    void setUnsupportsAndRelated( EInArgsMembers arg );
    /** \brief Must be called after the above function. */
    void setUnsupportsAndRelated( EOutArgsMembers arg );

#ifdef Thyra_BUILD_HESSIAN_SUPPORT
    /** \brief . */
    void setHessianSupports( const bool supports = true );
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT
   };

  //@}

  /** \brief constructor */
  //@{
 
  /** \brief . */
  ModelEvaluatorBase();

  //@}

private:
  // Not defined and not to be called
  ModelEvaluatorBase(const ModelEvaluatorBase&);
  ModelEvaluatorBase& operator=(const ModelEvaluatorBase&);

}; // ModelEvaluatorBase


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

// Extended InArgs
template<class Scalar>
template<typename ObjectType>
bool Thyra::ModelEvaluatorBase::InArgs<Scalar>::supports() const
{
  std::map<std::string,Teuchos::any>::const_iterator search =
    extended_inargs_.find(typeid(ObjectType).name());

  if (search == extended_inargs_.end())
    return false;

  return true;
} 

template<class Scalar>
template<typename ObjectType>
void Thyra::ModelEvaluatorBase::InArgs<Scalar>::set(const Teuchos::RCP<const ObjectType>& eo)
{
  std::map<std::string,Teuchos::any>::iterator search = extended_inargs_.find(typeid(ObjectType).name());
  TEUCHOS_TEST_FOR_EXCEPTION(search == extended_inargs_.end(),
                             std::runtime_error,
                             "ERROR: InArgs::set<ObjectType>() was called with unsupported extended data type \""
                             << typeid(ObjectType).name() << "\"!");

  search->second = Teuchos::any(eo);
}

template<class Scalar>
template<typename ObjectType>
Teuchos::RCP<const ObjectType>
Thyra::ModelEvaluatorBase::InArgs<Scalar>::get() const
{
  std::map<std::string,Teuchos::any>::const_iterator search = extended_inargs_.find(typeid(ObjectType).name());
  TEUCHOS_TEST_FOR_EXCEPTION(search == extended_inargs_.end(),
                             std::runtime_error,
                             "ERROR: InArgs::get<ObjectType>() was called with unsupported extended data type \""
                             << typeid(ObjectType).name() << "\"!");

  return Teuchos::any_cast<Teuchos::RCP<const ObjectType> >(search->second);
}

template<class Scalar>
template<class ObjectType>
void Thyra::ModelEvaluatorBase::InArgsSetup<Scalar>::
setSupports(const bool in_supports)
{
  this->template _setSupports<ObjectType>(in_supports);
}

template<class Scalar>
template<class ObjectType>
void Thyra::ModelEvaluatorBase::InArgs<Scalar>::
_setSupports(const bool in_supports)
{
  if (in_supports)
    // When supports() is called, the map is searched to check for
    // supporting a type. If we support the type, we will insert an
    // empty placholder for now so that the search is successful for
    // support checks.
    this->extended_inargs_[typeid(ObjectType).name()] = Teuchos::any();
  else {
    // if false, remove the entry
    std::map<std::string,Teuchos::any>::iterator search =
      this->extended_inargs_.find(typeid(ObjectType).name());

    if (search != this->extended_inargs_.end())
      this->extended_inargs_.erase(typeid(ObjectType).name());
  }
}

// Extended OutArgs
template<class Scalar>
template<typename ObjectType>
bool Thyra::ModelEvaluatorBase::OutArgs<Scalar>::supports() const
{
  std::map<std::string,Teuchos::any>::const_iterator search =
    extended_outargs_.find(typeid(ObjectType).name());

  if (search == extended_outargs_.end())
    return false;

  return true;
} 

template<class Scalar>
template<typename ObjectType>
void Thyra::ModelEvaluatorBase::OutArgs<Scalar>::set(const Teuchos::RCP<const ObjectType>& eo)
{
  std::map<std::string,Teuchos::any>::iterator search = extended_outargs_.find(typeid(ObjectType).name());
  TEUCHOS_TEST_FOR_EXCEPTION(search == extended_outargs_.end(),
                             std::runtime_error,
                             "ERROR: OutArgs::set<ObjectType>() was called with unsupported extended data type \""
                             << typeid(ObjectType).name() << "\"!");

  search->second = Teuchos::any(eo);
}

template<class Scalar>
template<typename ObjectType>
Teuchos::RCP<const ObjectType>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>::get() const
{
  std::map<std::string,Teuchos::any>::const_iterator search = extended_outargs_.find(typeid(ObjectType).name());
  TEUCHOS_TEST_FOR_EXCEPTION(search == extended_outargs_.end(),
                             std::runtime_error,
                             "ERROR: OutArgs::get<ObjectType>() was called with unsupported extended data type \""
                             << typeid(ObjectType).name() << "\"!");

  return Teuchos::any_cast<Teuchos::RCP<const ObjectType> >(search->second);
}

template<class Scalar>
template<class ObjectType>
void Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar>::
setSupports(const bool in_supports)
{
  this->template _setSupports<ObjectType>(in_supports);
}

template<class Scalar>
template<class ObjectType>
void Thyra::ModelEvaluatorBase::OutArgs<Scalar>::
_setSupports(const bool in_supports)
{
  if (in_supports)
    // When supports() is called, the map is searched to check for
    // supporting a type. If we support the type, we will insert an
    // empty placholder for now so that the search is successful for
    // support checks.
    this->extended_outargs_[typeid(ObjectType).name()] = Teuchos::any();
  else {
    // if false, remove the entry
    std::map<std::string,Teuchos::any>::iterator search =
      this->extended_outargs_.find(typeid(ObjectType).name());

    if (search != this->extended_outargs_.end())
      this->extended_outargs_.erase(typeid(ObjectType).name());
  }
}

//
// Thyra_MEB_helper_functions_grp
//


inline
std::string Thyra::toString(ModelEvaluatorBase::EInArgsMembers arg)
{
  switch(arg) {
    case ModelEvaluatorBase::IN_ARG_x_dot_dot:
      return "IN_ARG_x_dot_dot";
    case ModelEvaluatorBase::IN_ARG_x_dot:
      return "IN_ARG_x_dot";
    case ModelEvaluatorBase::IN_ARG_x:
      return "IN_ARG_x";
    case ModelEvaluatorBase::IN_ARG_x_dot_poly:
      return "IN_ARG_x_dot_poly";
    case ModelEvaluatorBase::IN_ARG_x_poly:
      return "IN_ARG_x_poly";
    case ModelEvaluatorBase::IN_ARG_x_dot_mp:
      return "IN_ARG_x_dot_mp";
    case ModelEvaluatorBase::IN_ARG_x_mp:
      return "IN_ARG_x_mp";
    case ModelEvaluatorBase::IN_ARG_t:
      return "IN_ARG_t";
    case ModelEvaluatorBase::IN_ARG_alpha:
      return "IN_ARG_alpha";
    case ModelEvaluatorBase::IN_ARG_beta:
      return "IN_ARG_beta";
    case ModelEvaluatorBase::IN_ARG_W_x_dot_dot_coeff:
      return "IN_ARG_W_x_dot_dot_coeff";
    case ModelEvaluatorBase::IN_ARG_step_size:
      return "IN_ARG_step_size";
    case ModelEvaluatorBase::IN_ARG_stage_number:
      return "IN_ARG_stage_number";
#ifdef TEUCHOS_DEBUG
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
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
    case ModelEvaluatorBase::OUT_ARG_f_mp:
      return "OUT_ARG_f_mp";
    case ModelEvaluatorBase::OUT_ARG_W_mp:
      return "OUT_ARG_W_mp";
    case ModelEvaluatorBase::OUT_ARG_W_op:
      return "OUT_ARG_W_op";
    case ModelEvaluatorBase::OUT_ARG_W_prec:
      return "OUT_ARG_W_prec";
    case ModelEvaluatorBase::OUT_ARG_f_poly:
      return "OUT_ARG_f_poly";
#ifdef TEUCHOS_DEBUG
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
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
      TEUCHOS_TEST_FOR_EXCEPT(true);
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
      TEUCHOS_TEST_FOR_EXCEPT(true);
#endif
  }
  return ModelEvaluatorBase::DERIV_MV_BY_COL; // Should never execute this!
}


#endif // THYRA_MODEL_EVALUATOR_BASE_DECL_HPP
