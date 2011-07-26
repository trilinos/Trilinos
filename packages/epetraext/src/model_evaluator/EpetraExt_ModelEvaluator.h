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

#include "EpetraExt_ConfigDefs.h"
#include "EpetraExt_PolynomialVectorTraits.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_Polynomial.hpp"
#include "Teuchos_Array.hpp"

#ifdef HAVE_PYTRILINOS
#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif
#endif

class Epetra_Map;
class Epetra_Vector;
class Epetra_Operator;

// Forward declaration of Stochastic Galerkin (SG) argument types
namespace Stokhos {
  class EpetraVectorOrthogPoly;
  class EpetraMultiVectorOrthogPoly;
  class EpetraOperatorOrthogPoly;
  template <typename ordinal_type, typename scalar_type> class OrthogPolyBasis;
  template <typename ordinal_type, typename scalar_type> class Quadrature;
  template <typename ordinal_type, typename scalar_type> class StandardStorage;
  template <typename ordinal_type, typename scalar_type, typename node_type> class OrthogPolyExpansion;

  class ProductEpetraVector;
  class ProductEpetraMultiVector;
  class ProductEpetraOperator;
}

namespace EpetraExt {

/** \brief Base interface for evaluating a stateless "model".
 *
 * ToDo: Finish Documentation!
 */
class ModelEvaluator : virtual public Teuchos::Describable {
public:

  /** \name Public types */
  //@{

  typedef Teuchos::RCP<const Stokhos::ProductEpetraVector> mp_const_vector_t;
  typedef Teuchos::RCP<const Stokhos::ProductEpetraMultiVector> mp_const_multivector_t;
  typedef Teuchos::RCP<const Stokhos::ProductEpetraOperator > mp_const_operator_t;
  typedef Teuchos::RCP<Stokhos::ProductEpetraVector> mp_vector_t;
  typedef Teuchos::RCP<Stokhos::ProductEpetraMultiVector> mp_multivector_t;
  typedef Teuchos::RCP<Stokhos::ProductEpetraOperator > mp_operator_t;

  /** \brief.  */
  enum EInArgsMembers {
    IN_ARG_x_dot
    ,IN_ARG_x
    ,IN_ARG_x_dot_poly ///< Time derivative vector Taylor polynomial
    ,IN_ARG_x_poly    ///< Solution vector Taylor polynomial
    ,IN_ARG_x_dot_sg  ///< Stochastic Galerkin time derivative vector polynomial
    ,IN_ARG_x_sg      ///< Stochastic Galerkin solution vector polynomial
    ,IN_ARG_x_dot_mp  ///< Multi-point time derivative vector
    ,IN_ARG_x_mp      ///< Multi-point solution vector
    ,IN_ARG_t
    ,IN_ARG_alpha
    ,IN_ARG_beta
    ,IN_ARG_sg_basis ///< Stochastic Galerkin basis
    ,IN_ARG_sg_quadrature ///< Stochastic Galerkin quadrature
    ,IN_ARG_sg_expansion ///< Stochastic Galerkin expansion
  };
  static const int NUM_E_IN_ARGS_MEMBERS=14;

  /** \brief . */
  enum EInArgs_p_sg {
    IN_ARG_p_sg   ///< .
  };

  /** \brief . */
  enum EInArgs_p_mp {
    IN_ARG_p_mp   ///< .
  };

  /** \brief . */
  class InArgs {
  public:

    //! Short-hand for stochastic Galerkin vector type
    typedef Teuchos::RefCountPtr<const Stokhos::EpetraVectorOrthogPoly> sg_const_vector_t;
        
    /** \brief. */
    InArgs();
    /** \brief . */
    std::string modelEvalDescription() const;
    /** \brief .  */
    int Np() const;
    /** \brief. */
    void set_x_dot( const Teuchos::RefCountPtr<const Epetra_Vector> &x_dot );
    /** \brief. */
    Teuchos::RefCountPtr<const Epetra_Vector> get_x_dot() const;
    /** \brief. */
    void set_x( const Teuchos::RefCountPtr<const Epetra_Vector> &x );
    /** \brief Set solution vector Taylor polynomial. */
    Teuchos::RefCountPtr<const Epetra_Vector> get_x() const;
    void set_x_poly(
      const Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > &x_poly
      );
    /** \brief Get solution vector Taylor polynomial.  */
    Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > get_x_poly() const;
    /** \brief Set time derivative vector Taylor polynomial.  */
    void set_x_dot_poly(
      const Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > &x_dot_poly
      );
    /** \brief Get time derivative vector Taylor polynomial.  */
    Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > get_x_dot_poly() const;
    /** \brief Set stochastic Galerkin solution vector polynomial.  */
    void set_x_sg(const sg_const_vector_t &x_sg);
    /** \brief Get stochastic Galerkin solution vector polynomial.  */
    sg_const_vector_t get_x_sg() const;
    /** \brief Set stochastic Galerkin time derivative vector polynomial.  */
    void set_x_dot_sg(const sg_const_vector_t &x_dot_sg);
    /** \brief Get stochastic Galerkin time derivative vector polynomial.  */
    sg_const_vector_t get_x_dot_sg() const;
    /** \brief Set multi-point solution vector.  */
    void set_x_mp(const mp_const_vector_t &x_mp);
    /** \brief Get multi-point solution vector.  */
    mp_const_vector_t get_x_mp() const;
    /** \brief Set multi-point time derivative vector.  */
    void set_x_dot_mp(const mp_const_vector_t &x_dot_mp);
    /** \brief Get multi-point time derivative vector.  */
    mp_const_vector_t get_x_dot_mp() const;
    /** \brief. */
    void set_p( int l, const Teuchos::RefCountPtr<const Epetra_Vector> &p_l );
    /** \brief. */
    Teuchos::RefCountPtr<const Epetra_Vector> get_p(int l) const;
    /** \brief Set stochastic Galerkin vector polynomial parameter. */
    void set_p_sg( int l, const sg_const_vector_t &p_sg_l );
    /** \brief Get stochastic Galerkin vector polynomial parameter. */
    sg_const_vector_t get_p_sg(int l) const;
    /** \brief Set multi-point parameter vector. */
    void set_p_mp( int l, const mp_const_vector_t &p_mp_l );
    /** \brief Get multi-point parameter vector. */
    mp_const_vector_t get_p_mp(int l) const;
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
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > get_sg_basis() const;
    /** \brief. */
    void set_sg_basis( const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis );
    /** \brief. */
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > get_sg_quadrature() const;
    /** \brief. */
    void set_sg_quadrature( const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& quad );
    /** \brief. */
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double,Stokhos::StandardStorage<int,double> > > get_sg_expansion() const;
    /** \brief. */
    void set_sg_expansion( const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double,Stokhos::StandardStorage<int,double> > >& exp );
    /** \brief. */
    bool supports(EInArgsMembers arg) const;
    /** Whether p_sg is supported for parameter vector l */
    bool supports(EInArgs_p_sg arg, int l) const;
    /** Whether p_mp is supported for parameter vector l */
    bool supports(EInArgs_p_mp arg, int l) const;
  protected:
    /** \brief . */
    void _setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void _set_Np(int Np);
    /** \brief . */
    void _setSupports( EInArgsMembers arg, bool supports );
    /** \brief . */
    void _setSupports( EInArgs_p_sg arg, int l, bool supports );
    /** \brief . */
    void _setSupports( EInArgs_p_mp arg, int l, bool supports );
  private:
    // types
    typedef Teuchos::Array<Teuchos::RefCountPtr<const Epetra_Vector> > p_t;
    typedef Teuchos::Array<sg_const_vector_t > p_sg_t;
    typedef Teuchos::Array<mp_const_vector_t > p_mp_t;
    typedef Teuchos::Array<bool> supports_p_sg_t;
    // data
    std::string modelEvalDescription_;
    Teuchos::RefCountPtr<const Epetra_Vector>  x_dot_;
    Teuchos::RefCountPtr<const Epetra_Vector>  x_;
    Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > x_dot_poly_;
    Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > x_poly_;
    sg_const_vector_t                          x_dot_sg_;
    sg_const_vector_t                          x_sg_;
    mp_const_vector_t                          x_dot_mp_;
    mp_const_vector_t                          x_mp_;
    p_t                                        p_;
    p_sg_t                                     p_sg_;
    p_mp_t                                     p_mp_;
    double                                     t_;
    double                                     alpha_;
    double                                     beta_;
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > sg_basis_;
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > sg_quad_;
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double,Stokhos::StandardStorage<int,double> > > sg_exp_;
    bool supports_[NUM_E_IN_ARGS_MEMBERS];
    supports_p_sg_t supports_p_sg_; // Np
    supports_p_sg_t supports_p_mp_; // Np
    // functions
    void assert_supports(EInArgsMembers arg) const;
    void assert_supports(EInArgs_p_sg arg, int l) const;
    void assert_supports(EInArgs_p_mp arg, int l) const;
    void assert_l(int l) const;
  };

  /** \brief. */
  enum EEvalType {
    EVAL_TYPE_EXACT                ///< Exact function evaluation
    ,EVAL_TYPE_APPROX_DERIV        ///< An approximate derivative (i.e. for a Jacobian)
    ,EVAL_TYPE_VERY_APPROX_DERIV   ///< A very approximate derivative (i.e. for a preconditioner)
  };

  /** \brief . */
  template<class ObjType>
  class Evaluation : public Teuchos::RefCountPtr<ObjType> {
  public:
    /** \brief . */
    Evaluation() : evalType_(EVAL_TYPE_EXACT) {}
    /** \brief . */
    Evaluation( const Teuchos::RefCountPtr<ObjType> &obj )
      : Teuchos::RefCountPtr<ObjType>(obj), evalType_(EVAL_TYPE_EXACT) {}
    /** \brief . */
    Evaluation( const Teuchos::RefCountPtr<ObjType> &obj, EEvalType evalType )
      : Teuchos::RefCountPtr<ObjType>(obj), evalType_(evalType) {}
    /** \brief . */
    EEvalType getType() const { return evalType_; }
    /** \brief . */
    void reset( const Teuchos::RefCountPtr<ObjType> &obj, EEvalType evalType ) 
    { this->operator=(obj); evalType_ = evalType; }
  private:
    EEvalType                      evalType_;
  };
  
  /** \brief . */
  enum EDerivativeMultiVectorOrientation {
    DERIV_MV_BY_COL           ///< .
    ,DERIV_TRANS_MV_BY_ROW    ///< .
  };

  /** \brief . */
  enum EDerivativeLinearOp { DERIV_LINEAR_OP };

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
    DerivativeSupport(
      EDerivativeLinearOp, EDerivativeMultiVectorOrientation mvOrientation )
      :supportsLinearOp_(true), supportsMVByCol_(mvOrientation==DERIV_MV_BY_COL)
      ,supportsTransMVByRow_(mvOrientation==DERIV_TRANS_MV_BY_ROW)
      {}
    /** \brief . */
    DerivativeSupport(
      EDerivativeMultiVectorOrientation mvOrientation1,
      EDerivativeMultiVectorOrientation mvOrientation2
      )
      :supportsLinearOp_(false)
      ,supportsMVByCol_(
        mvOrientation1==DERIV_MV_BY_COL||mvOrientation2==DERIV_MV_BY_COL )
      ,supportsTransMVByRow_(
        mvOrientation1==DERIV_TRANS_MV_BY_ROW||mvOrientation2==DERIV_TRANS_MV_BY_ROW )
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
    EDerivativeLinearity linearity;
    /** \brief . */
    ERankStatus rank;
    /** \brief . */
    bool supportsAdjoint;
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
  class DerivativeMultiVector {
  public:
    /** \brief . */
    DerivativeMultiVector() {}
    /** \brief . */
    DerivativeMultiVector(
      const Teuchos::RefCountPtr<Epetra_MultiVector> &mv
      ,const EDerivativeMultiVectorOrientation orientation = DERIV_MV_BY_COL
      ,const Teuchos::Array<int> &paramIndexes = Teuchos::Array<int>()
      ) : mv_(mv), orientation_(orientation), paramIndexes_(paramIndexes) {}
    /** \brief . */
    void changeOrientation( const EDerivativeMultiVectorOrientation orientation )
      { orientation_ = orientation; };
    /** \brief . */
    Teuchos::RefCountPtr<Epetra_MultiVector> getMultiVector() const
      { return mv_; }
    /** \brief . */
    EDerivativeMultiVectorOrientation getOrientation() const
      { return orientation_; }
    /** \brief . */
    const Teuchos::Array<int>& getParamIndexes() const
      { return paramIndexes_; }
  private:
    Teuchos::RefCountPtr<Epetra_MultiVector> mv_;
    EDerivativeMultiVectorOrientation orientation_;
    Teuchos::Array<int> paramIndexes_;
  };

  /** \brief Simple aggregate class that stores a derivative object
   * as a general linear operator or as a multi-vector.
   */
  class Derivative {
  public:
    /** \brief . */
    Derivative() {}
    /** \brief . */
    Derivative( const Teuchos::RefCountPtr<Epetra_Operator> &lo )
      : lo_(lo) {}
     /** \brief . */
    Derivative(
      const Teuchos::RefCountPtr<Epetra_MultiVector> &mv
      ,const EDerivativeMultiVectorOrientation orientation = DERIV_MV_BY_COL
      ) : dmv_(mv,orientation) {}
   /** \brief . */
    Derivative( const DerivativeMultiVector &dmv )
      : dmv_(dmv) {}
    /** \brief . */
    Teuchos::RefCountPtr<Epetra_Operator> getLinearOp() const
      { return lo_; }
    /** \brief . */
    Teuchos::RefCountPtr<Epetra_MultiVector> getMultiVector() const
      { return dmv_.getMultiVector(); }
    /** \brief . */
    EDerivativeMultiVectorOrientation getMultiVectorOrientation() const
      { return dmv_.getOrientation(); }
    /** \brief . */
    DerivativeMultiVector getDerivativeMultiVector() const
      { return dmv_; }
    /** \brief . */
    bool isEmpty() const
        { return !lo_.get() && !dmv_.getMultiVector().get(); }
  private:
    Teuchos::RefCountPtr<Epetra_Operator> lo_;
    DerivativeMultiVector dmv_;
  };

  /** \brief Simple aggregate struct that stores a preconditioner
   * as an Epetra_Operator and a bool, about whether it is inverted or not.
   */
  struct Preconditioner {
    /** \brief Default constructor of null Operatir */
    Preconditioner() : PrecOp(Teuchos::null), isAlreadyInverted(false) {}
    /** \brief Usable constructor to set the (Epetra_Operator,bool) pair */
    Preconditioner(const Teuchos::RCP<Epetra_Operator>& PrecOp_,
                     bool isAlreadyInverted_ )
      : PrecOp(PrecOp_), isAlreadyInverted(isAlreadyInverted_) {}
    /** \brief Accessor for the Epetra_Operator */
    Teuchos::RCP<Epetra_Operator> PrecOp;
    /** \brief Bool flag. When isAlreadyInverted==true, then the PrecOp
        is an approximation to the inverse of the Matrix. If false, then
        this Operator is an approximation to the matrix, and still must
        be sent to a preconditioning algorithm like ILU. */
    bool isAlreadyInverted;
  };

  /** \brief Simple aggregate class for a derivative object represented as a
   * column-wise multi-vector or its transpose as a row-wise multi-vector.
   */
  class SGDerivativeMultiVector {
  public:
    /** \brief . */
    SGDerivativeMultiVector() {}
    /** \brief . */
    SGDerivativeMultiVector(
      const Teuchos::RefCountPtr< Stokhos::EpetraMultiVectorOrthogPoly > &mv
      ,const EDerivativeMultiVectorOrientation orientation = DERIV_MV_BY_COL
      ,const Teuchos::Array<int> &paramIndexes = Teuchos::Array<int>()
      ) : mv_(mv), orientation_(orientation), paramIndexes_(paramIndexes) {}
    /** \brief . */
    void changeOrientation( const EDerivativeMultiVectorOrientation orientation )
      { orientation_ = orientation; };
    /** \brief . */
    Teuchos::RefCountPtr< Stokhos::EpetraMultiVectorOrthogPoly > getMultiVector() const
      { return mv_; }
    /** \brief . */
    EDerivativeMultiVectorOrientation getOrientation() const
      { return orientation_; }
    /** \brief . */
    const Teuchos::Array<int>& getParamIndexes() const
      { return paramIndexes_; }
  private:
    Teuchos::RefCountPtr< Stokhos::EpetraMultiVectorOrthogPoly > mv_;
    EDerivativeMultiVectorOrientation orientation_;
    Teuchos::Array<int> paramIndexes_;
  };

  /** \brief Simple aggregate class that stores a derivative object
   * as a general linear operator or as a multi-vector.
   */
  class SGDerivative {
  public:
    /** \brief . */
    SGDerivative() {}
    /** \brief . */
    SGDerivative( const Teuchos::RefCountPtr< Stokhos::EpetraOperatorOrthogPoly > &lo )
      : lo_(lo) {}
     /** \brief . */
    SGDerivative(
      const Teuchos::RefCountPtr< Stokhos::EpetraMultiVectorOrthogPoly > &mv
      ,const EDerivativeMultiVectorOrientation orientation = DERIV_MV_BY_COL
      ) : dmv_(mv,orientation) {}
   /** \brief . */
    SGDerivative( const SGDerivativeMultiVector &dmv )
      : dmv_(dmv) {}
    /** \brief . */
    Teuchos::RefCountPtr< Stokhos::EpetraOperatorOrthogPoly > getLinearOp() const
      { return lo_; }
    /** \brief . */
    Teuchos::RefCountPtr< Stokhos::EpetraMultiVectorOrthogPoly > getMultiVector() const
      { return dmv_.getMultiVector(); }
    /** \brief . */
    EDerivativeMultiVectorOrientation getMultiVectorOrientation() const
      { return dmv_.getOrientation(); }
    /** \brief . */
    SGDerivativeMultiVector getDerivativeMultiVector() const
      { return dmv_; }
    /** \brief . */
    bool isEmpty() const
        { return !lo_.get() && !dmv_.getMultiVector().get(); }
  private:
    Teuchos::RefCountPtr< Stokhos::EpetraOperatorOrthogPoly > lo_;
    SGDerivativeMultiVector dmv_;
  };

  /** \brief Simple aggregate class for a derivative object represented as a
   * column-wise multi-vector or its transpose as a row-wise multi-vector.
   */
  class MPDerivativeMultiVector {
  public:
    
    /** \brief . */
    MPDerivativeMultiVector() {}
    /** \brief . */
    MPDerivativeMultiVector(
      const mp_multivector_t &mv
      ,const EDerivativeMultiVectorOrientation orientation = DERIV_MV_BY_COL
      ,const Teuchos::Array<int> &paramIndexes = Teuchos::Array<int>()
      ) : mv_(mv), orientation_(orientation), paramIndexes_(paramIndexes) {}
    /** \brief . */
    void changeOrientation( const EDerivativeMultiVectorOrientation orientation )
      { orientation_ = orientation; };
    /** \brief . */
    mp_multivector_t getMultiVector() const
      { return mv_; }
    /** \brief . */
    EDerivativeMultiVectorOrientation getOrientation() const
      { return orientation_; }
    /** \brief . */
    const Teuchos::Array<int>& getParamIndexes() const
      { return paramIndexes_; }
  private:
    mp_multivector_t mv_;
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
    MPDerivative( const mp_operator_t &lo )
      : lo_(lo) {}
     /** \brief . */
    MPDerivative(
      const mp_multivector_t &mv
      ,const EDerivativeMultiVectorOrientation orientation = DERIV_MV_BY_COL
      ) : dmv_(mv,orientation) {}
    /** \brief . */
    MPDerivative( const MPDerivativeMultiVector &dmv )
      : dmv_(dmv) {}
    /** \brief . */
    mp_operator_t getLinearOp() const
      { return lo_; }
    /** \brief . */
    mp_multivector_t getMultiVector() const
      { return dmv_.getMultiVector(); }
    /** \brief . */
    EDerivativeMultiVectorOrientation getMultiVectorOrientation() const
      { return dmv_.getOrientation(); }
    /** \brief . */
    MPDerivativeMultiVector getDerivativeMultiVector() const
      { return dmv_; }
    /** \brief . */
    bool isEmpty() const
        { return !lo_.get() && !dmv_.getMultiVector().get(); }
  private:
    mp_operator_t lo_;
    MPDerivativeMultiVector dmv_;
  };

  /** \brief.  */
  enum EOutArgsMembers {
    OUT_ARG_f
    ,OUT_ARG_W
    ,OUT_ARG_f_poly   ///< Residual vector Taylor polynomial
    ,OUT_ARG_f_sg     ///< Stochastic Galerkin residual vector polynomial
    ,OUT_ARG_W_sg     ///< Stochastic Galerkin "W" operator polyomial
    ,OUT_ARG_f_mp     ///< Multi-point residual vector
    ,OUT_ARG_W_mp     ///< Multi-point "W" operator 
    ,OUT_ARG_WPrec   ///< Preconditioner operator (approx Jacobian)
  };
  static const int NUM_E_OUT_ARGS_MEMBERS=9;

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

  /** \brief . */
  enum EOutArgsDfDp_sg {
    OUT_ARG_DfDp_sg   ///< .
  };

  /** \brief . */
  enum EOutArgs_g_sg {
    OUT_ARG_g_sg   ///< .
  };

  /** \brief . */
  enum EOutArgsDgDx_dot_sg {
    OUT_ARG_DgDx_dot_sg   ///< .
  };

  /** \brief . */
  enum EOutArgsDgDx_sg {
    OUT_ARG_DgDx_sg   ///< .
  };

  /** \brief . */
  enum EOutArgsDgDp_sg {
    OUT_ARG_DgDp_sg   ///< .
  };

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

  /** \brief . */
  class OutArgs {
  public:

    //! Short-hand for stochastic Galerkin vector type
    typedef Teuchos::RefCountPtr<Stokhos::EpetraVectorOrthogPoly> sg_vector_t;

    //! Short-hand for stochastic Galerkin operator type
    typedef Teuchos::RefCountPtr<Stokhos::EpetraOperatorOrthogPoly > sg_operator_t;

    /** \brief. */
    OutArgs();
    /** \brief . */
    std::string modelEvalDescription() const;
    /** \brief .  */
    int Np() const;
    /** \brief .  */
    int Ng() const;
    /** \brief. */
    bool supports(EOutArgsMembers arg) const;
    /** \brief <tt>0 <= l && l < Np()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDfDp arg, int l) const;
    /** \brief <tt>0 <= j && j < Ng()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDx_dot arg, int j) const;
    /** \brief <tt>0 <= j && j < Ng()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDx arg, int j) const;
    /** \brief <tt>0 <= j && j < Ng()</tt> and <tt>0 <= l && l < Np()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDp arg, int j, int l) const;
    /** Whether g_sg is supported for response vector j */
    bool supports(EOutArgs_g_sg arg, int j) const;
    /** \brief <tt>0 <= l && l < Np()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDfDp_sg arg, int l) const;
    /** \brief <tt>0 <= j && j < Ng()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDx_dot_sg arg, int j) const;
    /** \brief <tt>0 <= j && j < Ng()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDx_sg arg, int j) const;
    /** \brief <tt>0 <= j && j < Ng()</tt> and <tt>0 <= l && l < Np()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDp_sg arg, int j, int l) const;
    /** \brief <tt>0 <= l && l < Np()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDfDp_mp arg, int l) const;
    /** Whether g_mp is supported for response vector j */
    bool supports(EOutArgs_g_mp arg, int j) const;
    /** \brief <tt>0 <= j && j < Ng()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDx_dot_mp arg, int j) const;
    /** \brief <tt>0 <= j && j < Ng()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDx_mp arg, int j) const;
    /** \brief <tt>0 <= j && j < Ng()</tt> and <tt>0 <= l && l < Np()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDp_mp arg, int j, int l) const;
    /** \brief. */
    void set_f( const Evaluation<Epetra_Vector> &f );
    /** \brief. */
    Evaluation<Epetra_Vector> get_f() const;
    /** \brief Set stochastic Galerkin residual vector polynomial.  */
    void set_f_sg( const sg_vector_t& f_sg );
    /** \brief Get stochastic Galerkin residual vector polynomial.  */
    sg_vector_t get_f_sg() const;
    /** \brief Set multi-point residual vector.  */
    void set_f_mp( const mp_vector_t& f_sg );
    /** \brief Get multi-point residual vector.  */
    mp_vector_t get_f_mp() const;
    /** \brief Set <tt>g(j)</tt> where <tt>0 <= j && j < this->Ng()</tt>.  */
    void set_g( int j, const Evaluation<Epetra_Vector> &g_j );
    /** \brief Get <tt>g(j)</tt> where <tt>0 <= j && j < this->Ng()</tt>.  */
    Evaluation<Epetra_Vector> get_g(int j) const;
    /** \brief Set stochastic Galerkin vector polynomial response. */
    /** <tt>0 <= j && j < this->Ng()</tt>.  */
    void set_g_sg( int j, const sg_vector_t &g_sg_j );
    /** \brief Get stochastic Galerkin vector polynomial response. */
    /** <tt>0 <= j && j < this->Ng()</tt>.  */
    sg_vector_t get_g_sg(int j) const;
    /** \brief Set multi-point response. */
    /** <tt>0 <= j && j < this->Ng()</tt>.  */
    void set_g_mp( int j, const mp_vector_t &g_mp_j );
    /** \brief Get multi-point response. */
    /** <tt>0 <= j && j < this->Ng()</tt>.  */
    mp_vector_t get_g_mp(int j) const;
    /** \brief. */
    void set_W( const Teuchos::RefCountPtr<Epetra_Operator> &W );
    void set_WPrec( const Teuchos::RefCountPtr<Epetra_Operator> &WPrec );
    /** \brief. */
    Teuchos::RefCountPtr<Epetra_Operator> get_W() const;
    Teuchos::RefCountPtr<Epetra_Operator> get_WPrec() const;
    /** \brief . */
    DerivativeProperties get_W_properties() const;
    DerivativeProperties get_WPrec_properties() const;
    /** \brief Set stochastic Galerkin W operator polynomial. */
    void set_W_sg( const sg_operator_t& W_sg );
    /** \brief Get stochastic Galerkin W operator polynomial. */
    sg_operator_t get_W_sg() const;
    /** \brief Set multi-point W. */
    void set_W_mp( const mp_operator_t& W_sg );
    /** \brief Get multi-point W. */
    mp_operator_t get_W_mp() const;
    /** \brief .  */
    void set_DfDp(int l,  const Derivative &DfDp_l);
    /** \brief .  */
    Derivative get_DfDp(int l) const;
    /** \brief . */
    DerivativeProperties get_DfDp_properties(int l) const;
    /** \brief .  */
    void set_DfDp_sg(int l,  const SGDerivative &DfDp_sg_l);
    /** \brief .  */
    SGDerivative get_DfDp_sg(int l) const;
    /** \brief . */
    DerivativeProperties get_DfDp_sg_properties(int l) const;
    /** \brief .  */
    void set_DfDp_mp(int l,  const MPDerivative &DfDp_mp_l);
    /** \brief .  */
    MPDerivative get_DfDp_mp(int l) const;
    /** \brief . */
    DerivativeProperties get_DfDp_mp_properties(int l) const;
    /** \brief .  */
    void set_DgDx_dot(int j, const Derivative &DgDx_dot_j);
    /** \brief .  */
    Derivative get_DgDx_dot(int j) const;
    /** \brief . */
    DerivativeProperties get_DgDx_dot_properties(int j) const;
    /** \brief .  */
    void set_DgDx_dot_sg(int j, const SGDerivative &DgDx_dot_j);
    /** \brief .  */
    SGDerivative get_DgDx_dot_sg(int j) const;
    /** \brief . */
    DerivativeProperties get_DgDx_dot_sg_properties(int j) const;
    /** \brief .  */
    void set_DgDx_dot_mp(int j, const MPDerivative &DgDx_dot_j);
    /** \brief .  */
    MPDerivative get_DgDx_dot_mp(int j) const;
    /** \brief . */
    DerivativeProperties get_DgDx_dot_mp_properties(int j) const;
    /** \brief .  */
    void set_DgDx(int j, const Derivative &DgDx_j);
    /** \brief .  */
    Derivative get_DgDx(int j) const;
    /** \brief . */
    DerivativeProperties get_DgDx_properties(int j) const;
    /** \brief .  */
    void set_DgDx_sg(int j, const SGDerivative &DgDx_j);
    /** \brief .  */
    SGDerivative get_DgDx_sg(int j) const;
    /** \brief . */
    DerivativeProperties get_DgDx_sg_properties(int j) const;
    /** \brief .  */
    void set_DgDx_mp(int j, const MPDerivative &DgDx_j);
    /** \brief .  */
    MPDerivative get_DgDx_mp(int j) const;
    /** \brief . */
    DerivativeProperties get_DgDx_mp_properties(int j) const;
    /** \brief .  */
    void set_DgDp( int j, int l, const Derivative &DgDp_j_l );
    /** \brief .  */
    Derivative get_DgDp(int j, int l) const;
    /** \brief . */
    DerivativeProperties get_DgDp_properties(int j, int l) const;
    /** \brief .  */
    void set_DgDp_sg( int j, int l, const SGDerivative &DgDp_sg_j_l );
    /** \brief .  */
    SGDerivative get_DgDp_sg(int j, int l) const;
    /** \brief . */
    DerivativeProperties get_DgDp_sg_properties(int j, int l) const;
    /** \brief .  */
    void set_DgDp_mp( int j, int l, const MPDerivative &DgDp_mp_j_l );
    /** \brief .  */
    MPDerivative get_DgDp_mp(int j, int l) const;
    /** \brief . */
    DerivativeProperties get_DgDp_mp_properties(int j, int l) const;

    /** \brief Set residual vector Taylor polynomial.  */
    void set_f_poly( const Teuchos::RefCountPtr<Teuchos::Polynomial<Epetra_Vector> > &f_poly );
    /** \brief Get residual vector Taylor polynomial.  */
    Teuchos::RefCountPtr<Teuchos::Polynomial<Epetra_Vector> > get_f_poly() const;
    
    
    /** \brief Return true if the function or its derivatives are set. */
    bool funcOrDerivesAreSet(EOutArgsMembers arg) const;
    
    /** \brief Set that the evaluation as a whole failed.
     *
     * Note that this function is declared as <tt>const</tt> even through it
     * technically changes the state of <tt>*this</tt> object.
     */
    void setFailed() const;
    /** \brief Return if the evaluation failed or not.
     *
     * If the evaluation failed, no assumptions should be made at all about
     * the state of the output objects.
     */
    bool isFailed() const;
    
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
    void _setSupports( EOutArgs_g_sg arg, int j, bool supports );
    /** \brief . */
    void _setSupports( EOutArgsDfDp_sg arg, int l, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDx_dot_sg arg, int j, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDx_sg arg, int j, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDp_sg arg, int j, int l, const DerivativeSupport& );

    /** \brief . */
    void _setSupports( EOutArgs_g_mp arg, int j, bool supports );
    /** \brief . */
    void _setSupports( EOutArgsDfDp_mp arg, int l, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDx_dot_mp arg, int j, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDx_mp arg, int j, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDp_mp arg, int j, int l, const DerivativeSupport& );
    /** \brief . */
    void _set_W_properties( const DerivativeProperties &W_properties );
    void _set_WPrec_properties( const DerivativeProperties &WPrec_properties );
    /** \brief . */
    void _set_DfDp_properties( int l, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDx_dot_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDx_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDp_properties( int j, int l, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DfDp_sg_properties( int l, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDx_dot_sg_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDx_sg_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDp_sg_properties( int j, int l, const DerivativeProperties &properties );

    /** \brief . */
    void _set_DfDp_mp_properties( int l, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDx_dot_mp_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDx_mp_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDp_mp_properties( int j, int l, const DerivativeProperties &properties );
  private:
    // types
    typedef Teuchos::Array<Evaluation<Epetra_Vector> > g_t;
    typedef Teuchos::Array<sg_vector_t > g_sg_t;
    typedef Teuchos::Array<mp_vector_t > g_mp_t;
    typedef Teuchos::Array<Derivative> deriv_t;
    typedef Teuchos::Array<SGDerivative> sg_deriv_t;
    typedef Teuchos::Array<MPDerivative> mp_deriv_t;
    typedef Teuchos::Array<DerivativeProperties> deriv_properties_t;
    typedef Teuchos::Array<DerivativeSupport> supports_t;
    typedef Teuchos::Array<bool> supports_g_sg_t;
    // data
    std::string modelEvalDescription_;
    mutable bool isFailed_;
    bool supports_[NUM_E_OUT_ARGS_MEMBERS];
    supports_t supports_DfDp_; // Np
    supports_t supports_DgDx_dot_; // Ng
    supports_t supports_DgDx_; // Ng
    supports_t supports_DgDp_; // Ng x Np
    supports_g_sg_t supports_g_sg_; // Ng
    supports_t supports_DfDp_sg_; // Np
    supports_t supports_DgDx_dot_sg_; // Ng
    supports_t supports_DgDx_sg_; // Ng
    supports_t supports_DgDp_sg_; // Ng x Np
    supports_g_sg_t supports_g_mp_; // Ng
    supports_t supports_DfDp_mp_; // Np_mp
    supports_t supports_DgDx_dot_mp_; // Ng_mp
    supports_t supports_DgDx_mp_; // Ng_mp
    supports_t supports_DgDp_mp_; // Ng_mp x Np_mp
    Evaluation<Epetra_Vector> f_;
    g_t g_;
    g_sg_t g_sg_;
    g_mp_t g_mp_;
    Teuchos::RefCountPtr<Epetra_Operator> W_;
    Teuchos::RefCountPtr<Epetra_Operator> WPrec_;
    DerivativeProperties W_properties_;
    DerivativeProperties WPrec_properties_;
    deriv_t DfDp_; // Np
    deriv_properties_t DfDp_properties_; // Np
    deriv_t DgDx_dot_; // Ng
    deriv_t DgDx_; // Ng
    deriv_properties_t DgDx_dot_properties_; // Ng
    deriv_properties_t DgDx_properties_; // Ng
    deriv_t DgDp_; // Ng x Np
    deriv_properties_t DgDp_properties_; // Ng x Np
    Teuchos::RefCountPtr<Teuchos::Polynomial<Epetra_Vector> > f_poly_;
    sg_vector_t f_sg_;
    sg_operator_t W_sg_;
    sg_deriv_t DfDp_sg_; // Np
    deriv_properties_t DfDp_sg_properties_; // Np
    sg_deriv_t DgDx_dot_sg_; // Ng
    sg_deriv_t DgDx_sg_; // Ng
    deriv_properties_t DgDx_dot_sg_properties_; // Ng
    deriv_properties_t DgDx_sg_properties_; // Ng
    sg_deriv_t DgDp_sg_; // Ng x Np
    deriv_properties_t DgDp_sg_properties_; // Ng x Np
    mp_vector_t f_mp_;
    mp_operator_t W_mp_;
    mp_deriv_t DfDp_mp_; // Np
    deriv_properties_t DfDp_mp_properties_; // Np
    mp_deriv_t DgDx_dot_mp_; // Ng
    mp_deriv_t DgDx_mp_; // Ng
    deriv_properties_t DgDx_dot_mp_properties_; // Ng
    deriv_properties_t DgDx_mp_properties_; // Ng
    mp_deriv_t DgDp_mp_; // Ng x Np
    deriv_properties_t DgDp_mp_properties_; // Ng x Np
    // functions
    void assert_supports(EOutArgsMembers arg) const;
    void assert_supports(EOutArgsDfDp arg, int l) const;
    void assert_supports(EOutArgsDgDx_dot arg, int j) const;
    void assert_supports(EOutArgsDgDx arg, int j) const;
    void assert_supports(EOutArgsDgDp arg, int j, int l) const;
    void assert_supports(EOutArgs_g_sg arg, int j) const;
    void assert_supports(EOutArgsDfDp_sg arg, int l) const;
    void assert_supports(EOutArgsDgDx_dot_sg arg, int j) const;
    void assert_supports(EOutArgsDgDx_sg arg, int j) const;
    void assert_supports(EOutArgsDgDp_sg arg, int j, int l) const;
    void assert_supports(EOutArgs_g_mp arg, int j) const;
    void assert_supports(EOutArgsDfDp_mp arg, int l) const;
    void assert_supports(EOutArgsDgDx_dot_mp arg, int j) const;
    void assert_supports(EOutArgsDgDx_mp arg, int j) const;
    void assert_supports(EOutArgsDgDp_mp arg, int j, int l) const;
    void assert_l(int l) const;
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

  /** \brief Get the names of the parameters associated with parameter
   * subvector l if available.
   *
   * \return Returns an RCP to a Teuchos::Array<std::string> object that
   * contains the names of the parameters.  If returnVal == Teuchos::null,
   * then there are no names available for the parameter subvector p(l).  If
   * returnVal->size() == 1, then a single name is given to the entire
   * parameter subvector.  If returnVal->size() ==
   * get_p_map(l)->GlobalNumElements(), then a name is given to every
   * parameter scalar entry.
   *
   * The default implementation return returnVal==Teuchos::null which means
   * by default, parameters have no names associated with them.
   */
  virtual Teuchos::RefCountPtr<const Teuchos::Array<std::string> > get_p_names(int l) const;

  /** \breif . */
  virtual Teuchos::RefCountPtr<const Epetra_Map> get_g_map(int j) const;

  //@}

  /** \name Initial guesses for variables/parameters */
  //@{

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_x_init() const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_x_dot_init() const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_p_init(int l) const;

  /** \brief . */
  virtual double get_t_init() const;

  //@}

  /** \name Bounds for variables/parameters */
  //@{

  /** \brief Return the value of an infinite bound.
   *
   * The default implementation returns 1e+50.
   */
  virtual double getInfBound() const;

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
   * <tt>W</tt> to be evaluated. Same for preconditioner <tt>M</tt>
   *
   * The default implementation returns <tt>return.get()==NULL</tt>
   * (i.e. implicit solvers are not supported by default).
   */
  virtual Teuchos::RefCountPtr<Epetra_Operator> create_W() const;
  virtual Teuchos::RefCountPtr<EpetraExt::ModelEvaluator::Preconditioner> create_WPrec() const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<Epetra_Operator> create_DfDp_op(int l) const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<Epetra_Operator> create_DgDx_dot_op(int j) const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<Epetra_Operator> create_DgDx_op(int j) const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<Epetra_Operator> create_DgDp_op( int j, int l ) const;

  //@}

  /** \name Computational functions */
  //@{

  /** \brief . */
  virtual InArgs createInArgs() const = 0;

  /** \brief . */
  virtual OutArgs createOutArgs() const = 0;

  /** \brief . */
  virtual void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const = 0;

#ifdef HAVE_PYTRILINOS
  /** \brief . */
  friend InArgs convertInArgsFromPython(PyObject * source);

  /** \brief . */
  friend OutArgs convertOutArgsFromPython(PyObject * source);
#endif
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
    /** \brief . */
    void setSupports( EInArgs_p_sg arg, int l, bool supports );
    /** \brief . */
    void setSupports( EInArgs_p_mp arg, int l, bool supports );
  };

  /** \brief . */
  class OutArgsSetup : public OutArgs {
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
    void setSupports(EOutArgsDgDx_dot arg, int j, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDx arg, int j, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDp arg, int j, int l, const DerivativeSupport& );
    /** \brief . */
    void setSupports( EOutArgs_g_sg arg, int j, bool supports );
    /** \brief . */
    void setSupports(EOutArgsDfDp_sg arg, int l, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDx_dot_sg arg, int j, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDx_sg arg, int j, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDp_sg arg, int j, int l, const DerivativeSupport& );
    /** \brief . */
    void setSupports( EOutArgs_g_mp arg, int j, bool supports );
    /** \brief . */
    void setSupports(EOutArgsDfDp_mp arg, int l, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDx_dot_mp arg, int j, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDx_mp arg, int j, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDp_mp arg, int j, int l, const DerivativeSupport& );
    /** \brief . */
    void set_W_properties( const DerivativeProperties &properties );
    void set_WPrec_properties( const DerivativeProperties &properties );
    /** \brief . */
    void set_DfDp_properties( int l, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDx_dot_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDx_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDp_properties( int j, int l, const DerivativeProperties &properties );
    /** \brief . */
    void set_DfDp_sg_properties( int l, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDx_dot_sg_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDx_sg_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDp_sg_properties( int j, int l, const DerivativeProperties &properties );
    /** \brief . */
    void set_DfDp_mp_properties( int l, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDx_dot_mp_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDx_mp_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDp_mp_properties( int j, int l, const DerivativeProperties &properties );
  };

  //@}

};

// ////////////////////////////
// Helper functions

/** \brief . */
std::string toString( ModelEvaluator::EDerivativeMultiVectorOrientation orientation );

/** \brief . */
std::string toString( ModelEvaluator::EInArgsMembers inArg );

/** \brief . */
std::string toString( ModelEvaluator::EOutArgsMembers outArg );

/** \brief . */
Teuchos::RefCountPtr<Epetra_Operator>
getLinearOp(
  const std::string &modelEvalDescription,
  const ModelEvaluator::Derivative &deriv,
  const std::string &derivName
  );

/** \brief . */
Teuchos::RefCountPtr<Epetra_MultiVector>
getMultiVector(
  const std::string &modelEvalDescription,
  const ModelEvaluator::Derivative &deriv,
  const std::string &derivName,
  const ModelEvaluator::EDerivativeMultiVectorOrientation mvOrientation
  );

/** \brief . */
Teuchos::RefCountPtr<Epetra_Operator>
get_DfDp_op(
  const int l
  ,const ModelEvaluator::OutArgs &outArgs
  );

/** \brief . */
Teuchos::RefCountPtr<Epetra_MultiVector>
get_DfDp_mv(
  const int l
  ,const ModelEvaluator::OutArgs &outArgs
  );

/** \brief . */
Teuchos::RefCountPtr<Epetra_MultiVector>
get_DgDx_dot_mv(
  const int j
  ,const ModelEvaluator::OutArgs &outArgs
  ,const ModelEvaluator::EDerivativeMultiVectorOrientation mvOrientation
  );

/** \brief . */
Teuchos::RefCountPtr<Epetra_MultiVector>
get_DgDx_mv(
  const int j
  ,const ModelEvaluator::OutArgs &outArgs
  ,const ModelEvaluator::EDerivativeMultiVectorOrientation mvOrientation
  );

/** \brief . */
Teuchos::RefCountPtr<Epetra_MultiVector>
get_DgDp_mv(
  const int j
  ,const int l
  ,const ModelEvaluator::OutArgs &outArgs
  ,const ModelEvaluator::EDerivativeMultiVectorOrientation mvOrientation
  );

// ///////////////////////////
// Inline Functions

//
// ModelEvaluator::InArgs
//

inline
std::string ModelEvaluator::InArgs::modelEvalDescription() const
{ return modelEvalDescription_; }

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
void ModelEvaluator::InArgs::set_x_dot_poly( const Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > &x_dot_poly )
{ assert_supports(IN_ARG_x_dot_poly); x_dot_poly_ = x_dot_poly; }

inline 
Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> >
ModelEvaluator::InArgs::get_x_dot_poly() const
{ assert_supports(IN_ARG_x_dot_poly); return x_dot_poly_; }

inline 
void ModelEvaluator::InArgs::set_x_poly( const Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > &x_poly )
{ assert_supports(IN_ARG_x_poly); x_poly_ = x_poly; }

inline 
Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> >
ModelEvaluator::InArgs::get_x_poly() const
{ assert_supports(IN_ARG_x_poly); return x_poly_; }

inline 
void ModelEvaluator::InArgs::set_x_dot_sg( const ModelEvaluator::InArgs::sg_const_vector_t &x_dot_sg )
{ assert_supports(IN_ARG_x_dot_sg); x_dot_sg_ = x_dot_sg; }

inline 
ModelEvaluator::InArgs::sg_const_vector_t
ModelEvaluator::InArgs::get_x_dot_sg() const
{ assert_supports(IN_ARG_x_dot_sg); return x_dot_sg_; }

inline 
void ModelEvaluator::InArgs::set_x_dot_mp( const ModelEvaluator::mp_const_vector_t &x_dot_mp )
{ assert_supports(IN_ARG_x_dot_mp); x_dot_mp_ = x_dot_mp; }

inline 
ModelEvaluator::mp_const_vector_t
ModelEvaluator::InArgs::get_x_dot_mp() const
{ assert_supports(IN_ARG_x_dot_mp); return x_dot_mp_; }

inline 
void ModelEvaluator::InArgs::set_x_sg( const ModelEvaluator::InArgs::sg_const_vector_t &x_sg )
{ assert_supports(IN_ARG_x_sg); x_sg_ = x_sg; }

inline 
ModelEvaluator::InArgs::sg_const_vector_t
ModelEvaluator::InArgs::get_x_sg() const
{ assert_supports(IN_ARG_x_sg); return x_sg_; }

inline 
void ModelEvaluator::InArgs::set_x_mp( const ModelEvaluator::mp_const_vector_t &x_mp )
{ assert_supports(IN_ARG_x_mp); x_mp_ = x_mp; }

inline 
ModelEvaluator::mp_const_vector_t
ModelEvaluator::InArgs::get_x_mp() const
{ assert_supports(IN_ARG_x_mp); return x_mp_; }

inline
void ModelEvaluator::InArgs::set_p( int l, const Teuchos::RefCountPtr<const Epetra_Vector> &p_l )
{ assert_l(l); p_[l] = p_l; }

inline
Teuchos::RefCountPtr<const Epetra_Vector> ModelEvaluator::InArgs::get_p(int l) const
{ assert_l(l); return p_[l]; }

inline
void ModelEvaluator::InArgs::set_p_sg( int l, 
				       const ModelEvaluator::InArgs::sg_const_vector_t &p_sg_l )
{ assert_supports(IN_ARG_p_sg, l); p_sg_[l] = p_sg_l; }

inline
ModelEvaluator::InArgs::sg_const_vector_t 
ModelEvaluator::InArgs::get_p_sg(int l) const
{ assert_supports(IN_ARG_p_sg, l); return p_sg_[l]; }

inline
void ModelEvaluator::InArgs::set_p_mp( int l, 
				       const ModelEvaluator::mp_const_vector_t &p_mp_l )
{ assert_supports(IN_ARG_p_mp, l); p_mp_[l] = p_mp_l; }

inline
ModelEvaluator::mp_const_vector_t 
ModelEvaluator::InArgs::get_p_mp(int l) const
{ assert_supports(IN_ARG_p_mp, l); return p_mp_[l]; }

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
void ModelEvaluator::InArgs::set_sg_basis( const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis )
{ assert_supports(IN_ARG_sg_basis); sg_basis_ = basis; }

inline
Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >
ModelEvaluator::InArgs::get_sg_basis() const
{ assert_supports(IN_ARG_sg_basis); return sg_basis_; }

inline
void ModelEvaluator::InArgs::set_sg_quadrature( const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& quad )
{ assert_supports(IN_ARG_sg_quadrature); sg_quad_ = quad; }

inline
Teuchos::RCP<const Stokhos::Quadrature<int,double> >
ModelEvaluator::InArgs::get_sg_quadrature() const
{ assert_supports(IN_ARG_sg_quadrature); return sg_quad_; }

inline
void ModelEvaluator::InArgs::set_sg_expansion( const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double,Stokhos::StandardStorage<int,double> > >& exp )
{ assert_supports(IN_ARG_sg_expansion); sg_exp_ = exp; }

inline
Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double,Stokhos::StandardStorage<int,double> > >
ModelEvaluator::InArgs::get_sg_expansion() const
{ assert_supports(IN_ARG_sg_expansion); return sg_exp_; }

inline
void ModelEvaluator::InArgs::_setModelEvalDescription( const std::string &new_modelEvalDescription )
{
  modelEvalDescription_ = new_modelEvalDescription;
}

inline
void ModelEvaluator::InArgs::_set_Np(int new_Np)
{
  p_.resize(new_Np);
  p_sg_.resize(new_Np);
  p_mp_.resize(new_Np);
  supports_p_sg_.resize(new_Np);
  supports_p_mp_.resize(new_Np);
}

//
// ModelEvaluator::OutArgs
//

inline
std::string ModelEvaluator::OutArgs::modelEvalDescription() const
{ return modelEvalDescription_; }

inline
int ModelEvaluator::OutArgs::Np() const
{
  return DfDp_.size();
}

inline
int ModelEvaluator::OutArgs::Ng() const
{ 
  return g_.size();
}

inline
void ModelEvaluator::OutArgs::set_f( const Evaluation<Epetra_Vector> &f ) { f_ = f; }

inline
ModelEvaluator::Evaluation<Epetra_Vector>
ModelEvaluator::OutArgs::get_f() const { return f_; }

inline
void ModelEvaluator::OutArgs::set_g( int j, const Evaluation<Epetra_Vector> &g_j )
{
  assert_j(j);
  g_[j] = g_j;
}

inline
ModelEvaluator::Evaluation<Epetra_Vector>
ModelEvaluator::OutArgs::get_g(int j) const
{
  assert_j(j);
  return g_[j];
}

inline
void ModelEvaluator::OutArgs::set_g_sg( int j, const sg_vector_t &g_sg_j )
{
  assert_supports(OUT_ARG_g_sg, j);
  g_sg_[j] = g_sg_j;
}

inline
ModelEvaluator::OutArgs::sg_vector_t
ModelEvaluator::OutArgs::get_g_sg(int j) const
{
  assert_supports(OUT_ARG_g_sg, j);
  return g_sg_[j];
}

inline
void ModelEvaluator::OutArgs::set_g_mp( int j, const mp_vector_t &g_mp_j )
{
  assert_supports(OUT_ARG_g_mp, j);
  g_mp_[j] = g_mp_j;
}

inline
ModelEvaluator::mp_vector_t
ModelEvaluator::OutArgs::get_g_mp(int j) const
{
  assert_supports(OUT_ARG_g_mp, j);
  return g_mp_[j];
}

inline
void ModelEvaluator::OutArgs::set_W( const Teuchos::RefCountPtr<Epetra_Operator> &W ) { W_ = W; }
inline
void ModelEvaluator::OutArgs::set_WPrec( const Teuchos::RefCountPtr<Epetra_Operator> &WPrec ) { WPrec_ = WPrec; }

inline
Teuchos::RefCountPtr<Epetra_Operator> ModelEvaluator::OutArgs::get_W() const { return W_; }
inline
Teuchos::RefCountPtr<Epetra_Operator> ModelEvaluator::OutArgs::get_WPrec() const { return WPrec_; }

inline
ModelEvaluator::DerivativeProperties ModelEvaluator::OutArgs::get_W_properties() const
{ return W_properties_; }
inline
ModelEvaluator::DerivativeProperties ModelEvaluator::OutArgs::get_WPrec_properties() const
{ return WPrec_properties_; }

inline
void ModelEvaluator::OutArgs::set_DfDp( int l, const Derivative &DfDp_l )
{
  assert_supports(OUT_ARG_DfDp,l);
  DfDp_[l] = DfDp_l;
}

inline
ModelEvaluator::Derivative
ModelEvaluator::OutArgs::get_DfDp(int l) const
{
  assert_supports(OUT_ARG_DfDp,l);
  return DfDp_[l];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DfDp_properties(int l) const
{
  assert_supports(OUT_ARG_DfDp,l);
  return DfDp_properties_[l];
}

inline
void ModelEvaluator::OutArgs::set_DfDp_sg( int l, const SGDerivative &DfDp_sg_l )
{
  assert_supports(OUT_ARG_DfDp_sg,l);
  DfDp_sg_[l] = DfDp_sg_l;
}

inline
ModelEvaluator::SGDerivative
ModelEvaluator::OutArgs::get_DfDp_sg(int l) const
{
  assert_supports(OUT_ARG_DfDp_sg,l);
  return DfDp_sg_[l];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DfDp_sg_properties(int l) const
{
  assert_supports(OUT_ARG_DfDp_sg,l);
  return DfDp_sg_properties_[l];
}

inline
void ModelEvaluator::OutArgs::set_DfDp_mp( int l, const MPDerivative &DfDp_mp_l )
{
  assert_supports(OUT_ARG_DfDp_mp,l);
  DfDp_mp_[l] = DfDp_mp_l;
}

inline
ModelEvaluator::MPDerivative
ModelEvaluator::OutArgs::get_DfDp_mp(int l) const
{
  assert_supports(OUT_ARG_DfDp_mp,l);
  return DfDp_mp_[l];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DfDp_mp_properties(int l) const
{
  assert_supports(OUT_ARG_DfDp_mp,l);
  return DfDp_mp_properties_[l];
}

inline
void ModelEvaluator::OutArgs::set_DgDx_dot( int j, const Derivative &DgDx_dot_j )
{
  assert_supports(OUT_ARG_DgDx_dot,j);
  DgDx_dot_[j] = DgDx_dot_j;
}

inline
ModelEvaluator::Derivative
ModelEvaluator::OutArgs::get_DgDx_dot(int j) const
{
  assert_supports(OUT_ARG_DgDx_dot,j);
  return DgDx_dot_[j];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DgDx_dot_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx_dot,j);
  return DgDx_dot_properties_[j];
}

inline
void ModelEvaluator::OutArgs::set_DgDx_dot_sg( int j, const SGDerivative &DgDx_dot_sg_j )
{
  assert_supports(OUT_ARG_DgDx_dot_sg,j);
  DgDx_dot_sg_[j] = DgDx_dot_sg_j;
}

inline
ModelEvaluator::SGDerivative
ModelEvaluator::OutArgs::get_DgDx_dot_sg(int j) const
{
  assert_supports(OUT_ARG_DgDx_dot_sg,j);
  return DgDx_dot_sg_[j];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DgDx_dot_sg_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx_dot_sg,j);
  return DgDx_dot_sg_properties_[j];
}

inline
void ModelEvaluator::OutArgs::set_DgDx_dot_mp( int j, const MPDerivative &DgDx_dot_mp_j )
{
  assert_supports(OUT_ARG_DgDx_dot_mp,j);
  DgDx_dot_mp_[j] = DgDx_dot_mp_j;
}

inline
ModelEvaluator::MPDerivative
ModelEvaluator::OutArgs::get_DgDx_dot_mp(int j) const
{
  assert_supports(OUT_ARG_DgDx_dot_mp,j);
  return DgDx_dot_mp_[j];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DgDx_dot_mp_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx_dot_mp,j);
  return DgDx_dot_mp_properties_[j];
}

inline
void ModelEvaluator::OutArgs::set_DgDx( int j, const Derivative &DgDx_j )
{
  assert_supports(OUT_ARG_DgDx,j);
  DgDx_[j] = DgDx_j;
}

inline
ModelEvaluator::Derivative
ModelEvaluator::OutArgs::get_DgDx(int j) const
{
  assert_supports(OUT_ARG_DgDx,j);
  return DgDx_[j];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DgDx_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx,j);
  return DgDx_properties_[j];
}

inline
void ModelEvaluator::OutArgs::set_DgDx_sg( int j, const SGDerivative &DgDx_sg_j )
{
  assert_supports(OUT_ARG_DgDx_sg,j);
  DgDx_sg_[j] = DgDx_sg_j;
}

inline
ModelEvaluator::SGDerivative
ModelEvaluator::OutArgs::get_DgDx_sg(int j) const
{
  assert_supports(OUT_ARG_DgDx_sg,j);
  return DgDx_sg_[j];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DgDx_sg_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx_sg,j);
  return DgDx_sg_properties_[j];
}

inline
void ModelEvaluator::OutArgs::set_DgDx_mp( int j, const MPDerivative &DgDx_mp_j )
{
  assert_supports(OUT_ARG_DgDx_mp,j);
  DgDx_mp_[j] = DgDx_mp_j;
}

inline
ModelEvaluator::MPDerivative
ModelEvaluator::OutArgs::get_DgDx_mp(int j) const
{
  assert_supports(OUT_ARG_DgDx_mp,j);
  return DgDx_mp_[j];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DgDx_mp_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx_mp,j);
  return DgDx_mp_properties_[j];
}

inline
void ModelEvaluator::OutArgs::set_DgDp( int j, int l, const Derivative &DgDp_j_l )
{
  assert_supports(OUT_ARG_DgDp,j,l);
  DgDp_[ j*Np() + l ] = DgDp_j_l;
}

inline
ModelEvaluator::Derivative
ModelEvaluator::OutArgs::get_DgDp(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp,j,l);
  return DgDp_[ j*Np() + l ];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DgDp_properties(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp,j,l);
  return DgDp_properties_[ j*Np() + l ];
}

inline
void ModelEvaluator::OutArgs::set_DgDp_sg( int j, int l, const SGDerivative &DgDp_sg_j_l )
{
  assert_supports(OUT_ARG_DgDp_sg,j,l);
  DgDp_sg_[ j*Np() + l ] = DgDp_sg_j_l;
}

inline
ModelEvaluator::SGDerivative
ModelEvaluator::OutArgs::get_DgDp_sg(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp_sg,j,l);
  return DgDp_sg_[ j*Np() + l ];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DgDp_sg_properties(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp_sg,j,l);
  return DgDp_sg_properties_[ j*Np() + l ];
}

inline
void ModelEvaluator::OutArgs::set_DgDp_mp( int j, int l, const MPDerivative &DgDp_mp_j_l )
{
  assert_supports(OUT_ARG_DgDp_mp,j,l);
  DgDp_mp_[ j*Np() + l ] = DgDp_mp_j_l;
}

inline
ModelEvaluator::MPDerivative
ModelEvaluator::OutArgs::get_DgDp_mp(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp_mp,j,l);
  return DgDp_mp_[ j*Np() + l ];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DgDp_mp_properties(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp_mp,j,l);
  return DgDp_mp_properties_[ j*Np() + l ];
}

inline
void ModelEvaluator::OutArgs::set_f_poly( const Teuchos::RefCountPtr<Teuchos::Polynomial<Epetra_Vector> > &f_poly )
{ f_poly_ = f_poly; }

inline
Teuchos::RefCountPtr<Teuchos::Polynomial<Epetra_Vector> >
ModelEvaluator::OutArgs::get_f_poly() const
{ return f_poly_; }

inline
void ModelEvaluator::OutArgs::set_f_sg( const ModelEvaluator::OutArgs::sg_vector_t& f_sg )
{ f_sg_ = f_sg; }

inline
ModelEvaluator::OutArgs::sg_vector_t
ModelEvaluator::OutArgs::get_f_sg() const
{ return f_sg_; }

inline
void ModelEvaluator::OutArgs::set_W_sg( const ModelEvaluator::OutArgs::sg_operator_t& W_sg ) { W_sg_ = W_sg; }

inline
ModelEvaluator::OutArgs::sg_operator_t ModelEvaluator::OutArgs::get_W_sg() const { return W_sg_; }

inline
void ModelEvaluator::OutArgs::set_f_mp( const ModelEvaluator::mp_vector_t& f_mp )
{ f_mp_ = f_mp; }

inline
ModelEvaluator::mp_vector_t
ModelEvaluator::OutArgs::get_f_mp() const
{ return f_mp_; }

inline
void ModelEvaluator::OutArgs::set_W_mp( const ModelEvaluator::mp_operator_t& W_mp ) { W_mp_ = W_mp; }

inline
ModelEvaluator::mp_operator_t ModelEvaluator::OutArgs::get_W_mp() const { return W_mp_; }

//
// ModelEvaluator::InArgsSetup
//

inline
void ModelEvaluator::InArgsSetup::setModelEvalDescription( const std::string &new_modelEvalDescription )
{
  this->_setModelEvalDescription(new_modelEvalDescription);
}

inline
void ModelEvaluator::InArgsSetup::set_Np(int new_Np)
{ this->_set_Np(new_Np); }

inline
void ModelEvaluator::InArgsSetup::setSupports( EInArgsMembers arg, bool new_supports )
{ this->_setSupports(arg,new_supports); }

inline
void ModelEvaluator::InArgsSetup::setSupports( EInArgs_p_sg arg, int l, bool new_supports )
{ this->_setSupports(arg,l,new_supports); }

inline
void ModelEvaluator::InArgsSetup::setSupports( EInArgs_p_mp arg, int l, bool new_supports )
{ this->_setSupports(arg,l,new_supports); }

//
// ModelEvaluator::OutArgsSetup
//

inline
void ModelEvaluator::OutArgsSetup::setModelEvalDescription( const std::string &new_modelEvalDescription )
{
  this->_setModelEvalDescription(new_modelEvalDescription);
}

inline
void ModelEvaluator::OutArgsSetup::set_Np_Ng(int new_Np, int new_Ng)
{ this->_set_Np_Ng(new_Np,new_Ng); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsMembers arg, bool new_supports )
{ this->_setSupports(arg,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDfDp arg, int l, const DerivativeSupport& new_supports )
{ this->_setSupports(arg,l,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDgDx_dot arg, int j, const DerivativeSupport& new_supports )
{ this->_setSupports(arg,j,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDgDx arg, int j, const DerivativeSupport& new_supports )
{ this->_setSupports(arg,j,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDgDp arg, int j, int l, const DerivativeSupport& new_supports )
{ this->_setSupports(arg,j,l,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgs_g_sg arg, int j, bool new_supports )
{ this->_setSupports(arg,j,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDfDp_sg arg, int l, const DerivativeSupport& new_supports )
{ this->_setSupports(arg,l,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDgDx_dot_sg arg, int j, const DerivativeSupport& new_supports )
{ this->_setSupports(arg,j,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDgDx_sg arg, int j, const DerivativeSupport& new_supports )
{ this->_setSupports(arg,j,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDgDp_sg arg, int j, int l, const DerivativeSupport& new_supports )
{ this->_setSupports(arg,j,l,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgs_g_mp arg, int j, bool new_supports )
{ this->_setSupports(arg,j,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDfDp_mp arg, int l, const DerivativeSupport& new_supports )
{ this->_setSupports(arg,l,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDgDx_dot_mp arg, int j, const DerivativeSupport& new_supports )
{ this->_setSupports(arg,j,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDgDx_mp arg, int j, const DerivativeSupport& new_supports )
{ this->_setSupports(arg,j,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDgDp_mp arg, int j, int l, const DerivativeSupport& new_supports )
{ this->_setSupports(arg,j,l,new_supports); }

inline
void ModelEvaluator::OutArgsSetup::set_W_properties( const DerivativeProperties &properties )
{ this->_set_W_properties(properties); }
inline
void ModelEvaluator::OutArgsSetup::set_WPrec_properties( const DerivativeProperties &properties )
{ this->_set_WPrec_properties(properties); }

inline
void ModelEvaluator::OutArgsSetup::set_DfDp_properties( int l, const DerivativeProperties &properties )
{
  this->_set_DfDp_properties(l,properties);
}

inline
void ModelEvaluator::OutArgsSetup::set_DgDx_dot_properties( int j, const DerivativeProperties &properties )
{
  this->_set_DgDx_dot_properties(j,properties);
}

inline
void ModelEvaluator::OutArgsSetup::set_DgDx_properties( int j, const DerivativeProperties &properties )
{
  this->_set_DgDx_properties(j,properties);
}

inline
void ModelEvaluator::OutArgsSetup::set_DgDp_properties( int j, int l, const DerivativeProperties &properties )
{
  this->_set_DgDp_properties(j,l,properties);
}

inline
void ModelEvaluator::OutArgsSetup::set_DfDp_sg_properties( int l, const DerivativeProperties &properties )
{
  this->_set_DfDp_sg_properties(l,properties);
}

inline
void ModelEvaluator::OutArgsSetup::set_DgDx_dot_sg_properties( int j, const DerivativeProperties &properties )
{
  this->_set_DgDx_dot_sg_properties(j,properties);
}

inline
void ModelEvaluator::OutArgsSetup::set_DgDx_sg_properties( int j, const DerivativeProperties &properties )
{
  this->_set_DgDx_sg_properties(j,properties);
}

inline
void ModelEvaluator::OutArgsSetup::set_DgDp_sg_properties( int j, int l, const DerivativeProperties &properties )
{
  this->_set_DgDp_sg_properties(j,l,properties);
}

inline
void ModelEvaluator::OutArgsSetup::set_DfDp_mp_properties( int l, const DerivativeProperties &properties )
{
  this->_set_DfDp_mp_properties(l,properties);
}

inline
void ModelEvaluator::OutArgsSetup::set_DgDx_dot_mp_properties( int j, const DerivativeProperties &properties )
{
  this->_set_DgDx_dot_mp_properties(j,properties);
}

inline
void ModelEvaluator::OutArgsSetup::set_DgDx_mp_properties( int j, const DerivativeProperties &properties )
{
  this->_set_DgDx_mp_properties(j,properties);
}

inline
void ModelEvaluator::OutArgsSetup::set_DgDp_mp_properties( int j, int l, const DerivativeProperties &properties )
{
  this->_set_DgDp_mp_properties(j,l,properties);
}

} // namespace EpetraExt

#endif // EPETRA_EXT_MODEL_EVALUATOR_HPP
