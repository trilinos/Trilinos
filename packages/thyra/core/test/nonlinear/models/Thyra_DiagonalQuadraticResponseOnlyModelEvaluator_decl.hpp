#ifndef THYRA_DIAGONAL_QUADRATIC_RESPONSE_ONLY_MODEL_EVALUATOR_DECL_HPP
#define THYRA_DIAGONAL_QUADRATIC_RESPONSE_ONLY_MODEL_EVALUATOR_DECL_HPP


#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"
#include "Teuchos_Comm.hpp"


namespace Thyra {


/** \brief Simple parallel response-only ModelEvaluator.
 *
 * The representation of the model in coefficient form is:

 \verbatim

   g(p) = 0.5 * sum( diag[i] * (p[i] - ps[i])^2, i=0...n-1 )
          + 0.5 * nonlinearTermFactor * sum( (p[i] - ps[i])^3, i=0...n-1 )
          + g_offset

 \endverbatim
 
 * In vector coefficient form it becomes:

 \verbatim

   g(p) = 0.5 * (p-ps)^T * D * (p-ps) + cubicTerm(p) + g_offset

 \endverbatim

 * where <tt>D = diag(diag)</tt> and <tt>cubitTerm(p)</tt> is given above.
 *
 * This test model implementation also supports the definition of diagonal
 * scalar product implementation.  The impact of this is the definition of new
 * abstract Thyra Euclidean vectors of the form:
 
 \verbatim

   p_bar = E_p * p

   g_grad_bar = E_p * g_grad

 \endverbatim

 * In this notation, we say that the Thyra vectors themselves <tt>p_bar</tt>
 * are vectors in a Euclidean space (from the standpoint of the ANA clients)
 * and that <tt>p</tt> are the coefficient vectors in the new non-Euclidean
 * space <tt>S_p</tt> with scalar product:
 
 \verbatim

   y_bar^T * x_bar = y^T * E_p^T * E_p * x
                   = y^T * S_p * x
                   = <y, x>_p

 \endverbatim

 * The ideas is that by providing this new vector space basis <tt>E_p</tt> and
 * inner product <tt><.,.>_p</tt>, we change the function that the ANA sees.
 *
 * This testing class allows a client to specify the desirable diagonal
 * <tt>D_bar</tt> matrix by passing in the vector <tt>diag_bar</tt>.  From
 * <tt>diag_bar</tt> and <tt>diag</tt>, the diagonal <tt>s_diag</tt> for
 * <tt>S_p = diag(s_diag)</tt> as:

 \verbatim

   s_diag[i] = diag_bar[i] / diag[i]

 \endverbatim

 * The inner product is therefore:

 \verbatim

   scalarProd(y_bar, x_bar) = sum( s_diag[i] * y[i] * x[i], i=0...n-1 )

 \endverbatim
 
 */
template<class Scalar>
class DiagonalQuadraticResponseOnlyModelEvaluator
   : public ResponseOnlyModelEvaluatorBase<Scalar>
{
public:

  //@}

  /** \name Constructors/Initializers/Accessors. */
  //@{

  /** \brief . */
  DiagonalQuadraticResponseOnlyModelEvaluator(
    const int localDim,
    const RCP<const Teuchos::Comm<Ordinal> > &comm = Teuchos::null
    );

  /** \brief Set the solution vector ps . */
  void setSolutionVector(const RCP<const VectorBase<Scalar> > &ps);

  /** \brief Get the solution vector ps . */
  const RCP<const VectorBase<Scalar> >
  getSolutionVector() const;
  
  /** \brief Set the diagonal vector diag. */
  void setDiagonalVector(const RCP<const VectorBase<Scalar> > &diag);
  
  /** \brief Set the diagonal vector diag_bar.
   *
   * NOTE: You must call setDiagonalVector(diag) first in order to set the
   * objective diagonal.
   */
  void setDiagonalBarVector(const RCP<const VectorBase<Scalar> > &diag_bar);

  /** \brief Set nonlinear term factory. */
  void setNonlinearTermFactor(const Scalar &nonlinearTermFactor);

  /** \brief Set offset scalar g_offset . */
  void setScalarOffset(const Scalar &g_offset);

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  int Np() const;
  /** \brief . */
  int Ng() const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > get_g_space(int j) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const ModelEvaluatorBase::OutArgs<Scalar>& outArgs
    ) const;

  //@}

private:

  // //////////////////////
  // Private data members

  int Np_;
  int Ng_;
  RCP<const Teuchos::Comm<Ordinal> > comm_;
  const int localDim_;

  RCP<const VectorSpaceBase<Scalar> > g_space_;

  // Declared non-const so we can change the space in place!
  RCP<VectorSpaceBase<Scalar> > p_space_;

  RCP<const VectorBase<Scalar> > ps_;
  RCP<const VectorBase<Scalar> > diag_;
  Scalar nonlinearTermFactor_;
  Scalar g_offset_;

  RCP<const VectorBase<Scalar> > diag_bar_;
  RCP<const VectorBase<Scalar> > s_bar_;

};


/** \brief Non-member constructor. */
template<class Scalar>
RCP<DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> >
diagonalQuadraticResponseOnlyModelEvaluator(
  const int localDim,
  const RCP<const Teuchos::Comm<Ordinal> > &comm = Teuchos::null
  )
{
  using Teuchos::rcp;
  return rcp(new DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>(localDim, comm));
}


} // namespace Thyra


#endif // THYRA_DIAGONAL_QUADRATIC_RESPONSE_ONLY_MODEL_EVALUATOR_DECL_HPP
