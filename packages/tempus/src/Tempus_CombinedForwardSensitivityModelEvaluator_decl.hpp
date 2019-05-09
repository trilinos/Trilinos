// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_CombinedForwardSensitivityModelEvaluator_decl_hpp
#define Tempus_CombinedForwardSensitivityModelEvaluator_decl_hpp

#include "Tempus_SensitivityModelEvaluatorBase.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"

namespace Tempus {

/** \brief Transform a ModelEvaluator's sensitivity equations to its residual
 *
 * This class wraps a given ModelEvalutor encapsulating f(x,p) and creates
 * a new "residual" for it and the forward sensitivity equations:
 *        F(X) = [                f(x,p)             ] = 0
 *               [ (df/dx)(x,p) * dx/dp + df/dp(x,p) ]
 * where X = [ x; dx/dp ] (transient terms supressed for simplicity). This model
 * evaluator can then be handed to a regular (non)linear solver to compute X.
 * Note that even though these equations are linear in X, it is not necessarily
 * the case that the underlying model evaluator accurately computes df/dx in its
 * evaluation of W.  Therefore this model evaluator can optionally reinterpret
 * the model's df/dp out-arg as (df/dx)(x,p) * dx/dp + df/dp(x,p) where dx/dp is
 * passed as another parameter (product) vector (encapsulated in the
 * Thyra::DefaultMultiVectorProductVector).  This is not standard model
 * evaluator behavior, but is useful for models where W is only an approximation
 * to df/dx and/or the model is capable of directly computing
 * (df/dx)(x,p) * dx/dp + df/dp(x,p).
 */
template <typename Scalar>
class CombinedForwardSensitivityModelEvaluator :
    public Thyra::StateFuncModelEvaluatorBase<Scalar>,
    public SensitivityModelEvaluatorBase<Scalar> {
public:
  typedef Thyra::VectorBase<Scalar>  Vector;
  typedef Thyra::MultiVectorBase<Scalar>  MultiVector;

  //! Constructor
  /*!
   * The optionally supplied parameter list supports the following options:
   * <ul>
   *   <li> "Use DfDp as Tangent" (default:  false) Reinterpret the df/dp
   *        out-arg as the tangent vector (df/dx)(x,p) * dx/dp + df/dp(x,p)
   *        as described above.  If it is false, this implementation will
   *        compute the tangent through several calls to the models
   *        evalModel().
   *   <li> "Sensitivity Parameter Index" (default: 0) Model evaluator
   *        parameter index for which sensitivities will be computed.
   *   <li> "Sensitivity X Tangent Index" (default: 1) If "Use DfDp as Tangent"
   *        is true, the model evaluator parameter index for passing dx/dp
   *        as a Thyra::DefaultMultiVectorProductVector.
   *   <li> "Sensitivity X-Dot Tangent Index" (default: 2) If
   *        "Use DfDp as Tangent" is true, the model evaluator parameter index
   *        for passing dx_dot/dp as a Thyra::DefaultMultiVectorProductVector.
   *   <li> "Sensitivity X-Dot-Dot Tangent Index" (default: 3) If
   *        "Use DfDp as Tangent" is true, the model evaluator parameter index
   *        for passing dx_dot_dot/dp as a
   *        Thyra::DefaultMultiVectorProductVector (if the model supports
   *        x_dot_dot).
   * </ul>
   */
  CombinedForwardSensitivityModelEvaluator(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & model,
    const Teuchos::RCP<const Teuchos::ParameterList>& pList = Teuchos::null,
    const Teuchos::RCP<MultiVector>& dxdp_init = Teuchos::null,
    const Teuchos::RCP<MultiVector>& dx_dotdp_init = Teuchos::null,
    const Teuchos::RCP<MultiVector>& dx_dotdot_dp_init = Teuchos::null);

  /** \name Public functions overridden from SensitivityModelEvaulator. */
  //@{

  //! Get the underlying model 'f'
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getForwardModel() const
  { return model_; }

  //@}

  //! Get sensitivity parameter index
  int getSensitivityParamIndex() const { return p_index_; }

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int p) const;

  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int p) const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;

  Teuchos::ArrayView<const std::string> get_g_names(int j) const;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_DgDx_dot_op(int j) const;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_DgDx_op(int j) const;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_DgDp_op(int j, int l) const;

  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
  get_W_factory() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;

  //@}

  static Teuchos::RCP<const Teuchos::ParameterList> getValidParameters();

private:

  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;


  Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model_;
  Teuchos::RCP<MultiVector> dxdp_init_;
  Teuchos::RCP<MultiVector> dx_dotdp_init_;
  Teuchos::RCP<MultiVector> dx_dotdotdp_init_;
  int p_index_;
  int x_tangent_index_;
  int xdot_tangent_index_;
  int xdotdot_tangent_index_;
  bool use_dfdp_as_tangent_;

  int num_param_;
  Teuchos::RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > dxdp_space_;
  Teuchos::RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > x_dxdp_space_;
  Teuchos::RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > dfdp_space_;
  Teuchos::RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > f_dfdp_space_;

  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar> > my_dfdx_;
  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar> > my_dfdxdot_;
  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar> > my_dfdxdotdot_;
};

} // namespace Tempus

#endif
