// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Thyra_TpetraExplicitAdjointModelEvaluator_hpp
#define Thyra_TpetraExplicitAdjointModelEvaluator_hpp

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"

namespace Thyra {

/** \brief A model evaluator decorator for computing an explicit adjoint. */
/**
 * This ModelEvaluator only supports computing W_op and forms the adjoint by
 * computing W_op from the underlying ModelEvaluator and then explicitly
 * transposes it.
 *
 * Since it derives from ModelEvaluatorDelegatorBase, it forwards all
 * ModelEvaluator methods to the supplied underlying ModelEvaluator except
 * those overloaded here.
 */
  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal=LocalOrdinal, typename Node=Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class TpetraExplicitAdjointModelEvaluator :
    public ModelEvaluatorDelegatorBase<Scalar>{
public:

  //! Constructor
  TpetraExplicitAdjointModelEvaluator(
    const RCP<const ModelEvaluator<Scalar> >& model) :
    ModelEvaluatorDelegatorBase<Scalar>(model) {}

  //! Constructor
  TpetraExplicitAdjointModelEvaluator(
    const RCP<ModelEvaluator<Scalar> >& model) :
    ModelEvaluatorDelegatorBase<Scalar>(model) {}

  //! Destructor
  virtual ~TpetraExplicitAdjointModelEvaluator() = default;

  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const
  {
    // This ME should use the same InArgs as the underlying model.  However
    // we can't just use it's InArgs directly because the description won't
    // match (which is checked in debug builds).  Instead create a new
    // InArgsSetup initialized by the underlying model's createInArgs() and
    // set the description appropriately.
    ModelEvaluatorBase::InArgsSetup<Scalar> inArgs =
      this->getUnderlyingModel()->createInArgs();
    inArgs.setModelEvalDescription(this->description());
    return inArgs;
  }

  RCP<LinearOpBase<Scalar> > create_W_op() const {
    typedef TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> TLO;
    typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> TO;
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> TCM;
    typedef Tpetra::RowMatrixTransposer<Scalar,LocalOrdinal,GlobalOrdinal,Node> TRMT;
    if (thyra_fwd_op == Teuchos::null)
      thyra_fwd_op = this->getUnderlyingModel()->create_W_op();
    RCP<TLO> thyra_tpetra_fwd_op =
      Teuchos::rcp_dynamic_cast<TLO>(thyra_fwd_op,true);
    RCP<TO> tpetra_fwd_op = thyra_tpetra_fwd_op->getTpetraOperator();
    RCP<TCM> tpetra_fwd_mat =
      Teuchos::rcp_dynamic_cast<TCM>(tpetra_fwd_op,true);
    TRMT transposer(tpetra_fwd_mat);
    RCP<TCM> tpetra_trans_mat = transposer.createTranspose();
    return tpetraLinearOp(thyra_op->range(), thyra_op->domain(),
                          tpetra_trans_mat);
  }

private:

  mutable RCP< Thyra::LinearOpBase<Scalar> > thyra_fwd_op;

  ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const
  {
    typedef ModelEvaluatorBase MEB;
    MEB::OutArgs<Scalar> model_outArgs =
      this->getUnderlyingModel()->createOutArgs();
    MEB::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(MEB::OUT_ARG_f); //All models must support f, apparently
    outArgs.setSupports(MEB::OUT_ARG_W_op);
    TEUCHOS_ASSERT(model_outArgs.supports(MEB::OUT_ARG_W_op));
    return outArgs;
  }

  void evalModelImpl(
    const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
  {
    typedef ModelEvaluatorBase MEB;
    typedef TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> TLO;
    typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> TO;
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> TCM;
    typedef Tpetra::RowMatrixTransposer<Scalar,LocalOrdinal,GlobalOrdinal,Node> TRMT;

    MEB::OutArgs<Scalar> model_outArgs =
      this->getUnderlyingModel()->createOutArgs();

    if (model_outArgs.supports(MEB::OUT_ARG_W_op) &&
        outArgs.get_W_op() != Teuchos::null) {
      // Compute underlying W_op.
      if (thyra_fwd_op == Teuchos::null)
        thyra_fwd_op = this->getUnderlyingModel()->create_W_op();
      model_outArgs->set_W_op(thyra_fwd_op);
      this->getUnderlyingModel()->evalModel(inArgs, model_outArgs);

      // Transpose W_op
      // Unfortunately, because of the horrendous design of the Tpetra
      // RowMatrixTransposer, this creates a new transposed matrix each time
      RCP<TLO> thyra_tpetra_fwd_op =
        Teuchos::rcp_dynamic_cast<TLO>(thyra_fwd_op,true);
      RCP<TO> tpetra_fwd_op = thyra_tpetra_fwd_op->getTpetraOperator();
      RCP<TCM> tpetra_fwd_mat =
        Teuchos::rcp_dynamic_cast<TCM>(tpetra_fwd_op,true);
      TRMT transposer(tpetra_fwd_mat);
      RCP<TCM> tpetra_trans_mat = transposer.createTranspose();

      // Copy transposed matrix into our outArg
      RCP<LOB> thyra_adj_op = outArgs.get_W_op();
      RCP<TLO> thyra_tpetra_adj_op =
        Teuchos::rcp_dynamic_cast<TLO>(thyra_adj_op,true);
      RCP<TO> tpetra_adj_op = thyra_tpetra_adj_op->getTpetraOperator();
      RCP<TCM> tpetra_adj_mat =
        Teuchos::rcp_dynamic_cast<TCM>(tpetra_adj_op,true);
      *tpetra_adj_mat = *tpetra_trans_mat;
    }
  }

};

template <typename Scalar>
RCP<TpetraExplicitAdjointModelEvaluator<Scalar> >
tpetraExplicitAdjointModelEvaluator(
  const RCP<const ModelEvaluator<Scalar> >& model)
{
  return Teuchos::rcp(new TpetraExplicitAdjointModelEvaluator<Scalar>(model));
}

template <typename Scalar>
RCP<TpetraExplicitAdjointModelEvaluator<Scalar> >
tpetraExplicitAdjointModelEvaluator(
  const RCP<ModelEvaluator<Scalar> >& model)
{
  return Teuchos::rcp(new TpetraExplicitAdjointModelEvaluator<Scalar>(model));
}

} // namespace Thyra

#endif
