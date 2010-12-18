#ifndef PANZER_MODEL_EVALUATOR_DEF_HPP
#define PANZER_MODEL_EVALUATOR_DEF_HPP

#include "Thyra_TpetraThyraWrappers.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_DefaultComm.hpp"
//#include "Panzer_LinearObject_Factory.hpp"

// Constructors/Initializers/Accessors


template<typename Scalar, typename LO, typename GO, typename NODE>
panzer::ModelEvaluator<Scalar,LO,GO,NODE>::
ModelEvaluator()
  : d_(0.0)
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using Thyra::VectorBase;
  using Thyra::createMember;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  //
  // A) Create the structure for the problem
  //
  
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  prototypeInArgs_ = inArgs;
  
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  prototypeOutArgs_ = outArgs;

  //
  // B) Create the Tpetra objects
  //

//   lo_factory_ = rcp(new panzer::LinearObjectFactory<T>(topology));

//   const RCP<const Tpetra::Map<LO,GO,NODE> > map = lo_factory_->localMap();

//   W_op_graph_ = lo_factory_->localGraph();

//   p_.resize(2, ST::zero());

//   x0_ = rcp(new Tpetra::Vector<Scalar,LO,GO,NODE>(map));
//   x0_->putScalar(ST::zero());

  //
  // C) Create the Thyra wrapped objects
  //

//   x_space_ = Thyra::createVectorSpace<Scalar>(map);
//   f_space_ = x_space_;

//   nominalValues_ = inArgs;
//   nominalValues_.set_x(Thyra::createVector(x0_, x_space_));

  //
  // D) Set initial values through interface functions
  //

//   set_d(10.0);
//   set_p(Teuchos::tuple<Scalar>(2.0, 0.0)());

}


template<typename Scalar, typename LO, typename GO, typename NODE>
void panzer::ModelEvaluator<Scalar,LO,GO,NODE>::set_d(const Scalar &d)
{
  d_ = d;
}


template<typename Scalar, typename LO, typename GO, typename NODE>
void panzer::ModelEvaluator<Scalar,LO,GO,NODE>::
set_p(const Teuchos::ArrayView<const Scalar> &p)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(p_.size(), p.size());
#endif
  p_().assign(p);
}


template<typename Scalar, typename LO, typename GO, typename NODE>
void panzer::ModelEvaluator<Scalar,LO,GO,NODE>::
set_x0(const Teuchos::ArrayView<const Scalar> &x0_in)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(x_space_->dim(), x0_in.size());
#endif
  x0_->get1dViewNonConst()().assign(x0_in);
}


// Public functions overridden from ModelEvaulator


template<typename Scalar, typename LO, typename GO, typename NODE>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
panzer::ModelEvaluator<Scalar,LO,GO,NODE>::get_x_space() const
{
  return x_space_;
}


template<typename Scalar, typename LO, typename GO, typename NODE>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
panzer::ModelEvaluator<Scalar,LO,GO,NODE>::get_f_space() const
{
  return f_space_;
}


template<typename Scalar, typename LO, typename GO, typename NODE>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
panzer::ModelEvaluator<Scalar,LO,GO,NODE>::getNominalValues() const
{
  return nominalValues_;
}


template<typename Scalar, typename LO, typename GO, typename NODE>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
panzer::ModelEvaluator<Scalar,LO,GO,NODE>::create_W_op() const
{
  return Thyra::createLinearOp(
      Teuchos::RCP<Tpetra::Operator<Scalar,LO,GO,NODE> >(
	Teuchos::rcp(new Tpetra::CrsMatrix<Scalar,LO,GO,NODE>(W_op_graph_))
      )
    );
}


template<typename Scalar, typename LO, typename GO, typename NODE>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
panzer::ModelEvaluator<Scalar,LO,GO,NODE>::createInArgs() const
{
  return prototypeInArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<typename Scalar, typename LO, typename GO, typename NODE>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
panzer::ModelEvaluator<Scalar,LO,GO,NODE>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}


template<typename Scalar, typename LO, typename GO, typename NODE>
void panzer::ModelEvaluator<Scalar,LO,GO,NODE>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::Array;
  using Teuchos::tuple;
  using Teuchos::rcp_dynamic_cast;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::TpetraOperatorVectorExtraction<Scalar,LO,GO,NODE> ConverterT;

  const RCP<const Tpetra::Vector<Scalar,LO,GO,NODE> > x_vec =
    ConverterT::getConstTpetraVector(inArgs.get_x());
  const ArrayRCP<const Scalar> x = x_vec->get1dView();

  const RCP<Tpetra::Vector<Scalar,LO,GO,NODE> > f_vec =
    ConverterT::getTpetraVector(outArgs.get_f());

  const RCP<Tpetra::CrsMatrix<Scalar,LO,GO,NODE> > W =
    rcp_dynamic_cast<Tpetra::CrsMatrix<Scalar,LO,GO,NODE> >(
      ConverterT::getTpetraOperator(outArgs.get_W_op()),
      true
      );

  if (nonnull(f_vec)) {
//     const ArrayRCP<Scalar> f = f_vec->get1dViewNonConst();
//     f[0] = x[0] + x[1]*x[1] - p_[0];
//     f[1] = d_ * (x[0]*x[0] -x[1] - p_[1]);
  }

  if (nonnull(W)) {
//     W->setAllToScalar(ST::zero());
//     W->sumIntoGlobalValues(0, tuple<int>(0, 1), tuple<Scalar>(1.0, 2.0*x[1]));
//     W->sumIntoGlobalValues(1, tuple<int>(0, 1), tuple<Scalar>(2.0*d_*x[0], -d_));
  }

}


#endif
