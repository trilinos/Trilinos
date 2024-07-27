// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SIMPLE_2D_TPETRA_MODEL_EVALUATOR_HPP
#define SIMPLE_2D_TPETRA_MODEL_EVALUATOR_HPP


#include "Simple2DTpetraModelEvaluator_decl.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_DefaultComm.hpp"


// Constructors/Initializers/Accessors


template<class Scalar>
Simple2DTpetraModelEvaluator<Scalar>::Simple2DTpetraModelEvaluator()
  : d_(0.0)
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using Thyra::VectorBase;
  using Thyra::createMember;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  const int dim = 2;

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

  const RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  const RCP<const Tpetra::Map<> > map = rcp(new Tpetra::Map<> (dim, 0, comm));

  typedef Tpetra::Map<>::global_ordinal_type GO;

  W_op_graph_ = rcp(new Tpetra::CrsGraph<> (map, dim));
  W_op_graph_->insertGlobalIndices(0, tuple<GO>(0, 1)());
  W_op_graph_->insertGlobalIndices(1, tuple<GO>(0, 1)());
  W_op_graph_->fillComplete();

  p_.resize(dim, ST::zero());

  x0_ = rcp(new Tpetra::Vector<Scalar>(map));
  x0_->putScalar(ST::zero());

  //
  // C) Create the Thyra wrapped objects
  //

  x_space_ = Thyra::createVectorSpace<Scalar>(map);
  f_space_ = x_space_;

  nominalValues_ = inArgs;
  nominalValues_.set_x(Thyra::createVector(x0_, x_space_));

  //
  // D) Set initial values through interface functions
  //

  set_d(10.0);
  set_p(Teuchos::tuple<Scalar>(2.0, 0.0)());
  set_x0(Teuchos::tuple<Scalar>(1.0, 1.0)());

}


template<class Scalar>
void Simple2DTpetraModelEvaluator<Scalar>::set_d(const Scalar &d)
{
  d_ = d;
}


template<class Scalar>
void Simple2DTpetraModelEvaluator<Scalar>::set_p(const Teuchos::ArrayView<const Scalar> &p)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(p_.size(), p.size());
#endif
  p_().assign(p);
}


template<class Scalar>
void Simple2DTpetraModelEvaluator<Scalar>::set_x0(const Teuchos::ArrayView<const Scalar> &x0_in)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(x_space_->dim(), x0_in.size());
#endif
  x0_->get1dViewNonConst()().assign(x0_in);
}


// Public functions overridden from ModelEvaulator


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Simple2DTpetraModelEvaluator<Scalar>::get_x_space() const
{
  return x_space_;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Simple2DTpetraModelEvaluator<Scalar>::get_f_space() const
{
  return f_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Simple2DTpetraModelEvaluator<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Simple2DTpetraModelEvaluator<Scalar>::create_W_op() const
{
  return Thyra::createLinearOp(
    Teuchos::RCP<Tpetra::Operator<Scalar,int> >(
      Teuchos::rcp(new Tpetra::CrsMatrix<Scalar,int>(W_op_graph_))
      )
    );
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Simple2DTpetraModelEvaluator<Scalar>::createInArgs() const
{
  return prototypeInArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Simple2DTpetraModelEvaluator<Scalar>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}


template<class Scalar>
void Simple2DTpetraModelEvaluator<Scalar>::evalModelImpl(
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
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Thyra::TpetraOperatorVectorExtraction<Scalar> ConverterT;

  const RCP<const Tpetra::Vector<Scalar, int> > x_vec =
    ConverterT::getConstTpetraVector(inArgs.get_x());
  const ArrayRCP<const Scalar> x = x_vec->get1dView();

  const RCP<Tpetra::Vector<Scalar, int> > f_vec =
    ConverterT::getTpetraVector(outArgs.get_f());

  const RCP<Tpetra::CrsMatrix<Scalar, int> > W =
    rcp_dynamic_cast<Tpetra::CrsMatrix<Scalar,int> >(
      ConverterT::getTpetraOperator(outArgs.get_W_op()),
      true
      );

  if (nonnull(f_vec)) {
    const ArrayRCP<Scalar> f = f_vec->get1dViewNonConst();
    f[0] = x[0] + x[1]*x[1] - p_[0];
    f[1] = d_ * (x[0]*x[0] -x[1] - p_[1]);
  }

  if (nonnull(W)) {
    W->setAllToScalar(ST::zero());
    W->sumIntoGlobalValues(0, tuple<GO>(0, 1), tuple<Scalar>(1.0, 2.0*x[1]));
    W->sumIntoGlobalValues(1, tuple<GO>(0, 1), tuple<Scalar>(2.0*d_*x[0], -d_));
  }

}


#endif // SIMPLE_2D_TPETRA_MODEL_EVALUATOR_HPP
