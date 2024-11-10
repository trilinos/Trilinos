// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef NOX_THYRA_MODEL_EVALUATOR_HEQ_DEF_HPP
#define NOX_THYRA_MODEL_EVALUATOR_HEQ_DEF_HPP

// Thyra support
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_PreconditionerBase.hpp"

// Epetra support
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

// Nonmember constuctors

template<class Scalar>
Teuchos::RCP<ModelEvaluatorHeq<Scalar> >
modelEvaluatorHeq(const Teuchos::RCP<const Epetra_Comm>& comm,
            const int num_global_elements,
            const Scalar paramC)
{
  return Teuchos::rcp(new ModelEvaluatorHeq<Scalar>(comm,num_global_elements,paramC));
}

// Constructor

template<class Scalar>
ModelEvaluatorHeq<Scalar>::
ModelEvaluatorHeq(const Teuchos::RCP<const Epetra_Comm>& comm,
            const int num_global_elements,
            const Scalar paramC) :
  comm_(comm),
  num_global_elements_(num_global_elements),
  paramC_(paramC),
  showGetInvalidArg_(false)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using ::Thyra::VectorBase;
  typedef ::Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEUCHOS_ASSERT(nonnull(comm_));

  // owned space
  x_map_ = rcp(new Epetra_Map(num_global_elements_,0,*comm_));
  x_space_ = ::Thyra::create_VectorSpace(x_map_);

  // residual space
  f_map_ = x_map_;
  f_space_ = x_space_;

  x0_ = ::Thyra::createMember(x_space_);
  V_S(x0_.ptr(), ST::zero());

  // Create the nodal coordinates
  {
    node_coordinates_.resize(num_global_elements_);
    Scalar dx = 1.0/static_cast<Scalar>(num_global_elements_);
    for (int i=0; i < num_global_elements_; i++)
      node_coordinates_[i] = 0.5*dx + dx*static_cast<Scalar>(i);
  }

  // Compute the kernel of the integral operator, A_heq_
  A_heq_ = rcp(new Epetra_CrsMatrix(::Copy,*x_map_,num_global_elements_));
  int num_local_elements = x_map_->NumMyElements();
  for (int i=0; i<num_local_elements; i++) { // Loop over rows on process
    int row = x_map_->GID(i);
    for(int j=0; j<num_global_elements_; j++) { // Loop over entries in row
      Scalar value = node_coordinates_[row]/(node_coordinates_[row]+node_coordinates_[j]);
      A_heq_->InsertGlobalValues(row,1,&value,&j);
    }
  }
  A_heq_->FillComplete();
  Scalar cc = 0.5*paramC_/static_cast<Scalar>(num_global_elements_);
  A_heq_->Scale(cc);

  ones_ = Teuchos::rcp(new Epetra_Vector(*x_map_));
  ones_->PutScalar(1.0);

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  prototypeInArgs_ = inArgs;

  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  outArgs.setSupports(MEB::OUT_ARG_W_prec);
  prototypeOutArgs_ = outArgs;

  nominalValues_ = inArgs;
  nominalValues_.set_x(x0_);
}

// Initializers/Accessors

template<class Scalar>
void ModelEvaluatorHeq<Scalar>::set_x0(const Teuchos::ArrayView<const Scalar> &x0_in)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(x_space_->dim(), x0_in.size());
#endif
  Thyra::DetachedVectorView<Scalar> x0(x0_);
  x0.sv().values()().assign(x0_in);
}


template<class Scalar>
void ModelEvaluatorHeq<Scalar>::setShowGetInvalidArgs(bool showGetInvalidArg)
{
  showGetInvalidArg_ = showGetInvalidArg;
}

template<class Scalar>
void ModelEvaluatorHeq<Scalar>::
set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> >& W_factory)
{
  W_factory_ = W_factory;
}

// Public functions overridden from ModelEvaulator


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
ModelEvaluatorHeq<Scalar>::get_x_space() const
{
  return x_space_;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
ModelEvaluatorHeq<Scalar>::get_f_space() const
{
  return f_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluatorHeq<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
ModelEvaluatorHeq<Scalar>::create_W_op() const
{
  Teuchos::RCP<Epetra_CrsMatrix> W_epetra =
    Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*x_map_,num_global_elements_));

  return Thyra::nonconstEpetraLinearOp(W_epetra);
}

template<class Scalar>
Teuchos::RCP< ::Thyra::PreconditionerBase<Scalar> >
ModelEvaluatorHeq<Scalar>::create_W_prec() const
{
  Teuchos::RCP<Epetra_CrsMatrix> W_epetra =
    Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*x_map_,1));

  const Teuchos::RCP<Thyra::LinearOpBase< Scalar > > W_op =
    Thyra::nonconstEpetraLinearOp(W_epetra);

  Teuchos::RCP<Thyra::DefaultPreconditioner<Scalar> > prec =
    Teuchos::rcp(new Thyra::DefaultPreconditioner<Scalar>);

  prec->initializeRight(W_op);

  return prec;
}

template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
ModelEvaluatorHeq<Scalar>::get_W_factory() const
{
  return W_factory_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluatorHeq<Scalar>::createInArgs() const
{
  return prototypeInArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
ModelEvaluatorHeq<Scalar>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}


template<class Scalar>
void ModelEvaluatorHeq<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  TEUCHOS_ASSERT(nonnull(inArgs.get_x()));

  const RCP<Thyra::VectorBase<Scalar> > f_out = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();
  const RCP<Thyra::PreconditionerBase<Scalar> > W_prec_out = outArgs.get_W_prec();

  if (is_null(x_ptr_)) {
    x_ptr_ = Teuchos::rcp(new Epetra_Vector(*x_map_));
  }
  x_ptr_ = Thyra::get_Epetra_Vector(*x_map_,inArgs.get_x());
  const Epetra_Vector& x = *x_ptr_;

  if ( nonnull(f_out) || nonnull(W_out) || nonnull(W_prec_out) ) {

    // ****************
    // Get the underlying epetra objects
    // ****************

    RCP<Epetra_Vector> f;
    if (nonnull(f_out)) {
      f = Thyra::get_Epetra_Vector(*f_map_,outArgs.get_f());
    }

    RCP<Epetra_CrsMatrix> J;
    if (nonnull(W_out)) {
      RCP<Epetra_Operator> W_epetra = Thyra::get_Epetra_Operator(*W_out);
      J = rcp_dynamic_cast<Epetra_CrsMatrix>(W_epetra);
      TEUCHOS_ASSERT(nonnull(J));
    }

    RCP<Epetra_CrsMatrix> M_inv;
    if (nonnull(W_prec_out)) {
      RCP<Epetra_Operator> M_epetra = Thyra::get_Epetra_Operator(*(W_prec_out->getNonconstRightPrecOp()));
      M_inv = rcp_dynamic_cast<Epetra_CrsMatrix>(M_epetra);
      TEUCHOS_ASSERT(nonnull(M_inv));
    }

    // Zero out the objects that will be filled
    if (nonnull(f))
      f->PutScalar(0.0);
    if (nonnull(J))
      J->PutScalar(0.0);
    if (nonnull(M_inv))
      M_inv->PutScalar(0.0);

    // Shared computation for f or J evaluation
    RCP<Epetra_Vector> workVec = rcp(new Epetra_Vector(*x_map_));
    if (nonnull(f) || nonnull(J)) {
      A_heq_->Multiply(false, x, *workVec);
      workVec->Update(1.0,*ones_,-1.0);
      workVec->Reciprocal(*workVec);
    }

    // Computing f
    if (nonnull(f)){
      *f = *workVec;
      f->Update(-1.0,x,1.0);
    }

    // Computing Jacobian
    if (nonnull(J)){
      *J = *A_heq_;
      workVec->Multiply(1.0,*workVec,*workVec,0.0);
      J->LeftScale(*workVec);
      // Subtract off identity
      int num_local_elements = x_map_->NumMyElements();
      for (int i=0; i<num_local_elements; i++) {
        int row = x_map_->GID(i);
        Scalar value = -1.0;
        J->SumIntoGlobalValues(row,1,&value,&row);
      }
      J->FillComplete();
    }

    // Preconditioner, set to identity
    if (nonnull(M_inv)){
      M_inv->ReplaceDiagonalValues(*ones_);
      M_inv->FillComplete();
    }
  }
}


#endif
