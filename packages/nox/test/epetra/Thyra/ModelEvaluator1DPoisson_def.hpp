// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef NOX_THYRA_MODEL_EVALUATOR_1DFEM_DEF_HPP
#define NOX_THYRA_MODEL_EVALUATOR_1DFEM_DEF_HPP

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
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

// Nonmember constuctors

template<class Scalar>
Teuchos::RCP<ModelEvaluator1DPoisson<Scalar> >
modelEvaluator1DPoisson(const Teuchos::RCP<const Epetra_Comm>& comm,
            const int num_global_elements,
            const Scalar a,
            const Scalar b)
{
  return Teuchos::rcp(new ModelEvaluator1DPoisson<Scalar>(comm,num_global_elements,a,b));
}

// Constructor

template<class Scalar>
ModelEvaluator1DPoisson<Scalar>::
ModelEvaluator1DPoisson(const Teuchos::RCP<const Epetra_Comm>& comm,
            const int num_global_elements,
            const Scalar a,
            const Scalar b) :
  comm_(comm),
  num_global_elements_(num_global_elements),
  a_(a),
  b_(b),
  showGetInvalidArg_(false)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using ::Thyra::VectorBase;
  typedef ::Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEUCHOS_ASSERT(nonnull(comm_));

  const int num_global_nodes = num_global_elements_ + 1;

  // owned space
  x_owned_map_ = rcp(new Epetra_Map(num_global_nodes,0,*comm_));
  x_space_ = ::Thyra::create_VectorSpace(x_owned_map_);

  // ghosted space
  if (comm_->NumProc() == 1) {
    x_ghosted_map_ = x_owned_map_;
  } else {

    int num_overlap_nodes;
    int min_overlap_gid;
    num_overlap_nodes = x_owned_map_->NumMyElements() + 2;
    if ( (comm_->MyPID() == 0) || (comm_->MyPID() == (comm_->NumProc() - 1)) )
      num_overlap_nodes --;

    if (comm_->MyPID() == 0)
      min_overlap_gid = x_owned_map_->MinMyGID();
    else
      min_overlap_gid = x_owned_map_->MinMyGID() - 1;

    int* overlap_gids = new int[num_overlap_nodes];

    for (int i = 0; i < num_overlap_nodes; i ++)
      overlap_gids[i] = min_overlap_gid + i;

    x_ghosted_map_ =
      Teuchos::rcp(new Epetra_Map(-1,
                  num_overlap_nodes,
                  overlap_gids,
                  0,
                  *comm_));

    delete [] overlap_gids;

  }

  importer_ = Teuchos::rcp(new Epetra_Import(*x_ghosted_map_, *x_owned_map_));

  // residual space
  f_owned_map_ = x_owned_map_;
  f_space_ = x_space_;

  x0_ = ::Thyra::createMember(x_space_);
  V_S(x0_.ptr(), ST::zero());

  // Initialize the graph for W CrsMatrix object
  W_graph_ = createGraph();

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  prototypeInArgs_ = inArgs;

  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  prototypeOutArgs_ = outArgs;

  nominalValues_ = inArgs;
  nominalValues_.set_x(x0_);
}

// Initializers/Accessors

template<class Scalar>
Teuchos::RCP<Epetra_CrsGraph>
ModelEvaluator1DPoisson<Scalar>::createGraph()
{
  Teuchos::RCP<Epetra_CrsGraph> W_graph;

  // Create the shell for the
  W_graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *x_owned_map_, 3));

  // Declare required variables
  int num_owned_nodes = x_owned_map_->NumMyElements();

  // Loop Over # of nodes on Processor
  for (int i=0; i < num_owned_nodes; i++) {
    int gid = x_owned_map_->GID(i);
    if (gid == 0) {
      int self[1] = {gid};
      W_graph->InsertGlobalIndices(gid, 1, self);
    } else if (gid == num_global_elements_) {
      int left[2] = {gid - 1, gid};
      W_graph->InsertGlobalIndices(gid, 2, left);
    } else {
      int all[3] = {gid - 1, gid, gid + 1};
      W_graph->InsertGlobalIndices(gid, 3, all);
    }
  }
  W_graph->FillComplete();
  return W_graph;
}

template<class Scalar>
void ModelEvaluator1DPoisson<Scalar>::set_x0(const Teuchos::ArrayView<const Scalar> &x0_in)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(x_space_->dim(), x0_in.size());
#endif
  Thyra::DetachedVectorView<Scalar> x0(x0_);
  x0.sv().values()().assign(x0_in);
}


template<class Scalar>
void ModelEvaluator1DPoisson<Scalar>::setShowGetInvalidArgs(bool showGetInvalidArg)
{
  showGetInvalidArg_ = showGetInvalidArg;
}

template<class Scalar>
void ModelEvaluator1DPoisson<Scalar>::
set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> >& W_factory)
{
  W_factory_ = W_factory;
}

// Public functions overridden from ModelEvaulator


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
ModelEvaluator1DPoisson<Scalar>::get_x_space() const
{
  return x_space_;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
ModelEvaluator1DPoisson<Scalar>::get_f_space() const
{
  return f_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluator1DPoisson<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
ModelEvaluator1DPoisson<Scalar>::create_W_op() const
{
  Teuchos::RCP<Epetra_CrsMatrix> W_epetra =
    Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));

  return Thyra::nonconstEpetraLinearOp(W_epetra);
}

template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
ModelEvaluator1DPoisson<Scalar>::get_W_factory() const
{
  return W_factory_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluator1DPoisson<Scalar>::createInArgs() const
{
  return prototypeInArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
ModelEvaluator1DPoisson<Scalar>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}


template<class Scalar>
void ModelEvaluator1DPoisson<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  TEUCHOS_ASSERT(nonnull(inArgs.get_x()));

  //const Thyra::ConstDetachedVectorView<Scalar> x(inArgs.get_x());

  const RCP<Thyra::VectorBase<Scalar> > f_out = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();

  if ( nonnull(f_out) || nonnull(W_out) ) {

    // ****************
    // Get the underlying epetra objects
    // ****************

    RCP<Epetra_Vector> f;
    if (nonnull(f_out)) {
      f = Thyra::get_Epetra_Vector(*f_owned_map_,outArgs.get_f());
    }

    RCP<Epetra_CrsMatrix> J;
    if (nonnull(W_out)) {
      RCP<Epetra_Operator> W_epetra = Thyra::get_Epetra_Operator(*W_out);
      J = rcp_dynamic_cast<Epetra_CrsMatrix>(W_epetra);
      TEUCHOS_ASSERT(nonnull(J));
    }

    // ****************
    // Create ghosted objects
    // ****************

    if (is_null(u_ptr))
      u_ptr = Teuchos::rcp(new Epetra_Vector(*x_ghosted_map_));

    u_ptr->Import(*(Thyra::get_Epetra_Vector(*x_owned_map_,inArgs.get_x())), *importer_, Insert);

    Epetra_Vector& u = *u_ptr;

    // Zero out the objects that will be filled
    if (nonnull(f))
      f->PutScalar(0.0);
    if (nonnull(J))
      J->PutScalar(0.0);

    int num_owned_nodes = x_owned_map_->NumMyElements();
    for (int i=0; i < num_owned_nodes; i++) {
      int gid = x_owned_map_->GID(i);
      if (nonnull(f)) {
        if (gid == 0) {
          (*f)[i] = u[0] - a_;
        } else if (gid == num_global_elements_) {
          (*f)[i] = u[x_ghosted_map_->LID(gid)] - b_;
        } else {
          (*f)[i] = - u[x_ghosted_map_->LID(gid)]
                    + u[x_ghosted_map_->LID(gid - 1)]/2.
                    + u[x_ghosted_map_->LID(gid + 1)]/2.;
        }
      }
      if (nonnull(J)) {
        if (gid == 0) {
          double vals[1] = {1.0};
          int columns[1] = {gid};
          J->SumIntoGlobalValues(gid, 1, vals, columns);
        } else if (gid == num_global_elements_) {
          double vals[1] = {1.0};
          int columns[1] = {gid};
          J->SumIntoGlobalValues(gid, 1, vals, columns);
        } else {
          double vals[3] = {0.5, -1.0, 0.5};
          int columns[3] = {gid - 1, gid, gid + 1};
          J->SumIntoGlobalValues(gid, 3, vals, columns);
        }
      }
    }

    if (nonnull(J))
      J->FillComplete();

  }

}

#endif
