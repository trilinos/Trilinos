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
#include "Epetra_Import.h"
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
  x_owned_map_ = rcp(new Epetra_Map(num_global_elements_,0,*comm_));
  x_space_ = ::Thyra::create_VectorSpace(x_owned_map_);

  // residual space
  f_owned_map_ = x_owned_map_;
  f_space_ = x_space_;

  x0_ = ::Thyra::createMember(x_space_);
  V_S(x0_.ptr(), ST::zero());

  // Initialize the graph for W CrsMatrix object
  W_graph_ = createGraph();

  // Create the nodal coordinates
  {
    node_coordinates_ = Teuchos::rcp(new Epetra_Vector(*x_owned_map_));
    Scalar dx = 1.0/((double) num_global_elements_);
    for (int i=0; i < x_owned_map_->NumMyElements(); i++) {
      (*node_coordinates_)[i] = 0.5*dx + dx*((double) x_owned_map_->MinMyGID() + i);
    }
  }

  // Compute the kernel of the integral operator, A_heq
  A_heq = rcp(new Epetra_MultiVector(*x_owned_map_,num_global_elements));
  for (int i=0; i<num_global_elements_; i++)
    (*A_heq)(i)->PutScalar((*node_coordinates_)[i]);
  Teuchos::RCP<Epetra_MultiVector> Temp = rcp(new Epetra_MultiVector(*A_heq));
  Teuchos::RCP<Epetra_MultiVector> Eye = rcp(new Epetra_MultiVector(*x_owned_map_,num_global_elements_));
  for (int i=0; i<num_global_elements_; i++)
    Eye->ReplaceGlobalValue(i,i,1.0);
  A_heq->Multiply('N','T',1.0,*Eye,*Temp,0.0); //Is there a better way to transpose A_heq?
  Temp->Multiply('N','N',1.0,*Eye,*A_heq,1.0);
  Temp->Reciprocal(*Temp);
  double cc = 0.5*paramC/((double) num_global_elements);
  A_heq->Multiply(cc,*A_heq,*Temp,0.0);

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
Teuchos::RCP<Epetra_CrsGraph>
ModelEvaluatorHeq<Scalar>::createGraph()
{
  Teuchos::RCP<Epetra_CrsGraph> W_graph;

  // Create the shell for the
  W_graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *x_owned_map_, 5));

  W_graph->FillComplete();
  return W_graph;
}

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
    Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));

  return Thyra::nonconstEpetraLinearOp(W_epetra);
}

template<class Scalar>
Teuchos::RCP< ::Thyra::PreconditionerBase<Scalar> >
ModelEvaluatorHeq<Scalar>::create_W_prec() const
{
  Teuchos::RCP<Epetra_CrsMatrix> W_epetra =
    Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));

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

//  RCP<const Epetra_Vector> x;
//  x = Thyra::get_Epetra_Vector(*x_owned_map_,inArgs.get_x());
  if (is_null(x_ptr)) {
    x_ptr = Teuchos::rcp(new Epetra_Vector(*x_owned_map_));
  }

  x_ptr = Thyra::get_Epetra_Vector(*x_owned_map_,inArgs.get_x());
  const Epetra_Vector& x = *x_ptr;

  if ( nonnull(f_out) || nonnull(W_out) || nonnull(W_prec_out) ) {

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

    RCP<Epetra_CrsMatrix> M_inv;
    if (nonnull(W_prec_out)) {
      RCP<Epetra_Operator> M_epetra = Thyra::get_Epetra_Operator(*(W_prec_out->getNonconstRightPrecOp()));
      M_inv = rcp_dynamic_cast<Epetra_CrsMatrix>(M_epetra);
      TEUCHOS_ASSERT(nonnull(M_inv));
      J_diagonal_ = Teuchos::rcp(new Epetra_Vector(*x_owned_map_));
    }

    int ierr = 0;

    // Zero out the objects that will be filled
    if (nonnull(f))
      f->PutScalar(1.0);
    if (nonnull(J))
      J->PutScalar(0.0);
    if (nonnull(M_inv))
      M_inv->PutScalar(0.0);

    // Computing f
    RCP<Epetra_Vector> Scale = rcp(new Epetra_Vector(*x_owned_map_));
    if (nonnull(f)){
//      x.Print(std::cout);
      ierr = f->Multiply('N','N',-1.0,*A_heq,x,1.0);
      ierr = f->Reciprocal(*f);
      *Scale = *f;
      ierr = f->Update(-1.0,x,1.0);
    }

    // Computing Jacobian
    if (nonnull(J)){
      Scale->Multiply(-1.0,*Scale,*Scale,0.0);
      for (int i=0; i<num_global_elements_; i++){
        for (int j=0; j<num_global_elements_; j++){
          double Value = (*A_heq)[i][j];
          J->ReplaceGlobalValues(i,1,&Value,&j);
        }
      }
      J->FillComplete();
      J->LeftScale(*Scale);
      double Value = 1.0;
      for (int i=0; i<num_global_elements_; i++)
        J->SumIntoGlobalValues(i,1,&Value,&i);
    }

    // Preconditioner
    if (nonnull(M_inv)){
      RCP<Epetra_Vector> Diag = Teuchos::rcp(new Epetra_Vector(*x_owned_map_));
      Diag->PutScalar(1.0);
      M_inv->ReplaceDiagonalValues(*Diag);
      M_inv->FillComplete();
    }

    TEUCHOS_ASSERT(ierr > -1);

  }

}


#endif
