

#ifndef NOX_THYRA_MODEL_EVALUATOR_2DSIM_DEF_HPP
#define NOX_THYRA_MODEL_EVALUATOR_2DSIM_DEF_HPP


//#include "Thyra_SimpleDenseLinearOp.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"

// Epetra support
#include "Thyra_EpetraThyraWrappers.cpp"
#include "Thyra_get_Epetra_Operator.hpp"

// Nonmember constuctors

template<class Scalar>
Teuchos::RCP<ModelEvaluator2DSim<Scalar> >
modelEvaluator2DSim(const Teuchos::RCP<const Epetra_Comm>& comm,
		    const Scalar d,
		    const Scalar p0,
		    const Scalar p1,
		    const Scalar x0,
		    const Scalar x1)
{
  return Teuchos::rcp(new ModelEvaluator2DSim<Scalar>(comm,d,p0,p1,x0,x1));
}

// Initializers/Accessors

template<class Scalar>
void ModelEvaluator2DSim<Scalar>::set_d(const Scalar &d)
{
  d_ = d;
}


template<class Scalar>
void ModelEvaluator2DSim<Scalar>::set_p(const Teuchos::ArrayView<const Scalar> &p)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(p_.size(), p.size());
#endif
  p_().assign(p);
}


template<class Scalar>
void ModelEvaluator2DSim<Scalar>::set_x0(const Teuchos::ArrayView<const Scalar> &x0_in)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(x_space_->dim(), x0_in.size());
#endif
  Thyra::DetachedVectorView<Scalar> x0(x0_);
  x0.sv().values()().assign(x0_in);
}


template<class Scalar>
void ModelEvaluator2DSim<Scalar>::setShowGetInvalidArgs(bool showGetInvalidArg)
{
  showGetInvalidArg_ = showGetInvalidArg;
}

template<class Scalar>
void ModelEvaluator2DSim<Scalar>::
set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> >& W_factory)
{
  W_factory_ = W_factory;
}

// Public functions overridden from ModelEvaulator


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
ModelEvaluator2DSim<Scalar>::get_x_space() const
{
  return x_space_;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
ModelEvaluator2DSim<Scalar>::get_f_space() const
{
  return f_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluator2DSim<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
ModelEvaluator2DSim<Scalar>::create_W_op() const
{
  Teuchos::RCP<Epetra_CrsMatrix> W_epetra =
    Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));

  return Thyra::nonconstEpetraLinearOp(W_epetra);
}


template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
ModelEvaluator2DSim<Scalar>::get_W_factory() const
{
  return W_factory_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluator2DSim<Scalar>::createInArgs() const
{
  return prototypeInArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
ModelEvaluator2DSim<Scalar>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}


template<class Scalar>
void ModelEvaluator2DSim<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  
  TEUCHOS_ASSERT(nonnull(inArgs.get_x()));
		 
  const Thyra::ConstDetachedVectorView<Scalar> x(inArgs.get_x());

  const RCP< Thyra::VectorBase<Scalar> > f_out = outArgs.get_f();
  const RCP< Thyra::LinearOpBase< Scalar > > W_out = outArgs.get_W_op();

  if (nonnull(f_out)) {
    NOX_FUNC_TIME_MONITOR("ModelEvaluator2DSim::eval f_out");
    const Thyra::DetachedVectorView<Scalar> f(f_out);
    f[0] = x[0] + x[1] * x[1] - p_[0];
    f[1] = d_ * (x[0] * x[0] - x[1] - p_[1]);
  }

  if (nonnull(W_out)) {
    NOX_FUNC_TIME_MONITOR("ModelEvaluator2DSim::eval W_op_out");
    RCP<Epetra_Operator> W_epetra= Thyra::get_Epetra_Operator(*W_out);
    RCP<Epetra_CrsMatrix> W_epetracrs = rcp_dynamic_cast<Epetra_CrsMatrix>(W_epetra);
    TEUCHOS_ASSERT(nonnull(W_epetracrs));
    Epetra_CrsMatrix& DfDx = *W_epetracrs;
    DfDx.PutScalar(0.0);
    //
    // Fill W = DfDx
    //
    // W = DfDx = [      1.0,  2*x[1] ]
    //            [ 2*d*x[0],     -d  ]
    //
    double values[2];
    int indexes[2];
    // Row [0]
    values[0] = 1.0;           indexes[0] = 0;
    values[1] = 2.0*x[1];      indexes[1] = 1;
    DfDx.SumIntoGlobalValues( 0, 2, values, indexes );
    // Row [1]
    values[0] = 2.0*d_*x[0];   indexes[0] = 0;
    values[1] = -d_;           indexes[1] = 1;
    DfDx.SumIntoGlobalValues( 1, 2, values, indexes );
  }

}


// private


template<class Scalar>
ModelEvaluator2DSim<Scalar>::ModelEvaluator2DSim(const Teuchos::RCP<const Epetra_Comm>& comm,
						 const Scalar d,
						 const Scalar p0,
						 const Scalar p1,
						 const Scalar x0,
						 const Scalar x1) :
  epetra_comm_(comm),
  d_(d),
  p_(2, Teuchos::ScalarTraits<Scalar>::zero()),
  showGetInvalidArg_(false)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using ::Thyra::VectorBase;
  typedef ::Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  
  TEUCHOS_ASSERT(nonnull(epetra_comm_));

  const int nx = 2;

  x_epetra_map_ = rcp(new Epetra_Map(nx,0,*epetra_comm_));
  x_space_ = ::Thyra::create_VectorSpace(x_epetra_map_);

  f_epetra_map_ = x_epetra_map_;
  f_space_ = x_space_;

  x0_ = ::Thyra::createMember(x_space_);
  V_S(x0_.ptr(), ST::zero());

  set_p(Teuchos::tuple<Scalar>(p0, p1)());
  set_x0(Teuchos::tuple<Scalar>(x0, x1)());

  // Initialize the graph for W CrsMatrix object
  W_graph_ = rcp(new Epetra_CrsGraph(::Copy,*x_epetra_map_,nx));
  {
    int indices[nx] = { 0, 1 };
    for( int i = 0; i < nx; ++i )
      W_graph_->InsertGlobalIndices(i,nx,indices);
  }
  W_graph_->FillComplete();

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  prototypeInArgs_ = inArgs;
  
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
//   outArgs.set_W_properties(DerivativeProperties(
// 			     DERIV_LINEARITY_NONCONST
// 			     ,DERIV_RANK_FULL
// 			     ,true // supportsAdjoint
// 			     ));
  prototypeOutArgs_ = outArgs;

  nominalValues_ = inArgs;
  nominalValues_.set_x(x0_);
}

#endif
