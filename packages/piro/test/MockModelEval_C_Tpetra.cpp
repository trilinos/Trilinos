// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MockModelEval_C_Tpetra.hpp"
#include "Thyra_LinearOpWithSolveBase_decl.hpp"


#include <iostream>

using Teuchos::RCP;
using Teuchos::rcp;

MockModelEval_C_Tpetra::MockModelEval_C_Tpetra(const Teuchos::RCP<const Teuchos::Comm<int> >  appComm, bool /*adjoint*/, const Teuchos::RCP<Teuchos::ParameterList>& problemList, bool hessianSupport) //problem is self-adjoint
 {
    comm = appComm;
    hessSupport = hessianSupport;

    //set up map and initial guess for solution vector
    const int vecLength = 1;
    x_map = rcp(new Tpetra_Map(vecLength, 0, comm));
    x_vec = rcp(new Tpetra_Vector(x_map));
    x_dot_vec = rcp(new Tpetra_Vector(x_map));
    x_dotdot_vec = rcp(new Tpetra_Vector(x_map));
    x_vec->putScalar(0.0);
    x_dot_vec->putScalar(1.0);
    x_dotdot_vec->putScalar(0.0);

    Teuchos::RCP<const Thyra::VectorSpaceBase<double>> x_space =
        Thyra::createVectorSpace<double>(x_map);

    //set up responses
    const int numResponses = 1;
    g_map = rcp(new const Tpetra_Map(numResponses , 0, comm, Tpetra::LocallyReplicated));

    //set up parameters
    const int numParameters= 1;
    p_map = rcp(new const Tpetra_Map(numParameters, 0, comm, Tpetra::LocallyReplicated));

    Teuchos::RCP<const Thyra::VectorSpaceBase<double>> p_space =
        Thyra::createVectorSpace<double>(p_map);


    Teuchos::RCP<Tpetra_Vector> p_init = rcp(new Tpetra_Vector(p_map));
    Teuchos::RCP<Tpetra_Vector> p_lo = rcp(new Tpetra_Vector(p_map));
    Teuchos::RCP<Tpetra_Vector> p_up = rcp(new Tpetra_Vector(p_map));
    p_init->putScalar(1.0);
    p_lo->putScalar(0.1);
    p_up->putScalar(10.0);

    p_vec = rcp(new Tpetra_Vector(p_map));
    p_vec->assign(*p_init);

    //set up jacobian graph
    crs_graph = rcp(new Tpetra_CrsGraph(x_map, vecLength));
    {
      std::vector<typename Tpetra_CrsGraph::global_ordinal_type> indices(vecLength);
      for (int i=0; i<vecLength; i++) indices[i]=i;
      const int nodeNumElements = x_map->getLocalNumElements();
      for (int i=0; i<nodeNumElements; i++)
        crs_graph->insertGlobalIndices(x_map->getGlobalElement(i), vecLength, &indices[0]);
    }
    crs_graph->fillComplete();

    //set up hessian graph
    hess_crs_graph = rcp(new Tpetra_CrsGraph(p_map, numParameters));
    if (comm->getRank() == 0)
    {
      std::vector<typename Tpetra_CrsGraph::global_ordinal_type> indices(numParameters);
      for (int i=0; i<numParameters; i++) indices[i]=i;
      const int nodeNumElements = p_map->getLocalNumElements();
      for (int i=0; i<nodeNumElements; i++)
        hess_crs_graph->insertGlobalIndices(p_map->getGlobalElement(i), numParameters, &indices[0]);
    }
    hess_crs_graph->fillComplete();

    // Setup nominal values, lower and upper bounds
    nominalValues = this->createInArgsImpl();
    lowerBounds = this->createInArgsImpl();
    upperBounds = this->createInArgsImpl();

    nominalValues.set_x(Thyra::createVector(x_vec, x_space));
    nominalValues.set_x_dot(Thyra::createVector(x_dot_vec, x_space));
    this->xDotDot = Thyra::createVector(x_dotdot_vec, x_space);
    nominalValues.set_x_dot_dot(this->xDotDot);
      
 
    nominalValues.set_p(0, Thyra::createVector(p_init, p_space));
    lowerBounds.set_p(0, Thyra::createVector(p_lo, p_space));
    upperBounds.set_p(0, Thyra::createVector(p_up, p_space));

    probList_ = problemList;
}

MockModelEval_C_Tpetra::~MockModelEval_C_Tpetra()
{
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_C_Tpetra::get_x_space() const
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> x_space =
      Thyra::createVectorSpace<double>(x_map);
  return x_space;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_C_Tpetra::get_f_space() const
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> f_space =
      Thyra::createVectorSpace<double>(x_map);
  return f_space;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_C_Tpetra::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> p_space =
        Thyra::createVectorSpace<double>(p_map);
  return p_space;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_C_Tpetra::get_g_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  MockModelEval_C_Tpetra::get_g_map() only " <<
                     " supports 1 response.  Supplied index l = " <<
                     l << std::endl);
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> g_space =
        Thyra::createVectorSpace<double>(g_map);
  return g_space;
}

RCP<const  Teuchos::Array<std::string> > MockModelEval_C_Tpetra::get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_names() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);

  Teuchos::Ordinal num_p = p_map->getLocalNumElements();
  RCP<Teuchos::Array<std::string> > p_names =
      rcp(new Teuchos::Array<std::string>(num_p) );
  for (int i=0; i<num_p; i++) {
    std::stringstream ss;
    ss << "Parameter " << i;
    const std::string name = ss.str();
    (*p_names)[i] = name;
  }
  return p_names;
}


Teuchos::RCP<Thyra::LinearOpBase<double>>
MockModelEval_C_Tpetra::create_W_op() const
{
  const Teuchos::RCP<Tpetra_Operator> W =
      Teuchos::rcp(new Tpetra_CrsMatrix(crs_graph));
  return Thyra::createLinearOp(W);
}

//! Create preconditioner operator
Teuchos::RCP<Thyra::PreconditionerBase<double>>
MockModelEval_C_Tpetra::create_W_prec() const
{
  return Teuchos::null;
}

Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double>>
MockModelEval_C_Tpetra::get_W_factory() const
{
  return Teuchos::null;
}

Teuchos::RCP<Thyra::LinearOpBase<double>>
MockModelEval_C_Tpetra::create_hess_g_pp( int j, int l1, int l2 ) const
{
  const Teuchos::RCP<Tpetra_Operator> H =
      Teuchos::rcp(new Tpetra_CrsMatrix(hess_crs_graph));
  return Teuchos::rcp(new MatrixBased_LOWS(Thyra::createLinearOp(H)));
}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_C_Tpetra::getNominalValues() const
{
  return nominalValues;
}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_C_Tpetra::getLowerBounds() const
{
  return lowerBounds;
}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_C_Tpetra::getUpperBounds() const
{
  return upperBounds;
}


Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_C_Tpetra::createInArgs() const
{
  return this->createInArgsImpl();
}

void
MockModelEval_C_Tpetra::reportFinalPoint(
    const Thyra::ModelEvaluatorBase::InArgs<double>& /* finalPoint */,
    const bool /* wasSolved */) {
  // Do nothing  
}

Thyra::ModelEvaluatorBase::OutArgs<double>
MockModelEval_C_Tpetra::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<double> result;
  result.setModelEvalDescription(this->description());
  result.set_Np_Ng(1, 1);

  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f, true);
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op, true);
  result.set_W_properties(Thyra::ModelEvaluatorBase::DerivativeProperties(
      Thyra::ModelEvaluatorBase::DERIV_LINEARITY_UNKNOWN,
      Thyra::ModelEvaluatorBase::DERIV_RANK_FULL,
      true));

  result.setSupports(
      Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, 0, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
  result.setSupports(
      Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, 0, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, 0, 0, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);

  result.setSupports(
        Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, 0, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);

  result.setHessianSupports(hessSupport);

  return result;
}

void MockModelEval_C_Tpetra::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<double>&  inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<double>& outArgs) const
{

  // Parse InArgs

  const Teuchos::RCP<const Tpetra_Vector> x_in =
      ConverterT::getConstTpetraVector(inArgs.get_x());
  if (!Teuchos::nonnull(x_in)) std::cerr << "ERROR: MockModelEval_C_Tpetra requires x as inargs\n";

  const Teuchos::RCP<const Tpetra_Vector> x_dot_in =
      Teuchos::nonnull(inArgs.get_x_dot()) ?
          ConverterT::getConstTpetraVector(inArgs.get_x_dot()) :
          Teuchos::null;

  const Teuchos::RCP<const Thyra::VectorBase<double>> p_in = inArgs.get_p(0);
  if (Teuchos::nonnull(p_in))
    p_vec->assign(*ConverterT::getConstTpetraVector(p_in));

  int myVecLength = x_in->getLocalLength();

  // Parse OutArgs

  const Teuchos::RCP<Tpetra_Vector> f_out =
      Teuchos::nonnull(outArgs.get_f()) ?
          ConverterT::getTpetraVector(outArgs.get_f()) :
          Teuchos::null;

  const Teuchos::RCP<Thyra::VectorBase<double>> g_base = outArgs.get_g(0);
  Teuchos::RCP<Tpetra_Vector> g_out = Teuchos::nonnull(g_base) ?
      ConverterT::getTpetraVector(g_base) :
      Teuchos::null;

  const Teuchos::RCP<Tpetra_Operator> W_out =
      Teuchos::nonnull(outArgs.get_W_op()) ?
          ConverterT::getTpetraOperator(outArgs.get_W_op()) :
          Teuchos::null;

  const Teuchos::RCP<Thyra::MultiVectorBase<double>> dfdp_base =
      outArgs.get_DfDp(0).getMultiVector();

  const Teuchos::RCP<Tpetra_MultiVector> dfdp_out =
      Teuchos::nonnull(dfdp_base) ?
          ConverterT::getTpetraMultiVector(dfdp_base) :
          Teuchos::null;

  const Teuchos::RCP<Thyra::MultiVectorBase<double>> dgdp_base =
      outArgs.get_DgDp(0, 0).getMultiVector();
  const Teuchos::RCP<Tpetra_MultiVector> dgdp_out =
      Teuchos::nonnull(dgdp_base) ?
          ConverterT::getTpetraMultiVector(dgdp_base) :
          Teuchos::null;

  const Teuchos::RCP<Thyra::MultiVectorBase<double>> dgdx_base =
        outArgs.get_DgDx(0).getMultiVector();
  const Teuchos::RCP<Tpetra_MultiVector> dgdx_out =
      Teuchos::nonnull(dgdx_base) ?
          ConverterT::getTpetraMultiVector(dgdx_base) :
          Teuchos::null;

  const Teuchos::RCP<const Tpetra_MultiVector> p_direction =
      Teuchos::nonnull(inArgs.get_p_direction(0)) ?
        ConverterT::getConstTpetraMultiVector(inArgs.get_p_direction(0)):
        Teuchos::null;


  const Teuchos::RCP<const Tpetra_Vector> lag_multiplier_f_in =
      Teuchos::nonnull(inArgs.get_f_multiplier()) ?
        ConverterT::getConstTpetraVector(inArgs.get_f_multiplier()) :
        Teuchos::null;

  auto f_hess_xx_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_f_xx) ? outArgs.get_hess_vec_prod_f_xx() : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_xx_v_out =
      Teuchos::nonnull(f_hess_xx_v) ?
        ConverterT::getTpetraMultiVector(f_hess_xx_v) :
        Teuchos::null;

  auto f_hess_xp_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_f_xp,0) ? outArgs.get_hess_vec_prod_f_xp(0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_xp_v_out =
      Teuchos::nonnull(f_hess_xp_v) ?
        ConverterT::getTpetraMultiVector(f_hess_xp_v) :
        Teuchos::null;

  auto f_hess_px_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_f_px,0) ? outArgs.get_hess_vec_prod_f_px(0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_px_v_out =
      Teuchos::nonnull(f_hess_px_v) ?
        ConverterT::getTpetraMultiVector(f_hess_px_v) :
        Teuchos::null;

  auto f_hess_pp_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_f_pp,0,0) ? outArgs.get_hess_vec_prod_f_pp(0,0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_pp_v_out =
      Teuchos::nonnull(f_hess_pp_v) ?
        ConverterT::getTpetraMultiVector(f_hess_pp_v) :
        Teuchos::null;

  const Teuchos::RCP<const Tpetra_MultiVector> x_direction =
      Teuchos::nonnull(inArgs.get_x_direction()) ?
        ConverterT::getConstTpetraMultiVector(inArgs.get_x_direction()):
        Teuchos::null;


  const Teuchos::RCP<const Tpetra_Vector> lag_multiplier_g_in =
      Teuchos::nonnull(inArgs.get_g_multiplier(0)) ?
        ConverterT::getConstTpetraVector(inArgs.get_g_multiplier(0)) :
        Teuchos::null;

  auto g_hess_xx_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_g_xx,0) ? outArgs.get_hess_vec_prod_g_xx(0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_xx_v_out =
      Teuchos::nonnull(g_hess_xx_v) ?
        ConverterT::getTpetraMultiVector(g_hess_xx_v) :
        Teuchos::null;

  auto g_hess_xp_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_g_xp,0,0) ? outArgs.get_hess_vec_prod_g_xp(0,0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_xp_v_out =
      Teuchos::nonnull(g_hess_xp_v) ?
        ConverterT::getTpetraMultiVector(g_hess_xp_v) :
        Teuchos::null;

  auto g_hess_px_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_g_px,0,0) ? outArgs.get_hess_vec_prod_g_px(0,0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_px_v_out =
      Teuchos::nonnull(g_hess_px_v) ?
        ConverterT::getTpetraMultiVector(g_hess_px_v) :
        Teuchos::null;

  auto g_hess_pp_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_g_pp,0,0,0) ? outArgs.get_hess_vec_prod_g_pp(0,0,0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_pp_v_out =
      Teuchos::nonnull(g_hess_pp_v) ?
        ConverterT::getTpetraMultiVector(g_hess_pp_v) :
        Teuchos::null;

  Teuchos::RCP<const Tpetra_Vector> x_dotdot_in;
  double omega2(0.0);
  // This assumes we are using Piro second order integrator 
  if (Teuchos::nonnull(this->get_x_dotdot())) {
    const Teuchos::RCP<const Thyra::VectorBase<double>> x_dotdot = this->get_x_dotdot();
    x_dotdot_in = ConverterT::getConstTpetraVector(x_dotdot);
    omega2    = this->get_omega();
  }

  // If using Tempus, we should have 
  //  x_dotdot_in = Tuchos::nonnull(inArgs.get_x_dot_dot()) ?
  //        ConverterT::getConstTpetraVector(inArgs.get_x_dot_dot()) :
  //        Teuchos::null;
  //omega2 = inArgs.get_W_x_dot_dot_coeff();

  const double damping = 0;
  if (Teuchos::nonnull(f_out)) {
    f_out->putScalar(0.0);
    auto p = p_vec->getLocalViewHost(Tpetra::Access::ReadOnly);
    {
      auto f_out_view = f_out->getLocalViewHost(Tpetra::Access::OverwriteAll);
      for (int i=0; i<myVecLength; i++)
        f_out_view(i,0) = -p(0,0);
    }

    if (Teuchos::nonnull(x_dotdot_in)) {
      // f(x, x_dot) = x_dotdot - f(x)
      auto f_out_view = f_out->getLocalViewHost(Tpetra::Access::ReadWrite);
      auto x_dotdot_in_view = x_dotdot_in->getLocalViewHost(Tpetra::Access::ReadOnly);
      for (int i=0; i<myVecLength; i++) {
        f_out_view(i,0) = x_dotdot_in_view(i,0) - f_out_view(i,0);
      }
    }

    if (Teuchos::nonnull(f_out) && Teuchos::nonnull(x_dot_in)) {
      // f(x, x_dot) = f(x) -  damping * x_dot
      auto f_out_view = f_out->getLocalViewHost(Tpetra::Access::ReadWrite);
      auto x_dot_in_view = x_dot_in->getLocalViewHost(Tpetra::Access::ReadOnly);
      for (int i=0; i<myVecLength; i++) {
        f_out_view(i,0) += damping * x_dot_in_view(i,0);
      }
    }
  }

  if (W_out != Teuchos::null) {
    Teuchos::RCP<Tpetra_CrsMatrix> W_out_crs =
      Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(W_out, true);
    W_out_crs->resumeFill();
    W_out_crs->setAllToScalar(0.0);

    double diag=0.0;
    for (int i=0; i<myVecLength; i++)
      W_out_crs->replaceLocalValues(i, 1, &diag, &i);

    double alpha = inArgs.get_alpha();
    double beta = inArgs.get_beta();

    // W(x, x_dot, x_dotdot) = omega2 * Id + alpha * damping * Id - beta * W(x)
    if (Teuchos::nonnull(x_dotdot_in)) {
      W_out_crs->resumeFill();
      W_out_crs->scale(-beta);
      diag = omega2;
      for (int i=0; i<myVecLength; i++) {
        W_out_crs->sumIntoLocalValues(i, 1, &diag, &i);
      }
    }

    if (Teuchos::nonnull(x_dot_in)) {
      W_out_crs->resumeFill();
      diag = alpha*damping;
      for (int i=0; i<myVecLength; i++) {
        W_out_crs->sumIntoLocalValues(i, 1, &diag, &i);
      }
    }
    
    W_out_crs->fillComplete();
  }

  if (Teuchos::nonnull(dfdp_out)) {
    if (Teuchos::nonnull(x_dotdot_in))
      dfdp_out->putScalar(1.0);
    else
      dfdp_out->putScalar(-1.0);
  }

  if (Teuchos::nonnull(g_out)) {
    double mean = x_in->meanValue();
    g_out->putScalar(mean);
  }

  if (dgdx_out != Teuchos::null) {
    dgdx_out->putScalar(1.0/x_in->getGlobalLength());
  }

  if (Teuchos::nonnull(f_hess_xx_v_out)) {
    f_hess_xx_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_xp_v_out)) {
    f_hess_xp_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_px_v_out)) {
    f_hess_px_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_pp_v_out)) {
    f_hess_pp_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_xx_v_out)) {
    g_hess_xx_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_xp_v_out)) {
    g_hess_xp_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_px_v_out)) {
    g_hess_px_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_pp_v_out)) {
    g_hess_pp_v_out->putScalar(0);
  }

  // Modify for time dependent (implicit time integration or eigensolves)

}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_C_Tpetra::createInArgsImpl() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<double> result;
  result.setModelEvalDescription(this->description());

  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot_dot, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_t, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_alpha, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_beta, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_W_x_dot_dot_coeff,true);

  result.set_Np_Ng(1,1);

  return result;
}

