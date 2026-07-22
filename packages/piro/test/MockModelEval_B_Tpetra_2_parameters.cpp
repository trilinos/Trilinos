// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MockModelEval_B_Tpetra_2_parameters.hpp"
#include "Thyra_LinearOpWithSolveBase_decl.hpp"


#include <iostream>

using Teuchos::RCP;
using Teuchos::rcp;

MockModelEval_B_Tpetra_2_parameters::MockModelEval_B_Tpetra_2_parameters(const Teuchos::RCP<const Teuchos::Comm<int> >  appComm, bool /*adjoint*/, const Teuchos::RCP<Teuchos::ParameterList>& problemList, bool hessianSupport) //problem is self-adjoint
 {
    comm = appComm;
    hessSupport = hessianSupport;

    //set up map and initial guess for solution vector
    const int vecLength = 4;
    x_map = rcp(new Tpetra_Map(vecLength, 0, comm));
    x_vec = rcp(new Tpetra_Vector(x_map));
    x_dot_vec = rcp(new Tpetra_Vector(x_map));
    x_vec->putScalar(3.0);
    x_dot_vec->putScalar(1.0);

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

    p_init->getDataNonConst()[0]= 1.0;
    p_lo->getDataNonConst()[0]= 0.1;
    p_up->getDataNonConst()[0]= 10.0;

    p_vec_0 = rcp(new Tpetra_Vector(p_map));
    p_vec_0->assign(*p_init);

    p_vec_1 = rcp(new Tpetra_Vector(p_map));
    p_vec_1->assign(*p_init);

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

    for (int i=0; i<2; i++) {
      nominalValues.set_p(i, Thyra::createVector(p_init, p_space));
      lowerBounds.set_p(i, Thyra::createVector(p_lo, p_space));
      upperBounds.set_p(i, Thyra::createVector(p_up, p_space));
    }

    probList_ = problemList;
}

MockModelEval_B_Tpetra_2_parameters::~MockModelEval_B_Tpetra_2_parameters()
{
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_B_Tpetra_2_parameters::get_x_space() const
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> x_space =
      Thyra::createVectorSpace<double>(x_map);
  return x_space;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_B_Tpetra_2_parameters::get_f_space() const
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> f_space =
      Thyra::createVectorSpace<double>(x_map);
  return f_space;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_B_Tpetra_2_parameters::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l > 1, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 2 parameters.  Supplied index l = " <<
                     l << std::endl);
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> p_space =
        Thyra::createVectorSpace<double>(p_map);
  return p_space;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_B_Tpetra_2_parameters::get_g_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  MockModelEval_B_Tpetra_2_parameters::get_g_map() only " <<
                     " supports 1 response.  Supplied index l = " <<
                     l << std::endl);
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> g_space =
        Thyra::createVectorSpace<double>(g_map);
  return g_space;
}

RCP<const  Teuchos::Array<std::string> > MockModelEval_B_Tpetra_2_parameters::get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l > 1, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_names() only " <<
                     " supports 2 parameters.  Supplied index l = " <<
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
MockModelEval_B_Tpetra_2_parameters::create_W_op() const
{
  const Teuchos::RCP<Tpetra_Operator> W =
      Teuchos::rcp(new Tpetra_CrsMatrix(crs_graph));
  return Thyra::createLinearOp(W);
}

//! Create preconditioner operator
Teuchos::RCP<Thyra::PreconditionerBase<double>>
MockModelEval_B_Tpetra_2_parameters::create_W_prec() const
{
  return Teuchos::null;
}

Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double>>
MockModelEval_B_Tpetra_2_parameters::get_W_factory() const
{
  return Teuchos::null;
}

Teuchos::RCP<Thyra::LinearOpBase<double>>
MockModelEval_B_Tpetra_2_parameters::create_hess_g_pp( int j, int l1, int l2 ) const
{
  const Teuchos::RCP<Tpetra_Operator> H =
      Teuchos::rcp(new Tpetra_CrsMatrix(hess_crs_graph));
  return Teuchos::rcp(new MatrixBased_LOWS(Thyra::createLinearOp(H)));
}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_B_Tpetra_2_parameters::getNominalValues() const
{
  return nominalValues;
}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_B_Tpetra_2_parameters::getLowerBounds() const
{
  return lowerBounds;
}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_B_Tpetra_2_parameters::getUpperBounds() const
{
  return upperBounds;
}


Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_B_Tpetra_2_parameters::createInArgs() const
{
  return this->createInArgsImpl();
}

void
MockModelEval_B_Tpetra_2_parameters::reportFinalPoint(
    const Thyra::ModelEvaluatorBase::InArgs<double>& /* finalPoint */,
    const bool /* wasSolved */) {
  // Do nothing  
}

Thyra::ModelEvaluatorBase::OutArgs<double>
MockModelEval_B_Tpetra_2_parameters::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<double> result;
  result.setModelEvalDescription(this->description());
  result.set_Np_Ng(2, 1);

  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f, true);
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op, true);
  result.set_W_properties(Thyra::ModelEvaluatorBase::DerivativeProperties(
      Thyra::ModelEvaluatorBase::DERIV_LINEARITY_UNKNOWN,
      Thyra::ModelEvaluatorBase::DERIV_RANK_FULL,
      true));

  result.setSupports(
      Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, 0, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
  result.setSupports(
      Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, 1, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
  result.setSupports(
      Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, 0, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, 0, 0, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, 0, 1, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);

  result.setSupports(
        Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, 0, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);

  result.setHessianSupports(hessSupport);

  return result;
}

void MockModelEval_B_Tpetra_2_parameters::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<double>&  inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<double>& outArgs) const
{

  // Parse InArgs

  const Teuchos::RCP<const Tpetra_Vector> x_in =
      ConverterT::getConstTpetraVector(inArgs.get_x());
  if (!Teuchos::nonnull(x_in)) std::cerr << "ERROR: MockModelEval_B_Tpetra_2_parameters requires x as inargs\n";

  const Teuchos::RCP<const Tpetra_Vector> x_dot_in =
      Teuchos::nonnull(inArgs.get_x_dot()) ?
          ConverterT::getConstTpetraVector(inArgs.get_x_dot()) :
          Teuchos::null;

  const Teuchos::RCP<const Thyra::VectorBase<double>> p_in_0 = inArgs.get_p(0);
  if (Teuchos::nonnull(p_in_0))
    p_vec_0->assign(*ConverterT::getConstTpetraVector(p_in_0));

  const Teuchos::RCP<const Thyra::VectorBase<double>> p_in_1 = inArgs.get_p(1);
  if (Teuchos::nonnull(p_in_1))
    p_vec_1->assign(*ConverterT::getConstTpetraVector(p_in_1));

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

  const Teuchos::RCP<Thyra::MultiVectorBase<double>> dfdp_0_base =
      outArgs.get_DfDp(0).getMultiVector();

  const Teuchos::RCP<Tpetra_MultiVector> dfdp_0_out =
      Teuchos::nonnull(dfdp_0_base) ?
          ConverterT::getTpetraMultiVector(dfdp_0_base) :
          Teuchos::null;

  const Teuchos::RCP<Thyra::MultiVectorBase<double>> dfdp_1_base =
      outArgs.get_DfDp(1).getMultiVector();

  const Teuchos::RCP<Tpetra_MultiVector> dfdp_1_out =
      Teuchos::nonnull(dfdp_1_base) ?
          ConverterT::getTpetraMultiVector(dfdp_1_base) :
          Teuchos::null;

  const Teuchos::RCP<Thyra::MultiVectorBase<double>> dgdp_0_base =
      outArgs.get_DgDp(0, 0).getMultiVector();
  const Teuchos::RCP<Tpetra_MultiVector> dgdp_0_out =
      Teuchos::nonnull(dgdp_0_base) ?
          ConverterT::getTpetraMultiVector(dgdp_0_base) :
          Teuchos::null;

  const Teuchos::RCP<Thyra::MultiVectorBase<double>> dgdp_1_base =
      outArgs.get_DgDp(0, 1).getMultiVector();
  const Teuchos::RCP<Tpetra_MultiVector> dgdp_1_out =
      Teuchos::nonnull(dgdp_1_base) ?
          ConverterT::getTpetraMultiVector(dgdp_1_base) :
          Teuchos::null;

  const Teuchos::RCP<Thyra::MultiVectorBase<double>> dgdx_base =
        outArgs.get_DgDx(0).getMultiVector();
  const Teuchos::RCP<Tpetra_MultiVector> dgdx_out =
      Teuchos::nonnull(dgdx_base) ?
          ConverterT::getTpetraMultiVector(dgdx_base) :
          Teuchos::null;

  const Teuchos::RCP<const Tpetra_MultiVector> p_direction_0 =
      Teuchos::nonnull(inArgs.get_p_direction(0)) ?
        ConverterT::getConstTpetraMultiVector(inArgs.get_p_direction(0)):
        Teuchos::null;
  const Teuchos::RCP<const Tpetra_MultiVector> p_direction_1 =
      Teuchos::nonnull(inArgs.get_p_direction(1)) ?
        ConverterT::getConstTpetraMultiVector(inArgs.get_p_direction(1)):
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

  auto f_hess_xp_0_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_f_xp,0) ? outArgs.get_hess_vec_prod_f_xp(0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_xp_0_v_out =
      Teuchos::nonnull(f_hess_xp_0_v) ?
        ConverterT::getTpetraMultiVector(f_hess_xp_0_v) :
        Teuchos::null;

  auto f_hess_xp_1_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_f_xp,1) ? outArgs.get_hess_vec_prod_f_xp(1) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_xp_1_v_out =
      Teuchos::nonnull(f_hess_xp_1_v) ?
        ConverterT::getTpetraMultiVector(f_hess_xp_1_v) :
        Teuchos::null;

  auto f_hess_px_0_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_f_px,0) ? outArgs.get_hess_vec_prod_f_px(0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_px_0_v_out =
      Teuchos::nonnull(f_hess_px_0_v) ?
        ConverterT::getTpetraMultiVector(f_hess_px_0_v) :
        Teuchos::null;

  auto f_hess_px_1_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_f_px,1) ? outArgs.get_hess_vec_prod_f_px(1) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_px_1_v_out =
      Teuchos::nonnull(f_hess_px_1_v) ?
        ConverterT::getTpetraMultiVector(f_hess_px_1_v) :
        Teuchos::null;

  auto f_hess_pp_00_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_f_pp,0,0) ? outArgs.get_hess_vec_prod_f_pp(0,0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_pp_00_v_out =
      Teuchos::nonnull(f_hess_pp_00_v) ?
        ConverterT::getTpetraMultiVector(f_hess_pp_00_v) :
        Teuchos::null;

  auto f_hess_pp_01_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_f_pp,0,1) ? outArgs.get_hess_vec_prod_f_pp(0,1) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_pp_01_v_out =
      Teuchos::nonnull(f_hess_pp_01_v) ?
        ConverterT::getTpetraMultiVector(f_hess_pp_01_v) :
        Teuchos::null;

  auto f_hess_pp_10_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_f_pp,1,0) ? outArgs.get_hess_vec_prod_f_pp(1,0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_pp_10_v_out =
      Teuchos::nonnull(f_hess_pp_10_v) ?
        ConverterT::getTpetraMultiVector(f_hess_pp_10_v) :
        Teuchos::null;

  auto f_hess_pp_11_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_f_pp,1,1) ? outArgs.get_hess_vec_prod_f_pp(1,1) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_pp_11_v_out =
      Teuchos::nonnull(f_hess_pp_11_v) ?
        ConverterT::getTpetraMultiVector(f_hess_pp_11_v) :
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

  auto g_hess_xp_0_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_g_xp,0,0) ? outArgs.get_hess_vec_prod_g_xp(0,0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_xp_0_v_out =
      Teuchos::nonnull(g_hess_xp_0_v) ?
        ConverterT::getTpetraMultiVector(g_hess_xp_0_v) :
        Teuchos::null;

  auto g_hess_xp_1_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_g_xp,0,1) ? outArgs.get_hess_vec_prod_g_xp(0,1) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_xp_1_v_out =
      Teuchos::nonnull(g_hess_xp_1_v) ?
        ConverterT::getTpetraMultiVector(g_hess_xp_1_v) :
        Teuchos::null;

  auto g_hess_px_0_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_g_px,0,0) ? outArgs.get_hess_vec_prod_g_px(0,0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_px_0_v_out =
      Teuchos::nonnull(g_hess_px_0_v) ?
        ConverterT::getTpetraMultiVector(g_hess_px_0_v) :
        Teuchos::null;

  auto g_hess_px_1_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_g_px,0,1) ? outArgs.get_hess_vec_prod_g_px(0,1) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_px_1_v_out =
      Teuchos::nonnull(g_hess_px_1_v) ?
        ConverterT::getTpetraMultiVector(g_hess_px_1_v) :
        Teuchos::null;

  auto g_hess_pp_00_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_g_pp,0,0,0) ? outArgs.get_hess_vec_prod_g_pp(0,0,0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_pp_00_v_out =
      Teuchos::nonnull(g_hess_pp_00_v) ?
        ConverterT::getTpetraMultiVector(g_hess_pp_00_v) :
        Teuchos::null;

  auto g_hess_pp_01_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_g_pp,0,0,1) ? outArgs.get_hess_vec_prod_g_pp(0,0,1) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_pp_01_v_out =
      Teuchos::nonnull(g_hess_pp_01_v) ?
        ConverterT::getTpetraMultiVector(g_hess_pp_01_v) :
        Teuchos::null;

  auto g_hess_pp_10_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_g_pp,0,1,0) ? outArgs.get_hess_vec_prod_g_pp(0,1,0) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_pp_10_v_out =
      Teuchos::nonnull(g_hess_pp_10_v) ?
        ConverterT::getTpetraMultiVector(g_hess_pp_10_v) :
        Teuchos::null;

  auto g_hess_pp_11_v = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_vec_prod_g_pp,0,1,1) ? outArgs.get_hess_vec_prod_g_pp(0,1,1) : Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_pp_11_v_out =
      Teuchos::nonnull(g_hess_pp_11_v) ?
        ConverterT::getTpetraMultiVector(g_hess_pp_11_v) :
        Teuchos::null;

  auto x = x_in->getData();
  auto p_0 = p_vec_0->getData();
  auto p_1 = p_vec_1->getData();

  if (f_out != Teuchos::null) {
    f_out->putScalar(0.0);
    auto f_out_data = f_out->getDataNonConst();
    for (int i=0; i<myVecLength; i++)
      f_out_data[i] = x[i];
  }
  if (W_out != Teuchos::null) {
    Teuchos::RCP<Tpetra_CrsMatrix> W_out_crs =
      Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(W_out, true);
    W_out_crs->resumeFill();
    W_out_crs->setAllToScalar(0.0);

    double diag=1.0;
    for (int i=0; i<myVecLength; i++)
      W_out_crs->replaceLocalValues(i, 1, &diag, &i);
    if(!Teuchos::nonnull(x_dot_in))
      W_out_crs->fillComplete();
  }

  auto hess_g_pp_00 = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_g_pp,0,0,0) ? outArgs.get_hess_g_pp(0,0,0) : Teuchos::null; 
  const Teuchos::RCP<MatrixBased_LOWS> H_pp_00_out =
    Teuchos::nonnull(hess_g_pp_00) ?
      Teuchos::rcp_dynamic_cast<MatrixBased_LOWS>(outArgs.get_hess_g_pp(0,0,0)):
      Teuchos::null;

  auto hess_g_pp_01 = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_g_pp,0,0,1) ? outArgs.get_hess_g_pp(0,0,1) : Teuchos::null; 
  const Teuchos::RCP<MatrixBased_LOWS> H_pp_01_out =
    Teuchos::nonnull(hess_g_pp_01) ?
      Teuchos::rcp_dynamic_cast<MatrixBased_LOWS>(outArgs.get_hess_g_pp(0,0,1)):
      Teuchos::null;

  auto hess_g_pp_10 = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_g_pp,0,1,0) ? outArgs.get_hess_g_pp(0,1,0) : Teuchos::null; 
  const Teuchos::RCP<MatrixBased_LOWS> H_pp_10_out =
    Teuchos::nonnull(hess_g_pp_10) ?
      Teuchos::rcp_dynamic_cast<MatrixBased_LOWS>(outArgs.get_hess_g_pp(0,1,0)):
      Teuchos::null;

  auto hess_g_pp_11 = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_g_pp,0,1,1) ? outArgs.get_hess_g_pp(0,1,1) : Teuchos::null; 
  const Teuchos::RCP<MatrixBased_LOWS> H_pp_11_out =
    Teuchos::nonnull(hess_g_pp_11) ?
      Teuchos::rcp_dynamic_cast<MatrixBased_LOWS>(outArgs.get_hess_g_pp(0,1,1)):
      Teuchos::null;

  // Response: g = 0.5*(p0-6)^2 + 0.5*c*(p1-4)^2 + 0.5*(p0+p1-10)^2
  // min g(x(p), p) s.t. f(x, p) = 0 reached for p0 = 6, p1 = 4

  double term1, term2, term3, c;
  term1 = p_0[0]-6;
  term2 = p_1[0]-4;
  term3 = p_0[0]+p_1[0]-10;
  c = 5;

  if (Teuchos::nonnull(H_pp_00_out)) {
    Teuchos::RCP<Tpetra_CrsMatrix> H_pp_out_crs =
      Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(ConverterT::getTpetraOperator(H_pp_00_out->getMatrix()), true);
    H_pp_out_crs->resumeFill();
    H_pp_out_crs->setAllToScalar(0.0);

    if (comm->getRank() == 0) {
      std::vector<double> vals = {2};
      std::vector<typename Tpetra_CrsGraph::global_ordinal_type> indices = {0};
      H_pp_out_crs->replaceGlobalValues(0, 1, &vals[0], &indices[0]);
    }
    H_pp_out_crs->fillComplete();

    if(probList_->sublist("Hessian").sublist("Response 0").sublist("Parameter 0").isSublist("H_pp Solver")) {
      auto pl = probList_->sublist("Hessian").sublist("Response 0").sublist("Parameter 0").sublist("H_pp Solver");
      H_pp_00_out->initializeSolver(Teuchos::rcpFromRef(pl));
    }
  }

  if (Teuchos::nonnull(H_pp_01_out)) {
    Teuchos::RCP<Tpetra_CrsMatrix> H_pp_out_crs =
      Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(ConverterT::getTpetraOperator(H_pp_01_out->getMatrix()), true);
    H_pp_out_crs->resumeFill();
    H_pp_out_crs->setAllToScalar(0.0);

    if (comm->getRank() == 0) {
      std::vector<double> vals = {1};
      std::vector<typename Tpetra_CrsGraph::global_ordinal_type> indices = {0};
      H_pp_out_crs->replaceGlobalValues(0, 1, &vals[0], &indices[0]);
    }
    H_pp_out_crs->fillComplete();

    if(probList_->sublist("Hessian").sublist("Response 0").sublist("Parameter 0").isSublist("H_pp Solver")) {
      auto pl = probList_->sublist("Hessian").sublist("Response 0").sublist("Parameter 0").sublist("H_pp Solver");
      H_pp_01_out->initializeSolver(Teuchos::rcpFromRef(pl));
    }
  }

  if (Teuchos::nonnull(H_pp_10_out)) {
    Teuchos::RCP<Tpetra_CrsMatrix> H_pp_out_crs =
      Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(ConverterT::getTpetraOperator(H_pp_10_out->getMatrix()), true);
    H_pp_out_crs->resumeFill();
    H_pp_out_crs->setAllToScalar(0.0);

    if (comm->getRank() == 0) {
      std::vector<double> vals = {1};
      std::vector<typename Tpetra_CrsGraph::global_ordinal_type> indices = {0};
      H_pp_out_crs->replaceGlobalValues(0, 1, &vals[0], &indices[0]);
    }
    H_pp_out_crs->fillComplete();

    if(probList_->sublist("Hessian").sublist("Response 0").sublist("Parameter 0").isSublist("H_pp Solver")) {
      auto pl = probList_->sublist("Hessian").sublist("Response 0").sublist("Parameter 0").sublist("H_pp Solver");
      H_pp_10_out->initializeSolver(Teuchos::rcpFromRef(pl));
    }
  }

  if (Teuchos::nonnull(H_pp_11_out)) {
    Teuchos::RCP<Tpetra_CrsMatrix> H_pp_out_crs =
      Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(ConverterT::getTpetraOperator(H_pp_11_out->getMatrix()), true);
    H_pp_out_crs->resumeFill();
    H_pp_out_crs->setAllToScalar(0.0);

    if (comm->getRank() == 0) {
      std::vector<double> vals = {1+c};
      std::vector<typename Tpetra_CrsGraph::global_ordinal_type> indices = {0};
      H_pp_out_crs->replaceGlobalValues(0, 1, &vals[0], &indices[0]);
    }
    H_pp_out_crs->fillComplete();

    if(probList_->sublist("Hessian").sublist("Response 0").sublist("Parameter 0").isSublist("H_pp Solver")) {
      auto pl = probList_->sublist("Hessian").sublist("Response 0").sublist("Parameter 0").sublist("H_pp Solver");
      H_pp_11_out->initializeSolver(Teuchos::rcpFromRef(pl));
    }
  }

  if (Teuchos::nonnull(dfdp_0_out)) {
    dfdp_0_out->putScalar(0.0);
  }

  if (Teuchos::nonnull(dfdp_1_out)) {
    dfdp_1_out->putScalar(0.0);
  }

  if (Teuchos::nonnull(g_out)) {
    g_out->getDataNonConst()[0] = 0.5*term1*term1 + 0.5*c*term2*term2 + 0.5*term3*term3;
  }

  if (dgdx_out != Teuchos::null) {
    dgdx_out->putScalar(0);
  }
  if (dgdp_0_out != Teuchos::null) {
    dgdp_0_out->putScalar(0.0);
    dgdp_0_out->getVectorNonConst(0)->getDataNonConst()[0] = term1+term3;
  }
  if (dgdp_1_out != Teuchos::null) {
    dgdp_1_out->putScalar(0.0);
    dgdp_1_out->getVectorNonConst(0)->getDataNonConst()[0] = c*term2+term3;
  }

  if (Teuchos::nonnull(f_hess_xx_v_out)) {
    f_hess_xx_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_xp_0_v_out)) {
    f_hess_xp_0_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_xp_1_v_out)) {
    f_hess_xp_1_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_px_0_v_out)) {
    f_hess_px_0_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_px_1_v_out)) {
    f_hess_px_1_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_pp_00_v_out)) {
    f_hess_pp_00_v_out->getVectorNonConst(0)->putScalar(0);
  }
  if (Teuchos::nonnull(f_hess_pp_01_v_out)) {
    f_hess_pp_01_v_out->getVectorNonConst(0)->putScalar(0);
  }
  if (Teuchos::nonnull(f_hess_pp_10_v_out)) {
    f_hess_pp_10_v_out->getVectorNonConst(0)->putScalar(0);
  }
  if (Teuchos::nonnull(f_hess_pp_11_v_out)) {
    f_hess_pp_11_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_xx_v_out)) {
    g_hess_xx_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_xp_0_v_out)) {
    g_hess_xp_0_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_xp_1_v_out)) {
    g_hess_xp_1_v_out->getVectorNonConst(0)->putScalar(0);
  }
 
  if (Teuchos::nonnull(g_hess_px_0_v_out)) {
    g_hess_px_0_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_px_1_v_out)) {
    g_hess_px_1_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_pp_00_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(p_direction_0));
    const auto direction_p_0 = p_direction_0->getVector(0)->getData();
    g_hess_pp_00_v_out->getVectorNonConst(0)->getDataNonConst()[0] = 2*direction_p_0[0];
  }
  if (Teuchos::nonnull(g_hess_pp_01_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(p_direction_1));
    const auto direction_p_1 = p_direction_1->getVector(0)->getData();
    g_hess_pp_01_v_out->getVectorNonConst(0)->getDataNonConst()[0] = direction_p_1[0];
  }
  if (Teuchos::nonnull(g_hess_pp_10_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(p_direction_0));
    const auto direction_p_0 = p_direction_0->getVector(0)->getData();
    g_hess_pp_10_v_out->getVectorNonConst(0)->getDataNonConst()[0] = direction_p_0[0];
  }
  if (Teuchos::nonnull(g_hess_pp_11_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(p_direction_1));
    const auto direction_p_1 = p_direction_1->getVector(0)->getData();
    g_hess_pp_11_v_out->getVectorNonConst(0)->getDataNonConst()[0] = (c+1)*direction_p_1[0];
  }  
  // Modify for time dependent (implicit time integration or eigensolves)
  if (Teuchos::nonnull(x_dot_in)) {
    // Velocity provided: Time dependent problem
    double alpha = inArgs.get_alpha();
    double beta = inArgs.get_beta();
    if (alpha==0.0 && beta==0.0) {
      std::cerr << "MockModelEval Warning: alpha=beta=0 -- setting beta=1\n";
      beta = 1.0;
    }

    if (f_out != Teuchos::null) {
      // f(x, x_dot) = f(x) - x_dot
      auto f_out_data = f_out->getDataNonConst();
      for (int i=0; i<myVecLength; i++) {
        f_out_data[i] = -x_dot_in->getData()[i] + f_out->getData()[i];
      }
    }
    if (W_out != Teuchos::null) {
      // W(x, x_dot) = beta * W(x) - alpha * Id
      const Teuchos::RCP<Tpetra_CrsMatrix> W_out_crs =
        Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(W_out, true);
      W_out_crs->scale(beta);

      const double diag = -alpha;
      for (int i=0; i<myVecLength; i++) {
        W_out_crs->sumIntoLocalValues(i, 1, &diag, &i);
      }
      W_out_crs->fillComplete();
    }
  }
}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_B_Tpetra_2_parameters::createInArgsImpl() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<double> result;
  result.setModelEvalDescription(this->description());

  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x, true);


  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_t, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_alpha, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_beta, true);

  result.set_Np_Ng(2,1);

  return result;
}

