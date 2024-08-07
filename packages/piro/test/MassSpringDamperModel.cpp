// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MassSpringDamperModel.hpp"
#include "Thyra_LinearOpWithSolveBase_decl.hpp"

#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"

#include <iostream>

using Teuchos::RCP;
using Teuchos::rcp;

MassSpringDamperModel::MassSpringDamperModel(const Teuchos::RCP<const Teuchos::Comm<int> >  appComm, bool adjoint, const Teuchos::RCP<Teuchos::ParameterList>& problemList, bool hessianSupport, bool use_x_dot_in_g)
 {
    comm = appComm;
    hessSupport = hessianSupport;
    adjoint_ = adjoint;
    use_x_dot_in_g_ = use_x_dot_in_g;

    target_x_ = 1;
    target_x_dot_ = 0;
    scaling_ = 0.02;

    scaling_g_x_ = 1.;
    scaling_g_p_ = 100.;

    target_k_ = 1;
    target_m_ = 0.5;

    //set up map and initial guess for solution vector
    const int vecLength = 2;
    x_map = rcp(new Tpetra_Map(vecLength, 0, comm, Tpetra::LocallyReplicated));
    x_vec = rcp(new Tpetra_Vector(x_map));
    x_dot_vec = rcp(new Tpetra_Vector(x_map));

    x_vec->putScalar(0.0);
    x_dot_vec->putScalar(0.0);
    x_dot_vec->getDataNonConst()[x_map->getLocalElement(1)]= 1.0;                                  // F/m with F == m == 1

    Teuchos::RCP<const Thyra::VectorSpaceBase<double>> x_space =
        Thyra::createVectorSpace<double>(x_map);

    //set up responses
    const int numResponses = 1;
    g_map = rcp(new const Tpetra_Map(numResponses , 0, comm, Tpetra::LocallyReplicated));

    //set up parameters
    const int numParameters= 2;
    p_map = rcp(new const Tpetra_Map(numParameters, 0, comm, Tpetra::LocallyReplicated));

    Teuchos::RCP<const Thyra::VectorSpaceBase<double>> p_space =
        Thyra::createVectorSpace<double>(p_map);

    Teuchos::RCP<Tpetra_Vector> p_init = rcp(new Tpetra_Vector(p_map));
    Teuchos::RCP<Tpetra_Vector> p_lo = rcp(new Tpetra_Vector(p_map));
    Teuchos::RCP<Tpetra_Vector> p_up = rcp(new Tpetra_Vector(p_map));
    for (int i=0; i<numParameters; i++) {
      p_init->getDataNonConst()[i]= 1.0;
      p_lo->getDataNonConst()[i]= 0.5;
      p_up->getDataNonConst()[i]= 1.5;
    }

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

    nominalValues.set_p(0, Thyra::createVector(p_init, p_space));
    lowerBounds.set_p(0, Thyra::createVector(p_lo, p_space));
    upperBounds.set_p(0, Thyra::createVector(p_up, p_space));

    probList_ = problemList;
}

MassSpringDamperModel::~MassSpringDamperModel()
{
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MassSpringDamperModel::get_x_space() const
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> x_space =
      Thyra::createVectorSpace<double>(x_map);
  return x_space;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MassSpringDamperModel::get_f_space() const
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> f_space =
      Thyra::createVectorSpace<double>(x_map);
  return f_space;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MassSpringDamperModel::get_p_space(int l) const
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
MassSpringDamperModel::get_g_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  MassSpringDamperModel::get_g_map() only " <<
                     " supports 1 response.  Supplied index l = " <<
                     l << std::endl);
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> g_space =
        Thyra::createVectorSpace<double>(g_map);
  return g_space;
}

RCP<const  Teuchos::Array<std::string> > MassSpringDamperModel::get_p_names(int l) const
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
MassSpringDamperModel::create_W_op() const
{
  const Teuchos::RCP<Tpetra_Operator> W =
      Teuchos::rcp(new Tpetra_CrsMatrix(crs_graph));
  return Thyra::createLinearOp(W);
}

//! Create preconditioner operator
Teuchos::RCP<Thyra::PreconditionerBase<double>>
MassSpringDamperModel::create_W_prec() const
{
  return Teuchos::null;
}

Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double>>
MassSpringDamperModel::get_W_factory() const
{
  return Teuchos::null;
}

Teuchos::RCP<Thyra::LinearOpBase<double>>
MassSpringDamperModel::create_hess_g_pp( int j, int l1, int l2 ) const
{
  const Teuchos::RCP<Tpetra_Operator> H =
      Teuchos::rcp(new Tpetra_CrsMatrix(hess_crs_graph));
  return Teuchos::rcp(new MatrixBased_LOWS(Thyra::createLinearOp(H)));
}

Thyra::ModelEvaluatorBase::InArgs<double>
MassSpringDamperModel::getNominalValues() const
{
  return nominalValues;
}

Thyra::ModelEvaluatorBase::InArgs<double>
MassSpringDamperModel::getLowerBounds() const
{
  return lowerBounds;
}

Thyra::ModelEvaluatorBase::InArgs<double>
MassSpringDamperModel::getUpperBounds() const
{
  return upperBounds;
}


Thyra::ModelEvaluatorBase::InArgs<double>
MassSpringDamperModel::createInArgs() const
{
  return this->createInArgsImpl();
}

void
MassSpringDamperModel::reportFinalPoint(
    const Thyra::ModelEvaluatorBase::InArgs<double>& /* finalPoint */,
    const bool /* wasSolved */) {
  // Do nothing  
}

Thyra::ModelEvaluatorBase::OutArgs<double>
MassSpringDamperModel::createOutArgsImpl() const
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
  if (use_x_dot_in_g_) {
    result.setSupports(
        Thyra::ModelEvaluatorBase::OUT_ARG_DgDx_dot, 0, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
  }
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, 0, 0, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);

  result.setHessianSupports(hessSupport);

  return result;
}

void MassSpringDamperModel::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<double>&  inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<double>& outArgs) const
{
  // Parse InArgs

  const Teuchos::RCP<const Tpetra_Vector> x_in =
      ConverterT::getConstTpetraVector(inArgs.get_x());
  if (!Teuchos::nonnull(x_in)) std::cerr << "ERROR: MassSpringDamperModel requires x as inargs\n";

  const Teuchos::RCP<const Tpetra_Vector> x_dot_in =
      Teuchos::nonnull(inArgs.get_x_dot()) ?
          ConverterT::getConstTpetraVector(inArgs.get_x_dot()) :
          Teuchos::null;

  const Teuchos::RCP<const Thyra::VectorBase<double>> p_in = inArgs.get_p(0);
  if (Teuchos::nonnull(p_in)) {
    Teuchos::RCP<const Thyra::ProductVectorBase<double>> p_prod_in =
      Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<double>>(p_in);
    if(Teuchos::nonnull(p_prod_in)) {
      if(p_prod_in->productSpace()->numBlocks() == 1) {
        p_vec->assign(*ConverterT::getConstTpetraVector(p_prod_in->getVectorBlock(0)));
      } else {
        std::cerr << "ERROR: MassSpringDamperModel has a parameter with " << p_prod_in->productSpace()->numBlocks() << " blocks \n";
      }
    } else {
      p_vec->assign(*ConverterT::getConstTpetraVector(p_in));
    }
  }

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

  const Teuchos::RCP<Thyra::MultiVectorBase<double>> dgdx_dot_base =
      use_x_dot_in_g_ ? 
          outArgs.get_DgDx_dot(0).getMultiVector() :
          Teuchos::null;
  const Teuchos::RCP<Tpetra_MultiVector> dgdx_dot_out =
      Teuchos::nonnull(dgdx_dot_base) ?
          ConverterT::getTpetraMultiVector(dgdx_dot_base) :
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

  auto x = x_in->getData();
  auto x_dot = x_dot_in->getData();
  auto p = p_vec->getData();

  auto k = p[0];
  auto m = p[1];

  double F = 1;

  if (f_out != Teuchos::null) {
    f_out->putScalar(0.0);

    f_out->getDataNonConst()[x_map->getLocalElement(0)]= x[x_map->getLocalElement(1)];
    f_out->getDataNonConst()[x_map->getLocalElement(1)]= - (2*sqrt(m*k)*x[x_map->getLocalElement(1)] + k*x[x_map->getLocalElement(0)] - F ) / m;
  }
  if (W_out != Teuchos::null) {
    Teuchos::RCP<Tpetra_CrsMatrix> W_out_crs =
      Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(W_out, true);
    W_out_crs->resumeFill();
    W_out_crs->setAllToScalar(0.0);

    double val;
    for (int row=0; row<myVecLength; ++row) {
      for (int col=0; col<myVecLength; ++col) {
        if ( row == 0 && col == 0)
          val = 0.0;                // d(f0)/d(x0_n)
        if ( (row == 0 && col == 1 && !adjoint_) || (row == 1 && col == 0 && adjoint_))
          val = 1.0;                // d(f0)/d(x1_n)
        if ( (row == 1 && col == 0 && !adjoint_) || (row == 0 && col == 1 && adjoint_))
          val = -(k/m);             // d(f1)/d(x0_n)
        if ( row == 1 && col == 1)
          val = -2*sqrt(k/m);       // d(f1)/d(x1_n)
        W_out_crs->replaceLocalValues(row, 1, &val, &col);
      }
    }
    W_out_crs->fillComplete();
  }

  auto hess_g_pp = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_g_pp,0,0,0) ? outArgs.get_hess_g_pp(0,0,0) : Teuchos::null; 
  const Teuchos::RCP<MatrixBased_LOWS> H_pp_out =
    Teuchos::nonnull(hess_g_pp) ?
      Teuchos::rcp_dynamic_cast<MatrixBased_LOWS>(outArgs.get_hess_g_pp(0,0,0)):
      Teuchos::null;

  // Response:  g = scaling_g_x * (( x - target_x )^2 + scaling ( x_dot - target_x_dot )^2) 
  //                 + scaling_g_p * (( k - target_k )^2 + ( m - target_m )^2)

  double diff_x = (x[0] - target_x_);
  double diff_x_dot;
  if (use_x_dot_in_g_)
    diff_x_dot = (x_dot[0] - target_x_dot_);
  else
    diff_x_dot = (x[1] - target_x_dot_);

  double diff_k = (p[0] - target_k_);
  double diff_m = (p[1] - target_m_);

  if (Teuchos::nonnull(H_pp_out)) {
    Teuchos::RCP<Tpetra_CrsMatrix> H_pp_out_crs =
      Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(ConverterT::getTpetraOperator(H_pp_out->getMatrix()), true);
    H_pp_out_crs->resumeFill();
    H_pp_out_crs->setAllToScalar(0.0);
    H_pp_out_crs->fillComplete();

    if(probList_->sublist("Hessian").sublist("Response 0").sublist("Parameter 0").isSublist("H_pp Solver")) {
      auto pl = probList_->sublist("Hessian").sublist("Response 0").sublist("Parameter 0").sublist("H_pp Solver");
      H_pp_out->initializeSolver(Teuchos::rcpFromRef(pl));
    }
  }

  if (Teuchos::nonnull(dfdp_out)) {
    dfdp_out->putScalar(0.0);
    auto dfdp_out_data_0 = dfdp_out->getVectorNonConst(0)->getDataNonConst();
    auto dfdp_out_data_1 = dfdp_out->getVectorNonConst(1)->getDataNonConst();

    dfdp_out_data_0[0] = 0.0;
    dfdp_out_data_1[0] = 0.0;
    dfdp_out_data_0[1] = -1/sqrt(m*k) * x[1] - 1/m * x[0];
    dfdp_out_data_1[1] = sqrt(m*k)/std::pow(m,2) * x[1] + k/std::pow(m,2) * x[0] - F/std::pow(m,2);
  }

  if (Teuchos::nonnull(g_out)) {
    g_out->getDataNonConst()[0] = scaling_g_x_ * (diff_x*diff_x + scaling_ * diff_x_dot*diff_x_dot)
    + scaling_g_p_ * (diff_k*diff_k + diff_m*diff_m);
  }

  if (dgdx_out != Teuchos::null) {
    dgdx_out->getVectorNonConst(0)->getDataNonConst()[0] = scaling_g_x_*2*diff_x;
    if (use_x_dot_in_g_)
      dgdx_out->getVectorNonConst(0)->getDataNonConst()[1] = 0.0;
    else
      dgdx_out->getVectorNonConst(0)->getDataNonConst()[1] = scaling_g_x_*scaling_*2*diff_x_dot;
  }

  if (dgdx_dot_out != Teuchos::null) {
    if (use_x_dot_in_g_)
      dgdx_dot_out->getVectorNonConst(0)->getDataNonConst()[0] = scaling_g_x_*scaling_*2*diff_x_dot;
    else
      dgdx_dot_out->getVectorNonConst(0)->getDataNonConst()[0] = 0.0;
    dgdx_dot_out->getVectorNonConst(0)->getDataNonConst()[1] = 0.0;
  }

  if (dgdp_out != Teuchos::null) {
    dgdp_out->putScalar(0.0);
    dgdp_out->getVectorNonConst(0)->getDataNonConst()[0] = scaling_g_p_*2*diff_k;
    dgdp_out->getVectorNonConst(0)->getDataNonConst()[1] = scaling_g_p_*2*diff_m;
  }

  if (Teuchos::nonnull(f_hess_xx_v_out)) {
    f_hess_xx_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_xp_v_out)) {
    f_hess_xp_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_px_v_out)) {
    f_hess_px_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_pp_v_out)) {
    f_hess_pp_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_xx_v_out)) {
    g_hess_xx_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_xp_v_out)) {
    g_hess_xp_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_px_v_out)) {
    g_hess_px_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_pp_v_out)) {
    g_hess_pp_v_out->getVectorNonConst(0)->putScalar(0);
  }

  // Modify for time dependent (implicit time integration or eigensolves)
  if (Teuchos::nonnull(x_dot_in)) {
    // Velocity provided: Time dependent problem

    if (f_out != Teuchos::null) {
      // f(x, x_dot) = f(x) - x_dot
      auto f_out_data = f_out->getDataNonConst();
      for (int i=0; i<myVecLength; i++) {
        f_out_data[i] -= x_dot_in->getData()[i];
      }
    }
    if (W_out != Teuchos::null) {
      double alpha = inArgs.get_alpha();
      double beta = inArgs.get_beta();
      if (alpha==0.0 && beta==0.0) {
        std::cerr << "MockModelEval Warning: alpha=beta=0 -- setting beta=1\n";
        beta = 1.0;
      }

      // W(x, x_dot) = beta * W(x) - alpha * Id
      const Teuchos::RCP<Tpetra_CrsMatrix> W_out_crs =
        Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(W_out, true);
      W_out_crs->resumeFill();
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
MassSpringDamperModel::createInArgsImpl() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<double> result;
  result.setModelEvalDescription(this->description());

  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x, true);


  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_t, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_alpha, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_beta, true);

  result.set_Np_Ng(1,1);

  return result;
}

