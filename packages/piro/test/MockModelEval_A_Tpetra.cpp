// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MockModelEval_A_Tpetra.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"

#include <iostream>

using Teuchos::RCP;
using Teuchos::rcp;

MockModelEval_A_Tpetra::MockModelEval_A_Tpetra(const Teuchos::RCP<const Teuchos::Comm<int> >  appComm, int paramVecDim, bool adjoint, const Teuchos::RCP<Teuchos::ParameterList>& problemList, bool hessianSupport)
{
    comm = appComm;
    adjointModel = adjoint;
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
    TEUCHOS_TEST_FOR_EXCEPTION((paramVecDim < 1) || (paramVecDim > 2), std::logic_error,
                     std::endl <<
                     "Error!  MockModelEval_A_Tpetra::dimension of parameer space should be 1 or 2. Dimension provided: " << paramVecDim << std::endl);

    p_map = rcp(new const Tpetra_Map(paramVecDim, 0, comm, Tpetra::LocallyReplicated));

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
    hess_crs_graph = rcp(new Tpetra_CrsGraph(p_map, paramVecDim));
    if (comm->getRank() == 0)
    {
      std::vector<typename Tpetra_CrsGraph::global_ordinal_type> indices(paramVecDim);
      for (int i=0; i<paramVecDim; i++) indices[i]=i;
      const int nodeNumElements = p_map->getLocalNumElements();
      for (int i=0; i<nodeNumElements; i++)
        hess_crs_graph->insertGlobalIndices(p_map->getGlobalElement(i), paramVecDim, &indices[0]);
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

MockModelEval_A_Tpetra::~MockModelEval_A_Tpetra()
{
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_A_Tpetra::get_x_space() const
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> x_space =
      Thyra::createVectorSpace<double>(x_map);
  return x_space;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_A_Tpetra::get_f_space() const
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> f_space =
      Thyra::createVectorSpace<double>(x_map);
  return f_space;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_A_Tpetra::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  MockModelEval_A_Tpetra::get_p_space() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> p_space =
        Thyra::createVectorSpace<double>(p_map);
  return p_space;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_A_Tpetra::get_g_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  MockModelEval_A_Tpetra::get_g_space() only " <<
                     " supports 1 response.  Supplied index l = " <<
                     l << std::endl);
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> g_space =
        Thyra::createVectorSpace<double>(g_map);
  return g_space;
}

RCP<const  Teuchos::Array<std::string> > MockModelEval_A_Tpetra::get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  MockModelEval_A_Tpetra::get_p_names() only " <<
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
MockModelEval_A_Tpetra::create_W_op() const
{
  const Teuchos::RCP<Tpetra_Operator> W =
      Teuchos::rcp(new Tpetra_CrsMatrix(crs_graph));
  return Thyra::createLinearOp(W);
}

//! Create preconditioner operator
Teuchos::RCP<Thyra::PreconditionerBase<double>>
MockModelEval_A_Tpetra::create_W_prec() const
{
  return Teuchos::null;
}

Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double>>
MockModelEval_A_Tpetra::get_W_factory() const
{
  return Teuchos::null;
}

Teuchos::RCP<Thyra::LinearOpBase<double>>
MockModelEval_A_Tpetra::create_hess_g_pp( int j, int l1, int l2 ) const
{
  const Teuchos::RCP<Tpetra_Operator> H =
      Teuchos::rcp(new Tpetra_CrsMatrix(hess_crs_graph));
  return Teuchos::rcp(new MatrixBased_LOWS(Thyra::createLinearOp(H)));
}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_A_Tpetra::getNominalValues() const
{
  return nominalValues;
}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_A_Tpetra::getLowerBounds() const
{
  return lowerBounds;
}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_A_Tpetra::getUpperBounds() const
{
  return upperBounds;
}


Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_A_Tpetra::createInArgs() const
{
  return this->createInArgsImpl();
}

void
MockModelEval_A_Tpetra::reportFinalPoint(
    const Thyra::ModelEvaluatorBase::InArgs<double>& /* finalPoint */,
    const bool /* wasSolved */) {
  // Do nothing  
}

Thyra::ModelEvaluatorBase::OutArgs<double>
MockModelEval_A_Tpetra::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<double> result;
  result.setModelEvalDescription(this->description());
  result.set_Np_Ng(1, 1);

  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f, true);
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op, true);
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W, true);
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

void MockModelEval_A_Tpetra::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<double>&  inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<double>& outArgs) const
{

  // Parse InArgs

  const Teuchos::RCP<const Tpetra_Vector> x_in =
      ConverterT::getConstTpetraVector(inArgs.get_x());
  if (!Teuchos::nonnull(x_in)) std::cerr << "ERROR: MockModelEval_A_Tpetra requires x as inargs\n";

  const Teuchos::RCP<const Tpetra_Vector> x_dot_in =
      Teuchos::nonnull(inArgs.get_x_dot()) ?
          ConverterT::getConstTpetraVector(inArgs.get_x_dot()) :
          Teuchos::null;


  const Teuchos::RCP<const Thyra::VectorBase<double>> p_in_thyra = inArgs.get_p(0);
  if (Teuchos::nonnull(p_in_thyra)) {
    const Teuchos::RCP<const Tpetra_Vector> p_in = ConverterT::getConstTpetraVector(p_in_thyra);
    p_vec->assign(*p_in);
  }

  int num_p = p_map->getLocalNumElements();
  int vecLength = x_in->getGlobalLength();
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

  {
    auto x_host = x_in->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto p_host = p_vec->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto x = Kokkos::subview(x_host,Kokkos::ALL(),0);
    auto p = Kokkos::subview(p_host,Kokkos::ALL(),0);

    double x0;
    for (int i=0; i<myVecLength; i++) {
      if( x_in->getMap()->getGlobalElement(i) == 0) {
        x0 = x(i);
        TEUCHOS_ASSERT(comm->getRank() == 0);
      }
    }
    comm->broadcast(0, sizeof(double), (char*)(&x0));
    
    if (f_out != Teuchos::null) {
      auto f_out_view = f_out->getLocalViewHost(Tpetra::Access::OverwriteAll);
      for (int i=0; i<myVecLength; i++) {
        int gid = x_in->getMap()->getGlobalElement(i);
        
        if (gid==0) { // f_0 = (x_0)^3 - p_0
          f_out_view(i,0) = x(i) * x(i) * x(i) -  p(0);
        }
        else{ // f_i = x_i * (1 + x_0 - p_0^(1/3)) - (i+p_j) - 0.5*(x_0 - p_0),  (for i != 0);   j=1 if num_p>1, j=0 otherwise
          int j = (num_p > 1) ? 1 : 0;
          f_out_view(i,0) = x(i) - (gid + p(j)) - 0.5*(x0 - p(0)) + x(i) * (x0 - std::cbrt(p(0)));
        }
      }
    }

    if (W_out != Teuchos::null) {
      Teuchos::RCP<Tpetra_CrsMatrix> W_out_crs =
        Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(W_out, true);
      W_out_crs->resumeFill();
      W_out_crs->setAllToScalar(0.0);
      
      double diag=0, extra_diag=0;
      for (int i=0; i<myVecLength; i++) {
        auto gid = x_in->getMap()->getGlobalElement(i);
        if(gid==0) {
          diag = 3.0 * x0 * x0;
          W_out_crs->replaceLocalValues(i, 1, &diag, &i);
        } 
        else {
          diag = 1.0 + x0 - std::cbrt(p(0));
          W_out_crs->replaceLocalValues(i, 1, &diag, &i);
          decltype(gid) col = 0;
          extra_diag = -0.5 + x(i);
          W_out_crs->replaceGlobalValues(gid, 1, &extra_diag, &col);
        }
      }

      if(!Teuchos::nonnull(x_dot_in))
        W_out_crs->fillComplete();

      if(adjointModel) {
        Tpetra::RowMatrixTransposer<double,LO,GO> transposer(W_out_crs);
        *W_out_crs = *transposer.createTranspose();
      }
    }
    
    auto hess_g_pp = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_g_pp,0,0,0) ? outArgs.get_hess_g_pp(0,0,0) : Teuchos::null;
    const Teuchos::RCP<MatrixBased_LOWS> H_pp_out =
      Teuchos::nonnull(hess_g_pp) ?
        Teuchos::rcp_dynamic_cast<MatrixBased_LOWS>(outArgs.get_hess_g_pp(0,0,0)):
        Teuchos::null;
    
    if (Teuchos::nonnull(H_pp_out)) {
      Teuchos::RCP<Tpetra_CrsMatrix> H_pp_out_crs =
        Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(ConverterT::getTpetraOperator(H_pp_out->getMatrix()), true);
      H_pp_out_crs->resumeFill();
      H_pp_out_crs->setAllToScalar(2.0);
      
      if ((comm->getRank() == 0) && (num_p > 1)) {
        std::vector<double> vals = {2, 1};
        std::vector<typename Tpetra_CrsGraph::global_ordinal_type> indices = {0, 1};
        H_pp_out_crs->replaceGlobalValues(0, 2, &vals[0], &indices[0]);
        vals[0] = 1;
        H_pp_out_crs->replaceGlobalValues(1, 2, &vals[0], &indices[0]);
      }
      H_pp_out_crs->fillComplete();

      if(probList_->sublist("Hessian").sublist("Response 0").sublist("Parameter 0").isSublist("H_pp Solver")) {
        auto pl = probList_->sublist("Hessian").sublist("Response 0").sublist("Parameter 0").sublist("H_pp Solver");
        H_pp_out->initializeSolver(Teuchos::rcpFromRef(pl));
      }
    }
    
    if (Teuchos::nonnull(dfdp_out)) {
      dfdp_out->putScalar(0.0);
      auto dfdp_out_view = dfdp_out->getLocalViewHost(Tpetra::Access::ReadWrite);
      for (int i=0; i<myVecLength; i++) {
        const int gid = x_in->getMap()->getGlobalElement(i);
        if  (gid==0) {
          dfdp_out_view(i,0) = -1.0;
        }
        else {
          if(num_p > 1) {
            dfdp_out_view(i,0) = 0.5 - x(i) / (3.0 * std::cbrt(p(0)*p(0)));
            dfdp_out_view(i,1) = -1.0;
          } else {
            dfdp_out_view(i,0) = -0.5 - x(i) / (3.0 * std::cbrt(p(0)*p(0)));
          }
        }
      }
    }
  }
  // Response: g = 0.5*(Sum(x)-Sum(p)-12)^2 + 0.5*(p0-1)^2
  // min g(x(p), p) s.t. f(x, p) = 0 reached for p_0 = 1, p_1 = 3

  double term1, term2;
  term1 = x_in->meanValue()*vecLength - p_vec->meanValue()*num_p -12;

  {
    auto p_host = p_vec->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto p = Kokkos::subview(p_host,Kokkos::ALL(),0);
    term2 = p(0) - 1.0;
  }

  if (Teuchos::nonnull(g_out)) {
    g_out->putScalar(0.5*term1*term1 + 0.5*term2*term2);
  }

  if (dgdx_out != Teuchos::null) {
    dgdx_out->putScalar(term1);
  }
  if (dgdp_out != Teuchos::null) {
    auto dgdp_out_view = dgdp_out->getLocalViewHost(Tpetra::Access::OverwriteAll);
    dgdp_out_view(0,0) = -term1 + term2;
    if(num_p  > 1)
      dgdp_out_view(1,0) = -term1;
  }

  if (Teuchos::nonnull(f_hess_xx_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(x_direction));
    TEUCHOS_TEST_FOR_EXCEPTION(adjointModel, std::logic_error,
                     std::endl << "Error!  MockModelEval_A_Tpetra::evalModelImpl " <<
                     " adjoint Hessian not implemented." << std::endl);

    double x_direction_0;
    for (int i=0; i<myVecLength; i++) {
      if( x_in->getMap()->getGlobalElement(i) == 0) {
        auto x_direction_0_view = x_direction->getLocalViewHost(Tpetra::Access::ReadOnly);
        x_direction_0 = x_direction_0_view(i,0);
        TEUCHOS_ASSERT(comm->getRank() == 0);
      }
    }

    double temp= lag_multiplier_f_in->dot(*x_direction->getVector(0));
    comm->broadcast(0, sizeof(double), (char*)(&x_direction_0));

    auto x_host = x_in->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto x = Kokkos::subview(x_host,Kokkos::ALL(),0);
    auto f_hess_xx_v_out_view = f_hess_xx_v_out->getLocalViewHost(Tpetra::Access::OverwriteAll);
    auto lag_multiplier_f_in_view = lag_multiplier_f_in->getLocalViewHost(Tpetra::Access::ReadOnly);
    for (int i=0; i<myVecLength; i++) {
      if (x_in->getMap()->getGlobalElement(i)==0){
        f_hess_xx_v_out_view(i,0) =
            (6.0* x(i) - 1.0) * lag_multiplier_f_in_view(i,0) * x_direction->getLocalViewHost(Tpetra::Access::ReadOnly)(i,0)  + temp;
      }
      else {
        f_hess_xx_v_out_view(i,0) = x_direction_0 * lag_multiplier_f_in_view(i,0);
      }
    }
  }

  if (Teuchos::nonnull(f_hess_xp_v_out)) {
    TEUCHOS_TEST_FOR_EXCEPTION(adjointModel, std::logic_error,
                  std::endl << "Error!  MockModelEval_A_Tpetra::evalModelImpl " <<
                  " adjoint Hessian not implemented." << std::endl);
    auto p_host = p_vec->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto p = Kokkos::subview(p_host,Kokkos::ALL(),0);
    TEUCHOS_ASSERT(Teuchos::nonnull(p_direction));
    auto p_direction_view = p_direction->getLocalViewHost(Tpetra::Access::ReadOnly);
    f_hess_xp_v_out->putScalar(0);
    auto f_hess_xp_v_out_view = f_hess_xp_v_out->getLocalViewHost(Tpetra::Access::ReadWrite);
    auto lag_multiplier_f_in_view = lag_multiplier_f_in->getLocalViewHost(Tpetra::Access::ReadOnly);
    for (int i=0; i<myVecLength; i++) {
      if (x_in->getMap()->getGlobalElement(i)!=0){
        f_hess_xp_v_out_view(i,0) = - p_direction_view(0,0)*lag_multiplier_f_in_view(i,0)/(3.0* std::cbrt(p(0)*p(0)));
      }
    }
  }

  if (Teuchos::nonnull(f_hess_px_v_out)) {
    TEUCHOS_TEST_FOR_EXCEPTION(adjointModel, std::logic_error,
                  std::endl << "Error!  MockModelEval_A_Tpetra::evalModelImpl " <<
                  " adjoint Hessian not implemented." << std::endl);
    TEUCHOS_ASSERT(Teuchos::nonnull(x_direction));
    f_hess_px_v_out->putScalar(0);
    Tpetra_Vector temp_vec(lag_multiplier_f_in->getMap());
    auto p_host = p_vec->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto p = Kokkos::subview(p_host,Kokkos::ALL(),0);
    {
      auto temp_vec_view = temp_vec.getLocalViewHost(Tpetra::Access::OverwriteAll);
      for (int i=0; i<myVecLength; i++) {
        if (x_in->getMap()->getGlobalElement(i)==0)
          temp_vec_view(i,0) = 0;
        else
          temp_vec_view(i,0) = lag_multiplier_f_in->getLocalViewHost(Tpetra::Access::ReadOnly)(i,0);
      }
    }
    double temp= temp_vec.dot(*x_direction->getVector(0));
    auto f_hess_px_v_out_view = f_hess_px_v_out->getLocalViewHost(Tpetra::Access::ReadWrite);
    f_hess_px_v_out_view(0,0) = -temp/(3.0* std::cbrt(p(0)*p(0)));
  }

  if (Teuchos::nonnull(f_hess_pp_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(p_direction));
    auto p_direction_view = p_direction->getLocalViewHost(Tpetra::Access::ReadOnly);
    f_hess_pp_v_out->putScalar(0);
    Tpetra_Vector temp_vec(lag_multiplier_f_in->getMap());
    auto p_host = p_vec->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto p = Kokkos::subview(p_host,Kokkos::ALL(),0);
    {
      auto temp_vec_view = temp_vec.getLocalViewHost(Tpetra::Access::OverwriteAll);
      for (int i=0; i<myVecLength; i++) {
        if (x_in->getMap()->getGlobalElement(i)==0)
          temp_vec_view(i,0) = 0;
        else
          temp_vec_view(i,0) = lag_multiplier_f_in->getLocalViewHost(Tpetra::Access::ReadOnly)(i,0);
      }
    }
    double temp = temp_vec.dot(*x_in);
    auto f_hess_pp_v_out_view = f_hess_pp_v_out->getLocalViewHost(Tpetra::Access::ReadWrite);
    f_hess_pp_v_out_view(0,0) = 2.0*p_direction_view(0,0)*temp/(9.0* std::cbrt(std::pow(p(0),5)));
  }

  double mult = Teuchos::nonnull(lag_multiplier_g_in) ? 
                lag_multiplier_g_in->getLocalViewHost(Tpetra::Access::ReadOnly)(0,0) : 1.0;
  if (Teuchos::nonnull(g_hess_xx_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(x_direction));
    term1 = x_direction->getVector(0)->meanValue() * vecLength;
    g_hess_xx_v_out->putScalar(mult*term1);
  }

  if (Teuchos::nonnull(g_hess_xp_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(p_direction));
    auto p_direction_view = p_direction->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto g_hess_xp_v_out_view = g_hess_xp_v_out->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (int j=0; j<myVecLength; j++){
      g_hess_xp_v_out_view(j,0) = - mult*p_direction_view(0,0);
      if(num_p > 1)
        g_hess_xp_v_out_view(j,0) -= mult*p_direction_view(1,0);
    }
  }

  if (Teuchos::nonnull(g_hess_px_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(x_direction));
    term1 = x_direction->getVector(0)->meanValue() * vecLength;
    auto g_hess_px_v_out_view = g_hess_px_v_out->getLocalViewHost(Tpetra::Access::OverwriteAll);
    g_hess_px_v_out_view(0,0) = - mult*term1;
    if(num_p > 1)
      g_hess_px_v_out_view(1,0) = - mult*term1;
  }

  if (Teuchos::nonnull(g_hess_pp_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(p_direction));
    auto p_direction_view = p_direction->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto g_hess_pp_v_out_view = g_hess_pp_v_out->getLocalViewHost(Tpetra::Access::OverwriteAll);
    if (num_p ==1) {
      g_hess_pp_v_out_view(0,0) = mult*2.0*p_direction_view(0,0);
    } else {
      g_hess_pp_v_out_view(0,0) = mult*(2.0*p_direction_view(0,0)+p_direction_view(1,0));
      g_hess_pp_v_out_view(1,0) = mult*(p_direction_view(0,0)+p_direction_view(1,0));
    }
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
      auto f_out_view = f_out->getLocalViewHost(Tpetra::Access::ReadWrite);
      auto x_dot_in_view = x_dot_in->getLocalViewHost(Tpetra::Access::ReadOnly);
      for (int i=0; i<myVecLength; i++) {
        f_out_view(i,0) = -x_dot_in_view(i,0) + f_out_view(i,0);
      }
    }

    if (W_out != Teuchos::null) {
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
MockModelEval_A_Tpetra::createInArgsImpl() const
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

