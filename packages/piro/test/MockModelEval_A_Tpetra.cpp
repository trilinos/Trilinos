// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "MockModelEval_A_Tpetra.hpp"

#include <iostream>

using Teuchos::RCP;
using Teuchos::rcp;

MockModelEval_A_Tpetra::MockModelEval_A_Tpetra(const Teuchos::RCP<const Teuchos::Comm<int> >  appComm)
{
    comm = appComm;

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
    const int numParameters= 2;
    p_map = rcp(new const Tpetra_Map(numParameters, 0, comm, Tpetra::LocallyReplicated));

    Teuchos::RCP<const Thyra::VectorSpaceBase<double>> p_space =
        Thyra::createVectorSpace<double>(p_map);


    Teuchos::RCP<Tpetra_Vector> p_init = rcp(new Tpetra_Vector(p_map));
    Teuchos::RCP<Tpetra_Vector> p_lo = rcp(new Tpetra_Vector(p_map));
    Teuchos::RCP<Tpetra_Vector> p_up = rcp(new Tpetra_Vector(p_map));
    for (int i=0; i<numParameters; i++) {
      p_init->getDataNonConst()[i]= 1.0;
      p_lo->getDataNonConst()[i]= 0.1;
      p_up->getDataNonConst()[i]= 10.0;
    }

    p_vec = rcp(new Tpetra_Vector(p_map));
    p_vec->assign(*p_init);

    //set up jacobian graph
    crs_graph = rcp(new Tpetra_CrsGraph(x_map, vecLength));
    std::vector<typename Tpetra_CrsGraph::global_ordinal_type> indices(vecLength);
    for (int i=0; i<vecLength; i++) indices[i]=i;
    const int nodeNumElements = x_map->getNodeNumElements();
    for (int i=0; i<nodeNumElements; i++)
      crs_graph->insertGlobalIndices(x_map->getGlobalElement(i), vecLength, &indices[0]);
    crs_graph->fillComplete();

    // Setup nominal values, lower and upper bounds
    nominalValues = this->createInArgsImpl();
    lowerBounds = this->createInArgsImpl();
    upperBounds = this->createInArgsImpl();

    nominalValues.set_x(Thyra::createVector(x_vec, x_space));
    nominalValues.set_x_dot(Thyra::createVector(x_dot_vec, x_space));

    nominalValues.set_p(0, Thyra::createVector(p_init, p_space));
    lowerBounds.set_p(0, Thyra::createVector(p_lo, p_space));
    upperBounds.set_p(0, Thyra::createVector(p_up, p_space));
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
                     "Error!  App::ModelEval::get_p_map() only " <<
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
                     "Error!  MockModelEval_A_Tpetra::get_g_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> g_space =
        Thyra::createVectorSpace<double>(g_map);
  return g_space;
}

RCP<const  Teuchos::Array<std::string> > MockModelEval_A_Tpetra::get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_names() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);

  Teuchos::Ordinal num_p = p_map->getNodeNumElements();
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

  result.setHessianSupports(true);

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


  const Teuchos::RCP<const Thyra::VectorBase<double>> p_in = inArgs.get_p(0);
  if (Teuchos::nonnull(p_in))
    p_vec->assign(*ConverterT::getConstTpetraVector(p_in));



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

  auto f_hess_xx_v = outArgs.get_hess_vec_prod_f_xx();
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_xx_v_out =
      Teuchos::nonnull(f_hess_xx_v) ?
        ConverterT::getTpetraMultiVector(f_hess_xx_v) :
        Teuchos::null;

  auto f_hess_xp_v = outArgs.get_hess_vec_prod_f_xp(0);
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_xp_v_out =
      Teuchos::nonnull(f_hess_xp_v) ?
        ConverterT::getTpetraMultiVector(f_hess_xp_v) :
        Teuchos::null;

  auto f_hess_px_v = outArgs.get_hess_vec_prod_f_px(0);
  const Teuchos::RCP<Tpetra_MultiVector> f_hess_px_v_out =
      Teuchos::nonnull(f_hess_px_v) ?
        ConverterT::getTpetraMultiVector(f_hess_px_v) :
        Teuchos::null;

  auto f_hess_pp_v = outArgs.get_hess_vec_prod_f_pp(0,0);
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

  auto g_hess_xx_v = outArgs.get_hess_vec_prod_g_xx(0);
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_xx_v_out =
      Teuchos::nonnull(g_hess_xx_v) ?
        ConverterT::getTpetraMultiVector(g_hess_xx_v) :
        Teuchos::null;

  auto g_hess_xp_v = outArgs.get_hess_vec_prod_g_xp(0,0);
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_xp_v_out =
      Teuchos::nonnull(g_hess_xp_v) ?
        ConverterT::getTpetraMultiVector(g_hess_xp_v) :
        Teuchos::null;

  auto g_hess_px_v = outArgs.get_hess_vec_prod_g_px(0,0);
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_px_v_out =
      Teuchos::nonnull(g_hess_px_v) ?
        ConverterT::getTpetraMultiVector(g_hess_px_v) :
        Teuchos::null;

  auto g_hess_pp_v = outArgs.get_hess_vec_prod_g_pp(0,0,0);
  const Teuchos::RCP<Tpetra_MultiVector> g_hess_pp_v_out =
      Teuchos::nonnull(g_hess_pp_v) ?
        ConverterT::getTpetraMultiVector(g_hess_pp_v) :
        Teuchos::null;

  auto x = x_in->getData();
  auto p = p_vec->getData();

  double x0;
  for (int i=0; i<myVecLength; i++) {
    if( x_in->getMap()->getGlobalElement(i) == 0) {
      x0 = x[i];
      TEUCHOS_ASSERT(comm->getRank() == 0);
    }
  }
  comm->broadcast(0, sizeof(double), (char*)(&x0));

  if (f_out != Teuchos::null) {
    f_out->putScalar(0.0);
    auto f_out_data = f_out->getDataNonConst();
    for (int i=0; i<myVecLength; i++) {
      int gid = x_in->getMap()->getGlobalElement(i);

      if (gid==0) { // f_0 = (x_0)^3 - p_0
        f_out_data[i] = x[i] * x[i] * x[i] -  p[0];
      }
      else // f_i = x_i * (1 + x_0 - p_0^(1/3)) - (i+p_1) - 0.5*(x_0 - p_0),  (for i != 0)
        f_out_data[i] = x[i] - (gid + p[1]) - 0.5*(x0 - p[0]) + x[i] * (x0 - std::cbrt(p[0]));
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
        diag = 1.0 + x0 - std::cbrt(p[0]);
        W_out_crs->replaceLocalValues(i, 1, &diag, &i);
        decltype(gid) col = 0;
        extra_diag = -0.5 + x[i];
        W_out_crs->replaceGlobalValues(gid, 1, &extra_diag, &col);
      }
    }
    if(!Teuchos::nonnull(x_dot_in))
      W_out_crs->fillComplete();
  }

  if (Teuchos::nonnull(dfdp_out)) {
    dfdp_out->putScalar(0.0);
    auto dfdp_out_data_0 = dfdp_out->getVectorNonConst(0)->getDataNonConst();
    auto dfdp_out_data_1 = dfdp_out->getVectorNonConst(1)->getDataNonConst();
    for (int i=0; i<myVecLength; i++) {
      const int gid = x_in->getMap()->getGlobalElement(i);
      if  (gid==0) {
        dfdp_out_data_0[i] = -1.0;
      }
      else {
        dfdp_out_data_0[i] = 0.5 - x[i] / (3.0 * std::cbrt(p[0]*p[0]));
        dfdp_out_data_1[i] = -1.0;
      }
    }
  }

  // Response: g = 0.5*(Sum(x)-Sum(p)-12)^2 + 0.5*(p0-1)^2
  // min g(x(p), p) s.t. f(x, p) = 0 reached for p_0 = 1, p_1 = 3

  double term1, term2;
  term1 = x_in->meanValue();
  term1 =  vecLength * term1 - (p[0] + p[1]) - 12.0;
  term2 = p[0] - 1.0;

  if (Teuchos::nonnull(g_out)) {
    g_out->getDataNonConst()[0] = 0.5*term1*term1 + 0.5*term2*term2;
  }

  if (dgdx_out != Teuchos::null) {
    dgdx_out->putScalar(term1);
    auto dgdx = dgdx_out->getVector(0)->getData();
  }
  if (dgdp_out != Teuchos::null) {
    dgdp_out->putScalar(0.0);
    dgdp_out->getVectorNonConst(0)->getDataNonConst()[0] = -term1 + term2;
    dgdp_out->getVectorNonConst(0)->getDataNonConst()[1] = -term1;
  }

  if (Teuchos::nonnull(f_hess_xx_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(x_direction));

    double x_direction_0;
    for (int i=0; i<myVecLength; i++) {
      if( x_in->getMap()->getGlobalElement(i) == 0) {
        x_direction_0 = x_direction->getVector(0)->getData()[i];
        TEUCHOS_ASSERT(comm->getRank() == 0);
      }
    }

    double temp= lag_multiplier_f_in->dot(*x_direction->getVector(0));
    comm->broadcast(0, sizeof(double), (char*)(&x_direction_0));

    f_hess_xx_v_out->getVectorNonConst(0)->putScalar(0);
    for (int i=0; i<myVecLength; i++) {
      if (x_in->getMap()->getGlobalElement(i)==0){
        f_hess_xx_v_out->getVectorNonConst(0)->getDataNonConst()[i] =
            (6.0* x[i] - 1.0) *lag_multiplier_f_in->getData()[i] * x_direction->getVector(0)->getData()[i]  + temp;
      }
      else {
        f_hess_xx_v_out->getVectorNonConst(0)->getDataNonConst()[i] = x_direction_0 * lag_multiplier_f_in->getData()[i];
      }
    }
  }

  if (Teuchos::nonnull(f_hess_xp_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(p_direction));
    f_hess_xp_v_out->getVectorNonConst(0)->putScalar(0);
    for (int i=0; i<myVecLength; i++) {
      if (x_in->getMap()->getGlobalElement(i)!=0){
        f_hess_xp_v_out->getVectorNonConst(0)->getDataNonConst()[i] =
            - p_direction->getVector(0)->getData()[0]*lag_multiplier_f_in->getData()[i]/(3.0* std::cbrt(p[0]*p[0]));
      }
    }
  }

  if (Teuchos::nonnull(f_hess_px_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(x_direction));
    f_hess_px_v_out->getVectorNonConst(0)->putScalar(0);
    Tpetra_Vector temp_vec(lag_multiplier_f_in->getMap());
    for (int i=0; i<myVecLength; i++) {
      if (x_in->getMap()->getGlobalElement(i)==0)
        temp_vec.getDataNonConst()[i] = 0;
      else
        temp_vec.getDataNonConst()[i] = lag_multiplier_f_in->getData()[i];
    }
    double temp= temp_vec.dot(*x_direction->getVector(0));
    f_hess_px_v_out->getVectorNonConst(0)->getDataNonConst()[0] = -temp/(3.0* std::cbrt(p[0]*p[0]));
  }

  if (Teuchos::nonnull(f_hess_pp_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(p_direction));
    f_hess_pp_v_out->getVectorNonConst(0)->putScalar(0);
    Tpetra_Vector temp_vec(lag_multiplier_f_in->getMap());
    for (int i=0; i<myVecLength; i++) {
      if (x_in->getMap()->getGlobalElement(i)==0)
        temp_vec.getDataNonConst()[i] = 0;
      else
        temp_vec.getDataNonConst()[i] = lag_multiplier_f_in->getData()[i];
    }
    double temp = temp_vec.dot(*x_in);
    f_hess_pp_v_out->getVectorNonConst(0)->getDataNonConst()[0] = 2.0*p_direction->getVector(0)->getData()[0]*temp/(9.0* std::cbrt(std::pow(p[0],5)));
  }

  double mult = Teuchos::nonnull(lag_multiplier_g_in) ? lag_multiplier_g_in->getData()[0] : 1.0;
  if (Teuchos::nonnull(g_hess_xx_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(x_direction));
    term1 = x_direction->getVector(0)->meanValue() * vecLength;
    for (int j=0; j<myVecLength; j++)
      g_hess_xx_v_out->getVectorNonConst(0)->getDataNonConst()[j] = mult*term1;
  }

  if (Teuchos::nonnull(g_hess_xp_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(p_direction));
    const auto direction_p = p_direction->getVector(0)->getData();
    for (int j=0; j<myVecLength; j++)
      g_hess_xp_v_out->getVectorNonConst(0)->getDataNonConst()[j] = - mult*(direction_p[0]+direction_p[1]);
  }

  if (Teuchos::nonnull(g_hess_px_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(x_direction));
    term1 = x_direction->getVector(0)->meanValue() * vecLength;
    g_hess_px_v_out->getVectorNonConst(0)->getDataNonConst()[0] = - mult*term1;
    g_hess_px_v_out->getVectorNonConst(0)->getDataNonConst()[1] = - mult*term1;
  }

  if (Teuchos::nonnull(g_hess_pp_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(p_direction));
    const auto direction_p = p_direction->getVector(0)->getData();
    g_hess_pp_v_out->getVectorNonConst(0)->getDataNonConst()[0] = mult*(2.0*direction_p[0]+direction_p[1]);
    g_hess_pp_v_out->getVectorNonConst(0)->getDataNonConst()[1] = mult*(direction_p[0]+direction_p[1]);
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

