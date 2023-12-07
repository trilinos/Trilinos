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

#include "MockModelEval_H_Tpetra.hpp"
#include "Thyra_LinearOpWithSolveBase_decl.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "MatrixMarket_Tpetra.hpp"


#include <iostream>

using Teuchos::RCP;
using Teuchos::rcp;


DfDpOp::DfDpOp(
      const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> x_space,
      const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> p_space) :
      x_space_(x_space),
      p_space_(p_space) {};

    //! Overrides Thyra::LinearOpBase purely virtual method
void 
DfDpOp::applyImpl (const Thyra::EOpTransp M_trans,
                    const Thyra::MultiVectorBase<double>& X,
                    const Teuchos::Ptr<Thyra::MultiVectorBase<double>>& Y,
                    const double /* alpha */,
                    const double /* beta */) const {
      (void) M_trans;
      Teuchos::RCP<const Tpetra_MultiVector> x_vec = ConverterT::getConstTpetraMultiVector(Teuchos::rcpFromRef(X));
      Teuchos::RCP<Tpetra_MultiVector> y_vec = ConverterT::getTpetraMultiVector(Teuchos::rcpFromPtr(Y));

     
      const int nodeNumElements = x_vec->getMap()->getLocalNumElements();
      y_vec->putScalar(0.0);
      for (int i=0; i<nodeNumElements; i++) {
        y_vec->getDataNonConst(0)[i] = -3*std::pow(p_vec_->getData()[i],2)*x_vec->getData(0)[i];
      }
    }

MockModelEval_H_Tpetra::MockModelEval_H_Tpetra(const Teuchos::RCP<const Teuchos::Comm<int> >  appComm, bool /*adjoint*/, const Teuchos::RCP<Teuchos::ParameterList>& problemList, bool hessianSupport) //problem is self-adjoint
 {
    comm = appComm;
    hessSupport = hessianSupport;

    //set up map and initial guess for solution vector
    const int vecLength = 51;
    h = 1.0/(vecLength-1);
    x_map = rcp(new Tpetra_Map(vecLength, 0, comm));
    x_vec = rcp(new Tpetra_Vector(x_map));
    x_dot_vec = rcp(new Tpetra_Vector(x_map));

    coords_vec = rcp(new Tpetra_Vector(x_map));
    for (int i=0; i< x_map->getLocalNumElements (); ++i)
      coords_vec->getDataNonConst()[i]=x_map->getGlobalElement(i)*h;

    x_space = Thyra::createVectorSpace<double>(x_map);

    //set up responses
    const int numResponses = 1;
    g_map = rcp(new const Tpetra_Map(numResponses , 0, comm, Tpetra::LocallyReplicated));

    //set up parameters
    const int numParameters = 2;
    p_map = rcp(new Tpetra_Map(vecLength, 0, comm));

    p_space = Thyra::createVectorSpace<double>(p_map);

    Teuchos::RCP<Tpetra_Vector> p_init = rcp(new Tpetra_Vector(p_map));
    Teuchos::RCP<Tpetra_Vector> p1_init = rcp(new Tpetra_Vector(p_map));
    Teuchos::RCP<Tpetra_Vector> p_lo = rcp(new Tpetra_Vector(p_map));
    Teuchos::RCP<Tpetra_Vector> p_up = rcp(new Tpetra_Vector(p_map));

    int nodeNumElements = p_map->getLocalNumElements();
    for (int i=0; i<nodeNumElements; i++) {
      p_init->getDataNonConst()[i] = 1+coords_vec->getData()[i];
    }

    nodeNumElements = x_map->getLocalNumElements();
    for (int i=0; i<nodeNumElements; i++) {
      p1_init->getDataNonConst()[i] = std::pow(1+coords_vec->getData()[i],3);
    }

    p_lo->putScalar(0.1);
    p_up->putScalar(10.0);

    p_vec_0 = rcp(new Tpetra_Vector(p_map));
    p_vec_0->assign(*p_init);

    p_vec_1 = rcp(new Tpetra_Vector(x_map)); //this is morally the state x
    p_vec_1->assign(*p1_init);


    //set up jacobian graph
    crs_graph = rcp(new Tpetra_CrsGraph(x_map, 3));
    {
      std::vector<typename Tpetra_CrsGraph::global_ordinal_type> indices(1);
      const int nodeNumElements = x_map->getLocalNumElements();
      for (int i=0; i<nodeNumElements; i++) {
        auto gid = x_map->getGlobalElement(i);
        indices[0] = gid;      
        crs_graph->insertGlobalIndices(gid, 1, &indices[0]);
      }
    }
    crs_graph->fillComplete();

    //set up hessian graph
    hess_crs_graph_p = rcp(new Tpetra_CrsGraph(p_map, 3));
    {
      std::vector<typename Tpetra_CrsGraph::global_ordinal_type> indices(3);
      const int nodeNumElements = p_map->getLocalNumElements();
      for (int i=0; i<nodeNumElements; i++) {
        auto gid = p_map->getGlobalElement(i);
        indices.clear();
        if(gid > 0)
          indices.push_back(gid-1);
        indices.push_back(gid);
        if(gid < vecLength-1)
          indices.push_back(gid+1);        
        hess_crs_graph_p->insertGlobalIndices(gid, indices.size(), &indices[0]);
      }
    }
    hess_crs_graph_p->fillComplete();

    hess_crs_graph_x = rcp(new Tpetra_CrsGraph(x_map, 3));
    {
      std::vector<typename Tpetra_CrsGraph::global_ordinal_type> indices(3);
      const int nodeNumElements = x_map->getLocalNumElements();
      for (int i=0; i<nodeNumElements; i++) {
        auto gid = x_map->getGlobalElement(i);
        indices.clear();
        if(gid > 0)
          indices.push_back(gid-1);
        indices.push_back(gid);
        if(gid < vecLength-1)
          indices.push_back(gid+1);        
        hess_crs_graph_x->insertGlobalIndices(gid, indices.size(), &indices[0]);
      }
    }
    hess_crs_graph_x->fillComplete();


    // Setup nominal values, lower and upper bounds
    nominalValues = this->createInArgsImpl();
    lowerBounds = this->createInArgsImpl();
    upperBounds = this->createInArgsImpl();

    nominalValues.set_x(Thyra::createVector(x_vec, x_space));
    nominalValues.set_x_dot(Thyra::createVector(x_dot_vec, x_space));

    nominalValues.set_p(0, Thyra::createVector(p_init, p_space));
    nominalValues.set_p(1, Thyra::createVector(p1_init, p_space));
    lowerBounds.set_p(0, Thyra::createVector(p_lo, p_space));
    upperBounds.set_p(0, Thyra::createVector(p_up, p_space));

    nominalValues.set_p(1, Thyra::createVector(x_vec, x_space));

    probList_ = problemList;
}

MockModelEval_H_Tpetra::~MockModelEval_H_Tpetra()
{
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_H_Tpetra::get_x_space() const
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> x_space =
      Thyra::createVectorSpace<double>(x_map);
  return x_space;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_H_Tpetra::get_f_space() const
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> f_space =
      Thyra::createVectorSpace<double>(x_map);
  return f_space;
}

Teuchos::RCP<Thyra::VectorBase<double>>
MockModelEval_H_Tpetra::get_p_opt(int j) const {
    if(j ==0) {
      Teuchos::RCP<Tpetra_Vector> p0_opt = rcp(new Tpetra_Vector(p_map));
      const int nodeNumElements = p_map->getLocalNumElements();
      for (int i=0; i<nodeNumElements; i++) {
        p0_opt->getDataNonConst()[i] = 1+coords_vec->getData()[i];
      }
      auto p0_out = Thyra::createVector(p0_opt, p_space);
      return p0_out;
    } else {
      Teuchos::RCP<Tpetra_Vector> p1_opt = rcp(new Tpetra_Vector(x_map));
      const int nodeNumElements = x_map->getLocalNumElements();
      for (int i=0; i<nodeNumElements; i++) {
        p1_opt->getDataNonConst()[i] = std::pow(1+coords_vec->getData()[i],3);
      }
      auto p1_out = Thyra::createVector(p1_opt, p_space);
      return p1_out;
    }
}

Teuchos::RCP<Thyra::VectorBase<double>>
MockModelEval_H_Tpetra::get_true_p_opt(int j) const {
    TEUCHOS_ASSERT(j==0);
    Teuchos::RCP<Tpetra_Vector> true_p0_opt  = rcp(new Tpetra_Vector(p_map));
    Teuchos::RCP<Tpetra_Vector> T = rcp(new Tpetra_Vector(p_map));
    Teuchos::RCP<Tpetra_Vector> true_f = rcp(new Tpetra_Vector(x_map));
    const int nodeNumElements = p_map->getLocalNumElements();
    for (int i=0; i<nodeNumElements; i++) {
      auto& p0 = true_p0_opt->getDataNonConst()[i];
      p0 = 1+coords_vec->getData()[i];
      T->getDataNonConst()[i] = std::pow(1+coords_vec->getData()[i],3);
      true_f->getDataNonConst()[i] =  p0*p0*(p0 + 0.2) - T->getData()[i];
    }
    // Newton method
    while (true_f->normInf() > 1e-12) {
      for (int i=0; i<nodeNumElements; i++) {
        auto& p0 = true_p0_opt->getDataNonConst()[i];
        auto& f = true_f->getDataNonConst()[i];
        auto df = p0*(3*p0+0.4);
        p0 -= f/df;
        f =  p0*p0*(p0 + 0.2) - T->getData()[i];
      }
    }
    auto true_p0_out = Thyra::createVector(true_p0_opt, p_space);
    return true_p0_out;
}

Teuchos::RCP<Thyra::VectorBase<double>>
MockModelEval_H_Tpetra::get_param_samples(int k) const {   
    Teuchos::RCP<Tpetra_Vector> sample = rcp(new Tpetra_Vector(p_map));
    const int nodeNumElements = p_map->getLocalNumElements();
    if(k ==0) {
      for (int i=0; i<nodeNumElements; i++)
        sample->getDataNonConst()[i] = 1+coords_vec->getData()[i];
    } else {
      for (int i=0; i<nodeNumElements; i++)
        sample->getDataNonConst()[i] = (1+coords_vec->getData()[i])*coords_vec->getData()[i];
    }

    // A Thyra::ProductVector is expected.
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> spaces(1);
    spaces[0] = p_space;
    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<double>> prd_space = Teuchos::rcp(new Thyra::DefaultProductVectorSpace<double>(spaces));
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<double>>> vecs(1);
    vecs[0] = Thyra::createVector(sample, p_space);
    Teuchos::RCP<Thyra::DefaultProductVector<double>> prd_sample = Thyra::defaultProductVector<double>(prd_space, vecs());
    return prd_sample;
}

Teuchos::RCP<Thyra::VectorBase<double>>
MockModelEval_H_Tpetra::get_solution_diff_at_samples(int k) const {
    Teuchos::RCP<Tpetra_Vector> diff = rcp(new Tpetra_Vector(x_map));
    const int nodeNumElements = x_map->getLocalNumElements();
    if(k ==0) {
      for (int i=0; i<nodeNumElements; i++)
        diff->getDataNonConst()[i] = 0.2*std::pow(1+coords_vec->getData()[i],2);
    } else {
      for (int i=0; i<nodeNumElements; i++)
        diff->getDataNonConst()[i] = 0.2*std::pow((1+coords_vec->getData()[i])*coords_vec->getData()[i],2);
    }
    auto diff_out = Thyra::createVector(diff, p_space);
    return diff_out;
}


Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_H_Tpetra::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l > 1, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 2 parameter vectors.  Supplied index l = " <<
                     l << std::endl);
  if (l==0)
    return p_space;
  else 
    return x_space;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MockModelEval_H_Tpetra::get_g_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l > 5, std::logic_error,
                     std::endl <<
                     "Error!  MockModelEval_H_Tpetra::get_g_map() only " <<
                     " supports 6 responses.  Supplied index l = " <<
                     l << std::endl);
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> g_space =
        Thyra::createVectorSpace<double>(g_map);
  return g_space;
}

RCP<const  Teuchos::Array<std::string> > MockModelEval_H_Tpetra::get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l > 1, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_names() only " <<
                     " supports 2 parameter vectors.  Supplied index l = " <<
                     l << std::endl);

  RCP<Teuchos::Array<std::string> > p_names =
      rcp(new Teuchos::Array<std::string>(1) );

  std::stringstream ss;
  ss << "Parameter " << l;
  const std::string name = ss.str();
  (*p_names)[0] = name;

  return p_names;
}


Teuchos::RCP<Thyra::LinearOpBase<double>>
MockModelEval_H_Tpetra::create_W_op() const
{
  const Teuchos::RCP<Tpetra_Operator> W =
      Teuchos::rcp(new Tpetra_CrsMatrix(crs_graph));
  return Thyra::createLinearOp(W);
}

//! Create preconditioner operator
Teuchos::RCP<Thyra::PreconditionerBase<double>>
MockModelEval_H_Tpetra::create_W_prec() const
{
  return Teuchos::null;
}

Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double>>
MockModelEval_H_Tpetra::get_W_factory() const
{
  return Teuchos::null;
}

Teuchos::RCP<Thyra::LinearOpBase<double>>
MockModelEval_H_Tpetra::create_hess_g_pp( int j, int l1, int l2 ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!((l1 == l2) && (j > 0) && (j < 6) && (l1 < 2)), std::logic_error,
                    std::endl <<
                    "Error!  MockModelEval_H_Tpetra::create_hess_g_pp() " <<
                    " supports only indices (1,0,0), (2,0,0), (3,0,0), (4,1,1) and (5,1,1). Supplied index l1 = " <<
                    l1 << ", l2 = " << l2 << std::endl);
  Teuchos::RCP<Tpetra_Operator> H;

  if(l1 == 0)
    H = Teuchos::rcp(new Tpetra_CrsMatrix(hess_crs_graph_p));
  else
    H = Teuchos::rcp(new Tpetra_CrsMatrix(hess_crs_graph_x));
  std::cout << __FILE__ << " " << __LINE__ << " " << j << " "<< l1 << " " << l2 << std::endl;

  return Teuchos::rcp(new MatrixBased_LOWS(Thyra::createLinearOp(H)));
}

Teuchos::RCP<Thyra::LinearOpBase<double>>
MockModelEval_H_Tpetra::create_DfDp_op_impl(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      j != 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl
          << "Error! MockModelEval_H_Tpetra::create_DfDp_op_impl():  "
          << "Supports only index j==0. Function called with index j = "
          << j
          << std::endl);

  return Teuchos::rcp( new DfDpOp(x_space, p_space) );
}


Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_H_Tpetra::getNominalValues() const
{
  return nominalValues;
}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_H_Tpetra::getLowerBounds() const
{
  return lowerBounds;
}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_H_Tpetra::getUpperBounds() const
{
  return upperBounds;
}


Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_H_Tpetra::createInArgs() const
{
  return this->createInArgsImpl();
}

void
MockModelEval_H_Tpetra::reportFinalPoint(
    const Thyra::ModelEvaluatorBase::InArgs<double>& /* finalPoint */,
    const bool /* wasSolved */) {
  // Do nothing  
}

Thyra::ModelEvaluatorBase::OutArgs<double>
MockModelEval_H_Tpetra::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<double> result;
  result.setModelEvalDescription(this->description());
  result.set_Np_Ng(2, 6);

  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f, true);
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op, true);
  result.set_W_properties(Thyra::ModelEvaluatorBase::DerivativeProperties(
      Thyra::ModelEvaluatorBase::DERIV_LINEARITY_UNKNOWN,
      Thyra::ModelEvaluatorBase::DERIV_RANK_FULL,
      true));

  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, 0, Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, 0, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
  result.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, 0, 0, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);

  result.setSupports(
        Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, 0, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);

  result.setHessianSupports(hessSupport);

  return result;
}

void MockModelEval_H_Tpetra::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<double>&  inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<double>& outArgs) const
{

  // Parse InArgs

  const Teuchos::RCP<const Tpetra_Vector> x_in =
      ConverterT::getConstTpetraVector(inArgs.get_x());
  if (!Teuchos::nonnull(x_in)) std::cerr << "ERROR: MockModelEval_H_Tpetra requires x as inargs\n";

  const Teuchos::RCP<const Tpetra_Vector> x_dot_in =
      Teuchos::nonnull(inArgs.get_x_dot()) ?
          ConverterT::getConstTpetraVector(inArgs.get_x_dot()) :
          Teuchos::null;

  const Teuchos::RCP<const Thyra::VectorBase<double>> p0_in = inArgs.get_p(0);
  if (Teuchos::nonnull(p0_in))
    p_vec_0->assign(*ConverterT::getConstTpetraVector(p0_in));

  const Teuchos::RCP<const Thyra::VectorBase<double>> p1_in = inArgs.get_p(1);
  if (Teuchos::nonnull(p1_in))
    p_vec_1->assign(*ConverterT::getConstTpetraVector(p1_in));

  int myVecLength = x_in->getLocalLength();

  // Parse OutArgs

  const Teuchos::RCP<Tpetra_Vector> f_out =
      Teuchos::nonnull(outArgs.get_f()) ?
          ConverterT::getTpetraVector(outArgs.get_f()) :
          Teuchos::null;

  const Teuchos::RCP<Thyra::VectorBase<double>> g_base0 = outArgs.get_g(0);
  Teuchos::RCP<Tpetra_Vector> g0_out = Teuchos::nonnull(g_base0) ?
      ConverterT::getTpetraVector(g_base0) :
      Teuchos::null;

  const Teuchos::RCP<Thyra::VectorBase<double>> g_base1 = outArgs.get_g(1);
  Teuchos::RCP<Tpetra_Vector> g1_out = Teuchos::nonnull(g_base1) ?
      ConverterT::getTpetraVector(g_base1) :
      Teuchos::null;

  const Teuchos::RCP<Tpetra_Operator> W_out =
      Teuchos::nonnull(outArgs.get_W_op()) ?
          ConverterT::getTpetraOperator(outArgs.get_W_op()) :
          Teuchos::null;

  const Teuchos::RCP<Thyra::LinearOpBase<double>> dfdp_0_out =
      outArgs.get_DfDp(0).getLinearOp();

  const Teuchos::RCP<Thyra::MultiVectorBase<double>> dg0dp_base =
      outArgs.get_DgDp(0, 0).getMultiVector();
  const Teuchos::RCP<Tpetra_MultiVector> dg0dp_out =
      Teuchos::nonnull(dg0dp_base) ?
          ConverterT::getTpetraMultiVector(dg0dp_base) :
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
  auto coords = coords_vec->getData();

  if (f_out != Teuchos::null) {
    f_out->putScalar(0.0);
    auto f_out_data = f_out->getDataNonConst();
    for (int i=0; i<myVecLength; i++)
      f_out_data[i] = x[i]-std::pow(p_0[i],3);
  }
  if (W_out != Teuchos::null) {
    Teuchos::RCP<Tpetra_CrsMatrix> W_out_crs =
      Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(W_out, true);
    W_out_crs->resumeFill();
    W_out_crs->setAllToScalar(0.0);

    double diag=1.0;
    for (int i=0; i<myVecLength; i++)
      W_out_crs->replaceLocalValues(i, 1, &diag, &i);
      W_out_crs->fillComplete();
  }

  {
    auto hess_g_pp = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_g_pp,1,0,0) ? outArgs.get_hess_g_pp(1,0,0) : Teuchos::null; 
    const Teuchos::RCP<MatrixBased_LOWS> H_pp_out =
      Teuchos::nonnull(hess_g_pp) ?
        Teuchos::rcp_dynamic_cast<MatrixBased_LOWS>(outArgs.get_hess_g_pp(1,0,0)):
        Teuchos::null;

    if (Teuchos::nonnull(H_pp_out)) {
      Teuchos::RCP<Tpetra_CrsMatrix> H_pp_out_crs =
        Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(ConverterT::getTpetraOperator(H_pp_out->getMatrix()), true);
      H_pp_out_crs->resumeFill();
      H_pp_out_crs->setAllToScalar(0.0);

      {
        Teuchos::Array<typename Tpetra_CrsGraph::global_ordinal_type> indices(3);
        Teuchos::Array<double> vals(3);
        const int nodeNumElements = p_map->getLocalNumElements();
        for (int i=0; i<nodeNumElements; i++) {
          auto gid = p_map->getGlobalElement(i);
          indices.clear();
          vals.clear();

          double s=0.01, m=1.0;
          
          if(gid > 0) {
            indices.push_back(gid-1);
            vals.push_back(-s/h+m*h/6);
          }

          indices.push_back(gid);
          if(gid == 0 || gid == x_map->getGlobalNumElements()-1)
            vals.push_back(s/h+m*h/3);
          else
            vals.push_back(2*s/h+2*m*h/3);

          if(gid < x_map->getGlobalNumElements()-1) {
            indices.push_back(gid+1);
            vals.push_back(-s/h+m*h/6);
          }
          H_pp_out_crs->replaceGlobalValues(gid, indices, vals);
        }
      }
      H_pp_out_crs->fillComplete();

      Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> >   matWriter;
      matWriter.writeSparseFile("H_100", H_pp_out_crs);

      if(probList_->sublist("Hessian").sublist("Response 1").sublist("Parameter 0").isSublist("H_pp Solver")) {
        auto pl = probList_->sublist("Hessian").sublist("Response 1").sublist("Parameter 0").sublist("H_pp Solver");
        H_pp_out->initializeSolver(Teuchos::rcpFromRef(pl));
      }
    }
  }

    {
    auto hess_g_pp = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_g_pp,2,0,0) ? outArgs.get_hess_g_pp(2,0,0) : Teuchos::null; 
    const Teuchos::RCP<MatrixBased_LOWS> H_pp_out =
      Teuchos::nonnull(hess_g_pp) ?
        Teuchos::rcp_dynamic_cast<MatrixBased_LOWS>(outArgs.get_hess_g_pp(2,0,0)):
        Teuchos::null;

    if (Teuchos::nonnull(H_pp_out)) {
      Teuchos::RCP<Tpetra_CrsMatrix> H_pp_out_crs =
        Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(ConverterT::getTpetraOperator(H_pp_out->getMatrix()), true);
      H_pp_out_crs->resumeFill();
      H_pp_out_crs->setAllToScalar(0.0);

      {
        Teuchos::Array<typename Tpetra_CrsGraph::global_ordinal_type> indices(3);
        Teuchos::Array<double> vals(3);
        const int nodeNumElements = p_map->getLocalNumElements();
        for (int i=0; i<nodeNumElements; i++) {
          auto gid = p_map->getGlobalElement(i);
          indices.clear();
          vals.clear();
          if(gid > 0) {
            indices.push_back(gid-1);
            vals.push_back(h/6);
          }

          indices.push_back(gid);
          if(gid == 0 || gid == x_map->getGlobalNumElements()-1)
            vals.push_back(h/3);
          else
            vals.push_back(2*h/3);

          if(gid < x_map->getGlobalNumElements()-1) {
            indices.push_back(gid+1);
            vals.push_back(h/6);
          }
          H_pp_out_crs->replaceGlobalValues(gid, indices, vals);
        }
      }
      H_pp_out_crs->fillComplete();

      Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> >   matWriter;
      matWriter.writeSparseFile("H_200", H_pp_out_crs);

      if(probList_->sublist("Hessian").sublist("Response 2").sublist("Parameter 0").isSublist("H_pp Solver")) {
        auto pl = probList_->sublist("Hessian").sublist("Response 2").sublist("Parameter 0").sublist("H_pp Solver");
        H_pp_out->initializeSolver(Teuchos::rcpFromRef(pl));
      }
    }
  }

  {
    auto hess_g_pp = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_g_pp,3,0,0) ? outArgs.get_hess_g_pp(3,0,0) : Teuchos::null; 
    const Teuchos::RCP<MatrixBased_LOWS> H_pp_out =
      Teuchos::nonnull(hess_g_pp) ?
        Teuchos::rcp_dynamic_cast<MatrixBased_LOWS>(outArgs.get_hess_g_pp(3,0,0)):
        Teuchos::null;

    if (Teuchos::nonnull(H_pp_out)) {
      Teuchos::RCP<Tpetra_CrsMatrix> H_pp_out_crs =
        Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(ConverterT::getTpetraOperator(H_pp_out->getMatrix()), true);
      H_pp_out_crs->resumeFill();
      H_pp_out_crs->setAllToScalar(0.0);
      
      //computing Cholesky coeff of mass matrix. Not scalable
      int M = x_map->getGlobalNumElements()-1;
      std::vector<double> a(M+1), b(M);
      a[0] = sqrt(h/3);
      for (size_t i=0; i< M-1; ++i) {
        b[i] = h/6.0/a[i];
        a[i+1] = sqrt(2.0*h/3 - b[i]*b[i]);
      }
      b[M-1] = h/6.0/a[M-1];
      a[M] = sqrt(1.0*h/3 - b[M-1]*b[M-1]);

      {
        Teuchos::Array<typename Tpetra_CrsGraph::global_ordinal_type> indices(3);
        Teuchos::Array<double> vals(3);
        const int nodeNumElements = p_map->getLocalNumElements();
        for (int i=0; i<nodeNumElements; i++) {
          auto gid = p_map->getGlobalElement(i);
          indices.clear();
          vals.clear();
          if(gid > 0) {
            indices.push_back(gid-1);
            vals.push_back(b[gid-1]);
          }

          indices.push_back(gid);
          vals.push_back(a[gid]);
       
          H_pp_out_crs->replaceGlobalValues(gid, indices, vals);
        }
      }
      H_pp_out_crs->fillComplete();

      Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> >   matWriter;
      matWriter.writeSparseFile("H_300", H_pp_out_crs);

      if(probList_->sublist("Hessian").sublist("Response 3").sublist("Parameter 0").isSublist("H_pp Solver")) {
        auto pl = probList_->sublist("Hessian").sublist("Response 3").sublist("Parameter 0").sublist("H_pp Solver");
        H_pp_out->initializeSolver(Teuchos::rcpFromRef(pl));
      }
    }
  }

  {
    int p_index = 1;
    auto hess_g_pp = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_g_pp,4,p_index,p_index) ? outArgs.get_hess_g_pp(4,p_index,p_index) : Teuchos::null; 
    const Teuchos::RCP<MatrixBased_LOWS> H_pp_out =
      Teuchos::nonnull(hess_g_pp) ?
        Teuchos::rcp_dynamic_cast<MatrixBased_LOWS>(outArgs.get_hess_g_pp(4,p_index,p_index)):
        Teuchos::null;

    if (Teuchos::nonnull(H_pp_out)) {
      Teuchos::RCP<Tpetra_CrsMatrix> H_pp_out_crs =
        Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(ConverterT::getTpetraOperator(H_pp_out->getMatrix()), true);
      H_pp_out_crs->resumeFill();
      H_pp_out_crs->setAllToScalar(0.0);

      {
        Teuchos::Array<typename Tpetra_CrsGraph::global_ordinal_type> indices(3);
        const int nodeNumElements = x_map->getLocalNumElements();
         Teuchos::Array<double> vals(3);
        
        for (int i=0; i<nodeNumElements; i++) {

          auto gid = x_map->getGlobalElement(i);
          indices.clear();
          vals.clear();

          double s=0.05, m=1.0;
          
          if(gid > 0) {
            indices.push_back(gid-1);
            vals.push_back(-s/h+m*h/6);
          }

          indices.push_back(gid);
          if(gid == 0 || gid == x_map->getGlobalNumElements()-1)
            vals.push_back(s/h+m*h/3);
          else
            vals.push_back(2*s/h+2*m*h/3);

          if(gid < x_map->getGlobalNumElements()-1) {
            indices.push_back(gid+1);
            vals.push_back(-s/h+m*h/6);
          }        
          
          H_pp_out_crs->replaceGlobalValues(gid, indices, vals);
        }
      }
      H_pp_out_crs->fillComplete();
      Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> >   matWriter;
      matWriter.writeSparseFile("H_411", H_pp_out_crs);

      if(probList_->sublist("Hessian").sublist("Response 4").sublist("Parameter 1").isSublist("H_pp Solver")) {
        auto pl = probList_->sublist("Hessian").sublist("Response 4").sublist("Parameter 1").sublist("H_pp Solver");
        H_pp_out->initializeSolver(Teuchos::rcpFromRef(pl));
      }
      
    }
  }

  {
    int p_index=1;
    auto hess_g_pp = outArgs.supports(Thyra::ModelEvaluator<double>::OUT_ARG_hess_g_pp,5,p_index,p_index) ? outArgs.get_hess_g_pp(5,p_index,p_index) : Teuchos::null; 
    const Teuchos::RCP<MatrixBased_LOWS> H_pp_out =
      Teuchos::nonnull(hess_g_pp) ?
        Teuchos::rcp_dynamic_cast<MatrixBased_LOWS>(outArgs.get_hess_g_pp(5,p_index,p_index)):
        Teuchos::null;

    if (Teuchos::nonnull(H_pp_out)) {
      Teuchos::RCP<Tpetra_CrsMatrix> H_pp_out_crs =
        Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(ConverterT::getTpetraOperator(H_pp_out->getMatrix()), true);
      H_pp_out_crs->resumeFill();
      H_pp_out_crs->setAllToScalar(0.0);

      {
        Teuchos::Array<typename Tpetra_CrsGraph::global_ordinal_type> indices(3);
        const int nodeNumElements = x_map->getLocalNumElements();
         Teuchos::Array<double> vals(3);
        
        for (int i=0; i<nodeNumElements; i++) {

          auto gid = x_map->getGlobalElement(i);
          indices.clear();
          vals.clear();

               
          if(gid > 0) {
            indices.push_back(gid-1);
            vals.push_back(h/6.0);
          }

          indices.push_back(gid);
          if(gid == 0 || gid == x_map->getGlobalNumElements()-1)
            vals.push_back(h/3);
          else
            vals.push_back(2*h/3);

          if(gid < x_map->getGlobalNumElements()-1) {
            indices.push_back(gid+1);
            vals.push_back(h/6);
          }  
          H_pp_out_crs->replaceGlobalValues(gid, indices, vals);
        }
      }
      H_pp_out_crs->fillComplete();
      Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> >   matWriter;
      matWriter.writeSparseFile("H_511", H_pp_out_crs);

      if(probList_->sublist("Hessian").sublist("Response 5").sublist("Parameter 1").isSublist("H_pp Solver")) {
        auto pl = probList_->sublist("Hessian").sublist("Response 5").sublist("Parameter 1").sublist("H_pp Solver");
        H_pp_out->initializeSolver(Teuchos::rcpFromRef(pl));
      }
    }
  }

  if (Teuchos::nonnull(dfdp_0_out)) {
    Teuchos::RCP<DfDpOp> dfdp_op =
        Teuchos::rcp_dynamic_cast<DfDpOp>(dfdp_0_out);
    dfdp_op->set(p_vec_0);
  }

  if (Teuchos::nonnull(g0_out)) {
    const int nodeNumElements = p_map->getLocalNumElements();
    auto temp_vec = rcp(new Tpetra_Vector(x_map));
    for (int i=0; i<nodeNumElements; i++) {
      temp_vec->getDataNonConst()[i] = x[i]-std::pow(1+coords[i],3);
    }
    double norm2 = temp_vec->dot(*temp_vec);
    g0_out->getDataNonConst(0)[0] = 0.5*norm2;
  }

  const Teuchos::RCP<Thyra::MultiVectorBase<double>> dg0dx_base =
      outArgs.get_DgDx(0).getMultiVector();
  const Teuchos::RCP<Tpetra_MultiVector> dg0dx_out =
      Teuchos::nonnull(dg0dx_base) ?
          ConverterT::getTpetraMultiVector(dg0dx_base) :
          Teuchos::null;

  if (Teuchos::nonnull(dg0dx_out)) {
    const int nodeNumElements = x_map->getLocalNumElements();
    for (int i=0; i<nodeNumElements; i++) {
      dg0dx_out->getDataNonConst(0)[i] = x[i]-std::pow(1+coords[i],3);
    }
  }

  const Teuchos::RCP<Thyra::MultiVectorBase<double>> dg0dp0_base =
      outArgs.get_DgDp(0, 0).getMultiVector();
  const Teuchos::RCP<Tpetra_MultiVector> dg0dp0_out =
      Teuchos::nonnull(dg0dp0_base) ?
          ConverterT::getTpetraMultiVector(dg0dp_base) :
          Teuchos::null;

  if (Teuchos::nonnull(dg0dp0_out)) {
    dg0dp0_out->putScalar(0.0);
  }

  if (Teuchos::nonnull(f_hess_xx_v_out)) {
    f_hess_xx_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_xp_0_v_out)) {
    f_hess_xp_0_v_out->getVectorNonConst(0)->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_px_1_v_out)) {
    f_hess_px_1_v_out->getVectorNonConst(0)->putScalar(0);
  }   

  if (Teuchos::nonnull(f_hess_pp_00_v_out)) {
    TEUCHOS_ASSERT(Teuchos::nonnull(p_direction_0));
    const auto direction_p_0 = p_direction_0->getVector(0)->getData();
    const auto lag_mult = lag_multiplier_f_in->getData();
    const int nodeNumElements = p_map->getLocalNumElements();
    for (int i=0; i<nodeNumElements; i++) {
      f_hess_pp_00_v_out->getDataNonConst(0)[i] = -6.0*p_0[i]*direction_p_0[i]*lag_mult[i];
    }
  }

  if (Teuchos::nonnull(f_hess_pp_01_v_out)) {
    f_hess_pp_01_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_pp_10_v_out)) {
    f_hess_pp_10_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(f_hess_pp_11_v_out)) {
    f_hess_pp_11_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_xx_v_out)) {
    g_hess_xx_v_out->putScalar(0);
        TEUCHOS_ASSERT(Teuchos::nonnull(x_direction));
    const auto direction_x = x_direction->getVector(0)->getData();
    const int nodeNumElements = x_map->getLocalNumElements();
    for (int i=0; i<nodeNumElements; i++) {
      g_hess_xx_v_out->getDataNonConst(0)[i] = direction_x[i];
    }
  }

  if (Teuchos::nonnull(g_hess_xp_0_v_out)) {
    g_hess_xp_0_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_px_0_v_out)) {
    g_hess_px_0_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_pp_00_v_out)) {
    g_hess_pp_00_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_pp_01_v_out)) {
    g_hess_pp_01_v_out->putScalar(0);
  }

    if (Teuchos::nonnull(g_hess_pp_10_v_out)) {
    g_hess_pp_10_v_out->putScalar(0);
  }

  if (Teuchos::nonnull(g_hess_pp_11_v_out)) {
    g_hess_pp_11_v_out->putScalar(0);
  }
}

Thyra::ModelEvaluatorBase::InArgs<double>
MockModelEval_H_Tpetra::createInArgsImpl() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<double> result;
  result.setModelEvalDescription(this->description());

  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x, true);


  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_t, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_alpha, true);
  result.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_beta, true);

  result.set_Np_Ng(2,6);

  return result;
}

