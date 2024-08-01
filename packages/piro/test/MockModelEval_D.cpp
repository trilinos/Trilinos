// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MockModelEval_D.hpp"
#include "Piro_ConfigDefs.hpp"

#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"

#ifdef HAVE_PIRO_STOKHOS
#include "Stokhos_Epetra.hpp"
#endif

using Teuchos::RCP;
using Teuchos::rcp;

MockModelEval_D::MockModelEval_D(const Teuchos::RCP<const Epetra_Comm>& comm_) :
  comm(comm_)
{
  //set up map and initial guess for solution vector
  x_map = rcp(new Epetra_Map(1, 0, *comm));
  x_init = rcp(new Epetra_Vector(*x_map));
  x_init->PutScalar(0.0);

  //set up responses
  g_map = rcp(new Epetra_LocalMap(1, 0, *comm));

  //set up parameters
  p1_map = rcp(new Epetra_LocalMap(1, 0, *comm));
  p1_init = rcp(new Epetra_Vector(*p1_map));
  (*p1_init)[0]= -0.41;

  p2_map = rcp(new Epetra_LocalMap(1, 0, *comm));
  p2_init = rcp(new Epetra_Vector(*p2_map));
  (*p2_init)[0]= 2.0;
  
  //setup Jacobian graph
  graph = rcp(new Epetra_CrsGraph(Copy, *x_map, 1));
  int index = 0;
  graph->InsertGlobalIndices(0,1,&index);
  graph->FillComplete();
}

MockModelEval_D::
~MockModelEval_D()
{
}

RCP<const Epetra_Map> 
MockModelEval_D::
get_x_map() const
{
  return x_map;
}

RCP<const Epetra_Map> 
MockModelEval_D::
get_f_map() const
{
  return x_map;
}

RCP<const Epetra_Map> 
MockModelEval_D::
get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l < 0 || l > 1, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map():  invalid  " <<
                     " parameter index l = " << l << std::endl);
  if (l==0) 
    return p1_map;
  else
    return p2_map;
}

RCP<const  Teuchos::Array<std::string> > 
MockModelEval_D::
get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l < 0 || l > 1, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_names():  invalid  " <<
                     " parameter index l = " << l << std::endl);

  RCP<Teuchos::Array<std::string> > p_names = 
      rcp(new Teuchos::Array<std::string>(1) );
  if (l == 0) 
    (*p_names)[0] = "p";
  else
    (*p_names)[0] = "xi";

  return p_names;
}

RCP<const Epetra_Map> 
MockModelEval_D::
get_g_map(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(j != 0, std::logic_error,
                     std::endl <<
                     "Error!  MockModelEval_D::get_g_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     j << std::endl);
  return g_map;
}

RCP<const Epetra_Vector> 
MockModelEval_D::
get_x_init() const
{
  return x_init;
}

RCP<const Epetra_Vector> 
MockModelEval_D::
get_x_dot_init() const
{
  return Teuchos::null;
}

RCP<const Epetra_Vector> 
MockModelEval_D::
get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l < 0 || l > 1, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_init():  invalid  " <<
                     " parameter index l = " << l << std::endl);
  if (l == 0)
    return p1_init;
  else
    return p2_init;
}

RCP<Epetra_Operator> 
MockModelEval_D::
create_W() const
{
  RCP<Epetra_CrsMatrix> W = rcp(new Epetra_CrsMatrix(Copy, *graph));
  W->FillComplete();
  return W;
}

EpetraExt::ModelEvaluator::InArgs 
MockModelEval_D::
createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(2);
  inArgs.setSupports(IN_ARG_x, true);

#ifdef HAVE_PIRO_STOKHOS
  inArgs.setSupports(IN_ARG_x_sg, true);
  inArgs.setSupports(IN_ARG_p_sg, 0, true);
  inArgs.setSupports(IN_ARG_p_sg, 1, true);
  inArgs.setSupports(IN_ARG_sg_basis, true);
  inArgs.setSupports(IN_ARG_sg_quadrature, true);
  inArgs.setSupports(IN_ARG_sg_expansion, true);
#endif

  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs 
MockModelEval_D::
createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(2, 1);

  outArgs.setSupports(OUT_ARG_f, true);
  outArgs.setSupports(OUT_ARG_W, true);
  outArgs.set_W_properties( DerivativeProperties(
      DERIV_LINEARITY_UNKNOWN, DERIV_RANK_FULL, true));
  outArgs.setSupports(OUT_ARG_DfDp, 0, DERIV_MV_BY_COL);
  outArgs.setSupports(OUT_ARG_DfDp, 1, DERIV_MV_BY_COL);
  outArgs.setSupports(OUT_ARG_DgDx, 0, DERIV_TRANS_MV_BY_ROW);
  outArgs.setSupports(OUT_ARG_DgDp, 0, 0, DERIV_MV_BY_COL);
  outArgs.setSupports(OUT_ARG_DgDp, 0, 1, DERIV_MV_BY_COL);

#ifdef HAVE_PIRO_STOKHOS
  outArgs.setSupports(OUT_ARG_f_sg, true);
  outArgs.setSupports(OUT_ARG_W_sg, true);
  outArgs.setSupports(OUT_ARG_g_sg, 0, true);
  outArgs.setSupports(OUT_ARG_DfDp_sg, 0, DERIV_MV_BY_COL);
  outArgs.setSupports(OUT_ARG_DfDp_sg, 1, DERIV_MV_BY_COL);
  outArgs.setSupports(OUT_ARG_DgDx_sg, 0, DERIV_TRANS_MV_BY_ROW);
  outArgs.setSupports(OUT_ARG_DgDp_sg, 0, 0, DERIV_MV_BY_COL);
  outArgs.setSupports(OUT_ARG_DgDp_sg, 0, 1, DERIV_MV_BY_COL);
#endif

  return outArgs;
}

void 
MockModelEval_D::
evalModel(const InArgs& inArgs, const OutArgs& outArgs) const
{
  int proc = comm->MyPID();

  // 
  // Deterministic calculation
  //

  // Parse InArgs
  RCP<const Epetra_Vector> p1_in = inArgs.get_p(0);
  if (p1_in == Teuchos::null)
    p1_in = p1_init;
  RCP<const Epetra_Vector> p2_in = inArgs.get_p(1);
  if (p2_in == Teuchos::null)
    p2_in = p2_init;

  RCP<const Epetra_Vector> x_in = inArgs.get_x();

  // Parse OutArgs
  RCP<Epetra_Vector> f_out = outArgs.get_f(); 
  if (f_out != Teuchos::null) {
    double p = (*p1_in)[0];
    double xi = (*p2_in)[0];
    if (proc == 0) {
      double x = (*x_in)[0];
      (*f_out)[0] = x - p + xi;
    }
  }

  RCP<Epetra_CrsMatrix> W_out = 
    Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(outArgs.get_W()); 
  if (W_out != Teuchos::null) {
    if (proc == 0) {
      double val = 1.0;
      int i = 0;
      W_out->ReplaceMyValues(i, 1, &val, &i);
    }
  }

  RCP<Epetra_MultiVector> dfdp1 = outArgs.get_DfDp(0).getMultiVector();
  if (dfdp1 != Teuchos::null) {
    if (proc == 0)
      (*dfdp1)[0][0] = -1.0;
  }
  RCP<Epetra_MultiVector> dfdp2 = outArgs.get_DfDp(1).getMultiVector();
  if (dfdp2 != Teuchos::null) {
    if (proc == 0)
      (*dfdp2)[0][0] = 1.0;
  }

  RCP<Epetra_Vector> g_out = outArgs.get_g(0); 
  if (g_out != Teuchos::null) {
    if (proc == 0) {
      double x = (*x_in)[0];
      (*g_out)[0] = 1.0 / x;
    }
  }
    

  RCP<Epetra_MultiVector> dgdx = outArgs.get_DgDx(0).getMultiVector();
  if (dgdx != Teuchos::null) {
    if (proc == 0) {
      double x = (*x_in)[0];
      (*dgdx)[0][0] = -1.0 / (x*x);
    }
  }

  RCP<Epetra_MultiVector> dgdp1 = outArgs.get_DgDp(0,0).getMultiVector();
  if (dgdp1 != Teuchos::null) {
    if (proc == 0) {
      (*dgdp1)[0][0] = 0.0;
    }
  }
  RCP<Epetra_MultiVector> dgdp2 = outArgs.get_DgDp(0,1).getMultiVector();
  if (dgdp2 != Teuchos::null) {
    if (proc == 0) {
      (*dgdp2)[0][0] = 0.0;
    }
  }

  // 
  // Stochastic calculation
  //

#ifdef HAVE_PIRO_STOKHOS
  // Parse InArgs
  RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
    inArgs.get_sg_basis();
  RCP<Stokhos::OrthogPolyExpansion<int,double> > expn = 
    inArgs.get_sg_expansion();
  InArgs::sg_const_vector_t x_sg = inArgs.get_x_sg();
  InArgs::sg_const_vector_t p1_sg = inArgs.get_p_sg(0);
  InArgs::sg_const_vector_t p2_sg = inArgs.get_p_sg(1);

  // Parse OutArgs
  OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
  if (f_sg != Teuchos::null && proc == 0) {
    for (int block=0; block<f_sg->size(); block++) {
      (*f_sg)[block][0] = 
	(*x_sg)[block][0] - (*p1_sg)[block][0] + (*p2_sg)[block][0];
    }
  }

  OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();
  if (W_sg != Teuchos::null) {
    W_sg->init(0.0);
    Teuchos::RCP<Epetra_CrsMatrix> W = 
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_sg->getCoeffPtr(0), 
						  true);
    if (proc == 0) {
      int i = 0;
      double val = 1.0;
      W->ReplaceMyValues(i, 1, &val, &i);
    }
  }

  RCP<Stokhos::EpetraMultiVectorOrthogPoly> dfdp1_sg = 
    outArgs.get_DfDp_sg(0).getMultiVector();
  if (dfdp1_sg != Teuchos::null) {
    dfdp1_sg->init(0.0);
    if (proc == 0) {
      (*dfdp1_sg)[0][0][0] = -1.0;
    }
  }
  RCP<Stokhos::EpetraMultiVectorOrthogPoly> dfdp2_sg = 
    outArgs.get_DfDp_sg(1).getMultiVector();
  if (dfdp2_sg != Teuchos::null) {
    dfdp2_sg->init(0.0);
    if (proc == 0) {
      (*dfdp2_sg)[0][0][0] = 1.0;
    }
  }

  Stokhos::OrthogPolyApprox<int,double> x(basis);
  if (x_sg != Teuchos::null && proc == 0) {
    for (int i=0; i<basis->size(); i++) {
      x[i] = (*x_sg)[i][0];
    }   
  }

  OutArgs::sg_vector_t g_sg = outArgs.get_g_sg(0); 
  if (g_sg != Teuchos::null && proc == 0) {
    Stokhos::OrthogPolyApprox<int,double> xinv(basis);
    expn->divide(xinv, 1.0, x);
    for (int block=0; block<g_sg->size(); block++) {
      (*g_sg)[block][0] = xinv[block];
    }
  }

  RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdx_sg = 
    outArgs.get_DgDx_sg(0).getMultiVector();
  if (dgdx_sg != Teuchos::null && proc == 0) {
    Stokhos::OrthogPolyApprox<int,double> x2(basis), x2inv(basis);
    expn->times(x2, x, x);
    expn->divide(x2inv, -1.0, x2);
    for (int block=0; block<dgdx_sg->size(); block++) {
      (*dgdx_sg)[block][0][0] = x2inv[block];
    }
  }

  RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdp1_sg = 
    outArgs.get_DgDp_sg(0,0).getMultiVector();
  if (dgdp1_sg != Teuchos::null) {
    dgdp1_sg->init(0.0);
  }
  RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdp2_sg = 
    outArgs.get_DgDp_sg(0,1).getMultiVector();
  if (dgdp2_sg != Teuchos::null) {
    dgdp2_sg->init(0.0);
  }
#endif
} 
