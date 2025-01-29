// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MockModelEval_C.hpp"
#include "Piro_ConfigDefs.hpp"

#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"

using Teuchos::RCP;
using Teuchos::rcp;

MockModelEval_C::MockModelEval_C(const Teuchos::RCP<const Epetra_Comm>& comm_) :
  comm(comm_)
{
  //set up map and initial guess for solution vector
  const int vecLength = 1;
  x_map = rcp(new Epetra_Map(vecLength, 0, *comm));
  x_init = rcp(new Epetra_Vector(*x_map));
  x_init->PutScalar(1.0);
  x_dot_init = rcp(new Epetra_Vector(*x_map));
  x_dot_init->PutScalar(0.0);

  //set up responses
  const int numResponses = 1;
  g_map = rcp(new Epetra_LocalMap(numResponses, 0, *comm));

  //set up parameters
  const int numParameters = 1;
  p_map = rcp(new Epetra_LocalMap(numParameters, 0, *comm));
  p_init = rcp(new Epetra_Vector(*p_map));
  for (int i=0; i<numParameters; i++) (*p_init)[i]= i+1;
  
  //setup Jacobian graph
  graph = rcp(new Epetra_CrsGraph(Copy, *x_map, 1));
  for (int i=0; i<vecLength; i++)
    graph->InsertGlobalIndices(i,1,&i);
  graph->FillComplete();
}

MockModelEval_C::~MockModelEval_C()
{
}

RCP<const Epetra_Map> 
MockModelEval_C::get_x_map() const
{
  return x_map;
}

RCP<const Epetra_Map> 
MockModelEval_C::get_f_map() const
{
  return x_map;
}

RCP<const Epetra_Map> 
MockModelEval_C::get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  return p_map;
}

RCP<const  Teuchos::Array<std::string> > 
MockModelEval_C::get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_names() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);

  Teuchos::Ordinal num_p = p_init->MyLength();
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

RCP<const Epetra_Map> 
MockModelEval_C::get_g_map(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(j != 0, std::logic_error,
                     std::endl <<
                     "Error!  MockModelEval_C::get_g_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     j << std::endl);
  return g_map;
}

RCP<const Epetra_Vector> 
MockModelEval_C::get_x_init() const
{
  return x_init;
}

RCP<const Epetra_Vector> 
MockModelEval_C::get_x_dot_init() const
{
  return x_dot_init;
}

RCP<const Epetra_Vector> 
MockModelEval_C::get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  return p_init;
}

RCP<Epetra_Operator> 
MockModelEval_C::create_W() const
{
  RCP<Epetra_CrsMatrix> W = rcp(new Epetra_CrsMatrix(Copy, *graph));
  W->FillComplete();
  return W;
}

EpetraExt::ModelEvaluator::InArgs 
MockModelEval_C::createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(1);
  inArgs.setSupports(IN_ARG_x, true);

  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs 
MockModelEval_C::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1, 1);

  outArgs.setSupports(OUT_ARG_f, true);
  outArgs.setSupports(OUT_ARG_W, true);
  outArgs.set_W_properties( DerivativeProperties(
      DERIV_LINEARITY_UNKNOWN, DERIV_RANK_FULL, true));
  outArgs.setSupports(OUT_ARG_DfDp, 0, DERIV_MV_BY_COL);
  outArgs.setSupports(OUT_ARG_DgDx, 0, DERIV_TRANS_MV_BY_ROW);
  outArgs.setSupports(OUT_ARG_DgDp, 0, 0, DERIV_MV_BY_COL);

  return outArgs;
}

void 
MockModelEval_C::evalModel(const InArgs& inArgs, const OutArgs& outArgs) const
{
  int proc = comm->MyPID();

  // 
  // Deterministic calculation
  //

  // Parse InArgs
  RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  if (p_in == Teuchos::null)
    p_in = p_init;

  RCP<const Epetra_Vector> x_in = inArgs.get_x();

  // Parse OutArgs
  RCP<Epetra_Vector> f_out = outArgs.get_f(); 
  if (f_out != Teuchos::null) {
    double p = (*p_in)[0];
    if (proc == 0) {
      double x = (*x_in)[0];
      (*f_out)[0] = 0.5*(x*x - p*p);
    }
  }

  RCP<Epetra_CrsMatrix> W_out = 
    Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(outArgs.get_W()); 
  if (W_out != Teuchos::null) {
    if (proc == 0) {
      double x = (*x_in)[0];
      int i = 0;
      W_out->ReplaceMyValues(i, 1, &x, &i);
    }
  }

  RCP<Epetra_MultiVector> dfdp = outArgs.get_DfDp(0).getMultiVector();
  if (dfdp != Teuchos::null) {
    double p = (*p_in)[0];
    if (proc == 0)
      (*dfdp)[0][0] = -p;
  }

  RCP<Epetra_Vector> g_out = outArgs.get_g(0); 
  if (g_out != Teuchos::null) {
    double p = (*p_in)[0];
    if (proc == 0) {
      double x = (*x_in)[0];
      (*g_out)[0] = 0.5*(x*x + p*p);
    }
  }
    

  RCP<Epetra_MultiVector> dgdx = outArgs.get_DgDx(0).getMultiVector();
  if (dgdx != Teuchos::null) {
    if (proc == 0) {
      double x = (*x_in)[0];
      (*dgdx)[0][0] = x;
    }
  }

  RCP<Epetra_MultiVector> dgdp = outArgs.get_DgDp(0,0).getMultiVector();
  if (dgdp != Teuchos::null) {
    double p = (*p_in)[0];
    if (proc == 0) {
      (*dgdp)[0][0] = p;
    }
  }
} 
