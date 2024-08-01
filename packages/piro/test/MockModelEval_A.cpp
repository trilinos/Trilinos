// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MockModelEval_A.hpp"

#include <iostream>

using Teuchos::RCP;
using Teuchos::rcp;

MockModelEval_A::MockModelEval_A(const MPI_Comm appComm)
{

#ifdef HAVE_MPI
    Comm = rcp(new Epetra_MpiComm(appComm));
#else
    Comm = rcp(new Epetra_SerialComm);
#endif

    //set up map and initial guess for solution vector
    const int vecLength = 4;
    x_map = rcp(new Epetra_Map(vecLength, 0, *Comm));
    x_vec = rcp(new Epetra_Vector(*x_map));
    x_dot_vec = rcp(new Epetra_Vector(*x_map));
    x_vec->PutScalar(3.0);
    x_dot_vec->PutScalar(1.0);

    //set up responses
    const int numResponses = 1;
    g_map = rcp(new Epetra_LocalMap(numResponses, 0, *Comm));

    //set up parameters
    const int numParameters= 2;
    p_map = rcp(new Epetra_LocalMap(numParameters, 0, *Comm));
    p_init = rcp(new Epetra_Vector(*p_map));
    p_lo = rcp(new Epetra_Vector(*p_map));
    p_up = rcp(new Epetra_Vector(*p_map));
    for (int i=0; i<numParameters; i++) {
      (*p_init)[i]= 1.0;
      (*p_lo)[i]= 0.1;
      (*p_up)[i]= 10.0;
    }

    //set up jacobian graph
    jacGraph = rcp(new Epetra_CrsGraph(Copy, *x_map, vecLength, true));
    std::vector<int> indices(vecLength);
    for (int i=0; i<vecLength; i++) indices[i]=i;
    for (int i=0; i<x_map->NumMyElements(); i++)
      jacGraph->InsertGlobalIndices(x_map->GID(i), vecLength, &indices[0]);
    jacGraph->FillComplete();
}

MockModelEval_A::~MockModelEval_A()
{
}

RCP<const Epetra_Map> MockModelEval_A::get_x_map() const
{
  return x_map;
}

RCP<const Epetra_Map> MockModelEval_A::get_f_map() const
{
  return x_map;
}

RCP<Epetra_Operator> MockModelEval_A::create_W() const
{
  return rcp(new Epetra_CrsMatrix(::Copy, *jacGraph));
}

RCP<const Epetra_Map> MockModelEval_A::get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  return p_map;
}

RCP<const  Teuchos::Array<std::string> > MockModelEval_A::get_p_names(int l) const
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


RCP<const Epetra_Map> MockModelEval_A::get_g_map(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(j != 0, std::logic_error,
                     std::endl <<
                     "Error!  MockModelEval_A::get_g_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     j << std::endl);
  return g_map;
}

RCP<const Epetra_Vector> MockModelEval_A::get_x_init() const
{
  return x_vec;
}

RCP<const Epetra_Vector> MockModelEval_A::get_x_dot_init() const
{
  return x_dot_vec;
}

RCP<const Epetra_Vector> MockModelEval_A::get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_init() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  return p_init;
}

RCP<const Epetra_Vector> MockModelEval_A::get_p_lower_bounds(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_lower_bounds() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  return p_lo;
}

RCP<const Epetra_Vector> MockModelEval_A::get_p_upper_bounds(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_upper_bounds() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  return p_up;
}

EpetraExt::ModelEvaluator::InArgs MockModelEval_A::createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(1);
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.setSupports(IN_ARG_t,true);

  // This ModelEvaluator only supports explicit time integration...
  inArgs.setSupports(IN_ARG_x_dot,true);
  inArgs.setSupports(IN_ARG_alpha,true);
  inArgs.setSupports(IN_ARG_beta,true);

  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs MockModelEval_A::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1, 1);
  outArgs.setSupports(OUT_ARG_DgDp, 0, 0, DerivativeSupport(DERIV_MV_BY_COL));
  outArgs.setSupports(OUT_ARG_DgDx, 0, DerivativeSupport(DERIV_TRANS_MV_BY_ROW));

  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_DfDp, 0, DerivativeSupport(DERIV_MV_BY_COL));
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties( DerivativeProperties(
      DERIV_LINEARITY_UNKNOWN, DERIV_RANK_FULL, true));


  return outArgs;
}

void MockModelEval_A::evalModel( const InArgs& inArgs,
                              const OutArgs& outArgs ) const
{

  // Parse InArgs
  RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  if (!p_in.get()) std::cerr << "ERROR: MockModelEval_A requires p as inargs\n";

  RCP<const Epetra_Vector> x_in = inArgs.get_x();
  if (!x_in.get()) std::cerr << "ERROR: MockModelEval_A requires x as inargs\n";
  int vecLength = x_in->GlobalLength();
  int myVecLength = x_in->MyLength();

  // Parse OutArgs

  RCP<Epetra_Vector> f_out = outArgs.get_f();
  RCP<Epetra_Vector> g_out = outArgs.get_g(0);
  Teuchos::RCP<Epetra_Operator> W_out = outArgs.get_W();
  Teuchos::RCP<Epetra_MultiVector> dfdp_out;
  if (outArgs.Np() > 0)
    dfdp_out = outArgs.get_DfDp(0).getMultiVector();
  RCP<Epetra_MultiVector> dgdp_out;
  dgdp_out = outArgs.get_DgDp(0,0).getMultiVector();
  RCP<Epetra_MultiVector> dgdx_out;
  dgdx_out = outArgs.get_DgDx(0).getMultiVector();

  if (f_out != Teuchos::null) {
    for (int i=0; i<myVecLength; i++) {
      int gid = x_in->Map().GID(i);
      if (gid==0) // f_0 = (x_0)^2 - p_0
       (*f_out)[i] = (*x_in)[i] * (*x_in)[i] -  (*p_in)[i];
      else // f_i = x_i^2 - (i+p_1)^2 (for i != 0)
       (*f_out)[i] = (*x_in)[i] * (*x_in)[i] - (gid + (*p_in)[1])*(gid + (*p_in)[1]);
    }
  }
  if (W_out != Teuchos::null) {
    Teuchos::RCP<Epetra_CrsMatrix> W_out_crs =
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_out, true);
    W_out_crs->PutScalar(0.0);

    double diag=0.0;
    for (int i=0; i<myVecLength; i++) {
      diag = 2.0 * (*x_in)[i];
      W_out_crs->ReplaceMyValues(i, 1, &diag, &i);
    }
  }

  if (dfdp_out != Teuchos::null) {
    dfdp_out->PutScalar(0.0);
    for (int i=0; i<myVecLength; i++) {
      const int gid = x_in->Map().GID(i);
      if   (gid==0) (*dfdp_out)[0][i] = -1.0;
      else          (*dfdp_out)[1][i] =  -2.0* (gid + (*p_in)[1]);
    }
  }

  // Response: g = 0.5*(Sum(x)-Sum(p)-12)^2 + 0.5*(p0-1)^2
  // min g(x(p), p) s.t. f(x, p) = 0 reached for p_0 = 1, p_1 = 3

  double term1, term2;
  x_in->MeanValue(&term1);
  term1 =  vecLength * term1 - ((*p_in)[0] + (*p_in)[1]) - 12.0;
  term2 = (*p_in)[0] - 1.0;

  if (!is_null(g_out)) {
    (*g_out)[0] = 0.5*term1*term1 + 0.5*term2*term2;
  }

  if (dgdx_out != Teuchos::null) {
    dgdx_out->PutScalar(term1);
  }
  if (dgdp_out != Teuchos::null) {
    dgdp_out->PutScalar(0.0);
    (*dgdp_out)[0][0] = -term1 + term2;
    (*dgdp_out)[0][1] = -term1;
  }

  // Modify for time dependent (implicit time integration or eigensolves)
  const RCP<const Epetra_Vector> x_dot = inArgs.get_x_dot();

  if (x_dot.get()) {
    // Velocity provided: Time dependent problem
    double alpha = inArgs.get_alpha();
    double beta = inArgs.get_beta();
    if (alpha==0.0 && beta==0.0) {
      std::cerr << "MockModelEval Warning: alpha=beta=0 -- setting beta=1\n";
      beta = 1.0;
    }

    if (f_out != Teuchos::null) {
      // f(x, x_dot) = f(x) - x_dot
      for (int i=0; i<myVecLength; i++) {
        (*f_out)[i] = -(*x_dot)[i] + (*f_out)[i];
      }
    }
    if (W_out != Teuchos::null) {
      // W(x, x_dot) = beta * W(x) - alpha * Id
      const Teuchos::RCP<Epetra_CrsMatrix> W_out_crs =
        Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_out, true);
      W_out_crs->Scale(beta);

      const double diag = -alpha;
      for (int i=0; i<myVecLength; i++) {
        W_out_crs->SumIntoMyValues(i, 1, &diag, &i);
      }
    }
  }
}
