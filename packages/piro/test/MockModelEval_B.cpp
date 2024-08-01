// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MockModelEval_B.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

MockModelEval_B::MockModelEval_B(const MPI_Comm appComm) 
{

#ifdef HAVE_MPI
    Comm = rcp(new Epetra_MpiComm(appComm));
#else
    Comm = rcp(new Epetra_SerialComm);
#endif

    //set up map and initial guess for solution vector
    const int vecLength = 1;
    x_map = rcp(new Epetra_Map(vecLength, 0, *Comm));
    x_vec = rcp(new Epetra_Vector(*x_map));
    x_vec->PutScalar(0.0);
    x_dot_vec = rcp(new Epetra_Vector(*x_map));
    x_dot_vec->PutScalar(1.0);
    x_dotdot_vec = rcp(new Epetra_Vector(*x_map));
    x_dotdot_vec->PutScalar(0.0);

    //set up responses
    const int numResponses = 1;
    g_map = rcp(new Epetra_LocalMap(numResponses, 0, *Comm));

    //set up parameters
    const int numParameters= 1;
    p_map = rcp(new Epetra_LocalMap(numParameters, 0, *Comm));
    p_init = rcp(new Epetra_Vector(*p_map));
    for (int i=0; i<numParameters; i++) (*p_init)[i]= 1.0;

    Epetra_CrsGraph graph(Copy, *x_map, 1);
    int z=0;
    graph.InsertGlobalIndices(0,1,&z);
    graph.FillComplete();

    W = Teuchos::rcp(new Epetra_CrsMatrix(Copy, graph));

}

MockModelEval_B::~MockModelEval_B()
{
}

RCP<const Epetra_Map> MockModelEval_B::get_x_map() const
{
  return x_map;
}

RCP<const Epetra_Map> MockModelEval_B::get_f_map() const
{
  return x_map;
}

RCP<Epetra_Operator> MockModelEval_B::create_W() const
{
  return W;
}

RCP<const Epetra_Map> MockModelEval_B::get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  return p_map;
}

RCP<const  Teuchos::Array<std::string> > MockModelEval_B::get_p_names(int l) const
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


RCP<const Epetra_Map> MockModelEval_B::get_g_map(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(j != 0, std::logic_error,
                     std::endl <<
                     "Error!  MockModelEval_B::get_g_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     j << std::endl);
  return g_map;
}

RCP<const Epetra_Vector> MockModelEval_B::get_x_init() const
{
  return x_vec;
}

RCP<const Epetra_Vector> MockModelEval_B::get_x_dot_init() const
{
  return x_dot_vec;
}

RCP<const Epetra_Vector> MockModelEval_B::get_x_dotdot_init() const
{
  return x_dotdot_vec;
}

RCP<const Epetra_Vector> MockModelEval_B::get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  return p_init;
}

EpetraExt::ModelEvaluator::InArgs MockModelEval_B::createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(1);
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.setSupports(IN_ARG_x_dot,true);
  inArgs.setSupports(IN_ARG_x_dotdot,true);
  inArgs.setSupports(IN_ARG_t,true);
  inArgs.setSupports(IN_ARG_omega,true);
  inArgs.setSupports(IN_ARG_alpha,true);
  inArgs.setSupports(IN_ARG_beta,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs MockModelEval_B::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1, 1);

  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties( DerivativeProperties(
      DERIV_LINEARITY_UNKNOWN, DERIV_RANK_FULL, true));

  return outArgs;
}

void MockModelEval_B::evalModel( const InArgs& inArgs,
                              const OutArgs& outArgs ) const
{

  // Parse InArgs
  RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  if (!p_in.get()) std::cout << "ERROR: MockModelEval_B requires p as inargs" << std::endl;

  RCP<const Epetra_Vector> x_in = inArgs.get_x();
  if (!x_in.get()) std::cout << "ERROR: MockModelEval_B requires x as inargs" << std::endl;
  int myVecLength = x_in->MyLength();

  RCP<const Epetra_Vector> x_dot_in = inArgs.get_x_dot();
  double alpha = inArgs.get_alpha();

  RCP<const Epetra_Vector> x_dotdot_in = inArgs.get_x_dotdot();
  double omega = inArgs.get_omega();

  // Parse OutArgs

  RCP<Epetra_Vector> f_out = outArgs.get_f(); 
  RCP<Epetra_Vector> g_out = outArgs.get_g(0); 
  RCP<Epetra_CrsMatrix> W_out = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(outArgs.get_W()); 

  double damping=0.0;
  if (f_out != Teuchos::null) {
    for (int i=0; i<myVecLength; i++) {
       (*f_out)[i] = -(*p_in)[0];
    }
    if (x_dotdot_in != Teuchos::null) {
       for (int i=0; i<myVecLength; i++) {
       (*f_out)[i] = (*x_dotdot_in)[i] - (*f_out)[i];
       }
    }
    if (x_dot_in != Teuchos::null) {
       for (int i=0; i<myVecLength; i++) {
       (*f_out)[i] += damping*(*x_dot_in)[i];
       }
    }
  }
  if (W_out != Teuchos::null) {
    int z=0;
    if (omega==0.0) throw "omega=0.0";
    W_out->ReplaceGlobalValues(0, 1, &omega, &z);

    if (x_dot_in != Teuchos::null) {
      double da = damping*alpha;
      W_out->SumIntoGlobalValues(0, 1, &da, &z);
    }
    W_out->FillComplete();
  }

  if (g_out != Teuchos::null) x_in->MeanValue(&(*g_out)[0]);
} 
