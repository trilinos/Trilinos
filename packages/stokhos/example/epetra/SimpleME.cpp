// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "SimpleME.hpp"
#include "Teuchos_Assert.hpp"
#include "Stokhos_Epetra.hpp"

SimpleME::SimpleME(const Teuchos::RCP<const Epetra_Comm>& comm) 
{
  // Solution vector map
  x_map = Teuchos::rcp(new Epetra_Map(2, 0, *comm));

  // Overlapped solution vector map
  x_overlapped_map = Teuchos::rcp(new Epetra_LocalMap(2, 0, *comm));

  // Importer
  importer = Teuchos::rcp(new Epetra_Import(*x_overlapped_map, *x_map));

  // Initial guess, initialized to 1.5
  x_init = Teuchos::rcp(new Epetra_Vector(*x_map));
  x_init->PutScalar(1.5);

  // Overlapped solution vector
  x_overlapped = Teuchos::rcp(new Epetra_Vector(*x_overlapped_map));

  // Parameter vector map
  p_map = Teuchos::rcp(new Epetra_LocalMap(1, 0, *comm));

  // Initial parameters
  p_init = Teuchos::rcp(new Epetra_Vector(*p_map));
  (*p_init)[0] = 2.0;

  // Parameter names
  p_names = Teuchos::rcp(new Teuchos::Array<std::string>(1));
  (*p_names)[0] = "alpha";

  // Jacobian graph (dense 2x2 matrix)
  graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *x_map, 2));
  int indices[2];
  indices[0] = 0; indices[1] = 1;
  if (x_map->MyGID(0))
    graph->InsertGlobalIndices(0, 2, indices);
  if (x_map->MyGID(1))
    graph->InsertGlobalIndices(1, 2, indices);
  graph->FillComplete();
  graph->OptimizeStorage();
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
SimpleME::get_x_map() const
{
  return x_map;
}

Teuchos::RCP<const Epetra_Map>
SimpleME::get_f_map() const
{
  return x_map;
}

Teuchos::RCP<const Epetra_Map>
SimpleME::get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  SimpleME::get_p_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_map;
}

Teuchos::RCP<const Teuchos::Array<std::string> >
SimpleME::get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  SimpleME::get_p_names():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_names;
}

Teuchos::RCP<const Epetra_Vector>
SimpleME::get_x_init() const
{
  return x_init;
}

Teuchos::RCP<const Epetra_Vector>
SimpleME::get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  SimpleME::get_p_init():  " <<
                     "Invalid parameter index l = " << l << std::endl);
  
  return p_init;
}

Teuchos::RCP<Epetra_Operator>
SimpleME::create_W() const
{
  Teuchos::RCP<Epetra_CrsMatrix> A = 
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
  A->FillComplete();
  A->OptimizeStorage();
  return A;
}

EpetraExt::ModelEvaluator::InArgs
SimpleME::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription("Simple Model Evaluator");

  // Deterministic InArgs
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.set_Np(1);    // 1 parameter vector

  // Stochastic InArgs
  inArgs.setSupports(IN_ARG_x_sg,true);
  inArgs.setSupports(IN_ARG_p_sg, 0, true); // 1 SG parameter vector
  inArgs.setSupports(IN_ARG_sg_basis,true);
  inArgs.setSupports(IN_ARG_sg_quadrature,true);
  inArgs.setSupports(IN_ARG_sg_expansion,true);
  
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
SimpleME::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription("Simple Model Evaluator");

  // Deterministic OutArgs
  outArgs.set_Np_Ng(1, 0);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  
  // Stochastic OutArgs
  outArgs.setSupports(OUT_ARG_f_sg,true);
  outArgs.setSupports(OUT_ARG_W_sg,true);

  return outArgs;
}

void 
SimpleME::evalModel(const InArgs& inArgs, const OutArgs& outArgs) const
{
  //
  // Determinisic calculation
  //

  // Solution vector
  Teuchos::RCP<const Epetra_Vector> x = inArgs.get_x();
  if (x != Teuchos::null) 
    x_overlapped->Import(*x, *importer, Insert);
  double x0 = (*x_overlapped)[0];
  double x1 = (*x_overlapped)[1];

  // Parameters
  Teuchos::RCP<const Epetra_Vector> p = inArgs.get_p(0);
  if (p == Teuchos::null)
    p = p_init;
  double a = (*p)[0];
  
  // Residual
  // f = |  a*a  - x0 |
  //     | x1*x1 - x0 |
  // where a = p[0].
  Teuchos::RCP<Epetra_Vector> f = outArgs.get_f();
  if (f != Teuchos::null) {
    int row;
    double v;
    if (x_map->MyGID(0)) {
      row = 0;
      v = a*a - x0;
      f->ReplaceGlobalValues(1, &v, &row);
    }
    if (x_map->MyGID(1)) {
      row = 1;
      v = x1*x1 - x0;
      f->ReplaceGlobalValues(1, &v, &row);
    }
  }

  // Jacobian
  // J = | -1   0   |
  //     | -1  2*x1 |
  Teuchos::RCP<Epetra_Operator> W = outArgs.get_W();
  if (W != Teuchos::null) {
    Teuchos::RCP<Epetra_CrsMatrix> jac = 
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W, true);
    int indices[2] = { 0, 1 };
    double values[2];
    if (x_map->MyGID(0)) {
      values[0] = -1.0; values[1] = 0.0;
      jac->ReplaceGlobalValues(0, 2, values, indices);
    }
    if (x_map->MyGID(1)) {
      values[0] = -1.0; values[1] = 2.0*x1;
      jac->ReplaceGlobalValues(1, 2, values, indices);
    }
  }

  //
  // Stochastic Galerkin calculation
  //

  // Get stochastic expansion data
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
    inArgs.get_sg_basis();
  Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expn = 
    inArgs.get_sg_expansion();

  // Stochastic solution vector
  InArgs::sg_const_vector_t x_sg = inArgs.get_x_sg();
  Stokhos::OrthogPolyApprox<int,double> x0_sg(basis), x1_sg(basis);
  if (x_sg != Teuchos::null) {
    for (int i=0; i<basis->size(); i++) {
      x_overlapped->Import((*x_sg)[i], *importer, Insert);
      x0_sg[i] = (*x_overlapped)[0];
      x1_sg[i] = (*x_overlapped)[1];
    }
  }

  // Stochastic parameters
  InArgs::sg_const_vector_t p_sg = inArgs.get_p_sg(0);
  Stokhos::OrthogPolyApprox<int,double> a_sg(basis);
  if (p_sg != Teuchos::null) {
    for (int i=0; i<basis->size(); i++) {
      a_sg[i] = (*p_sg)[i][0];
    }
  }

  // Stochastic residual
  // f[i] = | <a*a - x0, psi_i>/<psi_i^2>   |
  //        | <x1*x1 - x0, psi_i>/<psi_i^2> |
  OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
  Stokhos::OrthogPolyApprox<int,double> tmp0_sg(basis), tmp1_sg(basis);
  if (f_sg != Teuchos::null) {
    int row;
    if (x_map->MyGID(0)) {
      row = 0;
      expn->times(tmp0_sg, a_sg, a_sg);
      expn->minus(tmp1_sg, tmp0_sg, x0_sg);
      for (int i=0; i<basis->size(); i++)
	(*f_sg)[i].ReplaceGlobalValues(1, &tmp1_sg[i], &row);
    }
    if (x_map->MyGID(1)) {
      row = 1;
      expn->times(tmp0_sg, x1_sg, x1_sg);
      expn->minus(tmp1_sg, tmp0_sg, x0_sg);
      for (int i=0; i<basis->size(); i++)
	(*f_sg)[i].ReplaceGlobalValues(1, &tmp1_sg[i], &row);
    }
  }

  // Stochastic Jacobian
  // J[0] = | -1     0    |,   J[i] = | 0     0    |,  i > 0
  //        | -1  2*x0[0] |           | 0  2*x0[i] |
  OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();
  if (W_sg != Teuchos::null) {
    for (int i=0; i<basis->size(); i++) {
      Teuchos::RCP<Epetra_CrsMatrix> jac = 
	Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_sg->getCoeffPtr(i), true);
      int indices[2] = { 0, 1 };
      double values[2];
      if (x_map->MyGID(0)) {
	if (i == 0)
	  values[0] = -1.0;
	else
	  values[0] = 0.0;
	values[1] = 0.0;
	jac->ReplaceGlobalValues(0, 2, values, indices);
      }
      if (x_map->MyGID(1)) {
	if (i == 0)
	  values[0] = -1.0;
	else
	  values[0] = 0.0;
	values[1] = 2.0*x1_sg[i]; 
	jac->ReplaceGlobalValues(1, 2, values, indices);
      }
    }
  }
}
