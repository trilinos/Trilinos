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
#include "Stokhos_Sacado.hpp"
#include "Sacado.hpp"

namespace {

  template <typename ScalarA, typename ScalarX>
  void func(const ScalarA& a, const ScalarX x[2], ScalarX y[2]) {
    y[0] = a*a       - x[0];
    y[1] = x[1]*x[1] - x[0];
  }

}

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
  if (x != Teuchos::null) {
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
      double x[2] = { x0, x1 };
      double y[2];
      func(a, x, y);

      if (x_map->MyGID(0)) {
        int row = 0;
        f->ReplaceGlobalValues(1, &y[0], &row);
      }
      if (x_map->MyGID(1)) {
        int row = 1;
        f->ReplaceGlobalValues(1, &y[1], &row);
      }
    }

    // Jacobian
    // J = | -1   0   |
    //     | -1  2*x1 |
    Teuchos::RCP<Epetra_Operator> W = outArgs.get_W();
    if (W != Teuchos::null) {
      typedef Sacado::Fad::SFad<double,2> fad_type;
      fad_type x[2], y[2];
      x[0] = fad_type(2, 0, x0);
      x[1] = fad_type(2, 1, x1);
      func(a, x, y);

      Teuchos::RCP<Epetra_CrsMatrix> jac =
        Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W, true);
      int indices[2] = { 0, 1 };
      if (x_map->MyGID(0)) {
        int row = 0;
        jac->ReplaceGlobalValues(row, 2, y[0].dx(), indices);
      }
      if (x_map->MyGID(1)) {
        int row = 1;
        jac->ReplaceGlobalValues(row, 2, y[1].dx(), indices);
      }
    }
  }

  //
  // Stochastic Galerkin calculation
  //

  // Stochastic solution vector
  InArgs::sg_const_vector_t x_sg = inArgs.get_x_sg();
  if (x_sg != Teuchos::null) {

    // Get stochastic expansion data
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis =
      inArgs.get_sg_basis();
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expn =
      inArgs.get_sg_expansion();
    typedef Stokhos::StandardStorage<int,double> storage_type;
    typedef Sacado::PCE::OrthogPoly<double, storage_type> pce_type;

    pce_type x0(expn), x1(expn);
    for (int i=0; i<basis->size(); i++) {
      x_overlapped->Import((*x_sg)[i], *importer, Insert);
      x0.fastAccessCoeff(i) = (*x_overlapped)[0];
      x1.fastAccessCoeff(i) = (*x_overlapped)[1];
    }

    // Stochastic parameters
    InArgs::sg_const_vector_t p_sg = inArgs.get_p_sg(0);
    pce_type a(expn);
    if (p_sg != Teuchos::null) {
      for (int i=0; i<basis->size(); i++) {
        a.fastAccessCoeff(i) = (*p_sg)[i][0];
      }
    }

    // Stochastic residual
    // f[i] = | <a*a - x0, psi_i>/<psi_i^2>   |
    //        | <x1*x1 - x0, psi_i>/<psi_i^2> |
    OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
    if (f_sg != Teuchos::null) {
      pce_type x[2] = { x0, x1 };
      pce_type y[2];
      func(a, x, y);

      if (x_map->MyGID(0)) {
        int row = 0;
        for (int i=0; i<basis->size(); i++) {
          double c = y[0].coeff(i);
          (*f_sg)[i].ReplaceGlobalValues(1, &c, &row);
        }
      }
      if (x_map->MyGID(1)) {
        int row = 1;
        for (int i=0; i<basis->size(); i++) {
          double c = y[1].coeff(i);
          (*f_sg)[i].ReplaceGlobalValues(1, &c, &row);
        }
      }
    }

    // Stochastic Jacobian
    // J[0] = | -1     0    |,   J[i] = | 0     0    |,  i > 0
    //        | -1  2*x0[0] |           | 0  2*x0[i] |
    OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();
    if (W_sg != Teuchos::null) {
      typedef Sacado::Fad::SFad<pce_type,2> fad_type;
      fad_type x[2], y[2];
      x[0] = fad_type(2, 0, x0);
      x[1] = fad_type(2, 1, x1);
      func(a, x, y);

      for (int i=0; i<basis->size(); i++) {
        Teuchos::RCP<Epetra_CrsMatrix> jac =
          Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_sg->getCoeffPtr(i),
                                                      true);
        int indices[2] = { 0, 1 };
        if (x_map->MyGID(0)) {
          int row = 0;
          double values[2] = { y[0].dx(0).coeff(i), y[0].dx(1).coeff(i) };
          jac->ReplaceGlobalValues(row, 2, values, indices);
        }
        if (x_map->MyGID(1)) {
          int row = 1;
          double values[2] = { y[1].dx(0).coeff(i), y[1].dx(1).coeff(i) };
          jac->ReplaceGlobalValues(row, 2, values, indices);
        }
      }
    }
  }
}
