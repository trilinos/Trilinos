// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_Epetra_StokhosSolver.hpp"
#include "Stokhos_Epetra.hpp"

Piro::Epetra::StokhosSolver::
StokhosSolver(const Teuchos::RCP<Teuchos::ParameterList>& piroParams_,
	      const Teuchos::RCP<const Epetra_Comm>& globalComm) :
  piroParams(piroParams_),
  sg_solver_factory(piroParams_, globalComm)
{
}

Piro::Epetra::StokhosSolver::~StokhosSolver()
{
}

void
Piro::Epetra::StokhosSolver::
setup(const Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
      const Teuchos::RCP<NOX::Epetra::Observer>& noxObserver)
{
  sg_nonlin_model = sg_solver_factory.createSGModel(model);
  const Teuchos::RCP<NOX::Epetra::Observer> sg_observer =
    sg_solver_factory.createSGObserver(noxObserver);
  const Teuchos::RCP<EpetraExt::ModelEvaluator> sg_block_solver =
    sg_solver_factory.createSGSolver(sg_nonlin_model, sg_observer);
  sg_solver = sg_solver_factory.createSGSolverAdapter(sg_block_solver);
}

void
Piro::Epetra::StokhosSolver::
resetSolverParameters(const Teuchos::ParameterList& new_solver_params)
{
  sg_solver_factory.resetSolverParameters(new_solver_params);
}

Teuchos::RCP<const Epetra_Comm>
Piro::Epetra::StokhosSolver::
getSpatialComm() const
{
  return sg_solver_factory.getSpatialComm();
}

Teuchos::RCP<const Epetra_Comm>
Piro::Epetra::StokhosSolver::
getStochasticComm() const
{
  return sg_solver_factory.getStochasticComm();
}

Teuchos::RCP<const EpetraExt::MultiComm>
Piro::Epetra::StokhosSolver::
getGlobalMultiComm() const
{
  return sg_solver_factory.getGlobalMultiComm();
}

Teuchos::RCP<const Epetra_Map> 
Piro::Epetra::StokhosSolver::get_x_map() const
{
  return sg_solver->get_x_map();
}

Teuchos::RCP<const Epetra_Map> 
Piro::Epetra::StokhosSolver::get_f_map() const
{
  return sg_solver->get_f_map();
}

Teuchos::RCP<const Epetra_Map> 
Piro::Epetra::StokhosSolver::get_p_map(int l) const
{
  return sg_solver->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> 
Piro::Epetra::StokhosSolver::get_g_map(int j) const
{
  return sg_solver->get_g_map(j);
}

Teuchos::RCP<const Epetra_Vector> 
Piro::Epetra::StokhosSolver::get_x_init() const
{
  return sg_solver->get_x_init();
}

Teuchos::RCP<const Epetra_Vector> 
Piro::Epetra::StokhosSolver::get_p_init(int l) const
{
  return sg_solver->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs 
Piro::Epetra::StokhosSolver::createInArgs() const
{
  return sg_solver->createInArgs();
}

EpetraExt::ModelEvaluator::OutArgs 
Piro::Epetra::StokhosSolver::createOutArgs() const
{
  return sg_solver->createOutArgs();
}

void 
Piro::Epetra::StokhosSolver::evalModel(const InArgs& inArgs,
			     const OutArgs& outArgs ) const
{
  sg_solver->evalModel(inArgs, outArgs);
}

Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> 
Piro::Epetra::StokhosSolver::create_g_sg(int l, Epetra_DataAccess CV, 
					 const Epetra_Vector* v) const 
{
  OutArgs outargs = sg_nonlin_model->createOutArgs();
  int ng = outargs.Ng();
  Piro::Epetra::StokhosSolverFactory::SG_METHOD sg_method = 
    sg_solver_factory.getSGMethod();
  if (sg_method != Piro::Epetra::StokhosSolverFactory::SG_NI && 
      sg_method != Piro::Epetra::StokhosSolverFactory::SG_MPNI && 
      piroParams->get<std::string>("Solver Type") == "NOX" && l == ng) {
    return sg_nonlin_model->create_x_sg(CV, v);
  }
  else
    return sg_nonlin_model->create_g_sg(l, CV, v);
}

Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> 
Piro::Epetra::StokhosSolver::create_g_mv_sg(int l, int num_vecs, 
					    Epetra_DataAccess CV, 
					    const Epetra_MultiVector* v) const
{
  OutArgs outargs = sg_nonlin_model->createOutArgs();
  int ng = outargs.Ng();
  Piro::Epetra::StokhosSolverFactory::SG_METHOD sg_method = 
    sg_solver_factory.getSGMethod();
  if (sg_method != Piro::Epetra::StokhosSolverFactory::SG_NI && 
      sg_method != Piro::Epetra::StokhosSolverFactory::SG_MPNI && 
      piroParams->get<std::string>("Solver Type") == "NOX" && l == ng) {
    return sg_nonlin_model->create_x_mv_sg(num_vecs, CV, v);
  }
  else
    return sg_nonlin_model->create_g_mv_sg(l, num_vecs, CV, v);
}

