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

#include "Piro_Epetra_StokhosMPSolver.hpp"

#include "Piro_Epetra_SolverFactory.hpp"
#include "Piro_ExtensibleFactory.hpp"

#include "Stokhos_Epetra.hpp"
#include "NOX_Epetra_ModelEvaluatorInterface.H"
#include "NOX_Epetra_LinearSystem_Stratimikos.H"
#include "NOX_Epetra_LinearSystem_MPBD.hpp"

Piro::Epetra::StokhosMPSolver::
StokhosMPSolver(const Teuchos::RCP<Teuchos::ParameterList>& piroParams_,
		const Teuchos::RCP<Teuchos::ParameterList>& mpParams_,
		const Teuchos::RCP<const Epetra_Comm>& globalComm,
		int block_size, int num_spatial_procs) :
  piroParams(piroParams_),
  mpParams(mpParams_),
  num_mp(block_size)
{
  product_comm =
    Stokhos::buildMultiComm(*globalComm, block_size, num_spatial_procs);
}

Piro::Epetra::StokhosMPSolver::~StokhosMPSolver()
{
}

void
Piro::Epetra::StokhosMPSolver::
setup(const Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
      const Teuchos::RCP<NOX::Epetra::Observer>& noxObserver)
{
  Teuchos::RCP<const Epetra_Comm> mp_comm = 
    Stokhos::getStochasticComm(product_comm);
  Teuchos::RCP<const Epetra_Map> mp_block_map = 
    Teuchos::rcp(new Epetra_Map(num_mp, 0, *mp_comm));
  mp_model = model;

  // Turn mp_model into an MP-nonlinear problem
  mp_nonlin_model =
    Teuchos::rcp(new Stokhos::MPModelEvaluator(mp_model, product_comm,
					       mp_block_map, mpParams));

  bool use_mpbd_solver = mpParams->get("Use MPBD Solver", false);
  Teuchos::RCP<NOX::Epetra::LinearSystem> linsys;
  Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> nox_interface;
  if (use_mpbd_solver) {
    nox_interface = 
      Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(mp_nonlin_model));
    Teuchos::RCP<Epetra_Operator> A = 
      mp_nonlin_model->create_W();
    Teuchos::RCP<Epetra_Operator> M = 
      mp_nonlin_model->create_WPrec()->PrecOp;
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = 
      nox_interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = 
      nox_interface;
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = 
      nox_interface;

    Teuchos::ParameterList& noxParams = piroParams->sublist("NOX");
    Teuchos::ParameterList& printParams = noxParams.sublist("Printing");
    Teuchos::ParameterList& newtonParams = 
      noxParams.sublist("Direction").sublist("Newton");
    Teuchos::ParameterList& noxstratlsParams = 
      newtonParams.sublist("Stratimikos Linear Solver");
    Teuchos::ParameterList& mpbdParams = 
      mpParams->sublist("MPBD Linear Solver");
    mpbdParams.sublist("Deterministic Solver Parameters") = 
      noxstratlsParams;
    Teuchos::RCP<Epetra_Operator> inner_A = model->create_W();
    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> inner_nox_interface = 
      Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(model));
    Teuchos::RCP<NOX::Epetra::Interface::Required> inner_iReq = 
      inner_nox_interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> inner_iJac = 
      inner_nox_interface;
    Teuchos::RCP<const Epetra_Vector> inner_u = model->get_x_init();
    Teuchos::RCP<NOX::Epetra::LinearSystem> inner_linsys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
		     printParams, 
		     noxstratlsParams,
		     inner_iJac, inner_A, *inner_u));
    linsys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemMPBD(printParams, 
						     mpbdParams,
						     inner_linsys,
						     iReq, iJac, A,
						     model->get_x_map()));
  }

  Piro::Epetra::SolverFactory solverFactory;
  {
    const std::string token = "Interface";
    solverFactory.setNOXInterfaceProvider(token, nox_interface);
    Teuchos::sublist(piroParams, token)->set("Type", token);
  }
  {
    const std::string token = "Linear System";
    solverFactory.setNOXLinearSystemProvider(token, linsys);
    Teuchos::sublist(piroParams, token)->set("Type", token);
  }

  // Create solver to map p -> g
  mp_solver = solverFactory.createSolver(piroParams, mp_nonlin_model);

  // Create MP inverse model evaluator to map p_mp -> g_mp
  Teuchos::Array<int> mp_p_index_map = 
    mp_nonlin_model->get_p_mp_map_indices();
  Teuchos::Array<int> mp_g_index_map = 
    mp_nonlin_model->get_g_mp_map_indices();
  Teuchos::Array< Teuchos::RCP<const Epetra_Map> > base_g_maps = 
    mp_nonlin_model->get_g_mp_base_maps();
  mp_g_index_map.push_back(base_g_maps.size());
  base_g_maps.push_back(model->get_x_map());
  mp_inverse_solver =
    Teuchos::rcp(new Stokhos::MPInverseModelEvaluator(mp_solver,
						      mp_p_index_map,
						      mp_g_index_map,
						      base_g_maps));
}

Teuchos::RCP<const Epetra_Comm>
Piro::Epetra::StokhosMPSolver::
getSpatialComm() const
{
  return Stokhos::getSpatialComm(product_comm);
}

Teuchos::RCP<const Epetra_Comm>
Piro::Epetra::StokhosMPSolver::
getStochasticComm() const
{
  return Stokhos::getStochasticComm(product_comm);
}

Teuchos::RCP<const EpetraExt::MultiComm>
Piro::Epetra::StokhosMPSolver::
getGlobalMultiComm() const
{
  return product_comm;
}

Teuchos::RCP<const Epetra_Map> 
Piro::Epetra::StokhosMPSolver::get_x_map() const
{
  return mp_inverse_solver->get_x_map();
}

Teuchos::RCP<const Epetra_Map> 
Piro::Epetra::StokhosMPSolver::get_f_map() const
{
  return mp_inverse_solver->get_f_map();
}

Teuchos::RCP<const Epetra_Map> 
Piro::Epetra::StokhosMPSolver::get_p_map(int l) const
{
  return mp_inverse_solver->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> 
Piro::Epetra::StokhosMPSolver::get_g_map(int j) const
{
  return mp_inverse_solver->get_g_map(j);
}

Teuchos::RCP<const Epetra_Vector> 
Piro::Epetra::StokhosMPSolver::get_x_init() const
{
  return mp_inverse_solver->get_x_init();
}

Teuchos::RCP<const Epetra_Vector> 
Piro::Epetra::StokhosMPSolver::get_p_init(int l) const
{
  return mp_nonlin_model->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs 
Piro::Epetra::StokhosMPSolver::createInArgs() const
{
  return mp_inverse_solver->createInArgs();
}

EpetraExt::ModelEvaluator::OutArgs 
Piro::Epetra::StokhosMPSolver::createOutArgs() const
{
  return mp_inverse_solver->createOutArgs();
}

void 
Piro::Epetra::StokhosMPSolver::evalModel(const InArgs& inArgs,
			     const OutArgs& outArgs ) const
{
  mp_inverse_solver->evalModel(inArgs, outArgs);
}

Teuchos::RCP<Stokhos::ProductEpetraVector> 
Piro::Epetra::StokhosMPSolver::create_g_mp(int l, Epetra_DataAccess CV, 
					 const Epetra_Vector* v) const 
{
  OutArgs outargs = mp_nonlin_model->createOutArgs();
  int ng = outargs.Ng();
  if (piroParams->get<std::string>("Solver Type") == "NOX" && l == ng) {
    return mp_nonlin_model->create_x_mp(CV, v);
  }
  else
    return mp_nonlin_model->create_g_mp(l, CV, v);
}

Teuchos::RCP<Stokhos::ProductEpetraMultiVector> 
Piro::Epetra::StokhosMPSolver::create_g_mv_mp(int l, int num_vecs, 
					    Epetra_DataAccess CV, 
					    const Epetra_MultiVector* v) const
{
  OutArgs outargs = mp_nonlin_model->createOutArgs();
  int ng = outargs.Ng();
  if (piroParams->get<std::string>("Solver Type") == "NOX" && l == ng) {
    return mp_nonlin_model->create_x_mv_mp(num_vecs, CV, v);
  }
  else
    return mp_nonlin_model->create_g_mv_mp(l, num_vecs, CV, v);
}

