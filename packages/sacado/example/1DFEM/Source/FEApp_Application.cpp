// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "FEApp_Application.hpp"
#include "FEApp_ProblemFactory.hpp"
#include "FEApp_QuadratureFactory.hpp"
#include "FEApp_DiscretizationFactory.hpp"
#include "FEApp_InitPostOps.hpp"
#include "FEApp_GlobalFill.hpp"

FEApp::Application::Application(
		   const std::vector<double>& coords,
		   const Teuchos::RefCountPtr<const Epetra_Comm>& comm,
		   const Teuchos::RefCountPtr<Teuchos::ParameterList>& params) 
{
  // Create problem object
  Teuchos::RefCountPtr<Teuchos::ParameterList> problemParams = 
    Teuchos::rcp(&(params->sublist("Problem")),false);
  FEApp::ProblemFactory problemFactory(problemParams);
  Teuchos::RefCountPtr<FEApp::AbstractProblem> problem = 
    problemFactory.create();

  // Get number of equations
  unsigned int num_equations = problem->numEquations();

  // Create quadrature object
  Teuchos::RefCountPtr<Teuchos::ParameterList> quadParams = 
    Teuchos::rcp(&(params->sublist("Quadrature")),false);
  FEApp::QuadratureFactory quadFactory(quadParams);
  quad = quadFactory.create();

  // Create discretization object
  Teuchos::RefCountPtr<Teuchos::ParameterList> discParams = 
    Teuchos::rcp(&(params->sublist("Discretization")),false);
  FEApp::DiscretizationFactory discFactory(discParams);
  disc = discFactory.create(coords, num_equations, comm);
  disc->createMesh();
  disc->createMaps();
  disc->createJacobianGraphs();

  // Create Epetra objects
  importer = Teuchos::rcp(new Epetra_Import(*(disc->getOverlapMap()), 
					    *(disc->getMap())));
  exporter = Teuchos::rcp(new Epetra_Export(*(disc->getOverlapMap()), 
					    *(disc->getMap())));
  overlapped_x = Teuchos::rcp(new Epetra_Vector(*(disc->getOverlapMap())));
  overlapped_f = Teuchos::rcp(new Epetra_Vector(*(disc->getOverlapMap())));
  overlapped_jac = 
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, 
				      *(disc->getOverlapJacobianGraph())));

  // Initialize PDEs
  problem->buildPDEs(pdeTM);
  typedef FEApp::AbstractPDE_TemplateManager<ValidTypes>::iterator iterator;
  int nqp = quad->numPoints();
  int nn = disc->getNumNodesPerElement();
  for (iterator it = pdeTM.begin(); it != pdeTM.end(); ++it)
    it->init(nqp, nn);

  // Initialize BCs
  bc = problem->buildBCs(*(disc->getMap()));

  // Intialize initial solution
  initial_x = problem->buildInitialSolution(*(disc->getMap()));
}

FEApp::Application::~Application()
{
}

Teuchos::RefCountPtr<const Epetra_Map>
FEApp::Application::getMap() const
{
  return disc->getMap();
}

Teuchos::RefCountPtr<const Epetra_CrsGraph>
FEApp::Application::getJacobianGraph() const
{
  return disc->getJacobianGraph();
}

Teuchos::RefCountPtr<const Epetra_Vector>
FEApp::Application::getInitialSolution() const
{
  return initial_x;
}

void
FEApp::Application::computeGlobalResidual(const Epetra_Vector& x,
					  Epetra_Vector& f)
{
  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Zero out overlapped residual
  overlapped_f->PutScalar(0.0);
  f.PutScalar(0.0);

  // Create residual init/post op
  FEApp::ResidualOp op(overlapped_x, overlapped_f);

  // Get template PDE instantiation
  Teuchos::RefCountPtr< FEApp::AbstractPDE<ResidualOp::fill_type> > pde = 
    Teuchos::rcp(&(pdeTM.getAsObject<ResidualOp::fill_type>()),false);

  // Do global fill
  FEApp::GlobalFill<ResidualOp::fill_type> globalFill(disc->getMesh(), quad, 
						      pde);
  globalFill.computeGlobalFill(op);

  // Assemble global residual
  f.Export(*overlapped_f, *exporter, Add);

  // Apply boundary conditions
  for (unsigned int i=0; i<bc.size(); i++)
    bc[i]->applyResidual(x, f);
}

void
FEApp::Application::computeGlobalJacobian(
			 const Epetra_Vector& x,
			 Epetra_Vector& f,
			 Epetra_CrsMatrix& jac)
{
  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Zero out overlapped residual
  overlapped_f->PutScalar(0.0);
  f.PutScalar(0.0);

  // Zero out Jacobian
  overlapped_jac->PutScalar(0.0);
  jac.PutScalar(0.0);

  // Create Jacobian init/post op
  FEApp::JacobianOp op(overlapped_x, overlapped_f, overlapped_jac);

  // Get template PDE instantiation
  Teuchos::RefCountPtr< FEApp::AbstractPDE<JacobianOp::fill_type> > pde = 
    Teuchos::rcp(&(pdeTM.getAsObject<JacobianOp::fill_type>()),false);

  // Do global fill
  FEApp::GlobalFill<JacobianOp::fill_type> globalFill(disc->getMesh(), 
						      quad, pde);
  globalFill.computeGlobalFill(op);

  // Assemble global residual
  f.Export(*overlapped_f, *exporter, Add);

  // Assemble global Jacobian
  jac.Export(*overlapped_jac, *exporter, Add);

  // Apply boundary conditions
  for (unsigned int i=0; i<bc.size(); i++) {
    bc[i]->applyResidual(x, f);
    bc[i]->applyJacobian(x, f, jac);
  }
}
