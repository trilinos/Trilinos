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
#include "FEApp_InitPostOps.hpp"
#include "FEApp_GlobalFill.hpp"

FEApp::Application::~Application()
{
}

Teuchos::RefCountPtr<Epetra_Map>
FEApp::Application::getMap()
{
  return disc->getMap();
}

Teuchos::RefCountPtr<Epetra_CrsGraph>
FEApp::Application::getJacobianGraph()
{
  return disc->getJacobianGraph();
}

void
FEApp::Application::setBCs(
       const std::vector< Teuchos::RefCountPtr<const FEApp::AbstractBC> > bcs)
{
  bc = bcs;
}

void
FEApp::Application::computeGlobalResidual(const Epetra_Vector& x,
					  Epetra_Vector& f)
{
  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Zero out overlapped residual
  overlapped_f->PutScalar(0.0);

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

  // Zero out Jacobian
  overlapped_jac->PutScalar(0.0);

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
