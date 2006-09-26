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

#include "Application.hpp"
#include "InitPostOps.hpp"

#include "GlobalFill.hpp"

Application::Application(
	    const Teuchos::RefCountPtr<CZeroDiscretization> c0_disc,
	    const std::vector< Teuchos::RefCountPtr<const AbstractBC> > bcs,
	    const Teuchos::RefCountPtr<const AbstractQuadrature>& quadRule) :
  disc(c0_disc),
  bc(bcs),
  quad(quadRule),
  importer(*(disc->getOverlapMap()), *(disc->getMap())),
  exporter(*(disc->getOverlapMap()), *(disc->getMap())),
  overlapped_x(*(disc->getOverlapMap())),
  overlapped_f(*(disc->getOverlapMap())),
  overlapped_jac(Copy, *(disc->getOverlapJacobianGraph()))
{
}

Application::~Application()
{
}

void
Application::computeGlobalResidual(AbstractPDE<double>& pde,
				   const Epetra_Vector& x,
				   Epetra_Vector& f)
{
  // Scatter x to the overlapped distrbution
  overlapped_x.Import(x, importer, Insert);

  // Zero out overlapped residual
  overlapped_f.PutScalar(0.0);

  // Create residual init/post op
  ResidualOp op(Teuchos::rcp(&overlapped_x, false), 
		Teuchos::rcp(&overlapped_f, false));

  // Do global fill
  GlobalFill<ResidualOp::fill_type> globalFill(disc->getMesh(), quad, 
					       Teuchos::rcp(&pde, false));
  globalFill.computeGlobalFill(op);

  // Assemble global residual
  f.Export(overlapped_f, exporter, Add);

  // Apply boundary conditions
  for (unsigned int i=0; i<bc.size(); i++)
    bc[i]->applyResidual(x, f);
}

void
Application::computeGlobalJacobian(
			       AbstractPDE< Sacado::Fad::DFad<double> >& pde,
			       const Epetra_Vector& x,
			       Epetra_Vector& f,
			       Epetra_CrsMatrix& jac)
{
  // Scatter x to the overlapped distrbution
  overlapped_x.Import(x, importer, Insert);

  // Zero out overlapped residual
  overlapped_f.PutScalar(0.0);

  // Zero out Jacobian
  overlapped_jac.PutScalar(0.0);

  // Create Jacobian init/post op
  JacobianOp op(Teuchos::rcp(&overlapped_x, false), 
		Teuchos::rcp(&overlapped_f, false), 
		Teuchos::rcp(&overlapped_jac, false));

  // Do global fill
  GlobalFill<JacobianOp::fill_type> globalFill(disc->getMesh(), 
					       quad, 
					       Teuchos::rcp(&pde, false));
  globalFill.computeGlobalFill(op);

  // Assemble global residual
  f.Export(overlapped_f, exporter, Add);

  // Assemble global Jacobian
  jac.Export(overlapped_jac, exporter, Add);

  // Apply boundary conditions
  for (unsigned int i=0; i<bc.size(); i++) {
    bc[i]->applyResidual(x, f);
    bc[i]->applyJacobian(x, f, jac);
  }
}
