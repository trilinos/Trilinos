// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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

#include "FEApp_BrusselatorProblem.hpp"
#include "FEApp_ConstantDirichletBC.hpp"

FEApp::BrusselatorProblem::BrusselatorProblem(
		   const Teuchos::RefCountPtr<Teuchos::ParameterList>& params)
{
  alpha = params->get("alpha", 1.0);
  beta = params->get("beta", 1.0);
  D1 = params->get("D1", 1.0);
  D2 = params->get("D2", 1.0);
}

FEApp::BrusselatorProblem::~BrusselatorProblem()
{
}

unsigned int
FEApp::BrusselatorProblem::numEquations() const
{
  return 2;
}

void
FEApp::BrusselatorProblem:: buildPDEs(
		       FEApp::AbstractPDE_TemplateManager<ValidTypes>& pdeTM)
{
  FEApp::BrusselatorPDE_TemplateBuilder pdeBuilder(alpha, beta, D1, D2);
  pdeTM.buildObjects(pdeBuilder);
}

std::vector< Teuchos::RefCountPtr<const FEApp::AbstractBC> >
FEApp::BrusselatorProblem::buildBCs(const Epetra_Map& dofMap)
{
  std::vector< Teuchos::RefCountPtr<const FEApp::AbstractBC> > bc(4);
  bc[0] = Teuchos::rcp(new FEApp::ConstantDirichletBC(dofMap.MinAllGID(),
						      alpha));
  bc[1] = Teuchos::rcp(new FEApp::ConstantDirichletBC(dofMap.MinAllGID()+1,
						      beta/alpha));
  bc[2] = Teuchos::rcp(new FEApp::ConstantDirichletBC(dofMap.MaxAllGID()-1,
						      alpha));
  bc[3] = Teuchos::rcp(new FEApp::ConstantDirichletBC(dofMap.MaxAllGID(),
						      beta/alpha));

  return bc;
}

Teuchos::RefCountPtr<Epetra_Vector>
FEApp::BrusselatorProblem::buildInitialSolution(const Epetra_Map& dofMap)
{
  Teuchos::RefCountPtr<Epetra_Vector> u =
    Teuchos::rcp(new Epetra_Vector(dofMap, false));
//   for (int i=0; i<u->MyLength(); i++) {
//     (*u)[2*i]   = alpha;
//     (*u)[2*i+1] = beta/alpha;
//   }
  u->PutScalar(0.5);

  return u;
}
