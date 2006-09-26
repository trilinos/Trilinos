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

template <typename ScalarT>
GlobalFill<ScalarT>::GlobalFill(
	  const Teuchos::RefCountPtr<const Mesh>& elementMesh,
	  const Teuchos::RefCountPtr<const AbstractQuadrature>& quadRule,
	  const Teuchos::RefCountPtr< AbstractPDE<ScalarT> >& pdeEquations):
  mesh(elementMesh),
  quad(quadRule),
  pde(pdeEquations),
  nnode(0),
  neqn(pde->numEquations()),
  ndof(0),
  elem_x(),
  elem_f()
{
  Teuchos::RefCountPtr<const AbstractElement> e0 = *(mesh->begin());
  nnode = e0->numNodes();
  ndof = nnode*neqn;
  elem_x.resize(ndof);
  elem_f.resize(ndof);
}

template <typename ScalarT>
GlobalFill<ScalarT>::~GlobalFill()
{
}

template <typename ScalarT>
void
GlobalFill<ScalarT>::computeGlobalFill(AbstractInitPostOp<ScalarT>& initPostOp)
{
  // Loop over elements
  Teuchos::RefCountPtr<const AbstractElement> e;
  for (Mesh::const_iterator eit=mesh->begin(); eit!=mesh->end(); ++eit) {
    e = *eit;

    // Zero out element residual
    for (unsigned int i=0; i<ndof; i++)
      elem_f[i] = 0.0;

    // Initialize element solution
    initPostOp.evalInit(*e, neqn, elem_x);

    // Compute element residual
    pde->evaluateElementResidual(*quad, *e, elem_x, elem_f);

    // Post-process element residual
    initPostOp.evalPost(*e, neqn, elem_f);

  }

}
