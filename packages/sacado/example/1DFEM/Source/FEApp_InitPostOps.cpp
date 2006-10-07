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

#include "FEApp_InitPostOps.hpp"

#include "Epetra_Map.h"

FEApp::ResidualOp::ResidualOp(
	      const Teuchos::RefCountPtr<const Epetra_Vector>& overlapped_x,
	      const Teuchos::RefCountPtr<Epetra_Vector>& overlapped_f) :
  x(overlapped_x),
  f(overlapped_f)
{
}

FEApp::ResidualOp::~ResidualOp()
{
}

void
FEApp::ResidualOp::evalInit(const FEApp::AbstractElement& e,
			    unsigned int neqn,
			    std::vector<double>& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Copy element solution
  for (unsigned int i=0; i<e.numNodes(); i++) {
    node_GID = e.nodeGID(i);
    firstDOF = x->Map().LID(node_GID*neqn);
    for (unsigned int j=0; j<neqn; j++) {
      elem_x[neqn*i+j] = (*x)[firstDOF+j];
    }
  }

//   // Copy element solution
//   int row;
//   unsigned int lrow;
//   for (unsigned int node_row=0; node_row<nnode; node_row++) {
//     for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
//       lrow = neqn*node_row+eq_row;
//       row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);
//       if (!x->Map().MyGID(row)) {
// 	std::cout << "ResidualOp::evalInit:  invalid row " << row 
// 		  << " for node " << node_row 
// 		  << "and equation " << eq_row << std::endl;
//       }
//       elem_x[lrow] = (*x)[x->Map().LID(row)];
//     }
//   }
}

void
FEApp::ResidualOp::evalPost(const FEApp::AbstractElement& e,
			    unsigned int neqn,
			    std::vector<double>& elem_f)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

//   // Sum element residual into global residual
//   for (unsigned int i=0; i<e.numNodes(); i++) {
//     node_GID = e.nodeGID(i);
//     firstDOF = f->Map().LID(node_GID*neqn);
//     for (unsigned int j=0; j<neqn; j++) {
//       (*f)[firstDOF+j] += elem_f[neqn*i+j];
//     }
//   }

  // Sum element residual into global residual
  int row;
  unsigned int lrow;
  for (unsigned int node_row=0; node_row<nnode; node_row++) {
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);
      f->SumIntoGlobalValues(1, &(elem_f[lrow]), &row);
    }
  }
}

FEApp::JacobianOp::JacobianOp(
	     const Teuchos::RefCountPtr<const Epetra_Vector>& overlapped_x,
	     const Teuchos::RefCountPtr<Epetra_Vector>& overlapped_f,
	     const Teuchos::RefCountPtr<Epetra_CrsMatrix>& overlapped_jac) :
  x(overlapped_x),
  f(overlapped_f),
  jac(overlapped_jac)
{
}

FEApp::JacobianOp::~JacobianOp()
{
}

void
FEApp::JacobianOp::evalInit(const FEApp::AbstractElement& e,
			    unsigned int neqn,
			    std::vector< Sacado::Fad::DFad<double> >& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Number of dof
  unsigned int ndof = nnode*neqn;

  // Copy element solution
  for (unsigned int i=0; i<nnode; i++) {
    node_GID = e.nodeGID(i);
    firstDOF = x->Map().LID(node_GID*neqn);
    for (unsigned int j=0; j<neqn; j++) {
      elem_x[neqn*i+j] = 
	Sacado::Fad::DFad<double>(ndof, neqn*i+j, (*x)[firstDOF+j]);
    }
  }

//   // Copy element solution
//   int row;
//   unsigned int lrow;
//   for (unsigned int node_row=0; node_row<nnode; node_row++) {
//     for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
//       lrow = neqn*node_row+eq_row;
//       row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);
//       if (!x->Map().MyGID(row)) {
// 	std::cout << "JacobianOp::evalInit:  invalid row " << row 
// 		  << " for node " << node_row 
// 		  << "and equation " << eq_row << std::endl;
//       }
//       elem_x[lrow] = Sacado::Fad::DFad<double>(ndof, lrow, 
// 					       (*x)[x->Map().LID(row)]);
//     }
//   }
}

void
FEApp::JacobianOp::evalPost(const FEApp::AbstractElement& e,
			    unsigned int neqn,
			    std::vector< Sacado::Fad::DFad<double> >& elem_f)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

   // Number of nodes
  unsigned int nnode = e.numNodes();

  // Sum element residual and Jacobian into global residual, Jacobian
  // Loop over nodes in element
  int row, col;
  unsigned int lrow, lcol;
  for (unsigned int node_row=0; node_row<nnode; node_row++) {

    // Loop over equations per node
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;

      // Global row
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);

      // Sum residual
      f->SumIntoGlobalValues(1, &(elem_f[lrow].val()), &row);
	
      // Check derivative array is nonzero
      if (elem_f[lrow].hasFastAccess()) {

	// Loop over nodes in element
	for (unsigned int node_col=0; node_col<nnode; node_col++){
	    
	  // Loop over equations per node
	  for (unsigned int eq_col=0; eq_col<neqn; eq_col++) {
	    lcol = neqn*node_col+eq_col;
	      
	    // Global column
	    col = static_cast<int>(e.nodeGID(node_col)*neqn + eq_col);

	    // Sum Jacobian
	    jac->SumIntoGlobalValues(row, 1, 
				     &(elem_f[lrow].fastAccessDx(lcol)),
				     &col);

	  } // column equations
	  
	} // column nodes
	
      } // has fast access
      
    } // row equations

  } // row node

}
  
