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

#include "CZeroDiscretization.hpp"

#include "LinearElement.hpp"

CZeroDiscretization::CZeroDiscretization(
		      const std::vector<double>& coords,
		      unsigned int num_equations,
		      Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm) :
  n_elem(coords.size()-1),
  x(coords),
  comm(epetra_comm),
  mesh(),
  map(),
  overlap_map(),
  graph(),
  overlap_graph(),
  num_proc(comm->NumProc()),
  myPID(comm->MyPID()),
  elem_per_proc(n_elem / num_proc),
  numMyElements(elem_per_proc),
  nodes_per_element(0),
  neq(num_equations),
  node_GID_base(elem_per_proc*(nodes_per_element-1)*myPID)
{
  if (myPID == num_proc - 1)
    numMyElements = n_elem - elem_per_proc * (num_proc - 1);

  AbstractElement *base_element = new LinearElement();  // replace w/factory
  nodes_per_element = base_element->numNodes();
  delete base_element;
}

CZeroDiscretization::~CZeroDiscretization()
{
}

void
CZeroDiscretization::createMesh()
{
  mesh = Teuchos::rcp(new Mesh);

  // Create elements and node IDs
  Teuchos::RefCountPtr<AbstractElement> e;
  unsigned int elem_GID;
  for (unsigned int i=0; i<numMyElements; i++) {
    elem_GID = elem_per_proc*myPID + i;
    e = Teuchos::rcp(new LinearElement); // replace w/factory
    e->createNodes(x[elem_GID], x[elem_GID+1], 
		   node_GID_base + i*(nodes_per_element-1));
    mesh->addElement(e);
  }
}

void
CZeroDiscretization::createMaps()
{
  // Create overlap DOF map
  unsigned int overlapNumMyNodes = numMyElements*(nodes_per_element-1) + 1;
  unsigned int overlapNumMyDOF = overlapNumMyNodes*neq;
  unsigned int overlap_dof_GID_base = node_GID_base*neq;
  std::vector<int> overlap_dof_GID(overlapNumMyDOF);
  for (unsigned int i=0; i<overlapNumMyDOF; i++)
    overlap_dof_GID[i] = overlap_dof_GID_base + i;
  overlap_map = 
    Teuchos::rcp(new Epetra_Map(-1, overlapNumMyDOF, &(overlap_dof_GID[0]), 
				0, *comm));
  
  // Create non-overlap DOF map
  unsigned int numMyNodes = overlapNumMyNodes - 1;
  if (myPID == 0)
    numMyNodes = overlapNumMyNodes;
  unsigned int numMyDOF = numMyNodes*neq;
  unsigned int dof_GID_base;
  if (myPID == 0)
    dof_GID_base = overlap_dof_GID_base;
  else
    dof_GID_base = overlap_dof_GID_base + 1;
  std::vector<int> dof_GID(numMyDOF);
  for (unsigned int i=0; i<numMyDOF; i++)
    dof_GID[i] = dof_GID_base + i;
  map = 
    Teuchos::rcp(new Epetra_Map(-1, numMyDOF, &(dof_GID[0]), 0, *comm));
}

void
CZeroDiscretization::createJacobianGraphs()
{
  // Generate matrix graphs
  graph = 
    Teuchos::rcp(new Epetra_CrsGraph(Copy, *map, neq*nodes_per_element, 
				     false));
  overlap_graph = 
    Teuchos::rcp(new Epetra_CrsGraph(Copy, *overlap_map, 
				     neq*nodes_per_element, false));
  int row, col;
  
  // Loop over elements
  Teuchos::RefCountPtr<AbstractElement> e;
  for (Mesh::iterator eit=mesh->begin(); eit!=mesh->end(); ++eit) {
    e = *eit;

    // Loop over nodes in element
    for (unsigned int node_row=0; node_row<nodes_per_element; node_row++) {

      // Loop over equations per node
      for (unsigned int eq_row=0; eq_row<neq; eq_row++) {

	// Matrix row
	row = static_cast<int>(e->nodeGID(node_row)*neq + eq_row);

	// Loop over nodes in element
	for (unsigned int node_col=0; node_col<nodes_per_element; 
	     node_col++){
	    
	  // Loop over equations per node
	  for (unsigned int eq_col=0; eq_col<neq; eq_col++) {
	      
	    // Matrix column
	    col = static_cast<int>(e->nodeGID(node_col)*neq + eq_col);

	    // Add column indices
	    overlap_graph->InsertGlobalIndices(row, 1, &col);
	    if (map->MyGID(row))
	      graph->InsertGlobalIndices(row, 1, &col);

	  } // column equations

	} // column nodes

      } // row equations
      
    } // row node
    
  } // element

  graph->FillComplete();
  overlap_graph->FillComplete();
}
	    
Teuchos::RefCountPtr<Mesh>
CZeroDiscretization::getMesh() 
{
  return mesh;
}

Teuchos::RefCountPtr<Epetra_Map>
CZeroDiscretization::getMap()
{
  return map;
}

Teuchos::RefCountPtr<Epetra_Map>
CZeroDiscretization::getOverlapMap()
{
  return overlap_map;
}

Teuchos::RefCountPtr<Epetra_CrsGraph>
CZeroDiscretization::getJacobianGraph()
{
  return graph;
}

Teuchos::RefCountPtr<Epetra_CrsGraph>
CZeroDiscretization::getOverlapJacobianGraph()
{
  return overlap_graph;
}
