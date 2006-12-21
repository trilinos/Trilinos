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
#include "Epetra_Export.h"

#include "FEApp_CZeroDiscretization.hpp"

FEApp::CZeroDiscretization::CZeroDiscretization(
		  const std::vector<double>& coords,
		  unsigned int num_equations,
		  const Teuchos::RefCountPtr<const Epetra_Comm>& epetra_comm,
		  const Teuchos::RefCountPtr<Teuchos::ParameterList>& params) :
  x(coords),
  comm(epetra_comm),
  elemFactory(Teuchos::rcp(&(params->sublist("Element")),false)),
  mesh(),
  elem_map(),
  map(),
  overlap_map(),
  graph(),
  overlap_graph(),
  myPID(comm->MyPID()),
  numMyElements(0),
  nodes_per_element(0),
  neq(num_equations)
{
  // Distribute the elements equally among processors
  elem_map = Teuchos::rcp(new Epetra_Map(x.size()-1, 0, *comm));
  numMyElements = elem_map->NumMyElements();

  Teuchos::RefCountPtr<AbstractElement> base_element = 
    elemFactory.create();
  nodes_per_element = base_element->numNodes();
}

FEApp::CZeroDiscretization::~CZeroDiscretization()
{
}

void
FEApp::CZeroDiscretization::createMesh()
{
  mesh = Teuchos::rcp(new FEApp::Mesh);

  // Create elements and node IDs
  Teuchos::RefCountPtr<FEApp::AbstractElement> e;
  unsigned int elem_GID;
  for (unsigned int i=0; i<numMyElements; i++) {
    elem_GID = elem_map->GID(i);
    e = elemFactory.create();
    e->createNodes(x[elem_GID], x[elem_GID+1], elem_GID*(nodes_per_element-1));
    mesh->addElement(e);
  }
}

void
FEApp::CZeroDiscretization::createMaps()
{
  // Create overlap DOF map
  unsigned int overlapNumMyNodes = numMyElements*(nodes_per_element-1) + 1;
  unsigned int overlapNumMyDOF = overlapNumMyNodes*neq;
  unsigned int overlap_dof_GID_base = 
    elem_map->MinMyGID()*(nodes_per_element-1)*neq;
  std::vector<int> overlap_dof_GID(overlapNumMyDOF);
  for (unsigned int i=0; i<overlapNumMyDOF; i++)
    overlap_dof_GID[i] = overlap_dof_GID_base + i;
  overlap_map = 
    Teuchos::rcp(new Epetra_Map(-1, overlapNumMyDOF, &(overlap_dof_GID[0]), 
				0, *comm));
  
  // Create non-overlap DOF map
  if (myPID == 0) {
    int numMyDOF = overlapNumMyNodes*neq;
    map = 
      Teuchos::rcp(new Epetra_Map(-1, numMyDOF, &(overlap_dof_GID[0]), 
				  0, *comm));
  }
  else {
    int numMyDOF = (overlapNumMyNodes - 1)*neq;
    map = 
      Teuchos::rcp(new Epetra_Map(-1, numMyDOF, &(overlap_dof_GID[neq]), 
				  0, *comm));
  }
}

void
FEApp::CZeroDiscretization::createJacobianGraphs()
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
  Teuchos::RefCountPtr<FEApp::AbstractElement> e;
  for (FEApp::Mesh::iterator eit=mesh->begin(); eit!=mesh->end(); ++eit) {
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

	  } // column equations

	} // column nodes

      } // row equations
      
    } // row node
    
  } // element

  overlap_graph->FillComplete();

  Epetra_Export exporter(*overlap_map, *map);
  graph->Export(*overlap_graph, exporter, Insert);
  graph->FillComplete();
}
	    
Teuchos::RefCountPtr<const FEApp::Mesh>
FEApp::CZeroDiscretization::getMesh() const
{
  return mesh;
}

Teuchos::RefCountPtr<const Epetra_Map>
FEApp::CZeroDiscretization::getMap() const
{
  return map;
}

Teuchos::RefCountPtr<const Epetra_Map>
FEApp::CZeroDiscretization::getOverlapMap() const
{
  return overlap_map;
}

Teuchos::RefCountPtr<const Epetra_CrsGraph>
FEApp::CZeroDiscretization::getJacobianGraph() const
{
  return graph;
}

Teuchos::RefCountPtr<const Epetra_CrsGraph>
FEApp::CZeroDiscretization::getOverlapJacobianGraph() const
{
  return overlap_graph;
}

int
FEApp::CZeroDiscretization::getNumNodesPerElement() const
{
  return nodes_per_element;
}
