// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include <list>
#include <vector>
#include <algorithm>
#include "LinearObjectFactory.hpp"
#include "Element_Linear2D.hpp"
#include "Teuchos_TestForException.hpp"
#include "Epetra_Export.h"

//**********************************************************************
LinearObjectFactory::LinearObjectFactory(const MeshBuilder& mb, 
					 const Teuchos::RCP<Epetra_Comm>& comm,
					 int number_of_equations_per_node)
  :
  m_num_eq(number_of_equations_per_node)
{
  std::vector<Element_Linear2D>& cells = *(mb.myElements());
  std::vector<Element_Linear2D>::iterator cell;
  
  // Maps
  {
    std::list<int> list_owned_unknown_gids(0);
    std::list<int> list_overlapped_unknown_gids(0);


    // Map
    for (cell = cells.begin(); cell != cells.end(); ++cell) {
      for (int node = 0; node < cell->numNodes(); ++node) {
	for (int eq = 0; eq < m_num_eq; ++eq) {
	  
	  int index = cell->globalNodeId(node) * m_num_eq + eq;

	  if (cell->ownsNode(node))
	    list_owned_unknown_gids.push_back(index);

	  list_overlapped_unknown_gids.push_back(index);

	}	
      }
    }

    list_owned_unknown_gids.sort();
    list_owned_unknown_gids.unique();
    list_overlapped_unknown_gids.sort();
    list_overlapped_unknown_gids.unique();

    // Copy lists into vectors (need contiguous arrays)
    std::vector<int> owned_unknown_gids(list_owned_unknown_gids.size());
    std::vector<int> overlapped_unknown_gids(list_overlapped_unknown_gids.size());

    std::size_t index = 0;
    for (std::list<int>::iterator node = list_owned_unknown_gids.begin();
	 node != list_owned_unknown_gids.end(); ++node, ++index)
      owned_unknown_gids[index] = *node;
    index = 0;
    for (std::list<int>::iterator node = list_overlapped_unknown_gids.begin();
	 node != list_overlapped_unknown_gids.end(); ++node, ++index)
      overlapped_unknown_gids[index] = *node;


    m_owned_map = Teuchos::rcp(new Epetra_Map(-1, owned_unknown_gids.size(),
					      &owned_unknown_gids[0], 0, 
					      *comm));
    m_overlapped_map = 
      Teuchos::rcp(new Epetra_Map(-1, overlapped_unknown_gids.size(),
				  &overlapped_unknown_gids[0], 0, *comm));
  }


  
  // Graph
  {
     
    Epetra_DataAccess copy = ::Copy;
    const std::size_t approximate_indices_per_row = 9 * m_num_eq;

    m_owned_graph = 
      Teuchos::rcp(new Epetra_CrsGraph(copy, *m_owned_map, 
				       approximate_indices_per_row));

    m_overlapped_graph = 
      Teuchos::rcp(new Epetra_CrsGraph(copy, *m_overlapped_map, 
				       approximate_indices_per_row));

    std::vector<Element_Linear2D>::iterator cell = cells.begin();
    for (cell = cells.begin(); cell != cells.end(); ++cell) {
            
      std::vector<int> col_indices(0);
      for (int node = 0; node < cell->numNodes(); ++node)
	for (int eq = 0; eq < m_num_eq; ++eq)
	  col_indices.push_back(cell->globalNodeId(node) * m_num_eq + eq);
      
      for (int node = 0; node < cell->numNodes(); ++node) {
	for (int eq = 0; eq < m_num_eq; ++eq) {

	  int global_row = cell->globalNodeId(node) * m_num_eq + eq;

	  m_overlapped_graph->InsertGlobalIndices(global_row, 
						  col_indices.size(), 
						  &col_indices[0]);



	}
      }
    }
    m_overlapped_graph->FillComplete();
    
    // Export overlapped entires to assemble the owned graph
    Epetra_Export exporter(*m_overlapped_map, *m_owned_map);
    m_owned_graph->Export(*m_overlapped_graph, exporter, Insert);
    m_owned_graph->FillComplete();

  }

}

//**********************************************************************
Teuchos::RCP<const Epetra_Map> LinearObjectFactory::ownedMap() const
{
  return m_owned_map;
}

//**********************************************************************
Teuchos::RCP<const Epetra_Map> LinearObjectFactory::overlappedMap() const
{
  return m_overlapped_map;
}

//**********************************************************************
Teuchos::RCP<const Epetra_CrsGraph> LinearObjectFactory::ownedGraph() const
{
  return m_owned_graph;
}

//**********************************************************************
Teuchos::RCP<const Epetra_CrsGraph> 
LinearObjectFactory::overlappedGraph() const
{
  return m_overlapped_graph;
}

//**********************************************************************
void LinearObjectFactory::print(std::ostream& os) const
{
  os << "Inside LinearObjectFactory!" << std::endl;
}

//**********************************************************************
std::ostream& operator<<(std::ostream& os, const LinearObjectFactory& f)
{
  f.print(os);
  return os;
}

//**********************************************************************
