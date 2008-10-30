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

#include "LinearObjectFactoryVBR.hpp"
#include "Element_Linear2D.hpp"
#include "Teuchos_TestForException.hpp"
#include "Epetra_Export.h"

//**********************************************************************
LinearObjectFactoryVBR::LinearObjectFactoryVBR(const MeshBuilder& mb, 
					 const Teuchos::RCP<Epetra_Comm>& comm,
					 int number_of_equations_per_node)
  :
  m_num_eq(number_of_equations_per_node)
{
  std::vector<Element_Linear2D>& cells = *(mb.myElements());
  std::vector<Element_Linear2D>::iterator cell;
  
  // Maps
  {
    std::vector<int> owned_unknown_gids(0);
    std::vector<int> overlapped_unknown_gids(0);


    // Map
    for (cell = cells.begin(); cell != cells.end(); ++cell) {
      for (std::size_t node = 0; node < cell->numNodes(); ++node) {
	  
	int index = cell->globalNodeId(node);

	if (cell->ownsNode(node)) {
	  if ( std::find(owned_unknown_gids.begin(),
			 owned_unknown_gids.end(),
			 index) 
	       == owned_unknown_gids.end() ) {
	    owned_unknown_gids.push_back(index);
	  }
	}

	if ( std::find(overlapped_unknown_gids.begin(),
		       overlapped_unknown_gids.end(),
		       index) 
	     == overlapped_unknown_gids.end() )
	  overlapped_unknown_gids.push_back(index);
	
      }
    }
    
    m_owned_map = 
      Teuchos::rcp(new Epetra_BlockMap(-1, owned_unknown_gids.size(),
				       &owned_unknown_gids[0], m_num_eq, 
				       0, *comm));
    m_overlapped_map = 
      Teuchos::rcp(new Epetra_BlockMap(-1, overlapped_unknown_gids.size(),
				       &overlapped_unknown_gids[0], m_num_eq,
				       0, *comm));
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

      for (std::size_t row = 0; row < cell->numNodes(); ++row) {
	for (std::size_t col = 0; col < cell->numNodes(); ++col) {
	  int global_row = cell->globalNodeId(row);
	  int global_col = cell->globalNodeId(col);
	  m_overlapped_graph->InsertGlobalIndices(global_row, 1, &global_col);
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
Teuchos::RCP<const Epetra_BlockMap> LinearObjectFactoryVBR::ownedMap() const
{
  return m_owned_map;
}

//**********************************************************************
Teuchos::RCP<const Epetra_BlockMap> LinearObjectFactoryVBR::overlappedMap() const
{
  return m_overlapped_map;
}

//**********************************************************************
Teuchos::RCP<const Epetra_CrsGraph> LinearObjectFactoryVBR::ownedGraph() const
{
  return m_owned_graph;
}

//**********************************************************************
Teuchos::RCP<const Epetra_CrsGraph> 
LinearObjectFactoryVBR::overlappedGraph() const
{
  return m_overlapped_graph;
}

//**********************************************************************
void LinearObjectFactoryVBR::print(std::ostream& os) const
{
  os << "Inside LinearObjectFactoryVBR!" << std::endl;
}

//**********************************************************************
std::ostream& operator<<(std::ostream& os, const LinearObjectFactoryVBR& f)
{
  f.print(os);
  return os;
}

//**********************************************************************
