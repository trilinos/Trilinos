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

#ifndef PHX_EXAMPLE_LINEAR_OBJECT_FACTORY_HPP
#define PHX_EXAMPLE_LINEAR_OBJECT_FACTORY_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "MeshBuilder.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"



/** \brief Builds linear solver objects given the number of unknowns per node

*/
class LinearObjectFactory {
  
public:

  LinearObjectFactory(const MeshBuilder& mb, 
		      const Teuchos::RCP<Epetra_Comm>& comm, 
		      int number_of_equations_per_node);

  Teuchos::RCP<const Epetra_Map> ownedMap() const;

  Teuchos::RCP<const Epetra_Map> overlappedMap() const;

  Teuchos::RCP<const Epetra_CrsGraph> ownedGraph() const;

  Teuchos::RCP<const Epetra_CrsGraph> overlappedGraph() const;

  void print(std::ostream& os) const;

private:
  
  //! Number of equations per node
  int m_num_eq;

  Teuchos::RCP<Epetra_Map> m_owned_map;

  Teuchos::RCP<Epetra_Map> m_overlapped_map;

  Teuchos::RCP<Epetra_CrsGraph> m_owned_graph;

  Teuchos::RCP<Epetra_CrsGraph> m_overlapped_graph;

};

std::ostream& operator<<(std::ostream& os, const LinearObjectFactory& b);

#endif
