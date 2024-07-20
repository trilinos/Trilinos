// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_LINEAR_OBJECT_FACTORY_HPP
#define PHX_EXAMPLE_LINEAR_OBJECT_FACTORY_HPP

#include "Phalanx_config.hpp"
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
