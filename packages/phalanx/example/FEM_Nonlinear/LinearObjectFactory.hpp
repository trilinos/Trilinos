// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
