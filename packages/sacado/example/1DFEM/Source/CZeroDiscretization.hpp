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

#ifndef CZERODISCRETIZATION_HPP
#define CZERODISCRETIZATION_HPP

#include <vector>

#include "Teuchos_RefCountPtr.hpp"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"

#include "Mesh.hpp"

class CZeroDiscretization {
public:

  //! Constructor
  CZeroDiscretization(const std::vector<double>& coords,
		      unsigned int num_equations,
		      Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm);

  //! Destructor
  virtual ~CZeroDiscretization();

  //! Create element mesh
  virtual void createMesh();

  //! Create DOF maps
  virtual void createMaps();

  //! Create Jacobian graph
  virtual void createJacobianGraphs();

  //! Get element mesh
  virtual Teuchos::RefCountPtr<Mesh> getMesh(); 

  //! Get DOF map
  virtual Teuchos::RefCountPtr<Epetra_Map> getMap();

  //! Get overlapped DOF map
  virtual Teuchos::RefCountPtr<Epetra_Map> getOverlapMap();

  //! Get Jacobian graph
  virtual Teuchos::RefCountPtr<Epetra_CrsGraph> getJacobianGraph();

  //! Get overlap Jacobian graph
  virtual Teuchos::RefCountPtr<Epetra_CrsGraph> getOverlapJacobianGraph();


private:

  //! Private to prohibit copying
  CZeroDiscretization(const CZeroDiscretization&);

  //! Private to prohibit copying
  CZeroDiscretization& operator=(const CZeroDiscretization&);

protected:

  //! Total number of elements
  unsigned int n_elem;

  //! Coordinates of mesh nodes
  std::vector<double> x;

  //! Epetra communicator
  Teuchos::RefCountPtr<const Epetra_Comm> comm;

  //! Element mesh
  Teuchos::RefCountPtr<Mesh> mesh;

  //! Map
  Teuchos::RefCountPtr<Epetra_Map> map;

  //! Overlapped map
  Teuchos::RefCountPtr<Epetra_Map> overlap_map;

  //! Jacobian matrix graph
  Teuchos::RefCountPtr<Epetra_CrsGraph> graph;

  //! Overlapped Jacobian matrix graph
  Teuchos::RefCountPtr<Epetra_CrsGraph> overlap_graph;

  //! Number of processors
  unsigned int num_proc;

  //! Processor ID
  unsigned int myPID;

  //! Number of elements per processor
  unsigned int elem_per_proc;

  //! Number of elements on this processor
  unsigned int numMyElements;

  //! Number of nodes per element
  unsigned int nodes_per_element;

  //! Number of equations per node
  unsigned int neq;

  //! Global ID of first node
  unsigned int node_GID_base;

};

#endif // CZERODISCRETIZATION_HPP
