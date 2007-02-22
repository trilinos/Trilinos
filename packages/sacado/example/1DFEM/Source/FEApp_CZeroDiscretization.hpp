// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef FEAPP_CZERODISCRETIZATION_HPP
#define FEAPP_CZERODISCRETIZATION_HPP

#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Comm.h"

#include "FEApp_AbstractDiscretization.hpp"
#include "FEApp_ElementFactory.hpp"

namespace FEApp {

  class CZeroDiscretization : public FEApp::AbstractDiscretization {
  public:

    //! Constructor
    CZeroDiscretization(
		 const std::vector<double>& coords,
		 unsigned int num_equations,
		 const Teuchos::RefCountPtr<const Epetra_Comm>& epetra_comm,
		 const Teuchos::RefCountPtr<Teuchos::ParameterList>& params);

    //! Destructor
    virtual ~CZeroDiscretization();

    //! Create element mesh
    virtual void createMesh();

    //! Create DOF maps
    virtual void createMaps();

    //! Create Jacobian graph
    virtual void createJacobianGraphs();

    //! Get element mesh
    virtual Teuchos::RefCountPtr<const FEApp::Mesh> 
    getMesh() const; 

    //! Get DOF map
    virtual Teuchos::RefCountPtr<const Epetra_Map> 
    getMap() const;

    //! Get overlapped DOF map
    virtual Teuchos::RefCountPtr<const Epetra_Map> 
    getOverlapMap() const;

    //! Get Jacobian graph
    virtual Teuchos::RefCountPtr<const Epetra_CrsGraph> 
    getJacobianGraph() const;

    //! Get overlap Jacobian graph
    virtual Teuchos::RefCountPtr<const Epetra_CrsGraph> 
    getOverlapJacobianGraph() const;

    //! Get number of nodes per element
    virtual int getNumNodesPerElement() const;


  private:

    //! Private to prohibit copying
    CZeroDiscretization(const CZeroDiscretization&);

    //! Private to prohibit copying
    CZeroDiscretization& operator=(const CZeroDiscretization&);

  protected:
    
    //! Coordinates of mesh nodes
    std::vector<double> x;

    //! Epetra communicator
    Teuchos::RefCountPtr<const Epetra_Comm> comm;

    //! Element factory
    FEApp::ElementFactory elemFactory;

    //! Element mesh
    Teuchos::RefCountPtr<FEApp::Mesh> mesh;

    //! Element map
    Teuchos::RefCountPtr<Epetra_Map> elem_map;

    //! Unknown Map
    Teuchos::RefCountPtr<Epetra_Map> map;

    //! Overlapped unknown map
    Teuchos::RefCountPtr<Epetra_Map> overlap_map;

    //! Jacobian matrix graph
    Teuchos::RefCountPtr<Epetra_CrsGraph> graph;

    //! Overlapped Jacobian matrix graph
    Teuchos::RefCountPtr<Epetra_CrsGraph> overlap_graph;

    //! Processor ID
    unsigned int myPID;

    //! Number of elements on this processor
    unsigned int numMyElements;

    //! Number of nodes per element
    unsigned int nodes_per_element;

    //! Number of equations per node
    unsigned int neq;

  };

}

#endif // FEAPP_CZERODISCRETIZATION_HPP
