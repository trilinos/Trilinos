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

#ifndef FEAPP_ABSTRACTDISCRETIZATION_HPP
#define FEAPP_ABSTRACTDISCRETIZATION_HPP

#include <vector>

#include "Teuchos_RefCountPtr.hpp"

#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"

#include "FEApp_Mesh.hpp"

namespace FEApp {

  class AbstractDiscretization {
  public:

    //! Constructor
    AbstractDiscretization() {}

    //! Destructor
    virtual ~AbstractDiscretization() {}

    //! Create element mesh
    virtual void createMesh() = 0;

    //! Create DOF maps
    virtual void createMaps() = 0;

    //! Create Jacobian graph
    virtual void createJacobianGraphs() = 0;

    //! Get element mesh
    virtual Teuchos::RefCountPtr<const FEApp::Mesh> 
    getMesh() const = 0; 

    //! Get DOF map
    virtual Teuchos::RefCountPtr<const Epetra_Map> 
    getMap() const = 0;

    //! Get overlapped DOF map
    virtual Teuchos::RefCountPtr<const Epetra_Map> 
    getOverlapMap() const = 0;

    //! Get Jacobian graph
    virtual Teuchos::RefCountPtr<const Epetra_CrsGraph> 
    getJacobianGraph() const = 0;

    //! Get overlap Jacobian graph
    virtual Teuchos::RefCountPtr<const Epetra_CrsGraph> 
    getOverlapJacobianGraph() const = 0;

    //! Get number of nodes per element
    virtual int getNumNodesPerElement() const = 0;


  private:

    //! Private to prohibit copying
    AbstractDiscretization(const AbstractDiscretization&);

    //! Private to prohibit copying
    AbstractDiscretization& operator=(const AbstractDiscretization&);

  };

}

#endif // FEAPP_ABSTRACTDISCRETIZATION_HPP
