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

#ifndef FEAPP_APPLICATION_HPP
#define FEAPP_APPLICATION_HPP

#include <vector>

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "FEApp_AbstractBC.hpp"
#include "FEApp_AbstractPDE.hpp"
#include "FEApp_AbstractQuadrature.hpp"
#include "FEApp_AbstractDiscretization.hpp"
#include "FEApp_AbstractProblem.hpp"
#include "FEApp_TemplateTypes.hpp"

namespace FEApp {

  class Application {
  public:

    //! Constructor 
    Application(const std::vector<double>& coords,
		const Teuchos::RefCountPtr<const Epetra_Comm>& comm,
		const Teuchos::RefCountPtr<Teuchos::ParameterList>& params);

    //! Destructor
    ~Application();

    //! Get DOF map
    Teuchos::RefCountPtr<const Epetra_Map> getMap() const;

    //! Get Jacobian graph
    Teuchos::RefCountPtr<const Epetra_CrsGraph> getJacobianGraph() const;

    //! Get initial solution
    Teuchos::RefCountPtr<const Epetra_Vector> getInitialSolution() const;

    //! Compute global residual
    void computeGlobalResidual(const Epetra_Vector& x,
			       Epetra_Vector& f);

    //! Compute global Jacobian
    void computeGlobalJacobian(const Epetra_Vector& x,
			       Epetra_Vector& f,
			       Epetra_CrsMatrix& jac);

  private:
    
    //! Private to prohibit copying
    Application(const Application&);

    //! Private to prohibit copying
    Application& operator=(const Application&);

  protected:
    
    //! Element discretization
    Teuchos::RefCountPtr<FEApp::AbstractDiscretization> disc;

    //! Boundary conditions
    std::vector< Teuchos::RefCountPtr<const FEApp::AbstractBC> > bc;

    //! Quadrature rule
    Teuchos::RefCountPtr<const FEApp::AbstractQuadrature> quad;

    //! PDE equations
    FEApp::AbstractPDE_TemplateManager<ValidTypes> pdeTM;

    //! Initial solution vector
    Teuchos::RefCountPtr<const Epetra_Vector> initial_x;

    //! Importer for overlapped data
    Teuchos::RefCountPtr<Epetra_Import> importer;

    //! Exporter for overlapped data
    Teuchos::RefCountPtr<Epetra_Export> exporter;

    //! Overlapped solution vector
    Teuchos::RefCountPtr<Epetra_Vector> overlapped_x;

    //! Overlapped residual vector
    Teuchos::RefCountPtr<Epetra_Vector> overlapped_f;

    //! Overlapped Jacobian matrix
    Teuchos::RefCountPtr<Epetra_CrsMatrix> overlapped_jac;

  };

}

#endif // FEAPP_APPLICATION_HPP
