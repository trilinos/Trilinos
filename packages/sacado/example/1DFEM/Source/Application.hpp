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

#ifndef APPLICATION_HPP
#define APPLICATION_HPP

#include <vector>

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "CZeroDiscretization.hpp"
#include "AbstractBC.hpp"
#include "AbstractQuadrature.hpp"

#include "AbstractPDE.hpp"
#include "Sacado_Fad_DFad.hpp"

class Application {
public:

  //! Constructor 
  Application(const Teuchos::RefCountPtr<CZeroDiscretization> c0_disc,
	      const std::vector< Teuchos::RefCountPtr<const AbstractBC> > bcs,
	      const Teuchos::RefCountPtr<const AbstractQuadrature>& quadRule);

  //! Destructor
  ~Application();

  //! Compute global residual
  void computeGlobalResidual(AbstractPDE<double>& pde,
			     const Epetra_Vector& x,
			     Epetra_Vector& f);

  //! Compute global Jacobian
  void computeGlobalJacobian(AbstractPDE< Sacado::Fad::DFad<double> >& pde,
			     const Epetra_Vector& x,
			     Epetra_Vector& f,
			     Epetra_CrsMatrix& jac);

private:

  //! Private to prohibit copying
  Application(const Application&);

  //! Private to prohibit copying
  Application& operator=(const Application&);

protected:

  //! Element discretization
  Teuchos::RefCountPtr<CZeroDiscretization> disc;

  //! Boundary conditions
  std::vector< Teuchos::RefCountPtr<const AbstractBC> > bc;

   //! Quadrature rule
  Teuchos::RefCountPtr<const AbstractQuadrature> quad;

  //! Importer for overlapped data
  Epetra_Import importer;

  //! Exporter for overlapped data
  Epetra_Export exporter;

  //! Overlapped solution vector
  Epetra_Vector overlapped_x;

  //! Overlapped residual vector
  Epetra_Vector overlapped_f;

  //! Overlapped Jacobian matrix
  Epetra_CrsMatrix overlapped_jac;

};

#endif // APPLICATION_HPP
