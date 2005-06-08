//
// @HEADER
// ***********************************************************************
// 
//                           Rythmos Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef Rythmos_EXAMPLE_APPLICATION_H
#define Rythmos_EXAMPLE_APPLICATION_H

//#include "Thyra_VectorBase.hpp"
#include "Epetra_Vector.h"
#include "Teuchos_RefCountPtr.hpp"

//-----------------------------------------------------------------------------
// Class         : ExampleApplication
// Purpose       : Example application code with DE: \dot{x}=\lambda x
// Special Notes :
// Creator       : Todd Coffey, SNL
// Creation Date : 05/05/05
//-----------------------------------------------------------------------------
class ExampleApplication
{
  public:
    
    // Destructor
    ~ExampleApplication() {};
    
    // Default Constructor
    ExampleApplication() {};

    // Constructor
    ExampleApplication(double lambda, int numElements);

    // Evaluate residual:
    int evalResidual(Epetra_Vector *y, const Epetra_Vector &x, double t);
    
    // return ODE decay coefficient
    double getCoeff();
    
    // Return nominal x0 vector
    const Teuchos::RefCountPtr<Epetra_Vector> get_x0();

    // Return epetra_map 
    const Teuchos::RefCountPtr<const Epetra_Map> get_epetra_map()
      { return(epetra_map_); };

    // Return epetra_comm
    const Teuchos::RefCountPtr<const Epetra_Comm> get_epetra_comm()
      { return(epetra_comm_); };

  protected:

    // Coefficient for ODE
    double lambda_;
    // Epetra Comm:
    Teuchos::RefCountPtr<Epetra_Comm> epetra_comm_;
    // Epetra Map:
    Teuchos::RefCountPtr<Epetra_Map> epetra_map_;
    // Number of unknowns:
    int numElements_;

};


#endif // Rythmos_EXAMPLE_APPLICATION_H
