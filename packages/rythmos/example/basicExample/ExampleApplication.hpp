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
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif // HAVE_MPI

#include "Rythmos_ConfigDefs.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

#include <vector>

class ExampleApplication
{
  public:
    
    // Destructor
    ~ExampleApplication() {};
    
    // Default Constructor
    ExampleApplication() {};

    // Constructor
    ExampleApplication(Teuchos::ParameterList &params);

    // Evaluate residual:
    int evalResidual(Epetra_Vector *y, const Epetra_Vector &x, double t);
    
    // return ODE decay coefficient
    Teuchos::RefCountPtr<const Epetra_Vector> get_coeff() const;
    
    // Return nominal x0 vector
    Teuchos::RefCountPtr<Epetra_Vector> get_x0() const;

    // Return epetra_map 
    Teuchos::RefCountPtr<Epetra_Map> get_epetra_map() const;

    // Return epetra_comm
    Teuchos::RefCountPtr<Epetra_Comm> get_epetra_comm() const;

  private:

    // Epetra Comm:
    Teuchos::RefCountPtr<Epetra_Comm> epetra_comm_ptr_;
    // Epetra Map:
    Teuchos::RefCountPtr<Epetra_Map> epetra_map_ptr_;
    
    // Global number of unknowns:
    int numElements_;
    // Coefficients for ODE
    double lambda_min_;
    double lambda_max_;
    Teuchos::RefCountPtr<Epetra_Vector> lambda_ptr_;
    // Constant initial condition for the problem:
    double x0_;

};


#endif // Rythmos_EXAMPLE_APPLICATION_H
