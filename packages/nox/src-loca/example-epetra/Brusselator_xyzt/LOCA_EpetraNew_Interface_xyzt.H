//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_EPETRANEW_INTERFACE_XYZT_H
#define LOCA_EPETRANEW_INTERFACE_XYZT_H

#ifdef HAVE_MPI

#include <iostream>
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_DistObject.h"
#include "Epetra_MpiComm.h"

#include "NOX_Common.H"
#include "NOX_EpetraNew_Interface_Required.H"
#include "NOX_EpetraNew_Interface_Jacobian.H"
  
namespace LOCA {
namespace EpetraNew {

  /*! 
    \brief Provides a set of interfaces for users to provide information about the nonlinear problem to NOX.  

    Contains interfaces for the user to supply (1) the evaluation of the nonlinear equations, (2) the Jacobian, and (3) any preconditioning if required.
  */
namespace Interface {

  /*!
    \brief Supplies NOX with the set nonlinear equations.
 
    This is the minimum required information to solve a nonlinear
    problem using the NOX::EpetraNew objects for the linear algebra
    implementation.  Used by NOX::EpetraNew::Group to provide a link
    to the external code for residual fills.
  */
class xyzt : public NOX::EpetraNew::Interface::Required,
             public NOX::EpetraNew::Interface::Jacobian
{

public:
 
  //! Code that takes several replicas of problem interface, each working on a different
  //! time step, and assembles a single space-time system

  //! Constructor
  xyzt(NOX::EpetraNew::Interface::Required &iReq,  NOX::EpetraNew::Interface::Jacobian &iJac,
     Epetra_Vector &splitVec, Epetra_RowMatrix &splitJac, Epetra_MpiComm &globalComm,
     int replica);

  //! Destructor
  virtual ~xyzt();

  //! Compute the function, F, given the specified input vector x.  Returns true if computation was successful.
  virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, 
			const FillType fillFlag);
  
  virtual bool computeJacobian(const Epetra_Vector& x);

  virtual Epetra_Vector& getSolution();
  virtual Epetra_RowMatrix& getJacobian();

private:
  
  // Info from split (3D) problem
  NOX::EpetraNew::Interface::Required &iReq;
  NOX::EpetraNew::Interface::Jacobian &iJac;
  Epetra_Vector& splitVec;
  Epetra_RowMatrix& splitJac;
  Epetra_MpiComm& globalComm;
  int replica;

  // xyzt objects created in constructor
  Epetra_RowMatrix* jacobian;
  Epetra_Vector* solution;

};
} // namespace Interface
} // namespace EpetraNew
} // namespace LOCA

#endif
#endif
