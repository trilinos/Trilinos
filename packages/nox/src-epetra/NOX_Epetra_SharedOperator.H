// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef NOX_PETRA_SHARED_OPERATOR_H
#define NOX_PETRA_SHARED_OPERATOR_H

#include "NOX_Common.H"

class Epetra_Operator;

namespace NOX {
namespace Epetra {

class Group;			// forward declaration

  /*! \brief Shared operator object for NOX::Epetra::Group objects.
   * 
   * Due to the size of certain operators (i.e. the Jacobian and 
   * Preconditioning Matrices), we cannot afford to have
   * multiple copies. Instead we implement a shared operator class that
   * that all groups use. 
   */

class SharedOperator {

public:

  //! Constructor
  SharedOperator(Epetra_Operator& m);

  //! Destructor
  ~SharedOperator();

  //! Get a non-const reference to the underlying matrix (an epetra operator)
  Epetra_Operator& getOperator(const Group* newowner);

  //! Return a const reference to the Jacobian
  const Epetra_Operator& getOperator() const;

  //! Return true if grp is the owner of the Jacobian
  bool isOwner(const Group* grp) const;

  //! RPP: This is a hack that is needed to get MPSalsa running.  I will remove this when the Epetra_Vbr is fixed to get a matrix "View" with Aztec2Petra.  Essentially, this function allows the underlying Epetra_Operator to be changed.  We need this if the operator is destroyed and a new one is created during the NOX::Epetra::Interface::computeJacobian() call. 
  bool setOperator(Epetra_Operator& m);

protected:

  //! Pointer to Epetra Jacobian matrix 
  Epetra_Operator* matrix;

  //! Pointer to const group that owns Jacobian
  const Group* owner;
};
} // namespace Epetra
} // namespace NOX

#endif
