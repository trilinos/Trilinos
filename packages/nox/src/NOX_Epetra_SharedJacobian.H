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

#ifndef NOX_PETRA_SHARED_JACOBIAN_H
#define NOX_PETRA_SHARED_JACOBIAN_H

class Epetra_CrsMatrix;
class Epetra_RowMatrix;

namespace NOX {
namespace Epetra {

class Group;			// forward declaration

  /*! \brief Shared Jacobian for NOX::Epetra::Group objects.
   * 
   * Due to the size of the Jacobian, we cannot afford to have
   * multiple copies. Instead we implement a shared Jacobian that
   * that all groups use. 
   *
   * NOTE: I have added a preconditioner Matrix to the SharedJacobian
   * in anticipation of having a separate Matrix for the Jacobian 
   * and the Preconditioner.  This option will be in Trilinos by the end 
   * of Feb 02. 
   */

class SharedJacobian {

public:

  //! Constructor
  SharedJacobian(Epetra_RowMatrix& j);

  //! Constructor with separate Matrices for the Jacobian and preconditioner 
  SharedJacobian(Epetra_RowMatrix& j, Epetra_RowMatrix& p);

  //! Destructor
  ~SharedJacobian();

  //! Get a reference to the Jacobian AND take ownership
  Epetra_RowMatrix& getJacobian(const Group* newowner);

  //! Return a const pointer to the Jacobian
  const Epetra_RowMatrix& getJacobian() const;

  //! Return true if grp is the owner of the Jacobian
  bool isOwner(const Group* grp) const;

  //! Get a reference to the Jacobian AND take ownership
  Epetra_RowMatrix& getPrec(const Group* newowner);

  //! Return a const pointer to the Jacobian
  const Epetra_RowMatrix& getPrec() const;

  //! This is a hack that is needed to get MPSalsa running.  I will remove this when the Epetra_Vbr is fixed to get a view with Aztec2Petra.
  bool setJacobian(Epetra_RowMatrix& j);

protected:

  //! Pointer to Epetra Jacobian matrix 
  Epetra_RowMatrix* jacobian;

  //! Pointer to Epetra Preconditioning Matrix 
  Epetra_RowMatrix* prec;

  //! Pointer to const group that owns Jacobian
  const Group* owner;
};
} // namespace Epetra
} // namespace NOX

#endif
