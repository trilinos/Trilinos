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

#ifndef NOX_PETSC_SHARED_JACOBIAN_H
#define NOX_PETSC_SHARED_JACOBIAN_H

#include "petscmat.h"

namespace NOX {
namespace Petsc {

class Group;  // forward declaration

  /*! \brief Shared Jacobian for NOX::Petsc::Group objects.
   * 
   * Due to the size of the Jacobian, we cannot afford to have
   * multiple copies. Instead we implement a shared Jacobian that
   * all groups use. 
   *
   * NOTE: An additional matrix is included with the SharedJacobian
   * in anticipation of having a separate Matrix for the Jacobian 
   * and the Preconditioner.  This option will be available in the future.
   */

class SharedJacobian {

public:

  //! Constructor
  SharedJacobian(Mat& j);

  //! Constructor with separate Matrices for the Jacobian and preconditioner 
  SharedJacobian(Mat& j, Mat& p);

  //! Destructor
  ~SharedJacobian();

  //! Get a reference to the Jacobian AND take ownership
  Mat& getJacobian(const Group* newowner);

  //! Return a const pointer to the Jacobian
  const Mat& getJacobian() const;

  //! Return true if grp is the owner of the Jacobian
  bool isOwner(const Group* grp) const;

  //! Get a reference to the Jacobian AND take ownership
  Mat& getPrec(const Group* newowner);

  //! Return a const pointer to the Jacobian
  const Mat& getPrec() const;

protected:

  //! Pointer to Petsc Jacobian matrix 
  Mat* jacobian;

  //! Pointer to Petsc Preconditioning Matrix 
  Mat* prec;

  //! Pointer to const group that owns Jacobian
  const Group* owner;
};
} // namespace Petsc
} // namespace NOX

#endif
