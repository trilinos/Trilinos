// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
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
