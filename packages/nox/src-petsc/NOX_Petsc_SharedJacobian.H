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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
