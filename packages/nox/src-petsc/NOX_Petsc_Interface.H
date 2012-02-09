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

#ifndef NOX_PETSC_INTERFACE_H
#define NOX_PETSC_INTERFACE_H
  
// Forward declarations

namespace NOX {
namespace Petsc {

  /*! \brief Used by NOX::Petsc::Group to provide a link to the
   *  external code for F and Jacobian evaluations (if needed). */
  /*! At a minimum the user must supply a computeF().  
   * If the user is supplying their own Petsc Mat for the Jacobian,
   * then they must supply a computeJacobian().  Preconditioning 
   * functions are only needed if the user wishes to supply their own 
   * preconditioner.  Otherwise Petsc SLES has built-in preconditioners 
   * to use.  */
class Interface {

public:
 
  //! Constructor
  Interface() {};

  //! Destructor
  virtual ~Interface() {};

  //! Compute F given the specified input vector, x.  Returns true if computation was successful.
  virtual bool computeF(const Vec& x, Vec& RHS) = 0;
  
  //! Compute Jacobian given the specified input vector, x.  Returns true if computation was successful.
  virtual bool computeJacobian(const Vec& x, Mat& Jac) = 0;

  //! Compute the matrix M that will be used as the preconditioner.  Returns true if computation was successful.
  virtual bool computePreconditioner(Mat& M) = 0;

  //! Return the action of the preconditioner on a vector.  Returns true if computation was successful.
  virtual bool preconditionVector(Vec& y) = 0;

};
} // namespace Petsc
} // namespace NOX

#endif
