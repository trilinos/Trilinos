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

#ifndef NOX_EPETRA_INTERFACE_JACOBIAN_H
#define NOX_EPETRA_INTERFACE_JACOBIAN_H

#include "NOX_Common.H"
  
// Forward declarations
class Epetra_Vector;
class Epetra_Operator;

namespace NOX {
namespace Epetra {
namespace Interface {

  /*! \brief Used by NOX::Epetra to provide a link to the
    external code for Jacobian fills.

    This is only required if the user wishes to supply their own Jacobian 
    operator.   
  */
class Jacobian {

public:
 
  //! Constructor.
  Jacobian() {};

  //! Destructor.
  virtual ~Jacobian() {};

  /*! Compute Jacobian given the specified input vector x.  Returns
    true if computation was successful.
   */
  virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) = 0;

};
} // namespace Interface
} // namespace Epetra
} // namespace NOX

#endif
