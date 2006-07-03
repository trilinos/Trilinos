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

#ifndef NOX_EPETRA_INTERFACE_PRECONDITIONER_H
#define NOX_EPETRA_INTERFACE_PRECONDITIONER_H

#include "NOX_Common.H"
  
// Forward declarations
class Epetra_Vector;
class Epetra_Operator;
namespace Teuchos {
  class ParameterList;
}

namespace NOX {
namespace Epetra {
namespace Interface {

  /*! \brief Used by NOX::Epetra to provide a link to the
    external code for Precondtioner fills.
    
    This is only required if the user wishes to supply their own 
    preconditioner operator.   
  */
class Preconditioner {

public:
 
  //! Constructor
  Preconditioner() {};

  //! Destructor
  virtual ~Preconditioner() {};

  //! Computes a user defined preconditioner.
  virtual bool computePreconditioner(const Epetra_Vector& x, 
				     Epetra_Operator& M,
				     Teuchos::ParameterList* precParams = 0) = 0;
 
};
} // namespace Interface
} // namespace Epetra
} // namespace NOX

#endif
