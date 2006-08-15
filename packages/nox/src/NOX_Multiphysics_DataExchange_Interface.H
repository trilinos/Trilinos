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

#ifndef NOX_MULTIPHYSICS_DATAEXCHANGE_INTERFACE_H
#define NOX_MULTIPHYSICS_DATAEXCHANGE_INTERFACE_H

#include "NOX_Common.H"
  
// Forward declarations
class Epetra_Vector;

namespace NOX {
namespace Multiphysics {
namespace DataExchange {

  /*! 
    \brief Provides a set of interfaces for users to provide information about exchanging data between registered NOX solvers
  */

  /*!
    \brief Supplies NOX with the set nonlinear equations.
 
    This is the minimum required information to solve a nonlinear
    problem using the NOX::Epetra objects for the linear algebra
    implementation.  Used by NOX::Epetra::Group to provide a link
    to the external code for residual fills.
  */
class Interface {

public:
 
  //! Type of fill that a computeF() method is used for.  
  /*! computeF() can be called for a variety of reasons:

  - To evaluate the function residuals.
  - To be used in an approximation to the Jacobian (finite difference or directional derivative).
  - To be used in an approximation to the preconditioner.
   */

  //! Constructor
  Interface() {};

  //! Destructor
  virtual ~Interface() {};

  //! Exchange data for all registered problems
  virtual bool exchangeAllData() = 0;
  
  //! Exchange data for a specified problem - brings needed data from others to this problem
  virtual bool exchangeDataTo(int solverId) = 0;
  
};
} // namespace DataExchange
} // namespace Multiphysics
} // namespace NOX

#endif
