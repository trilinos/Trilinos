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

#ifndef NOX_LAPACK_INTERFACE_H
#define NOX_LAPACK_INTERFACE_H

// Forward declarations
namespace NOX {
  namespace LAPACK {
    class Vector;
    class Matrix;
  }
}

namespace NOX {
  namespace LAPACK {

    /*! \brief Virtual interface for NOX::BLAS::Group to provide a
      link to the external code for RHS and Jacobian fills. */
    class Interface {

    public:
 
      //! Constructor
      Interface() {};

      //! Destructor
      virtual ~Interface() {};

      //! Get initial guess
      virtual const NOX::LAPACK::Vector& getInitialGuess() = 0;

      //! Compute RHS given the specified input vector x.  Returns true if computation was successful.
      virtual bool computeF(NOX::LAPACK::Vector& rhs, const NOX::LAPACK::Vector &x) = 0;
  
      //! Compute Jacobian given the specified input vector x.  Returns true if computation was successful.
      virtual bool computeJacobian(NOX::LAPACK::Matrix& J, const NOX::LAPACK::Vector & x) = 0;

    };
  } // namespace LAPACK
} // namespace NOX

#endif
