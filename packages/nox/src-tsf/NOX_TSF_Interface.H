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

#ifndef NOX_TSF_INTERFACE_H
#define NOX_TSF_INTERFACE_H


// Forward declarations
namespace NOX 
{
  namespace TSF 
  {
    class Vector;
  }
}

namespace TSFExtended
{
  class LinearOp;
  class NonlinearOp;
}
namespace NOX {
  namespace TSF {

    /*! \brief Virtual interface for NOX::BLAS::Group to provide a
     *  link to the external code for RHS and Jacobian fills.
     *
     * This class is not needed, as TSFExtended::NonlinearOperator
     * provides all these capabilities. - KL 13 July 2004. 
     */
    class Interface {

    public:
 
      //! Constructor
      Interface() {};

      //! Destructor
      virtual ~Interface() {};

      //! Get initial guess
      virtual const NOX::TSF::Vector& getInitialGuess() = 0;

      //! Compute RHS given the specified input vector x.  Returns true if computation was successful.
      virtual bool computeF(NOX::TSF::Vector& rhs, const NOX::TSF::Vector &x) = 0;
  
      //! Compute Jacobian given the specified input vector x.  Returns true if computation was successful.
      //      virtual bool computeJacobian(TSFExtended::LinearOp& J, const NOX::TSF::Vector & x) = 0;
      virtual bool computeJacobian(NOX::TSF::SharedOperator& J, const NOX::TSF::Vector & x) = 0;
    };
  } // namespace TSF
} // namespace NOX

#endif
