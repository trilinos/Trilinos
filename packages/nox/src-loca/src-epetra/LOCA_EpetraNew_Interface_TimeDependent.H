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

#ifndef LOCA_EPETRANEW_INTERFACE_TIMEDEPENDENT_H
#define LOCA_EPETRANEW_INTERFACE_TIMEDEPENDENT_H

#include "NOX_Common.H"
#include "LOCA_EpetraNew_Interface_Required.H"
#include "NOX_Abstract_Group.H"
  
// Forward declarations
class Epetra_Vector;

namespace LOCA {
  class ParameterVector;
}

namespace LOCA {
  namespace EpetraNew {
    namespace Interface {

      /*! \brief 
       * Used by LOCA::EpetraNew::Group to provide a link to the
       * external code for RHS and Jacobian fills if needed. 
      */
      /*!
       * This interface is derived from the 
       * NOX::EpetraNew::Interface::Required and additionally provides a 
       * method for accessing a mass matrix. This is needed for linear
       * stability analysis and Hopf tracking.
      */
      class TimeDependent : public LOCA::EpetraNew::Interface::Required {

      public:
 
	//! Constructor
	TimeDependent() {};

	//! Destructor
	virtual ~TimeDependent() {};

	/*! 
	 * \brief Call user routine for applying the Mass Matrix to a 
	 * vector, if using matrix-free mode for mass matrix
	 */
	virtual NOX::Abstract::Group::ReturnType 
        applyMassMatrix(const NOX::Epetra::Vector& input,
                        NOX::Epetra::Vector& result) const = 0;	   

      };
    } // namespace Interface
  } // namespace Epetra
} // namespace LOCA

#endif
