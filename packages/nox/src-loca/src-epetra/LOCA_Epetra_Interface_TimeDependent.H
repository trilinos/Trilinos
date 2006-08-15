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

#ifndef LOCA_EPETRA_INTERFACE_TIMEDEPENDENT_H
#define LOCA_EPETRA_INTERFACE_TIMEDEPENDENT_H

#include "NOX_Common.H"
#include "LOCA_Epetra_Interface_Required.H"
#include "NOX_Abstract_Group.H"
  
// Forward declarations
class Epetra_Vector;

namespace LOCA {
  class ParameterVector;
}

namespace LOCA {
  namespace Epetra {
    namespace Interface {

      /*!
       * \brief Used by LOCA::Epetra::Group to provide a link to the
       * external code for computing the shifted matrix.
       */
      /*!
       * This interface is derived from the 
       * NOX::Epetra::Interface::Jacobian and additionally provides a 
       * method for computing the shifted matrix \f$\alpha J + \beta M\f$. 
       * This is needed for linear stability analysis and Hopf tracking.
       */
      class TimeDependent : public virtual LOCA::Epetra::Interface::Required {

      public:
 
	//! Constructor
	TimeDependent() {};

	//! Destructor
	virtual ~TimeDependent() {};

	/*! 
	 * \brief Call user routine for computing the shifted matrix
	 * \f$\alpha J + \beta M\f$ where \f$J\f$ is the Jacobian matrix
	 * and \f$M\f$ is the mass matrix.
	 */
	virtual bool
        computeShiftedMatrix(double alpha, double beta, 
			     const Epetra_Vector& x,
			     Epetra_Operator& A) = 0;

      };
    } // namespace Interface
  } // namespace Epetra
} // namespace LOCA

#endif
