// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef LOCA_EPETRA_INTERFACE_TIMEDEPENDENT_H
#define LOCA_EPETRA_INTERFACE_TIMEDEPENDENT_H

#include "NOX_Common.H"
#include "LOCA_Epetra_Interface_Required.H"
#include "NOX_Abstract_Group.H"
  
// Forward declarations
class Epetra_Operator;

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

        /*!
          \brief Routine used in XYZT to set x_dot and time in the interface.

          The computeF() routine for XYZT problems needs to be a function
          of x_dot, but th NOX/LOCA computeF() does not take x_dot as an 
          argument. This is used to set x_dot in the application interface
          so the subsequent call to computeF has the correct x_dot value.
          The  timeStep  argument is sent so the use can set the global time,
          for cases when computeF, computeJacobian, computeMassMatrix fills
          are functions of time (nonautonomous systems).
        */
        virtual void setXdot(const Epetra_Vector& xdot, const double time) {
	  std::cout << "WARNING: "
		    << "LOCA::Epetra::Interface::TimeDependent::setXdot"
		    << "\n\tempty default implementation not overloaded!" << std::endl;
        };

      };
    } // namespace Interface
  } // namespace Epetra
} // namespace LOCA

#endif
