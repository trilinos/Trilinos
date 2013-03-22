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

#ifndef LOCA_EPETRA_INTERFACE_MASSMATRIX_H
#define LOCA_EPETRA_INTERFACE_MASSMATRIX_H

#include "NOX_Common.H"
#include "NOX_Abstract_Group.H"
  
// Forward declarations
class Epetra_Vector;

namespace LOCA {
  namespace Epetra {
    namespace Interface {
      
      /*! 
	\brief 
	Used by LOCA::Epetra::Group to provide a link to the external code 
	for the MassMatrix (coefficients of time dependent terms). 

	This is used for Hopf bifurcation tracking, linear stability analysis, 
	and space-time solutions (xyzt).
      */
      class MassMatrix  {
	
      public:
	
	//! Constructor
	MassMatrix() {};
	
	//! Destructor
	virtual ~MassMatrix() {};

        /*! 
	  \brief 
	  Compute MassMatrix given the specified input vector x.  
	  Returns true if computation was successful.
	*/
        virtual bool computeMassMatrix(const Epetra_Vector& x) = 0;

        /*! 
	  \brief Routines used in XYZT to set the old solution, 
	  the one from the previous time step. 

	  There is a different routine for first step,
	  where the old solution is set by the user and not part of the 
	  solution vector. These routines are used by space-time (xyzt) 
	  problems, where the residual vector is a function of the previous
	  solution, which is also being solved for, and where the MassMatrix 
	  is calculated as a function of a different solution vector then 
	  the Jacobian (that is, the previous time step).
	  The  timeStep  argument is sent so the use can set the global time,
	  for cases when computeF, computeJacobian, computeMassMatrix fills
	  are functions of time (nonautonomous systems).
	*/
        virtual void setOldSolution(const Epetra_Vector& x,const int timeStep) {
	  cout << "WARNING: " 
	       << "LOCA::Epetra::Interface::MassMatrix::setOldSolution"
	       << "\n\tempty default implementation not overloaded!" << std::endl;
	};
	
	//! See setOldSolution description.
        virtual void setOldSolutionFirstStep() { 
	  cout << "WARNING: " 
	       << "LOCA::Epetra::Interface::MassMatrix::setOldSolutionFirstStep"
	       << "\n\tempty default implementation not overloaded!" << std::endl;
	};

        /*! 
	  \brief Provides data to application for output files.

	  This routine is called from Interface::xyzt::printSolution() just
	  before the call to Interface::Required::printSolution(x,param),
	  and gives the application some indices that can be used for
	  creating a unique name/index for the output files.
	*/
        virtual void dataForPrintSolution(const int conStep, const int timeStep,
                                          const int totalTimeSteps) {};

      };
    } // namespace Interface
  } // namespace Epetra
} // namespace LOCA

#endif
