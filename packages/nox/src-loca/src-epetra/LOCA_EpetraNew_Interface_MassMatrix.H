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

#ifndef LOCA_EPETRANEW_INTERFACE_MASSMATRIX_H
#define LOCA_EPETRANEW_INTERFACE_MASSMATRIX_H

#include "NOX_Common.H"
#include "NOX_Abstract_Group.H"
  
// Forward declarations
class Epetra_Vector;

namespace LOCA {
  namespace EpetraNew {
    namespace Interface {

      /*! \brief 
       * Used by LOCA::EpetraNew::Group to provide a link to the
       * external code for the MassMatrix (coefficient of time dependent
       * terms).  This is used for: Hopf bifurcation tracking, linear
       * stability analysis, and space-time solutions (xyzt).
      */
      class MassMatrix  {

      public:
 
	//! Constructor
	MassMatrix() {};

	//! Destructor
	virtual ~MassMatrix() {};

        /*! Compute MassMatrix given the specified input vector x.  Returns
         *  true if computation was successful.
         *  */
        virtual bool computeMassMatrix(const Epetra_Vector& x) = 0;

        /*! Routine to reset the old solution, the one from the previous
	 * time step. Different routine for first step, where the old
	 * solution is set by the user and not part of the solution vector.
	 * These routines are used by space-time (xyzt) problems, where the
	 * MassMatrix is calculated as a function of a different solution
	 * vector then the Jacobian (that is, the previous time step).
         */
        virtual void setOldSolution(const Epetra_Vector& x) { cout <<
           "WARNING: LOCA::EpetraNew::Interface::MassMatrix::setOldSolution"
           << "\n\tempty default implementation not overloaded!" << endl;
	};
	//! See above
        virtual void setOldSolutionFirstStep() { cout << 
           "WARNING: LOCA::EpetraNew::Interface::MassMatrix::setOldSolutionFirstStep"
           << "\n\tempty default implementation not overloaded!" << endl;
	};

        /*! This routine is called from Interface::xyzt::printSolution() just
	 *  before the call to Interface::Required::printSolution(x,param),
	 *  and gives the application some indices that can be used for
	 *  creating a unique name/index for the output files.
         *  */
        virtual void dataForPrintSolution(const int conStep, const int timeStep,
                                          const int totalTimeSteps) {};

      };
    } // namespace Interface
  } // namespace Epetra
} // namespace LOCA

#endif
