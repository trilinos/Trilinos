// $Id: LOCA_Epetra_Interface_FreeEnergy.H,v 1.8 2007/04/13 19:28:25 agsalin Exp $ 
// $Source: /space/CVS/Trilinos/packages/nox/src-loca/src-epetra/LOCA_Epetra_Interface_FreeEnergy.H,v $ 

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
//  $Source: /space/CVS/Trilinos/packages/nox/src-loca/src-epetra/LOCA_Epetra_Interface_FreeEnergy.H,v $
//  $Author: agsalin $
//  $Date: 2007/04/13 19:28:25 $
//  $Revision: 1.8 $
// ************************************************************************
//@HEADER

#ifndef LOCA_EPETRA_INTERFACE_FREEENERGY_H
#define LOCA_EPETRA_INTERFACE_FREEENERGY_H

#include "NOX_Common.H"
  
// Forward declarations
class Epetra_Vector;

namespace LOCA {
  namespace Epetra {
    namespace Interface {

      /*!
       * \brief Used by LOCA::Epetra::Group to provide a link to the
       * external code for computing the free energy.
       */

      class FreeEnergy {

      public:
 
	//! Constructor
	FreeEnergy() {};

	//! Destructor
	virtual ~FreeEnergy() {};

	/*! 
	 * \brief Call user routine for computing the FreeEnergy
         * of a system, for use in Phase Transition Tracking alg
	 */

	virtual double
        computeFreeEnergy(const Epetra_Vector& x) = 0;

      };
    } // namespace Interface
  } // namespace Epetra
} // namespace LOCA

#endif
