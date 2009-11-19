// $Id: LOCA_PhaseTransition_MooreSpence_ExtendedVector.H,v 1.6 2007/06/21 16:22:52 rhoope Exp $ 
// $Source: /space/CVS/Trilinos/packages/nox/src-loca/src/LOCA_PhaseTransition_MooreSpence_ExtendedVector.H,v $ 

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
//  $Source: /space/CVS/Trilinos/packages/nox/src-loca/src/LOCA_PhaseTransition_MooreSpence_ExtendedVector.H,v $
//  $Author: rhoope $
//  $Date: 2007/06/21 16:22:52 $
//  $Revision: 1.6 $
// ************************************************************************
//@HEADER

#ifndef LOCA_PHASE_TRANSITION_EXTENDEDVECTOR_H
#define LOCA_PHASE_TRANSITION_EXTENDEDVECTOR_H

#include "LOCA_Extended_Vector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace PhaseTransition {
    class ExtendedMultiVector;
  }
}

namespace LOCA {

  namespace PhaseTransition {

      /*! 
       * \brief %Vector class to hold solution vectors, %Newton vectors, etc.
       * for the Phase Transition tracking formulation.
       */
      /*! 
       * This class uses the LOCA::Extended::Vector implementation to store the
       * solution1, solution2, and parameter components of the phase transition vector
       * and merely provides an interface for naming which components of the
       * vector these quantities correspond to.
      */
      class ExtendedVector : public LOCA::Extended::Vector {

	/*! 
	 * \brief Declare LOCA::PhaseTransition::MooreSpence::ExtendedMultiVector 
	 * as a friend class so it can call protected methods.
	 */
	friend class ExtendedMultiVector;

      public:
  
	//! %Constructor.
	ExtendedVector(
		   const Teuchos::RCP<LOCA::GlobalData>& global_data,
		   const NOX::Abstract::Vector& x1Vec,
		   const NOX::Abstract::Vector& x2Vec,
		   double ptp);

	//! Copy constructor.
	ExtendedVector(const ExtendedVector& source, 
		       NOX::CopyType type = NOX::DeepCopy);

	//! Destructor.
	virtual ~ExtendedVector();

	//! Assignment operator
	virtual NOX::Abstract::Vector& 
	operator=(const NOX::Abstract::Vector& y);

	//! Assignment operator
	virtual LOCA::Extended::Vector& 
	operator=(const LOCA::Extended::Vector& y);
  
	//! Assignment operator
	virtual ExtendedVector& operator=(const ExtendedVector& y);

	//! Cloning function
	virtual Teuchos::RCP<NOX::Abstract::Vector>
	clone(NOX::CopyType type = NOX::DeepCopy) const;

	//! Sets the %Vector by setting its three components 
	virtual void setVec(const NOX::Abstract::Vector& xVec,
			    const NOX::Abstract::Vector& nullVec,
			    double bifPar);

	//! Returns the solution1 vector component of extended vector
	virtual Teuchos::RCP<const NOX::Abstract::Vector> 
	X1() const;

	//! Returns the solution2 vector component of extended vector
	virtual Teuchos::RCP<const NOX::Abstract::Vector> 
	X2() const;

	//! Get %Bifurcation parameter
	virtual double PTP() const;

	//! Returns the solution vector component of extended vector
	virtual Teuchos::RCP<NOX::Abstract::Vector> X1();

	//! Returns the null vector component of extended vector
	virtual Teuchos::RCP<NOX::Abstract::Vector> X2();

	//! Get %Bifurcation parameter
	virtual double& PTP();

      protected:

	//! Default Constructor
	ExtendedVector(
		  const Teuchos::RCP<LOCA::GlobalData>& global_data);

	//! generate an extended multi-vector
	/*!
	 * Returns an empty multivector of type 
	 * LOCA::PhaseTransition::ExtendedMultiVector
	 */
	virtual Teuchos::RCP<LOCA::Extended::MultiVector>
	generateMultiVector(int nColumns, int nVectorRows, 
			    int nScalarRows) const;

      }; // class ExtendedVector
  } // namespace PhaseTransition
} // namespace LOCA

#endif
