// $Id: LOCA_PhaseTransition_MooreSpence_AbstractGroup.H,v 1.7 2008/07/22 23:29:57 etphipp Exp $ 
// $Source: /space/CVS/Trilinos/packages/nox/src-loca/src/LOCA_PhaseTransition_MooreSpence_AbstractGroup.H,v $ 

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
//  $Source: /space/CVS/Trilinos/packages/nox/src-loca/src/LOCA_PhaseTransition_MooreSpence_AbstractGroup.H,v $
//  $Author: etphipp $
//  $Date: 2008/07/22 23:29:57 $
//  $Revision: 1.7 $
// ************************************************************************
//@HEADER

#ifndef LOCA_PHASE_TRANSITION_ABSTRACTGROUP_H
#define LOCA_PHASE_TRANSITION_ABSTRACTGROUP_H

#include "LOCA_MultiContinuation_AbstractGroup.H"    // base class

namespace LOCA { 
  
  //! Groups and vectors for turning point bifurcations
  namespace PhaseTransition { 

    /*! 
     * \brief Interface to underlying groups for phase transition calculations
     */
    /*!
     * This abstract class provides the required interface for underlying
     * groups to locate phase transitions using the bordering
     * algorithm from the Salinger&Frink (2003) paper.
     *
     * This class is derived from the LOCA::MultiContinuation::AbstractGroup 
     * and declares a pure virtual method for computing the free energy.
     */
     class AbstractGroup : 
       public virtual LOCA::MultiContinuation::AbstractGroup {

    public:

	//! Default constructor.
	AbstractGroup() {}

	//! Destructor
	virtual ~AbstractGroup() {}

	/*! 
	 * @name Pure virtual methods 
	 * These methods must be defined by any concrete implementation
	 */
	//@{

	/*! 
	 * \brief Computes the free energy at the current solution and parameter values
	 */
	virtual double computeFreeEnergy() = 0;

    }; // class AbstractGroup

  } // namespace PhaseTransition

} // namespace LOCA

#endif
