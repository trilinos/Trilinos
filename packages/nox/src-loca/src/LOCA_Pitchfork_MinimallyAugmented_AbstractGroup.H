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

#ifndef LOCA_PITCHFORK_MINIMALLYAUGMENTED_ABSTRACTGROUP_H
#define LOCA_PITCHFORK_MINIMALLYAUGMENTED_ABSTRACTGROUP_H

#include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H" // base class
#include "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"    // base class

namespace LOCA { 
  
  namespace Pitchfork { 

    /*! 
     * \brief Groups and vectors for locating pitchfork bifurcations
     * using the minimally augmented turning point formulation.
     */
    namespace MinimallyAugmented {

      /*! 
       * \brief Interface to underlying groups for pitchfork calculations
       * using the minimally augmented formulation.
       */
      /*!
       * This abstract class provides the required interface for underlying
       * groups to locate pitchforks using the minimally augmented pitchfork
       * formulation (see 
       * LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup for a 
       * description of the governing equations).  
       *
       * This class is derived from the 
       * LOCA::Pitchfork::MooreSpence::AbstractGroup and
       * LOCA::TurningPoint::MinimallyAugmented::AbstractGroup and does not
       * declare any new virtual methods.
       */
      class AbstractGroup : 
	public virtual LOCA::TurningPoint::MinimallyAugmented::AbstractGroup,
	public virtual LOCA::Pitchfork::MooreSpence::AbstractGroup {

      public:
  
	//! Default constructor.
	AbstractGroup() {}

	//! Destructor
	virtual ~AbstractGroup() {}

	//! Bring NOX::Abstract::Group::operator=() into scope
	using NOX::Abstract::Group::operator=;

      }; // class AbstractGroup

    } // namespace MinimallyAugmented

  } // namespace Pitchfork

} // namespace LOCA

#endif
