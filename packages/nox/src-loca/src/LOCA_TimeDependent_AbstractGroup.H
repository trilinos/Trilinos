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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_TIMEDEPENDENT_ABSTRACTGROUP_H
#define LOCA_TIMEDEPENDENT_ABSTRACTGROUP_H

#include "LOCA_Continuation_AbstractGroup.H"    // base class

namespace LOCA { 

  namespace TimeDependent {

    /*!  
     * \brief Interface to underlying groups for time dependent systems
     */
    /*!
     * This abstract class provides an interface for time dependent problems,
     * i.e., problems with a mass matrix (typically used in eignvalue or
     * Hopf calculations).  It provides two pure virtual methods, one to
     * compute the mass matrix and the second to apply the mass matrix to 
     * a vector.
     */
    class AbstractGroup : 
      public virtual LOCA::Continuation::AbstractGroup {

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
      
      //! Assignment operator
      virtual AbstractGroup& 
      operator=(const AbstractGroup& source) = 0;
      
      //! Compute the mass matrix of the system
      virtual NOX::Abstract::Group::ReturnType
      computeMassMatrix() = 0;

      //! Multiply the mass matrix by a vector.
      virtual NOX::Abstract::Group::ReturnType
      applyMassMatrix(const NOX::Abstract::Vector& input,
		      NOX::Abstract::Vector& result) const = 0;

      //@}

    };
  } // namespace TimeDependent
} // namespace LOCA

#endif
