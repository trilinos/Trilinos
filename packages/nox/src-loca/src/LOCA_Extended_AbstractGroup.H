// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_EXTENDED_ABSTRACTGROUP_H
#define LOCA_EXTENDED_ABSTRACTGROUP_H

#include "NOX_Abstract_Group.H"    // base class

namespace LOCA {
  namespace Continuation {
    class AbstractGroup;
  }
}

namespace LOCA { 

  //! %Extended groups and vectors
  namespace Extended { 

    /*! 
     * \brief %LOCA abstract interface for extended groups, derived from the 
     * NOX::Abstract::Group, i.e., an abstract interface for "super" groups
     * that have an underlying group component.
    */
    /*! 
     * Concrete implemenations of this interface must provide implementations 
     * of all of the methods in the NOX::Abstract::Group interface as well as 
     * the additional interface defined here.
    */
    class AbstractGroup : public virtual NOX::Abstract::Group {

    public:
  
      //! Default constructor.
      AbstractGroup() {}
 
      //! Destructor.
      virtual ~AbstractGroup() {}

      /*! 
       * @name Pure virtual methods 
       * These methods must be defined by any concrete implementation
       */
      //@{

      //! Assignment operator
      virtual AbstractGroup& operator=(const AbstractGroup& source) = 0;

      //! Return underlying group
      /*!
       * This method should the underlying group data member.
       */
      virtual const LOCA::Continuation::AbstractGroup& 
      getUnderlyingGroup() const = 0;
      
      //! Return underlying group
      /*!
       * This method should the underlying group data member.
       */
      virtual LOCA::Continuation::AbstractGroup& 
      getUnderlyingGroup() = 0;

      //@}

      
      /*! 
       * @name Virtual methods with default implementations 
       * These methods should be overloaded in a concrete implementation
       * if more appropriate/efficient approaches are available.
       */
      //@{

      //! Return base-level underlying group
      /*!
       * This method is intended
       * for composite groups (such as extended  %bifurcation groups) which
       * have an underlying group as a data member.  This method is supposed 
       * to return the base level group and has a default recursive 
       * implementation that should work in most cases.
       */
      virtual const LOCA::Continuation::AbstractGroup& 
      getBaseLevelUnderlyingGroup() const;
      
      //! Return base-level underlying group
       /*!
       * This method is intended
       * for composite groups (such as extended  %bifurcation groups) which
       * have an underlying group as a data member.  This method is supposed 
       * to return the base level group and has a default recursive 
       * implementation that should work in most cases.
       */
      virtual LOCA::Continuation::AbstractGroup& 
      getBaseLevelUnderlyingGroup();

      //@}

    };
  } // namespace Extended
} // namespace LOCA

#endif
