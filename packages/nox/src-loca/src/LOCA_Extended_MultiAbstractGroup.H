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

#ifndef LOCA_EXTENDED_MULTIABSTRACTGROUP_H
#define LOCA_EXTENDED_MULTIABSTRACTGROUP_H

#include "Teuchos_RCP.hpp"

#include "NOX_Abstract_Group.H"    // base class

namespace LOCA {
  namespace MultiContinuation {
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
    class MultiAbstractGroup : public virtual NOX::Abstract::Group {

    public:

      //! Default constructor.
      MultiAbstractGroup() {}

      //! Destructor.
      virtual ~MultiAbstractGroup() {}

      /*!
       * @name Pure virtual methods
       * These methods must be defined by any concrete implementation
       */
      //@{

      //! Return underlying group
      /*!
       * This method should the underlying group data member.
       */
      virtual
      Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
      getUnderlyingGroup() const = 0;

      //! Return underlying group
      /*!
       * This method should the underlying group data member.
       */
      virtual
      Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
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
      virtual
      Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
      getBaseLevelUnderlyingGroup() const;

      //! Return base-level underlying group
       /*!
       * This method is intended
       * for composite groups (such as extended  %bifurcation groups) which
       * have an underlying group as a data member.  This method is supposed
       * to return the base level group and has a default recursive
       * implementation that should work in most cases.
       */
      virtual
      Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
      getBaseLevelUnderlyingGroup();

      //@}

      //! Bring NOX::Abstract::Group::operator=() into scope
      using NOX::Abstract::Group::operator=;

    };
  } // namespace Extended
} // namespace LOCA

#endif
