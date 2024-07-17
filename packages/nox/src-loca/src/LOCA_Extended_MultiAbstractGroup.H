// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

      //! Override from NOX::Abstract::Group base class. Calls getUnderlyingGroup() from this class.
      virtual Teuchos::RCP<NOX::Abstract::Group> getNestedGroup();

      //! Override from NOX::Abstract::Group base class. Calls getUnderlyingGroup() from this class.
      virtual Teuchos::RCP<const NOX::Abstract::Group> getNestedGroup() const;
    };
  } // namespace Extended
} // namespace LOCA

#endif
