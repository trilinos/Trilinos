// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
