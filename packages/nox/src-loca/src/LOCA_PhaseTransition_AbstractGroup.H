// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
