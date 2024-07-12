// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_TURNINGPOINT_MOORESPENCE_ABSTRACTGROUP_H
#define LOCA_TURNINGPOINT_MOORESPENCE_ABSTRACTGROUP_H

#include "LOCA_MultiContinuation_AbstractGroup.H"    // base class

namespace LOCA {

  //! Groups and vectors for turning point bifurcations
  namespace TurningPoint {

    /*!
     * \brief Groups and vectors for locating turning point bifurcations
     * using the Moore-Spence formulation.
     */
    namespace MooreSpence {

      /*!
       * \brief Interface to underlying groups for turning point calculations
       * using the Moore-Spence formulation.
       */
      /*!
       * This abstract class provides the required interface for underlying
       * groups to locate turning points using the bordering
       * algorithm for the Moore-Spence turning point formulation (see
       * LOCA::TurningPoint::MooreSpence::ExtendedGroup for a description of
       * the governing equations).
       *
       * This class is derived from the LOCA::MultiContinuation::AbstractGroup
       * and declares several pure virtual methods compute various
       * derivatives of \f$ Jn\f$ for a given \f$ n\f$.  Default implementations
       * for the derivatives using finite differencing are implemented in the
       * LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup.
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
     * \brief Computes the derivative \f$\partial Jn/\partial p\f$
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDJnDpMulti(const std::vector<int>& paramIDs,
              const NOX::Abstract::Vector& nullVector,
              NOX::Abstract::MultiVector& result,
              bool isValid) = 0;

    /*!
     * \brief Computes the directional derivative
     * \f$\frac{\partial Jn}{\partial x} a\f$ for the given direction
     * \f$ a\f$.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDJnDxaMulti(const NOX::Abstract::Vector& nullVector,
               const NOX::Abstract::MultiVector& aVector,
               NOX::Abstract::MultiVector& result) = 0;

    /*!
     * \brief Computes the directional derivative
     * \f$\frac{\partial Jn}{\partial x} a\f$ for the given direction
     * \f$ a\f$.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDJnDxaMulti(const NOX::Abstract::Vector& nullVector,
               const NOX::Abstract::Vector& JnVector,
               const NOX::Abstract::MultiVector& aVector,
               NOX::Abstract::MultiVector& result) = 0;

    /*!
     * \brief Computes the derivative
     * \f$\frac{\partial w^TJn}{\partial x}\f$.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDwtJnDxMulti(const NOX::Abstract::MultiVector& w,
                const NOX::Abstract::Vector& nullVector,
                NOX::Abstract::MultiVector& result) = 0;

    //@}

    //! Bring NOX::Abstract::Group::operator=() into scope
    using NOX::Abstract::Group::operator=;

      }; // class AbstractGroup

    } // namespace MooreSpence

  } // namespace TurningPoint

} // namespace LOCA

#endif
