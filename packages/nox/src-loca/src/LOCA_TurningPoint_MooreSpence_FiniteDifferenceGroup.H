// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_TURNINGPOINT_MOORESPENCE_FINITEDIFFERENCEGROUP_H
#define LOCA_TURNINGPOINT_MOORESPENCE_FINITEDIFFERENCEGROUP_H

#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"   // base class
#include "LOCA_MultiContinuation_FiniteDifferenceGroup.H" // base class

namespace LOCA {

  namespace TurningPoint {

    namespace MooreSpence {

      /*!
       * \brief Concrete class that provides concrete implementations of the
       * derivative computation methods of the
       * LOCA::TurningPoint::MooreSpence::AbstractGroup using first-order
       * finite differencing.
       */
      /*!
       * The finite-differencing calculations are actually implemented by
       * the LOCA::DerivUtils class, and a custom DerivUtils object can be
       * passed through the constructor of this class.  However, in the future
       * the calculations encapsulated in the DerivUtils class may be
       * incorporated directly into this class and other finite-differencing
       * child classes.
       */
      class FiniteDifferenceGroup :
    public virtual LOCA::TurningPoint::MooreSpence::AbstractGroup,
    public virtual LOCA::MultiContinuation::FiniteDifferenceGroup {

      public:

    //! Constructor.
    FiniteDifferenceGroup();

    //! Copy constructor
    FiniteDifferenceGroup(const FiniteDifferenceGroup& source,
                  NOX::CopyType type = NOX::DeepCopy);

    //! Destructor
    virtual ~FiniteDifferenceGroup();

    /*!
     * \brief Computes the derivative \f$\partial Jn/\partial p\f$
     */
    /*!
     * The calculation is implemented by calling the corresponding
     * LOCA::DerivUtils::computeDJnDp() method of the passed
     * LOCA::DerivUtils object.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDJnDpMulti(const std::vector<int>& paramIDs,
              const NOX::Abstract::Vector& nullVector,
              NOX::Abstract::MultiVector& result,
              bool isValid);

    /*!
     * \brief Computes the directional derivative
     * \f$\frac{\partial Jn}{\partial x} a\f$ for the given direction
     * \f$ a\f$.
     */
    /*!
     * The calculation is implemented by calling the corresponding
     * LOCA::DerivUtils::computeDJnDxa() method of the passed
     * LOCA::DerivUtils object.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDJnDxaMulti(const NOX::Abstract::Vector& nullVector,
               const NOX::Abstract::MultiVector& aVector,
               NOX::Abstract::MultiVector& result);

    /*!
     * \brief Computes the directional derivative
     * \f$\frac{\partial Jn}{\partial x} a\f$ for the given direction
     * \f$ a\f$.
     */
    /*!
     * The calculation is implemented by calling the corresponding
     * LOCA::DerivUtils::computeDJnDxa() method of the passed
     * LOCA::DerivUtils object.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDJnDxaMulti(const NOX::Abstract::Vector& nullVector,
               const NOX::Abstract::Vector& JnVector,
               const NOX::Abstract::MultiVector& aVector,
               NOX::Abstract::MultiVector& result);

    /*!
     * \brief Computes the derivative
     * \f$\frac{\partial w^TJn}{\partial x}\f$.
     */
    /*!
     * The calculation is implemented by calling the corresponding
     * LOCA::DerivUtils::computeDwtJnDx() method of the passed
     * LOCA::DerivUtils object.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDwtJnDxMulti(const NOX::Abstract::MultiVector& w,
                const NOX::Abstract::Vector& nullVector,
                NOX::Abstract::MultiVector& result);

    //! Bring NOX::Abstract::Group::operator=() into scope
    using NOX::Abstract::Group::operator=;

      }; // class FiniteDifferenceGroup
    } // namespace MooreSpence
  } // namespace TurningPoint
} // namespace LOCA

#endif
