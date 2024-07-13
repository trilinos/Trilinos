// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_TURNINGPOINT_MINIMALLYAUGMENTED_FINITEDIFFERENCEGROUP_H
#define LOCA_TURNINGPOINT_MINIMALLYAUGMENTED_FINITEDIFFERENCEGROUP_H

#include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"   // base class
#include "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H" // base class

namespace LOCA {

  namespace TurningPoint {

    namespace MinimallyAugmented {

      /*!
       * \brief Concrete class that provides concrete implementations of the
       * derivative computation methods of the
       * LOCA::TurningPoint::MinimallyAugmented::AbstractGroup using
       * first-order finite differencing.
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
    public virtual LOCA::TurningPoint::MinimallyAugmented::AbstractGroup,
    public virtual LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup {

      public:

    //! Constructor.
    FiniteDifferenceGroup();

    //! Copy constructor
    FiniteDifferenceGroup(const FiniteDifferenceGroup& source,
                  NOX::CopyType type = NOX::DeepCopy);

    //! Destructor
    virtual ~FiniteDifferenceGroup();

    /*!
     * \brief Computes the derivative \f$\partial w^TJn/\partial p\f$
     */
    /*!
     * The calculation is implemented by calling the corresponding
     * LOCA::DerivUtils::computeDwtJnDp() method of the passed
     * LOCA::DerivUtils object.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDwtJnDp(const std::vector<int>& paramIDs,
               const NOX::Abstract::Vector& w,
               const NOX::Abstract::Vector& nullVector,
               NOX::Abstract::MultiVector::DenseMatrix& result,
               bool isValid);

    /*!
     * \brief Computes the derivative \f$\partial w^TJ/\partial p\f$
     */
    /*!
     * The calculation is implemented by calling the corresponding
     * LOCA::DerivUtils::computeDwtJDp() method of the passed
     * LOCA::DerivUtils object.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDwtJDp(const std::vector<int>& paramIDs,
              const NOX::Abstract::Vector& w,
              NOX::Abstract::MultiVector& result,
              bool isValid);

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
    computeDwtJnDx(const NOX::Abstract::Vector& w,
               const NOX::Abstract::Vector& nullVector,
               NOX::Abstract::Vector& result);

    //! Bring NOX::Abstract::Group::operator=() into scope
    using NOX::Abstract::Group::operator=;

      }; // class FiniteDifferenceGroup
    } // namespace MinimallyAugmented
  } // namespace TurningPoint
} // namespace LOCA

#endif
