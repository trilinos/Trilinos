// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_TURNINGPOINT_MINIMALLYAUGMENTED_ABSTRACTGROUP_H
#define LOCA_TURNINGPOINT_MINIMALLYAUGMENTED_ABSTRACTGROUP_H

#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"    // base class

namespace LOCA {

  namespace TurningPoint {

    /*!
     * \brief Groups and vectors for locating turning point bifurcations
     * using the minimally augmented turning point formulation.
     */
    namespace MinimallyAugmented {

      /*!
       * \brief Interface to underlying groups for turning point calculations
       * using the minimally augmented formulation.
       */
      /*!
       * This abstract class provides the required interface for underlying
       * groups to locate turning points using the minimally augmented turning
       * point formulation (see
       * LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup for a
       * description of the governing equations).
       *
       * This class is derived from the
       * LOCA::TurningPoint::MooreSpence::AbstractGroup
       * and declares several pure virtual methods compute various
       * derivatives of \f$ w^TJn\f$ for a given \f$ w\f$ and \f$ n\f$.  Default
       * implementations for the derivatives using finite differencing are
       * implemented in the
       * LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup.
       */
      class AbstractGroup :
    public virtual LOCA::TurningPoint::MooreSpence::AbstractGroup {

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
     * \brief Computes the derivative \f$\partial w^TJn/\partial p\f$
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDwtJnDp(const std::vector<int>& paramIDs,
               const NOX::Abstract::Vector& w,
               const NOX::Abstract::Vector& nullVector,
               NOX::Abstract::MultiVector::DenseMatrix& result,
               bool isValid) = 0;

    /*!
     * \brief Computes the derivative \f$\partial w^TJ/\partial p\f$
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDwtJDp(const std::vector<int>& paramIDs,
              const NOX::Abstract::Vector& w,
              NOX::Abstract::MultiVector& result,
              bool isValid) = 0;

    /*!
     * \brief Computes the derivative
     * \f$\frac{\partial w^TJn}{\partial x}\f$.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDwtJnDx(const NOX::Abstract::Vector& w,
               const NOX::Abstract::Vector& nullVector,
               NOX::Abstract::Vector& result) = 0;

    //@}

    //! Bring NOX::Abstract::Group::operator=() into scope
    using NOX::Abstract::Group::operator=;

      }; // class AbstractGroup

    } // namespace MinimallyAugmented

  } // namespace TurningPoint

} // namespace LOCA

#endif
