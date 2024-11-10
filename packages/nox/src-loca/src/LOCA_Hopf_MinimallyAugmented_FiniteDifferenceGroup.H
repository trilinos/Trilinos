// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_HOPF_MINIMALLYAUGMENTED_FINITEDIFFERENCEGROUP_H
#define LOCA_HOPF_MINIMALLYAUGMENTED_FINITEDIFFERENCEGROUP_H

#include "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"   // base class
#include "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"  // base class

namespace LOCA {

  namespace Hopf {

    namespace MinimallyAugmented {

      /*!
       * \brief Concrete class that provides concrete implementations of the
       * derivative computation methods of the
       * LOCA::Hopf::MinimallyAugmented::AbstractGroup using
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
    public virtual LOCA::Hopf::MinimallyAugmented::AbstractGroup,
    public virtual LOCA::Hopf::MooreSpence::FiniteDifferenceGroup {

      public:

    //! Constructor.
    FiniteDifferenceGroup();

    //! Copy constructor
    FiniteDifferenceGroup(const FiniteDifferenceGroup& source,
                  NOX::CopyType type = NOX::DeepCopy);

    //! Destructor
    virtual ~FiniteDifferenceGroup();

    /*!
     * \brief Computes the derivative \f$\partial w^TCe/\partial p\f$
     */
    /*!
     * The calculation is implemented by calling the corresponding
     * LOCA::DerivUtils::computeDwtCeDp() method of the passed
     * LOCA::DerivUtils object.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDwtCeDp(const std::vector<int>& paramIDs,
               const NOX::Abstract::Vector& w1,
               const NOX::Abstract::Vector& w2,
               const NOX::Abstract::Vector& y,
               const NOX::Abstract::Vector& x,
               double omega,
               NOX::Abstract::MultiVector::DenseMatrix& result_real,
               NOX::Abstract::MultiVector::DenseMatrix& result_imag,
               bool isValid);

    /*!
     * \brief Computes the derivative
     * \f$\frac{\partial w^TCe}{\partial x}\f$.
     */
    /*!
     * The calculation is implemented by calling the corresponding
     * LOCA::DerivUtils::computeDwtCeDx() method of the passed
     * LOCA::DerivUtils object.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDwtCeDx(const NOX::Abstract::Vector& w1,
               const NOX::Abstract::Vector& w2,
               const NOX::Abstract::Vector& y,
               const NOX::Abstract::Vector& z,
               double omega,
               NOX::Abstract::Vector& result_real,
               NOX::Abstract::Vector& result_imag);

    //! Bring NOX::Abstract::Group::operator=() into scope
    using NOX::Abstract::Group::operator=;

      }; // class FiniteDifferenceGroup
    } // namespace MinimallyAugmented
  } // namespace Hopf
} // namespace LOCA

#endif
