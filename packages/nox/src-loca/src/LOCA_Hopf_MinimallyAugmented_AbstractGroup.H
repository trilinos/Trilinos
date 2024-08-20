// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_HOPF_MINIMALLYAUGMENTED_ABSTRACTGROUP_H
#define LOCA_HOPF_MINIMALLYAUGMENTED_ABSTRACTGROUP_H

#include "LOCA_Hopf_MooreSpence_AbstractGroup.H"                // base class

namespace LOCA {

  namespace Hopf {

    /*!
     * \brief Groups and vectors for locating Hopf bifurcations
     * using the minimally augmented Hopf formulation.
     */
    namespace MinimallyAugmented {

      /*!
       * \brief Interface to underlying groups for Hopf calculations
       * using the minimally augmented formulation.
       */
      /*!
       * This abstract class provides the required interface for underlying
       * groups to locate Hopfs using the minimally augmented Hopf (see
       * LOCA::Hopf::MinimallyAugmented::ExtendedGroup for a
       * description of the governing equations).
       *
       * This class is derived from the
       * LOCA::Hopf::MooreSpence::AbstractGroup
       * and declares several pure virtual methods to compute various
       * derivatives of \f$ w^H C e\f$ for a given \f$ w = w_1 + i w_2\f$ and
       * \f$ e = y + i z\f$ where \f$ C = J + i \omega M\f$.  Default
       * implementations for the derivatives using finite differencing are
       * implemented in the
       * LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup.
       */
      class AbstractGroup :
    public virtual LOCA::Hopf::MooreSpence::AbstractGroup {

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
     * Computes conjugate-tranpose matrix vector product
     * \f$ (J+i\omega M)^H (x + iy) \f$.
     */
    virtual NOX::Abstract::Group::ReturnType
    applyComplexTranspose(const NOX::Abstract::Vector& input_real,
                  const NOX::Abstract::Vector& input_imag,
                  NOX::Abstract::Vector& result_real,
                  NOX::Abstract::Vector& result_imag) const = 0;

    /*!
     * Computes conjugate-tranpose matrix vector product
     * \f$ (J+i\omega M)^H (x + iy) \f$.
     */
    virtual NOX::Abstract::Group::ReturnType
    applyComplexTransposeMultiVector(
               const NOX::Abstract::MultiVector& input_real,
               const NOX::Abstract::MultiVector& input_imag,
               NOX::Abstract::MultiVector& result_real,
               NOX::Abstract::MultiVector& result_imag) const = 0;

    //! Solve \f$(J+i\omega M)^H (x + iy) = a+ib\f$
    virtual NOX::Abstract::Group::ReturnType
    applyComplexTransposeInverseMultiVector(
                Teuchos::ParameterList& params,
                const NOX::Abstract::MultiVector& input_real,
                const NOX::Abstract::MultiVector& input_imag,
                NOX::Abstract::MultiVector& result_real,
                NOX::Abstract::MultiVector& result_imag) const = 0;

    /*!
     * \brief Computes the derivative \f$\partial w^H C e/\partial p\f$,
     * \f$ w = w_1 + i w_2\f$, \f$ e = y + i z\f$, \f$ C = J + i \omega M\f$.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDwtCeDp(const std::vector<int>& paramIDs,
               const NOX::Abstract::Vector& w1,
               const NOX::Abstract::Vector& w2,
               const NOX::Abstract::Vector& y,
               const NOX::Abstract::Vector& z,
               double omega,
               NOX::Abstract::MultiVector::DenseMatrix& result_real,
               NOX::Abstract::MultiVector::DenseMatrix& result_imag,
               bool isValid) = 0;

    /*!
     * \brief Computes the derivative
     * \f$\frac{\partial w^H C e}{\partial x}\f$,
     * \f$ w = w_1 + i w_2\f$, \f$ e = y + i z\f$, \f$ C = J + i \omega M\f$.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDwtCeDx(const NOX::Abstract::Vector& w1,
               const NOX::Abstract::Vector& w2,
               const NOX::Abstract::Vector& y,
               const NOX::Abstract::Vector& z,
               double omega,
               NOX::Abstract::Vector& result_real,
               NOX::Abstract::Vector& result_imag) = 0;

    //@}

    //! Bring NOX::Abstract::Group::operator=() into scope
    using NOX::Abstract::Group::operator=;

      }; // class AbstractGroup

    } // namespace MinimallyAugmented

  } // namespace Hopf

} // namespace LOCA

#endif
