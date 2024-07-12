// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_HOPF_MOORESPENCE_ABSTRACTGROUP_H
#define LOCA_HOPF_MOORESPENCE_ABSTRACTGROUP_H

#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"    // base class
#include "LOCA_TimeDependent_AbstractGroup.H"               // base class

namespace LOCA {

  /*!
   * \brief Groups and vectors for Hopf bifurcations
   */
  namespace Hopf {

    /*!
     * \brief Groups and vectors for locating Hopf bifurcations using the
     * Moore-Spence formulation.
     */
    namespace MooreSpence {

      /*!
       * \brief Interface to underlying groups for Hopf point calculations
       * using the Moore-Spence formulation.
       */
      /*!
       * This abstract class provides the required interface for underlying
       * groups to locate Hopf bifurcations using the bordering
       * algorithm for the Moore-Spence Hopf fomulation (see
       * LOCA::Hopf::MooreSpence::ExtendedGroup for a description of the
       * governing equations).
       *
       * This class is derived from the
       * LOCA::TurningPoint::MooreSpence::AbstractGroup and
       * LOCa::TimeDependent::AbstractGroup and declares several
       * pure virtual methods for applying, solving, and computing derivatives
       * of the complex matrix \f$ J+i\omega B\f$ where \f$ J\f$ is the Jacobian
       * matrix, \f$ B\f$ is the mass matrix, and \f$\omega\f$ is a (real)
       * scalar.
       */
      class AbstractGroup :
    public virtual LOCA::TurningPoint::MooreSpence::AbstractGroup,
    public virtual LOCA::TimeDependent::AbstractGroup {

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

    //! Is  \f$ J+i\omega B\f$ valid
    virtual bool isComplex() const = 0;

    //! Compute \f$ J+i\omega B\f$
    /*!
     * The argument \b frequency stores \f$\omega\f$.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeComplex(double frequency) = 0;

    //! Compute \f$(J+i\omega B)(y+iz)\f$
    virtual NOX::Abstract::Group::ReturnType
    applyComplex(const NOX::Abstract::Vector& input_real,
             const NOX::Abstract::Vector& input_imag,
             NOX::Abstract::Vector& result_real,
             NOX::Abstract::Vector& result_imag) const = 0;

    //! Compute \f$(J+i\omega B)(y+iz)\f$
    virtual NOX::Abstract::Group::ReturnType
    applyComplexMultiVector(
                const NOX::Abstract::MultiVector& input_real,
                const NOX::Abstract::MultiVector& input_imag,
                NOX::Abstract::MultiVector& result_real,
                NOX::Abstract::MultiVector& result_imag) const = 0;

    //! Solve \f$(J+i\omega B)(y+iz) = a+ib\f$
    virtual NOX::Abstract::Group::ReturnType
    applyComplexInverseMultiVector(
                Teuchos::ParameterList& params,
                const NOX::Abstract::MultiVector& input_real,
                const NOX::Abstract::MultiVector& input_imag,
                NOX::Abstract::MultiVector& result_real,
                NOX::Abstract::MultiVector& result_imag) const = 0;

    /*!
     * \brief Computes the derivative
     * \f$\frac{\partial (J+i\omega B)(y+iz)}{\partial p}\f$
     * where \f$ p\f$ is the parameter indexed by \b paramIDs.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDCeDp(const std::vector<int>& paramIDs,
             const NOX::Abstract::Vector& yVector,
             const NOX::Abstract::Vector& zVector,
             double w,
             NOX::Abstract::MultiVector& result_real,
             NOX::Abstract::MultiVector& result_imag,
             bool isValid) = 0;

    /*!
     * \brief Computes the directional derivative
     * \f$\frac{\partial (J+i\omega B)(y+iz)}{\partial x} a\f$ for the
     * given direction \f$ a\f$.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDCeDxa(const NOX::Abstract::Vector& yVector,
              const NOX::Abstract::Vector& zVector,
              double w,
              const NOX::Abstract::MultiVector& aVector,
              NOX::Abstract::MultiVector& result_real,
              NOX::Abstract::MultiVector& result_imag) = 0;

    /*!
     * \brief Computes the directional derivative
     * \f$\frac{\partial (J+i\omega B)(y+iz)}{\partial x} a\f$ for the
     * given direction \f$ a\f$. The
     * arguments \b Ce_real and \b Ce_imag hold the real and imaginary
     * components of \f$(J+i\omega B)(y+iz)\f$.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDCeDxa(const NOX::Abstract::Vector& yVector,
              const NOX::Abstract::Vector& zVector,
              double w,
              const NOX::Abstract::MultiVector& aVector,
              const NOX::Abstract::Vector& Ce_real,
              const NOX::Abstract::Vector& Ce_imag,
              NOX::Abstract::MultiVector& result_real,
              NOX::Abstract::MultiVector& result_imag) = 0;

    //@}

    //! Bring NOX::Abstract::Group::operator=() into scope
    using NOX::Abstract::Group::operator=;

      };
    } // namespace MooreSpence
  } // namespace Hopf
} // namespace LOCA

#endif
