// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_TIMEDEPENDENT_ABSTRACTGROUP_H
#define LOCA_TIMEDEPENDENT_ABSTRACTGROUP_H

#include "LOCA_MultiContinuation_AbstractGroup.H"    // base class

namespace LOCA {

  namespace TimeDependent {

    /*!
     * \brief Interface to underlying groups for time dependent systems
     */
    /*!
     * This abstract class provides an interface for time dependent problems,
     * i.e., problems with a mass matrix (typically used in eignvalue or
     * Hopf calculations).  It provides pure virtual methods for computing
     * and manipulating the shifted matrix \f$\alpha J + \beta M\f$ where
     * \f$ J\f$ is the Jacobian matrix and \f$ M\f$ is the mass matrix.
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

      //! Compute the shifted matrix
      virtual NOX::Abstract::Group::ReturnType
      computeShiftedMatrix(double alpha, double beta) = 0;

      //! Multiply the shifted matrix by a vector.
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrix(const NOX::Abstract::Vector& input,
             NOX::Abstract::Vector& result) const = 0;

      //! Multiply the shifted matrix by a multi-vector.
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixMultiVector(
                const NOX::Abstract::MultiVector& input,
                NOX::Abstract::MultiVector& result) const = 0;

      /*!
       * \brief Apply the inverse of the shifted matrix by a multi-vector, as
       * needed by the shift-and-invert and generalized Cayley transformations.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixInverseMultiVector(
                    Teuchos::ParameterList& params,
                const NOX::Abstract::MultiVector& input,
                NOX::Abstract::MultiVector& result) const = 0;


      //@}

      //! Compute the second shifted matrix. Can avoid recomputing if two are stored.
      /*!
       * Implementation here prints an error message and returns
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeSecondShiftedMatrix(double alpha, double beta) = 0;

      //! Multiply the shifted matrix by a vector.
      /*!
       * Implementation here prints an error message and returns
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applySecondShiftedMatrix(const NOX::Abstract::Vector& input,
                         NOX::Abstract::Vector& result) const = 0;

      //! Multiply the shifted matrix by a multi-vector.
      /*!
       * Implementation here prints an error message and returns
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applySecondShiftedMatrixMultiVector(
                                const NOX::Abstract::MultiVector& input,
                                NOX::Abstract::MultiVector& result) const = 0;

      //@}


      //! Bring NOX::Abstract::Group::operator=() into scope
      using NOX::Abstract::Group::operator=;

    };
  } // namespace TimeDependent
} // namespace LOCA

#endif
