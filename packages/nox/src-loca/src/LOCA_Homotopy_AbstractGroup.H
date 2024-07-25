// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_HOMOTOPY_ABSTRACTGROUP_H
#define LOCA_HOMOTOPY_ABSTRACTGROUP_H

#include "LOCA_MultiContinuation_AbstractGroup.H"    // base class

namespace LOCA {

  namespace Homotopy {

    /*!
     * \brief Interface to underlying groups for homotopy calculations
     */
    /*!
     * This abstract class provides an interface for a homotopy technique
     * for solving nonlinear equations.  See LOCA::Homotopy::Group for a
     * description of the technique used.  This class provides a single
     * pure virtual method, augmentJacobianForHomotopy(), which scales the
     * diagonal of the Jacobian by a constant times the identity matrix.
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
       * \brief Replace Jacobian \f$ J\f$ by \f$ aJ+bI\f$ where \f$ I\f$ is
       * the identity matrix
       */
      virtual NOX::Abstract::Group::ReturnType
      augmentJacobianForHomotopy(double a, double b) = 0;

      //@}

      //! Bring NOX::Abstract::Group::operator=() into scope
      using NOX::Abstract::Group::operator=;

    };
  } // namespace Homotopy
} // namespace LOCA

#endif
