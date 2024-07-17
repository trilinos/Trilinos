// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_SINGULARJACOBIANSOLVE_GENERIC_H
#define LOCA_SINGULARJACOBIANSOLVE_GENERIC_H

#include "NOX_Abstract_Group.H"  // for ReturnType, Teuchos::ParameterList

// Forward declarations
namespace LOCA {
  namespace Continuation {
    class AbstractGroup;
  }
}

namespace LOCA {

  //! %Strategies for solving \f$Jx=b\f$ when \f$J\f$ is (nearly) singular
  namespace SingularJacobianSolve {

    //! %Generic singular jacobian solve interface.
    /*!
     * %Generic interface for solving \f$Jx=b\f$ when \f$J\f$ is (nearly)
     * singular.  All classes the implement a method for computing solutions
     * to nearly singular systems should be derived from this class.
    */
    class Generic {

    public:

      //! Constructor.
      /*!
    Constructors of derived objects should look like reset.
      */
      Generic() {};

      //! Copy constructor
      Generic(const Generic& source) {}

      //! Destructor
      virtual ~Generic() {};

      //! Clone function
      virtual Generic* clone() const = 0;

      //! Assignment operator
      virtual Generic& operator = (const Generic& source) = 0;

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType
      reset(Teuchos::ParameterList& params) = 0;

      /*!
    \brief Compute solution to singular system.

    The passed parameters are assumed be the (nonsingular) linear solver
    parameters.
      */
      virtual NOX::Abstract::Group::ReturnType
      compute(Teuchos::ParameterList& params,
          LOCA::Continuation::AbstractGroup& grp,
          const NOX::Abstract::Vector& input,
          const NOX::Abstract::Vector& approxNullVec,
          const NOX::Abstract::Vector& jacApproxNullVec,
          NOX::Abstract::Vector& result) = 0;

      /*!
    \brief Compute solution to singular system with multiple RHS.

    The passed parameters are assumed be the (nonsingular) linear solver
    parameters.
      */
      virtual NOX::Abstract::Group::ReturnType
      computeMulti(Teuchos::ParameterList& params,
           LOCA::Continuation::AbstractGroup& grp,
           const NOX::Abstract::Vector*const* inputs,
           const NOX::Abstract::Vector& approxNullVec,
           const NOX::Abstract::Vector& jacApproxNullVec,
           NOX::Abstract::Vector** results,
           int nVecs) = 0;

    };
  } // namespace SingularJacobianSolve
} // namespace LOCA

#endif
