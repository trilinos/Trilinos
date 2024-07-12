// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_SINGULARJACOBIANSOLVE_DEFAULT_H
#define LOCA_SINGULARJACOBIANSOLVE_DEFAULT_H

#include "LOCA_SingularJacobianSolve_Generic.H"  // base class

namespace LOCA {

  namespace SingularJacobianSolve {

    //! %Default singular Jacobian solve computation class
    /*!
     * This class computes the solution to \f$J x = b\f$ using the
     * applyJacobianInverse method of the underlying group ignoring the
     * null vector data.
    */
    class Default : public LOCA::SingularJacobianSolve::Generic {

    public:

      //! Constructor.
      Default(Teuchos::ParameterList& params);

      //! Copy constructor
      Default(const Default& source);

      //! Destructor
      virtual ~Default();

      //! Clone function
      virtual LOCA::SingularJacobianSolve::Generic* clone() const;

      //! Assignment operator
      virtual LOCA::SingularJacobianSolve::Generic&
      operator = (const LOCA::SingularJacobianSolve::Generic& source);

      //! Assignment operator
      virtual Default&
      operator = (const Default& source);

      //! Reset parameters
      /*!
    There are no additional parameters for the default calculation.
      */
      virtual NOX::Abstract::Group::ReturnType
      reset(Teuchos::ParameterList& params);

      //! Computes the solution as described above.
      virtual NOX::Abstract::Group::ReturnType
      compute(Teuchos::ParameterList& params,
          LOCA::Continuation::AbstractGroup& grp,
          const NOX::Abstract::Vector& input,
          const NOX::Abstract::Vector& approxNullVec,
          const NOX::Abstract::Vector& jacApproxNullVec,
          NOX::Abstract::Vector& result);

      //! Computes solution for multiple RHS using applyJacobianInverseMulti
      virtual NOX::Abstract::Group::ReturnType
      computeMulti(Teuchos::ParameterList& params,
           LOCA::Continuation::AbstractGroup& grp,
           const NOX::Abstract::Vector*const* inputs,
           const NOX::Abstract::Vector& approxNullVec,
           const NOX::Abstract::Vector& jacApproxNullVec,
           NOX::Abstract::Vector** results,
           int nVecs);

    };
  } // namespace SingularJacobianSolve
} // namespace LOCA

#endif
