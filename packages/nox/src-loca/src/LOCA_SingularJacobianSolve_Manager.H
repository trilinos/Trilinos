// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_SINGULARJACOBIANSOLVE_MANAGER_H
#define LOCA_SINGULARJACOBIANSOLVE_MANAGER_H

#include "LOCA_SingularJacobianSolve_Generic.H"  // base class
#include "Teuchos_ParameterList.hpp"

namespace LOCA {

  namespace SingularJacobianSolve {

    //! %Manager for all singular Jacobian solve computations
    /*!
      The parameters passed to the constructor or reset should specify
      the "Method", as described below, as well as any additional
      parameters for that particular method.

      <ul>
      <li> "Method" - Name of the singular jacobian solve method.
           Valid choices are
      <ul>
      <li> "Default" (LOCA::SingularJacobianSolve::Default) [Default]
      <li> "Nic" (LOCA::SingularJacobianSolve::Nic)
      <li> "Nic-Day" (LOCA::SingularJacobianSolve::NicDay)
      <li> "Iterative Refinement" (LOCA::SingularJacobianSolve::ItRef)
      </ul>
      </ul>
    */
    class Manager : public LOCA::SingularJacobianSolve::Generic {

    public:

      //! Constructor.
      Manager(Teuchos::ParameterList& params);

      //! Constructor.
      Manager(const Teuchos::ParameterList& params = Teuchos::ParameterList());

      //! Copy constructor
      Manager(const Manager& source);

      //! Destructor
      virtual ~Manager();

      //! Clone function
      virtual LOCA::SingularJacobianSolve::Generic* clone() const;

      //! Assignment operator
      virtual LOCA::SingularJacobianSolve::Generic&
      operator = (const LOCA::SingularJacobianSolve::Generic& source);

      //! Assignment operator
      virtual Manager&
      operator = (const Manager& source);

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType
      reset(Teuchos::ParameterList& params);

      //! Computes solution based on method parameter
      virtual NOX::Abstract::Group::ReturnType
      compute(Teuchos::ParameterList& params,
          LOCA::Continuation::AbstractGroup& grp,
          const NOX::Abstract::Vector& input,
          const NOX::Abstract::Vector& approxNullVec,
          const NOX::Abstract::Vector& jacApproxNullVec,
          NOX::Abstract::Vector& result);

      //! Computes solution based on method parameter for multiple RHS
      virtual NOX::Abstract::Group::ReturnType
      computeMulti(Teuchos::ParameterList& params,
           LOCA::Continuation::AbstractGroup& grp,
           const NOX::Abstract::Vector*const* inputs,
           const NOX::Abstract::Vector& approxNullVec,
           const NOX::Abstract::Vector& jacApproxNullVec,
           NOX::Abstract::Vector** results,
           int nVecs);

      //! Returns current method
      const std::string& getMethod() const;

    protected:
      //! Stores predictor method name
      std::string method;

      //! Stores pointer to current singular solver object
      LOCA::SingularJacobianSolve::Generic *singularSolverPtr;

    };
  } // namespace SingularJacobianSolve
} // namespace LOCA

#endif
