// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_SINGULARJACOBIANSOLVE_ITREF_H
#define LOCA_SINGULARJACOBIANSOLVE_ITREF_H

#include "LOCA_SingularJacobianSolve_Generic.H"  // base class

namespace LOCA {

  namespace SingularJacobianSolve {

    /*!
     * \brief This class computes the solution to \f$J x = b\f$ using
     * one step of iterative refinement.
     */
    /*!
     * This singular solve method uses one step of iterative refinement to
     * improve the accuracy of the solution to the linear system
     * \f$J x = b\f$.  In particular, the algorithm used here is
      * \f[
     *   \begin{aligned}
     *     &\text{Solve}\; Jx_1 = b \\
     *     &r = b - Jx_1 \\
     *     &\text{Solve}\; Jx_2 = r \\
     *     &x = x_1 + x_2
     *   \end{aligned}
     * \f]
     * Both solves use the underlying group's applyJacobianInverse() method
     * and therefore this is a generic technique for computing solutions to
     * nearly singular system since it uses any supplied linear solver.

     * This algorithm is selected by setting the "Method" parameter of the
     * "Singular Solve" sublist of the NOX linear solver parameter list to
     * "Iterative Refinement".
     */
    class ItRef : public LOCA::SingularJacobianSolve::Generic {

    public:

      //! Constructor.
      ItRef(Teuchos::ParameterList& params);

      //! Copy constructor
      ItRef(const ItRef& source);

      //! Destructor
      virtual ~ItRef();

      //! Clone function
      virtual LOCA::SingularJacobianSolve::Generic* clone() const;

      //! Assignment operator
      virtual LOCA::SingularJacobianSolve::Generic&
      operator = (const LOCA::SingularJacobianSolve::Generic& source);

      //! Assignment operator
      virtual ItRef&
      operator = (const ItRef& source);

      //! Reset parameters
      /*!
    There are no additional parameters for the Nic calculation.
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

      //! Computes solution for multiple RHS
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
