// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_SINGULARJACOBIANSOLVE_NIC_H
#define LOCA_SINGULARJACOBIANSOLVE_NIC_H

#include "LOCA_SingularJacobianSolve_Generic.H"  // base class

namespace LOCA {

  namespace SingularJacobianSolve {

    /*!
     * \brief This class computes the solution to \f$J x = b\f$ using the
     * %Nic method.
     */
    /*!
     * The idea here is to use deflation of the right hand side to improve
     * the conditioning of the linear system.  Typically a solution to
     * \f$J x = b\f$ when \f$J\f$ is nearly singular will have a large
     * component in the direction of the null vector \f$v\f$.  The idea
     * then is to deflate \f$Jv\f$ out of the right hand side \f$b\f$.
     * The complete algorithm used here is:
     * \f[
     *   \begin{aligned}
     *     &\tilde{b} = b - \frac{v^T b}{v^T J v} Jv \\
     *     &\text{Solve}\; J\tilde{x} = \tilde{b} \\
     *     &x = \tilde{x} + \frac{v^T b}{v^T J v} v
     *   \end{aligned}
     * \f]
     * The solve \f$J\tilde{x} = \tilde{b}\f$ uses the underlying group's
     * applyJacobianInverse() method and therefore this is a generic
     * technique for computing solutions to nearly singular system since
     * it uses any supplied linear solver.

     * This algorithm is selected by
     * setting the "Method" parameter of the "Singular Solve" sublist of
     * the NOX linear solver parameter list to "Nic".  The idea for this
     * algorithm is taken from:  R. A. Nicolaides, "Deflation of Conjugate
     * Gradients With Applications to Boundary Value Problems," SIAM J.
     * Numer. Anal., 24(2), 1987.
     */
    class Nic : public LOCA::SingularJacobianSolve::Generic {

    public:

      //! Constructor.
      Nic(Teuchos::ParameterList& params);

      //! Copy constructor
      Nic(const Nic& source);

      //! Destructor
      virtual ~Nic();

      //! Clone function
      virtual LOCA::SingularJacobianSolve::Generic* clone() const;

      //! Assignment operator
      virtual LOCA::SingularJacobianSolve::Generic&
      operator = (const LOCA::SingularJacobianSolve::Generic& source);

      //! Assignment operator
      virtual Nic&
      operator = (const Nic& source);

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
