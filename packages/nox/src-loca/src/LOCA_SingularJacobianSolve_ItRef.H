// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
