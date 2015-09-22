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

#ifndef LOCA_PITCHFORK_MOORESPENCE_SOLVERSTRATEGY_H
#define LOCA_PITCHFORK_MOORESPENCE_SOLVERSTRATEGY_H

#include "Teuchos_RCP.hpp"

#include "NOX_Abstract_Group.H"  // for ReturnType

// forward declarations
namespace LOCA {
  namespace Pitchfork {
    namespace MooreSpence {
      class AbstractGroup;
      class ExtendedGroup;
      class ExtendedMultiVector;
    }
  }
}

namespace LOCA {

  namespace Pitchfork {

    namespace MooreSpence {

      /*!
       * \brief Abstract strategy for solving the Moore-Spence pitchfork
       * equations.
       */
      /*!
       * This class provides an abstract interface for solver strategies to
       * solve the Moore-Spence pitchfork Newton system:
       * \f[
       *     \begin{bmatrix}
       *          J   &    0   & \psi &  f_p   \\
       *       (Jv)_x &    J   &   0  & (Jv)_p \\
       *       \psi^T &    0   &   0  &   0    \\
       *          0   & \phi^T &   0  &   0
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        X \\
       *        Y \\
       *        w \\
       *        z
       *     \end{bmatrix} =
       *     \begin{bmatrix}
       *        F \\
       *        G \\
       *        s \\
       *        h
       *     \end{bmatrix}.
       *  \f]
       * After instantiating a solver Solvers (via
       * LOCA::Pitchfork::MooreSpence::SolverFactory), the linear
       * system is set up by setBlocks() and can then be solved by solve().
       */
      class SolverStrategy {

      public:

    //! Constructor
    SolverStrategy() {}

    //! Destructor
    virtual ~SolverStrategy() {}

    //! Set blocks in extended linear system
    /*!
     * \param group [in] Underlying group representing J
     * \param pfGroup [in] Pitchfork group representing the pitchfork
     *                     equations.
     * \param asymMultiVector [in] Multivector representing the
     *                             asymmetric vector
     * \param nullVector [in] Vector representing v
     * \param JnVector [in] Vector representing Jv
     * \param dfdp [in] Vector representing df/dp
     * \param dJndp [in] Vector representing d(Jv)/dp
     */
    virtual void setBlocks(
      const Teuchos::RCP<LOCA::Pitchfork::MooreSpence::AbstractGroup>& group,
      const Teuchos::RCP<LOCA::Pitchfork::MooreSpence::ExtendedGroup>& pfGroup,
      const Teuchos::RCP<const NOX::Abstract::MultiVector>& asymMultiVector,
      const Teuchos::RCP<const NOX::Abstract::Vector>& nullVector,
      const Teuchos::RCP<const NOX::Abstract::Vector>& JnVector,
      const Teuchos::RCP<const NOX::Abstract::Vector>& dfdp,
      const Teuchos::RCP<const NOX::Abstract::Vector>& dJndp) = 0;

    /*!
     * \brief Solves the extended system as defined above
     */
    /*!
     * The \em params argument is the linear solver parameters.
     */
    virtual NOX::Abstract::Group::ReturnType
    solve(
      Teuchos::ParameterList& params,
      const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& input,
          LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& result) const = 0;

      private:

    //! Private to prohibit copying
    SolverStrategy(const SolverStrategy&);

    //! Private to prohibit copying
    SolverStrategy& operator = (const SolverStrategy&);

      }; // Class SolverStrategy

    } // Namespace MooreSpence

  } // Namespace TurningPont

} // Namespace LOCA

#endif // LOCA_PITCHFORK_MOORESPENCE_SOLVERSTRATEGY_H
