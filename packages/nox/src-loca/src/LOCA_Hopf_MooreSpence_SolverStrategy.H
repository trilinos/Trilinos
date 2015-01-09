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

#ifndef LOCA_HOPF_MOORESPENCE_SOLVERSTRATEGY_H
#define LOCA_HOPF_MOORESPENCE_SOLVERSTRATEGY_H

#include "Teuchos_RCP.hpp"

#include "NOX_Abstract_Group.H"  // for ReturnType

// forward declarations
namespace LOCA {
  namespace Hopf {
    namespace MooreSpence {
      class AbstractGroup;
      class ExtendedGroup;
      class ExtendedMultiVector;
    }
  }
}

namespace LOCA {

  namespace Hopf {

    namespace MooreSpence {

      /*!
       * \brief Abstract strategy for solving the Moore-Spence Hopf
       * equations.
       */
      /*!
       * This class provides an abstract interface for solver strategies to
       * solve the Moore-Spence Hopf Newton system:
       * \f[
       *     \begin{bmatrix}
       *          J       &    0   &    0   &   0  &    f_p     \\
       *       (Jy-wBz)_x &    J   &  -wB   & -Bz  & (Jy-wBz)_p \\
       *       (Jz+wBy)_x &   wB   &    J   &  By  & (Jz+wBy)_p \\
       *          0       & \phi^T &    0   &   0  &     0      \\
       *          0       &    0   & \phi^T &   0  &     0
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        X \\
       *        Y \\
       *        Z \\
       *        \Omega^T \\
       *        \Lambda^T
       *     \end{bmatrix} =
       *     \begin{bmatrix}
       *        F \\
       *        G \\
       *        H \\
       *        U^T \\
       *        V^T
       *     \end{bmatrix}.
       *  \f]
       * After instantiating a solver (via
       * LOCA::Hopf::MooreSpence::SolverFactory), the linear
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
     * \param hopfGroup [in] Hopf group representing the Hopf equations.
     * \param yVector [in] Vector representing y
     * \param zVector [in] Vector representing z
     * \param CeRealVector [in] Vector representing Jy-wBz
     * \param CeImagVector [in] Vector representing Jz+wBy
     * \param dfdp [in] Vector representing df/dp
     * \param dCedpReal [in] Vector representing d(Jy-wBz)/dp
     * \param dCedpImag [in] Vector representing d(Jz+wBy)/dp
     * \param ByVector [in] Vector representing By
     * \param mBzVector [in] Vector representing -Bz
     * \param w [in] Bifurcation frequency w
     */
    virtual void setBlocks(
     const Teuchos::RCP<LOCA::Hopf::MooreSpence::AbstractGroup>& group,
     const Teuchos::RCP<LOCA::Hopf::MooreSpence::ExtendedGroup>& hopfGroup,
     const Teuchos::RCP<const NOX::Abstract::Vector>& yVector,
     const Teuchos::RCP<const NOX::Abstract::Vector>& zVector,
     const Teuchos::RCP<const NOX::Abstract::Vector>& CeRealVector,
     const Teuchos::RCP<const NOX::Abstract::Vector>& CeImagVector,
     const Teuchos::RCP<const NOX::Abstract::Vector>& dfdp,
     const Teuchos::RCP<const NOX::Abstract::Vector>& dCedpReal,
     const Teuchos::RCP<const NOX::Abstract::Vector>& dCedpImag,
     const Teuchos::RCP<const NOX::Abstract::Vector>& ByVector,
     const Teuchos::RCP<const NOX::Abstract::Vector>& mBzVector,
     double w) = 0;

    /*!
     * \brief Solves the extended system as defined above
     */
    /*!
     * The \em params argument is the linear solver parameters.
     */
    virtual NOX::Abstract::Group::ReturnType
    solve(
      Teuchos::ParameterList& params,
      const LOCA::Hopf::MooreSpence::ExtendedMultiVector& input,
      LOCA::Hopf::MooreSpence::ExtendedMultiVector& result) const = 0;

      private:

    //! Private to prohibit copying
    SolverStrategy(const SolverStrategy&);

    //! Private to prohibit copying
    SolverStrategy& operator = (const SolverStrategy&);

      }; // Class SolverStrategy

    } // Namespace MooreSpence

  } // Namespace Hopf

} // Namespace LOCA

#endif // LOCA_HOPF_MOORESPENCE_SOLVERSTRATEGY_H
