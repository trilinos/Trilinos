// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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

#ifndef NOX_SOLVER_FACTORY_H
#define NOX_SOLVER_FACTORY_H

#include "Teuchos_RCP.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace NOX {

  namespace Abstract {
    class Group;
  }
  namespace StatusTest {
    class Generic;
  }

  namespace Solver{

    class Generic;

    /*! \brief %Factory class to control the creation of solvers derived from the NOX::Solver::Generic object.

    <B>Parameters</B>

    The following entries may be specified in the parameter list.

    <ul>
    <li> "Nonlinear Solver" <std::string> Name of the solver method. Valid choices are
    <ul>
    <li> "Line Search Based" or "Newton" (NOX::Solver::LineSearchBased) [Default]
    <li> "Trust Region Based" (NOX::Solver::TrustRegionBased)
    <li> "Inexact Trust Region Based" (NOX::Solver::InexactTrustRegionBased)
    <li> "Tensor Based" (NOX::Solver::TensorBased)
    </ul>
    <li> "Line Search Based" <Teuchos::ParameterList> Sublist for the NOX::Solver::LineSearchBased solver.
    <li> "Trust Region Based" <Teuchos::ParameterList> Sublist for the NOX::Solver::TrustRegionBased solver.
    <li> "Inexact Trust Region Based" <Teuchos::ParameterList> Sublist for the NOX::Solver::InexactTrustRegionBased solver.
    <li> "Tensor Based" <Teuchos::ParameterList> Sublist for the NOX::Solver::TensorBased solver.
    <li> "Pseudo-Transient" <Teuchos::ParameterList> Sublist for the NOX::Solver::PseudoTransient solver. (Requires Thyra support enabled)
    <li> "Tensor-Krylov Based" <Teuchos::ParameterList> Sublist for the NOX::Solver::TensorBasedTest solver. (Prerelease only)
    </ul>

    \note Solvers can be constructed using a nonmember function NOX::Solver::buildSolver instead of using this object directly.

    \author Roger Pawlowski (SNL 1416)
    */

    class Factory {

    public:

      //! Constructor.
      Factory();

      //! Destructor.
      ~Factory();

      Teuchos::RCP<NOX::Solver::Generic>
      buildSolver(const Teuchos::RCP<NOX::Abstract::Group>& grp,
          const Teuchos::RCP<NOX::StatusTest::Generic>& tests,
          const Teuchos::RCP<Teuchos::ParameterList>& params);

    };

    /*! \brief Nonmember helper function for the NOX::Solver::Factory.

    \relates NOX::Solver::Factory

    */

    Teuchos::RCP<NOX::Solver::Generic>
    buildSolver(const Teuchos::RCP<NOX::Abstract::Group>& grp,
        const Teuchos::RCP<NOX::StatusTest::Generic>& tests,
        const Teuchos::RCP<Teuchos::ParameterList>& params);

  }  // namespace Solver
}  // namespace NOX
#endif
