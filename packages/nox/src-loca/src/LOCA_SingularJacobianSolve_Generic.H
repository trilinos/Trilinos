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
