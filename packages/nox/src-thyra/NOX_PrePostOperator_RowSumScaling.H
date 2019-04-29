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

#ifndef NOX_PREPOSTOPERATOR_THYRA_ROWSUMSCALING_H
#define NOX_PREPOSTOPERATOR_THYRA_ROWSUMSCALING_H

#include "NOX_Common.H"  // for NOX_Config.h
#include "NOX_Abstract_PrePostOperator.H"
#include "Thyra_VectorBase.hpp"

// Forward Declarations
namespace NOX {
  namespace Solver {
    class Generic;
  }
}

namespace NOX {

  /** \brief Updates an inverse row sum scaling vector at a user specified time.

  */
  class RowSumScaling : public NOX::Abstract::PrePostOperator {

  public:

    enum ENOX_WhenToUpdateScaling {
      UpdateInvRowSumVectorAtBeginningOfSolve,
      UpdateInvRowSumVectorAtBeginningOfIteration,
      UpdateInvRowSumVectorAtEndOfIteration,
      UpdateInvRowSumVectorAtEndOfSolve
    };

    /** \brief Ctor
        
        @param[in] <tt>inv_row_sum_vec</tt> Inverse row sum vector.
        @param[in] <tt>s</tt> When to update the scaling
        @param[in] <tt>useSimilarityTransformJacobian</tt> If true, this algorithm will use the similarity transformed row sum scaled matrix, \$D_{x}^{-1}JD_{x}\$, when computing the row sums. Additionally the inverse row sum vector will be the product of the inverse row sums and the inverse right scaling vector. If false, it uses \$JD_x\$ to compute the row sums.  
    */
    RowSumScaling(const Teuchos::RCP< ::Thyra::VectorBase<double> >& inv_row_sum_vec,
                  const ENOX_WhenToUpdateScaling s,
                  const bool useSimilarityTransformJacobian);

    void runPreSolve(const NOX::Solver::Generic& solver) override;
    void runPreIterate(const NOX::Solver::Generic& solver) override;
    void runPostIterate(const NOX::Solver::Generic& solver) override;
    void runPostSolve(const NOX::Solver::Generic& solver) override;

    Teuchos::RCP<const ::Thyra::VectorBase<double> >
    getInvRowSumScalingVector() const;

  private:

    void computeScaling(const NOX::Solver::Generic& solver);

    //! Inverse row sum vector
    Teuchos::RCP< ::Thyra::VectorBase<double> > inv_row_sum_vec_;

    ENOX_WhenToUpdateScaling when_to_update;

    bool use_similarity_transform_jacobian_;
  };
}

#endif
