// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_SINGLE_RHS_LINEAR_OP_WITH_SOLVE_BASE_HPP
#define THYRA_SINGLE_RHS_LINEAR_OP_WITH_SOLVE_BASE_HPP

#include "Thyra_SingleRhsLinearOpWithSolveBaseDecl.hpp"
#include "Thyra_SingleScalarLinearOpWithSolveBase.hpp"

namespace Thyra {
  
// Overridden from SingleScalarLinearOpWithSolveBase

template <class Scalar>
void SingleRhsLinearOpWithSolveBase<Scalar>::solve(
  const ETransp                         M_trans
  ,const MultiVectorBase<Scalar>        &B
  ,MultiVectorBase<Scalar>              *X
  ,const int                            numBlocks
  ,const BlockSolveCriteria<Scalar>     blockSolveCriteria[]
  ,SolveStatus<Scalar>                  blockSolveStatus[]
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( numBlocks < 0 );
  TEST_FOR_EXCEPT( numBlocks > 0 && blockSolveCriteria==NULL );
#endif
  const VectorSpaceBase<Scalar> &space_mv_rows = *B.domain();
  const Index num_mv_cols = space_mv_rows.dim();
  if(numBlocks) {
    // There is client-requested solve criteria so we need to keep track of this
    Index j = 0;
    for( int block_i = 0; block_i < numBlocks; ++block_i ) {
      const BlockSolveCriteria<Scalar> &solveCriteria = blockSolveCriteria[block_i];
      SolveStatus<Scalar> overallSolveStatus;
      overallSolveStatus.solveStatus = SOLVE_STATUS_CONVERGED; // Initialized for accumulateSolveStatus()
      for( Index block_col_i = 0; block_col_i < solveCriteria.numRhs; ++block_col_i, ++j ) {
        SolveStatus<Scalar>
          solveStatus = this->solve(M_trans,*B.col(j),&*X->col(j),&solveCriteria.solveCriteria);
        accumulateSolveStatus( solveCriteria.solveCriteria, solveStatus, &overallSolveStatus );
      }
      if(blockSolveStatus) blockSolveStatus[block_i] = overallSolveStatus;
    }
  }
  else {
    // There is not client-requested solve criteria so just solve the systems
    // with the default tolerenaces.
    for( Index j = 0; j < num_mv_cols; ++j )
      this->solve(M_trans,*B.col(j),&*X->col(j),NULL);
  }
}

} // namespace Thyra

#endif // THYRA_SINGLE_RHS_LINEAR_OP_WITH_SOLVE_BASE_HPP
