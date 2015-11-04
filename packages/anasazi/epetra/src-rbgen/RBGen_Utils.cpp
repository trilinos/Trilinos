// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "RBGen_Utils.h"

#include "Epetra_LAPACK.h"
#include "Epetra_LocalMap.h"
#include "Epetra_MultiVector.h"
#include <vector>

namespace RBGen {

  double BasisAngle( const Epetra_MultiVector& S, const Epetra_MultiVector& T )
  {
    //
    //  Computes the largest acute angle between the two orthogonal basis
    //
    if (S.NumVectors() != T.NumVectors()) {
      return acos( 0.0 );
    } else { 
      int lwork, info = 0;
      int m = S.NumVectors(), n = T.NumVectors();
      int num_vecs = ( m < n ? m : n );
      double U[ 1 ], Vt[ 1 ];
      lwork = 3*m*n;
      std::vector<double> work( lwork );
      std::vector<double> theta( num_vecs );
      Epetra_LAPACK lapack;
      Epetra_LocalMap localMap( S.NumVectors(), 0, S.Map().Comm() );
      Epetra_MultiVector Pvec( localMap, T.NumVectors() );
      info = Pvec.Multiply( 'T', 'N', 1.0, S, T, 0.0 );
      //
      // Perform SVD on S^T*T
      //
      lapack.GESVD( 'N', 'N', num_vecs, num_vecs, Pvec.Values(), Pvec.Stride(), 
		    &theta[0], U, 1, Vt, 1, &work[0], &lwork, &info );
      assert( info == 0 );
      return (acos( theta[num_vecs-1] ) );
    }
    //
    // Default return statement, should never be executed.
    //
    return acos( 0.0 );
  } 

} // end namespace RBGen
