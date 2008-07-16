// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Robert Kirby (robert.c.kirby@ttu.edu) 
// ************************************************************************
// @HEADER

/** \file FIAT_AffineMapDef.hpp
    \brief  Contains implementation of affine mappings between simplices
    \author Created by R. Kirby.
*/

#include "FIAT_AffineMap.hpp"

namespace FIAT
{
  template<class Scalar>
  AffineMap<Scalar>::AffineMap( RefCountPtr<SerialDenseMatrix<int,Scalar> > baseVerts ,
				RefCountPtr<SerialDenseMatrix<int,Scalar> > targetVerts ) 
  {
    int m = baseVerts->numRows();
    int dim_x = baseVerts->numCols();
    int mprime = targetVerts->numRows();
    int dim_y = targetVerts->numCols();
    
    TEST_FOR_EXCEPTION( m != mprime ,
			std::invalid_argument ,
			">>> ERROR (FIAT::AffineMap::AffineMap): transformation will not be well-defined" );

    SerialDenseMatrix<int,Scalar> mat( dim_x*dim_y+dim_y , dim_x*dim_y+dim_y);
    SerialDenseMatrix<int,Scalar> rhs( dim_x*dim_y+dim_y , 1);
    SerialDenseMatrix<int,Scalar> sol( dim_x*dim_y+dim_y , 1); 

    for (int i=0;i<m;i++) {
      for (int j=0;j<dim_y;j++) {
  	    int row_cur = i*dim_y+j;
	    int col_start = dim_x * j;

	    for (int k=0;k<dim_x;k++) {
	      mat(row_cur,k+col_start) = (*baseVerts)(i,k);
	    }
	    rhs(row_cur,0) = (*targetVerts)(i,j);
	    mat(row_cur,dim_y*dim_x+j) = 1.0;	
      }
    }

    int info;
    SerialDenseSolver<int,Scalar> solver;
    solver.setMatrix( Teuchos::rcp( &mat , false ) );
    solver.setVectors( Teuchos::rcp( &sol , false ) , Teuchos::rcp( &rhs , false ) );

    info = solver.factor( );
    TEST_FOR_EXCEPTION( info ,
			std::runtime_error ,
			">>> ERROR (FIAT::AffineMap::AffineMap): matrix factorization failed" );
    info = solver.solve( );
    TEST_FOR_EXCEPTION( info ,
			std::runtime_error ,
			">>> ERROR (FIAT::AffineMap::AffineMap): back substitution failed" );
    
    A_ = Teuchos::rcp( new SerialDenseMatrix<int,Scalar>( Teuchos::Copy ,
							  sol.values( ) , dim_y , dim_y , dim_x ) );
    b_ = Teuchos::rcp( new SerialDenseMatrix<int,Scalar>( Teuchos::Copy ,
							  sol.values( ) + dim_y * dim_x ,
							  dim_y , dim_y , 1 ) );

  }
}
