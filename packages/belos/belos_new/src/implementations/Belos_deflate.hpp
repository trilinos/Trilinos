// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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

#ifndef BELOS_DEFLATE_HPP
#define BELOS_DEFLATE_HPP

#include "Belos_Types.hpp"
#include "TSFCoreMultiVector.hpp"

namespace Belos {

///
/** Deflate out columns from a MultiVector.
 *
 * @param  numToRemove      [in] The number of columns to deflate out.
 * @param  indexesToRemove  [in] Array (length <tt>numToRemove</tt>) of
 *                          (1-based) indexes for columns to deflate out.
 * @param  size             [in] The number of columns in the multi-vectors <tt>mv[]</tt> to be considered
 *                          for deflation.
 * @param  numMv            [in] Number of multi-vectors
 * @param  mv               [out] Array (length <tt>numMv</tt>) of pointers to multi-vectors to
 *                          have columns deflated out.
 *
 * Preconditions:<ul>
 * <li><tt>size > 0</tt>
 * <li><tt>0 < indexesToRemove[i] <= size</tt>, for <tt>i=0...numToRemove-1</tt>
 * <li><tt>indexesToRemove[i] < indexesToRemove[i+1], for <tt>i=0...numToRemove-2</tt>
 * <li><tt>mv != NULL</tt>
 * <li><tt>mv[j] != NULL</tt>, for <tt>j=0...numMv-1</tt>
 * <li><tt>mv[j]->domain()->dim() >= size</tt>, for <tt>j=0...numMv-1</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li>The input <tt>size</tt> columns in <tt>*mv</tt> will have been compressed into
 *     <tt>size-numToRemove</tt> columns by deflating out the columns specified by
 *     <tt>indexesToRemove[]</tt>.
 * </ul>
 */
template<class Scalar>
void deflate(
	const int                               numToRemove
	,const int                              indexesToRemove[]
	,const int                              size
	,const int                              numMv[]
	,TSFCore::MultiVector<Scalar>*          mvs[]
	);

///
/** Deflate out entries from list.
 *
 * @param  numToRemove      [in] The number of entries to deflate out.
 * @param  indexesToRemove  [in] Array (length <tt>numToRemove</tt>) of
 *                          (1-based) indexes for entries to deflate out.
 * @param  size             [in] The number of total entries in <tt>a[]</tt> to be considered
 *                          for deflation.
 * @param  m                [in] Number of arrays to be deflated
 * @param  a                [out] Array (length <tt>m</tt>) of pointers to array (length <tt>size</tt>)
 *                          that will have entries deflated out.
 *
 * Preconditions:<ul>
 * <li><tt>0 < indexesToRemove[i] <= size</tt>, for <tt>i=0...numToRemove-1</tt>
 * <li><tt>indexesToRemove[i] < indexesToRemove[i+1], for <tt>i=0...numToRemove-2</tt>
 * <li><tt>size > 0</tt>
 * <li><tt>a != NULL</tt>
 * <li><tt>a[j] != NULL</tt>, for <tt>j=0...m-1</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li>The input <tt>size</tt> entries in <tt>a[]</tt> will have been compressed into
 *     <tt>size-numToRemove</tt> entries by deflating out the entries specified by
 *     <tt>indexesToRemove[]</tt>.
 * </ul>
 */
template<class T>
void deflate(
	const int                               numToRemove
	,const int                              indexesToRemove[]
	,const int                              size
	,const int                              m
	,T*                                     a[]
	);

// /////////////////////////////////////////
// Implementations

template<class Scalar>
void deflate(
	const int                               numToRemove
	,const int                              indexesToRemove[]
	,const int                              size
	,const int                              numMv
	,TSFCore::MultiVector<Scalar>*          mv[]
	)
{
	int i_lower = 0;
	for( int k = 0; k < numToRemove; ++k ) {
#ifdef TEUCHOS_DEBUG
		TEST_FOR_EXCEPT( !( 0 < indexesToRemove[k] && indexesToRemove[k] <= size ) );
		TEST_FOR_EXCEPT( k < numToRemove -1 && !( 0 < indexesToRemove[k+1] && indexesToRemove[k+1] <= size ) );
		TEST_FOR_EXCEPT( k < numToRemove -1 && !( 0 < indexesToRemove[k+1] && indexesToRemove[k+1] <= size ) );
#endif
		const int
			i_start = indexesToRemove[k]+1,
			i_end = ( k < numToRemove -1 ? indexesToRemove[k+1]-1 : size );
		if( k==0 ) i_lower = i_start-1;
		for( int i = i_start; i <= i_end; ++i ) {
			for( int j = 0; j < numMv; ++j )
				assign( &*mv[j]->col(i_lower), *mv[j]->col(i) );
			++i_lower;
		}
	}
}

template<class T>
void deflate(
	const int                               numToRemove
	,const int                              indexesToRemove[]
	,const int                              size
	,const int                              m
	,T*                                     a[]
	)
{
	int i_lower = 0;
	for( int k = 0; k < numToRemove; ++k ) {
#ifdef TEUCHOS_DEBUG
		TEST_FOR_EXCEPT( !( 0 < indexesToRemove[k] && indexesToRemove[k] <= size ) );
		TEST_FOR_EXCEPT( k < numToRemove -1 && !( 0 < indexesToRemove[k+1] && indexesToRemove[k+1] <= size ) );
		TEST_FOR_EXCEPT( k < numToRemove -1 && !( 0 < indexesToRemove[k+1] && indexesToRemove[k+1] <= size ) );
#endif
		const int
			i_start = indexesToRemove[k]+1,
			i_end = ( k < numToRemove -1 ? indexesToRemove[k+1]-1 : size );
		if( k==0 ) i_lower = i_start-1;
		for( int i = i_start; i <= i_end; ++i ) {
			for( int j = 0; j < m; ++j )
				a[j][i_lower-1] = a[j][i-1];
			++i_lower;
		}
	}
}

} // namespace Belos

#endif  // BELOS_DEFLATE_HPP
