/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef _IFPACK_SCALINGTYPE_H_
#define _IFPACK_SCALINGTYPE_H_
/*! \file Ifpack_ScalingType.h
    \brief Ifpack_ScalingType enumerable type
 */

//! Ifpack scaling type selector.
/*! Selects the type of scaling used (if any) for Ifpack preconditioners.
*/
enum Ifpack_ScalingType {None, LeftDiagonal, RightDiagonal, 
                     SymmetricDiagonal, RowSum, ColSum, 
		           RowAndColSum};

#endif /* _IFPACK_SCALINGTYPE_H_ */
