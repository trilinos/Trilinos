/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER
*/

#ifndef TIFPACK_CONDEST_HPP
#define TIFPACK_CONDEST_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_CondestType.hpp"
class Tifpack_Preconditioner;
class Tpetra_RowMatrix;

double Tifpack_Condest(const Tifpack_Preconditioner& IFP,
		      const Tifpack_CondestType CT,
		      const int MaxIters = 1550,
		      const double Tol = 1e-9,
		      Tpetra_RowMatrix* Matrix = 0);

#endif // TIFPACK_CONDEST_HPP
