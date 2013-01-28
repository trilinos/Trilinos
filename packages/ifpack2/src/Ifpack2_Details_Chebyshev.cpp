/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

#include "Ifpack2_Details_Chebyshev_decl.hpp"

#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION

#include "Ifpack2_Details_Chebyshev_def.hpp"
#include "Ifpack2_ETIHelperMacros.h"

// Expand this macro only in the Ifpack2::Details namespace!  It
// explicitly instantiates the Ifpack2::Details::Chebyshev class for
// the given S (scalar), LO (local ordinal), and GO (global ordinal)
// types.
#define IFPACK2_DETAILS_CHEBYSHEV(S,LO,GO) \
  template<> \
  class Chebyshev<S, \
		  Tpetra::MultiVector<S,LO,GO,Kokkos::DefaultNode::DefaultNodeType>, \
		  Tpetra::CrsMatrix<S,LO,GO,Kokkos::DefaultNode::DefaultNodeType,Kokkos::DefaultKernels<S,LO,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> >;

namespace Ifpack2 {
namespace Details {

  IFPACK2_ETI_MANGLING_TYPEDEFS()

  // Our Chebyshev implementation currently only makes sense for real
  // scalar types.
  IFPACK2_INSTANTIATE_SLG_REAL(IFPACK2_DETAILS_CHEBYSHEV)

} // namespace Details
} // namespace Ifpack2

#endif // HAVE_IFPACK2_EXPLICIT_INSTANTIATION
