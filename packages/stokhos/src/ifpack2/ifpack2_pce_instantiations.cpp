// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_Ifpack2.hpp"

typedef Stokhos::StandardStorage<int,double> Storage;
typedef Sacado::PCE::OrthogPoly<double,Storage> pce_type;
typedef pce_type Scalar;
typedef int LocalOrdinal;
typedef int GlobalOrdinal;

#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
#include "Ifpack2_RILUK_def.hpp"
#include "Ifpack2_ILUT_def.hpp"
#include "Ifpack2_Chebyshev_def.hpp"
#include "Ifpack2_Diagonal_def.hpp"
#include "Ifpack2_Relaxation_def.hpp"
#include "Ifpack2_Krylov_def.hpp"
#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION
#include "Tpetra_CrsGraph_def.hpp"
#include "Tpetra_CrsMatrix_def.hpp"
#endif

#include "Ifpack2_ExplicitInstantiationHelpers.hpp"

// Needed to copy these from Ifpack2_Krylov.cpp since they are defined in the 
// cpp file
#define IFPACK2_INST_SPARSE_RELAX(S,LO,GO) \
  template class Krylov<Tpetra::CrsMatrix<S,LO,GO,Kokkos::DefaultNode::DefaultNodeType,Kokkos::DefaultKernels<S,LO,Kokkos::DefaultNode::DefaultNodeType>::SparseOps>, \
			Ifpack2::Relaxation<Tpetra::CrsMatrix<S,LO,LO,Kokkos::DefaultNode::DefaultNodeType,Kokkos::DefaultKernels<S,LO,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> > >;

#define IFPACK2_INST_SPARSE_ILUT(S,LO,GO) \
  template class Krylov<Tpetra::CrsMatrix<S,LO,GO,Kokkos::DefaultNode::DefaultNodeType,Kokkos::DefaultKernels<S,LO,Kokkos::DefaultNode::DefaultNodeType>::SparseOps>, \
			Ifpack2::ILUT<Tpetra::CrsMatrix<S,LO,LO,Kokkos::DefaultNode::DefaultNodeType,Kokkos::DefaultKernels<S,LO,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> > >;

#define IFPACK2_INST_SPARSE_CHEBY(S,LO,GO) \
  template class Krylov<Tpetra::CrsMatrix<S,LO,GO,Kokkos::DefaultNode::DefaultNodeType,Kokkos::DefaultKernels<S,LO,Kokkos::DefaultNode::DefaultNodeType>::SparseOps>, \
			Ifpack2::Chebyshev<Tpetra::CrsMatrix<S,LO,LO,Kokkos::DefaultNode::DefaultNodeType,Kokkos::DefaultKernels<S,LO,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> > >;

namespace Ifpack2 {
IFPACK2_INST(RILUK,Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST(ILUT,Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST(Chebyshev,Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST(Diagonal,Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST(Relaxation,Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST_SPARSE_RELAX(Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST_SPARSE_ILUT(Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST_SPARSE_CHEBY(Scalar,LocalOrdinal,GlobalOrdinal)
}

#endif
