// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_VECTORTRAITS_HPP
#define AMESOS2_VECTORTRAITS_HPP

//#include "Amesos2_config.h"

#include <Tpetra_MultiVector.hpp>


#ifdef HAVE_AMESOS2_EPETRA
#  include <Epetra_MultiVector.h>
// and perhaps some others later...
#endif

namespace Amesos2 {

  // The declaration
  template <class Vector>
  struct VectorTraits {};

  /*******************
   * Specializations *
   *******************/

  template < typename Scalar,
             typename LocalOrdinal,
             typename GlobalOrdinal,
             typename Node >
  struct VectorTraits<
    Tpetra::MultiVector<Scalar,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Node> > {
    typedef Scalar scalar_t;
    typedef LocalOrdinal local_ordinal_t;
    typedef GlobalOrdinal global_ordinal_t;
    typedef Node node_t;

    typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>  multivector_type;
    typedef typename multivector_type::impl_scalar_type  ptr_scalar_type; // TODO Make this a pointer
  };

  template < typename Scalar,
             typename ExecutionSpace >
  struct VectorTraits<
    Kokkos::View<Scalar**,Kokkos::LayoutLeft,ExecutionSpace> > {
    typedef Scalar scalar_t;
    typedef int local_ordinal_t;
    typedef Tpetra::Map<>::global_ordinal_type global_ordinal_t;
    typedef Tpetra::Map<>::node_type node_t;

    typedef Kokkos::View<Scalar**,Kokkos::LayoutLeft,ExecutionSpace>  multivector_type;
    typedef Scalar  ptr_scalar_type; // TODO Make this a pointer
  };


#ifdef HAVE_AMESOS2_EPETRA

  template <>
  struct VectorTraits<Epetra_MultiVector> {
    typedef double scalar_t;
    typedef int local_ordinal_t;
    typedef Tpetra::Map<>::global_ordinal_type global_ordinal_t;
    typedef Tpetra::Map<>::node_type node_t;

    typedef Epetra_MultiVector multivector_type;
    typedef double ptr_scalar_type; // TODO Make this a pointer
  };

#endif

}

#endif  // AMESOS2_VECTORTRAITS_HPP
