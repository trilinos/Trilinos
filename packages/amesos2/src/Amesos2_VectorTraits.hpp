// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER


#ifndef AMESOS2_VECTORTRAITS_HPP
#define AMESOS2_VECTORTRAITS_HPP

//#include "Amesos2_config.h"

#include <Tpetra_MultiVector.hpp>


#ifdef HAVE_TPETRA_INST_INT_INT
#ifdef HAVE_AMESOS2_EPETRA
#  include <Epetra_MultiVector.h>
// and perhaps some others later...
#endif
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

#ifdef HAVE_TPETRA_INST_INT_INT
#ifdef HAVE_AMESOS2_EPETRA

  template <>
  struct VectorTraits<Epetra_MultiVector> {
    typedef double scalar_t;
    typedef int local_ordinal_t;
    typedef int global_ordinal_t;
    typedef Tpetra::Map<>::node_type node_t;

    typedef Epetra_MultiVector multivector_type;
    typedef double ptr_scalar_type; // TODO Make this a pointer
  };

#endif
#endif

}

#endif  // AMESOS2_VECTORTRAITS_HPP
