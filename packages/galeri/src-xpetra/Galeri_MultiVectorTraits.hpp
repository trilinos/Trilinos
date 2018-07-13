// @HEADER
//
// ***********************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef GALERI_MULTIVECTORTRAITS_HPP
#define GALERI_MULTIVECTORTRAITS_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>

#include "Galeri_ConfigDefs.h"

#ifdef HAVE_GALERI_TPETRA
#include <Tpetra_Map.hpp>
#endif

#ifdef HAVE_GALERI_XPETRA
#include <Xpetra_MapFactory.hpp>
#endif

#ifdef HAVE_GALERI_XPETRA
#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Exceptions.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MapExtractor.hpp>
#endif // HAVE_GALERI_XPETRA

namespace Galeri {

  namespace Xpetra {

    typedef size_t global_size_t;

    // TODO: Epetra_Map trait not implemented

    template <typename T>
    struct UndefinedMultiVectorTraits
    {
      static inline T notDefined() { return T::this_type_is_missing_a_specialization(); }
    };

    /* Default traits (not implemented) */
    template <class Map, class MultiVector>
    class MultiVectorTraits {
    public:
      typedef void type;

      static Teuchos::RCP<MultiVector> Build(const Teuchos::RCP<const Map>& map, size_t num) { return UndefinedMultiVectorTraits<MultiVector>::notDefined(); }
    };

#ifdef HAVE_GALERI_TPETRA
    /* Specialized traits for Map = Tpetra::Map<...> */
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    class MultiVectorTraits<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>,Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > {
    public:
      typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> type;

      static Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Build(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& map, size_t num) {
        return Teuchos::rcp(new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map, num));
      }
    };
#endif // HAVE_GALERI_TPETRA

#ifdef HAVE_GALERI_XPETRA
    /* Specialized traits for Map = Xpetra::TpetraMap<...> */
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class Map>
    class MultiVectorTraits<Map,::Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > {
    public:
      typedef ::Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> type;

      static Teuchos::RCP< ::Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Build(const Teuchos::RCP<const Map>& map, size_t num) {
        return ::Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(map, num);
      }
    };
#endif // HAVE_GALERI_XPETRA

  } // namespace Xpetra

} // namespace Galeri

#endif // GALERI_MULTIVECTORTRAITS_HPP
