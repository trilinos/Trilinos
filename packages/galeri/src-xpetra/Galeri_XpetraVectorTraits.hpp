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
/*
  Support for vectors.
*/

// TODO: rename variables (camelCase)
#include "Galeri_ConfigDefs.hpp"

#ifndef GALERI_XPETRAVECTORTRAITS_HPP
#define GALERI_XPETRAVECTORTRAITS_HPP

#ifdef HAVE_GALERI_XPETRA
#  include "Xpetra_Map.hpp"  // needed for specialized traits
#endif

namespace Galeri {
  
  namespace Xpetra {
    
    // Default traits
    template <class Map, class Vector>
    class VectorTraits 
    {
    public:
      static RCP<Vector> Build(const RCP<const Map> &map, size_t numVectors, bool zeroOut) {
        return rcp( new Vector(map, numVectors, zeroOut) );
      }
    };

#ifdef HAVE_GALERI_XPETRA

/*
    // Specialized traits for     Map = Xpetra::Map<...>, Vector = Xpetra::Vector<...>
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node >
    class VectorTraits <Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node>, Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal, Node > >
    {
    public:
      static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal, Node> >
        Build(const RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node> > &map, bool zeroOut)
      { return Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal, Node>::Build(map,  zeroOut); };
    };
*/

    // Specialized traits for     Map = Xpetra::Map<...>, Vector = Xpetra::MultiVector<...>
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    class VectorTraits <Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node>, Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal, Node> >
    {
    public:
      static RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal, Node> >
        Build(const RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node> > &map, size_t numVectors, bool zeroOut)
      { return Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal, Node>::Build(map, numVectors, zeroOut);};
    };

#endif

  } // namespace Xpetra

} // namespace Galeri

#endif //ifndef GALERI_XPETRAVECTORTRAITS_HPP
