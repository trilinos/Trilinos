// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef XPETRA_MAPEXTRACTORFACTORY_HPP_
#define XPETRA_MAPEXTRACTORFACTORY_HPP_

#include <Kokkos_DefaultNode.hpp>

#include <Xpetra_MapExtractor.hpp>

namespace Xpetra {

  // factory class
  template <class Scalar = MapExtractor<>::scalar_type,
            class LocalOrdinal = typename MapExtractor<Scalar>::local_ordinal_type
            class GlobalOrdinal = typename MapExtractor<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node = typename MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class MapExtractorFactory {
    typedef void LocalMatOps;
#undef XPETRA_MAPEXTRACTORFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

  private:
    //! Private construtor, since this is a static class
    MapExtractorFactory() {}

  public:
    /// \brief Constructor specifying the Maps.
    ///
    /// The Maps indirectly specify the linear algebra library to use
    /// (Tpetra or Epetra).
    static Teuchos::RCP<MapExtractor>
    Build (const Teuchos::RCP<const Map>& fullmap,
           const std::vector<Teuchos::RCP<const Map> >& maps)
    {
      return rcp (new MapExtractor (fullmap, maps));
    }
  };

}

#define XPETRA_MAPEXTRACTORFACTORY_SHORT
#endif /* XPETRA_MAPEXTRACTORFACTORY_HPP_ */
