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

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_STRIDEDMAPFACTORY_DECL_HPP
#define XPETRA_STRIDEDMAPFACTORY_DECL_HPP

#include <Kokkos_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_StridedMap_decl.hpp"

// This factory creates Xpetra::Map. User have to specify the exact class of
// object that he want to create (ie: a Xpetra::TpetraMap or a Xpetra::EpetraMap).

namespace Xpetra {

template<class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node = KokkosClassic::DefaultNode::DefaultNodeType>
class StridedMapFactory
{

#undef XPETRA_STRIDEDMAPFACTORY_SHORT
#include "Xpetra_UseShortNamesOrdinal.hpp"

  private:

    //! Private constructor. This is a static class.
    StridedMapFactory();

  public:

    //! Map constructor with Xpetra-defined contiguous uniform distribution.


    static RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>> 
    Build(UnderlyingLib                                 lib,
          global_size_t                                 numGlobalElements,
          GlobalOrdinal                                 indexBase,
          std::vector<size_t>&                          stridingInfo,
          const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
          LocalOrdinal                                  stridedBlockId = -1,
          GlobalOrdinal                                 offset         = 0,
          LocalGlobal                                   lg = Xpetra::GloballyDistributed);


    //! Map constructor with a user-defined contiguous distribution.


    static RCP<StridedMap> 
    Build(UnderlyingLib                                 lib,
          global_size_t                                 numGlobalElements,
          size_t                                        numLocalElements,
          GlobalOrdinal                                 indexBase,
          std::vector<size_t>&                          stridingInfo,
          const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
          LocalOrdinal                                  stridedBlockId = -1,
          GlobalOrdinal                                 offset         = 0);


    static RCP<StridedMap>
    Build(const RCP<const Map>& map, std::vector<size_t>& stridingInfo, LocalOrdinal stridedBlockId = -1, GlobalOrdinal offset = 0);


    // special constructor for generating a given subblock of a strided map
    static RCP<StridedMap> 
    Build(const RCP<const StridedMap>& map, LocalOrdinal stridedBlockId);


    //! Create copy of existing map (this just creates a copy of your map, it's not a clone in the sense of Tpetra)
    static RCP<StridedMap> 
    Build(const StridedMap& map);


    //! Map constructor with a user-defined contiguous distribution. 
    //! (for experts only. There is no special check whether the generated strided maps are valid)


    static RCP<StridedMap> 
    Build(UnderlyingLib                                  lib,
          global_size_t                                  numGlobalElements,
          const Teuchos::ArrayView<const GlobalOrdinal>& elementList,
          GlobalOrdinal                                  indexBase,
          std::vector<size_t>&                           stridingInfo,
          const Teuchos::RCP<const Teuchos::Comm<int>>&  comm,
          LocalOrdinal                                   stridedBlockId = -1,      // FIXME (mfh 03 Sep 2014) This breaks if LocalOrdinal is unsigned
          GlobalOrdinal                                  /* offset */  = 0);

};     // class StridedMapFactory


}      // namespace Xpetra


#define XPETRA_STRIDEDMAPFACTORY_SHORT

#endif  // XPETRA_STRIDEDMAPFACTORY_DECL_HPP__

// TODO: removed unused methods



