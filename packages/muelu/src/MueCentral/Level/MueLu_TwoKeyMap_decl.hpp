// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_TWOKEYMAP_DECL_HPP
#define MUELU_TWOKEYMAP_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_TwoKeyMap_fwd.hpp"

namespace MueLu {
  namespace UTILS {

    //! TwoKeyMap is an associative container combining two keys with a mapped value.

    //  TwoKeyMap is currently implemented as a map of map (ie: map< Key1, map< Key2, Value > >)

    template <class Key1, class Key2, class Value>
    class TwoKeyMap : public BaseClass {

    private:

      //! Sub-map container (Key2 -> Value)
      typedef Teuchos::map<Key2, Value> SubMap;

      //! Map of a map (Key1 -> SubMap)
      typedef Teuchos::map<Key1, SubMap> Map;

      //!
      Map map_;

    public:

      TwoKeyMap();

      void Set(const Key1 & key1, const Key2 & key2, const Value & entry);

      const Value & Get(const Key1 & key1, const Key2 & key2) const;

      Value & Get(const Key1 & key1, const Key2 & key2);

      void Remove(const Key1 & key1, const Key2 & key2);

      bool IsKey(const Key1 & key1, const Key2 & key2) const;

      //TODO: GetKeyList and GetKey2List looks expensive. Also return by value a vector...

      std::vector<Key1> GetKeyList() const;

      std::vector<Key2> GetKey2List(const Key1 & key1) const;

      void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

      const Teuchos::map<Key2, Value> & Get(const Key1 & key1) const;

      bool IsKey(const Key1 & key1) const;
    };

  } // namespace UTILS

} // namespace MueLu

#define MUELU_TWOKEYMAP_SHORT
#endif // MUELU_TWOKEYMAP_DECL_HPP
