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
#ifndef MUELU_TWOKEYMAP_DEF_HPP
#define MUELU_TWOKEYMAP_DEF_HPP

#include <Teuchos_TabularOutputter.hpp>

#include <Xpetra_Matrix.hpp>

#include "MueLu_TwoKeyMap_decl.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {
  namespace UTILS {

    template <class Key1, class Key2, class Value>
    TwoKeyMap<Key1, Key2, Value>::TwoKeyMap() { }

    template <class Key1, class Key2, class Value>
    void TwoKeyMap<Key1, Key2, Value>::Set(const Key1 & key1, const Key2 & key2, const Value & entry) {
      SubMap & subMap = map_[key1]; // operator [] will add an entry in map_ for key 'key1' if it does not exists

      // if (subMap.count(key2) != 0)
      //  GetOStream(Warnings0, 0) << "Warning: MueLu:TwoKeyMap::Set(): Value already exist for key (" << key1 << ", " << key2 << ")" << std::endl;

      subMap[key2] = entry; // operator [] will add an entry in subMap for key 'key2' if it does not exists
    }

    template <class Key1, class Key2, class Value>
    const Value & TwoKeyMap<Key1, Key2, Value>::Get(const Key1 & key1, const Key2 & key2) const {
      // Get SubMap
      typename Map::const_iterator it = map_.find(key1);
      TEUCHOS_TEST_FOR_EXCEPTION(it == map_.end(), Exceptions::RuntimeError, "MueLu::TwoKeyMap::Get(): Key (" << key1 << ", *) does not exist.");
      const SubMap & subMap = map_.find(key1)->second;

      // Get Value
      typename SubMap::const_iterator subit = it->second.find(key2);
      TEUCHOS_TEST_FOR_EXCEPTION(subit == it->second.end(), Exceptions::RuntimeError, "MueLu::TwoKeyMap::Get(): Key (" << key1 << ", " << key2 << ") does not exist.");
      return subit->second;
    }

    template <class Key1, class Key2, class Value>
    Value & TwoKeyMap<Key1, Key2, Value>::Get(const Key1 & key1, const Key2 & key2) {
      return const_cast<Value &>(const_cast<const TwoKeyMap &>(*this).Get(key1, key2)); // Valid cast. See Effective C++, Item 3.
    }

    template <class Key1, class Key2, class Value>
    const Teuchos::map<Key2, Value> & TwoKeyMap<Key1, Key2, Value>::Get(const Key1 & key1) const {
      // Get SubMap
      typename Map::const_iterator it = map_.find(key1);
      TEUCHOS_TEST_FOR_EXCEPTION(it == map_.end(), Exceptions::RuntimeError, "MueLu::TwoKeyMap::Get(): Key (" << key1 << ", *) does not exist.");
      return it->second;
    }

    template <class Key1, class Key2, class Value>
    void TwoKeyMap<Key1, Key2, Value>::Remove(const Key1 & key1, const Key2 & key2) {
      // Get SubMap
      typename Map::iterator it = map_.find(key1);
      TEUCHOS_TEST_FOR_EXCEPTION(it == map_.end(), Exceptions::RuntimeError, "MueLu::TwoKeyMap::Remove(): Key (" << key1 << ", *) does not exist.");
      SubMap & subMap = it->second;

      // Erase SubMap. If subMap[key2] do not exist, then subMap.erase() returns 0.
      int nElementErased = subMap.erase(key2);
      TEUCHOS_TEST_FOR_EXCEPTION(nElementErased != 1, Exceptions::RuntimeError, "MueLu::TwoKeyMap::Remove(): Key (" << key1 << ", " << key2 << ") does not exist.");

      // Check if there exist other instances of 'key1' (with another key2 than 'key2')
      if (subMap.size() == 0)
        map_.erase(it); // remove SubMap from the main Map.
    }

    template <class Key1, class Key2, class Value>
    bool TwoKeyMap<Key1, Key2, Value>::IsKey(const Key1 & key1, const Key2 & key2) const {
      typename Map::const_iterator it = map_.find(key1);
      if (it == map_.end())         // return False if (key1, *) does not exist
        return false;
      const SubMap & subMap = it->second;

      return subMap.count(key2);    // True or False according if entry (key1, key2) exists.
    }

    template <class Key1, class Key2, class Value>
    bool TwoKeyMap<Key1, Key2, Value>::IsKey(const Key1 & key1) const {
      return map_.count(key1);
    }

    //TODO: GetKeyList and GetKey2List looks expensive. Also return by value a vector...
    template <class Key1, class Key2, class Value>
    std::vector<Key1> TwoKeyMap<Key1, Key2, Value>::GetKeyList() const {
      std::vector<Key1> v;
      v.reserve(map_.size());
      for(typename Map::const_iterator it = map_.begin(); it != map_.end(); ++it)
        v.push_back(it->first);

      return v;
    }

    template <class Key1, class Key2, class Value>
    std::vector<Key2> TwoKeyMap<Key1, Key2, Value>::GetKey2List(const Key1 & key1) const {
      typename Map::const_iterator it = map_.find(key1);
      TEUCHOS_TEST_FOR_EXCEPTION(it == map_.end(), Exceptions::RuntimeError, "MueLu::TwoKeyMap::GetKey2List(): Key (" << key1 << ", *) does not exist.");
      const SubMap & subMap = it->second;

      std::vector<Key2> v;
      v.reserve(subMap.size());
      for(typename SubMap::const_iterator it = subMap.begin(); it != subMap.end(); ++it)
        v.push_back(it->first);

      return v;
    }

    template <class Key1, class Key2, class Value>
    void TwoKeyMap<Key1, Key2, Value>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
      Teuchos::TabularOutputter outputter(out);
      outputter.pushFieldSpec("key1",  Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 12);
      outputter.pushFieldSpec("key2",  Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 18);
      outputter.pushFieldSpec("value", Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 18);
      outputter.outputHeader();

      std::vector<Key1> ekeys = GetKeyList();
      for(typename Map::const_iterator it1 = map_.begin(); it1 != map_.end(); ++it1) {
        const Key1   & key1   = it1->first;
        const SubMap & subMap = it1->second;

        for(typename SubMap::const_iterator it2 = subMap.begin(); it2 != subMap.end(); ++it2) {
          Key2 key2 = it2->first;

          outputter.outputField(key1);
          outputter.outputField(key2);
          //outputter.outputField(GetType(key1, key2));

          outputter.nextRow();
        }
      }
    }

  } // namespace UTILS

} // namespace MueLu

#endif // MUELU_TWOKEYMAP_DEF_HPP
