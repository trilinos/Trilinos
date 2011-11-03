#ifndef MUELU_TWOKEYMAP_DEF_HPP
#define MUELU_TWOKEYMAP_DEF_HPP

#include <MueLu_TwoKeyMap_decl.hpp>

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
      TEUCHOS_TEST_FOR_EXCEPTION(map_.count(key1) != 1, Exceptions::RuntimeError, "MueLu::TwoKeyMap::Get(): Key (" << key1 << ", *) does not exist.");
      const SubMap & subMap = map_.find(key1)->second;
        
      // Get Value
      TEUCHOS_TEST_FOR_EXCEPTION(subMap.count(key2) != 1, Exceptions::RuntimeError, "MueLu::TwoKeyMap::Get(): Key (" << key1 << ", " << key2 << ") does not exist.");
      return subMap.find(key2)->second;
    }

    template <class Key1, class Key2, class Value>
    Value & TwoKeyMap<Key1, Key2, Value>::Get(const Key1 & key1, const Key2 & key2) {
      return const_cast<Value &>(const_cast<const TwoKeyMap &>(*this).Get(key1, key2)); // Valid cast. See Effective C++, Item 3.
    }

    template <class Key1, class Key2, class Value>
    void TwoKeyMap<Key1, Key2, Value>::Remove(const Key1 & key1, const Key2 & key2) {
      // Get SubMap
      TEUCHOS_TEST_FOR_EXCEPTION(map_.count(key1) != 1, Exceptions::RuntimeError, "MueLu::TwoKeyMap::Remove(): Key (" << key1 << ", *) does not exist.");
      SubMap & subMap = map_.find(key1)->second;

      // Erase SubMap. If subMap[key2] do not exist, then subMap.erase() returns 0.
      int nElementErased = subMap.erase(key2);
      TEUCHOS_TEST_FOR_EXCEPTION(nElementErased != 1, Exceptions::RuntimeError, "MueLu::TwoKeyMap::Remove(): Key (" << key1 << ", " << key2 << ") does not exist.");

      // Check if there exist other instances of 'key1' (with another key2 than 'key2')
      if (subMap.size() == 0)
        map_.erase(key1); // remove SubMap from the main Map.
    }

    template <class Key1, class Key2, class Value>
    bool TwoKeyMap<Key1, Key2, Value>::IsKey(const Key1 & key1, const Key2 & key2) const {
      // Return False if (key1, *) does not exist
      if (map_.count(key1) != 1) return false;
      const SubMap & subMap = map_.find(key1)->second;

      // True or False according if entry (key1, key2) exists.
      return (subMap.count(key2) == 1);
    }

    //TODO: GetKeyList and GetKey2List looks expensive. Also return by value a vector...
    template <class Key1, class Key2, class Value>
    std::vector<Key1> TwoKeyMap<Key1, Key2, Value>::GetKeyList() const {
      std::vector<Key1> v;
      for(typename Map::const_iterator it = map_.begin(); it != map_.end(); ++it) {
        v.push_back(it->first);
      }
      return v;
    }

    template <class Key1, class Key2, class Value>
    std::vector<Key2> TwoKeyMap<Key1, Key2, Value>::GetKey2List(const Key1 & key1) const {
      TEUCHOS_TEST_FOR_EXCEPTION(map_.count(key1) != 1, Exceptions::RuntimeError, "MueLu::TwoKeyMap::GetKey2List(): Key (" << key1 << ", *) does not exist.");

      std::vector<Key2> v;
      const SubMap & subMap = map_.find(key1)->second;

      for(typename SubMap::const_iterator it = subMap.begin(); it != subMap.end(); ++it) {
        v.push_back(it->first);
      }

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
