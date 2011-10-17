#ifndef MUELU_TWOKEYMAP_HPP
#define MUELU_TWOKEYMAP_HPP

#include <Teuchos_TabularOutputter.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_map.hpp>

#include "MueLu_Exceptions.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_FactoryBase.hpp"

namespace MueLu {
  namespace UTILS {

    using std::string;
    
    //! TwoKeyMap is an associative container combining two keys with a mapped value.
    //  Key1  = std::string
    //  Key2  = FactoryBase*
    //  Value = ParameterEntry
    //
    //  TwoKeyMap is currently implemented as a map of map (map< Key1, map< Key2, Value > >)
    class TwoKeyMap : public BaseClass {

    private:
     
      //! Sub-map container (Key2 -> Value)
      typedef Teuchos::map<const MueLu::FactoryBase*, Teuchos::ParameterEntry> SubMap;

      //! Map of a map (Key1 -> SubMap)
      typedef Teuchos::map<const string, SubMap > Map;

      //! 
      Map map_;

    public:

      TwoKeyMap() { };

      template<class Value> 
      void Set(const string & ename, const FactoryBase* factory, const Value & evalue) {
        SetEntry(ename, factory, Teuchos::ParameterEntry(evalue));
      }

      template<class Value>
      const Value & Get(const string & ename, const FactoryBase* factory) const {
        return Teuchos::getValue<Value>(GetEntry(ename, factory));
      }

      template<class Value>
      Value & Get(const string & ename, const FactoryBase* factory) {
        return const_cast<Value &>(const_cast<const TwoKeyMap &>(*this).Get<Value>(ename, factory)); // Valid cast. See Effective C++, Item 3.
      }

      void Remove(const string & ename, const FactoryBase* factory) {
        // Get SubMap
        TEST_FOR_EXCEPTION(map_.count(ename) != 1, Exceptions::RuntimeError, "MueLu::TwoKeyMap::Get(): Key (" << ename << ", *) does not exist.");
        SubMap & subMap = map_.find(ename)->second;

        // Erase SubMap. If subMap[factory] do not exist, then subMap.erase() returns 0.
        int nElementErased = subMap.erase(factory);
        TEST_FOR_EXCEPTION(nElementErased != 1, Exceptions::RuntimeError, "MueLu::TwoKeyMap::Get(): Key (" << ename << ", " << factory << ") does not exist.");

        // Check if there exist other instances of 'ename' (with another key2 than 'factory')
        if (subMap.size() == 0)
          map_.erase(ename); // remove SubMap from the main Map.
      }

      std::string GetType(const string & ename, const FactoryBase* factory) const {
        return GetEntry(ename, factory).getAny(true).typeName();
      }

      bool IsKey(const string & ename, const FactoryBase* factory) const {
        // Return False if (ename, *) does not exist
        if (map_.count(ename) != 1) return false;
        const SubMap & subMap = map_.find(ename)->second;

        // True or False according if entry (ename, factory) exists.
        return (subMap.count(factory) == 1);
      }

      //TODO: GetKeyList and GetFactoryList looks expensive. Also return by value a vector...

      std::vector<string> GetKeyList() const {
        std::vector<string> v;
        for(Map::const_iterator it = map_.begin(); it != map_.end(); ++it) {
          v.push_back(it->first);
        }
        return v;
      }

      std::vector<const MueLu::FactoryBase*> GetFactoryList(const string & ename) const {
        TEST_FOR_EXCEPTION(map_.count(ename) != 1, Exceptions::RuntimeError, "MueLu::TwoKeyMap::Get(): Key (" << ename << ", *) does not exist.");

        std::vector<const MueLu::FactoryBase*> v;
        const SubMap & subMap = map_.find(ename)->second;

        for(SubMap::const_iterator it = subMap.begin(); it != subMap.end(); ++it) {
          v.push_back(it->first);
        }

        return v;
      }

      void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
        Teuchos::TabularOutputter outputter(out);
        outputter.pushFieldSpec("name",               Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 12);
        outputter.pushFieldSpec("gen. factory addr.", Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 18);
        outputter.pushFieldSpec("type",               Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 18);
        outputter.outputHeader();

        std::vector<std::string> ekeys = GetKeyList();
        for(Map::const_iterator it1 = map_.begin(); it1 != map_.end(); ++it1) {
          const string & ename  = it1->first;
          const SubMap & subMap = it1->second;

          for(SubMap::const_iterator it2 = subMap.begin(); it2 != subMap.end(); ++it2) {
            const MueLu::FactoryBase* factory = it2->first;

            outputter.outputField(ename);
            outputter.outputField(factory);
            outputter.outputField(GetType(ename, factory));

            outputter.nextRow();
          }
        }
      }

    private:

      void SetEntry(const string & ename, const FactoryBase* factory, const Teuchos::ParameterEntry & entry) {
        SubMap & subMap = map_[ename]; // operator [] will add an entry in map_ for key 'ename' if it does not exists
        
        // if (subMap.count(factory) != 0)
        //  GetOStream(Warnings0, 0) << "Warning: MueLu:TwoKeyMap::Set(): Value already exist for key (" << ename << ", " << factory << ")" << std::endl;
        
        subMap[factory] = entry; // operator [] will add an entry in subMap for key 'factory' if it does not exists
      }

       const Teuchos::ParameterEntry & GetEntry(const string & ename, const FactoryBase* factory) const {
        // Get SubMap
        TEST_FOR_EXCEPTION(map_.count(ename) != 1, Exceptions::RuntimeError, "MueLu::TwoKeyMap::GetEntry(): Key (" << ename << ", *) does not exist.");
        const SubMap & subMap = map_.find(ename)->second;
        
        // Get ParameterEntry
        TEST_FOR_EXCEPTION(subMap.count(factory) != 1, Exceptions::RuntimeError, "MueLu::TwoKeyMap::GetEntry(): Key (" << ename << ", " << factory << ") does not exist.");
        return subMap.find(factory)->second;
      }

      // Not used
      Teuchos::ParameterEntry & GetEntry(const string & ename, const FactoryBase* factory) {
        return const_cast<Teuchos::ParameterEntry &>(const_cast<const TwoKeyMap &>(*this).GetEntry(ename, factory)); // Valid cast. See Effective C++, Item 3.
      }

    };

  } // namespace UTILS

} // namespace MueLu

#endif /* MUELU_TWOKEYMAP_HPP */

// JG: Some notes from code review 2011-10-12:
//
// TODO: 
// - Use factory ID instead of pointer to distinghuish factories
// - Use Teuchos::Array instead of std::vector
//
// QUESTIONS:
// - Can we use an std::map<Tuple<const std::string, const MueLu::FactoryBase*>, ... > instead?
// - Teuchos::any vs. Teuchos::ParameterEntry?
// - Teuchos::ConstNonconstObjectContainer? Avoid const and non-cosnt version of same method
// - Can we use an std::map<... , Tuple<counter, factory*> >  instead?
// - Can be more generic (template type for key1, key2)
// - Current implementation allows to get efficiently the list of entries sharing Key1. Is it useful? Is it useful to switch Key1 and Key2?
//
// - if templated Key1,Key2, Value, we can remove Set/Get with template and only keep SetEntry/GetEntry (as public methods)
//   we can then deal with ParameterEntry wrapper in Needs.
