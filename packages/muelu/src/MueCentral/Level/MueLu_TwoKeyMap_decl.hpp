#ifndef MUELU_TWOKEYMAP_DECL_HPP
#define MUELU_TWOKEYMAP_DECL_HPP

#include <Teuchos_TabularOutputter.hpp>
#include <Teuchos_map.hpp>

#include "MueLu_Exceptions.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_FactoryBase.hpp"

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

      TwoKeyMap() ;

      void Set(const Key1 & key1, const Key2 & key2, const Value & entry) ;

       const Value & Get(const Key1 & key1, const Key2 & key2) const ;

      Value & Get(const Key1 & key1, const Key2 & key2) ;

      void Remove(const Key1 & key1, const Key2 & key2) ;

      bool IsKey(const Key1 & key1, const Key2 & key2) const ;

      //TODO: GetKeyList and GetKey2List looks expensive. Also return by value a vector...

      std::vector<Key1> GetKeyList() const ;

      std::vector<Key2> GetKey2List(const Key1 & key1) const ;

      void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const ;

    };

  } // namespace UTILS

} // namespace MueLu

#endif // MUELU_TWOKEYMAP_DECL_HPP
