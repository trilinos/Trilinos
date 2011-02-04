#include <stk_adapt/NodeRegistry.hpp>

#if 0 && NODE_REGISTRY_MAP_TYPE_TEUCHOS_HASHTABLE
namespace Teuchos
{
  //template <class T> int hashCode(const T& x);

  template <> int hashCode(const SDCell_HashTable_Key& x)
  {
    return (int)x.getHash();
  }

}
#endif

namespace stk {
  namespace adapt {

#if !NODE_REGISTRY_MAP_ACCESSORS_INLINED
    SubDimCellData& NodeRegistry::getFromMap(SubDimCell_EntityId& subDimEntity)
      {
        //ftest(subDimEntity);
        return m_cell_2_data_map[subDimEntity];
      }
    void NodeRegistry::putInMap(SubDimCell_EntityId& subDimEntity, SubDimCellData& data)
      {
        //m_cell_2_data_map.insert(std::subDimEntity, data);
        //SubDimCellData& dataInMap = m_cell_2_data_map[subDimEntity];
        //dataInMap = data;
        //ftest(subDimEntity);
        m_cell_2_data_map[subDimEntity] = data;
      }
#endif

  }
}
