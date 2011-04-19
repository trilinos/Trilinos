#ifndef CTHULHU_OPERATOR_FACTORY_HPP
#define CTHULHU_OPERATOR_FACTORY_HPP

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Operator.hpp"
#include "Cthulhu_CrsOperator.hpp"
#include "Cthulhu_Map.hpp"
#include "Cthulhu_Exceptions.hpp"

#include "Cthulhu_Debug.hpp"

namespace Cthulhu {
  
  template <class ScalarType, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<ScalarType,LocalOrdinal,Node>::SparseOps>

  class OperatorFactory {
    
    typedef Cthulhu::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
    typedef Cthulhu::Operator<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Operator;
    typedef Cthulhu::CrsOperator<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOperator;

  private:
    //! Private constructor. This is a static class.
    OperatorFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Operator> Build(const RCP<const Map> &rowMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) {
      // if const block size && blocksize == 1

      return rcp( new CrsOperator(rowMap, maxNumEntriesPerRow, pftype) );

      // elseif
      
      // return vbr

      // else

      // TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"?");
    }
    
  };

}

#define CTHULHU_OPERATOR_FACTORY_SHORT
#endif
