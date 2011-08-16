#ifndef CTHULHU_OPERATORFACTORY_HPP
#define CTHULHU_OPERATORFACTORY_HPP

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Operator.hpp"
#include "Cthulhu_CrsOperator.hpp"
#include "Cthulhu_Map.hpp"
#include "Cthulhu_Exceptions.hpp"

namespace Cthulhu {
  
  template <class Scalar, class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps>
  class OperatorFactory {
#include "Cthulhu_UseShortNames.hpp"

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

#define CTHULHU_OPERATORFACTORY_SHORT
#endif
