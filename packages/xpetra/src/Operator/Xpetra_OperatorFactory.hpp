#ifndef XPETRA_OPERATORFACTORY_HPP
#define XPETRA_OPERATORFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Operator.hpp"
#include "Xpetra_CrsOperator.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {
  
  template <class Scalar, class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps>
  class OperatorFactory {
#include "Xpetra_UseShortNames.hpp"

  private:
    //! Private constructor. This is a static class.
    OperatorFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Operator> Build(const RCP<const Map> &rowMap, size_t maxNumEntriesPerRow, Xpetra::ProfileType pftype = Xpetra::DynamicProfile) {
      // if const block size && blocksize == 1

      return rcp( new CrsOperator(rowMap, maxNumEntriesPerRow, pftype) );

      // elseif
      
      // return vbr

      // else

      // TEST_FOR_EXCEPTION(1,Xpetra::Exceptions::BadCast,"?");
    }
    
  };

}

#define XPETRA_OPERATORFACTORY_SHORT
#endif
