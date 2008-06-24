#ifndef PHX_ALLOCATOR_NEW_HPP
#define PHX_ALLOCATOR_NEW_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"

namespace PHX {
  

  /*! \brief Allocator that uses the default "new" in the Teuchos::ArrayRCP object to allocate field data.

      This object just uses new on the fly to allocate field data.  The memory for separate fields is therefore non-contiguous.  This is a very safe allocator, but may be inefficient due to cache issues/page thrashing.
   */
  class NewAllocator {
    
  public:
    
    NewAllocator() {}
    
    ~NewAllocator() {}
    
    void reset() {}

    void addRequiredBytes(std::size_t num_bytes) {}
    
    void setup() {}

    template<class DataT> 
    Teuchos::ArrayRCP<DataT> allocate(std::size_t num_elements)
    { 
      using namespace Teuchos;
      ArrayRCP<DataT> ptr = Teuchos::arcp<DataT>(num_elements);
      return ptr;
    }

  };

} 
#endif
