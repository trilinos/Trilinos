#ifndef PHX_ALLOCATOR_DEFAULT_HPP
#define PHX_ALLOCATOR_DEFAULT_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"

namespace PHX {
  
  class DefaultAllocator {
    
  public:
    
    DefaultAllocator() {}
    
    ~DefaultAllocator() {}
    
    void addRequiredBytes(std::size_t num_bytes) {}
    
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
