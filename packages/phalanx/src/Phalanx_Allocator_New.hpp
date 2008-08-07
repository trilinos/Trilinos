#ifndef PHX_ALLOCATOR_NEW_HPP
#define PHX_ALLOCATOR_NEW_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"

namespace PHX {
  

  /*! \brief Allocator that uses "new" to allocate each array separately.

      This object just uses new on the fly to allocate field data.  The memory for different fields is therefore non-contiguous.  This is a very safe allocator, but may be inefficient due to cache issues/page thrashing.
   */
  class NewAllocator {
    
  public:
    
    NewAllocator() {}
    
    ~NewAllocator() {}
    
    void reset() {}

    //! data_type_size is the size of a single element of the data type and num_elements is the number of elements of the data type that need to be allocated.
    void addRequiredChunk(std::size_t size_of_data_type, 
			  std::size_t num_elements) {}
    
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
