#ifndef PHX_ALLOCATOR_CONTIGUOUS_HPP
#define PHX_ALLOCATOR_CONTIGUOUS_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_TestForException.hpp"

namespace PHX {
  
  class ContiguousAllocator {
    
  public:
    
    ContiguousAllocator() :
      allocate_called_(false),
      required_bytes_(0),
      offset_(0),
      chunk_(NULL)
    {}
    
    ~ContiguousAllocator() 
    {
      free(chunk_);
    }
    
    void addRequiredBytes(std::size_t num_bytes)
    { 
      if (allocate_called_) {
	std::string msg = "ERROR - PHX::ContiguousAllocator::addRequiredByte() - The method addRequiredByte() has been called after an allocate() has been called!";
	TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }
      required_bytes_ += num_bytes;
    }
    
    template<class DataT> Teuchos::RCP<DataT> allocate(std::size_t num_bytes)
    { 
      allocate_called_ = true;
      
      return Teuchos::rcp(new DataT);
      
      TEST_FOR_EXCEPTION(true, std::logic_error, "Roger - need to implement the allocate() method for the contiguous allocator!");
    }
    
  private:

    bool allocate_called_;

    std::size_t required_bytes_;

    std::size_t offset_;
    
    unsigned char* chunk_;
    
  };

} 
#endif
