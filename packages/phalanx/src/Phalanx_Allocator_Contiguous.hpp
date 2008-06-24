#ifndef PHX_ALLOCATOR_CONTIGUOUS_HPP
#define PHX_ALLOCATOR_CONTIGUOUS_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"

namespace PHX {
  
  class ContiguousAllocator {
    
  public:
    
    ContiguousAllocator() :
      allocate_called_(false),
      total_bytes_(0),
      offset_(0)
    {
      this->reset();
    }
    
    ~ContiguousAllocator() 
    { }
    
    void reset() {
      allocate_called_ = false;
      total_bytes_ = 0;
      offset_ = 0;
      chunk_ = Teuchos::null;
    }

    void addRequiredBytes(int num_bytes)
    { 
      if (allocate_called_) {
	std::string msg = "ERROR - PHX::ContiguousAllocator::addRequiredByte() - The method addRequiredByte() has been called after an allocate() has been called!";
	TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }
      total_bytes_ += num_bytes;
    }
    
    template<class DataT> 
    Teuchos::RCP<DataT> allocate(std::size_t num_elements)
    { 
      if (!allocate_called_) {
	chunk_ = Teuchos::arcp<char>(total_bytes_);
	allocate_called_ = true;
      }
      
      int required_bytes = num_elements * sizeof(DataT);

      TEST_FOR_EXCEPTION(offset_ + required_bytes >= total_bytes_, 
			 std::logic_error, 
			 "The requested number of bytes is larger than the size allocated!");

      Teuchos::ArrayRCP<DataT> array = 
	Teuchos::arcp<DataT>(reinterpret_cast<DataT*>(chunk_.get()[offset_]), 0, num_elements);
      offset_ += required_bytes;
      return array;
    }
    
    int getTotalBytes() const
    {
      return total_bytes_;
    }

  private:

    bool allocate_called_;

    int total_bytes_;

    int offset_;
    
    Teuchos::ArrayRCP<char> chunk_;
    
  };

} 
#endif
