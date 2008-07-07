#ifndef PHX_ALLOCATOR_CONTIGUOUS_HPP
#define PHX_ALLOCATOR_CONTIGUOUS_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"

namespace PHX {
  
  class ContiguousAllocator {
    
  public:
    
    ContiguousAllocator()
    {
      this->reset();
    }
    
    ~ContiguousAllocator() 
    { }
    
    void reset() {
      setup_called_ = false;
      total_bytes_ = 0;
      offset_ = 0;
      chunk_ = Teuchos::null;
    }

    void addRequiredBytes(int num_bytes)
    { 
      if (setup_called_) {
	std::string msg = "ERROR - PHX::ContiguousAllocator::addRequiredByte() - The method addRequiredBytes() has been called after the setup() has been called!  Please fix your logic.";
	TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }
      total_bytes_ += num_bytes;
    }
    
    //! Called after all byte requirements are registered.  Allocates the contiguous array.
    void setup()
    {
	chunk_ = Teuchos::arcp<char>(total_bytes_);
	setup_called_ = true;
    }

    template<class DataT> 
    Teuchos::ArrayRCP<DataT> allocate(std::size_t num_elements)
    { 
      
      TEST_FOR_EXCEPTION(!setup_called_, std::logic_error, 
			 "setup() has not been called.  The memory block has therefore not been allocated yet!  Please call setup before calling allocate().");

      int required_bytes = num_elements * sizeof(DataT);

      TEST_FOR_EXCEPTION(offset_ + required_bytes > total_bytes_, 
			 std::logic_error, 
			 "The requested number of bytes is larger than the size of the allocated contiguous block!");

      char* raw_data = chunk_.get();
      char* offset_data = &raw_data[offset_];
      DataT* data = reinterpret_cast<DataT*>(offset_data);

      //! \todo Need to figure out how to set extra data in ARCP so that the bulk object is deleted when all references to it are gone.

      Teuchos::ArrayRCP<DataT> array = Teuchos::arcp(data, 0, num_elements, false);
      offset_ += required_bytes;

      // call ctor on each element
      for (std::size_t i=0; i < num_elements; ++i)
	new (&data[i]) DataT;

      return array;
    }
    
    int getTotalBytes() const
    {
      return total_bytes_;
    }

  private:

    bool setup_called_;

    long int total_bytes_;

    int offset_;
    
    Teuchos::ArrayRCP<char> chunk_;
    
  };

} 
#endif
