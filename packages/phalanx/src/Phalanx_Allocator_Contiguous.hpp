// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PHX_ALLOCATOR_CONTIGUOUS_HPP
#define PHX_ALLOCATOR_CONTIGUOUS_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"

namespace PHX {
  
  /*! \brief Class that allocates a contiguous chunk of memory for all fields.
    
      This class will allocate all fields for all data types in one contiguous chunk of memory.  It is templated on AlignmentT which is the size that all variables should be aligned on. 

  */
  template<typename AlignmentT>
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

    //! data_type_size is the size of a single element of the data type and num_elements is the number of elements of the data type that need to be allocated.
    void addRequiredChunk(std::size_t size_of_data_type, 
			  std::size_t num_elements)
    { 
      if (setup_called_) {
	std::string msg = "ERROR - PHX::ContiguousAllocator::addRequiredByte() - The method addRequiredBytes() has been called after the setup() has been called!  Please fix your logic.";
	TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }

      std::size_t alignment_size = sizeof(AlignmentT);
      std::size_t residual = size_of_data_type % alignment_size;
      std::size_t element_size = size_of_data_type + residual;

      total_bytes_ += num_elements * element_size;
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

      std::size_t size_of_data_type = sizeof(DataT);
      std::size_t alignment_size = sizeof(AlignmentT);
      std::size_t residual = size_of_data_type % alignment_size;
      std::size_t element_size = size_of_data_type + residual;

      int required_bytes = num_elements * element_size;

      TEST_FOR_EXCEPTION(offset_ + required_bytes > total_bytes_, 
			 std::logic_error, 
			 "The requested number of bytes is larger than the size of the allocated contiguous block!");

      char* raw_data = chunk_.get();
      char* offset_data = &raw_data[offset_];
      DataT* data = reinterpret_cast<DataT*>(offset_data);

      Teuchos::ArrayRCP<DataT> array = 
	Teuchos::arcpWithEmbeddedObjPostDestroy(data, 0, num_elements, 
						chunk_, false);
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
