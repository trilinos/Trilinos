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

#include <cstddef>
#include <string>
#include <exception>
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
      m_setup_called = false;
      m_total_bytes = 0;
      m_offset = 0;
      m_memory = Teuchos::null;
    }

    //! data_type_size is the size of a single element of the data type and num_elements is the number of elements of the data type that need to be allocated.
    void addRequiredChunk(std::size_t size_of_data_type, 
			  std::size_t num_elements)
    { 
      if (m_setup_called) {
	std::string msg = "ERROR - PHX::ContiguousAllocator::addRequiredByte() - The method addRequiredBytes() has been called after the setup() has been called!  Please fix your logic.";
	TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }

      std::size_t field_size = size_of_data_type * num_elements;
      std::size_t alignment_size = sizeof(AlignmentT);
      std::size_t remainder = field_size % alignment_size;
      std::size_t padding = 0;
      if (remainder != 0)
	padding = alignment_size - remainder;
      std::size_t field_size_with_padding = field_size + padding;

//       std::cout << "\nSize of alignment type = " 
// 		<< alignment_size << std::endl;
//       std::cout << "Size of field no padding = " << field_size << std::endl;
//       std::cout << "Size of padding = " << padding << std::endl;
//       std::cout << "Size of field with padding = " 
// 		<< field_size_with_padding << std::endl;      

      m_total_bytes += field_size_with_padding;
    }
    
    //! Called after all byte requirements are registered.  Allocates the contiguous array.
    void setup()
    {
      if (m_total_bytes > 0) {

// 	std::cout << "Size of total_bytes = " << m_total_bytes << std::endl;
// 	std::cout << "Size of alginment type = " 
// 		  << sizeof(AlignmentT) << std::endl;

	// Make sure the alignment is correct
	TEST_FOR_EXCEPTION(m_total_bytes % sizeof(AlignmentT) != 0,
			   std::logic_error,
			   "Error - a remainder of " 
			   << m_total_bytes % sizeof(AlignmentT) 
			   << " exists for the total array size divided by the alignment type size.  The total memory is not aligned with the data type of the contiguous allocator.  This should never happen.  Please contact the Phalanx development team.");
	
	// Allocate using alignment type so that the boundaries are correct
	Teuchos::ArrayRCP<AlignmentT> aligned_memory = 
	  Teuchos::arcp<AlignmentT>(m_total_bytes / sizeof(AlignmentT));

	m_memory = Teuchos::arcp_reinterpret_cast<char>(aligned_memory);

        //m_memory = Teuchos::arcp<char>(m_total_bytes);
      }
      else {
        m_memory = Teuchos::null;
      }
      m_setup_called = true;
    }

    template<class DataT> 
    Teuchos::ArrayRCP<DataT> allocate(std::size_t num_elements)
    {       
      TEST_FOR_EXCEPTION(!m_setup_called, std::logic_error, 
			 "setup() has not been called.  The memory block has therefore not been allocated yet!  Please call setup before calling allocate().");

      std::size_t size_of_data_type = sizeof(DataT);
      std::size_t field_size = size_of_data_type * num_elements;
      std::size_t alignment_size = sizeof(AlignmentT);
      std::size_t remainder = field_size % alignment_size;
      std::size_t padding = 0;
      if (remainder != 0)
	padding = alignment_size - remainder;
      std::size_t field_size_with_padding = field_size + padding;

      TEST_FOR_EXCEPTION(m_offset + field_size_with_padding > m_total_bytes, 
			 std::logic_error, 
			 "The requested number of bytes is larger than the size of the allocated contiguous block!");

      Teuchos::ArrayRCP<DataT> array;

      if (field_size > 0) {
	
	Teuchos::ArrayRCP<char> chunk = 
	  m_memory.persistingView(m_offset, field_size);
	
	// This call sets up the arcp to call the ctor and dtor for the
	// DataT when the last rcp goes away.
	// nonpod = non-plain old data
	// DO NOT include padding - The ArrayRCP will be too big
	array = Teuchos::arcp_reinterpret_cast_nonpod<DataT>(chunk);
	
	m_offset += field_size_with_padding;
      }

      return array;
    }
    
    int getTotalBytes() const
    {
      return m_total_bytes;
    }

  private:

    //! True if setup() has been called
    bool m_setup_called;

    //! Total size of memory to allocate
    std::size_t m_total_bytes;

    //! Current offset into the allocator
    std::size_t m_offset;
    
    //! Pointer to block of all contiguous memory
    Teuchos::ArrayRCP<char> m_memory;
    
  };

} 
#endif
