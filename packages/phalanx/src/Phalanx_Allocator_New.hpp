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
      Teuchos::ArrayRCP<DataT> ptr = Teuchos::arcp<DataT>(num_elements);
      return ptr;
    }

  };

} 
#endif
