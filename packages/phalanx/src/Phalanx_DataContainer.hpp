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

#ifndef PHX_DATA_CONTAINER_HPP
#define PHX_DATA_CONTAINER_HPP

#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_FieldTag_STL_Functors.hpp"
#include "Phalanx_DataContainer_Base.hpp"

namespace PHX {

  /*! \brief Container that holds all fields associated with a specific DataT.
    
      One DataContainer is instantiated for each data type in each
      evaluation type.

  */
  template <typename DataT, typename Traits>
  class DataContainer : public PHX::DataContainerBase<Traits> {
    
  public:

    DataContainer() {}
    
    ~DataContainer() {}
    
    Teuchos::ArrayRCP<DataT> getFieldData(const PHX::FieldTag& t);
    
    void allocateField(const Teuchos::RCP<PHX::FieldTag>& t,
		       std::size_t max_num_cells,
		       typename Traits::Allocator& a);

    const std::type_info& dataTypeInfo() const;

    std::size_t getSizeOfDataType() const;

    void print(std::ostream& os) const;
    
  private:
    
    typedef std::map< Teuchos::RCP<const PHX::FieldTag>, 
                      Teuchos::ArrayRCP<DataT>, 
                      FTComp > m_data_t;

    m_data_t m_data;
    
  };
  
  template <typename DataT, typename Traits>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::DataContainer<DataT, Traits>& dc)
  {
    dc.print(os);
    return os;
  }

} 

#include "Phalanx_DataContainer_Def.hpp"

#endif 
