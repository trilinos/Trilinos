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

#ifndef PHX_DATA_CONTAINER_BASE_HPP
#define PHX_DATA_CONTAINER_BASE_HPP

#include <typeinfo>
#include "Teuchos_RCP.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_Evaluator_Manager.hpp"

namespace PHX {

  template<typename Traits>
  class DataContainerBase {

  public:

    DataContainerBase();

    virtual ~DataContainerBase();

    virtual void allocateField(const Teuchos::RCP<PHX::FieldTag>& t,
			       typename Traits::Allocator& a) = 0;
    
    virtual const std::type_info& dataTypeInfo() const = 0; 

    virtual std::size_t getSizeOfDataType() const = 0;

    virtual void print(std::ostream& os) const = 0;
    
  };

  template<typename Traits>
  std::ostream& 
  operator<<(std::ostream& os, 
	     const PHX::DataContainerBase<Traits>& dc);

}

#include "Phalanx_DataContainer_Base_Def.hpp"

#endif 
