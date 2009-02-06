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

#ifndef PHX_DATA_CONTAINER_DEF_HPP
#define PHX_DATA_CONTAINER_DEF_HPP

#include "Teuchos_TestForException.hpp"
#include <iostream>
#include <typeinfo>
#include <sstream>
#include <iterator>
#include "boost/mpl/at.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_TypeStrings.hpp"

// ************************************************************************
template <typename DataT, typename Traits>
Teuchos::ArrayRCP<DataT> PHX::DataContainer<DataT, Traits>::
getFieldData(const PHX::FieldTag& t)
{
  typename std::map< Teuchos::RCP<const PHX::FieldTag>, 
    Teuchos::ArrayRCP<DataT>, PHX::FTComp >::iterator it;
  it = m_data.find(Teuchos::rcp(&t, false));

  if (it == m_data.end()) {
    std::string type = PHX::typeAsString<DataT>();
    std::ostringstream msg;
    msg << "The field:\n\n" << t
	<< "\n\ndoes not exist in DataContainer of type: " 
	<< type << std::endl;
    TEST_FOR_EXCEPTION(it == m_data.end(), std::logic_error, msg.str());
  }

  return it->second;
}

// ************************************************************************
template <typename DataT, typename Traits>
void PHX::DataContainer<DataT, Traits>::
allocateField(const Teuchos::RCP<PHX::FieldTag>& t, 
	      std::size_t max_num_cells,
	      typename Traits::Allocator& a)
{
  std::size_t num_elements = t->dataLayout().size() * max_num_cells;
  m_data[t] = a.template allocate<DataT>(num_elements);
}

// ************************************************************************
template <typename DataT, typename Traits>
const std::type_info& PHX::DataContainer<DataT, Traits>::
dataTypeInfo() const
{
  return typeid(DataT);
}

// ************************************************************************
template <typename DataT, typename Traits>
std::size_t PHX::DataContainer<DataT, Traits>::
getSizeOfDataType() const
{
  return sizeof(DataT);
}

// ************************************************************************
template <typename DataT, typename Traits>
void PHX::DataContainer<DataT, Traits>::
print(std::ostream& os) const
{
  std::string type = PHX::typeAsString<DataT>();
  
  os << "********************************************" << std::endl;
  os << "PHX::DataContainer Output" << std::endl;
  os << "********************************************" << std::endl;
  os << "  Data Type = " << type << std::endl;
  os << "  My FieldTags:";

  if (m_data.size() == 0)
    os << " None!" << std::endl;
  else {
    os << std::endl;
    typename m_data_t::const_iterator it = m_data.begin();
    for (; it != m_data.end(); ++it)
      os << "    " << *(it->first) << std::endl;
  }

  os << "********************************************" << std::endl;
}

// ************************************************************************
// ************************************************************************
#endif 
