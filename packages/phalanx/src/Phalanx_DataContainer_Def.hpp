// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
	      typename Traits::Allocator& a)
{
  m_data[t] = a.template allocate<DataT>(t->dataLayout().size());
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
