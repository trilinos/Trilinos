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


#ifndef PHX_FIELDTAG_TAG_DEF_HPP
#define PHX_FIELDTAG_TAG_DEF_HPP

#include <sstream>
#include "Phalanx_TypeStrings.hpp"
#include "Teuchos_TestForException.hpp"

//**********************************************************************
template<typename DataT>
PHX::Tag<DataT>::Tag(const std::string& name,
		     const Teuchos::RCP<PHX::DataLayout>& dl) :
  m_name(name),
  m_data_layout(dl)
{ }

//**********************************************************************
template<typename DataT>
PHX::Tag<DataT>::~Tag()
{ }

//**********************************************************************
template<typename DataT>
Teuchos::RCP<PHX::FieldTag> PHX::Tag<DataT>::clone() const
{
  Teuchos::RCP<PHX::FieldTag> tag = 
    Teuchos::rcp(new PHX::Tag<DataT>(m_name, m_data_layout));
  return (tag);
}

//**********************************************************************
template<typename DataT>
void PHX::Tag<DataT>::operator=(const PHX::Tag<DataT>& t)
{
  m_name = t.m_name;
  m_data_layout = t.m_data_layout;
}

//**********************************************************************
template<typename DataT>
bool PHX::Tag<DataT>::operator==(const PHX::FieldTag& t) const
{
  return (  (this->name() == t.name()) &&
	    (this->dataLayout() == t.dataLayout()) &&
	    (this->dataTypeInfo() == t.dataTypeInfo()) );
}

//**********************************************************************
template<typename DataT>
const std::string& PHX::Tag<DataT>::name() const
{ return m_name; }

//**********************************************************************
template<typename DataT>
const PHX::DataLayout& PHX::Tag<DataT>::dataLayout() const
{ return *m_data_layout; }

//**********************************************************************
template<typename DataT>
const std::type_info& PHX::Tag<DataT>::dataTypeInfo() const
{ 
  DataT tmp;
  return typeid(tmp);
}

//**********************************************************************
template<typename DataT>
const std::string PHX::Tag<DataT>::identifier() const
{
  std::ostringstream ost;

  ost << this->name() 
      << this->dataTypeInfo().name() 
      << this->dataLayout().identifier();
  return ost.str(); 
}

//**********************************************************************
template<typename DataT>
void PHX::Tag<DataT>::print(std::ostream& os) const
{
//   DataT tmp;
//   os << "Tag: " << m_name << ", " << typeid(tmp).name()
//      << ", DataLayout: " << *m_data_layout;
  os << "Tag: " << m_name << ", " << PHX::typeAsString<DataT>()
     << ", DataLayout: " << *m_data_layout;

}

//**********************************************************************
template<typename DataT>
std::ostream& PHX::operator<<(std::ostream& os, const PHX::Tag<DataT>& v)
{
  v.print(os);
  return os;
}

//**********************************************************************

#endif
