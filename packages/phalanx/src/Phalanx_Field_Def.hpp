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

#ifndef PHX_FIELD_DEF_H
#define PHX_FIELD_DEF_H

#include "Teuchos_TestForException.hpp"

//**********************************************************************
#ifdef PHX_DEBUG
template<typename DataT>
const std::string PHX::Field<DataT>::m_field_tag_error_msg = 
    "Error - PHX::Field::fieldTag() - No tag has been set!";
template<typename DataT>
const std::string PHX::Field<DataT>::m_field_data_error_msg = "Error - PHX::Field::operator[] - No data has been set!  Please call getFieldData(this) on all PHX::Field objects in providers!";
#endif

//**********************************************************************
template<typename DataT>
PHX::Field<DataT>::Field(const std::string& name, 
			 const Teuchos::RCP<PHX::DataLayout>& t) :
  m_tag(name,t)
#ifdef PHX_DEBUG
  , m_tag_set(true),
  m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT>
PHX::Field<DataT>::Field(const PHX::Tag<DataT>& v) :
  m_tag(v)
#ifdef PHX_DEBUG
  ,m_tag_set(true),
  m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT>
PHX::Field<DataT>::Field() :
  m_tag("???", Teuchos::null)
#ifdef PHX_DEBUG
  ,m_tag_set(false),
  m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT>
PHX::Field<DataT>::~Field()
{ }

//**********************************************************************
template<typename DataT>
inline
const PHX::FieldTag& PHX::Field<DataT>::fieldTag() const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_tag_error_msg);
#endif
  return m_tag;
}

//**********************************************************************
template<typename DataT>
inline
DataT& PHX::Field<DataT>::operator[](int index)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data[index];
}

//**********************************************************************
template<typename DataT>
inline
typename Teuchos::ArrayRCP<DataT>::Ordinal PHX::Field<DataT>::size() const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data.size();
}

//**********************************************************************
template<typename DataT>
void PHX::Field<DataT>::setFieldTag(const PHX::Tag<DataT>& v)
{  
#ifdef PHX_DEBUG
  m_tag_set = true;
#endif
  m_tag = v;
}

//**********************************************************************
template<typename DataT>
void PHX::Field<DataT>::setFieldData(const Teuchos::ArrayRCP<DataT>& d)
{ 
#ifdef PHX_DEBUG
  m_data_set = true;
#endif
  m_field_data = d;
}

//**********************************************************************
template<typename DataT>
void PHX::Field<DataT>::print(std::ostream& os) const
{
  os << "Printing Field: \n" << m_tag << std::endl;
  typedef typename Teuchos::ArrayRCP<DataT>::Ordinal size_type;
  for (size_type i = 0; i < m_field_data.size(); ++i)
    os << "value[" << i << "] = " << m_field_data[i] << std::endl;
}

//**********************************************************************
template<typename DataT>
std::ostream& PHX::operator<<(std::ostream& os, const PHX::Field<DataT>& f)
{
  f.print(os);
  return os;
}

//**********************************************************************

#endif
