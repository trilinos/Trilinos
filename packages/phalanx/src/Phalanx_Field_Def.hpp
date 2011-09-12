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
