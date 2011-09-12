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


#ifndef PHX_MD_FIELD_DEF_H
#define PHX_MD_FIELD_DEF_H

#include <algorithm>
#include <sstream>
#include <vector>
#include "Teuchos_TestForException.hpp"
#include "Phalanx_Print_Utilities.hpp"

//**********************************************************************
#ifdef PHX_DEBUG
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
const std::string PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::m_field_tag_error_msg = 
  "Error - PHX::MDField::fieldTag() - No tag has been set!";

template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
const std::string PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::m_field_data_error_msg = 
  "Error - PHX::MDField::operator[] - No data has been set!  Please call getFieldData(this) on all PHX::MDField objects in providers!";
#endif

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDField(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& t) :
  m_tag(name,t)
#ifdef PHX_DEBUG
  , m_tag_set(true),
  m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDField(const PHX::Tag<DataT>& v) :
  m_tag(v)
#ifdef PHX_DEBUG
  ,m_tag_set(true),
  m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDField() :
  m_tag("???", Teuchos::null)
#ifdef PHX_DEBUG
  ,m_tag_set(false),
  m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
~MDField()
{ }

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
const PHX::FieldTag& 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
fieldTag() const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_tag_error_msg);
#endif
  return m_tag;
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5, size_type index6,
	   size_type index7, size_type index8)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6,index7,index8);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5, size_type index6,
	   size_type index7)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6,index7);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5, size_type index6)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2, size_type index3)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator[](size_type index)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data[index];
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
const DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5, size_type index6,
	   size_type index7, size_type index8) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6,index7,index8);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
const DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5, size_type index6,
	   size_type index7) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6,index7);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
const DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5, size_type index6) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
const DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
const DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
const DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2, size_type index3) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
const DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1, size_type index2) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
const DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(size_type index1) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
const DataT& PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator[](size_type index) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data[index];
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size_type 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
rank() const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_data_error_msg);
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data.rank();
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size_type 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
dimension(size_type ord) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_data_error_msg);
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data.dimension(ord);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
void PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
dimensions(std::vector<size_type>& dims)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_data_error_msg);
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  
  dims.resize(m_field_data.rank());
  for ( size_type i = 0 ; i <  m_field_data.rank(); ++i ) 
	dims[i] = m_field_data.dimension(i);
  
  // Doesn't work with shards????
  //m_field_data.template dimensions<size_type>(dims);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
inline
typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size_type 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size() const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data.size();
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
void PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
setFieldTag(const PHX::Tag<DataT>& v)
{  
#ifdef PHX_DEBUG
  m_tag_set = true;
#endif
  m_tag = v;
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
void PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
setFieldData(const Teuchos::ArrayRCP<DataT>& d)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_tag_error_msg);
  m_data_set = true;
#endif

  m_array_rcp = d;

  std::vector<size_type> array_dim;
  m_tag.dataLayout().dimensions(array_dim);

  TEST_FOR_EXCEPTION(m_array_rcp.size() != m_tag.dataLayout().size(),
		     std::runtime_error, "RCP Array size is not equal to the DataLayout size!");

#ifdef PHX_USE_COMPILETIME_ARRAY
  if (array_type::Rank != (array_dim.size()) ) {
    std::ostringstream os;
    os << "MDField rank must be equal to the DataLayout rank.\n"
       << "MDField rank = " << array_type::Rank 
       << ", data layout rank = " << array_dim.size() << " = " 
       << array_dim.size() << std::endl
       << "Offending MDField:\n"
       << *this << std::endl;
    TEST_FOR_EXCEPTION(array_type::Rank != (array_dim.size()), 
		       std::logic_error, os.str());
  }
  
  typename 
    shards::Array<DataT,shards::NaturalOrder,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> 
    array(d.get(), &array_dim[0]);
  
  m_field_data = array;
#else

  size_type rank = array_dim.size();

  Teuchos::ArrayRCP<shards::ArrayDimTag*> tags = 
    Teuchos::arcp<shards::ArrayDimTag*>(rank);

  typename 
      shards::Array<DataT,shards::NaturalOrder> array(d.get(),
						      rank,
						      &array_dim[0],
						      tags.get());

  m_field_data = array;
#endif

}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
void PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
print(std::ostream& os,	bool printValues) const
{

  std::vector<const char*> dim_names;

  PHX::PrintDimension<Tag0,array_type> pd0;
  pd0.addName(dim_names);
  PHX::PrintDimension<Tag1,array_type> pd1;
  pd1.addName(dim_names);
  PHX::PrintDimension<Tag2,array_type> pd2;
  pd2.addName(dim_names);
  PHX::PrintDimension<Tag3,array_type> pd3;
  pd3.addName(dim_names);
  PHX::PrintDimension<Tag4,array_type> pd4;
  pd4.addName(dim_names);
  PHX::PrintDimension<Tag5,array_type> pd5;
  pd5.addName(dim_names);
  PHX::PrintDimension<Tag6,array_type> pd6;
  pd6.addName(dim_names);
  PHX::PrintDimension<Tag7,array_type> pd7;
  pd7.addName(dim_names);

  os << "MDField<";

  for (std::size_t i=0; i < dim_names.size(); ++i) {
    if (i > 0)
      os << ",";
    os << std::string(dim_names[i]);
  }
  os << ">(";
  for (std::size_t i=0; i < dim_names.size(); ++i) {
    if (i > 0)
      os << ",";
    os << m_field_data.dimension(i);
  }
  os << "): ";
  
  os << m_tag;

  if (printValues) {
    os << std::endl;
    for (typename array_type::size_type i = 0; i < m_field_data.size(); ++i)
      os << "value[" << i << "] = " << m_field_data[i] << std::endl;
  }

}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
std::ostream& PHX::operator<<(std::ostream& os, 
			      const PHX::MDField<DataT,Tag0,Tag1,
			      Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>& f)
{
  f.print(os, false);
  return os;
}

//**********************************************************************




//**********************************************************************
//**********************************************************************
// Runtime Version
//**********************************************************************
//**********************************************************************




//**********************************************************************
template<typename DataT>
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
MDField(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& t) :
  m_tag(name,t)
#ifdef PHX_DEBUG
  , m_tag_set(true),
  m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT>
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
MDField(const PHX::Tag<DataT>& v) :
  m_tag(v)
#ifdef PHX_DEBUG
  ,m_tag_set(true),
  m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT>
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
MDField() :
  m_tag("???", Teuchos::null)
#ifdef PHX_DEBUG
  ,m_tag_set(false),
  m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT>
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
~MDField()
{ }

//**********************************************************************
template<typename DataT>
inline
const PHX::FieldTag& 
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
fieldTag() const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_tag_error_msg);
#endif
  return m_tag;
}

//**********************************************************************
template<typename DataT>
inline
DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5, size_type index6,
	   size_type index7, size_type index8)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6,index7,index8);
}

//**********************************************************************
template<typename DataT>
inline
DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5, size_type index6,
	   size_type index7)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6,index7);
}

//**********************************************************************
template<typename DataT>
inline
DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5, size_type index6)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6);
}

//**********************************************************************
template<typename DataT>
inline
DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5);
}

//**********************************************************************
template<typename DataT>
inline
DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4);
}

//**********************************************************************
template<typename DataT>
inline
DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2, size_type index3)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3);
}

//**********************************************************************
template<typename DataT>
inline
DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2);
}

//**********************************************************************
template<typename DataT>
inline
DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1);
}

//**********************************************************************
template<typename DataT>
inline
DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator[](size_type index)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data[index];
}

//**********************************************************************
template<typename DataT>
inline
const DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5, size_type index6,
	   size_type index7, size_type index8) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6,index7,index8);
}

//**********************************************************************
template<typename DataT>
inline
const DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5, size_type index6,
	   size_type index7) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6,index7);
}

//**********************************************************************
template<typename DataT>
inline
const DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5, size_type index6) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5,index6);
}

//**********************************************************************
template<typename DataT>
inline
const DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4, size_type index5) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4,index5);
}

//**********************************************************************
template<typename DataT>
inline
const DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2, size_type index3, 
	   size_type index4) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3,index4);
}

//**********************************************************************
template<typename DataT>
inline
const DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2, size_type index3) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2,index3);
}

//**********************************************************************
template<typename DataT>
inline
const DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1, size_type index2) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1,index2);
}

//**********************************************************************
template<typename DataT>
inline
const DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(size_type index1) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1);
}

//**********************************************************************
template<typename DataT>
inline
const DataT& PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator[](size_type index) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data[index];
}

//**********************************************************************
template<typename DataT>
inline
typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::size_type 
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
rank() const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_data_error_msg);
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data.rank();
}

//**********************************************************************
template<typename DataT>
inline
typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::size_type 
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
dimension(size_type ord) const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_data_error_msg);
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data.dimension(ord);
}

//**********************************************************************
template<typename DataT>
inline
void PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
dimensions(std::vector<size_type>& dims)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_data_error_msg);
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  
  dims.resize(m_field_data.rank());
  for ( size_type i = 0 ; i <  m_field_data.rank(); ++i ) 
	dims[i] = m_field_data.dimension(i);
  
  // Doesn't work with shards????
  //m_field_data.template dimensions<size_type>(dims);
}

//**********************************************************************
template<typename DataT>
inline
typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::size_type 
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::size() const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data.size();
}

//**********************************************************************
template<typename DataT>
void PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
setFieldTag(const PHX::Tag<DataT>& v)
{  
#ifdef PHX_DEBUG
  m_tag_set = true;
#endif
  m_tag = v;
}

//**********************************************************************
template<typename DataT>
void PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
setFieldData(const Teuchos::ArrayRCP<DataT>& d)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_tag_error_msg);
  m_data_set = true;
#endif

  m_array_rcp = d;

  std::vector<size_type> array_dim;
  m_tag.dataLayout().dimensions(array_dim);

  Teuchos::ArrayRCP<shards::ArrayDimTag*> dim_tags = 
    Teuchos::arcp<shards::ArrayDimTag*>(array_dim.size());

  typename 
    shards::Array<DataT,shards::NaturalOrder,void,void,void,void,void,void,void,void> 
    array(d.get(), (int) array_dim.size(), &array_dim[0], dim_tags.get());
  
  m_field_data = array;

}

//**********************************************************************
template<typename DataT>
void PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
print(std::ostream& os,	bool printValues) const
{

  os << "MDField(";
  for (size_type i=0; i < m_field_data.rank(); ++i) {
    if (i > 0)
      os << ",";
    os << m_field_data.dimension(i);
  }
  os << "): ";
  
  os << m_tag;

  if (printValues) {
    os << std::endl;
    for (typename array_type::size_type i = 0; i < m_field_data.size(); ++i)
      os << "value[" << i << "] = " << m_field_data[i] << std::endl;
  }

}

//**********************************************************************
template<typename DataT>
std::ostream& PHX::operator<<(std::ostream& os, 
			      const PHX::MDField<DataT,void,void,
			      void,void,void,void,void,void>& f)
{
  f.print(os, false);
  return os;
}

//**********************************************************************


#endif
