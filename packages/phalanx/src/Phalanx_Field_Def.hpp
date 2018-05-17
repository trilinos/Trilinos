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


#ifndef PHX_FIELD_DEF_HPP
#define PHX_FIELD_DEF_HPP

#include <algorithm>
#include <sstream>
#include <vector>
#include "Teuchos_Assert.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Kokkos_DynRankView_Fad.hpp" // for copy/assignment specializations
#include "Phalanx_FieldTag_Tag.hpp"

// **********************************************************************
#ifdef PHX_DEBUG
template<typename DataT,int Rank>
const std::string PHX::Field<DataT,Rank>::m_field_tag_error_msg =
  "Error - PHX::Field - No tag has been set!";

template<typename DataT,int Rank>
const std::string PHX::Field<DataT,Rank>::m_field_data_error_msg =
  "Error - PHX::Field - No data has been set!  Please bind memory (call getFieldData()) to MDField!";
#endif

// **********************************************************************
template<typename DataT,int Rank>
PHX::Field<DataT,Rank>::
Field(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& dl)
#ifdef PHX_DEBUG
  : m_data_set(false)
#endif
{
  m_tag = Teuchos::rcp(new PHX::Tag<DataT>(name,dl));
}

// **********************************************************************
template<typename DataT,int Rank>
PHX::Field<DataT,Rank>::Field(const PHX::FieldTag& t) :
  m_tag(t.clone())
#ifdef PHX_DEBUG
  , m_data_set(false)
#endif
{ }

// **********************************************************************
template<typename DataT,int Rank>
PHX::Field<DataT,Rank>::Field(const Teuchos::RCP<const PHX::FieldTag>& t) :
  m_tag(t)
#ifdef PHX_DEBUG
  , m_data_set(false)
#endif
{ }

// **********************************************************************
template<typename DataT,int Rank>
PHX::Field<DataT,Rank>::Field()
#ifdef PHX_DEBUG
  : m_data_set(false)
#endif
{ }

// **********************************************************************
template<typename DataT,int Rank>
template<typename CopyDataT>
PHX::Field<DataT,Rank>::Field(const Field<CopyDataT,Rank>& source) :
  m_tag(source.m_tag),
  m_field_data(source.m_field_data)
#ifdef PHX_DEBUG
  ,m_data_set(source.m_data_set)
#endif
{
  static_assert(std::is_same<typename std::decay<DataT>::type, typename std::decay<CopyDataT>::type>::value,
                "ERROR: Compiletime MDField copy ctor requires scalar types to be the same!");
}

// **********************************************************************
template<typename DataT,int Rank>
PHX::Field<DataT,Rank>::~Field()
{ }

// **********************************************************************
template<typename DataT,int Rank>
const PHX::FieldTag&
PHX::Field<DataT,Rank>::fieldTag() const
{
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, m_field_tag_error_msg);
#endif
  return *m_tag;
}

// **********************************************************************
template<typename DataT,int Rank>
Teuchos::RCP<const PHX::FieldTag>
PHX::Field<DataT,Rank>::fieldTagPtr() const
{
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, m_field_tag_error_msg);
#endif
  return m_tag;
}

// **********************************************************************
template<typename DataT,int Rank>
template<typename CopyDataT>
PHX::Field<DataT,Rank>&
PHX::Field<DataT,Rank>::operator=(const Field<CopyDataT,Rank>& source)
{
  m_tag = source.m_tag;
  m_field_data = source.m_field_data;
#ifdef PHX_DEBUG
  m_data_set = source.m_data_set;
#endif
  static_assert(std::is_same<typename std::decay<DataT>::type, typename std::decay<CopyDataT>::type>::value,
                "ERROR: Compiletime MDField assignment operator requires scalar types to be the same!");
  return *this;
}

// **********************************************************************
template<typename DataT,int Rank>
template<typename... index_pack>
KOKKOS_INLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::Field<DataT,Rank>::array_type>::return_type
PHX::Field<DataT,Rank>::operator()(const index_pack&... indices) const
{
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  static_assert(Rank == sizeof...(indices), "PHX::MDField::operator(const index_pack&... indices) : must have number of indices equal to rank!");
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif

  return m_field_data(indices...);
}

// **********************************************************************
template<typename DataT,int Rank>
KOKKOS_INLINE_FUNCTION
typename PHX::Field<DataT,Rank>::size_type
PHX::Field<DataT,Rank>::rank() const
{
  return m_field_data.Rank;
}

// **********************************************************************
template<typename DataT,int Rank>
KOKKOS_INLINE_FUNCTION
typename PHX::Field<DataT,Rank>::size_type
PHX::Field<DataT,Rank>::size() const
{
  return m_field_data.size();
}

// **********************************************************************
template<typename DataT,int Rank>
void PHX::Field<DataT,Rank>::setFieldTag(const PHX::FieldTag& t)
{
  m_tag = t.clone();
}

// **********************************************************************
template<typename DataT,int Rank>
void PHX::Field<DataT,Rank>::
setFieldTag(const Teuchos::RCP<const PHX::FieldTag>& t)
{
  m_tag = t;
}

// **********************************************************************
template<typename DataT,int Rank>
void PHX::Field<DataT,Rank>::setFieldData(const PHX::any& a)
{
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, m_field_tag_error_msg);
  m_data_set = true;
#endif

  // any object is always the non-const data type.  To correctly cast
  // the any object to the Kokkos::View, need to pull the const off
  // the scalar type if this Field has a const scalar type.
  typedef PHX::View<typename array_type::non_const_data_type> non_const_view;
  try {
    non_const_view tmp = PHX::any_cast<non_const_view>(a);
    m_field_data = tmp;
  }
  catch (std::exception& e) {
    std::cout << "\n\nError in compiletime PHX::Field::setFieldData() in PHX::any_cast. Tried to cast the field \""
	      << this->fieldTag().name()  << "\" with the identifier \"" << this->fieldTag().identifier()
	      << "\" to a type of \"" << Teuchos::demangleName(typeid(non_const_view).name())
	      << "\" from a PHX::any object containing a type of \""
	      << Teuchos::demangleName(a.type().name()) << "\"." << std::endl;
    throw;
  }
}

// **********************************************************************
template<typename DataT,int Rank>
void PHX::Field<DataT,Rank>::print(std::ostream& os, bool printValues) const
{
  os << "Field<" << Rank << ">(";
  for (int i=0; i < Rank; ++i) {
    if (i > 0)
      os << ",";
    os << m_field_data.extent(i);
  }
  os << "): ";

  if (nonnull(m_tag))
    os << *m_tag;

  if (printValues)
    os << "Error - Field no longer supports the \"printValues\" member of the MDField::print() method. Values may be on a device that does not support printing (e.g. GPU).  Please disconstinue the use of this call!" << std::endl;
}

// *********************************************************************
template<typename DataT,int Rank>
KOKKOS_INLINE_FUNCTION
Kokkos::DynRankView<DataT,typename PHX::DevLayout<DataT>::type,PHX::Device>
PHX::Field<DataT,Rank>::get_view()
{
  return m_field_data;
}

// *********************************************************************
template<typename DataT,int Rank>
KOKKOS_INLINE_FUNCTION
const Kokkos::DynRankView<DataT,typename PHX::DevLayout<DataT>::type,PHX::Device>
PHX::Field<DataT,Rank>::get_view() const
{
  return m_field_data;
}

// *********************************************************************
template<typename DataT,int Rank>
KOKKOS_INLINE_FUNCTION
typename PHX::Field<DataT,Rank>::array_type
PHX::Field<DataT,Rank>::get_static_view()
{
  return m_field_data;
}

// *********************************************************************
template<typename DataT,int Rank>
KOKKOS_INLINE_FUNCTION
const typename PHX::Field<DataT,Rank>::array_type
PHX::Field<DataT,Rank>::get_static_view() const
{
  return m_field_data;
}

// *********************************************************************
template<typename DataT,int Rank>
template<typename SrcDataT>
void
PHX::Field<DataT,Rank>::deep_copy(const PHX::Field<SrcDataT,Rank>& source)
{
  Kokkos::deep_copy(m_field_data, source.get_static_view());
}

// *************************************************************************
template<typename DataT,int Rank>
void
PHX::Field<DataT,Rank>::deep_copy(const DataT source)
{
  Kokkos::deep_copy(m_field_data, source);
}

// **********************************************************************
template<typename DataT,int Rank>
std::ostream& PHX::operator<<(std::ostream& os,
			      const PHX::Field<DataT,Rank>& f)
{
  f.print(os, false);
  return os;
}

// **********************************************************************

#endif
