// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
template<typename DataT,int Rank,typename Layout>
const std::string PHX::Field<DataT,Rank,Layout>::m_field_tag_error_msg =
  "Error - PHX::Field - No tag has been set!";

template<typename DataT,int Rank,typename Layout>
const std::string PHX::Field<DataT,Rank,Layout>::m_field_data_error_msg =
  "Error - PHX::Field - No data has been set!  Please bind memory to Field (call getFieldData())!";
#endif

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
PHX::Field<DataT,Rank,Layout>::
Field(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& dl)
#ifdef PHX_DEBUG
  : m_data_set(false)
#endif
{
  m_tag = Teuchos::rcp(new PHX::Tag<DataT>(name,dl));
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
PHX::Field<DataT,Rank,Layout>::Field(const PHX::FieldTag& t) :
  m_tag(t.clone())
#ifdef PHX_DEBUG
  , m_data_set(false)
#endif
{ }

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
PHX::Field<DataT,Rank,Layout>::Field(const Teuchos::RCP<const PHX::FieldTag>& t) :
  m_tag(t)
#ifdef PHX_DEBUG
  , m_data_set(false)
#endif
{ }

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
PHX::Field<DataT,Rank,Layout>::Field()
#ifdef PHX_DEBUG
  : m_data_set(false)
#endif
{ }

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
template<typename CopyDataT>
PHX::Field<DataT,Rank,Layout>::Field(const PHX::Field<CopyDataT,Rank,Layout>& source) :
  m_tag(source.m_tag),
  m_field_data(source.m_field_data)
#ifdef PHX_DEBUG
  ,m_data_set(source.m_data_set)
#endif
{
  static_assert(std::is_same<typename std::decay<DataT>::type, typename std::decay<CopyDataT>::type>::value,
                "ERROR: PXH::Field copy ctor requires scalar types to be the same!");
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
PHX::Field<DataT,Rank,Layout>::~Field()
{ }

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
const PHX::FieldTag&
PHX::Field<DataT,Rank,Layout>::fieldTag() const
{
#if defined( PHX_DEBUG)
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, m_field_tag_error_msg);
#endif
  return *m_tag;
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
Teuchos::RCP<const PHX::FieldTag>
PHX::Field<DataT,Rank,Layout>::fieldTagPtr() const
{
#if defined( PHX_DEBUG)
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, m_field_tag_error_msg);
#endif
  return m_tag;
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
template<typename CopyDataT>
PHX::Field<DataT,Rank,Layout>&
PHX::Field<DataT,Rank,Layout>::operator=(const Field<CopyDataT,Rank,Layout>& source)
{
  m_tag = source.m_tag;
  m_field_data = source.m_field_data;
#ifdef PHX_DEBUG
  m_data_set = source.m_data_set;
#endif
  static_assert(std::is_same<typename std::decay<DataT>::type, typename std::decay<CopyDataT>::type>::value,
                "ERROR: PHX::Field assignment operator requires scalar types to be the same!");
  return *this;
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
template<typename... index_pack>
KOKKOS_INLINE_FUNCTION
typename PHX::FieldReturnType<typename PHX::Field<DataT,Rank,Layout>::array_type>::return_type
PHX::Field<DataT,Rank,Layout>::operator()(const index_pack&... indices) const
{
#if defined( PHX_DEBUG)
  KOKKOS_IF_ON_HOST((
  static_assert(Rank == sizeof...(indices), "PHX::Field::operator(const index_pack&... indices) : must have number of indices equal to rank!");
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);))
#endif

  return m_field_data(indices...);
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
template<typename... index_pack>
KOKKOS_INLINE_FUNCTION
typename PHX::FieldReturnType<typename PHX::Field<DataT,Rank,Layout>::array_type>::return_type
PHX::Field<DataT,Rank,Layout>::access(const index_pack&... indices) const
{
#if defined( PHX_DEBUG)
  KOKKOS_IF_ON_HOST((
  static_assert(Rank == sizeof...(indices), "PHX::Field::operator(const index_pack&... indices) : must have number of indices equal to rank!");
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);))
#endif

  return m_field_data.access(indices...);
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
KOKKOS_INLINE_FUNCTION
typename PHX::Field<DataT,Rank,Layout>::size_type
PHX::Field<DataT,Rank,Layout>::rank() const
{
  return m_field_data.rank;
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
KOKKOS_INLINE_FUNCTION
typename PHX::Field<DataT,Rank,Layout>::size_type
PHX::Field<DataT,Rank,Layout>::size() const
{
  return m_field_data.size();
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
void PHX::Field<DataT,Rank,Layout>::setFieldTag(const PHX::FieldTag& t)
{
  m_tag = t.clone();
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
void PHX::Field<DataT,Rank,Layout>::
setFieldTag(const Teuchos::RCP<const PHX::FieldTag>& t)
{
  m_tag = t;
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
void PHX::Field<DataT,Rank,Layout>::setFieldData(const std::any& a)
{
#if defined( PHX_DEBUG)
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, m_field_tag_error_msg);
  m_data_set = true;
#endif

  // any object is always the non-const data type.  To correctly cast
  // the any object to the Kokkos::View, need to pull the const off
  // the scalar type if this Field has a const scalar type.
  typedef Kokkos::View<typename array_type::non_const_data_type,Layout,PHX::Device> non_const_view;
  try {
    non_const_view tmp = std::any_cast<non_const_view>(a);
    m_field_data = tmp;
  }
  catch (std::exception& ) {
    std::cout << "\n\nError in compiletime PHX::Field::setFieldData() in std::any_cast. Tried to cast the field \""
	      << this->fieldTag().name()  << "\" with the identifier \"" << this->fieldTag().identifier()
	      << "\" to a type of \"" << Teuchos::demangleName(typeid(non_const_view).name())
	      << "\" from a std::any object containing a type of \""
	      << Teuchos::demangleName(a.type().name()) << "\"." << std::endl;
    throw;
  }
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
void PHX::Field<DataT,Rank,Layout>::releaseFieldData()
{
#if defined( PHX_DEBUG)
  m_data_set = false;
#endif
  m_field_data = array_type();
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
void PHX::Field<DataT,Rank,Layout>::print(std::ostream& os, bool printValues) const
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
    os << "Error - Field no longer supports the \"printValues\" member of the Field::print() method. Values may be on a device that does not support printing (e.g. GPU).  Please disconstinue the use of this call!" << std::endl;
}

// *********************************************************************
template<typename DataT,int Rank,typename Layout>
KOKKOS_INLINE_FUNCTION
Kokkos::DynRankView<DataT,Layout,PHX::Device>
PHX::Field<DataT,Rank,Layout>::get_view()
{
  return m_field_data;
}

// *********************************************************************
template<typename DataT,int Rank,typename Layout>
KOKKOS_INLINE_FUNCTION
const Kokkos::DynRankView<DataT,Layout,PHX::Device>
PHX::Field<DataT,Rank,Layout>::get_view() const
{
  return m_field_data;
}

// *********************************************************************
template<typename DataT,int Rank,typename Layout>
KOKKOS_INLINE_FUNCTION
typename PHX::Field<DataT,Rank,Layout>::array_type
PHX::Field<DataT,Rank,Layout>::get_static_view()
{
  return m_field_data;
}

// *********************************************************************
template<typename DataT,int Rank,typename Layout>
KOKKOS_INLINE_FUNCTION
const typename PHX::Field<DataT,Rank,Layout>::array_type
PHX::Field<DataT,Rank,Layout>::get_static_view() const
{
  return m_field_data;
}

// *********************************************************************
template<typename DataT,int Rank,typename Layout>
template<typename SrcDataT>
void
PHX::Field<DataT,Rank,Layout>::deep_copy(const PHX::Field<SrcDataT,Rank,Layout>& source)
{
  Kokkos::deep_copy(m_field_data, source.get_static_view());
}

// *************************************************************************
template<typename DataT,int Rank,typename Layout>
void
PHX::Field<DataT,Rank,Layout>::deep_copy(const DataT source)
{
  Kokkos::deep_copy(m_field_data, source);
}

// **********************************************************************
template<typename DataT,int Rank,typename Layout>
std::ostream& PHX::operator<<(std::ostream& os,
			      const PHX::Field<DataT,Rank,Layout>& f)
{
  f.print(os, false);
  return os;
}

// **********************************************************************

#endif
