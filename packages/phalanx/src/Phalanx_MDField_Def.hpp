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


#ifndef PHX_MDFIELD_DEF_HPP
#define PHX_MDFIELD_DEF_HPP

#include <algorithm>
#include <sstream>
#include <vector>
#include "Teuchos_Assert.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Phalanx_config.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_Print_Utilities.hpp"

//**********************************************************************
#ifdef PHX_DEBUG
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
const std::string PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::m_field_tag_error_msg = 
  "Error - PHX::MDField - No tag has been set!";

template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
const std::string PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::m_field_data_error_msg = 
  "Error - PHX::MDField - No data has been set!  Please bind memory (call getFieldData()) to MDField!";

template<typename DataT>
const std::string PHX::MDField<DataT>::m_field_tag_error_msg = 
  "Error - PHX::MDField - No tag has been set!";

template<typename DataT>
const std::string PHX::MDField<DataT>::m_field_data_error_msg = 
  "Error - PHX::MDField - No data has been set!  Please bind memory (call getFieldData()) to MDField!";
#endif

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDField(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& dl)
#ifdef PHX_DEBUG
  : m_data_set(false)
#endif
{
  m_tag = Teuchos::rcp(new PHX::Tag<DataT>(name,dl));
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDField(const PHX::FieldTag& t) :
  m_tag(t.clone())
#ifdef PHX_DEBUG
  , m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDField(const Teuchos::RCP<const PHX::FieldTag>& t) :
  m_tag(t)
#ifdef PHX_DEBUG
  , m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDField()
#ifdef PHX_DEBUG
  : m_data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename CopyDataT,
         typename T0, typename T1, typename T2, 
         typename T3, typename T4, typename T5,
         typename T6, typename T7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDField(const MDField<CopyDataT,T0,T1,T2,T3,T4,T5,T6,T7>& source) :
  m_tag(source.m_tag),
  m_field_data(source.m_field_data)
#ifdef PHX_DEBUG
  ,m_data_set(source.m_data_set)
#endif  
{
  static_assert(ArrayRank == std::decay<decltype(source)>::type::ArrayRank,
                "ERROR: Compiletime MDField copy ctor requires rank to be the same!");
  static_assert(std::is_same<typename std::decay<DataT>::type, typename std::decay<CopyDataT>::type>::value,
                "ERROR: Compiletime MDField copy ctor requires scalar types to be the same!");
}

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
const PHX::FieldTag& 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
fieldTag() const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, m_field_tag_error_msg);
#endif
  return *m_tag;
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
Teuchos::RCP<const PHX::FieldTag> 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
fieldTagPtr() const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, m_field_tag_error_msg);
#endif
  return m_tag;
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename CopyDataT,
         typename T0, typename T1, typename T2, 
         typename T3, typename T4, typename T5,
         typename T6, typename T7>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>&
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator=(const MDField<CopyDataT,T0,T1,T2,T3,T4,T5,T6,T7>& source)
{
  m_tag = source.m_tag;
  m_field_data = source.m_field_data;
#ifdef PHX_DEBUG
  m_data_set = source.m_data_set;
#endif
  static_assert(ArrayRank == std::decay<decltype(source)>::type::ArrayRank,
                "ERROR: Compiletime MDField assignment operator requires rank to be the same!");
  static_assert(std::is_same<typename std::decay<DataT>::type, typename std::decay<CopyDataT>::type>::value,
                "ERROR: Compiletime MDField assignment operator requires scalar types to be the same!");

  return *this;
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename... index_pack>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type>::return_type
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(const index_pack&... indices) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  static_assert(ArrayRank == sizeof...(indices), "PHX::MDField::operator(const index_pack&... indices) : must have number of indices equal to rank!");
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(indices...);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size_type 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
rank() const 
{
  return m_field_data.Rank;
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename iType>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size_type 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
dimension(const iType& ord) const
{ 
  return m_field_data.extent(ord);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename iType>
void PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
dimensions(std::vector<iType>& dims)
{
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, m_field_tag_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif

  dims.resize(m_field_data.Rank);
  for ( size_type i = 0 ; i <  m_field_data.Rank; ++i ) 
    dims[i] = static_cast<iType>(m_field_data.extent(i));  // dangerous
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size_type 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size() const
{
  return m_field_data.size();
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
void PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
setFieldTag(const PHX::FieldTag& t)
{  
  m_tag = t.clone();
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
void PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
setFieldTag(const Teuchos::RCP<const PHX::FieldTag>& t)
{  
  m_tag = t;
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
void PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
setFieldData(const PHX::any& a)
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, m_field_tag_error_msg);
  m_data_set = true;
#endif

  // PHX::any object is always the non-const data type.  To correctly
  // cast the any object to the Kokkos::View, need to pull the const
  // off the scalar type if this MDField has a const scalar type.
  typedef PHX::View<typename array_type::non_const_data_type> non_const_view;
  try {
    non_const_view tmp = PHX::any_cast<non_const_view>(a);
    m_field_data = tmp;
  }
  catch (std::exception& e) {
    std::cout << "\n\nError in compiletime PHX::MDField::setFieldData() in PHX::any_cast. Tried to cast the field \"" 
	      << this->fieldTag().name()  << "\" with the identifier \"" << this->fieldTag().identifier() 
	      << "\" to a type of \"" << Teuchos::demangleName(typeid(non_const_view).name()) 
	      << "\" from a PHX::any object containing a type of \"" 
	      << Teuchos::demangleName(a.type().name()) << "\"." << std::endl;
    throw;
  }
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
    os << m_field_data.extent(i);
  }
  os << "): ";
  
  if (nonnull(m_tag))
    os << *m_tag;

  if (printValues)
    os << "Error - MDField no longer supports the \"printValues\" member of the MDField::print() method. Values may be on a device that does not support printing (e.g. GPU).  Please disconstinue the use of this call!" << std::endl;  
}

//*********************************************************************
template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7>
KOKKOS_FORCEINLINE_FUNCTION
Kokkos::DynRankView<DataT,typename PHX::DevLayout<DataT>::type,PHX::Device> 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
get_view() 
{
  return m_field_data;
}

//*********************************************************************
template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7>
KOKKOS_FORCEINLINE_FUNCTION
const Kokkos::DynRankView<DataT,typename PHX::DevLayout<DataT>::type,PHX::Device>
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
get_view() const
{
  return m_field_data;
}

//*********************************************************************
template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
get_static_view() 
{
  return m_field_data;
}

//*********************************************************************
template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7>
KOKKOS_FORCEINLINE_FUNCTION
const typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
get_static_view() const
{
  return m_field_data;
}

//*********************************************************************
template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename SrcDataT,
         typename SrcTag0,typename SrcTag1, typename SrcTag2, typename SrcTag3,
         typename SrcTag4,typename SrcTag5, typename SrcTag6, typename SrcTag7>
void
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
deep_copy(const PHX::MDField<SrcDataT,SrcTag0,SrcTag1,SrcTag2,SrcTag3,
          SrcTag4,SrcTag5,SrcTag6,SrcTag7>& source)
{
  Kokkos::deep_copy(m_field_data, source.get_static_view());
}

//*************************************************************************
template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7>
void
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
deep_copy(const DataT source) 
{
  Kokkos::deep_copy(m_field_data, source);
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

#include "Phalanx_MDField_DynRank_Def.hpp"

#endif
