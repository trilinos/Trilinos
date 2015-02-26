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
#include "Phalanx_Print_Utilities.hpp"

#ifdef Phalanx_ENABLE_Intrepid
#include "Intrepid_config.h" // for HAVE_INTREPID_KOKKOSCORE define
#ifdef HAVE_INTREPID_KOKKOSCORE
#include "KokkosRank.hpp"
#endif
#endif

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

template<typename DataT>
const std::string PHX::MDField<DataT>::m_field_tag_error_msg = 
  "Error - PHX::MDField::fieldTag() - No tag has been set!";

template<typename DataT>
const std::string PHX::MDField<DataT>::m_field_data_error_msg = 
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
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_tag_error_msg);
#endif
  return m_tag;
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename iType0, typename iType1, typename iType2, typename iType3,
	 typename iType4, typename iType5, typename iType6, typename iType7>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type>::return_type
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(iType0 index0, iType1 index1, iType2 index2, 
	   iType3 index3, iType4 index4, iType5 index5,
	   iType6 index6, iType7 index7) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index0,index1,index2,index3,index4,index5,index6,index7);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename iType0, typename iType1, typename iType2, typename iType3,
	 typename iType4, typename iType5, typename iType6>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type>::return_type
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(iType0 index0, iType1 index1, iType2 index2, 
	   iType3 index3, iType4 index4, iType5 index5,
	   iType6 index6)const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index0,index1,index2,index3,index4,index5,index6);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename iType0, typename iType1, typename iType2, typename iType3,
	 typename iType4, typename iType5>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type>::return_type
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(iType0 index0, iType1 index1, iType2 index2, 
	   iType3 index3, iType4 index4, iType5 index5)const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index0,index1,index2,index3,index4,index5);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename iType0, typename iType1, typename iType2, typename iType3,
	 typename iType4>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type>::return_type
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(iType0 index0, iType1 index1, iType2 index2, 
	   iType3 index3, iType4 index4)const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index0,index1,index2,index3,index4);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename iType0, typename iType1, typename iType2, typename iType3>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type>::return_type
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(iType0 index0, iType1 index1, iType2 index2, 
	   iType3 index3)const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index0,index1,index2,index3);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename iType0, typename iType1, typename iType2>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type>::return_type
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(iType0 index0, iType1 index1, iType2 index2) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index0,index1,index2);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename iType0, typename iType1>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type>::return_type
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(iType0 index0, iType1 index1) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index0,index1);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename iType0>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type>::return_type
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator()(iType0 index1) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data(index1);
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
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif

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
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_data.dimension(ord);
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename iType>
inline
void PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
dimensions(std::vector<iType>& dims)
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  
  dims.resize(m_field_data.Rank);
  for ( size_type i = 0 ; i <  m_field_data.Rank; ++i ) 
    dims[i] = static_cast<iType>(m_field_data.dimension(i));  // dangerous
}

//**********************************************************************
template<typename DataT,
	 typename Tag0,typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4,typename Tag5, typename Tag6, typename Tag7>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size_type 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size() const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
   return m_field_data.size();
  
  // loop over the dimensions and compute the size (this gets around an issue with
  // Kokkos-Fad 10/21/2014 where size() was working incorrectly)
/*  typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::size_type sz = 1;
  for ( size_type i = 0 ; i <  m_field_data.Rank; ++i ) 
    sz *= m_field_data.dimension(i);
  return sz;*/
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
setFieldData(const boost::any& a)
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_tag_error_msg);
  m_data_set = true;
#endif

  // Boost any object is always the non-const data type.  To correctly
  // cast the boost any object to the Kokkos::View, need to pull the
  // const off the scalar type if this MDField has a const scalar type.
  typedef Kokkos::View<typename array_type::non_const_array_intrinsic_type,PHX::Device> non_const_view;
  try {
    non_const_view tmp = boost::any_cast<non_const_view>(a);
    m_field_data = tmp;
  }
  catch (std::exception& e) {
    std::cout << "\n\nError in compiletime PHX::MDField::setFieldData() in boost::any_cast. Tried to cast the field \"" 
	      << this->fieldTag().name()  << "\" with the identifier \"" << this->fieldTag().identifier() 
	      << "\" to a type of \"" << Teuchos::demangleName(typeid(non_const_view).name()) 
	      << "\" from a boost::any object containing a type of \"" 
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
    os << m_field_data.dimension(i);
  }
  os << "): ";
  
  os << m_tag;

  if (printValues)
    os << "Error - MDField no longer supports the \"printValues\" member of the MDField::print() method. Values may be on a device that does not support printing (e.g. GPU).  Please disconstinue the use of this call!" << std::endl;  
}
//*********************************************************************
template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7>
typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type 
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
get_kokkos_view() 
{
  return m_field_data;
}
//*********************************************************************
template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7>
const typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
get_kokkos_view() const
{
  return m_field_data;
}
//*********************************************************************
template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename MDFieldType>
void
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
deep_copy(const MDFieldType& source)      
{
  Kokkos::deep_copy(m_field_data, source);
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
//*********************************************************************
/*template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7>
template<typename MDFieldType>
void
PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
V_Multiply(const MDFieldType& source)
{
  Kokkos::deep_copy(m_field_data
}
*/
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
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_tag_error_msg);
#endif
  return m_tag;
}

//**********************************************************************
template<typename DataT>
template<typename iType0, typename iType1, typename iType2, typename iType3,
	 typename iType4, typename iType5, typename iType6, typename iType7>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::array_type>::return_type
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(iType0 index0, iType1 index1, iType2 index2, 
	   iType3 index3, iType4 index4, iType5 index5,
	   iType6 index6, iType7 index7) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.dataLayout().rank() != 8,std::logic_error,"Error: The MDField \"" << m_tag.name() << "\" is of rank " << m_tag.dataLayout().rank() << " but was accessed with the rank 8 operator().");
#endif

  TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"Runtime PHX::MDField does not support 8 ranks!  Limited to 7.");
  return m_field_data7(index0,index1,index2,index3,index4,index5,index6,index7);
}

//**********************************************************************
template<typename DataT>
template<typename iType0, typename iType1, typename iType2, typename iType3,
	 typename iType4, typename iType5, typename iType6>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::array_type>::return_type
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(iType0 index0, iType1 index1, iType2 index2, 
	   iType3 index3, iType4 index4, iType5 index5,
	   iType6 index6) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.dataLayout().rank() != 7,std::logic_error,"Error: The MDField \"" << m_tag.name() << "\" is of rank " << m_tag.dataLayout().rank() << " but was accessed with the rank 7 operator().");
#endif
  return m_field_data7(index0,index1,index2,index3,index4,index5,index6);
}

//**********************************************************************
template<typename DataT>
template<typename iType0, typename iType1, typename iType2, typename iType3,
	 typename iType4, typename iType5>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::array_type>::return_type
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(iType0 index0, iType1 index1, iType2 index2, 
	   iType3 index3, iType4 index4, iType5 index5) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.dataLayout().rank() != 6,std::logic_error,"Error: The MDField \"" << m_tag.name() << "\" is of rank " << m_tag.dataLayout().rank() << " but was accessed with the rank 6 operator().");
#endif
  return m_field_data6(index0,index1,index2,index3,index4,index5);
}

//**********************************************************************
template<typename DataT>
template<typename iType0, typename iType1, typename iType2, typename iType3,
	 typename iType4>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::array_type>::return_type
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(iType0 index0, iType1 index1, iType2 index2, 
	   iType3 index3, iType4 index4) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.dataLayout().rank() != 5,std::logic_error,"Error: The MDField \"" << m_tag.name() << "\" is of rank " << m_tag.dataLayout().rank() << " but was accessed with the rank 5 operator().");
#endif
  return m_field_data5(index0,index1,index2,index3,index4);
}

//**********************************************************************
template<typename DataT>
template<typename iType0, typename iType1, typename iType2, typename iType3>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::array_type>::return_type
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(iType0 index0, iType1 index1, iType2 index2, 
	   iType3 index3)const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.dataLayout().rank() != 4,std::logic_error,"Error: The MDField \"" << m_tag.name() << "\" is of rank " << m_tag.dataLayout().rank() << " but was accessed with the rank 4 operator().");
#endif
  return m_field_data4(index0,index1,index2,index3);
}

//**********************************************************************
template<typename DataT>
template<typename iType0, typename iType1, typename iType2>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::array_type>::return_type
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(iType0 index0, iType1 index1, iType2 index2) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ ) 
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.dataLayout().rank() != 3,std::logic_error,"Error: The MDField \"" << m_tag.name() << "\" is of rank " << m_tag.dataLayout().rank() << " but was accessed with the rank 3 operator().");
#endif
  return m_field_data3(index0,index1,index2);
}

//**********************************************************************
template<typename DataT>
template<typename iType0, typename iType1>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::array_type>::return_type
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(iType0 index0, iType1 index1) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.dataLayout().rank() != 2,std::logic_error,"Error: The MDField \"" << m_tag.name() << "\" is of rank " << m_tag.dataLayout().rank() << " but was accessed with the rank 2 operator().");
#endif
  return m_field_data2(index0,index1);
}

//**********************************************************************
template<typename DataT>
template<typename iType0>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::array_type>::return_type
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator()(iType0 index0) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(m_tag.dataLayout().rank() != 1,std::logic_error,"Error: The MDField \"" << m_tag.name() << "\" is of rank " << m_tag.dataLayout().rank() << " but was accessed with the rank 1 operator().");
#endif
  return m_field_data1(index0);
}

//**********************************************************************
template<typename DataT>
template<typename iType0>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::array_type>::return_type
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
operator[](iType0 index0) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_field_oned_view(index0);
}

//**********************************************************************
template<typename DataT>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::size_type 
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
rank() const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_dimension_rank_size(7);
}

//**********************************************************************
template<typename DataT>
template<typename iType>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::size_type 
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
dimension(const iType& ord) const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_dimension_rank_size(ord);
}

//**********************************************************************
template<typename DataT>
template<typename iType>
inline
void PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
dimensions(std::vector<iType>& dims)
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_data_error_msg);
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  
  dims.resize(m_tag.dataLayout().rank());
  for ( size_type i = 0 ; i <  m_tag.dataLayout().rank(); ++i ) 
    dims[i] = static_cast<iType>(m_tag.dataLayout().dimension(i)); // dangerous
}

//**********************************************************************
template<typename DataT>
KOKKOS_FORCEINLINE_FUNCTION
typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::size_type 
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::size() const
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, m_field_data_error_msg);
#endif
  return m_dimension_rank_size(8);
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
template <typename T, typename L, typename D, typename M, typename S>
unsigned PHX::getSacadoSize(const Kokkos::View<T,L,D,M,S>& view) {
  return 1;
}

//**********************************************************************
template <typename T, typename L, typename D, typename M>
unsigned PHX::getSacadoSize(const Kokkos::View<T,L,D,M,Kokkos::Impl::ViewSpecializeSacadoFad>& view) {
  return view.storage_size();
}

//**********************************************************************
template<typename DataT>
void PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
setFieldData(const boost::any& a)
{ 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
  TEUCHOS_TEST_FOR_EXCEPTION(!m_tag_set, std::logic_error, m_field_tag_error_msg);
  m_data_set = true;
#endif

  m_dimension_rank_size = Kokkos::View<PHX::index_size_type*, PHX::Device>("m_dimension_rank_size",9);

  m_dimension_rank_size(0) = 1;  // dim0
  m_dimension_rank_size(1) = 1;  // dim1
  m_dimension_rank_size(2) = 1;  // dim2
  m_dimension_rank_size(3) = 1;  // dim3
  m_dimension_rank_size(4) = 1;  // dim4
  m_dimension_rank_size(5) = 1;  // dim5
  m_dimension_rank_size(6) = 1;  // dim6
  m_dimension_rank_size(7) = 0;  // rank
  m_dimension_rank_size(8) = 0;  // size

  typedef Kokkos::View<typename array_type1::non_const_array_intrinsic_type,PHX::Device> non_const_view1;
  typedef Kokkos::View<typename array_type2::non_const_array_intrinsic_type,PHX::Device> non_const_view2;
  typedef Kokkos::View<typename array_type3::non_const_array_intrinsic_type,PHX::Device> non_const_view3;
  typedef Kokkos::View<typename array_type4::non_const_array_intrinsic_type,PHX::Device> non_const_view4;
  typedef Kokkos::View<typename array_type5::non_const_array_intrinsic_type,PHX::Device> non_const_view5;
  typedef Kokkos::View<typename array_type6::non_const_array_intrinsic_type,PHX::Device> non_const_view6;
  typedef Kokkos::View<typename array_type7::non_const_array_intrinsic_type,PHX::Device> non_const_view7;
  try {

    if (m_tag.dataLayout().rank() == 1) {
      m_field_data1 = boost::any_cast<non_const_view1>(a);
      m_field_oned_view = array_oned_type(m_field_data1.ptr_on_device(),m_field_data1.size(),PHX::getSacadoSize(m_field_data1));
      m_dimension_rank_size(0) = m_field_data1.dimension_0();
      m_dimension_rank_size(7) = 1;
      m_dimension_rank_size(8) = m_field_data1.size();
    }
    else if (m_tag.dataLayout().rank() == 2) {
      m_field_data2 = boost::any_cast<non_const_view2>(a);
      m_field_oned_view = array_oned_type(m_field_data2.ptr_on_device(),m_field_data2.size(),PHX::getSacadoSize(m_field_data2));
      m_dimension_rank_size(0) = m_field_data2.dimension_0();
      m_dimension_rank_size(1) = m_field_data2.dimension_1();
      m_dimension_rank_size(7) = 2;
      m_dimension_rank_size(8) = m_field_data2.size();
    }
    else if (m_tag.dataLayout().rank() == 3) {
      m_field_data3 = boost::any_cast<non_const_view3>(a);
      m_field_oned_view = array_oned_type(m_field_data3.ptr_on_device(),m_field_data3.size(),PHX::getSacadoSize(m_field_data3));
      m_dimension_rank_size(0) = m_field_data3.dimension_0();
      m_dimension_rank_size(1) = m_field_data3.dimension_1();
      m_dimension_rank_size(2) = m_field_data3.dimension_2();
      m_dimension_rank_size(7) = 3;
      m_dimension_rank_size(8) = m_field_data3.size();
    }
    else if (m_tag.dataLayout().rank() == 4) {
      m_field_data4 = boost::any_cast<non_const_view4>(a);
      m_field_oned_view = array_oned_type(m_field_data4.ptr_on_device(),m_field_data4.size(),PHX::getSacadoSize(m_field_data4));
      m_dimension_rank_size(0) = m_field_data4.dimension_0();
      m_dimension_rank_size(1) = m_field_data4.dimension_1();
      m_dimension_rank_size(2) = m_field_data4.dimension_2();
      m_dimension_rank_size(3) = m_field_data4.dimension_3();
      m_dimension_rank_size(7) = 4;
      m_dimension_rank_size(8) = m_field_data4.size();
    }
    else if (m_tag.dataLayout().rank() == 5) {
      m_field_data5 = boost::any_cast<non_const_view5>(a);
      m_field_oned_view = array_oned_type(m_field_data5.ptr_on_device(),m_field_data5.size(),PHX::getSacadoSize(m_field_data5));
      m_dimension_rank_size(0) = m_field_data5.dimension_0();
      m_dimension_rank_size(1) = m_field_data5.dimension_1();
      m_dimension_rank_size(2) = m_field_data5.dimension_2();
      m_dimension_rank_size(3) = m_field_data5.dimension_3();
      m_dimension_rank_size(4) = m_field_data5.dimension_4();
      m_dimension_rank_size(7) = 5;
      m_dimension_rank_size(8) = m_field_data5.size();
    }
    else if (m_tag.dataLayout().rank() == 6) {
      m_field_data6 = boost::any_cast<non_const_view6>(a);
      m_field_oned_view = array_oned_type(m_field_data6.ptr_on_device(),m_field_data6.size(),PHX::getSacadoSize(m_field_data6));
      m_dimension_rank_size(0) = m_field_data6.dimension_0();
      m_dimension_rank_size(1) = m_field_data6.dimension_1();
      m_dimension_rank_size(2) = m_field_data6.dimension_2();
      m_dimension_rank_size(3) = m_field_data6.dimension_3();
      m_dimension_rank_size(4) = m_field_data6.dimension_4();
      m_dimension_rank_size(5) = m_field_data6.dimension_5();
      m_dimension_rank_size(7) = 6;
      m_dimension_rank_size(8) = m_field_data6.size();
    }
    else if (m_tag.dataLayout().rank() == 7) {
      m_field_data7 = boost::any_cast<non_const_view7>(a);
      m_field_oned_view = array_oned_type(m_field_data7.ptr_on_device(),m_field_data7.size(),PHX::getSacadoSize(m_field_data7));
      m_dimension_rank_size(0) = m_field_data7.dimension_0();
      m_dimension_rank_size(1) = m_field_data7.dimension_1();
      m_dimension_rank_size(2) = m_field_data7.dimension_2();
      m_dimension_rank_size(3) = m_field_data7.dimension_3();
      m_dimension_rank_size(4) = m_field_data7.dimension_4();
      m_dimension_rank_size(5) = m_field_data7.dimension_5(); 
      m_dimension_rank_size(6) = m_field_data7.dimension_6(); 
      m_dimension_rank_size(7) = 7;
      m_dimension_rank_size(8) = m_field_data7.size();   
    }

    // temporary size calculation to avoid bug in Kokkos 10/22/2014
    typename PHX::MDField<DataT,void,void,void,void,void,void,void,void>::size_type sz = 1;
    for ( size_type i = 0 ; i <  m_dimension_rank_size(7); ++i ) 
      sz *= m_dimension_rank_size(i);

    m_dimension_rank_size(8) = sz; 

  }
  catch (std::exception& e) {
    
    std::string type_cast_name = "???";
    if (m_tag.dataLayout().rank() == 1)
      type_cast_name = Teuchos::demangleName(typeid(non_const_view1).name());
    else if (m_tag.dataLayout().rank() == 2)
      type_cast_name = Teuchos::demangleName(typeid(non_const_view2).name());
    else if (m_tag.dataLayout().rank() == 3)
      type_cast_name = Teuchos::demangleName(typeid(non_const_view3).name());
    else if (m_tag.dataLayout().rank() == 4)
      type_cast_name = Teuchos::demangleName(typeid(non_const_view4).name());
    else if (m_tag.dataLayout().rank() == 5)
      type_cast_name = Teuchos::demangleName(typeid(non_const_view5).name());
    else if (m_tag.dataLayout().rank() == 6)
      type_cast_name = Teuchos::demangleName(typeid(non_const_view6).name());
    else if (m_tag.dataLayout().rank() == 7)
      type_cast_name = Teuchos::demangleName(typeid(non_const_view7).name());

    std::cout << "\n\nError in runtime PHX::MDField::setFieldData() in boost::any_cast. Tried to cast the field \"" 
	      << this->fieldTag().name()  << "\" with the identifier \"" << this->fieldTag().identifier() 
	      << "\" to a type of \"" << type_cast_name
	      << "\" from a boost::any object containing a type of \"" 
	      << Teuchos::demangleName(a.type().name()) << "\"." << std::endl;
    throw;
  }

}

//**********************************************************************
template<typename DataT>
void PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
print(std::ostream& os,	bool printValues) const
{

  os << "MDField(";
  for (size_type i=0; i < m_tag.dataLayout().rank(); ++i) {
    if (i > 0)
      os << ",";
    os << m_tag.dataLayout().dimension(i);
  }
  os << "): ";
  
  os << m_tag;

  if (printValues)
    os << "Error - MDField no longer supports the \"printValues\" member of the MDField::print() method. Values may be on a device that does not support printing (e.g. GPU).  Please disconstinue the use of this call!" << std::endl;  

}
//**********************************************************************
template<typename DataT>
template<typename MDFieldType>
void
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
deep_copy(const MDFieldType& source)
{
  if (m_tag.dataLayout().rank() == 1){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      m_field_data1(ind1) = source(ind1);
  }
  else if (m_tag.dataLayout().rank() == 2){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      for (int ind2=0; ind2<m_tag.dataLayout().dimension(1); ind2++)   
        m_field_data2(ind1,ind2) = source(ind1,ind2);
  }
  else if (m_tag.dataLayout().rank() == 3){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      for (int ind2=0; ind2<m_tag.dataLayout().dimension(1); ind2++)
         for (int ind3=0; ind3<m_tag.dataLayout().dimension(2); ind3++)
            m_field_data3(ind1,ind2,ind3) = source(ind1,ind2,ind3);
  }
  else if (m_tag.dataLayout().rank() == 4){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      for (int ind2=0; ind2<m_tag.dataLayout().dimension(1); ind2++)
         for (int ind3=0; ind3<m_tag.dataLayout().dimension(2); ind3++)
            for (int ind4=0; ind4<m_tag.dataLayout().dimension(3); ind4++)  
                m_field_data4(ind1,ind2,ind3,ind4) = source(ind1,ind2,ind3,ind4);   
  }
  else if (m_tag.dataLayout().rank() == 5){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      for (int ind2=0; ind2<m_tag.dataLayout().dimension(1); ind2++)
         for (int ind3=0; ind3<m_tag.dataLayout().dimension(2); ind3++)
            for (int ind4=0; ind4<m_tag.dataLayout().dimension(3); ind4++)
                for (int ind5=0; ind5<m_tag.dataLayout().dimension(4); ind5++)
                   m_field_data5(ind1,ind2,ind3,ind4,ind5) = source(ind1,ind2,ind3,ind4,ind5);
  }
  else if (m_tag.dataLayout().rank() == 6){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      for (int ind2=0; ind2<m_tag.dataLayout().dimension(1); ind2++)
         for (int ind3=0; ind3<m_tag.dataLayout().dimension(2); ind3++)
            for (int ind4=0; ind4<m_tag.dataLayout().dimension(3); ind4++)
                for (int ind5=0; ind5<m_tag.dataLayout().dimension(4); ind5++)
                   for (int ind6=0; ind6<m_tag.dataLayout().dimension(5); ind6++)
                      m_field_data6(ind1,ind2,ind3,ind4,ind5,ind6) = source(ind1,ind2,ind3,ind4,ind5,ind6);
  }
  else if (m_tag.dataLayout().rank() == 7){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      for (int ind2=0; ind2<m_tag.dataLayout().dimension(1); ind2++)
         for (int ind3=0; ind3<m_tag.dataLayout().dimension(2); ind3++)
            for (int ind4=0; ind4<m_tag.dataLayout().dimension(3); ind4++)
                for (int ind5=0; ind5<m_tag.dataLayout().dimension(4); ind5++)
                   for (int ind6=0; ind6<m_tag.dataLayout().dimension(5); ind6++)
                      for (int ind7=0; ind7<m_tag.dataLayout().dimension(6); ind7++)
                          m_field_data7(ind1,ind2,ind3,ind4,ind5,ind6,ind7) = source(ind1,ind2,ind3,ind4,ind5,ind6,ind7);
  }
}
//******************************************************************************
template<typename DataT>
template<typename MDFieldType>
void
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
V_Multiply(const MDFieldType& source)
{
  if (m_tag.dataLayout().rank() == 1){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      m_field_data1(ind1) = m_field_data1(ind1)*source(ind1);
  }
  else if (m_tag.dataLayout().rank() == 2){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      for (int ind2=0; ind2<m_tag.dataLayout().dimension(1); ind2++)
        m_field_data2(ind1,ind2) = m_field_data2(ind1,ind2)*source(ind1,ind2);
  }
  else if (m_tag.dataLayout().rank() == 3){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      for (int ind2=0; ind2<m_tag.dataLayout().dimension(1); ind2++)
         for (int ind3=0; ind3<m_tag.dataLayout().dimension(2); ind3++)
            m_field_data3(ind1,ind2,ind3) = m_field_data3(ind1,ind2,ind3)*source(ind1,ind2,ind3);
  }
  else if (m_tag.dataLayout().rank() == 4){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      for (int ind2=0; ind2<m_tag.dataLayout().dimension(1); ind2++)
         for (int ind3=0; ind3<m_tag.dataLayout().dimension(2); ind3++)
            for (int ind4=0; ind4<m_tag.dataLayout().dimension(3); ind4++)
                m_field_data4(ind1,ind2,ind3,ind4) = m_field_data4(ind1,ind2,ind3,ind4)*source(ind1,ind2,ind3,ind4);
  }
  else if (m_tag.dataLayout().rank() == 5){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      for (int ind2=0; ind2<m_tag.dataLayout().dimension(1); ind2++)
         for (int ind3=0; ind3<m_tag.dataLayout().dimension(2); ind3++)
            for (int ind4=0; ind4<m_tag.dataLayout().dimension(3); ind4++)
                for (int ind5=0; ind5<m_tag.dataLayout().dimension(4); ind5++)
                   m_field_data5(ind1,ind2,ind3,ind4,ind5) = m_field_data5(ind1,ind2,ind3,ind4,ind5)*source(ind1,ind2,ind3,ind4,ind5);
  }
  else if (m_tag.dataLayout().rank() == 6){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      for (int ind2=0; ind2<m_tag.dataLayout().dimension(1); ind2++)
         for (int ind3=0; ind3<m_tag.dataLayout().dimension(2); ind3++)
            for (int ind4=0; ind4<m_tag.dataLayout().dimension(3); ind4++)
                for (int ind5=0; ind5<m_tag.dataLayout().dimension(4); ind5++)
                   for (int ind6=0; ind6<m_tag.dataLayout().dimension(5); ind6++)
                      m_field_data6(ind1,ind2,ind3,ind4,ind5,ind6) = m_field_data6(ind1,ind2,ind3,ind4,ind5,ind6)*source(ind1,ind2,ind3,ind4,ind5,ind6);
  }
  else if (m_tag.dataLayout().rank() == 7){
   for (int ind1=0; ind1<m_tag.dataLayout().dimension(0); ind1++)
      for (int ind2=0; ind2<m_tag.dataLayout().dimension(1); ind2++)
         for (int ind3=0; ind3<m_tag.dataLayout().dimension(2); ind3++)
            for (int ind4=0; ind4<m_tag.dataLayout().dimension(3); ind4++)
                for (int ind5=0; ind5<m_tag.dataLayout().dimension(4); ind5++)
                   for (int ind6=0; ind6<m_tag.dataLayout().dimension(5); ind6++)
                      for (int ind7=0; ind7<m_tag.dataLayout().dimension(6); ind7++)
                          m_field_data7(ind1,ind2,ind3,ind4,ind5,ind6,ind7) = m_field_data7(ind1,ind2,ind3,ind4,ind5,ind6,ind7)*source(ind1,ind2,ind3,ind4,ind5,ind6,ind7);
  }
}

//******************************************************************************
template<typename DataT>
void
PHX::MDField<DataT,void,void,void,void,void,void,void,void>::
deep_copy(const DataT source)
{
  if (m_tag.dataLayout().rank() == 1){
    Kokkos::deep_copy(m_field_data1, source);  
  }
  else if (m_tag.dataLayout().rank() == 2){
    Kokkos::deep_copy(m_field_data2, source);
  }
  else if (m_tag.dataLayout().rank() == 3){
    Kokkos::deep_copy(m_field_data3, source);
  }
  else if (m_tag.dataLayout().rank() == 4){
   Kokkos::deep_copy(m_field_data4, source);
  }
  else if (m_tag.dataLayout().rank() == 5){
   Kokkos::deep_copy(m_field_data5, source);
  }
  else if (m_tag.dataLayout().rank() == 6){
   Kokkos::deep_copy(m_field_data6, source);
  }
  else if (m_tag.dataLayout().rank() == 7){
   Kokkos::deep_copy(m_field_data7, source);
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
// For interoperability with Intrepid+Kokkos

//template<class A>
//struct Rank{static const int value = -1;};
#ifdef HAVE_INTREPID_KOKKOSCORE
#ifdef Phalanx_ENABLE_Intrepid

template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7>
struct Rank <PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> > {
 static const int value=PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::ArrayRank;

};

template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7>
struct Rank <const PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> > {
 static const int value=PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::ArrayRank;

};

template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7, class ScalarT>
struct Return_Type <const PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> ,ScalarT> {
// static const int value=PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::ArrayRank;
typedef typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type>::return_type return_type;
typedef typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type>::return_type const_return_type;
};

template<typename DataT,
         typename Tag0,typename Tag1, typename Tag2, typename Tag3,
         typename Tag4,typename Tag5, typename Tag6, typename Tag7, class ScalarT>
struct Return_Type < PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>, ScalarT> {
// static const int value=PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::ArrayRank;
 typedef typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type>::return_type return_type;
 typedef typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::array_type>::return_type const_return_type;
 };


// ***********************
// Runtime
// ***********************
template<typename DataT>
struct Rank <PHX::MDField<DataT> > {
 static const int value = -1;
};

template<typename DataT>
struct Rank <const PHX::MDField<DataT> > {
  static const int value = -1;
};

template<typename DataT, class ScalarT>
struct Return_Type <const PHX::MDField<DataT> ,ScalarT> {
typedef typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT>::array_type>::return_type return_type;
typedef typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT>::array_type>::return_type const_return_type;
};

template<typename DataT, class ScalarT>
struct Return_Type < PHX::MDField<DataT>, ScalarT> {
 typedef typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT>::array_type>::return_type return_type;
 typedef typename PHX::MDFieldTypeTraits<typename PHX::MDField<DataT>::array_type>::return_type const_return_type;
 };

#endif
#endif
//********************************************************************************************




#endif
