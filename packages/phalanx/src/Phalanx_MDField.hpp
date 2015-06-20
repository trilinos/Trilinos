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


#ifndef PHX_MDFIELD_HPP
#define PHX_MDFIELD_HPP

#include <iostream>
#include <string>
#include <boost/any.hpp>
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Kokkos_View.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_MDFieldToKokkos.hpp"
#include "Phalanx_MDField_TypeTraits.hpp"
#include "Sacado.hpp"

namespace PHX {

  template<typename DataT,
	   typename Tag0 = void, typename Tag1 = void, typename Tag2 = void, 
	   typename Tag3 = void, typename Tag4 = void, typename Tag5 = void,
	   typename Tag6 = void, typename Tag7 = void>
  class MDField;

  template<typename DataT,
           typename Tag0 = void, typename Tag1 = void, typename Tag2 = void,
           typename Tag3 = void, typename Tag4 = void, typename Tag5 = void,
           typename Tag6 = void, typename Tag7 = void>
  struct KokkosDimentionType{
    typedef DataT******** type;
 };

  template<typename DataT,
           typename Tag0> 
  struct KokkosDimentionType<DataT, Tag0, void, void, void, void, void, void, void> {
    typedef DataT* type;
  };
  template<typename DataT,
           typename Tag0,
           typename Tag1> 
  struct KokkosDimentionType<DataT, Tag0, Tag1, void, void, void, void, void, void> {
    typedef DataT** type;
  };
   template<typename DataT,
           typename Tag0,
           typename Tag1,
           typename Tag2>
  struct KokkosDimentionType<DataT, Tag0, Tag1, Tag2, void, void, void, void, void> {
    typedef DataT*** type;
  };
  template<typename DataT,
           typename Tag0,
           typename Tag1,
           typename Tag2,
           typename Tag3>
  struct KokkosDimentionType<DataT, Tag0, Tag1, Tag2, Tag3, void, void, void, void> {
    typedef DataT**** type;
  };
 template<typename DataT,
           typename Tag0,
           typename Tag1,
           typename Tag2,
           typename Tag3,
           typename Tag4>
  struct KokkosDimentionType<DataT, Tag0, Tag1, Tag2, Tag3, Tag4, void, void, void> {
    typedef DataT***** type;
  };
  template<typename DataT,
           typename Tag0,
           typename Tag1,
           typename Tag2,
           typename Tag3,
           typename Tag4,
           typename Tag5>
  struct KokkosDimentionType<DataT, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, void, void> {
    typedef DataT****** type;
 };
 template<typename DataT,
           typename Tag0,
           typename Tag1,
           typename Tag2,
           typename Tag3,
           typename Tag4,
           typename Tag5,
           typename Tag6>
  struct KokkosDimentionType<DataT, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, void> {
    typedef DataT******* type;
 };

  // *************************************
  // Compile time checked MDField
  // *************************************

  template<typename DataT,
	   typename Tag0, typename Tag1, typename Tag2, 
	   typename Tag3, typename Tag4, typename Tag5,
	   typename Tag6, typename Tag7>
  class MDField {
    
  public:

    typedef DataT value_type;

    typedef typename KokkosDimentionType<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::type kokkos_data_type;
    typedef typename Kokkos::View <kokkos_data_type, PHX::Device> array_type;
    typedef typename array_type::array_layout layout_type;
    typedef typename array_type::device_type device_type;
    typedef typename PHX::Device::size_type size_type;
 
    KOKKOS_FORCEINLINE_FUNCTION
    MDField(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& t);
    
    KOKKOS_FORCEINLINE_FUNCTION
    MDField(const PHX::Tag<DataT>& v);
    
    KOKKOS_FORCEINLINE_FUNCTION
    MDField();
    
    KOKKOS_FORCEINLINE_FUNCTION
    ~MDField();

    static const int ArrayRank=array_type::Rank;
    
    KOKKOS_FORCEINLINE_FUNCTION
    const PHX::FieldTag& fieldTag() const;

    template<typename iType0, typename iType1, typename iType2, typename iType3,
	     typename iType4, typename iType5, typename iType6, typename iType7>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1, iType2 index2, 
	       iType3 index3, iType4 index4, iType5 index5,
	       iType6 index6, iType7 index7) const;

    template<typename iType0, typename iType1, typename iType2, typename iType3,
	     typename iType4, typename iType5, typename iType6>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1, iType2 index2, 
	       iType3 index3, iType4 index4, iType5 index5,
	       iType6 index6) const;

    template<typename iType0, typename iType1, typename iType2, typename iType3,
	     typename iType4, typename iType5>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1, iType2 index2, 
	       iType3 index3, iType4 index4, iType5 index5)const;
    
    template<typename iType0, typename iType1, typename iType2, typename iType3,
	     typename iType4>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1, iType2 index2, 
	       iType3 index3, iType4 index4)const;
    
    template<typename iType0, typename iType1, typename iType2, typename iType3>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1, iType2 index2, 
	       iType3 index3)const;

    template<typename iType0, typename iType1, typename iType2>
    KOKKOS_FORCEINLINE_FUNCTION    
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1, iType2 index2)const;
    
    template<typename iType0, typename iType1>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1)const;
    
    template<typename iType0>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0) const;

    KOKKOS_FORCEINLINE_FUNCTION
    size_type rank() const;

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_0() const 
    {return m_field_data.dimension_0();}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_1() const 
    {return m_field_data.dimension_1();}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_2() const 
    {return m_field_data.dimension_2();}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_3() const 
    {return m_field_data.dimension_3();}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_4() const 
    {return m_field_data.dimension_4();}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_5() const 
    {return m_field_data.dimension_5();}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_6() const 
    {return m_field_data.dimension_6();}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_7() const 
    {return m_field_data.dimension_7();}

    template<typename iType>
    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension(const iType& ord) const;

    KOKKOS_FORCEINLINE_FUNCTION
    size_type size() const;

    void setFieldTag(const PHX::Tag<DataT>& t);
    
    void setFieldData(const boost::any& a);
    
    void print(std::ostream& os, bool printValues = false) const;

    /** WARNING: The vector data in this method should be a "size_type" to be consistent with Kokkos, but for backwards compatibility during the transition, needs to be templated in the index type.

     void dimensions(std::vector<size_type>& dims);
    */
    template<typename iType>
    void dimensions(std::vector<iType>& dims);
   
    KOKKOS_FORCEINLINE_FUNCTION 
    array_type get_kokkos_view();

    KOKKOS_FORCEINLINE_FUNCTION
    const array_type get_kokkos_view()const;

    template<typename MDFieldType>
    void deep_copy(const MDFieldType& source);

    void deep_copy(const DataT source);

  private:
    
    PHX::Tag<DataT> m_tag;
    
    array_type m_field_data;

#ifdef PHX_DEBUG
    bool m_tag_set;
    bool m_data_set;
    static const std::string m_field_tag_error_msg;
    static const std::string m_field_data_error_msg;
#endif

  };
  
  template<typename DataT,
	   typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	   typename Tag4, typename Tag5, typename Tag6, typename Tag7>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::MDField<DataT, Tag0, Tag1, 
			   Tag2, Tag3, Tag4, Tag5, Tag6, Tag7>& h);
  
  // *************************************
  // Runtime time checked MDField
  // *************************************

  // temporary for bracket op support
  template <typename T, typename L, typename D, typename M, typename S>
  KOKKOS_FORCEINLINE_FUNCTION 
  unsigned getSacadoSize(const Kokkos::View<T,L,D,M,S>& view);
  
  template <typename T, typename L, typename D, typename M>
  KOKKOS_FORCEINLINE_FUNCTION 
  unsigned getSacadoSize(const Kokkos::View<T,L,D,M,Kokkos::Impl::ViewSpecializeSacadoFad>& view);

  template<typename DataT>
  class MDField<DataT,void,void,void,void,void,void,void,void> {
    
  public:

    typedef DataT value_type;
 
    typedef typename Kokkos::View <DataT*******, PHX::Device> array_type;
      
    typedef typename PHX::Device::size_type size_type;

    KOKKOS_FORCEINLINE_FUNCTION
    MDField(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& t);
    
    KOKKOS_FORCEINLINE_FUNCTION
    MDField(const PHX::Tag<DataT>& v);
    
    KOKKOS_FORCEINLINE_FUNCTION
    MDField();
    
    KOKKOS_FORCEINLINE_FUNCTION
    ~MDField();
    
    const PHX::FieldTag& fieldTag() const;
    
    template<typename iType0, typename iType1, typename iType2, typename iType3,
	     typename iType4, typename iType5, typename iType6, typename iType7>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1, iType2 index2, 
	       iType3 index3, iType4 index4, iType5 index5,
	       iType6 index6, iType7 index7) const;

    template<typename iType0, typename iType1, typename iType2, typename iType3,
	     typename iType4, typename iType5, typename iType6>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1, iType2 index2, 
	       iType3 index3, iType4 index4, iType5 index5,
	       iType6 index6) const;

    template<typename iType0, typename iType1, typename iType2, typename iType3,
	     typename iType4, typename iType5>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1, iType2 index2, 
	       iType3 index3, iType4 index4, iType5 index5) const;
    
    template<typename iType0, typename iType1, typename iType2, typename iType3,
	     typename iType4>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1, iType2 index2, 
	       iType3 index3, iType4 index4) const;
    
    template<typename iType0, typename iType1, typename iType2, typename iType3>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1, iType2 index2, 
	       iType3 index3) const;
    
    template<typename iType0, typename iType1, typename iType2>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1, iType2 index2) const;
   
    template<typename iType0, typename iType1>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0, iType1 index1) const;
    
    template<typename iType0>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(iType0 index0) const;

    template<typename iType0>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator[](iType0 index0) const;

    KOKKOS_FORCEINLINE_FUNCTION
    size_type rank() const;

    template<typename iType>
    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension(const iType& ord) const;

    /** WARNING: The vector data in this method should be a "size_type" to be consistent with Kokkos, but for backwards compatibility during the transition, needs to be templated in the index type.

     void dimensions(std::vector<size_type>& dims);
    */
    template<typename iType>
    KOKKOS_FORCEINLINE_FUNCTION
    void dimensions(std::vector<iType>& dims);

    KOKKOS_FORCEINLINE_FUNCTION
    size_type size() const;

    void setFieldTag(const PHX::Tag<DataT>& t);
    
    void setFieldData(const boost::any& a);
    
    void print(std::ostream& os, bool printValues = false) const;

    template<typename MDFieldType>
    void deep_copy(const MDFieldType& source);

    void deep_copy(const DataT source);
    
    template<typename MDFieldType>
    void V_Multiply(const MDFieldType& source);
 
    private:
   
    PHX::Tag<DataT> m_tag;  
    typedef Kokkos::View<DataT*, PHX::Device> array_type1;
    typedef Kokkos::View<DataT**, PHX::Device> array_type2;
    typedef Kokkos::View<DataT***, PHX::Device> array_type3;
    typedef Kokkos::View<DataT****, PHX::Device> array_type4;
    typedef Kokkos::View<DataT*****, PHX::Device> array_type5;
    typedef Kokkos::View<DataT******, PHX::Device> array_type6;
    typedef Kokkos::View<DataT*******, PHX::Device> array_type7;
 
    Kokkos::View<DataT*, PHX::Device> m_field_data1;
    Kokkos::View<DataT**, PHX::Device> m_field_data2;
    Kokkos::View<DataT***, PHX::Device> m_field_data3;
    Kokkos::View<DataT****, PHX::Device> m_field_data4;
    Kokkos::View<DataT*****, PHX::Device> m_field_data5;
    Kokkos::View<DataT******, PHX::Device> m_field_data6;
    Kokkos::View<DataT*******, PHX::Device> m_field_data7;

    typedef Kokkos::View<DataT*, PHX::Device, Kokkos::MemoryUnmanaged> array_oned_type; 
    array_oned_type m_field_oned_view;

    //! For fast access to rank/dimension/size data.  Entries 0-6 are dimensions, entry 7 is rank, entry 8 is size. 
    Kokkos::View<PHX::index_size_type*, PHX::Device> m_dimension_rank_size;

#ifdef PHX_DEBUG
    bool m_tag_set;
    bool m_data_set;
    static const std::string m_field_tag_error_msg;
    static const std::string m_field_data_error_msg;
#endif

  };
  
  template<typename DataT>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::MDField<DataT, void, void, 
			   void, void, void, void, void, void>& h);
  
} 

#include "Phalanx_MDField_Def.hpp"

#endif 
