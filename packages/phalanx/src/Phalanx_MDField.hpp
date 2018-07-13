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
#include <type_traits>
#include "Phalanx_any.hpp"
#include "Teuchos_RCP.hpp"
#include "Kokkos_View.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_MDField_TypeTraits.hpp"
#include "Sacado.hpp"

namespace PHX {

  class DataLayout;
  class FieldTag;

  template<typename DataT,
	   typename Tag0 = void, typename Tag1 = void, typename Tag2 = void, 
	   typename Tag3 = void, typename Tag4 = void, typename Tag5 = void,
	   typename Tag6 = void, typename Tag7 = void>
  class MDField;

  template<typename DataT,
           typename Tag0 = void, typename Tag1 = void, typename Tag2 = void,
           typename Tag3 = void, typename Tag4 = void, typename Tag5 = void,
           typename Tag6 = void, typename Tag7 = void>
  struct KokkosDimensionType{
    typedef DataT******** type;
 };

  template<typename DataT,
           typename Tag0> 
  struct KokkosDimensionType<DataT, Tag0, void, void, void, void, void, void, void> {
    typedef DataT* type;
  };
  template<typename DataT,
           typename Tag0,
           typename Tag1> 
  struct KokkosDimensionType<DataT, Tag0, Tag1, void, void, void, void, void, void> {
    typedef DataT** type;
  };
   template<typename DataT,
           typename Tag0,
           typename Tag1,
           typename Tag2>
  struct KokkosDimensionType<DataT, Tag0, Tag1, Tag2, void, void, void, void, void> {
    typedef DataT*** type;
  };
  template<typename DataT,
           typename Tag0,
           typename Tag1,
           typename Tag2,
           typename Tag3>
  struct KokkosDimensionType<DataT, Tag0, Tag1, Tag2, Tag3, void, void, void, void> {
    typedef DataT**** type;
  };
 template<typename DataT,
           typename Tag0,
           typename Tag1,
           typename Tag2,
           typename Tag3,
           typename Tag4>
  struct KokkosDimensionType<DataT, Tag0, Tag1, Tag2, Tag3, Tag4, void, void, void> {
    typedef DataT***** type;
  };
  template<typename DataT,
           typename Tag0,
           typename Tag1,
           typename Tag2,
           typename Tag3,
           typename Tag4,
           typename Tag5>
  struct KokkosDimensionType<DataT, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, void, void> {
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
  struct KokkosDimensionType<DataT, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, void> {
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
    typedef DataT& reference_type;

    typedef typename KokkosDimensionType<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::type kokkos_data_type;
    typedef typename PHX::View<kokkos_data_type> array_type;
    typedef typename array_type::array_layout layout_type;
    typedef typename array_type::device_type device_type;
    typedef typename PHX::Device::size_type size_type;
    typedef typename array_type::execution_space execution_space;
 
    MDField(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& dl);
    
    MDField(const PHX::FieldTag& t);

    MDField(const Teuchos::RCP<const PHX::FieldTag>& t);

    MDField();

    template<typename CopyDataT,
             typename T0 = void, typename T1 = void, typename T2 = void, 
             typename T3 = void, typename T4 = void, typename T5 = void,
             typename T6 = void, typename T7 = void>
    MDField(const MDField<CopyDataT,T0,T1,T2,T3,T4,T5,T6,T7>& source);

    ~MDField();

    static const int ArrayRank=array_type::Rank;
    
    const PHX::FieldTag& fieldTag() const;

    Teuchos::RCP<const PHX::FieldTag> fieldTagPtr() const;

    template<typename CopyDataT,
             typename T0 = void, typename T1 = void, typename T2 = void, 
             typename T3 = void, typename T4 = void, typename T5 = void,
             typename T6 = void, typename T7 = void>
    PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>&
    operator=(const MDField<CopyDataT,T0,T1,T2,T3,T4,T5,T6,T7>& source);    
    
    template<typename... index_pack>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(const index_pack&...) const;

    KOKKOS_FORCEINLINE_FUNCTION
    size_type rank() const;

    template< typename iType >
    KOKKOS_INLINE_FUNCTION constexpr
    typename std::enable_if< std::is_integral<iType>::value , size_t >::type
    extent( const iType & r ) const
    {return m_field_data.extent(r);}
    
    template< typename iType >
    KOKKOS_INLINE_FUNCTION constexpr
    typename std::enable_if< std::is_integral<iType>::value , int >::type
    extent_int( const iType & r ) const
    {return m_field_data.extent_int(r);}

    //TODO: Consider removing the following dimension* functions
    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_0() const 
    {return m_field_data.extent(0);}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_1() const 
    {return m_field_data.extent(1);}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_2() const 
    {return m_field_data.extent(2);}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_3() const 
    {return m_field_data.extent(3);}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_4() const 
    {return m_field_data.extent(4);}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_5() const 
    {return m_field_data.extent(5);}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_6() const 
    {return m_field_data.extent(6);}

    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension_7() const 
    {return m_field_data.extent(7);}

    template<typename iType>
    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension(const iType& ord) const;

    KOKKOS_FORCEINLINE_FUNCTION
    size_type size() const;

    void setFieldTag(const PHX::FieldTag& t);

    void setFieldTag(const Teuchos::RCP<const PHX::FieldTag>& t);
    
    void setFieldData(const PHX::any& a);
    
    void print(std::ostream& os, bool printValues = false) const;

    /** WARNING: The vector data in this method should be a "size_type" to be consistent with Kokkos, but for backwards compatibility during the transition, needs to be templated in the index type.

     void dimensions(std::vector<size_type>& dims);
    */
    template<typename iType>
    void dimensions(std::vector<iType>& dims);
   
    KOKKOS_FORCEINLINE_FUNCTION 
    Kokkos::DynRankView<DataT,typename PHX::DevLayout<DataT>::type,PHX::Device> get_view();

    KOKKOS_FORCEINLINE_FUNCTION
    const Kokkos::DynRankView<DataT,typename PHX::DevLayout<DataT>::type,PHX::Device> get_view() const;

    //! Returns a static view of the underlying kokkos static view.
    KOKKOS_FORCEINLINE_FUNCTION 
    array_type get_static_view();

    //! Returns a static view of the underlying kokkos static view.
    KOKKOS_FORCEINLINE_FUNCTION
    const array_type get_static_view() const;

    template<typename SrcDataT,
             typename SrcTag0,typename SrcTag1, typename SrcTag2, typename SrcTag3,
             typename SrcTag4,typename SrcTag5, typename SrcTag6, typename SrcTag7>
    void deep_copy(const PHX::MDField<SrcDataT,SrcTag0,SrcTag1,SrcTag2,SrcTag3,
                   SrcTag4,SrcTag5,SrcTag6,SrcTag7>& source);

    void deep_copy(const DataT source);

  private:
    
    Teuchos::RCP<const PHX::FieldTag> m_tag;
    
    array_type m_field_data;

#ifdef PHX_DEBUG
    bool m_data_set;
    static const std::string m_field_tag_error_msg;
    static const std::string m_field_data_error_msg;
#endif

    template<typename ScalarT,
             typename T0, typename T1, typename T2, 
             typename T3, typename T4, typename T5,
             typename T6, typename T7>
    friend class PHX::MDField;
    
  };
  
  template<typename DataT,
	   typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	   typename Tag4, typename Tag5, typename Tag6, typename Tag7>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::MDField<DataT, Tag0, Tag1, 
			   Tag2, Tag3, Tag4, Tag5, Tag6, Tag7>& h);
}

// *************************************
// Runtime time checked MDField
// *************************************

#include "Phalanx_MDField_DynRank.hpp"

// Definition file (this will also bring in the DynRank version definitions as
// well)
#include "Phalanx_MDField_Def.hpp"

#endif 
