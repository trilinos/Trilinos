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


#ifndef PHX_MDFIELD_DYN_RANK_HPP
#define PHX_MDFIELD_DYN_RANK_HPP

#include <iostream>
#include <string>
#include "Sacado.hpp"
#include "Kokkos_DynRankView_Fad.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Teuchos_RCP.hpp"

#include "Phalanx_config.hpp"
#include "Phalanx_any.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_MDField_TypeTraits.hpp"

namespace PHX {

  class DataLayout;
  class FieldTag;
  
  // *************************************
  // Runtime time checked MDField
  // *************************************

  // temporary for bracket op support
  template <typename ViewType>
  KOKKOS_FORCEINLINE_FUNCTION
  unsigned getSacadoSize(const ViewType& view);

  template<typename DataT>
  class MDField<DataT,void,void,void,void,void,void,void,void> {
    
  public:

    typedef DataT value_type;
    typedef DataT& reference_type;
 
    typedef typename Kokkos::DynRankView<DataT,typename PHX::DevLayout<DataT>::type,PHX::Device> array_type;
      
    typedef typename PHX::Device::size_type size_type;

    typedef typename array_type::execution_space execution_space;

    MDField(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& t);
    
    MDField(const PHX::FieldTag& v);

    MDField(const Teuchos::RCP<const PHX::FieldTag>& v);
    
    MDField();

    template<typename CopyDataT,
             typename T0, typename T1, typename T2, 
             typename T3, typename T4, typename T5,
             typename T6, typename T7>
    MDField(const MDField<CopyDataT,T0,T1,T2,T3,T4,T5,T6,T7>& source);
    
    ~MDField();
    
    const PHX::FieldTag& fieldTag() const;

    Teuchos::RCP<const PHX::FieldTag> fieldTagPtr() const;

    template<typename CopyDataT,
             typename T0, typename T1, typename T2, 
             typename T3, typename T4, typename T5,
             typename T6, typename T7>
    PHX::MDField<DataT,void,void,void,void,void,void,void,void>&
    operator=(const MDField<CopyDataT,T0,T1,T2,T3,T4,T5,T6,T7>& source);
    
    // template<typename iType0, typename iType1, typename iType2, typename iType3,
    // 	     typename iType4, typename iType5, typename iType6, typename iType7>
    // KOKKOS_FORCEINLINE_FUNCTION
    // typename PHX::MDFieldTypeTraits<array_type>::return_type
    // operator()(iType0 index0, iType1 index1, iType2 index2, 
    // 	       iType3 index3, iType4 index4, iType5 index5,
    // 	       iType6 index6, iType7 index7) const;

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

    template<typename iType>
    KOKKOS_FORCEINLINE_FUNCTION
    size_type dimension(const iType& ord) const;

    /** WARNING: The vector data in this method should be a "size_type" to be consistent with Kokkos, but for backwards compatibility during the transition, needs to be templated in the index type.

     void dimensions(std::vector<size_type>& dims);
    */
    template<typename iType>
    void dimensions(std::vector<iType>& dims);

    KOKKOS_FORCEINLINE_FUNCTION
    size_type size() const;

    void setFieldTag(const PHX::FieldTag& t);

    void setFieldTag(const Teuchos::RCP<const PHX::FieldTag>& t);
    
    void setFieldData(const PHX::any& a);
    
    void print(std::ostream& os, bool printValues = false) const;

    template<typename MDFieldType>
    void deep_copy(const MDFieldType& source);

    void deep_copy(const DataT source);
    
    template<typename MDFieldTypeA, typename MDFieldTypeB, unsigned int RANK>
    struct V_MultiplyFunctor{
      V_MultiplyFunctor(const MDFieldTypeA &base, const MDFieldTypeB &source) :base_(base), source_(source){}
      KOKKOS_INLINE_FUNCTION
      void operator()(const PHX::index_t& i) const;
      MDFieldTypeA base_;
      MDFieldTypeB source_;
    };

  public:

    template<typename MDFieldType>
    void V_Multiply(const MDFieldType& source);
 
    KOKKOS_FORCEINLINE_FUNCTION
    array_type get_view()
    {return m_field_data;}

    KOKKOS_FORCEINLINE_FUNCTION
    const array_type get_view() const
    {return m_field_data;}

    PHX::any& get_static_any_view()
    {return m_any;}

    const PHX::any& get_static_any_view() const
    {return m_any;}
   
  private:
   
    Teuchos::RCP<const PHX::FieldTag> m_tag;  
    PHX::any m_any; //! Store RCP to Kokkos::View
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
  
  template<typename DataT>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::MDField<DataT, void, void, 
			   void, void, void, void, void, void>& h);
  
} 

#endif 
