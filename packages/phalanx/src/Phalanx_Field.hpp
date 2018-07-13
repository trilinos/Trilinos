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


#ifndef PHX_FIELD_HPP
#define PHX_FIELD_HPP

#include <iostream>
#include <string>
#include <type_traits>
#include "Phalanx_config.hpp"
#include "Phalanx_any.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Kokkos_View.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_MDField_TypeTraits.hpp"
#include "Sacado.hpp"

namespace PHX {

  class DataLayout;
  class FieldTag;

  template<typename DataT, int Rank>
  struct KokkosDimType;

  template<typename DataT>
  struct KokkosDimType<DataT,1> {
    typedef DataT* type;
  };
  template<typename DataT>
  struct KokkosDimType<DataT,2> {
    typedef DataT** type;
  };
  template<typename DataT>
  struct KokkosDimType<DataT,3> {
    typedef DataT*** type;
  };
  template<typename DataT>
  struct KokkosDimType<DataT,4> {
    typedef DataT**** type;
  };
  template<typename DataT>
  struct KokkosDimType<DataT,5> {
    typedef DataT***** type;
  };
  template<typename DataT>
  struct KokkosDimType<DataT,6> {
    typedef DataT****** type;
  };
  template<typename DataT>
  struct KokkosDimType<DataT,7> {
    typedef DataT******* type;
  };
  template<typename DataT>
  struct KokkosDimType<DataT,8> {
    typedef DataT******** type;
  };

  // *************************************
  // Compile time checked MDField
  // *************************************

  template<typename DataT,int Rank>
  class Field {

  public:

    typedef DataT value_type;
    typedef DataT& reference_type;

    typedef typename KokkosDimType<DataT,Rank>::type kokkos_data_type;
    typedef typename PHX::View<kokkos_data_type> array_type;
    typedef typename array_type::array_layout layout_type;
    typedef typename array_type::device_type device_type;
    typedef typename PHX::Device::size_type size_type;
    typedef typename array_type::execution_space execution_space;

    Field(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& dl);

    Field(const PHX::FieldTag& t);

    Field(const Teuchos::RCP<const PHX::FieldTag>& t);

    Field();

    //! For const/non-const compatibility
    template<typename CopyDataT>
    Field(const Field<CopyDataT,Rank>& source);

    ~Field();

    static const int ArrayRank=array_type::Rank;

    const PHX::FieldTag& fieldTag() const;

    Teuchos::RCP<const PHX::FieldTag> fieldTagPtr() const;

    //! For const/non-const compatibility
    template<typename CopyDataT>
    PHX::Field<DataT,Rank>&
    operator=(const Field<CopyDataT,Rank>& source);

    template<typename... index_pack>
    KOKKOS_INLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(const index_pack&...) const;

    KOKKOS_INLINE_FUNCTION
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

    KOKKOS_INLINE_FUNCTION
    size_type size() const;

    void setFieldTag(const PHX::FieldTag& t);

    void setFieldTag(const Teuchos::RCP<const PHX::FieldTag>& t);

    void setFieldData(const PHX::any& a);

    void print(std::ostream& os, bool printValues = false) const;

    KOKKOS_INLINE_FUNCTION
    Kokkos::DynRankView<DataT,typename PHX::DevLayout<DataT>::type,PHX::Device> get_view();

    KOKKOS_INLINE_FUNCTION
    const Kokkos::DynRankView<DataT,typename PHX::DevLayout<DataT>::type,PHX::Device> get_view() const;

    //! Returns a static view of the underlying kokkos static view.
    KOKKOS_INLINE_FUNCTION
    array_type get_static_view();

    //! Returns a static view of the underlying kokkos static view.
    KOKKOS_INLINE_FUNCTION
    const array_type get_static_view() const;

    template<typename SrcDataT>
    void deep_copy(const PHX::Field<SrcDataT,Rank>& source);

    void deep_copy(const DataT source);

  private:

    Teuchos::RCP<const PHX::FieldTag> m_tag;

    array_type m_field_data;

#ifdef PHX_DEBUG
    bool m_data_set;
    static const std::string m_field_tag_error_msg;
    static const std::string m_field_data_error_msg;
#endif

    //! For copy/assignment between const/non-const
    template<typename ScalarT,int FriendRank> friend class PHX::Field;
  };

  template<typename DataT,int Rank>
  std::ostream& operator<<(std::ostream& os,
			   const PHX::Field<DataT,Rank>& h);
}

#include "Phalanx_Field_Def.hpp"

#endif
