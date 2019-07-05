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
#include "Kokkos_DynRankView_Fad.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_MDField_TypeTraits.hpp"
#include "Sacado.hpp"
#include "Phalanx_FieldTag_Tag.hpp"

namespace PHX {

  class DataLayout;
  class FieldTag;

  // ****************************
  // Layouts
  // ****************************
  template<typename L> struct is_layout : std::false_type {};
  template<> struct is_layout<Kokkos::LayoutRight> : std::true_type {};
  template<> struct is_layout<Kokkos::LayoutLeft> : std::true_type {};

  // ****************************
  // Extents
  // ****************************
  template<typename I> struct is_extent : std::false_type {};

  // ****************************
  // Devices
  // ****************************
  template<typename T> struct is_device : std::false_type {};
  template<> struct is_device<PHX::Device> : std::true_type {};

  // ****************************
  // Rank count
  // ****************************
  template<typename...Props> struct RankCount;

  template<typename...Props> struct RankCount : std::integral_constant<int,RankCount<void,Props...>::value> {};

  template<> struct RankCount<void> : std::integral_constant<int,0>{};

  template<typename Extent, typename...Props>
  struct RankCount<typename std::enable_if<is_extent<Extent>::value>::type,Extent,Props...>
    : std::integral_constant<int,1+RankCount<void,Props...>::value> {};

  template<typename NonExtent, typename...Props>
  struct RankCount<typename std::enable_if<!is_extent<NonExtent>::value>::type,NonExtent,Props...>
    : std::integral_constant<int,RankCount<void,Props...>::value> {};

  // ****************************
  // Add pointer (used to construct the static data type)
  // ****************************
  template<typename Data,int N> struct add_pointer;

  template<typename Data,int N> struct add_pointer
  { using type = typename add_pointer<Data*,N-1>::type; };

  template<typename Data> struct add_pointer<Data,0>
  { using type = Data; };

  // ****************************
  // ArrayType
  // ****************************
  template<typename Scalar,int Rank,typename...Props> struct ArrayType;

  // Static rank default
  template<typename Scalar,int Rank,typename...Props> struct ArrayType
  {
    static_assert(Rank > 0,"Error: Rank of static MDField must be greater than zero!");
    static_assert(Rank < 8,"Error: Rank of static MDField must be less than eight!");
    using data_type = typename add_pointer<Scalar,Rank>::type;
    using array_type = Kokkos::View<data_type,Props...>;
  };

  // Dynamic rank specialization
  template<typename Scalar,typename...Props> struct ArrayType<Scalar,0,Props...>
  {
    using data_type = Scalar;
    using array_type = Kokkos::DynRankView<Scalar,Props...>;
  };

  // ****************************
  // FieldTraits
  // ****************************
  template<typename Scalar, typename... Props> struct FieldTraits;

  template<>
  struct FieldTraits<void>
  {
    using layout = void;
    using device = void;
  };

  template<typename Extent, typename... Props>
  struct FieldTraits< typename std::enable_if<is_extent<Extent>::value>::type, Extent, Props...>
  {
    using layout = typename FieldTraits<void,Props...>::layout;
    using device = typename FieldTraits<void,Props...>::device;
  };

  template<typename Layout, typename... Props>
  struct FieldTraits< typename std::enable_if<is_layout<Layout>::value>::type, Layout, Props...>
  {
    using layout = Layout;
    using device = typename FieldTraits<void,Props...>::device;
  };

  template<typename Device, typename... Props>
  struct FieldTraits< typename std::enable_if<is_device<Device>::value>::type, Device, Props...>
  {
    using layout = void;
    using device = Device;
  };

  template<typename Scalar, typename... Props>
  struct FieldTraits {
    using prop = FieldTraits<void,Props...>;
    static constexpr int rank = RankCount<Props...>::value;
    // This sets defaults if not specified
    using layout = typename std::conditional< !std::is_same<typename prop::layout, void>::value,typename prop::layout, typename PHX::DevLayout<Scalar>::type>::type;
    using device = typename std::conditional< !std::is_same<typename prop::device, void>::value,typename prop::device, PHX::Device>::type;
    using data_type = typename ArrayType<Scalar,rank,layout,device>::data_type;
    using array_type = typename ArrayType<Scalar,rank,layout,device>::array_type;
  };

  // ****************************
  // MDField
  // ****************************

  template<typename Scalar,typename... Props>
  class MDField
  {
    // Trick to allow for member method partial specialization on View type (static or dynamic)
    template<int R> struct ViewSpecialization{};

  public:
    using traits = FieldTraits<Scalar,Props...>;
    using layout_type = typename traits::layout;
    using device_type = typename traits::device;
    using data_type = typename traits::data_type;
    using array_type = typename traits::array_type;
    using size_type = typename device_type::size_type;
    using execution_space = typename array_type::execution_space;

    typedef Scalar value_type;
    typedef Scalar& reference_type;


    // Not allowed - too easy to forget to bind memory!
    // template<typename...Extents>
    // MDField(const std::string name,Extents... e)
    //   : view(name,e...)
    // {}

    MDField(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& dl)
#ifdef PHX_DEBUG
      : m_data_set(false)
#endif
    {m_tag = Teuchos::rcp(new PHX::Tag<Scalar>(name,dl));}

    MDField(const PHX::FieldTag& t)
      : m_tag(t.clone())
#ifdef PHX_DEBUG
      , m_data_set(false)
#endif
    {}

    MDField(const Teuchos::RCP<const PHX::FieldTag>& t)
      : m_tag(t)
#ifdef PHX_DEBUG
      , m_data_set(false)
#endif
    {}

    MDField()
#ifdef PHX_DEBUG
      : m_data_set(false)
#endif
    {}

    template<typename SourceScalar,typename...SourceProps>
    MDField(const MDField<SourceScalar,SourceProps...>& source)
      : m_tag(source.m_tag),
        m_view(source.m_view)
#ifdef PHX_DEBUG
      , m_data_set(source.m_data_set)
#endif
    {
      static_assert(std::is_same<typename std::decay<Scalar>::type, typename std::decay<SourceScalar>::type>::value,
                    "ERROR: Compiletime MDField copy ctor requires scalar types to be the same!");
    }

    ~MDField() {}

    constexpr bool is_static() const {return (traits::rank != 0);}
    constexpr bool is_dynamic() const {return (traits::rank == 0);}

    KOKKOS_FORCEINLINE_FUNCTION
    constexpr size_type rank() const {return rank(ViewSpecialization<traits::rank>());}

    KOKKOS_FORCEINLINE_FUNCTION
    constexpr size_t size() const {return m_view.size();}

    const PHX::FieldTag& fieldTag() const
    {
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
      TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, fieldTagErrorMsg());
#endif
      return *m_tag;
    }

    Teuchos::RCP<const PHX::FieldTag> fieldTagPtr() const
    {
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
      TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, fieldTagErrorMsg());
#endif
      return m_tag;
    }

    template<typename SrcScalar,typename...SrcProps>
    PHX::MDField<Scalar,Props...>&
    operator=(const MDField<SrcScalar,SrcProps...>& source)
    {
      m_tag = source.m_tag;
      m_view = source.m_view;
#ifdef PHX_DEBUG
      m_data_set = source.m_data_set;
#endif
      static_assert(std::is_same<typename std::decay<Scalar>::type, typename std::decay<SrcScalar>::type>::value,
                    "ERROR: MDField assignment operator requires scalar types to be the same!");

      return *this;
    }

    template<typename... index_pack>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator()(const index_pack&... indices) const
    {
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
      TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, fieldDataErrorMsg());
#endif
      return m_view(indices...);
    }

    template<typename iType0>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldTypeTraits<array_type>::return_type
    operator[](iType0 index0) const
    {return m_view[index0];}
    
    template< typename iType >
    KOKKOS_INLINE_FUNCTION constexpr
    typename std::enable_if< std::is_integral<iType>::value , size_t >::type
    extent( const iType & r ) const
    {return m_view.extent(r);}

    template< typename iType >
    KOKKOS_INLINE_FUNCTION constexpr
    typename std::enable_if< std::is_integral<iType>::value , int >::type
    extent_int( const iType & r ) const
    {return static_cast<int>(m_view.extent(r));}

    //TODO: Consider removing the following dimension* functions
    KOKKOS_FORCEINLINE_FUNCTION
    constexpr size_type dimension_0() const
    {return m_view.extent(0);}

    KOKKOS_FORCEINLINE_FUNCTION
    constexpr size_type dimension_1() const
    {return m_view.extent(1);}

    KOKKOS_FORCEINLINE_FUNCTION
    constexpr size_type dimension_2() const
    {return m_view.extent(2);}

    KOKKOS_FORCEINLINE_FUNCTION
    constexpr size_type dimension_3() const
    {return m_view.extent(3);}

    KOKKOS_FORCEINLINE_FUNCTION
    constexpr size_type dimension_4() const
    {return m_view.extent(4);}

    KOKKOS_FORCEINLINE_FUNCTION
    constexpr size_type dimension_5() const
    {return m_view.extent(5);}

    KOKKOS_FORCEINLINE_FUNCTION
    constexpr size_type dimension_6() const
    {return m_view.extent(6);}

    KOKKOS_FORCEINLINE_FUNCTION
    constexpr size_type dimension_7() const
    {return m_view.extent(7);}

    template<typename iType>
    KOKKOS_FORCEINLINE_FUNCTION
    constexpr size_type dimension(const iType& ord) const
    {return m_view.extent(ord);}

    void setFieldTag(const PHX::FieldTag& t)
    {m_tag = t.clone();}

    void setFieldTag(const Teuchos::RCP<const PHX::FieldTag>& t)
    {m_tag = t;}

    void setFieldData(const PHX::any& a)
    {setFieldData(ViewSpecialization<traits::rank>(),a);}
    
    void print(std::ostream& os, bool printValues = false) const
    {print(ViewSpecialization<traits::rank>(),os,printValues);}

    /** WARNING: The vector data in this method should be a
        "size_type" to be consistent with Kokkos, but for backwards
        compatibility during the transition, needs to be templated in
        the index type.

        void dimensions(std::vector<size_type>& dims);
    */
    template<typename iType>
    void dimensions(std::vector<iType>& dims)
    {
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
      TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, fieldTagErrorMsg());
      TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, fieldDataErrorMsg());
#endif 
      dims.resize(this->rank());
      for (size_type i=0; i <  this->rank(); ++i) 
        dims[i] = static_cast<iType>(m_view.extent(i));  // dangerous
    }

    KOKKOS_FORCEINLINE_FUNCTION
    Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,PHX::Device> get_view()
    {return m_view;}
    
    KOKKOS_FORCEINLINE_FUNCTION
    const Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,PHX::Device> get_view() const
    {return m_view;}

    /// Returns a static view of the underlying kokkos static view.
    KOKKOS_FORCEINLINE_FUNCTION
    array_type get_static_view()
    {return m_view;}

    /// Returns a static view of the underlying kokkos static view.
    KOKKOS_FORCEINLINE_FUNCTION
    const array_type get_static_view() const
    {return m_view;}

    template<typename SrcScalar,typename...SrcProps>
    void deep_copy(const PHX::MDField<SrcScalar,SrcProps...>& source)
    {Kokkos::deep_copy(m_view, source.get_static_view());}
    
    void deep_copy(const Scalar source)
    {Kokkos::deep_copy(m_view, source);}

  private:
    template<int R> constexpr size_type rank(ViewSpecialization<R>) const {return traits::rank;}
    constexpr size_type rank(ViewSpecialization<0>) const {return m_view.rank();}

    template<int R> void setFieldData(ViewSpecialization<R>,const PHX::any& a)
    { 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
      TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, fieldTagErrorMsg());
      m_data_set = true;
#endif
      
      // PHX::any object is always the non-const data type.  To correctly
      // cast the any object to the Kokkos::View, need to pull the const
      // off the scalar type if this MDField has a const scalar type.
      typedef PHX::View<typename array_type::non_const_data_type> non_const_view;
      try {
        non_const_view tmp = PHX::any_cast<non_const_view>(a);
        m_view = tmp;
      }
      catch (std::exception& ) {
        std::cout << "\n\nError in compiletime PHX::MDField::setFieldData() in PHX::any_cast. Tried to cast the field \"" 
                  << this->fieldTag().name()  << "\" with the identifier \"" << this->fieldTag().identifier() 
                  << "\" to a type of \"" << Teuchos::demangleName(typeid(non_const_view).name()) 
                  << "\" from a PHX::any object containing a type of \"" 
                  << Teuchos::demangleName(a.type().name()) << "\"." << std::endl;
        throw;
      }
    }

    void setFieldData(ViewSpecialization<0>,const PHX::any& a)
    { 
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ )
      TEUCHOS_TEST_FOR_EXCEPTION(m_tag.is_null(), std::logic_error, fieldTagErrorMsg());
      m_data_set = true;
#endif
      //m_any = a;
      using NonConstDataT = typename std::remove_const<Scalar>::type;
      try {
        if (m_tag->dataLayout().rank() == 1)
          m_view =  PHX::any_cast<Kokkos::View<NonConstDataT*,typename PHX::DevLayout<NonConstDataT>::type,PHX::Device>>(a);
        else if (m_tag->dataLayout().rank() == 2)
          m_view =  PHX::any_cast<Kokkos::View<NonConstDataT**,typename PHX::DevLayout<NonConstDataT>::type,PHX::Device>>(a);
        else if (m_tag->dataLayout().rank() == 3)
          m_view =  PHX::any_cast<Kokkos::View<NonConstDataT***,typename PHX::DevLayout<NonConstDataT>::type,PHX::Device>>(a);
        else if (m_tag->dataLayout().rank() == 4)
          m_view =  PHX::any_cast<Kokkos::View<NonConstDataT****,typename PHX::DevLayout<NonConstDataT>::type,PHX::Device>>(a);
        else if (m_tag->dataLayout().rank() == 5)
          m_view =  PHX::any_cast<Kokkos::View<NonConstDataT*****,typename PHX::DevLayout<NonConstDataT>::type,PHX::Device>>(a);
        else if (m_tag->dataLayout().rank() == 6)
          m_view =  PHX::any_cast<Kokkos::View<NonConstDataT******,typename PHX::DevLayout<NonConstDataT>::type,PHX::Device>>(a);
        else if (m_tag->dataLayout().rank() == 7)
          m_view =  PHX::any_cast<Kokkos::View<NonConstDataT*******,typename PHX::DevLayout<NonConstDataT>::type,PHX::Device>>(a);
        else {
          throw std::runtime_error("ERROR - PHX::MDField::setFieldData (DynRank) - Invalid rank!");
        }
      }
      catch (std::exception& ) { 
        //std::string type_cast_name = Teuchos::demangleName(typeid(non_const_view).name());
        std::string type_cast_name = "???";
        std::cout << "\n\nError in runtime PHX::MDField::setFieldData() in PHX::any_cast. Tried to cast the field \"" 
                  << this->fieldTag().name()  << "\" with the identifier \"" << this->fieldTag().identifier() 
                  << "\" to a type of \"" << type_cast_name
                  << "\" from a PHX::any object containing a type of \"" 
                  << Teuchos::demangleName(a.type().name()) << "\"." << std::endl;
        throw;
      } 
    }
    
    template<int R> void print(ViewSpecialization<R>, std::ostream& os, bool printValues) const
    {
      std::vector<const char*> dim_names;

      // ROGER TODO: fix printing for static dims!!!!
      
      // PHX::PrintDimension<Tag0,array_type> pd0;
      // pd0.addName(dim_names);
      // PHX::PrintDimension<Tag1,array_type> pd1;
      // pd1.addName(dim_names);
      // PHX::PrintDimension<Tag2,array_type> pd2;
      // pd2.addName(dim_names);
      // PHX::PrintDimension<Tag3,array_type> pd3;
      // pd3.addName(dim_names);
      // PHX::PrintDimension<Tag4,array_type> pd4;
      // pd4.addName(dim_names);
      // PHX::PrintDimension<Tag5,array_type> pd5;
      // pd5.addName(dim_names);
      // PHX::PrintDimension<Tag6,array_type> pd6;
      // pd6.addName(dim_names);
      // PHX::PrintDimension<Tag7,array_type> pd7;
      // pd7.addName(dim_names);

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
        os << m_view.extent(i);
      }
      os << "): ";
      
      if (nonnull(m_tag))
        os << *m_tag;
      
      if (printValues)
        os << "Error - MDField no longer supports the \"printValues\" member of the MDField::print() method. Values may be on a device that does not support printing (e.g. GPU).  Please disconstinue the use of this call!" << std::endl;  
    }
    
    void print(ViewSpecialization<0>, std::ostream& os, bool printValues) const
    {
      os << "MDField(";
      for (size_type i=0; i < m_tag->dataLayout().rank(); ++i) {
        if (i > 0)
          os << ",";
        os << m_tag->dataLayout().dimension(i);
      }
      os << "): ";
      
      if (nonnull(m_tag))
        os << *m_tag;
      
      if (printValues)
        os << "Error - MDField no longer supports the \"printValues\" member of the MDField::print() method. Values may be on a device that does not support printing (e.g. GPU).  Please disconstinue the use of this call!" << std::endl;  
    }

#ifdef PHX_DEBUG
    std::string fieldTagErrorMsg() const
    {return "Error - PHX::MDField - No tag has been set!";}

    std::string fieldDataErrorMsg() const
    {return "Error - PHX::MDField - No data has been set!  Please bind memory (call getFieldData()) to MDField!";}
#endif
    
  private:
    Teuchos::RCP<const PHX::FieldTag> m_tag;
    array_type m_view;
#ifdef PHX_DEBUG
    bool m_data_set;
#endif

    template<typename FScalar,typename...FProps> friend class PHX::MDField;
  };

  template<typename DataT,typename...Props>
  std::ostream& operator<<(std::ostream& os, const PHX::MDField<DataT,Props...>& f)
  {
    f.print(os, false);
    return os;
  }

}

#endif
