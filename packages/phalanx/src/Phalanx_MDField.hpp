// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_MDFIELD_HPP
#define PHX_MDFIELD_HPP

#include <any>
#include <iostream>
#include <string>
#include <type_traits>
#include "Phalanx_config.hpp"
#include "Teuchos_RCP.hpp"
#include "Kokkos_DynRankView_Fad.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_ExtentTraits.hpp"
#include "Sacado.hpp"
#include "Sacado_mpl_vector.hpp"
#include "Sacado_mpl_for_each.hpp"
#include "Sacado_mpl_push_back.hpp"
#include "Phalanx_DataLayout_DynamicLayout.hpp"
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
  // Now declared in Phalanx_MDField_ExtentTraits.hpp.
  // Saving here for convenient user reference.
  //template<typename I> struct is_extent : std::false_type {};

  // ****************************
  // Devices
  // ****************************
  template<typename T> struct is_device : std::false_type {};
  template<> struct is_device<PHX::Device> : std::true_type {};

  // ****************************
  // Rank count
  // ****************************
  template<typename...Props> struct RankCount;

  template<typename...Props> struct RankCount : std::integral_constant<int,RankCount<void,Props...>::value>
  {
    using vector_type = typename RankCount<void,Props...>::vector_type;
  };

  template<> struct RankCount<void> : std::integral_constant<int,0>
  {
    using vector_type = Sacado::mpl::vector<>;
  };

  template<typename Extent, typename...Props>
  struct RankCount<typename std::enable_if<is_extent<Extent>::value>::type,Extent,Props...>
    : std::integral_constant<int,1+RankCount<void,Props...>::value>
  {
    using vector_type = typename Sacado::mpl::push_back<typename RankCount<Props...>::vector_type,Extent>::type;
  };

  template<typename NonExtent, typename...Props>
  struct RankCount<typename std::enable_if<!is_extent<NonExtent>::value>::type,NonExtent,Props...>
    : std::integral_constant<int,RankCount<void,Props...>::value>
  {
    using vector_type = typename RankCount<void,Props...>::vector_type;
  };

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
  // AnyType - used to avoid memory use in static rank MDFields
  // ****************************
  template<int Rank> struct AnyType;

  // Static rank default
  template<int Rank> struct AnyType
  {
    void set(const std::any& ){}
    std::any get() const {return std::any(nullptr);}
  };

  // Dynamic rank specialization
  template<> struct AnyType<0>
  {
    std::any m_any;
    void set(const std::any& a){m_any = a;}
    std::any get() const {return m_any;}
  };

  // ****************************
  // ReturnType
  // ****************************
  template<typename ViewType>
  struct MDFieldReturnType {
    using return_type = typename ViewType::reference_type;
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

  template<typename TemplateArg, typename... Props>
  struct FieldTraits< typename std::enable_if<!is_extent<TemplateArg>::value &&
                                              !is_layout<TemplateArg>::value &&
                                              !is_device<TemplateArg>::value >::type, TemplateArg, Props...>
  {
    static_assert(is_extent<TemplateArg>::value || is_layout<TemplateArg>::value || is_device<TemplateArg>::value,
                  "\n\nERROR: Invalid template argument on MDField!\nERROR: Please specialize all MDField template arguments to be an extent, a layout or a device!\n");
  };

  template<typename Scalar, typename... Props>
  struct FieldTraits {
    using prop = FieldTraits<void,Props...>;
    static constexpr int rank = RankCount<Props...>::value;
    // using extent_vector = typename RankCount<Props...>::vector_type;
    // This sets defaults if not specified
    using layout = typename std::conditional< !std::is_same<typename prop::layout, void>::value,typename prop::layout, typename PHX::DevLayout<Scalar>::type>::type;
    using device = typename std::conditional< !std::is_same<typename prop::device, void>::value,typename prop::device, PHX::Device>::type;
    using data_type = typename ArrayType<Scalar,rank,layout,device>::data_type;
    using array_type = typename ArrayType<Scalar,rank,layout,device>::array_type;
  };

  // ****************************
  // MDField
  // ****************************

  /*! A multidimensional array with optional compile time rank tags for self documentation.

    This class currently wraps a Kokkos::View as the underlying data
    structure for performance portability. It also carries along a
    field tag with identifier and a data layout for sizing the
    multidimensional array.

    Design Notes:

    - There are essentially two versions of MDField, a runtime and
      compile time version that switched based on class template
      parameters.

    - Tag dispatch is used to switch between the runtime rank and
      compile time rank implementations. The
      ViewSpecialization<traits::rank> is the tag. For the runtime
      rank version, the traits::rank is 0. For the compile time
      version, the rank is a positive integer greater than zero. Since
      c++17, we can now use "if consexpr" to replace some of the tag
      dispatch when appropriate.

    - The private member m_host_data is a pointer to all data that is
      only accessed on host. This data is not marked with device ctors
      and generally creates warnings with cuda compilers if a copy
      constructor for an MDField is called on device. By using a
      pointer, the warnings are removed and these objects are not
      created on device.

    - The runtime version has an extra data member in the m_host_data
      struct. The member is an any object to hold the true static
      type. The FieldManager always allocates all arrays using the
      static rank Kokkos::View. This ensures performance when we can
      use the static view. The runtime versions can always be created
      from the static versions, but tend to be less performant due to
      runtime indexing.

    - We can assign runtime and static MDFields to each other. This
      means that the rank comparison can be a runtime check.
   */
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
#ifdef PHX_DEBUG
    enum { rank_value = traits::rank }; // for printing in debug mode
#endif
    typedef Scalar value_type;
    typedef Scalar& reference_type;

  private:

    struct HostData {
      Teuchos::RCP<const PHX::FieldTag> m_tag;
      PHX::AnyType<traits::rank> m_any_holder; // only used for dynamic rank (empty otherwise)
    };

    /// Kokkos View or DynRankView of the field
    array_type m_view;

    /// Host data. Stored as a raw pointer to eliminate warnings in device code for ctors and dtors of member RCPs and std::any.
    HostData* m_host_data;

#ifdef PHX_DEBUG
    /// Debug flag that only exists if phalanx debug configure flag is set to true.
    bool m_data_set;
#endif

  public:

    /// ONLY USE THIS CTOR FOR UNMANAGED FIELDS!!!! It will allocate memory unassociated with the DAG!
    template<typename...Extents>
    MDField(const std::string name,const std::string layout_name,Extents... e)
      : m_view(name,e...)
#ifdef PHX_DEBUG
      , m_data_set(true)
#endif
    {
      m_host_data = new HostData;
      Teuchos::RCP<PHX::Layout> layout = Teuchos::rcp(new PHX::Layout(layout_name,e...));
      m_host_data->m_tag = Teuchos::rcp(new PHX::Tag<value_type>(name,layout));

      if constexpr (traits::rank==0) {
        m_host_data->m_any_holder.set(std::any(m_view));
      }
    }

    MDField(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& dl)
#ifdef PHX_DEBUG
      : m_data_set(false)
#endif
    {
      m_host_data = new HostData;
      m_host_data->m_tag = Teuchos::rcp(new PHX::Tag<Scalar>(name,dl));
    }

    MDField(const PHX::FieldTag& t)
#ifdef PHX_DEBUG
      : m_data_set(false)
#endif
    {
      m_host_data = new HostData;
      m_host_data->m_tag = t.clone();
    }

    MDField(const Teuchos::RCP<const PHX::FieldTag>& t)
#ifdef PHX_DEBUG
      : m_data_set(false)
#endif
    {
      m_host_data = new HostData;
      m_host_data->m_tag = t;
    }

    /// Default empty constructor
    MDField()
#ifdef PHX_DEBUG
      : m_data_set(false)
#endif
    {m_host_data = new HostData;}

    /// Copy ctor
    KOKKOS_FUNCTION
    MDField(const MDField& source)
      : m_view(source.m_view),
        m_host_data(nullptr)
#ifdef PHX_DEBUG
      , m_data_set(source.m_data_set)
#endif
    {
      // Do not create host data if on device
      KOKKOS_IF_ON_HOST( m_host_data = new HostData; )
      KOKKOS_IF_ON_HOST( m_host_data->m_tag = source.m_host_data->m_tag; )
      KOKKOS_IF_ON_HOST( m_host_data->m_any_holder.set(source.get_static_view_as_any()); )
    }

    /** \brief Templated copy ctor

        The templated version allows for different extra template
        parameters/dimension tags. For example, one evaluator could a
        use "Point" dim tag and be copied to an evaluator that uses a
        "QuadraturePoint" dim tag for the same field in a different
        evaluator. Another example is for assigning a compiletime
        MDFields to runtime MDFields and vice versa. Examples:

        MDField<double,Cell,Point> a("a",data_layout);
        a.setFieldData(memory);
        MDField<double,Cell,QuadPoint> b;
        b = a;

        MDField<double> c; // dynamic rank
        c = a;

        Another example could be for atomic access flags.
    */
    template<typename SourceScalar,typename...SourceProps>
    KOKKOS_FUNCTION
    MDField(const MDField<SourceScalar,SourceProps...>& source)
      : m_view(source.m_view),
        m_host_data(nullptr)
#ifdef PHX_DEBUG
      , m_data_set(source.m_data_set)
#endif
    {
      static_assert(std::is_same<typename std::decay<Scalar>::type, typename std::decay<SourceScalar>::type>::value,
                    "ERROR: Compiletime MDField copy ctor requires scalar types to be the same!");

      // Do not create host data if on device
      KOKKOS_IF_ON_HOST( m_host_data = new HostData; )
      KOKKOS_IF_ON_HOST( m_host_data->m_tag = source.m_host_data->m_tag; )
      KOKKOS_IF_ON_HOST( m_host_data->m_any_holder.set(source.get_static_view_as_any()); )
    }

    KOKKOS_FUNCTION
    ~MDField()
    {
      // Only delete if on host (will not ba allocated on device)
      KOKKOS_IF_ON_HOST( delete m_host_data; )
#ifdef PHX_DEBUG
      // Should not be needed, but helps in debug builds to catch ctor issues
      KOKKOS_IF_ON_HOST( m_host_data = nullptr; )
#endif
    }

    constexpr bool is_static() const {return (traits::rank != 0);}
    constexpr bool is_dynamic() const {return (traits::rank == 0);}

    KOKKOS_INLINE_FUNCTION
    constexpr size_type rank() const {return rank(ViewSpecialization<traits::rank>());}

    KOKKOS_INLINE_FUNCTION
    constexpr size_t size() const {return m_view.size();}

    KOKKOS_INLINE_FUNCTION
    constexpr size_t span() const {return m_view.span();}

    const PHX::FieldTag& fieldTag() const
    {
#if defined(PHX_DEBUG)
      TEUCHOS_TEST_FOR_EXCEPTION(m_host_data == nullptr, std::runtime_error, hostDataErrorMsg());
      TEUCHOS_TEST_FOR_EXCEPTION(m_host_data->m_tag.is_null(), std::logic_error, fieldTagErrorMsg());
#endif
      return *(m_host_data->m_tag);
    }

    Teuchos::RCP<const PHX::FieldTag> fieldTagPtr() const
    {
#if defined(PHX_DEBUG)
      TEUCHOS_TEST_FOR_EXCEPTION(m_host_data->m_tag.is_null(), std::logic_error, fieldTagErrorMsg());
#endif
      return m_host_data->m_tag;
    }

    PHX::MDField<Scalar,Props...>&
    operator=(const MDField<Scalar,Props...>& source)
    {
      m_view = source.m_view;
      m_host_data->m_tag = source.m_host_data->m_tag;
      m_host_data->m_any_holder.set(source.get_static_view_as_any());
#ifdef PHX_DEBUG
      m_data_set = source.m_data_set;
#endif
      return *this;
    }

    template<typename SrcScalar,typename...SrcProps>
    PHX::MDField<Scalar,Props...>&
    operator=(const MDField<SrcScalar,SrcProps...>& source)
    {
      m_view = source.m_view;
      m_host_data->m_tag = source.m_host_data->m_tag;
      m_host_data->m_any_holder.set(source.get_static_view_as_any());
#ifdef PHX_DEBUG
      m_data_set = source.m_data_set;
#endif
      static_assert(std::is_same<typename std::decay<Scalar>::type, typename std::decay<SrcScalar>::type>::value,
                    "ERROR: MDField assignment operator requires scalar types to be the same!");

      return *this;
    }

    template<typename... index_pack>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldReturnType<array_type>::return_type
    operator()(const index_pack&... indices) const
    {
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ ) && !defined(__HIP_DEVICE_COMPILE__)
      TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, fieldDataErrorMsg());
#endif
      return m_view(indices...);
    }

    template<typename... index_pack>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldReturnType<array_type>::return_type
    access(const index_pack&... indices) const
    {
#if defined( PHX_DEBUG) && !defined (__CUDA_ARCH__ ) && !defined(__HIP_DEVICE_COMPILE__)
      TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, fieldDataErrorMsg());
#endif
      return m_view.access(indices...);
    }

    template<typename iType0>
    KOKKOS_FORCEINLINE_FUNCTION
    typename PHX::MDFieldReturnType<array_type>::return_type
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
    {m_host_data->m_tag = t.clone();}

    void setFieldTag(const Teuchos::RCP<const PHX::FieldTag>& t)
    {m_host_data->m_tag = t;}

    void setFieldData(const std::any& a)
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
#if defined( PHX_DEBUG)
      TEUCHOS_TEST_FOR_EXCEPTION(m_host_data->m_tag.is_null(), std::logic_error, fieldTagErrorMsg());
      TEUCHOS_TEST_FOR_EXCEPTION(!m_data_set, std::logic_error, fieldDataErrorMsg());
#endif
      dims.resize(this->rank());
      for (size_type i=0; i <  this->rank(); ++i)
        dims[i] = static_cast<iType>(m_view.extent(i));  // dangerous
    }

    KOKKOS_FORCEINLINE_FUNCTION
    operator array_type () const { return get_static_view();}

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

    void deep_copy(const Scalar& source)
    {Kokkos::deep_copy(m_view, source);}

    std::any get_static_view_as_any() const
    {return get_static_view_as_any(ViewSpecialization<traits::rank>());}

    /// Resets the underlying view ptr to null.
    void releaseFieldData()
    {
#if defined(PHX_DEBUG)
      m_data_set = false;
#endif
      m_view = array_type();
    }

  private:
    template<int R> KOKKOS_INLINE_FUNCTION constexpr size_type rank(ViewSpecialization<R>) const {return traits::rank;}
    KOKKOS_INLINE_FUNCTION constexpr size_type rank(ViewSpecialization<0>) const {return m_view.rank();}

    template<int R> void setFieldData(ViewSpecialization<R>,const std::any& a)
    {
#if defined(PHX_DEBUG)
      TEUCHOS_TEST_FOR_EXCEPTION(m_host_data == nullptr, std::logic_error, hostDataErrorMsg());
      TEUCHOS_TEST_FOR_EXCEPTION(m_host_data->m_tag.is_null(), std::logic_error, fieldTagErrorMsg());
      m_data_set = true;
#endif

      // std::any object is always the non-const data type.  To correctly
      // cast the any object to the Kokkos::View, need to pull the const
      // off the scalar type if this MDField has a const scalar type.
      using non_const_view = PHX::View<typename array_type::non_const_data_type>;
      try {
        non_const_view tmp = std::any_cast<non_const_view>(a);
        m_view = tmp;
      }
      catch (std::exception& ) {
        std::cout << "\n\nError in compiletime PHX::MDField::setFieldData() in std::any_cast. Tried to cast the field \""
                  << this->fieldTag().name()  << "\" with the identifier \"" << this->fieldTag().identifier()
                  << "\" to a type of \"" << Teuchos::demangleName(typeid(non_const_view).name())
                  << "\" from a std::any object containing a type of \""
                  << Teuchos::demangleName(a.type().name()) << "\"." << std::endl;
        throw;
      }
    }

    void setFieldData(ViewSpecialization<0>,const std::any& a)
    {
#if defined(PHX_DEBUG)
      TEUCHOS_TEST_FOR_EXCEPTION(m_host_data == nullptr, std::logic_error, hostDataErrorMsg());
      TEUCHOS_TEST_FOR_EXCEPTION(m_host_data->m_tag.is_null(), std::logic_error, fieldTagErrorMsg());
      m_data_set = true;
#endif

      m_host_data->m_any_holder.set(a);

      using NonConstDataT = typename std::remove_const<Scalar>::type;
      try {
        if (m_host_data->m_tag->dataLayout().rank() == 1)
          m_view =  std::any_cast<Kokkos::View<NonConstDataT*,layout_type,device_type>>(a);
        else if (m_host_data->m_tag->dataLayout().rank() == 2)
          m_view =  std::any_cast<Kokkos::View<NonConstDataT**,layout_type,device_type>>(a);
        else if (m_host_data->m_tag->dataLayout().rank() == 3)
          m_view =  std::any_cast<Kokkos::View<NonConstDataT***,layout_type,device_type>>(a);
        else if (m_host_data->m_tag->dataLayout().rank() == 4)
          m_view =  std::any_cast<Kokkos::View<NonConstDataT****,layout_type,device_type>>(a);
        else if (m_host_data->m_tag->dataLayout().rank() == 5)
          m_view =  std::any_cast<Kokkos::View<NonConstDataT*****,layout_type,device_type>>(a);
        else if (m_host_data->m_tag->dataLayout().rank() == 6)
          m_view =  std::any_cast<Kokkos::View<NonConstDataT******,layout_type,device_type>>(a);
        else if (m_host_data->m_tag->dataLayout().rank() == 7)
          m_view =  std::any_cast<Kokkos::View<NonConstDataT*******,layout_type,device_type>>(a);
        else {
          throw std::runtime_error("ERROR - PHX::MDField::setFieldData (DynRank) - Invalid rank!");
        }
      }
      catch (std::exception& ) {
        //std::string type_cast_name = Teuchos::demangleName(typeid(non_const_view).name());
        std::string type_cast_name = "???";
        std::cout << "\n\nError in runtime PHX::MDField::setFieldData() in std::any_cast. Tried to cast the field \""
                  << this->fieldTag().name()  << "\" with the identifier \"" << this->fieldTag().identifier()
                  << "\" to a type of \"" << type_cast_name
                  << "\" from a std::any object containing a type of \""
                  << Teuchos::demangleName(a.type().name()) << "\"." << std::endl;
        throw;
      }
    }

    // Functor for mpl::for_each for printing dimensions
    struct AddString {
      std::vector<std::string>& names;
      AddString(std::vector<std::string>& extent_names) : names(extent_names) {}
      template<typename Extent> void operator()(Extent ) const
      {names.push_back(PHX::print<Extent>());}
    };

    template<int R> void print(ViewSpecialization<R>, std::ostream& os, bool printValues) const
    {
      std::vector<std::string> dim_names;
      using extent_vector = typename RankCount<Props...>::vector_type;
      Sacado::mpl::for_each_no_kokkos<extent_vector>(AddString(dim_names));

      os << "MDField<";

      int count = 0;
      for (auto i=dim_names.rbegin(); i < dim_names.rend(); ++i,++count) {
        if (count > 0)
          os << ",";
        os << std::string(*i);
      }
      os << ">(";
      for (std::size_t i=0; i < dim_names.size(); ++i) {
        if (i > 0)
          os << ",";
        os << m_view.extent(i);
      }
      os << "): ";

      if (nonnull(m_host_data->m_tag))
        os << *(m_host_data->m_tag);

      if (printValues)
        os << "Error - MDField no longer supports the \"printValues\" member of the MDField::print() method. Values may be on a device that does not support printing (e.g. GPU).  Please discontinue the use of this call!" << std::endl;
    }

    void print(ViewSpecialization<0>, std::ostream& os, bool printValues) const
    {
      os << "MDField(";
      for (size_type i=0; i < m_host_data->m_tag->dataLayout().rank(); ++i) {
        if (i > 0)
          os << ",";
        os << m_host_data->m_tag->dataLayout().dimension(i);
      }
      os << "): ";

      if (nonnull(m_host_data->m_tag))
        os << *(m_host_data->m_tag);

      if (printValues)
        os << "Error - MDField no longer supports the \"printValues\" member of the MDField::print() method. Values may be on a device that does not support printing (e.g. GPU).  Please discontinue the use of this call!" << std::endl;
    }

    template<int R> std::any get_static_view_as_any(ViewSpecialization<R>) const
    {return std::any(m_view);}

    std::any get_static_view_as_any(ViewSpecialization<0>) const
    {return m_host_data->m_any_holder.get();}

#ifdef PHX_DEBUG
    std::string hostDataErrorMsg() const
    {return "Error - PHX::MDField - m_host_data is null!";}

    std::string fieldTagErrorMsg() const
    {return "Error - PHX::MDField - No tag has been set!";}

    std::string fieldDataErrorMsg() const
    {return "Error - PHX::MDField - No data has been set!  Please bind memory (call getFieldData()) to MDField!";}
#endif

  private:

    template<typename FScalar,typename...FProps> friend class PHX::MDField;
  };

  template<typename DataT,typename...Props>
  std::ostream& operator<<(std::ostream& os, const PHX::MDField<DataT,Props...>& f)
  {
    f.print(os, false);
    return os;
  }

  /// \brief free function to allow one to pass in a kokkos view or MDField and get out a view
  template <typename ...Args>
  const auto as_view(const Kokkos::View<Args...> &a) { return a; }
  template <typename ...Args>
  const auto as_view(const Kokkos::DynRankView<Args...> &a) { return a; }
  template <typename ...Args>
  const auto as_view(const PHX::MDField<Args...> &a) { return a.get_static_view(); }

}

#endif
