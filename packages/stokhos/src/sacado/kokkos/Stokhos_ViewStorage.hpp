// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_VIEW_STORAGE_HPP
#define STOKHOS_VIEW_STORAGE_HPP

#include "Stokhos_DynArrayTraits.hpp"

#include "Kokkos_Macros.hpp"
#include "impl/Kokkos_Traits.hpp"

#include "Sacado_Traits.hpp"
#include "Stokhos_KokkosTraits.hpp"
#include <sstream>

namespace Stokhos {

  struct error_storage_type_is_not_allocateable {};
  struct error_storage_type_is_not_resizeable {};

  //! Dynamic storage with view semantics and contiguous access
  template < typename ordinal_t ,
             typename value_t ,
             unsigned static_length ,
             unsigned static_stride ,
             typename device_t >
  class ViewStorage {
  private:
    // Enumerated flag so logic is evaluated at compile-time
    enum { stride_one = 1 == static_stride };
  public:

    static const bool is_static       = static_length != 0 ;
    static const int  static_size     = static_length ;
    static const bool supports_reset  = false ;
    static const bool supports_resize = ! is_static ;
    static const bool supports_view   = true ;

    typedef ordinal_t          ordinal_type;
    typedef value_t            value_type;
    typedef device_t           execution_space;
    typedef value_type       & reference;
    typedef const value_type & const_reference;
    typedef value_type       * pointer;
    typedef const value_type * const_pointer;
    // typedef Stokhos::DynArrayTraits<value_type,execution_space> ds;

    //! Turn ViewStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t = value_t , typename dev_t = device_t >
    struct apply {
      typedef ViewStorage<ord_t,val_t,static_length,static_stride,dev_t> type;
    };

    //! Constructor to satisfy Sacado::MP::Vector, disabled via error type.
    KOKKOS_INLINE_FUNCTION
    ViewStorage( const error_storage_type_is_not_allocateable & z =
                   error_storage_type_is_not_allocateable(),
                 const value_type& x = value_type(0) );

    //! Constructor
    KOKKOS_INLINE_FUNCTION
    ViewStorage( pointer v ,
                 const ordinal_type & arg_size = 0 ,
                 const ordinal_type & arg_stride = 0 ) :
      coeff_(v), size_(arg_size), stride_(arg_stride) {}

    //! Constructor from array
    KOKKOS_INLINE_FUNCTION
    ViewStorage(const ordinal_type& sz, const value_type* x) {}

    //! Constructor for creating a view
    KOKKOS_INLINE_FUNCTION
    ViewStorage(const ordinal_type& sz, pointer v, bool owned) :
      coeff_(v), size_(sz) {}

    //! Constructor
    KOKKOS_INLINE_FUNCTION
    ViewStorage(const ViewStorage& s) :
      coeff_(s.coeff_), size_(s.size_), stride_(s.stride_) {}

    KOKKOS_INLINE_FUNCTION
    ViewStorage( const ViewStorage& s , const ordinal_type & arg_begin , const ordinal_type & arg_end )
      : coeff_( ( ordinal_type(size_.value) < arg_end ) ? pointer(0) : s.coeff_ + arg_begin * ordinal_type(s.stride_.value) )
      , size_(  ( ordinal_type(size_.value) < arg_end ) ? ordinal_type(0) : arg_end - arg_begin )
      , stride_( s.stride_ )
      {}

    //! Destructor
    KOKKOS_INLINE_FUNCTION
    ~ViewStorage() {}

private:
    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    ViewStorage& operator=(const ViewStorage& );
public:

    //! Initialize values to a constant value
    KOKKOS_INLINE_FUNCTION
    void init(const_reference v) {
      if ( stride_one ) {
        for ( ordinal_type i = 0 ; i < size_.value ; ++i ) { coeff_[i] = v ; }
      }
      else {
        for ( ordinal_type i = 0 ; i < size_.value * stride_.value ; i += stride_.value ) {
          coeff_[i] = v ;
        }
      }
    }

    //! Initialize values to an array of values
    KOKKOS_INLINE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      const ordinal_type n = sz ? sz : size_.value ;
      if ( stride_one ) {
        for ( ordinal_type i = 0 ; i < n ; ++i ) { coeff_[i] = v[i] ; }
      }
      else {
        for ( ordinal_type i = 0 ; i < n ; ++i ) { coeff_[i*stride_.value] = v[i] ; }
      }
    }

    //! Load values to an array of values
    KOKKOS_INLINE_FUNCTION
    void load(pointer v) {
      if ( stride_one ) {
        for ( ordinal_type i = 0 ; i < size_.value ; ++i ) { coeff_[i] = v[i] ; }
      }
      else {
        for ( ordinal_type i = 0 ; i < size_.value ; ++i ) { coeff_[i*stride_.value] = v[i] ; }
      }
    }

    //! Resize function disabled.
    KOKKOS_INLINE_FUNCTION
    void resize( ordinal_type s ) {}
    //void resize( const error_storage_type_is_not_resizeable & );

    //! Reset function disabled.
    KOKKOS_INLINE_FUNCTION
    void shallowReset( pointer v ,
                       const error_storage_type_is_not_resizeable & ,
                       const error_storage_type_is_not_resizeable & ,
                       const bool owned );

    //! Return size
    KOKKOS_INLINE_FUNCTION
    ordinal_type size() const { return size_.value ; }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    reference operator[] (const ordinal_type& i) const
    {
      return coeff_[ stride_one ? i : i * stride_.value ];
    }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    reference getCoeff() { return coeff_[ stride_one ? i : i * stride_.value ]; }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    const_reference getCoeff() const { return coeff_[ stride_one ? i : i * stride_.value ]; }

    //! Get coefficients
    KOKKOS_INLINE_FUNCTION
    const_pointer coeff() const { return coeff_; }

    //! Get coefficients
    KOKKOS_INLINE_FUNCTION
    pointer coeff() { return coeff_; }

  private:

    //! Coefficient values
    const pointer coeff_;

    //! Size of array
    const Kokkos::Impl::integral_nonzero_constant< ordinal_t , static_length > size_ ;

    //! Stride of array
    const Kokkos::Impl::integral_nonzero_constant< ordinal_t , static_stride > stride_ ;
  };

  template< class Storage >
  struct is_ViewStorage { enum { value = false }; };

  template < typename ordinal_t ,
             typename value_t ,
             unsigned static_length ,
             unsigned static_stride ,
             typename device_t >
  struct is_ViewStorage< ViewStorage< ordinal_t , value_t , static_length , static_stride , device_t > >
    { enum { value = true }; };


} /* namespace Stokhos */

namespace Sacado {
  template <typename ordinal_t, typename value_t, unsigned static_length,
            unsigned static_stride, typename device_t>
  struct StringName< Stokhos::ViewStorage<ordinal_t,
                                          value_t,
                                          static_length,
                                          static_stride,
                                          device_t> > {
    static std::string eval() {
      std::stringstream ss;
      ss << "Stokhos::ViewStorage<"
         << StringName<ordinal_t>::eval() << ","
         << StringName<value_t>::eval() << ","
         << static_length << ","
         << static_stride << ","
         << StringName<device_t>::eval() << ">";
      return ss.str();
    }
  };
}

#endif /* #ifndef STOKHOS_VIEW_STORAGE_HPP */
