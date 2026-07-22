// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if defined( __CUDA_ARCH__ )

namespace Stokhos {

  template <typename ordinal_t, typename value_t>
  class DynamicThreadedStorage<ordinal_t, value_t, Kokkos::Cuda> {
  public:

    static const bool is_static = false;
    static const int static_size = 0;
    static const bool supports_reset = true;

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef Kokkos::Cuda execution_space;
    typedef value_type& reference;
    typedef volatile value_type& volatile_reference;
    typedef const value_type& const_reference;
    typedef const volatile value_type& const_volatile_reference;
    typedef value_type* pointer;
    typedef volatile value_type* volatile_pointer;
    typedef const value_type* const_pointer;
    typedef const volatile value_type* const_volatile_pointer;
    typedef Stokhos::DynArrayTraits<value_type,execution_space> ds;

    //! Turn DynamicThreadedStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t = value_t , typename dev_t = Kokkos::Cuda >
    struct apply {
      typedef DynamicThreadedStorage<ord_t,val_t,dev_t> type;
    };

    //! Constructor
    __device__
    DynamicThreadedStorage(const ordinal_type& sz = 1,
                           const value_type& x = value_type(0.0)) :
      sz_(sz), stride_(num_threads()), total_sz_(sz_*stride_) {
      allocate_coeff_array(coeff_, is_owned_, total_sz_, x);
    }

    //! Constructor from array
    __device__
    DynamicThreadedStorage(const ordinal_type& sz, const value_type* x) :
      sz_(sz), stride_(num_threads()), total_sz_(sz_*stride_) {
      allocate_coeff_array(coeff_, is_owned_, total_sz_, x);
    }

    //! Constructor for creating a view
    __device__
    DynamicThreadedStorage(const ordinal_type& sz, pointer v, bool owned) :
      coeff_(v), sz_(sz), stride_(num_threads()), total_sz_(sz_*stride_),
      is_owned_(owned) {}

    //! Constructor
    __device__
    DynamicThreadedStorage(const DynamicThreadedStorage& s) :
    sz_(s.sz_), stride_(s.stride_), total_sz_(s.total_sz_) {
      allocate_coeff_array(coeff_, is_owned_, total_sz_);
      for (ordinal_type i=0; i<total_sz_; i+=stride_)
        coeff_[i] = s.coeff_[i];
    }

    //! Constructor
    __device__
    DynamicThreadedStorage(const volatile DynamicThreadedStorage& s) :
    sz_(s.sz_), stride_(s.stride_), total_sz_(s.total_sz_) {
      allocate_coeff_array(coeff_, is_owned_, total_sz_);
      for (ordinal_type i=0; i<total_sz_; i+=stride_)
        coeff_[i] = s.coeff_[i];
    }

    //! Destructor
    __device__
    ~DynamicThreadedStorage() {
      destroy_coeff_array(coeff_, is_owned_, total_sz_);
    }

    //! Assignment operator
    __device__
    DynamicThreadedStorage& operator=(const DynamicThreadedStorage& s) {
      if (&s != this) {
        if (s.sz_ != sz_) {
          destroy_coeff_array(coeff_, is_owned_, total_sz_);
          sz_ = s.sz_;
          stride_ = s.stride_;
          total_sz_ = sz_*stride_;
          allocate_coeff_array(coeff_, is_owned_, total_sz_);
          for (ordinal_type i=0; i<total_sz_; i+=stride_)
            coeff_[i] = s.coeff_[i];
        }
        else {
          for (ordinal_type i=0; i<total_sz_; i+=stride_)
            coeff_[i] = s.coeff_[i];
        }
      }
      return *this;
    }

    //! Assignment operator
    __device__
    DynamicThreadedStorage&
    operator=(const volatile DynamicThreadedStorage& s) {
      if (&s != this) {
        if (s.sz_ != sz_) {
          destroy_coeff_array(coeff_, is_owned_, total_sz_);
          sz_ = s.sz_;
          stride_ = s.stride_;
          total_sz_ = sz_*stride_;
          allocate_coeff_array(coeff_, is_owned_, total_sz_);
          for (ordinal_type i=0; i<total_sz_; i+=stride_)
            coeff_[i] = s.coeff_[i];
        }
        else {
          for (ordinal_type i=0; i<total_sz_; i+=stride_)
            coeff_[i] = s.coeff_[i];
        }
      }
      return *this;
    }

    //! Assignment operator
    __device__
    volatile DynamicThreadedStorage&
    operator=(const DynamicThreadedStorage& s) volatile {
      if (&s != this) {
        if (s.sz_ != sz_) {
          destroy_coeff_array(coeff_, is_owned_, total_sz_);
          sz_ = s.sz_;
          stride_ = s.stride_;
          total_sz_ = sz_*stride_;
          allocate_coeff_array(coeff_, is_owned_, total_sz_);
          for (ordinal_type i=0; i<total_sz_; i+=stride_)
            coeff_[i] = s.coeff_[i];
        }
        else {
          for (ordinal_type i=0; i<total_sz_; i+=stride_)
            coeff_[i] = s.coeff_[i];
        }
      }
      return *this;
    }

    //! Assignment operator
    __device__
    volatile DynamicThreadedStorage&
    operator=(const volatile DynamicThreadedStorage& s) volatile {
      if (&s != this) {
        if (s.sz_ != sz_) {
          destroy_coeff_array(coeff_, is_owned_, total_sz_);
          sz_ = s.sz_;
          stride_ = s.stride_;
          total_sz_ = sz_*stride_;
          allocate_coeff_array(coeff_, is_owned_, total_sz_);
          for (ordinal_type i=0; i<total_sz_; i+=stride_)
            coeff_[i] = s.coeff_[i];
        }
        else {
          for (ordinal_type i=0; i<total_sz_; i+=stride_)
            coeff_[i] = s.coeff_[i];
        }
      }
      return *this;
    }

    //! Initialize values to a constant value
    __device__
    void init(const_reference v) {
      for (ordinal_type i=0; i<total_sz_; i+=stride_)
        coeff_[i] = v;
    }

    //! Initialize values to a constant value
    __device__
    void init(const_reference v) volatile {
      for (ordinal_type i=0; i<total_sz_; i+=stride_)
        coeff_[i] = v;
    }

    //! Initialize values to an array of values
    __device__
    void init(const_pointer v, const ordinal_type& sz = 0) {
      ordinal_type my_sz = stride_*sz;
      if (sz == 0)
        my_sz = total_sz_;
      for (ordinal_type i=0; i<my_sz; i+=stride_)
        coeff_[i] = v[i];
    }

    //! Initialize values to an array of values
    __device__
    void init(const_pointer v, const ordinal_type& sz = 0) volatile {
      ordinal_type my_sz = stride_*sz;
      if (sz == 0)
        my_sz = total_sz_;
      for (ordinal_type i=0; i<my_sz; i+=stride_)
        coeff_[i] = v[i];
    }

    //! Load values to an array of values
    __device__
    void load(pointer v) {
      for (ordinal_type i=0; i<total_sz_; i+=stride_)
        coeff_[i] = v[i];
    }

    //! Load values to an array of values
    __device__
    void load(pointer v) volatile {
      for (ordinal_type i=0; i<total_sz_; i+=stride_)
        coeff_[i] = v[i];
    }

    //! Resize to new size (values are preserved)
    __device__
    void resize(const ordinal_type& sz) {
      if (sz != sz_) {
        value_type *coeff_new;
        bool owned_new;
        ordinal_type total_sz_new = sz*stride_;
        allocate_coeff_array(coeff_new, owned_new, total_sz_new);
        ordinal_type my_tsz = total_sz_;
        if (total_sz_ > total_sz_new)
          my_tsz = total_sz_new;
        for (ordinal_type i=0; i<my_tsz; i+=stride_)
          coeff_new[i] = coeff_[i];
        destroy_coeff_array(coeff_, is_owned_, total_sz_);
        coeff_ = coeff_new;
        sz_ = sz;
        total_sz_ = total_sz_new;
        is_owned_ = owned_new;
      }
    }

    //! Resize to new size (values are preserved)
    __device__
    void resize(const ordinal_type& sz) volatile {
      if (sz != sz_) {
        value_type *coeff_new;
        bool owned_new;
        ordinal_type total_sz_new = sz*stride_;
        allocate_coeff_array(coeff_new, owned_new, total_sz_new);
        ordinal_type my_tsz = total_sz_;
        if (total_sz_ > total_sz_new)
          my_tsz = total_sz_new;
        for (ordinal_type i=0; i<my_tsz; i+=stride_)
          coeff_new[i] = coeff_[i];
        destroy_coeff_array(coeff_, is_owned_, total_sz_);
        coeff_ = coeff_new;
        sz_ = sz;
        total_sz_ = total_sz_new;
        is_owned_ = owned_new;
      }
    }

    //! Reset storage to given array, size, and stride
    __device__
    void shallowReset(pointer v, const ordinal_type& sz,
                      const ordinal_type& stride, bool owned) {
      destroy_coeff_array(coeff_, is_owned_, total_sz_);
      coeff_ = v;
      sz_ = sz;
      stride_ = stride;
      total_sz_ = sz_*stride_;
      is_owned_ = owned;
    }

    //! Reset storage to given array, size, and stride
    __device__
    void shallowReset(pointer v, const ordinal_type& sz,
                      const ordinal_type& stride, bool owned) volatile {
      destroy_coeff_array(coeff_, is_owned_, total_sz_);
      coeff_ = v;
      sz_ = sz;
      stride_ = stride;
      total_sz_ = sz_*stride_;
      is_owned_ = owned;
    }

    //! Return size
    __device__
    ordinal_type size() const { return sz_; }

    //! Return size
    __device__
    ordinal_type size() const volatile { return sz_; }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    const_reference operator[] (const ordinal_type& i) const {
      return coeff_[i*stride_];
    }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    const_volatile_reference operator[] (const ordinal_type& i) const volatile {
      return coeff_[i*stride_];
    }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    reference operator[] (const ordinal_type& i) {
      return coeff_[i*stride_];
    }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    volatile_reference operator[] (const ordinal_type& i) volatile {
      return coeff_[i*stride_];
    }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    reference getCoeff() { return coeff_[i*stride_]; }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    volatile_reference getCoeff() volatile { return coeff_[i*stride_]; }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    const_reference getCoeff() const { return coeff_[i*stride_]; }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    const_volatile_reference getCoeff() const volatile { return coeff_[i*stride_]; }

    //! Get coefficients
    KOKKOS_INLINE_FUNCTION
    const_volatile_pointer coeff() const volatile { return coeff_; }

    //! Get coefficients
    KOKKOS_INLINE_FUNCTION
    const_pointer coeff() const { return coeff_; }

    //! Get coefficients
    KOKKOS_INLINE_FUNCTION
    volatile_pointer coeff() volatile { return coeff_; }

    //! Get coefficients
    KOKKOS_INLINE_FUNCTION
    pointer coeff() { return coeff_; }

  protected:

    //! Compute number of threads in each block
    __device__
    ordinal_type num_threads() const {
      return blockDim.x*blockDim.y*blockDim.z;
    }

    //! Compute number of threads in each block
    __device__
    ordinal_type num_threads() const volatile {
      return blockDim.x*blockDim.y*blockDim.z;
    }

    //! Compute thread index within a block
    __device__
    ordinal_type thread_index() const {
      return threadIdx.x + (threadIdx.y + threadIdx.z*blockDim.y)*blockDim.x;
    }

    //! Compute thread index within a block
    __device__
    ordinal_type thread_index() const volatile {
      return threadIdx.x + (threadIdx.y + threadIdx.z*blockDim.y)*blockDim.x;
    }

    //! Allocate coefficient array
    __device__
    void allocate_coeff_array(pointer& c, bool& owned,
                              ordinal_type total_size,
                              const value_type& x = value_type(0.0)) {

      // Allocate coefficient array on thread 0
      __shared__ pointer ptr;
      ordinal_type tidx = thread_index();
      if (tidx == 0) {
        ptr = ds::get_and_fill(total_size,x);
        owned = true;
      }
      else
        owned = false;
      __syncthreads();

      // Give each thread its portion of the array
      c = ptr + tidx;
    }

    //! Allocate coefficient array
    __device__
    void allocate_coeff_array(pointer& c, bool& owned,
                              ordinal_type total_size,
                              const value_type& x = value_type(0.0)) volatile {

      // Allocate coefficient array on thread 0
      __shared__ pointer ptr;
      ordinal_type tidx = thread_index();
      if (tidx == 0) {
        ptr = ds::get_and_fill(total_size,x);
        owned = true;
      }
      else
        owned = false;
      __syncthreads();

      // Give each thread its portion of the array
      c = ptr + tidx;
    }

    //! Allocate coefficient array
    __device__
    void allocate_coeff_array(pointer& c, bool& owned,
                              ordinal_type total_size,
                              const value_type* x) {

      // Allocate coefficient array on thread 0
      __shared__ pointer ptr;
      ordinal_type tidx = thread_index();
      if (tidx == 0) {
        ptr = ds::get_and_fill(x, total_size);
        owned = true;
      }
      else
        owned = false;
      __syncthreads();

      // Give each thread its portion of the array
      c = ptr + tidx;
    }

    //! Destroy coefficient array
    __device__
    void destroy_coeff_array(pointer c, bool owned, ordinal_type total_size) {
      __syncthreads();
      if (owned)
        ds::destroy_and_release(c, total_size);
    }

    //! Destroy coefficient array
    __device__
    void destroy_coeff_array(pointer c, bool owned, ordinal_type total_size) volatile {
      __syncthreads();
      if (owned)
        ds::destroy_and_release(c, total_size);
    }

  private:

    //! Coefficient values
    pointer coeff_;

    //! Size of array used
    ordinal_type sz_;

    //! Stride of array
    ordinal_type stride_;

    //! Total size of array
    ordinal_type total_sz_;

    //! Do we own the array
    bool is_owned_;

  };

}

#endif
