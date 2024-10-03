// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef Intrepid2_FunctorIterator_h
#define Intrepid2_FunctorIterator_h

/** \file   Intrepid2_FunctorIterator.hpp
    \brief  Defines the Intrepid2::FunctorIterator class, as well as the Intrepid2::functor_returns_ref SFINAE helper.
    \author Created by Nathan V. Roberts.
*/

namespace Intrepid2
{
  //! SFINAE helper to detect whether a functor returns a reference type.
  template<typename FunctorType, typename ScalarType, int rank>
  class functor_returns_ref{};

  //! SFINAE helper to detect whether rank-0 functor returns a reference type.
  template<typename FunctorType, typename ScalarType>
  class functor_returns_ref<FunctorType,ScalarType,0>
  {
    using return_type        = decltype(std::declval<FunctorType>().operator()());
    using return_type_no_ref = typename std::remove_reference<return_type>::type;
  public:
    static constexpr bool value = !std::is_same<return_type, return_type_no_ref>::value;
  };

  //! SFINAE helper to detect whether rank-1 functor returns a reference type.
  template<typename FunctorType, typename ScalarType>
  class functor_returns_ref<FunctorType,ScalarType,1>
  {
    using return_type        = decltype(std::declval<FunctorType>().operator()(0));
    using return_type_no_ref = typename std::remove_reference<return_type>::type;
  public:
    static constexpr bool value = !std::is_same<return_type, return_type_no_ref>::value;
  };

  //! SFINAE helper to detect whether rank-2 functor returns a reference type.
  template<typename FunctorType, typename ScalarType>
  class functor_returns_ref<FunctorType,ScalarType,2>
  {
    using return_type        = decltype(std::declval<FunctorType>().operator()(0,0));
    using return_type_no_ref = typename std::remove_reference<return_type>::type;
  public:
    static constexpr bool value = !std::is_same<return_type, return_type_no_ref>::value;
  };

  //! SFINAE helper to detect whether rank-3 functor returns a reference type.
  template<typename FunctorType, typename ScalarType>
  class functor_returns_ref<FunctorType,ScalarType,3>
  {
    using return_type        = decltype(std::declval<FunctorType>().operator()(0,0,0));
    using return_type_no_ref = typename std::remove_reference<return_type>::type;
  public:
    static constexpr bool value = !std::is_same<return_type, return_type_no_ref>::value;
  };

  //! SFINAE helper to detect whether rank-4 functor returns a reference type.
  template<typename FunctorType, typename ScalarType>
  class functor_returns_ref<FunctorType,ScalarType,4>
  {
    using return_type        = decltype(std::declval<FunctorType>().operator()(0,0,0,0));
    using return_type_no_ref = typename std::remove_reference<return_type>::type;
  public:
    static constexpr bool value = !std::is_same<return_type, return_type_no_ref>::value;
  };

  //! SFINAE helper to detect whether rank-5 functor returns a reference type.
  template<typename FunctorType, typename ScalarType>
  class functor_returns_ref<FunctorType,ScalarType,5>
  {
    using return_type        = decltype(std::declval<FunctorType>().operator()(0,0,0,0,0));
    using return_type_no_ref = typename std::remove_reference<return_type>::type;
  public:
    static constexpr bool value = !std::is_same<return_type, return_type_no_ref>::value;
  };

  //! SFINAE helper to detect whether rank-6 functor returns a reference type.
  template<typename FunctorType, typename ScalarType>
  class functor_returns_ref<FunctorType,ScalarType,6>
  {
    using return_type        = decltype(std::declval<FunctorType>().operator()(0,0,0,0,0,0));
    using return_type_no_ref = typename std::remove_reference<return_type>::type;
  public:
    static constexpr bool value = !std::is_same<return_type, return_type_no_ref>::value;
  };

  //! SFINAE helper to detect whether rank-7 functor returns a reference type.
  template<typename FunctorType, typename ScalarType>
  class functor_returns_ref<FunctorType,ScalarType,7>
  {
    using return_type        = decltype(std::declval<FunctorType>().operator()(0,0,0,0,0,0,0));
    using return_type_no_ref = typename std::remove_reference<return_type>::type;
  public:
    static constexpr bool value = !std::is_same<return_type, return_type_no_ref>::value;
  };

  //! essentially, a read-only variant of ViewIterator, for a general functor (extent_int() and rank() support required)
  template<class FunctorType, typename ScalarType, int rank>
  class FunctorIterator
  {
    const FunctorType &functor_;
    Kokkos::Array<int,7> dims_; // 7 is the maximum rank of a Kokkos view
    Kokkos::Array<int,7> index_;
  public:
    //! Constructor.  A reference to the functor is stored.  This means that FunctorIterators should be constructed where they will be used — on device, or on host, e.g., but not copied from host to device.
    //! \param [in] functor - the functor to iterate over
    KOKKOS_INLINE_FUNCTION
    FunctorIterator(const FunctorType &functor)
    :
    functor_(functor)
    {
      for (int d=0; d<rank; d++)
      {
        dims_[d] = functor.extent_int(d);
        index_[d] = 0;
      }
      for (int d=rank; d<7; d++)
      {
        dims_[d] = 1;
        index_[d] = 0;
      }
    }
    
    using return_type = typename std::conditional< functor_returns_ref<FunctorType, ScalarType, rank>::value, const ScalarType &, const ScalarType>::type;
    
    template< bool B, class T = return_type >
    using enable_if_t = typename std::enable_if<B,T>::type;
    
    //! Getter
    //! \return the value at the current location
    template<int M = rank>
    enable_if_t<M==0>
    KOKKOS_INLINE_FUNCTION
    get() const
    {
      return functor_();
    }
    
    //! Getter
    //! \return the value at the current location
    template<int M = rank>
    enable_if_t<M==1>
    KOKKOS_INLINE_FUNCTION
    get() const
    {
      return functor_(index_[0]);
    }
    
    //! Getter
    //! \return the value at the current location
    template<int M = rank>
    enable_if_t<M==2>
    KOKKOS_INLINE_FUNCTION
    get() const
    {
      return functor_(index_[0], index_[1]);
    }
    
    //! Getter
    //! \return the value at the current location
    template<int M = rank>
    enable_if_t<M==3>
    KOKKOS_INLINE_FUNCTION
    get() const
    {
      return functor_(index_[0], index_[1], index_[2]);
    }
    
    //! Getter
    //! \return the value at the current location
    template<int M = rank>
    enable_if_t<M==4>
    KOKKOS_INLINE_FUNCTION
    get() const
    {
      return functor_(index_[0], index_[1], index_[2], index_[3]);
    }
    
    //! Getter
    //! \return the value at the current location
    template<int M = rank>
    enable_if_t<M==5>
    KOKKOS_INLINE_FUNCTION
    get() const
    {
      return functor_(index_[0], index_[1], index_[2], index_[3], index_[4]);
    }
    
    //! Getter
    //! \return the value at the current location
    
    template<int M = rank>
    enable_if_t<M==6>
    KOKKOS_INLINE_FUNCTION
    get() const
    {
      return functor_(index_[0], index_[1], index_[2], index_[3], index_[4], index_[5]);
    }
    
    //! Getter
    //! \return the value at the current location
    template<int M = rank>
    enable_if_t<M==7>
    KOKKOS_INLINE_FUNCTION
    get() const
    {
      return functor_(index_[0], index_[1], index_[2], index_[3], index_[4], index_[5], index_[6]);
    }
    
    //!
    //! \return the leftmost rank ordinal whose index will change on next increment (-1 if next increment will go out of bounds)
    KOKKOS_INLINE_FUNCTION
    int nextIncrementRank()
    {
      for (int r=rank-1; r>=0; r--)
      {
        if (index_[r]+1 < dims_[r]) // can increment without going out of bounds in this dimension
        {
          return r;
        }
      }
      // next increment will take us out of bounds
      return -1;
    }
    
    //! move to the next location
    //! \return the rank of the leftmost index that was changed; -1 if increment reached the end of the view
    KOKKOS_INLINE_FUNCTION
    int increment()
    {
      for (int r=rank-1; r>=0; r--)
      {
        if (index_[r]+1 < dims_[r]) // can increment without going out of bounds in this dimension
        {
          index_[r]++;
          // we've completed the increment
          return r;
        }
        else
        {
          // next rank should be incremented -- this one should reset to 0
          index_[r] = 0;
        }
      }
      // if we get here, we have run through all ranks, setting them to 0 -- we've cycled around
      // and in that sense have not completed the increment
      return -1;
    }
    
    //! move to the previous location
    //! \return the rank of the leftmost index that was changed
    KOKKOS_INLINE_FUNCTION
    bool decrement()
    {
      for (int r=rank-1; r>=0; r--)
      {
        if (index_[r]-1 >= 0) // can decrement without going out of bounds in this dimension
        {
          index_[r]--;
          return true; // we've completed the decrement
        }
        else
        {
          // next rank should be decremented -- this one should cycle round to dim_[r]-1
          index_[r] = dims_[r]-1;
        }
      }
      // if we get here, we've gone past 0 in every dimension, so we should return false
      // -- we have not completed the decrement in an in-bounds fashion, but have cycled round to the last value
      // to maintain a clean state, let's reset
      reset();
      return false;
    }
    
    //! Enumeration index refers to a 1D enumeration of the entries in the View, with dimensions in order of their significance (dimension 0 is the slowest-moving).
    //! \return the enumeration index at current location.
    KOKKOS_INLINE_FUNCTION
    int getEnumerationIndex()
    {
      int index_1D = 0;
      for (int d=0; d<7; d++)
      {
        if (d>0) index_1D *= dims_[d-1];
        index_1D += index_[d];
      }
      
      return index_1D;
    }
    
    //! The enumeration index refers to a 1D enumeration of the entries in the View, with dimensions in order of their significance (dimension 0 is the slowest-moving).
    //! \param [in] enumerationIndex - the index to which the location should be set
    KOKKOS_INLINE_FUNCTION
    void setEnumerationIndex(const int &enumerationIndex)
    {
      Kokkos::Array<int,7> location;
      int remainder = enumerationIndex;
      for (int d=6; d>=0; d--)
      {
        location[d] = remainder % dims_[d];
        remainder  /= dims_[d];
      }
      
      setLocation(location);
    }
    
    //! The index of the current location in the specified dimension.  (Indices in dimensions beyond the rank of the View, but less than 7, are defined to be 0.)
    //! \param [in] dimension - the dimension for which the current index should be returned.
    //! \return index in the specified dimension.
    KOKKOS_INLINE_FUNCTION
    int getIndex(int dimension)
    {
      return index_[dimension];
    }
    
    //! The extent of the View in the specified dimension.  (Extents in dimensions beyond the rank of the View, but less than 7, are defined to be 1.)
    //! \param [in] dimension - the dimension for which the extent should be returned.
    //! \return extent of the View in the specified dimension.
    KOKKOS_INLINE_FUNCTION
    int getExtent(int dimension)
    {
      return dims_[dimension];
    }
    
    //! Resets the location to index 0 in each dimension, starting from the specified dimension.
    //! \param [in] from_rank_number - the first dimension in which to set the index to 0.
    KOKKOS_INLINE_FUNCTION
    void reset(unsigned from_rank_number=0)
    {
      for (unsigned d=from_rank_number; d<functor_.rank(); d++)
      {
        index_[d] = 0;
      }
    }
    
    //! Sets the current location.
    //! \param [in] location - the location as a 7-element array value.
    KOKKOS_INLINE_FUNCTION
    void setLocation(const Kokkos::Array<int,7> &location)
    {
      index_ = location;
    }
    
    //! Sets the current location in the specified dimension.
    //! \param [in] dim - which dimension to set the location in
    //! \param [in] i - the index to use in dimension <b>dim</b>
    KOKKOS_INLINE_FUNCTION
    void setLocationInDim(const int &dim, const int &i)
    {
      index_[dim] = i;
    }
    
    //! Sets the current location.
    //! \param [in] location - the location as a 7-element array value.
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<int,7> & getLocation()
    {
      return index_;
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_FunctorIterator_h */
