// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ViewIterator.hpp
    \brief  Iterator allows linear traversal of (part of) a Kokkos View in a manner that is agnostic to its rank.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_ViewIterator_h
#define Intrepid2_ViewIterator_h

#include "Intrepid2_Utils.hpp"

namespace Intrepid2
{
  /** \class  Intrepid2::ViewIterator
      \brief  A helper class that allows iteration over some part of a Kokkos View, while allowing
              the calling code to remain agnostic as to the rank of the view.
   
   Usage is as follows:
   1. Construct ViewIterator from View.
   2. Use the setLocation() method to set the location appropriate to each thread.
   3. Use get() to access the current value, set() to set it, and increment() to move to the next entry.
   4. Stop condition: if only part of the View is traversed on a thread, nextIncrementRank(), or the returned rank value from increment(), may be used to determine when the View locations for a given dimension have been traversed.
   
   For example, if a rank-4 View has dimensions (D1,D2,D3,D4) = (10,15,20,25), and logical threads are assigned to each (D1,D2) pair,
   then each thread would call setLocation({D1,D2,0,0,0,0,0}), and would stop traversal when increment() returned 1, indicating that
   the D2 rank had been incremented, which in turn would imply that all (D3,D4) pairs had been visited.
   
      @see Intrepid2::TensorViewIterator

      \remark This class is under active development; its interface should be understood to be in flux.
              Furthermore, this may eventually move to another Trilinos package.
  */
  template<class ViewType,typename ScalarType>
  class ViewIterator
  {
    ViewType view_;
    Kokkos::Array<int,7> dims_; // 7 is the maximum rank of a Kokkos view
    Kokkos::Array<int,7> index_;
  public:
    //! Constructor
    //! \param [in] view - the View to iterate over
    KOKKOS_INLINE_FUNCTION
    ViewIterator(ViewType view)
    :
    view_(view)
    {
      for (unsigned d=0; d<getFunctorRank(view); d++)
      {
        dims_[d] = view.extent_int(d);
        index_[d] = 0;
      }
      for (unsigned d=getFunctorRank(view); d<7; d++)
      {
        dims_[d] = 1;
        index_[d] = 0;
      }
    }
    
    //! Getter
    //! \return the value at the current location
    KOKKOS_INLINE_FUNCTION
    ScalarType get()
    {
      return view_.access(index_[0],index_[1],index_[2],index_[3],index_[4],index_[5],index_[6]);
    }
    
    //! Setter
    //! \param [in] value - the value to set at the current location
    KOKKOS_INLINE_FUNCTION
    void set(const ScalarType &value)
    {
      view_.access(index_[0],index_[1],index_[2],index_[3],index_[4],index_[5],index_[6]) = value;
    }
    
    //!
    //! \return the leftmost rank ordinal whose index will change on next increment (-1 if next increment will go out of bounds)
    KOKKOS_INLINE_FUNCTION
    int nextIncrementRank()
    {
      const int rank = getFunctorRank(view_);
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
      const int rank = getFunctorRank(view_);
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
      const auto rank = view_.rank();
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
      for (unsigned d=from_rank_number; d<view_.rank(); d++)
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
    
    //! Sets the current location.
    //! \param [in] location - the location as a 7-element array value.
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<int,7> & getLocation()
    {
      return index_;
    }
  };
} // namespace Intrepid2

#endif /* Intrepid2_ViewIterator_h */
