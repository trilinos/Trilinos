// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_ViewIterator.hpp
    \brief  Iterator allows linear traversal of (part of) a Kokkos View in a manner that is agnostic to its rank.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_ViewIterator_h
#define Intrepid2_ViewIterator_h

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
      for (unsigned d=0; d<view.rank(); d++)
      {
        dims_[d] = view.extent_int(d);
        index_[d] = 0;
      }
      for (unsigned d=view.rank(); d<7; d++)
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
    void set(ScalarType &value)
    {
      view_.access(index_[0],index_[1],index_[2],index_[3],index_[4],index_[5],index_[6]) = value;
    }
    
    //!
    //! \return the leftmost rank ordinal whose index will change on next increment (-1 if next increment will go out of bounds)
    KOKKOS_INLINE_FUNCTION
    int nextIncrementRank()
    {
      const auto rank = view_.rank();
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
      const auto rank = view_.rank();
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
    void reset(int from_rank_number=0)
    {
      for (unsigned d=from_rank_number; d<view_.rank(); d++)
      {
        index_[d] = 0;
      }
    }
    
    //! Sets the current location.
    //! \param [in] location - the location as a 7-element array value.
    KOKKOS_INLINE_FUNCTION
    void setLocation(const Kokkos::Array<int,7> location)
    {
      index_ = location;
    }
  };
} // namespace Intrepid2

#endif /* Intrepid2_ViewIterator_h */
