// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_TensorArgumentIterator.hpp
    \brief  Defines TensorArgumentIterator, which allows systematic enumeration of a TensorData object.

    \author Nathan V. Roberts
*/
#ifndef Intrepid2_TensorArgumentIterator_h
#define Intrepid2_TensorArgumentIterator_h

#include "Intrepid2_TensorData.hpp"
#include "Intrepid2_Types.hpp"

namespace Intrepid2
{
/** \class Intrepid2::TensorArgumentIterator
    \brief Allows systematic enumeration of all entries in a TensorData object, tracking indices for each tensor component.
 
   For example, you might have a 3D basis that is the product of 3 1D bases.  The components might have shape (F_d,P_d), with N_F_d and N_P_d as the extents for each.  If we think of the tensor container as having shape (F,P), then F = N_F_0 * N_F_1 * N_F_2, and P = N_P_0 * N_P_1 * N_P_2.  TensorArgumentIterator lets you say that you want to iterate through all the point indices.  When you want to move to the next point, call increment(); this will return the index of the most significant component whose argument changed.  (Component indices that are greater than this one may also have had their argument changed; lower-indexed components will have stayed the same.)
*/
  class TensorArgumentIterator {
    Kokkos::Array<ordinal_type,Parameters::MaxTensorComponents> arguments_;
    Kokkos::Array<ordinal_type,Parameters::MaxTensorComponents> bounds_;
    ordinal_type numTensorComponents_;
  public:
    template<class Scalar, typename ExecSpaceType>
    KOKKOS_INLINE_FUNCTION
    TensorArgumentIterator(const TensorData<Scalar,ExecSpaceType> &tensorData, const ordinal_type argumentOrdinal)
    :
    numTensorComponents_(tensorData.numTensorComponents())
    {
      for (ordinal_type r=0; r<numTensorComponents_; r++)
      {
        arguments_[r] = 0;
        bounds_[r]    = tensorData.getTensorComponent(r).extent_int(argumentOrdinal);
      }
    }
    
    //! Variant that allows truncation of the tensor components at the specified number of components.
    template<class Scalar, typename ExecSpaceType>
    KOKKOS_INLINE_FUNCTION
    TensorArgumentIterator(const TensorData<Scalar,ExecSpaceType> &tensorData, const ordinal_type argumentOrdinal, const ordinal_type numTensorComponents)
    :
    numTensorComponents_(numTensorComponents)
    {
      for (ordinal_type r=0; r<numTensorComponents_; r++)
      {
        arguments_[r] = 0;
        bounds_[r]    = tensorData.getTensorComponent(r).extent_int(argumentOrdinal);
      }
    }
    
    //! Basic constructor in which only the bounds of the tensor components are required.
    TensorArgumentIterator(const std::vector<ordinal_type> tensorComponentBounds)
    :
    numTensorComponents_(tensorComponentBounds.size())
    {
      for (ordinal_type r=0; r<numTensorComponents_; r++)
      {
        arguments_[r] = 0;
        bounds_[r]    = tensorComponentBounds[r];
      }
    }
    
    //! Basic constructor in which only the bounds of the tensor components are required.
    template<size_t rank>
    KOKKOS_INLINE_FUNCTION
    TensorArgumentIterator(const Kokkos::Array<ordinal_type,rank> &tensorComponentBounds)
    :
    numTensorComponents_(rank)
    {
      for (ordinal_type r=0; r<rank; r++)
      {
        arguments_[r] = 0;
        bounds_[r]    = tensorComponentBounds[r];
      }
    }
    
    //! Proceed to next entry.
    KOKKOS_INLINE_FUNCTION ordinal_type increment()
    {
      ordinal_type r = numTensorComponents_ - 1;
      while (arguments_[r] + 1 >= bounds_[r])
      {
        arguments_[r] = 0; // reset
        r--;
        if (r < 0) break;
      }
      if (r >= 0) ++arguments_[r];
      return r;
    }
    
    //!
    //! \return the least tensor component ordinal whose index will change on next increment (-1 if next increment will go out of bounds)
    KOKKOS_INLINE_FUNCTION
    ordinal_type nextIncrementResult() const
    {
      ordinal_type r = numTensorComponents_ - 1;
      while (arguments_[r] + 1 >= bounds_[r])
      {
        r--;
        if (r < 0) break;
      }
      return r;
    }
    
    //! \param [in] r - component whose argument is being requested.
    //! \return the current index into component <b><var>r</var></b>.
    KOKKOS_INLINE_FUNCTION const ordinal_type & argument(const ordinal_type &r) const
    {
      return arguments_[r];
    }
    
    //! Note: enumerationIndex() matches the ordering in TensorData.  This is different from the order in which this iterator proceeds through the tensor arguments.
    KOKKOS_INLINE_FUNCTION ordinal_type enumerationIndex() const
    {
      // commented-out code belongs to implementation with rightmost argument as the fastest-moving.  We may want to support this as an option.
//      ordinal_type i = 0;
//      for (ordinal_type r=0; r<numTensorComponents_-1; r++)
//      {
//        i += arguments_[r];
//        i *= bounds_[r+1];
//      }
//      i += arguments_[numTensorComponents_-1];
//      return i;
      
      // TensorData's numbering has the leftmost argument as the fastest-moving
      // We return that numbering here.
      ordinal_type i = 0;
      for (ordinal_type r=numTensorComponents_-1; r>0; r--)
      {
        i += arguments_[r];
        i *= bounds_[r-1];
      }
      i += arguments_[0];
      return i;
    }
    
    //! Note: relativeEnumerationIndex() matches the ordering in TensorData.  This is different from the order in which this iterator proceeds through the tensor arguments.
    KOKKOS_INLINE_FUNCTION ordinal_type relativeEnumerationIndex(const ordinal_type &startingComponent) const
    {
      // commented-out code belongs to implementation with rightmost argument as the fastest-moving.  We may want to support this as an option.
//      ordinal_type i = 0;
//      for (ordinal_type r=startingComponent; r<numTensorComponents_-1; r++)
//      {
//        i += arguments_[r];
//        i *= bounds_[r+1];
//      }
//      i += arguments_[numTensorComponents_-1];
//      return i;
      
      // TensorData's numbering has the leftmost argument as the fastest-moving
      // We return that numbering here.
      ordinal_type i = 0;
      for (ordinal_type r=numTensorComponents_-1; r>startingComponent; r--)
      {
        i += arguments_[r];
        i *= bounds_[r-1];
      }
      i += arguments_[startingComponent];
      return i;
    }
    
    //! total number of enumeration indices with arguments prior to the startingComponent fixed
    KOKKOS_INLINE_FUNCTION ordinal_type relativeEnumerationSpan(const ordinal_type &startingComponent) const
    {
      ordinal_type i = 1;
      for (ordinal_type r=startingComponent; r<numTensorComponents_; r++)
      {
        i *= bounds_[r];
      }
      return i;
    }
    
    //! Resets the location to index 0 in each dimension, starting from the specified dimension.
    //! \param [in] from_rank_number - the first dimension in which to set the index to 0.
    KOKKOS_INLINE_FUNCTION
    void reset(ordinal_type from_component_number=0)
    {
      for (ordinal_type r=from_component_number; r<numTensorComponents_; r++)
      {
        arguments_[r] = 0;
      }
    }
    
    //! Sets the current argument in the specified component.
    //! \param [in] r - the component index
    //! \param [in] i - the argument to set for component <b>r<\b>
    KOKKOS_INLINE_FUNCTION
    void setArgumentForComponent(const ordinal_type &r, const ordinal_type &i)
    {
      arguments_[r] = i;
    }
    
    /** \brief  Sets the enumeration index; this refers to a 1D enumeration of the possible in-bound arguments.
        \param [in] enumerationIndex - the index to corresponding to the arguments to set.
        \see enumerationIndex()
     \note WARNING: this method does not have any tests against it.  We should add a test.
    */
    KOKKOS_INLINE_FUNCTION
    void setEnumerationIndex(const ordinal_type &enumerationIndex)
    {
      ordinal_type remainder = enumerationIndex;
      for (ordinal_type d=0; d<numTensorComponents_; d++)
      {
        arguments_[d] = remainder % bounds_[d];
        remainder  /= bounds_[d];
      }
    }
    
    /** \brief Sets a subset of this iterator's component arguments to match the component arguments from <var>otherArgumentIterator</var>.
     */
    KOKKOS_INLINE_FUNCTION
    void copyArguments(TensorArgumentIterator &otherArgumentIterator, const ordinal_type &r0_from, const ordinal_type &r0_to, const ordinal_type &numArguments)
    {
      for (ordinal_type i=0; i<numArguments; i++)
      {
        arguments_[r0_to + i] = otherArgumentIterator.argument(r0_from + i);
      }
    }
    
  };
} // namespace Intrepid2

#endif /* Intrepid2_TensorArgumentIterator_h */
