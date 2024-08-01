// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_TensorViewIterator.hpp
    \brief  Implementation of support for traversing component views alongside a view that represents a combination of those views; support is provided for rank-preserving, rank-increasing, and rank-reducing operations.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_TensorViewIterator_h
#define Intrepid2_TensorViewIterator_h

#include "Intrepid2_DeviceAssert.hpp"
#include "Intrepid2_ViewIterator.hpp"

#include <Kokkos_Core.hpp>

#include <vector>

namespace Intrepid2
{
  /** \class  Intrepid2::TensorViewIterator
      \brief  A helper class that allows iteration over three Kokkos Views simultaneously, according to
              tensor combination rules:
              - component view 1,
              - component view 2, and
              - the combination tensor view
   
      @see Intrepid2::ViewIterator
   
      \remark This class is under active development; its interface should be understood to be in flux.
              Furthermore, this may eventually move to another Trilinos package.
   */
  template<class TensorViewType, class ViewType1, class ViewType2 ,typename ScalarType>
  class TensorViewIterator
  {
  public:
    enum RankCombinationType : int
    {
      DIMENSION_MATCH,
      TENSOR_PRODUCT,
      TENSOR_CONTRACTION
    };
    using RankCombinationViewType = Kokkos::View<RankCombinationType*, typename TensorViewType::device_type>;
  protected:
    
    ViewIterator<TensorViewType, ScalarType> tensor_view_iterator_;
    ViewIterator<ViewType1, ScalarType> view1_iterator_;
    ViewIterator<ViewType2, ScalarType> view2_iterator_;
    
    RankCombinationViewType rank_combination_types_;
  public:
    /** \brief Constructor
        \param [in] tensor_view            - the view that stores the tensor combination
        \param [in] view1                  - first component view
        \param [in] view2                  - second component view
        \param [in] rank_combination_types - vector with length equal to the maximum rank of the views provided, specifying how component views should be combined to produce tensor view.
     rank_combination_types entries can be as follows:
     - DIMENSION_MATCH:    component views and tensor view all have same extent in this rank, and should iterate in tandem.
     - TENSOR_PRODUCT:     if component views have extents a and b in this rank, the tensor view has extent a*b.
     - TENSOR_CONTRACTION: component views have matching extents in this rank; tensor view has extent 1.
     
     For TENSOR_PRODUCT combinations in a given rank, view1 increments first.  For example, if view1 has entries [x,y]
     and view2 has entries [0,1] in a given TENSOR_PRODUCT rank, there will be four entries in the tensor_view, in the following order:
     - (x,0)
     - (y,0)
     - (x,1)
     - (y,1)
     
     Tensor contractions are only allowed in the final dimension(s); if rank d has a TENSOR_CONTRACTION entry, then all subsequent ranks must also contract.
    */
    KOKKOS_INLINE_FUNCTION
    TensorViewIterator(TensorViewType tensor_view, ViewType1 view1, ViewType2 view2,
        RankCombinationViewType rank_combination_types)
    :
    tensor_view_iterator_(tensor_view),
    view1_iterator_(view1),
    view2_iterator_(view2),
    rank_combination_types_(rank_combination_types)
    {
      // rank_combination_type should have length equal to the maximum rank of the views provided
      /*
       Examples:
       1. vector dot product in third dimension: {DIMENSION_MATCH, DIMENSION_MATCH, TENSOR_CONTRACTION}
       - view1 and view2 should both be rank 3, and should match in all dimensions
       - tensor_view should be rank 2, and should match view1 and view2 in first two dimensions
       2. vector outer product in third dimension: {DIMENSION_MATCH, DIMENSION_MATCH, TENSOR_PRODUCT}
       - view1 and view2 should both be rank 3, and should match in first two dimensions
       - tensor_view should be rank 3, and should match view1 and view2 in first two dimensions
       - in third dimension, tensor_view should have dimension equal to the product of the third dimension of view1 and the third dimension of view2
       3. rank-3 view1 treated as vector times scalar rank-2 view2: {DIMENSION_MATCH, DIMENSION_MATCH, TENSOR_PRODUCT}
       - here, the rank-2 view2 is interpreted as having an extent 1 third dimension
       
       We only allow TENSOR_CONTRACTION in final dimension(s)
       */
      // check that the above rules are satisfied:
      unsigned max_component_rank = (view1.rank() > view2.rank()) ? view1.rank() : view2.rank();
      unsigned max_rank           = (tensor_view.rank() > max_component_rank) ? tensor_view.rank() : max_component_rank;
      
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(rank_combination_types.extent(0) != max_rank, std::invalid_argument, "need to provide RankCombinationType for the largest-rank View");
      
      unsigned expected_rank = 0;
      bool contracting = false;
      for (unsigned d=0; d<rank_combination_types.extent(0); d++)
      {
        if (rank_combination_types[d] == TENSOR_CONTRACTION)
        {
          // check that view1 and view2 agree on the length of this dimension
          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(view1.extent_int(d) != view2.extent_int(d), std::invalid_argument, "Contractions can only occur along ranks of equal length");
          contracting = true;
        }
        else
        {
          if (!contracting)
          {
            expected_rank++;
            if (rank_combination_types[d] == TENSOR_PRODUCT)
            {
              INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(tensor_view.extent_int(d) != view1.extent_int(d) * view2.extent_int(d), std::invalid_argument, "For TENSOR_PRODUCT rank combination, the tensor View must have length in that dimension equal to the product of the two component views in that dimension");
            }
            else // matching
            {
              INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(view1.extent_int(d) != view2.extent_int(d), std::invalid_argument, "For DIMENSION_MATCH rank combination, all three views must have length equal to each other in that rank");
              INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(tensor_view.extent_int(d) != view1.extent_int(d), std::invalid_argument, "For DIMENSION_MATCH rank combination, all three views must have length equal to each other in that rank");
            }
          }
          else
          {
            INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(contracting, std::invalid_argument, "encountered a non-contraction rank combination after a contraction; contractions can only go at the end");
          }
        }
      }
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(expected_rank != tensor_view.rank(), std::invalid_argument, "Tensor view does not match expected rank");
    }
    
    //!
    //! \return the leftmost rank ordinal whose index will change on next increment in any of the views (-1 if next increment will go out of bounds)
    KOKKOS_INLINE_FUNCTION
    int nextIncrementRank()
    {
      int view2_next_increment_rank = view2_iterator_.nextIncrementRank();
      int view1_next_increment_rank = view1_iterator_.nextIncrementRank();
      if (view2_next_increment_rank > view1_next_increment_rank) return view2_next_increment_rank;
      else return view1_next_increment_rank;
    }
    
    //! Move to the next location.  Note that during contractions, the tensor view location will not change; the returned rank corresponds to the other two views in that case.
    //! \return the rank of the leftmost index that was changed in any of the views; -1 if increment reached the end of the views
    KOKKOS_INLINE_FUNCTION
    int increment()
    {
      // proceed to the next view1/view2 combination
      // where we're doing a dimension match, then all three iterators should increment in tandem
      // where we're doing a contraction, view1/view2 should increment in tandem, while tensor_view should be fixed
      // where we're doing a tensor product, view1 and tensor_view increment in tandem, while view2 is fixed
      
      // note that regardless of the choice, view1 should be incremented, with one exception:
      //   If we are doing a tensor product, then view1 can be understood to be in an interior for loop, and it should loop around.
      //   We can detect this by checking which the least rank that would be updated -- if view2's least rank exceeds view1's, then:
      //   - view1 should be reset, AND
      //   - view2 should be incremented (as should the tensor view)
      int view2_next_increment_rank = view2_iterator_.nextIncrementRank();
      int view1_next_increment_rank = view1_iterator_.nextIncrementRank();
      if (view2_next_increment_rank > view1_next_increment_rank)
      {
        // if we get here, we should be doing a tensor product in the view2 rank that will change
        device_assert(rank_combination_types_[view2_next_increment_rank]==TENSOR_PRODUCT);
        view1_iterator_.reset(view2_next_increment_rank); // set to 0 from the tensor product rank inward -- this is "looping around"
        view2_iterator_.increment();
        tensor_view_iterator_.increment();
        return view2_next_increment_rank;
      }
      else
      {
        int view1_rank_change = view1_iterator_.increment();
        if (view1_rank_change >= 0)
        {
          switch (rank_combination_types_[view1_rank_change])
          {
            case DIMENSION_MATCH:
              view2_iterator_.increment();
              tensor_view_iterator_.increment();
              break;
            case TENSOR_PRODUCT:
              // view1 increments fastest; the only time we increment view2 is when view1 loops around; we handle that above
              tensor_view_iterator_.increment();
              break;
            case TENSOR_CONTRACTION:
              // view1 and view2 increment in tandem; we don't increment tensor_view while contraction is taking place
              view2_iterator_.increment();
          }
        }
        return view1_rank_change;
      }
    }
    
    //! Sets the current location in all three views to the specified location.
    //! \param [in] location - the location as a 7-element array value.
    KOKKOS_INLINE_FUNCTION
    void setLocation(const Kokkos::Array<int,7> location)
    {
      view1_iterator_.setLocation(location);
      view2_iterator_.setLocation(location);
      tensor_view_iterator_.setLocation(location);
    }
    
    //! Sets the current location in the two component views to the specified locations, and sets the tensor view location to the location corresponding to the two component locations.
    //! \param [in] location1 - the location in view1 as a 7-element array value.
    //! \param [in] location2 - the location in view2 as a 7-element array value.
    KOKKOS_INLINE_FUNCTION
    void setLocation(Kokkos::Array<int,7> location1, Kokkos::Array<int,7> location2)
    {
      view1_iterator_.setLocation(location1);
      view2_iterator_.setLocation(location2);
      Kokkos::Array<int,7> tensor_location = location1;
      for (unsigned d=0; d<rank_combination_types_.extent(0); d++)
      {
        switch (rank_combination_types_[d])
        {
          case TENSOR_PRODUCT:
            // view1 index is fastest-moving:
            tensor_location[d] = location2[d] * view1_iterator_.getExtent(d) + location1[d];
            break;
          case DIMENSION_MATCH:
            // we copied location1 into tensor_location to initialize -- that's correct in this dimension
            break;
          case TENSOR_CONTRACTION:
            tensor_location[d] = 0;
            break;
        }
      }
#ifdef HAVE_INTREPID2_DEBUG
      // check that the location makes sense
      for (unsigned d=0; d<rank_combination_types_.extent(0); d++)
      {
        switch (rank_combination_types_[d])
        {
          case TENSOR_PRODUCT:
            // in this case, the two indices are independent
            break;
          case DIMENSION_MATCH:
          case TENSOR_CONTRACTION:
            device_assert(location1[d] == location2[d]);
            break;
        }
        // let's check that the indices are in bounds:
        device_assert(location1[d] < view1_iterator_.getExtent(d));
        device_assert(location2[d] < view2_iterator_.getExtent(d));
        device_assert(tensor_location[d] < tensor_view_iterator_.getExtent(d));
      }
#endif
      tensor_view_iterator_.setLocation(tensor_location);
    }
    
    //! Getter for view1.
    //! \return the value at the current location in view1
    KOKKOS_INLINE_FUNCTION
    ScalarType getView1Entry()
    {
      return view1_iterator_.get();
    }
    
    //! Getter for view1.
    //! \return the value at the current location in view2
    KOKKOS_INLINE_FUNCTION
    ScalarType getView2Entry()
    {
      return view2_iterator_.get();
    }
    
    //! Setter for tensor view.
    //! \param [in] value - the value to set in the tensor view
    KOKKOS_INLINE_FUNCTION
    void set(ScalarType value)
    {
      tensor_view_iterator_.set(value);
    }
  };

} // namespace Intrepid2

#endif /* Intrepid2_TensorViewIterator_h */
