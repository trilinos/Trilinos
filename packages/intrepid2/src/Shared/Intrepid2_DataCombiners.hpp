// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  Intrepid2_DataCombiners.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 5/31/23.
//

#ifndef Intrepid2_DataCombiners_hpp
#define Intrepid2_DataCombiners_hpp

/** \file  Intrepid2_DataCombiners.hpp
   \brief  Defines functors that help with combinations of Data objects, such as in-place sums and products.
   \author Created by N.V. Roberts.
*/

#include "Intrepid2_ArgExtractor.hpp"
#include "Intrepid2_Data.hpp"
#include "Intrepid2_DataFunctors.hpp"
#include "Intrepid2_DataVariationType.hpp"
#include "Intrepid2_ScalarView.hpp"

namespace Intrepid2 {
  template<class DataScalar,typename DeviceType>
  class Data;

  template<class BinaryOperator, class ThisUnderlyingViewType, class AUnderlyingViewType, class BUnderlyingViewType,
           class ArgExtractorThis, class ArgExtractorA, class ArgExtractorB, bool includeInnerLoop=false>
  struct InPlaceCombinationFunctor
  {
  private:
    ThisUnderlyingViewType this_underlying_;
    AUnderlyingViewType A_underlying_;
    BUnderlyingViewType B_underlying_;
    BinaryOperator binaryOperator_;
    int innerLoopSize_;
  public:
    InPlaceCombinationFunctor(ThisUnderlyingViewType this_underlying, AUnderlyingViewType A_underlying, BUnderlyingViewType B_underlying,
                              BinaryOperator binaryOperator)
    :
    this_underlying_(this_underlying),
    A_underlying_(A_underlying),
    B_underlying_(B_underlying),
    binaryOperator_(binaryOperator)
    {
      INTREPID2_TEST_FOR_EXCEPTION(includeInnerLoop,std::invalid_argument,"If includeInnerLoop is true, must specify the size of the inner loop");
    }
    
    InPlaceCombinationFunctor(ThisUnderlyingViewType this_underlying, AUnderlyingViewType A_underlying, BUnderlyingViewType B_underlying,
                              BinaryOperator binaryOperator, int innerLoopSize)
    :
    this_underlying_(this_underlying),
    A_underlying_(A_underlying),
    B_underlying_(B_underlying),
    binaryOperator_(binaryOperator),
    innerLoopSize_(innerLoopSize)
    {
      INTREPID2_TEST_FOR_EXCEPTION(includeInnerLoop,std::invalid_argument,"If includeInnerLoop is true, must specify the size of the inner loop");
    }
    
    template<class ...IntArgs, bool M=includeInnerLoop>
    KOKKOS_INLINE_FUNCTION
    enable_if_t<!M, void>
    operator()(const IntArgs&... args) const
    {
      auto      & result = ArgExtractorThis::get( this_underlying_, args... );
      const auto & A_val =    ArgExtractorA::get(    A_underlying_, args... );
      const auto & B_val =    ArgExtractorB::get(    B_underlying_, args... );
      
      result = binaryOperator_(A_val,B_val);
    }
    
    template<class ...IntArgs, bool M=includeInnerLoop>
    KOKKOS_INLINE_FUNCTION
    enable_if_t<M, void>
    operator()(const IntArgs&... args) const
    {
      using int_type = std::tuple_element_t<0, std::tuple<IntArgs...>>;
      for (int_type iFinal=0; iFinal<static_cast<int_type>(innerLoopSize_); iFinal++)
      {
        auto      & result = ArgExtractorThis::get( this_underlying_, args..., iFinal );
        const auto & A_val =    ArgExtractorA::get(    A_underlying_, args..., iFinal );
        const auto & B_val =    ArgExtractorB::get(    B_underlying_, args..., iFinal );
        
        result = binaryOperator_(A_val,B_val);
      }
    }
  };

  //! functor definition for the constant-data case.
  template<class BinaryOperator, class ThisUnderlyingViewType, class AUnderlyingViewType, class BUnderlyingViewType>
  struct InPlaceCombinationFunctorConstantCase
  {
  private:
    ThisUnderlyingViewType this_underlying_;
    AUnderlyingViewType A_underlying_;
    BUnderlyingViewType B_underlying_;
    BinaryOperator binaryOperator_;
  public:
    InPlaceCombinationFunctorConstantCase(ThisUnderlyingViewType this_underlying,
                                          AUnderlyingViewType A_underlying,
                                          BUnderlyingViewType B_underlying,
                                          BinaryOperator binaryOperator)
    :
    this_underlying_(this_underlying),
    A_underlying_(A_underlying),
    B_underlying_(B_underlying),
    binaryOperator_(binaryOperator)
    {
      INTREPID2_TEST_FOR_EXCEPTION(this_underlying.extent(0) != 1,std::invalid_argument,"all views for InPlaceCombinationFunctorConstantCase should have rank 1 and extent 1");
      INTREPID2_TEST_FOR_EXCEPTION(A_underlying.extent(0) != 1,std::invalid_argument,"all views for InPlaceCombinationFunctorConstantCase should have rank 1 and extent 1");
      INTREPID2_TEST_FOR_EXCEPTION(B_underlying.extent(0) != 1,std::invalid_argument,"all views for InPlaceCombinationFunctorConstantCase should have rank 1 and extent 1");
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const int arg0) const
    {
      auto & result      = this_underlying_(0);
      const auto & A_val = A_underlying_(0);
      const auto & B_val = B_underlying_(0);
      
      result = binaryOperator_(A_val,B_val);
    }
  };

  //! For use with Data object into which a value will be stored.  We use passThroughBlockDiagonalArgs = true for storeInPlaceCombination().
  template<bool passThroughBlockDiagonalArgs>
  struct FullArgExtractorWritableData
  {
    template<class ViewType, class ...IntArgs>
    static KOKKOS_INLINE_FUNCTION typename ViewType::reference_type get(const ViewType &view, const IntArgs&... intArgs)
    {
      return view.getWritableEntryWithPassThroughOption(passThroughBlockDiagonalArgs, intArgs...);
    }
  };

  //! For use with Data object into which a value will be stored.  We use passThroughBlockDiagonalArgs = true for storeInPlaceCombination().
  template<bool passThroughBlockDiagonalArgs>
  struct FullArgExtractorData
  {
    template<class ViewType, class ...IntArgs>
    static KOKKOS_INLINE_FUNCTION typename ViewType::const_reference_type get(const ViewType &view, const IntArgs&... intArgs)
    {
      return view.getEntryWithPassThroughOption(passThroughBlockDiagonalArgs, intArgs...);
    }
  };

// static class for combining two Data objects using a specified binary operator
  template <class DataScalar,typename DeviceType, class BinaryOperator>
  class DataCombiner
{
  using reference_type       = typename ScalarView<DataScalar,DeviceType>::reference_type;
  using const_reference_type = typename ScalarView<const DataScalar,DeviceType>::reference_type;
public:
  //! storeInPlaceCombination implementation for rank < 7, with compile-time underlying views and argument interpretation.  Intended for internal and expert use.
  template<class PolicyType, class ThisUnderlyingViewType, class AUnderlyingViewType, class BUnderlyingViewType,
           class ArgExtractorThis, class ArgExtractorA, class ArgExtractorB>
  static void storeInPlaceCombination(PolicyType &policy, ThisUnderlyingViewType &this_underlying,
                                      AUnderlyingViewType &A_underlying, BUnderlyingViewType &B_underlying,
                                      BinaryOperator &binaryOperator, ArgExtractorThis argThis, ArgExtractorA argA, ArgExtractorB argB)
  {
    using Functor = InPlaceCombinationFunctor<BinaryOperator, ThisUnderlyingViewType, AUnderlyingViewType, BUnderlyingViewType, ArgExtractorThis, ArgExtractorA, ArgExtractorB>;
    Functor functor(this_underlying, A_underlying, B_underlying, binaryOperator);
    Kokkos::parallel_for("compute in-place", policy, functor);
  }
  
  //! storeInPlaceCombination with compile-time rank -- implementation for rank < 7.
  template<int rank>
  static
  enable_if_t<rank != 7, void>
  storeInPlaceCombination(Data<DataScalar,DeviceType> &thisData, const Data<DataScalar,DeviceType> &A, const Data<DataScalar,DeviceType> &B, BinaryOperator binaryOperator)
  {
    auto policy = thisData.template dataExtentRangePolicy<rank>();
    
    const bool A_1D          = A.getUnderlyingViewRank() == 1;
    const bool B_1D          = B.getUnderlyingViewRank() == 1;
    const bool this_1D       = thisData.getUnderlyingViewRank() == 1;
    const bool A_constant    = A_1D && (A.getUnderlyingViewSize() == 1);
    const bool B_constant    = B_1D && (B.getUnderlyingViewSize() == 1);
    const bool this_constant = this_1D && (thisData.getUnderlyingViewSize() == 1);
    const bool A_full        = A.underlyingMatchesLogical();
    const bool B_full        = B.underlyingMatchesLogical();
    const bool this_full     = thisData.underlyingMatchesLogical();
    
    const ConstantArgExtractor<reference_type> constArg;
    
    const FullArgExtractor<reference_type> fullArgs;
    const FullArgExtractorData<true> fullArgsData; // true: pass through block diagonal args.  This is due to the behavior of dataExtentRangePolicy() for block diagonal args.
    const FullArgExtractorWritableData<true> fullArgsWritable; // true: pass through block diagonal args.  This is due to the behavior of dataExtentRangePolicy() for block diagonal args.
    
    const SingleArgExtractor<reference_type,0> arg0;
    const SingleArgExtractor<reference_type,1> arg1;
    const SingleArgExtractor<reference_type,2> arg2;
    const SingleArgExtractor<reference_type,3> arg3;
    const SingleArgExtractor<reference_type,4> arg4;
    const SingleArgExtractor<reference_type,5> arg5;
    
    // this lambda returns -1 if there is not a rank-1 underlying view whose data extent matches the logical extent in the corresponding dimension;
    // otherwise, it returns the logical index of the corresponding dimension.
    auto get1DArgIndex = [](const Data<DataScalar,DeviceType> &data) -> int
    {
      const auto & variationTypes = data.getVariationTypes();
      for (int d=0; d<rank; d++)
      {
        if (variationTypes[d] == GENERAL)
        {
          return d;
        }
      }
      return -1;
    };
    if (this_constant)
    {
      // then A, B are constant, too
      auto thisAE = constArg;
      auto AAE    = constArg;
      auto BAE    = constArg;
      auto & this_underlying = thisData.template getUnderlyingView<1>();
      auto & A_underlying    = A.template getUnderlyingView<1>();
      auto & B_underlying    = B.template getUnderlyingView<1>();
      storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, AAE, BAE);
    }
    else if (this_full && A_full && B_full)
    {
      auto thisAE = fullArgs;
      auto AAE    = fullArgs;
      auto BAE    = fullArgs;
      
      auto & this_underlying = thisData.template getUnderlyingView<rank>();
      auto & A_underlying    = A.template getUnderlyingView<rank>();
      auto & B_underlying    = B.template getUnderlyingView<rank>();
      
      storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, AAE, BAE);
    }
    else if (A_constant)
    {
      auto AAE = constArg;
      auto & A_underlying = A.template getUnderlyingView<1>();
      if (this_full)
      {
        auto thisAE = fullArgs;
        auto & this_underlying = thisData.template getUnderlyingView<rank>();
        
        if (B_full)
        {
          auto BAE = fullArgs;
          auto & B_underlying = B.template getUnderlyingView<rank>();
          storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, AAE, BAE);
        }
        else // this_full, not B_full: B may have modular data, etc.
        {
          auto BAE = fullArgsData;
          storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, thisAE, AAE, BAE);
        }
      }
      else // this is not full
      {
        // below, we optimize for the case of 1D data in B, when A is constant.  Still need to handle other cases…
        if (B_1D && (get1DArgIndex(B) != -1) )
        {
          // since A is constant, that implies that this_1D is true, and has the same 1DArgIndex
          const int argIndex = get1DArgIndex(B);
          auto & B_underlying    = B.template getUnderlyingView<1>();
          auto & this_underlying = thisData.template getUnderlyingView<1>();
          switch (argIndex)
          {
            case 0: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg0, AAE, arg0); break;
            case 1: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg1, AAE, arg1); break;
            case 2: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg2, AAE, arg2); break;
            case 3: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg3, AAE, arg3); break;
            case 4: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg4, AAE, arg4); break;
            case 5: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg5, AAE, arg5); break;
            default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid/unexpected arg index");
          }
        }
        else
        {
          // since storing to Data object requires a call to getWritableEntry(), we use FullArgExtractorWritableData
          auto thisAE = fullArgsWritable;
          auto BAE    = fullArgsData;
          storeInPlaceCombination(policy, thisData, A_underlying, B, binaryOperator, thisAE, AAE, BAE);
        }
      }
    }
    else if (B_constant)
    {
      auto BAE = constArg;
      auto & B_underlying = B.template getUnderlyingView<1>();
      if (this_full)
      {
        auto thisAE = fullArgs;
        auto & this_underlying = thisData.template getUnderlyingView<rank>();
        if (A_full)
        {
          auto AAE = fullArgs;
          auto & A_underlying = A.template getUnderlyingView<rank>();
          
          storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, AAE, BAE);
        }
        else  // this_full, not A_full: A may have modular data, etc.
        {
          // use A (the Data object).  This could be further optimized by using A's underlying View and an appropriately-defined ArgExtractor.
          auto AAE = fullArgsData;
          storeInPlaceCombination(policy, this_underlying, A, B_underlying, binaryOperator, thisAE, AAE, BAE);
        }
      }
      else // this is not full
      {
        // below, we optimize for the case of 1D data in A, when B is constant.  Still need to handle other cases…
        if (A_1D && (get1DArgIndex(A) != -1) )
        {
          // since B is constant, that implies that this_1D is true, and has the same 1DArgIndex as A
          const int argIndex = get1DArgIndex(A);
          auto & A_underlying    = A.template getUnderlyingView<1>();
          auto & this_underlying = thisData.template getUnderlyingView<1>();
          switch (argIndex)
          {
            case 0: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg0, arg0, BAE); break;
            case 1: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg1, arg1, BAE); break;
            case 2: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg2, arg2, BAE); break;
            case 3: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg3, arg3, BAE); break;
            case 4: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg4, arg4, BAE); break;
            case 5: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg5, arg5, BAE); break;
            default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid/unexpected arg index");
          }
        }
        else
        {
          // since storing to Data object requires a call to getWritableEntry(), we use FullArgExtractorWritableData
          auto thisAE = fullArgsWritable;
          auto AAE    = fullArgsData;
          storeInPlaceCombination(policy, thisData, A, B_underlying, binaryOperator, thisAE, AAE, BAE);
        }
      }
    }
    else // neither A nor B constant
    {
      if (this_1D && (get1DArgIndex(thisData) != -1))
      {
        // possible ways that "this" could have full-extent, 1D data
        // 1. A constant, B 1D
        // 2. A 1D, B constant
        // 3. A 1D, B 1D
        // The constant possibilities are already addressed above, leaving us with (3).  Note that A and B don't have to be full-extent, however
        const int argThis = get1DArgIndex(thisData);
        const int argA    = get1DArgIndex(A); // if not full-extent, will be -1
        const int argB    = get1DArgIndex(B); // ditto
        
        auto & A_underlying    = A.template getUnderlyingView<1>();
        auto & B_underlying    = B.template getUnderlyingView<1>();
        auto & this_underlying = thisData.template getUnderlyingView<1>();
        if ((argA != -1) && (argB != -1))
        {
#ifdef INTREPID2_HAVE_DEBUG
          INTREPID2_TEST_FOR_EXCEPTION(argA != argThis, std::logic_error, "Unexpected 1D arg combination.");
          INTREPID2_TEST_FOR_EXCEPTION(argB != argThis, std::logic_error, "Unexpected 1D arg combination.");
#endif
          switch (argThis)
          {
            case 0: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg0, arg0, arg0); break;
            case 1: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg1, arg1, arg1); break;
            case 2: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg2, arg2, arg2); break;
            case 3: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg3, arg3, arg3); break;
            case 4: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg4, arg4, arg4); break;
            case 5: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, arg5, arg5, arg5); break;
            default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid/unexpected arg index");
          }
        }
        else if (argA != -1)
        {
          // B is not full-extent in dimension argThis; use the Data object
          switch (argThis)
          {
            case 0: storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, arg0, arg0, fullArgsData); break;
            case 1: storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, arg1, arg1, fullArgsData); break;
            case 2: storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, arg2, arg2, fullArgsData); break;
            case 3: storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, arg3, arg3, fullArgsData); break;
            case 4: storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, arg4, arg4, fullArgsData); break;
            case 5: storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, arg5, arg5, fullArgsData); break;
            default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid/unexpected arg index");
          }
        }
        else
        {
          // A is not full-extent in dimension argThis; use the Data object
          switch (argThis)
          {
            case 0: storeInPlaceCombination(policy, this_underlying, A, B_underlying, binaryOperator, arg0, fullArgsData, arg0); break;
            case 1: storeInPlaceCombination(policy, this_underlying, A, B_underlying, binaryOperator, arg1, fullArgsData, arg1); break;
            case 2: storeInPlaceCombination(policy, this_underlying, A, B_underlying, binaryOperator, arg2, fullArgsData, arg2); break;
            case 3: storeInPlaceCombination(policy, this_underlying, A, B_underlying, binaryOperator, arg3, fullArgsData, arg3); break;
            case 4: storeInPlaceCombination(policy, this_underlying, A, B_underlying, binaryOperator, arg4, fullArgsData, arg4); break;
            case 5: storeInPlaceCombination(policy, this_underlying, A, B_underlying, binaryOperator, arg5, fullArgsData, arg5); break;
            default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid/unexpected arg index");
          }
        }
      }
      else if (this_full)
      {
        // This case uses A,B Data objects; could be optimized by dividing into subcases and using underlying Views with appropriate ArgExtractors.
        auto & this_underlying = thisData.template getUnderlyingView<rank>();
        auto thisAE = fullArgs;
        
        if (A_full)
        {
          auto & A_underlying = A.template getUnderlyingView<rank>();
          auto AAE = fullArgs;
          
          if (B_1D && (get1DArgIndex(B) != -1))
          {
            const int argIndex = get1DArgIndex(B);
            auto & B_underlying = B.template getUnderlyingView<1>();
            switch (argIndex)
            {
              case 0: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, AAE, arg0); break;
              case 1: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, AAE, arg1); break;
              case 2: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, AAE, arg2); break;
              case 3: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, AAE, arg3); break;
              case 4: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, AAE, arg4); break;
              case 5: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, AAE, arg5); break;
              default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid/unexpected arg index");
            }
          }
          else
          {
            // A is full; B is not full, but not constant or full-extent 1D
            // unoptimized in B access:
            FullArgExtractor<const_reference_type> BAE;
            storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, thisAE, AAE, BAE);
          }
        }
        else // A is not full
        {
          if (A_1D && (get1DArgIndex(A) != -1))
          {
            const int argIndex = get1DArgIndex(A);
            auto & A_underlying  = A.template getUnderlyingView<1>();
            if (B_full)
            {
              auto & B_underlying = B.template getUnderlyingView<rank>();
              auto BAE = fullArgs;
              switch (argIndex)
              {
                case 0: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, arg0, BAE); break;
                case 1: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, arg1, BAE); break;
                case 2: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, arg2, BAE); break;
                case 3: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, arg3, BAE); break;
                case 4: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, arg4, BAE); break;
                case 5: storeInPlaceCombination(policy, this_underlying, A_underlying, B_underlying, binaryOperator, thisAE, arg5, BAE); break;
                default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid/unexpected arg index");
              }
            }
            else
            {
              auto BAE = fullArgsData;
              switch (argIndex)
              {
                case 0: storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, thisAE, arg0, BAE); break;
                case 1: storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, thisAE, arg1, BAE); break;
                case 2: storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, thisAE, arg2, BAE); break;
                case 3: storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, thisAE, arg3, BAE); break;
                case 4: storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, thisAE, arg4, BAE); break;
                case 5: storeInPlaceCombination(policy, this_underlying, A_underlying, B, binaryOperator, thisAE, arg5, BAE); break;
                default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid/unexpected arg index");
              }
            }
          }
          else // A not full, and not full-extent 1D
          {
            // unoptimized in A, B accesses.
            auto AAE    = fullArgsData;
            auto BAE    = fullArgsData;
            storeInPlaceCombination(policy, this_underlying, A, B, binaryOperator, thisAE, AAE, BAE);
          }
        }
      }
      else
      {
        // completely un-optimized case: we use Data objects for this, A, B.
        auto thisAE = fullArgsWritable;
        auto AAE    = fullArgsData;
        auto BAE    = fullArgsData;
        storeInPlaceCombination(policy, thisData, A, B, binaryOperator, thisAE, AAE, BAE);
      }
    }
  }
  
  //! storeInPlaceCombination with compile-time rank -- implementation for rank of 7.  (Not optimized; expectation is this case will be rarely used.)
  template<int rank>
  static
  enable_if_t<rank == 7, void>
  storeInPlaceCombination(Data<DataScalar,DeviceType> &thisData, const Data<DataScalar,DeviceType> &A, const Data<DataScalar,DeviceType> &B, BinaryOperator binaryOperator)
  {
    auto policy = thisData.template dataExtentRangePolicy<rank>();
    
    using DataType = Data<DataScalar,DeviceType>;
    using ThisAE = FullArgExtractorWritableData<true>;
    using AAE    = FullArgExtractor<const_reference_type>;
    using BAE    = FullArgExtractor<const_reference_type>;
    
    const ordinal_type dim6 = thisData.getDataExtent(6);
    const bool includeInnerLoop = true;
    using Functor = InPlaceCombinationFunctor<BinaryOperator, DataType, DataType, DataType, ThisAE, AAE, BAE, includeInnerLoop>;
    Functor functor(thisData, A, B, binaryOperator, dim6);
    Kokkos::parallel_for("compute in-place", policy, functor);
  }
  
  static void storeInPlaceCombination(Data<DataScalar,DeviceType> &thisData, const Data<DataScalar,DeviceType> &A, const Data<DataScalar,DeviceType> &B, BinaryOperator binaryOperator)
  {
    using ExecutionSpace = typename DeviceType::execution_space;

#ifdef INTREPID2_HAVE_DEBUG
    // check logical extents
    for (int d=0; d<rank_; d++)
    {
      INTREPID2_TEST_FOR_EXCEPTION(A.extent_int(d) != thisData.extent_int(d), std::invalid_argument, "A, B, and this must agree on all logical extents");
      INTREPID2_TEST_FOR_EXCEPTION(B.extent_int(d) != thisData.extent_int(d), std::invalid_argument, "A, B, and this must agree on all logical extents");
    }
    // TODO: add some checks that data extent of this suffices to accept combined A + B data.
#endif
    
    const bool this_constant = (thisData.getUnderlyingViewRank() == 1) && (thisData.getUnderlyingViewSize() == 1);

    // we special-case for constant output here; since the constant case is essentially all overhead, we want to avoid as much of the overhead of storeInPlaceCombination() as possible…
    if (this_constant)
    {
      // constant data
      Kokkos::RangePolicy<ExecutionSpace> policy(ExecutionSpace(),0,1); // just 1 entry
      
      auto this_underlying = thisData.template getUnderlyingView<1>();
      auto A_underlying = A.template getUnderlyingView<1>();
      auto B_underlying = B.template getUnderlyingView<1>();
      
      using ConstantCaseFunctor = InPlaceCombinationFunctorConstantCase<decltype(binaryOperator), decltype(this_underlying),
                                                                        decltype(A_underlying), decltype(B_underlying)>;
      
      ConstantCaseFunctor functor(this_underlying, A_underlying, B_underlying, binaryOperator);
      Kokkos::parallel_for("compute in-place", policy,functor);
    }
    else
    {
      switch (thisData.rank())
      {
        case 1: storeInPlaceCombination<1>(thisData, A, B, binaryOperator); break;
        case 2: storeInPlaceCombination<2>(thisData, A, B, binaryOperator); break;
        case 3: storeInPlaceCombination<3>(thisData, A, B, binaryOperator); break;
        case 4: storeInPlaceCombination<4>(thisData, A, B, binaryOperator); break;
        case 5: storeInPlaceCombination<5>(thisData, A, B, binaryOperator); break;
        case 6: storeInPlaceCombination<6>(thisData, A, B, binaryOperator); break;
        case 7: storeInPlaceCombination<7>(thisData, A, B, binaryOperator); break;
        default:
          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::logic_error, "unhandled rank in switch");
      }
    }
  }
};

} // end namespace Intrepid2

// We do ETI for basic double arithmetic on default device.
//template<class Scalar> struct ScalarSumFunctor;
//template<class Scalar> struct ScalarDifferenceFunctor;
//template<class Scalar> struct ScalarProductFunctor;
//template<class Scalar> struct ScalarQuotientFunctor;

extern template class Intrepid2::DataCombiner<double,Kokkos::DefaultExecutionSpace::device_type, Intrepid2::ScalarSumFunctor<double> >;
extern template class Intrepid2::DataCombiner<double,Kokkos::DefaultExecutionSpace::device_type, Intrepid2::ScalarDifferenceFunctor<double> >;
extern template class Intrepid2::DataCombiner<double,Kokkos::DefaultExecutionSpace::device_type, Intrepid2::ScalarProductFunctor<double> >;
extern template class Intrepid2::DataCombiner<double,Kokkos::DefaultExecutionSpace::device_type, Intrepid2::ScalarQuotientFunctor<double> >;

#endif /* Intrepid2_DataCombiners_hpp */
