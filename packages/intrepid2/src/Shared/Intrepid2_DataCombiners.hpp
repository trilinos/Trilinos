//
//  Intrepid2_DataCombiners.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 5/31/23.
//

#ifndef Intrepid2_DataCombiners_h
#define Intrepid2_DataCombiners_h

/** \file  Intrepid2_DataCombiners.hpp
   \brief  Defines functors that help with combinations of Data objects, such as in-place sums and products.
   \author Created by N.V. Roberts.
*/

namespace Intrepid2 {
  template<class DataScalar,typename DeviceType>
  class Data;

/** \enum  Intrepid2::DataVariationType
    \brief Enumeration to indicate how data varies in a particular dimension of an Intrepid2::Data object.
 CONSTANT indicates that the data does not vary; MODULAR indicates that it varies according to some (separately specified) modulus; BLOCK_PLUS_DIAGONAL allows specification of a matrix that has a non-diagonal block followed by a diagonal block; GENERAL indicates arbitrary variation.
 
 To give some examples for a Data object containing the Jacobian for reference-to-physical space mappings:
 - CONSTANT could be used in the point dimension for an affine transformation
 - MODULAR could be used for the cell dimension for a uniform grid that has been subdivided into simplices
 - BLOCK_PLUS_DIAGONAL could be used for the coordinate dimensions for an arbitrary 2D mesh that has been orthogonally extruded in the z dimension (resulting in diagonal entries in the final row and column of the Jacobian matrix)
 - GENERAL should be used in any dimension in which the data varies in a way not captured by the other options
*/
  enum DataVariationType
  {
    CONSTANT            /// does not vary
  , MODULAR             /// varies according to modulus of the index
  , BLOCK_PLUS_DIAGONAL /// one of two dimensions in a matrix; bottom-right part of matrix is diagonal
  , GENERAL             /// arbitrary variation
  };

/** \struct  Intrepid2::DimensionInfo
    \brief Struct expressing all variation information about a Data object in a single dimension, including its logical extent and storage extent.
*/
  struct DimensionInfo
  {
    int logicalExtent;
    DataVariationType variationType;
    int dataExtent;
    int variationModulus; // should be equal to dataExtent variationType other than MODULAR and CONSTANT
    int blockPlusDiagonalLastNonDiagonal = -1; // only relevant for variationType == BLOCK_PLUS_DIAGONAL
  };

  //! Returns DimensionInfo for a Data container that combines (through multiplication, say, or addition) the two specified DimensionInfo specifications in one of its dimensions.
  KOKKOS_INLINE_FUNCTION
  DimensionInfo combinedDimensionInfo(const DimensionInfo &myData, const DimensionInfo &otherData)
  {
    const int myNominalExtent    = myData.logicalExtent;
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(myNominalExtent != otherData.logicalExtent, std::invalid_argument, "both Data objects must match in their logical extent in the specified dimension");
#endif
    const DataVariationType & myVariation    = myData.variationType;
    const DataVariationType & otherVariation = otherData.variationType;
    
    const int & myVariationModulus    = myData.variationModulus;
    const int & otherVariationModulus = otherData.variationModulus;
    
    int myDataExtent    = myData.dataExtent;
    int otherDataExtent = otherData.dataExtent;
    
    DimensionInfo combinedDimensionInfo;
    combinedDimensionInfo.logicalExtent = myNominalExtent;
    
    switch (myVariation)
    {
      case CONSTANT:
        switch (otherVariation)
        {
          case CONSTANT:
          case MODULAR:
          case GENERAL:
          case BLOCK_PLUS_DIAGONAL:
            combinedDimensionInfo = otherData;
        }
        break;
      case MODULAR:
        switch (otherVariation)
        {
          case CONSTANT:
            combinedDimensionInfo = myData;
            break;
          case MODULAR:
            if (myVariationModulus == otherVariationModulus)
            {
              // in this case, expect both to have the same data extent
              INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(myDataExtent != otherDataExtent, std::logic_error, "Unexpected data extent/modulus combination");
              combinedDimensionInfo.variationType    = MODULAR;
              combinedDimensionInfo.dataExtent       = myDataExtent;
              combinedDimensionInfo.variationModulus = myVariationModulus;
            }
            else
            {
              // both modular with two different moduli
              // we could do something clever with e.g. least common multiples, but for now we just use GENERAL
              // (this is not a use case we anticipate being a common one)
              combinedDimensionInfo.variationType    = GENERAL;
              combinedDimensionInfo.dataExtent       = myNominalExtent;
              combinedDimensionInfo.variationModulus = myNominalExtent;
            }
            break;
          case BLOCK_PLUS_DIAGONAL:
            combinedDimensionInfo.variationType    = GENERAL;
            combinedDimensionInfo.dataExtent       = myNominalExtent;
            combinedDimensionInfo.variationModulus = myNominalExtent;
            break;
          case GENERAL:
            // otherData is GENERAL: its info dominates
            combinedDimensionInfo = otherData;
            break;
        }
        break;
      case BLOCK_PLUS_DIAGONAL:
        switch (otherVariation)
        {
          case CONSTANT:
            combinedDimensionInfo = myData;
            break;
          case MODULAR:
            combinedDimensionInfo.variationType    = GENERAL;
            combinedDimensionInfo.dataExtent       = myNominalExtent;
            combinedDimensionInfo.variationModulus = myNominalExtent;
            break;
          case GENERAL:
            // otherData is GENERAL: its info dominates
            combinedDimensionInfo = otherData;
            break;
          case BLOCK_PLUS_DIAGONAL:
            combinedDimensionInfo.variationType    = GENERAL;
            combinedDimensionInfo.dataExtent       = max(myDataExtent,otherDataExtent);
            combinedDimensionInfo.variationModulus = combinedDimensionInfo.dataExtent;
            // for this case, we want to take the minimum of the two Data objects' blockPlusDiagonalLastNonDiagonal as the combined object's blockPlusDiagonalLastNonDiagonal
            combinedDimensionInfo.blockPlusDiagonalLastNonDiagonal = min(myData.blockPlusDiagonalLastNonDiagonal, otherData.blockPlusDiagonalLastNonDiagonal);
        }
        break;
      case GENERAL:
        switch (otherVariation)
        {
          case CONSTANT:
          case MODULAR:
          case GENERAL:
          case BLOCK_PLUS_DIAGONAL:
            combinedDimensionInfo = myData;
        }
    }
    return combinedDimensionInfo;
  }

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
  //        std::cout << "ThisUnderlyingViewType: " << typeid(ThisUnderlyingViewType).name() << std::endl;
  //        std::cout << "AUnderlyingViewType:    " << typeid(AUnderlyingViewType).name() << std::endl;
  //        std::cout << "BUnderlyingViewType:    " << typeid(BUnderlyingViewType).name() << std::endl;
      
  //        std::cout << "ThisUnderlyingViewType: " << TypeParseTraits<ThisUnderlyingViewType>::name << std::endl;
  //        std::cout << "AUnderlyingViewType:    " << TypeParseTraits<AUnderlyingViewType>::name    << std::endl;
  //        std::cout << "BUnderlyingViewType:    " << TypeParseTraits<BUnderlyingViewType>::name    << std::endl;
//      std::cout << std::flush;
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

// static class for combining two Data objects in various ways
//  template <class DataScalar,typename DeviceType>
//  class DataCombiner
//{
//public:
//  //! storeInPlaceCombination implementation for rank < 7, with compile-time underlying views and argument interpretation.  Intended for internal and expert use.
//  template<class BinaryOperator, class PolicyType, class ThisUnderlyingViewType, class AUnderlyingViewType, class BUnderlyingViewType,
//           class ArgExtractorThis, class ArgExtractorA, class ArgExtractorB>
//  static void storeInPlaceCombination(PolicyType &policy, ThisUnderlyingViewType &this_underlying,
//                                      AUnderlyingViewType &A_underlying, BUnderlyingViewType &B_underlying,
//                                      BinaryOperator &binaryOperator, ArgExtractorThis argThis, ArgExtractorA argA, ArgExtractorB argB)
//  {
//    using Functor = InPlaceCombinationFunctor<BinaryOperator, ThisUnderlyingViewType, AUnderlyingViewType, BUnderlyingViewType, ArgExtractorThis, ArgExtractorA, ArgExtractorB>;
//    Functor functor(this_underlying, A_underlying, B_underlying, binaryOperator);
//    Kokkos::parallel_for("compute in-place", policy, functor);
//  }
//};

template<class Scalar>
struct ScalarProductFunctor
{
  KOKKOS_INLINE_FUNCTION
  Scalar operator()(const Scalar &a, const Scalar &b) const
  {
    return a * b;
  }
};

template<class Scalar>
struct ScalarSumFunctor
{
  KOKKOS_INLINE_FUNCTION
  Scalar operator()(const Scalar &a, const Scalar &b) const
  {
    return a + b;
  }
};



} // end namespace Intrepid2

#endif /* Intrepid2_DataCombiners_h */
