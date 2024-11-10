// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  Intrepid2_Data.hpp
//  QuadraturePerformance
//
//  Created by Roberts, Nathan V on 8/24/20.
//

#ifndef Intrepid2_Data_h
#define Intrepid2_Data_h

#include "Intrepid2_ArgExtractor.hpp"
#include "Intrepid2_DataDimensionInfo.hpp"
#include "Intrepid2_DataFunctors.hpp"
#include "Intrepid2_DataVariationType.hpp"
#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_Utils.hpp"

/** \file  Intrepid2_Data.hpp
   \brief  Defines the Data class, a wrapper around a Kokkos::View that allows data that is constant or repeating in various logical dimensions to be stored just once, while providing a similar interface to that of View.
   \author Created by N.V. Roberts.
*/

namespace Intrepid2 {

/**
\class  Intrepid2::ZeroView
\brief  A singleton class for a DynRankView containing exactly one zero entry.  (Technically, the entry is DataScalar(), the default value for the scalar type.)  This allows View-wrapping classes to return a reference to zero, even when that zero is not explicitly stored in the wrapped views.
 
This is used by Intrepid2::Data for its getEntry() and getWritableEntry() methods.
 
 \note There is no protection against the zero value being overwritten; perhaps we should add some (i.e., const-qualify DataScalar).  Because of implementation details in Intrepid2::Data, we don't do so yet.
 */
template<class DataScalar,typename DeviceType>
class ZeroView {
public:
  static ScalarView<DataScalar,DeviceType> zeroView()
  {
    static ScalarView<DataScalar,DeviceType> zeroView = ScalarView<DataScalar,DeviceType>("zero",1);
    static bool havePushedFinalizeHook = false;
    if (!havePushedFinalizeHook)
    {
      Kokkos::push_finalize_hook( [=] {
        zeroView = ScalarView<DataScalar,DeviceType>();
      });
      havePushedFinalizeHook = true;
    }
    return zeroView;
  }
};

    /**
      \class  Intrepid2::Data
      \brief  Wrapper around a Kokkos::View that allows data that is constant or repeating in various logical dimensions to be stored just once, while providing a similar interface to that of View.
     
      The Data class distinguishes between the logical extent and the data extent.  For example, one could construct a data container corresponding to constant (cell, point) data with 100 cells
     and 25 points per cell as follows:
          auto cpData = Data(value, Kokkos::Array<int>{100,25});
     The data extent of the container is 1 in every dimension, while the logical extent is 100 in the first dimension, and 25 in the second.  Similarly, the logical rank of the container is 2, but the rank of the
     underlying View is 1.
     
     There are four possible variation types in a logical dimension:
     - GENERAL: the data varies arbitrarily.  The underlying View will have the same extent in its corresponding dimension (which may be distinct from the logical dimension).
     - CONSTANT: the data does not vary.  The underlying View will not have a dimension corresponding to this dimension.
     - MODULAR: the data varies with a modulus.  The underlying View will have a corresponding dimension with extent corresponding to the modulus.
     - BLOCK_PLUS_DIAGONAL: the data varies in this logical dimension and one other, corresponding to a square matrix that has some (possibly trivial) full block, with diagonal entries in the remaining dimensions.  The underlying View will have one dimension corresponding to the two logical dimensions, with extent corresponding to the number of nonzeros in the matrix.
     
  */
  template<class DataScalar,typename DeviceType>
  class Data {
  public:
    using value_type      = DataScalar;
    using execution_space = typename DeviceType::execution_space;
    
    using reference_type       = typename ScalarView<DataScalar,DeviceType>::reference_type;
    using const_reference_type = typename ScalarView<const DataScalar,DeviceType>::reference_type;
  private:
    ordinal_type dataRank_;
    Kokkos::View<DataScalar*,       DeviceType> data1_;  // the rank 1 data that is explicitly stored
    Kokkos::View<DataScalar**,      DeviceType> data2_;  // the rank 2 data that is explicitly stored
    Kokkos::View<DataScalar***,     DeviceType> data3_;  // the rank 3 data that is explicitly stored
    Kokkos::View<DataScalar****,    DeviceType> data4_;  // the rank 4 data that is explicitly stored
    Kokkos::View<DataScalar*****,   DeviceType> data5_;  // the rank 5 data that is explicitly stored
    Kokkos::View<DataScalar******,  DeviceType> data6_;  // the rank 6 data that is explicitly stored
    Kokkos::View<DataScalar*******, DeviceType> data7_;  // the rank 7 data that is explicitly stored
    Kokkos::Array<int,7> extents_;                     // logical extents in each dimension
    Kokkos::Array<DataVariationType,7> variationType_; // for each dimension, whether the data varies in that dimension
    Kokkos::Array<int,7> variationModulus_;            // for each dimension, a value by which indices should be modulused (only used when variationType_ is MODULAR)
    int blockPlusDiagonalLastNonDiagonal_ = -1;        // last row/column that is part of the non-diagonal part of the matrix indicated by BLOCK_PLUS_DIAGONAL (if any dimensions are thus marked)
    
    bool hasNontrivialModulusUNUSED_;  // this is a little nutty, but having this UNUSED member variable improves performance, probably by shifting the alignment of underlyingMatchesLogical_.  This is true with nvcc; it may also be true with Apple clang
    bool underlyingMatchesLogical_;   // if true, this Data object has the same rank and extent as the underlying view
    Kokkos::Array<ordinal_type,7> activeDims_;
    int numActiveDims_; // how many of the 7 entries are actually filled in
    
    ordinal_type rank_;
    
    // we use (const_)reference_type as the return for operator() for performance reasons, especially significant when using Sacado types
    using return_type = const_reference_type;
    
    ScalarView<DataScalar,DeviceType> zeroView_; // one-entry (zero); used to allow getEntry() to return 0 for off-diagonal entries in BLOCK_PLUS_DIAGONAL
    
    //! Returns the number of non-diagonal entries based on the last non-diagonal.  Only applicable for BLOCK_PLUS_DIAGONAL DataVariationType.
    KOKKOS_INLINE_FUNCTION
    static int blockPlusDiagonalNumNondiagonalEntries(const int &lastNondiagonal)
    {
      return (lastNondiagonal + 1) * (lastNondiagonal + 1);
    }
    
    //! //! Returns flattened index of the specified (i,j) matrix entry, assuming that i,j ≤ lastNondiagonal.  Only applicable for BLOCK_PLUS_DIAGONAL DataVariationType.
    KOKKOS_INLINE_FUNCTION
    static int blockPlusDiagonalBlockEntryIndex(const int &lastNondiagonal, const int &numNondiagonalEntries, const int &i, const int &j)
    {
      return i * (lastNondiagonal + 1) + j;
    }
    
    //! Returns flattened index of the specified (i,i) matrix entry, assuming that i > lastNondiagonal.  Only applicable for BLOCK_PLUS_DIAGONAL DataVariationType.
    KOKKOS_INLINE_FUNCTION
    static int blockPlusDiagonalDiagonalEntryIndex(const int &lastNondiagonal, const int &numNondiagonalEntries, const int &i)
    {
      return i - (lastNondiagonal + 1) + numNondiagonalEntries;
    }
    
    //! Returns the extent of the underlying view in the specified dimension.
    KOKKOS_INLINE_FUNCTION
    int getUnderlyingViewExtent(const int &dim) const
    {
      switch (dataRank_)
      {
        case 1: return data1_.extent_int(dim);
        case 2: return data2_.extent_int(dim);
        case 3: return data3_.extent_int(dim);
        case 4: return data4_.extent_int(dim);
        case 5: return data5_.extent_int(dim);
        case 6: return data6_.extent_int(dim);
        case 7: return data7_.extent_int(dim);
        default: return -1;
      }
    }
    
    //! class initialization method.  Called by constructors.
    void setActiveDims()
    {
      zeroView_ = ZeroView<DataScalar,DeviceType>::zeroView(); // one-entry (zero); used to allow getEntry() to return 0 for off-diagonal entries in BLOCK_PLUS_DIAGONAL
      // check that rank is compatible with the claimed extents:
      for (int d=rank_; d<7; d++)
      {
        INTREPID2_TEST_FOR_EXCEPTION(extents_[d] > 1, std::invalid_argument, "Nominal extents may not be > 1 in dimensions beyond the rank of the container");
      }
      
      numActiveDims_ = 0;
      int blockPlusDiagonalCount = 0;
      underlyingMatchesLogical_ = true;
      for (ordinal_type i=0; i<7; i++)
      {
        if (variationType_[i] == GENERAL)
        {
          if (extents_[i] != 0)
          {
            variationModulus_[i] = extents_[i];
          }
          else
          {
            variationModulus_[i] = 1;
          }
          activeDims_[numActiveDims_] = i;
          numActiveDims_++;
        }
        else if (variationType_[i] == MODULAR)
        {
          underlyingMatchesLogical_ = false;
          if (extents_[i] != getUnderlyingViewExtent(numActiveDims_))
          {
            const int dataExtent = getUnderlyingViewExtent(numActiveDims_);
            const int logicalExtent = extents_[i];
            const int modulus = dataExtent;
            
            INTREPID2_TEST_FOR_EXCEPTION( dataExtent * (logicalExtent / dataExtent) != logicalExtent, std::invalid_argument, "data extent must evenly divide logical extent");
            
            variationModulus_[i] = modulus;
          }
          else
          {
            variationModulus_[i] = extents_[i];
          }
          activeDims_[numActiveDims_] = i;
          numActiveDims_++;
        }
        else if (variationType_[i] == BLOCK_PLUS_DIAGONAL)
        {
          underlyingMatchesLogical_ = false;
          blockPlusDiagonalCount++;
          if (blockPlusDiagonalCount == 1) // first dimension thus marked --> active
          {
            
#ifdef HAVE_INTREPID2_DEBUG
            const int numNondiagonalEntries = blockPlusDiagonalNumNondiagonalEntries(blockPlusDiagonalLastNonDiagonal_);
            const int dataExtent = getUnderlyingViewExtent(numActiveDims_); // flat storage of all matrix entries
            const int logicalExtent = extents_[i];
            const int numDiagonalEntries    = logicalExtent - (blockPlusDiagonalLastNonDiagonal_ + 1);
            const int expectedDataExtent = numNondiagonalEntries + numDiagonalEntries;
            INTREPID2_TEST_FOR_EXCEPTION(dataExtent != expectedDataExtent, std::invalid_argument, ("BLOCK_PLUS_DIAGONAL data extent of " + std::to_string(dataExtent) + " does not match expected based on blockPlusDiagonalLastNonDiagonal setting of " + std::to_string(blockPlusDiagonalLastNonDiagonal_)).c_str());
#endif
            
            activeDims_[numActiveDims_] = i;
            numActiveDims_++;
          }
          variationModulus_[i] = getUnderlyingViewExtent(numActiveDims_);
          INTREPID2_TEST_FOR_EXCEPTION(variationType_[i+1] != BLOCK_PLUS_DIAGONAL, std::invalid_argument, "BLOCK_PLUS_DIAGONAL ranks must be contiguous");
          i++; // skip over the next BLOCK_PLUS_DIAGONAL
          variationModulus_[i] = 1; // trivial modulus (should not ever be used)
          INTREPID2_TEST_FOR_EXCEPTION(blockPlusDiagonalCount > 1, std::invalid_argument, "BLOCK_PLUS_DIAGONAL can only apply to two ranks");
        }
        else // CONSTANT
        {
          if (i < rank_)
          {
            underlyingMatchesLogical_ = false;
          }
          variationModulus_[i] = 1; // trivial modulus
        }
      }
      
      if (rank_ != dataRank_)
      {
        underlyingMatchesLogical_ = false;
      }
      
      for (int d=numActiveDims_; d<7; d++)
      {
        // for *inactive* dims, the activeDims_ map just is the identity
        // (this allows getEntry() to work even when the logical rank of the Data object is lower than that of the underlying View.  This can happen for gradients in 1D.)
        activeDims_[d] = d;
      }
      for (int d=0; d<7; d++)
      {
        INTREPID2_TEST_FOR_EXCEPTION(variationModulus_[d] == 0, std::logic_error, "variationModulus should not ever be 0");
      }
    }

    public:
    //! applies the specified unary operator to each entry
    template<class UnaryOperator>
    void applyOperator(UnaryOperator unaryOperator)
    {
      using ExecutionSpace = typename DeviceType::execution_space;
      
      switch (dataRank_)
      {
        case 1:
        {
          const int dataRank = 1;
          auto view = getUnderlyingView<dataRank>();
          
          const int dataExtent = this->getDataExtent(0);
          Kokkos::RangePolicy<ExecutionSpace> policy(ExecutionSpace(),0,dataExtent);
          Kokkos::parallel_for("apply operator in-place", policy,
          KOKKOS_LAMBDA (const int &i0) {
            view(i0) = unaryOperator(view(i0));
          });
          
        }
        break;
        case 2:
        {
          const int dataRank = 2;
          auto policy = dataExtentRangePolicy<dataRank>();
          auto view = getUnderlyingView<dataRank>();
          
          Kokkos::parallel_for("apply operator in-place", policy,
          KOKKOS_LAMBDA (const int &i0, const int &i1) {
            view(i0,i1) = unaryOperator(view(i0,i1));
          });
        }
        break;
        case 3:
        {
          const int dataRank = 3;
          auto policy = dataExtentRangePolicy<dataRank>();
          auto view = getUnderlyingView<dataRank>();
          
          Kokkos::parallel_for("apply operator in-place", policy,
          KOKKOS_LAMBDA (const int &i0, const int &i1, const int &i2) {
            view(i0,i1,i2) = unaryOperator(view(i0,i1,i2));
          });
        }
        break;
        case 4:
        {
          const int dataRank = 4;
          auto policy = dataExtentRangePolicy<dataRank>();
          auto view = getUnderlyingView<dataRank>();
          
          Kokkos::parallel_for("apply operator in-place", policy,
          KOKKOS_LAMBDA (const int &i0, const int &i1, const int &i2, const int &i3) {
            view(i0,i1,i2,i3) = unaryOperator(view(i0,i1,i2,i3));
          });
        }
        break;
        case 5:
        {
          const int dataRank = 5;
          auto policy = dataExtentRangePolicy<dataRank>();
          auto view = getUnderlyingView<dataRank>();
          
          Kokkos::parallel_for("apply operator in-place", policy,
          KOKKOS_LAMBDA (const int &i0, const int &i1, const int &i2, const int &i3, const int &i4) {
            view(i0,i1,i2,i3,i4) = unaryOperator(view(i0,i1,i2,i3,i4));
          });
        }
        break;
        case 6:
        {
          const int dataRank = 6;
          auto policy = dataExtentRangePolicy<dataRank>();
          auto view = getUnderlyingView<dataRank>();
          
          Kokkos::parallel_for("apply operator in-place", policy,
          KOKKOS_LAMBDA (const int &i0, const int &i1, const int &i2, const int &i3, const int &i4, const int &i5) {
            view(i0,i1,i2,i3,i4,i5) = unaryOperator(view(i0,i1,i2,i3,i4,i5));
          });
        }
        break;
        case 7:
        {
          const int dataRank = 7;
          auto policy6 = dataExtentRangePolicy<6>();
          auto view = getUnderlyingView<dataRank>();
          
          const int dim_i6 = view.extent_int(6);
          
          Kokkos::parallel_for("apply operator in-place", policy6,
          KOKKOS_LAMBDA (const int &i0, const int &i1, const int &i2, const int &i3, const int &i4, const int &i5) {
            for (int i6=0; i6<dim_i6; i6++)
            {
              view(i0,i1,i2,i3,i4,i5,i6) = unaryOperator(view(i0,i1,i2,i3,i4,i5,i6));
            }
          });
        }
        break;
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Unsupported data rank");
      }
    }

    //! Returns an l-value reference to the specified logical entry in the underlying view.  Note that for variation types other than GENERAL, multiple valid argument sets will refer to the same memory location.  Intended for Intrepid2 developers and expert users only.  If passThroughBlockDiagonalArgs is TRUE, the corresponding arguments are interpreted as entries in the 1D packed matrix rather than as logical 2D matrix row and column.
    template<class ...IntArgs>
    KOKKOS_INLINE_FUNCTION
    reference_type getWritableEntryWithPassThroughOption(const bool &passThroughBlockDiagonalArgs, const IntArgs... intArgs) const
    {
  #ifdef INTREPID2_HAVE_DEBUG
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(numArgs != rank_, std::invalid_argument, "getWritableEntry() should have the same number of arguments as the logical rank.");
  #endif
        constexpr int numArgs = sizeof...(intArgs);
        if (underlyingMatchesLogical_)
        {
          // in this case, we require that numArgs == dataRank_
          return getUnderlyingView<numArgs>()(intArgs...);
        }
        
        // extract the type of the first argument; use that for the arrays below
        using int_type = std::tuple_element_t<0, std::tuple<IntArgs...>>;
        
        const Kokkos::Array<int_type, numArgs+1> args {intArgs...,0}; // we pad with one extra entry (0) to avoid gcc compiler warnings about references beyond the bounds of the array (the [d+1]'s below)
        Kokkos::Array<int_type, 7> refEntry;
        for (int d=0; d<numArgs; d++)
        {
          switch (variationType_[d])
          {
            case CONSTANT: refEntry[d] = 0;                              break;
            case GENERAL:  refEntry[d] = args[d];                        break;
            case MODULAR:  refEntry[d] = args[d] % variationModulus_[d]; break;
            case BLOCK_PLUS_DIAGONAL:
            {
              if (passThroughBlockDiagonalArgs)
              {
                refEntry[d]   = args[d];
                refEntry[d+1] = args[d+1]; // this had better be == 0
                INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(args[d+1] != 0, std::invalid_argument, "getWritableEntry() called with passThroughBlockDiagonalArgs = true, but nonzero second matrix argument.");
              }
              else
              {
                const int numNondiagonalEntries = blockPlusDiagonalNumNondiagonalEntries(blockPlusDiagonalLastNonDiagonal_);
                
                const int_type &i = args[d];
                if (d+1 >= numArgs)
                {
                  INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "BLOCK_PLUS_DIAGONAL must be present for two dimensions; here, encountered only one");
                }
                else
                {
                  const int_type &j = args[d+1];
                  
                  if ((i > static_cast<int_type>(blockPlusDiagonalLastNonDiagonal_)) || (j > static_cast<int_type>(blockPlusDiagonalLastNonDiagonal_)))
                  {
                    if (i != j)
                    {
                      // off diagonal: zero
                      return zeroView_(0); // NOTE: this branches in an argument-dependent way; this is not great for CUDA performance.  When using BLOCK_PLUS_DIAGONAL, should generally avoid calls to this getEntry() method.  (Use methods that directly take advantage of the data packing instead.)
                    }
                    else
                    {
                      refEntry[d] = blockPlusDiagonalDiagonalEntryIndex(blockPlusDiagonalLastNonDiagonal_, numNondiagonalEntries, i);
                    }
                  }
                  else
                  {
                    refEntry[d] = blockPlusDiagonalBlockEntryIndex(blockPlusDiagonalLastNonDiagonal_, numNondiagonalEntries, i, j);
                  }

                  // skip next d (this is required also to be BLOCK_PLUS_DIAGONAL, and we've consumed its arg as j above)
                  refEntry[d+1] = 0;
                }
              }
              d++;
            }
          }
        }
         // refEntry should be zero-filled beyond numArgs, for cases when rank_ < dataRank_ (this only is allowed if the extra dimensions each has extent 1).
        for (int d=numArgs; d<7; d++)
        {
          refEntry[d] = 0;
        }
        
        if (dataRank_ == 1)
        {
          return data1_(refEntry[activeDims_[0]]);
        }
        else if (dataRank_ == 2)
        {
          return data2_(refEntry[activeDims_[0]],refEntry[activeDims_[1]]);
        }
        else if (dataRank_ == 3)
        {
          return data3_(refEntry[activeDims_[0]],refEntry[activeDims_[1]],refEntry[activeDims_[2]]);
        }
        else if (dataRank_ == 4)
        {
          return data4_(refEntry[activeDims_[0]],refEntry[activeDims_[1]],refEntry[activeDims_[2]],refEntry[activeDims_[3]]);
        }
        else if (dataRank_ == 5)
        {
          return data5_(refEntry[activeDims_[0]],refEntry[activeDims_[1]],refEntry[activeDims_[2]],refEntry[activeDims_[3]],
                        refEntry[activeDims_[4]]);
        }
        else if (dataRank_ == 6)
        {
          return data6_(refEntry[activeDims_[0]],refEntry[activeDims_[1]],refEntry[activeDims_[2]],refEntry[activeDims_[3]],
                        refEntry[activeDims_[4]],refEntry[activeDims_[5]]);
        }
        else // dataRank_ == 7
        {
          return data7_(refEntry[activeDims_[0]],refEntry[activeDims_[1]],refEntry[activeDims_[2]],refEntry[activeDims_[3]],
                        refEntry[activeDims_[4]],refEntry[activeDims_[5]],refEntry[activeDims_[6]]);
        }
      
    }
    
    //! Returns an l-value reference to the specified logical entry in the underlying view.  Note that for variation types other than GENERAL, multiple valid argument sets will refer to the same memory location.  Intended for Intrepid2 developers and expert users only.  If passThroughBlockDiagonalArgs is TRUE, the corresponding arguments are interpreted as entries in the 1D packed matrix rather than as logical 2D matrix row and column.
    template<class ...IntArgs>
    KOKKOS_INLINE_FUNCTION
    reference_type getWritableEntry(const IntArgs... intArgs) const
    {
      return getWritableEntryWithPassThroughOption(false, intArgs...);
    }
  public:
    //! Generic data copying method to allow construction of Data object from DynRankViews for which deep_copy() to the underlying view would be disallowed.  This method made public to allow CUDA compilation (because it contains a Kokkos lambda).
    template<class ToContainer, class FromContainer>
    static void copyContainer(ToContainer to, FromContainer from)
    {
//      std::cout << "Entered copyContainer().\n";
      auto policy = Kokkos::MDRangePolicy<execution_space,Kokkos::Rank<6>>({0,0,0,0,0,0},{from.extent_int(0),from.extent_int(1),from.extent_int(2), from.extent_int(3), from.extent_int(4), from.extent_int(5)});
      
      Kokkos::parallel_for("copyContainer", policy,
      KOKKOS_LAMBDA (const int &i0, const int &i1, const int &i2, const int &i3, const int &i4, const int &i5) {
        for (int i6=0; i6<from.extent_int(6); i6++)
        {
          to.access(i0,i1,i2,i3,i4,i5,i6) = from.access(i0,i1,i2,i3,i4,i5,i6);
        }
      });
    }
    
    //! allocate an underlying View that matches the provided DynRankView in dimensions, and copy.  Called by constructors that accept a DynRankView as argument.
    void allocateAndCopyFromDynRankView(ScalarView<DataScalar,DeviceType> data)
    {
//      std::cout << "Entered allocateAndCopyFromDynRankView().\n";
      switch (dataRank_)
      {
        case 1: data1_ = Kokkos::View<DataScalar*,       DeviceType>("Intrepid2 Data", data.extent_int(0)); break;
        case 2: data2_ = Kokkos::View<DataScalar**,      DeviceType>("Intrepid2 Data", data.extent_int(0), data.extent_int(1)); break;
        case 3: data3_ = Kokkos::View<DataScalar***,     DeviceType>("Intrepid2 Data", data.extent_int(0), data.extent_int(1), data.extent_int(2)); break;
        case 4: data4_ = Kokkos::View<DataScalar****,    DeviceType>("Intrepid2 Data", data.extent_int(0), data.extent_int(1), data.extent_int(2), data.extent_int(3)); break;
        case 5: data5_ = Kokkos::View<DataScalar*****,   DeviceType>("Intrepid2 Data", data.extent_int(0), data.extent_int(1), data.extent_int(2), data.extent_int(3), data.extent_int(4)); break;
        case 6: data6_ = Kokkos::View<DataScalar******,  DeviceType>("Intrepid2 Data", data.extent_int(0), data.extent_int(1), data.extent_int(2), data.extent_int(3), data.extent_int(4), data.extent_int(5)); break;
        case 7: data7_ = Kokkos::View<DataScalar*******, DeviceType>("Intrepid2 Data", data.extent_int(0), data.extent_int(1), data.extent_int(2), data.extent_int(3), data.extent_int(4), data.extent_int(5), data.extent_int(6)); break;
        default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
      }
      
      switch (dataRank_)
      {
        case 1: copyContainer(data1_,data); break;
        case 2: copyContainer(data2_,data); break;
        case 3: copyContainer(data3_,data); break;
        case 4: copyContainer(data4_,data); break;
        case 5: copyContainer(data5_,data); break;
        case 6: copyContainer(data6_,data); break;
        case 7: copyContainer(data7_,data); break;
        default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
      }
    }
    
    //! Constructor in terms of DimensionInfo for each logical dimension; does not require a View to be specified.  Will allocate a View of appropriate rank, zero-filled.
    Data(std::vector<DimensionInfo> dimInfoVector)
    :
    // initialize member variables as if default constructor; if dimInfoVector is empty, we want default constructor behavior.
    dataRank_(0), extents_({0,0,0,0,0,0,0}), variationType_({CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT}), blockPlusDiagonalLastNonDiagonal_(-1), rank_(dimInfoVector.size())
    {
      // If dimInfoVector is empty, the member initialization above is correct; otherwise, we set as below.
      // Either way, once members are initialized, we must call setActiveDims().
      if (dimInfoVector.size() != 0)
      {
        std::vector<int> dataExtents;

        bool blockPlusDiagonalEncountered = false;
        for (int d=0; d<rank_; d++)
        {
          const DimensionInfo & dimInfo = dimInfoVector[d];
          extents_[d] = dimInfo.logicalExtent;
          variationType_[d] = dimInfo.variationType;
          const bool isBlockPlusDiagonal = (variationType_[d] == BLOCK_PLUS_DIAGONAL);
          const bool isSecondBlockPlusDiagonal = isBlockPlusDiagonal && blockPlusDiagonalEncountered;
          if (isBlockPlusDiagonal)
          {
            blockPlusDiagonalEncountered = true;
            blockPlusDiagonalLastNonDiagonal_ = dimInfo.blockPlusDiagonalLastNonDiagonal;
          }
          if ((variationType_[d] != CONSTANT) && (!isSecondBlockPlusDiagonal))
          {
            dataExtents.push_back(dimInfo.dataExtent);
          }
        }
        if (dataExtents.size() == 0)
        {
          // constant data
          dataExtents.push_back(1);
        }
        dataRank_ = dataExtents.size();
        switch (dataRank_)
        {
          case 1: data1_ = Kokkos::View<DataScalar*,       DeviceType>("Intrepid2 Data", dataExtents[0]); break;
          case 2: data2_ = Kokkos::View<DataScalar**,      DeviceType>("Intrepid2 Data", dataExtents[0], dataExtents[1]); break;
          case 3: data3_ = Kokkos::View<DataScalar***,     DeviceType>("Intrepid2 Data", dataExtents[0], dataExtents[1], dataExtents[2]); break;
          case 4: data4_ = Kokkos::View<DataScalar****,    DeviceType>("Intrepid2 Data", dataExtents[0], dataExtents[1], dataExtents[2], dataExtents[3]); break;
          case 5: data5_ = Kokkos::View<DataScalar*****,   DeviceType>("Intrepid2 Data", dataExtents[0], dataExtents[1], dataExtents[2], dataExtents[3], dataExtents[4]); break;
          case 6: data6_ = Kokkos::View<DataScalar******,  DeviceType>("Intrepid2 Data", dataExtents[0], dataExtents[1], dataExtents[2], dataExtents[3], dataExtents[4], dataExtents[5]); break;
          case 7: data7_ = Kokkos::View<DataScalar*******, DeviceType>("Intrepid2 Data", dataExtents[0], dataExtents[1], dataExtents[2], dataExtents[3], dataExtents[4], dataExtents[5], dataExtents[6]); break;
          default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
        }
      }
      setActiveDims();
    }
    
    //! DynRankView constructor.  Will copy to a View of appropriate rank.
    Data(const ScalarView<DataScalar,DeviceType> &data, int rank, Kokkos::Array<int,7> extents, Kokkos::Array<DataVariationType,7> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
    :
    dataRank_(data.rank()), extents_(extents), variationType_(variationType), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
    {
      allocateAndCopyFromDynRankView(data);
      setActiveDims();
    }
    
    //! copy-like constructor for differing device type, but same memory space.  This does a shallow copy of the underlying view.
    template<typename OtherDeviceType, class = typename std::enable_if< std::is_same<typename DeviceType::memory_space, typename OtherDeviceType::memory_space>::value>::type,
                                       class = typename std::enable_if<!std::is_same<DeviceType,OtherDeviceType>::value>::type>
    Data(const Data<DataScalar,OtherDeviceType> &data)
    :
    dataRank_(data.getUnderlyingViewRank()), extents_(data.getExtents()), variationType_(data.getVariationTypes()), blockPlusDiagonalLastNonDiagonal_(data.blockPlusDiagonalLastNonDiagonal()), rank_(data.rank())
    {
//      std::cout << "Entered copy-like Data constructor.\n";
      if (dataRank_ != 0) // dataRank_ == 0 indicates an invalid Data object (a placeholder, can indicate zero value)
      {
        const auto view = data.getUnderlyingView();
        switch (dataRank_)
        {
          case 1: data1_ = data.getUnderlyingView1(); break;
          case 2: data2_ = data.getUnderlyingView2(); break;
          case 3: data3_ = data.getUnderlyingView3(); break;
          case 4: data4_ = data.getUnderlyingView4(); break;
          case 5: data5_ = data.getUnderlyingView5(); break;
          case 6: data6_ = data.getUnderlyingView6(); break;
          case 7: data7_ = data.getUnderlyingView7(); break;
          default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
        }
      }
      setActiveDims();
    }
    
    //! copy-like constructor for differing execution spaces.  This does a deep_copy of the underlying view.
    template<typename OtherDeviceType, class = typename std::enable_if<!std::is_same<typename DeviceType::memory_space, typename OtherDeviceType::memory_space>::value>::type>
    Data(const Data<DataScalar,OtherDeviceType> &data)
    :
    dataRank_(data.getUnderlyingViewRank()), extents_(data.getExtents()), variationType_(data.getVariationTypes()), blockPlusDiagonalLastNonDiagonal_(data.blockPlusDiagonalLastNonDiagonal()), rank_(data.rank())
    {
//      std::cout << "Entered copy-like Data constructor.\n";
      if (dataRank_ != 0) // dataRank_ == 0 indicates an invalid Data object (a placeholder, can indicate zero value)
      {
        const auto view = data.getUnderlyingView();
        switch (dataRank_)
        {
          case 1: data1_ = Kokkos::View<DataScalar*,       DeviceType>("Intrepid2 Data", view.extent_int(0)); break;
          case 2: data2_ = Kokkos::View<DataScalar**,      DeviceType>("Intrepid2 Data", view.extent_int(0), view.extent_int(1)); break;
          case 3: data3_ = Kokkos::View<DataScalar***,     DeviceType>("Intrepid2 Data", view.extent_int(0), view.extent_int(1), view.extent_int(2)); break;
          case 4: data4_ = Kokkos::View<DataScalar****,    DeviceType>("Intrepid2 Data", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3)); break;
          case 5: data5_ = Kokkos::View<DataScalar*****,   DeviceType>("Intrepid2 Data", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3), view.extent_int(4)); break;
          case 6: data6_ = Kokkos::View<DataScalar******,  DeviceType>("Intrepid2 Data", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3), view.extent_int(4), view.extent_int(5)); break;
          case 7: data7_ = Kokkos::View<DataScalar*******, DeviceType>("Intrepid2 Data", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3), view.extent_int(4), view.extent_int(5), view.extent_int(6)); break;
          default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
        }
        
        // copy
        // (Note: Kokkos::deep_copy() will not generally work if the layouts are different; that's why we do a manual copy here once we have the data on the host):
        // first, mirror and copy dataView; then copy to the appropriate data_ member
        using MemorySpace = typename DeviceType::memory_space;
        switch (dataRank_)
        {
          case 1: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView1()); copyContainer(data1_, dataViewMirror);} break;
          case 2: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView2()); copyContainer(data2_, dataViewMirror);} break;
          case 3: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView3()); copyContainer(data3_, dataViewMirror);} break;
          case 4: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView4()); copyContainer(data4_, dataViewMirror);} break;
          case 5: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView5()); copyContainer(data5_, dataViewMirror);} break;
          case 6: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView6()); copyContainer(data6_, dataViewMirror);} break;
          case 7: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView7()); copyContainer(data7_, dataViewMirror);} break;
          default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
        }
      }
      setActiveDims();
    }
    
    //! copy constructor modeled after the copy-like constructor above.  Not as efficient as the implicit copy constructor, so this should only be uncommented to emulate the copy-like constructor above in situations where host and device space are the same, and the copy-like constructor does not exist, or to provide a debugging breakpoint to assess when copies are being constructed.
//    Data(const Data<DataScalar,ExecSpaceType> &data)
//    :
//    dataRank_(data.getUnderlyingViewRank()), extents_(data.getExtents()), variationType_(data.getVariationTypes()), blockPlusDiagonalLastNonDiagonal_(data.blockPlusDiagonalLastNonDiagonal()), rank_(data.rank())
//    {
//      std::cout << "Entered Data copy constructor.\n";
//      if (dataRank_ != 0) // dataRank_ == 0 indicates an invalid Data object (a placeholder, can indicate zero value)
//      {
//        const auto view = data.getUnderlyingView();
//        switch (dataRank_)
//        {
//          case 1: data1_ = Kokkos::View<DataScalar*,       DeviceType>("Intrepid2 Data - explicit copy constructor(for debugging)", view.extent_int(0)); break;
//          case 2: data2_ = Kokkos::View<DataScalar**,      DeviceType>("Intrepid2 Data - explicit copy constructor(for debugging)", view.extent_int(0), view.extent_int(1)); break;
//          case 3: data3_ = Kokkos::View<DataScalar***,     DeviceType>("Intrepid2 Data - explicit copy constructor(for debugging)", view.extent_int(0), view.extent_int(1), view.extent_int(2)); break;
//          case 4: data4_ = Kokkos::View<DataScalar****,    DeviceType>("Intrepid2 Data - explicit copy constructor(for debugging)", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3)); break;
//          case 5: data5_ = Kokkos::View<DataScalar*****,   DeviceType>("Intrepid2 Data - explicit copy constructor(for debugging)", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3), view.extent_int(4)); break;
//          case 6: data6_ = Kokkos::View<DataScalar******,  DeviceType>("Intrepid2 Data - explicit copy constructor(for debugging)", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3), view.extent_int(4), view.extent_int(5)); break;
//          case 7: data7_ = Kokkos::View<DataScalar*******, DeviceType>("Intrepid2 Data - explicit copy constructor(for debugging)", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3), view.extent_int(4), view.extent_int(5), view.extent_int(6)); break;
//          default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
//        }
//
//        // copy
//        // (Note: Kokkos::deep_copy() will not generally work if the layouts are different; that's why we do a manual copy here once we have the data on the host):
//        // first, mirror and copy dataView; then copy to the appropriate data_ member
//        using MemorySpace = typename DeviceType::memory_space;
//        switch (dataRank_)
//        {
//          case 1: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView1()); copyContainer(data1_, dataViewMirror);} break;
//          case 2: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView2()); copyContainer(data2_, dataViewMirror);} break;
//          case 3: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView3()); copyContainer(data3_, dataViewMirror);} break;
//          case 4: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView4()); copyContainer(data4_, dataViewMirror);} break;
//          case 5: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView5()); copyContainer(data5_, dataViewMirror);} break;
//          case 6: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView6()); copyContainer(data6_, dataViewMirror);} break;
//          case 7: {auto dataViewMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), data.getUnderlyingView7()); copyContainer(data7_, dataViewMirror);} break;
//          default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
//        }
//      }
//
//      setActiveDims();
//    }
    
    //! constructor for fully varying data container, no expressed redundancy/repetition.  Copies the data to a new Kokkos::View of matching rank.
    Data(ScalarView<DataScalar,DeviceType> data)
    :
    Data(data,
         data.rank(),
         Kokkos::Array<int,7> {data.extent_int(0),data.extent_int(1),data.extent_int(2),data.extent_int(3),data.extent_int(4),data.extent_int(5),data.extent_int(6)},
         Kokkos::Array<DataVariationType,7> {GENERAL,GENERAL,GENERAL,GENERAL,GENERAL,GENERAL,GENERAL}, -1)
    {}
    
    //! Constructor that accepts a DynRankView as an argument.  The data belonging to the DynRankView will be copied into a new View of matching dimensions.
    template<size_t rank, class ...DynRankViewProperties>
    Data(const Kokkos::DynRankView<DataScalar,DeviceType, DynRankViewProperties...> &data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
    :
    dataRank_(data.rank()), extents_({1,1,1,1,1,1,1}), variationType_({CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
    {
//      std::cout << "Entered a DynRankView Data() constructor.\n";
      allocateAndCopyFromDynRankView(data);
      for (unsigned d=0; d<rank; d++)
      {
        extents_[d]       = extents[d];
        variationType_[d] = variationType[d];
      }
      setActiveDims();
    }
        
    template<size_t rank, class ...ViewProperties>
    Data(Kokkos::View<DataScalar*,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
    :
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
    {
      data1_ = data;
      for (unsigned d=0; d<rank; d++)
      {
        extents_[d]       = extents[d];
        variationType_[d] = variationType[d];
      }
      setActiveDims();
    }
    
    template<size_t rank, class ...ViewProperties>
    Data(Kokkos::View<DataScalar**,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
    :
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
    {
      data2_ = data;
      for (unsigned d=0; d<rank; d++)
      {
        extents_[d]       = extents[d];
        variationType_[d] = variationType[d];
      }
      setActiveDims();
    }
    
    template<size_t rank, class ...ViewProperties>
    Data(Kokkos::View<DataScalar***,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
    :
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
    {
      data3_ = data;
      for (unsigned d=0; d<rank; d++)
      {
        extents_[d]       = extents[d];
        variationType_[d] = variationType[d];
      }
      setActiveDims();
    }
    
    template<size_t rank, class ...ViewProperties>
    Data(Kokkos::View<DataScalar****,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
    :
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
    {
      data4_ = data;
      for (unsigned d=0; d<rank; d++)
      {
        extents_[d]       = extents[d];
        variationType_[d] = variationType[d];
      }
      setActiveDims();
    }
    
    template<size_t rank, class ...ViewProperties>
    Data(Kokkos::View<DataScalar*****,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
    :
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
    {
      data5_ = data;
      for (unsigned d=0; d<rank; d++)
      {
        extents_[d]       = extents[d];
        variationType_[d] = variationType[d];
      }
      setActiveDims();
    }
    
    template<size_t rank, class ...ViewProperties>
    Data(Kokkos::View<DataScalar******,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
    :
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
    {
      data6_ = data;
      for (unsigned d=0; d<rank; d++)
      {
        extents_[d]       = extents[d];
        variationType_[d] = variationType[d];
      }
      setActiveDims();
    }
    
    template<size_t rank, class ...ViewProperties>
    Data(Kokkos::View<DataScalar*******,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
    :
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
    {
      data7_ = data;
      for (unsigned d=0; d<rank; d++)
      {
        extents_[d]       = extents[d];
        variationType_[d] = variationType[d];
      }
      setActiveDims();
    }
    
    //! constructor with run-time rank (requires full-length extents, variationType inputs; those beyond the rank will be ignored).
    template<class ViewScalar, class ...ViewProperties>
    Data(const unsigned rank, Kokkos::View<ViewScalar,DeviceType, ViewProperties...> data, Kokkos::Array<int,7> extents, Kokkos::Array<DataVariationType,7> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
    :
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
    {
      setUnderlyingView<data.rank>(data);
      for (unsigned d=0; d<rank; d++)
      {
        extents_[d]       = extents[d];
        variationType_[d] = variationType[d];
      }
      setActiveDims();
    }
    
    //! constructor for everywhere-constant data
    template<size_t rank>
    Data(DataScalar constantValue, Kokkos::Array<int,rank> extents)
    :
    dataRank_(1), extents_({1,1,1,1,1,1,1}), variationType_({CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT}), blockPlusDiagonalLastNonDiagonal_(-1), rank_(rank)
    {
      data1_ = Kokkos::View<DataScalar*,DeviceType>("Constant Data",1);
      Kokkos::deep_copy(data1_, constantValue);
      for (unsigned d=0; d<rank; d++)
      {
        extents_[d] = extents[d];
      }
      setActiveDims();
    }
            
    //! default constructor (empty data)
    Data()
    :
    dataRank_(0), extents_({0,0,0,0,0,0,0}), variationType_({CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT,CONSTANT}), blockPlusDiagonalLastNonDiagonal_(-1), rank_(0)
    {
      setActiveDims();
    }
    
    //! For a Data object containing data with variation type BLOCK_PLUS_DIAGONAL, returns the row and column (which must match) of the last non-diagonal entry.  For fully diagonal matrices, this is -1.
    KOKKOS_INLINE_FUNCTION
    const int & blockPlusDiagonalLastNonDiagonal() const
    {
      return blockPlusDiagonalLastNonDiagonal_;
    }
    
    //! Returns an array containing the logical extents in each dimension.
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<int,7> getExtents() const
    {
      return extents_;
    }
    
    //! Returns an object fully specifying the indicated dimension.  This is used in determining appropriate sizing of Data objects that depend on other Data objects (e.g., when two matrix Data objects are multiplied together).
    KOKKOS_INLINE_FUNCTION
    DimensionInfo getDimensionInfo(const int &dim) const
    {
      DimensionInfo dimInfo;
      
      dimInfo.logicalExtent = extent_int(dim);
      dimInfo.variationType = variationType_[dim];
      dimInfo.dataExtent    = getDataExtent(dim);
      dimInfo.variationModulus = variationModulus_[dim];
      
      if (dimInfo.variationType == BLOCK_PLUS_DIAGONAL)
      {
        dimInfo.blockPlusDiagonalLastNonDiagonal = blockPlusDiagonalLastNonDiagonal_;
      }
      return dimInfo;
    }
    
    //! Returns (DataVariationType, data extent) in the specified dimension for a Data container that combines (through multiplication, say, or addition) this container with otherData.
    KOKKOS_INLINE_FUNCTION
    DimensionInfo combinedDataDimensionInfo(const Data &otherData, const int &dim) const
    {
      const DimensionInfo myDimInfo    = getDimensionInfo(dim);
      const DimensionInfo otherDimInfo = otherData.getDimensionInfo(dim);
      
      return combinedDimensionInfo(myDimInfo, otherDimInfo);
    }
    
    //! Returns the underlying view.  Throws an exception if the underlying view is not rank 1.
    template<int rank>
    KOKKOS_INLINE_FUNCTION
    enable_if_t<rank==1, const Kokkos::View<typename RankExpander<DataScalar, rank>::value_type, DeviceType> &>
    getUnderlyingView() const
    {
      #ifdef HAVE_INTREPID2_DEBUG
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(dataRank_ != rank, std::invalid_argument, "getUnderlyingView() called for rank that does not match dataRank_");
      #endif
      return data1_;
    }
    
    //! Returns the underlying view.  Throws an exception if the underlying view is not rank 2.
    template<int rank>
    KOKKOS_INLINE_FUNCTION
    enable_if_t<rank==2, const Kokkos::View<typename RankExpander<DataScalar, rank>::value_type, DeviceType> &>
    getUnderlyingView() const
    {
      #ifdef HAVE_INTREPID2_DEBUG
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(dataRank_ != rank, std::invalid_argument, "getUnderlyingView() called for rank that does not match dataRank_");
      #endif
      return data2_;
    }
    
    //! Returns the underlying view.  Throws an exception if the underlying view is not rank 3.
    template<int rank>
    KOKKOS_INLINE_FUNCTION
    enable_if_t<rank==3, const Kokkos::View<typename RankExpander<DataScalar, rank>::value_type, DeviceType> &>
    getUnderlyingView() const
    {
      #ifdef HAVE_INTREPID2_DEBUG
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(dataRank_ != rank, std::invalid_argument, "getUnderlyingView() called for rank that does not match dataRank_");
      #endif
      return data3_;
    }
    
    //! Returns the underlying view.  Throws an exception if the underlying view is not rank 4.
    template<int rank>
    KOKKOS_INLINE_FUNCTION
    enable_if_t<rank==4, const Kokkos::View<typename RankExpander<DataScalar, rank>::value_type, DeviceType> &>
    getUnderlyingView() const
    {
      #ifdef HAVE_INTREPID2_DEBUG
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(dataRank_ != rank, std::invalid_argument, "getUnderlyingView() called for rank that does not match dataRank_");
      #endif
      return data4_;
    }
    
    //! Returns the underlying view.  Throws an exception if the underlying view is not rank 5.
    template<int rank>
    KOKKOS_INLINE_FUNCTION
    enable_if_t<rank==5, const Kokkos::View<typename RankExpander<DataScalar, rank>::value_type, DeviceType> &>
    getUnderlyingView() const
    {
      #ifdef HAVE_INTREPID2_DEBUG
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(dataRank_ != rank, std::invalid_argument, "getUnderlyingView() called for rank that does not match dataRank_");
      #endif
      return data5_;
    }
    
    //! Returns the underlying view.  Throws an exception if the underlying view is not rank 6.
    template<int rank>
    KOKKOS_INLINE_FUNCTION
    enable_if_t<rank==6, const Kokkos::View<typename RankExpander<DataScalar, rank>::value_type, DeviceType> &>
    getUnderlyingView() const
    {
      #ifdef HAVE_INTREPID2_DEBUG
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(dataRank_ != rank, std::invalid_argument, "getUnderlyingView() called for rank that does not match dataRank_");
      #endif
      return data6_;
    }
    
    //! Returns the underlying view.  Throws an exception if the underlying view is not rank 7.
    template<int rank>
    KOKKOS_INLINE_FUNCTION
    enable_if_t<rank==7, const Kokkos::View<typename RankExpander<DataScalar, rank>::value_type, DeviceType> &>
    getUnderlyingView() const
    {
      #ifdef HAVE_INTREPID2_DEBUG
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(dataRank_ != rank, std::invalid_argument, "getUnderlyingView() called for rank that does not match dataRank_");
      #endif
      return data7_;
    }
    
    //! returns the View that stores the unique data.  For rank-1 underlying containers.
    KOKKOS_INLINE_FUNCTION
    const Kokkos::View<DataScalar*, DeviceType> & getUnderlyingView1() const
    {
      return getUnderlyingView<1>();
    }
    
    //! returns the View that stores the unique data.  For rank-2 underlying containers.
    KOKKOS_INLINE_FUNCTION
    const Kokkos::View<DataScalar**, DeviceType> & getUnderlyingView2() const
    {
      return getUnderlyingView<2>();
    }
    
    //! returns the View that stores the unique data.  For rank-3 underlying containers.
    KOKKOS_INLINE_FUNCTION
    const Kokkos::View<DataScalar***, DeviceType> & getUnderlyingView3() const
    {
      return getUnderlyingView<3>();
    }
    
    //! returns the View that stores the unique data.  For rank-4 underlying containers.
    KOKKOS_INLINE_FUNCTION
    const Kokkos::View<DataScalar****, DeviceType> & getUnderlyingView4() const
    {
      return getUnderlyingView<4>();
    }
    
    //! returns the View that stores the unique data.  For rank-5 underlying containers.
    KOKKOS_INLINE_FUNCTION
    const Kokkos::View<DataScalar*****, DeviceType> & getUnderlyingView5() const
    {
      return getUnderlyingView<5>();
    }
    
    //! returns the View that stores the unique data.  For rank-6 underlying containers.
    KOKKOS_INLINE_FUNCTION
    const Kokkos::View<DataScalar******, DeviceType> & getUnderlyingView6() const
    {
      return getUnderlyingView<6>();
    }
    
    //! returns the View that stores the unique data.  For rank-7 underlying containers.
    KOKKOS_INLINE_FUNCTION
    const Kokkos::View<DataScalar*******, DeviceType> & getUnderlyingView7() const
    {
      return getUnderlyingView<7>();
    }
    
    //! sets the View that stores the unique data.  For rank-1 underlying containers.
    KOKKOS_INLINE_FUNCTION
    void setUnderlyingView1(const Kokkos::View<DataScalar*, DeviceType> & view)
    {
      data1_ = view;
    }
    
    //! sets the View that stores the unique data.  For rank-2 underlying containers.
    KOKKOS_INLINE_FUNCTION
    void setUnderlyingView2(const Kokkos::View<DataScalar**, DeviceType> & view)
    {
      data2_ = view;
    }
    
    //! sets the View that stores the unique data.  For rank-3 underlying containers.
    KOKKOS_INLINE_FUNCTION
    void setUnderlyingView3(const Kokkos::View<DataScalar***, DeviceType> & view)
    {
      data3_ = view;
    }
    
    //! sets the View that stores the unique data.  For rank-4 underlying containers.
    KOKKOS_INLINE_FUNCTION
    void setUnderlyingView4(const Kokkos::View<DataScalar****, DeviceType> & view)
    {
      data4_ = view;
    }
    
    //! sets the View that stores the unique data.  For rank-5 underlying containers.
    KOKKOS_INLINE_FUNCTION
    void setUnderlyingView5(const Kokkos::View<DataScalar*****, DeviceType> & view)
    {
      data5_ = view;
    }
    
    //! sets the View that stores the unique data.  For rank-6 underlying containers.
    KOKKOS_INLINE_FUNCTION
    void setUnderlyingView6(const Kokkos::View<DataScalar******, DeviceType> & view)
    {
      data6_ = view;
    }
    
    //! sets the View that stores the unique data.  For rank-7 underlying containers.
    KOKKOS_INLINE_FUNCTION
    void setUnderlyingView7(const Kokkos::View<DataScalar*******, DeviceType> & view)
    {
      data7_ = view;
    }
    
    template<int underlyingRank, class ViewScalar>
    KOKKOS_INLINE_FUNCTION
    void setUnderlyingView(const Kokkos::View<ViewScalar, DeviceType> & view)
    {
      if constexpr (underlyingRank == 1)
      {
        setUnderlyingView1(view);
      }
      else if constexpr (underlyingRank == 2)
      {
        setUnderlyingView2(view);
      }
      else if constexpr (underlyingRank == 3)
      {
        setUnderlyingView3(view);
      }
      else if constexpr (underlyingRank == 4)
      {
        setUnderlyingView4(view);
      }
      else if constexpr (underlyingRank == 5)
      {
        setUnderlyingView5(view);
      }
      else if constexpr (underlyingRank == 6)
      {
        setUnderlyingView6(view);
      }
      else if constexpr (underlyingRank == 7)
      {
        setUnderlyingView7(view);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "implementation for specialization missing");
      }
    }
    
    //! Returns a DynRankView constructed atop the same underlying data as the fixed-rank Kokkos::View used internally.
    ScalarView<DataScalar,DeviceType> getUnderlyingView() const
    {
      switch (dataRank_)
      {
        case 1: return data1_;
        case 2: return data2_;
        case 3: return data3_;
        case 4: return data4_;
        case 5: return data5_;
        case 6: return data6_;
        case 7: return data7_;
        default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
      }
    }
    
    //! returns the rank of the View that stores the unique data
    KOKKOS_INLINE_FUNCTION
    ordinal_type getUnderlyingViewRank() const
    {
      return dataRank_;
    }
    
    //! returns the number of entries in the View that stores the unique data
    KOKKOS_INLINE_FUNCTION
    ordinal_type getUnderlyingViewSize() const
    {
      ordinal_type size = 1;
      for (ordinal_type r=0; r<dataRank_; r++)
      {
        size *= getUnderlyingViewExtent(r);
      }
      return size;
    }
    
    //! Returns a DynRankView that matches the underlying Kokkos::View object in value_type, layout, and dimension.
    ScalarView<DataScalar,DeviceType> allocateDynRankViewMatchingUnderlying() const
    {
      switch (dataRank_)
      {
        case 1: return getMatchingViewWithLabel(data1_, "Intrepid2 Data", data1_.extent_int(0));
        case 2: return getMatchingViewWithLabel(data2_, "Intrepid2 Data", data2_.extent_int(0), data2_.extent_int(1));
        case 3: return getMatchingViewWithLabel(data3_, "Intrepid2 Data", data3_.extent_int(0), data3_.extent_int(1), data3_.extent_int(2));
        case 4: return getMatchingViewWithLabel(data4_, "Intrepid2 Data", data4_.extent_int(0), data4_.extent_int(1), data4_.extent_int(2), data4_.extent_int(3));
        case 5: return getMatchingViewWithLabel(data5_, "Intrepid2 Data", data5_.extent_int(0), data5_.extent_int(1), data5_.extent_int(2), data5_.extent_int(3), data5_.extent_int(4));
        case 6: return getMatchingViewWithLabel(data6_, "Intrepid2 Data", data6_.extent_int(0), data6_.extent_int(1), data6_.extent_int(2), data6_.extent_int(3), data6_.extent_int(4), data6_.extent_int(5));
        case 7: return getMatchingViewWithLabel(data7_, "Intrepid2 Data", data7_.extent_int(0), data7_.extent_int(1), data7_.extent_int(2), data7_.extent_int(3), data7_.extent_int(4), data7_.extent_int(5), data7_.extent_int(6));
        default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
      }
    }
    
    //! Returns a DynRankView that matches the underlying Kokkos::View object value_type and layout, but with the specified dimensions.
    template<class ... DimArgs>
    ScalarView<DataScalar,DeviceType> allocateDynRankViewMatchingUnderlying(DimArgs... dims) const
    {
      switch (dataRank_)
      {
        case 1: return getMatchingViewWithLabel(data1_, "Intrepid2 Data", dims...);
        case 2: return getMatchingViewWithLabel(data2_, "Intrepid2 Data", dims...);
        case 3: return getMatchingViewWithLabel(data3_, "Intrepid2 Data", dims...);
        case 4: return getMatchingViewWithLabel(data4_, "Intrepid2 Data", dims...);
        case 5: return getMatchingViewWithLabel(data5_, "Intrepid2 Data", dims...);
        case 6: return getMatchingViewWithLabel(data6_, "Intrepid2 Data", dims...);
        case 7: return getMatchingViewWithLabel(data7_, "Intrepid2 Data", dims...);
        default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
      }
    }
  
    //! Copies 0.0 to the underlying View.
    void clear() const
    {
      switch (dataRank_)
      {
        case 1: Kokkos::deep_copy(data1_, 0.0); break;
        case 2: Kokkos::deep_copy(data2_, 0.0); break;
        case 3: Kokkos::deep_copy(data3_, 0.0); break;
        case 4: Kokkos::deep_copy(data4_, 0.0); break;
        case 5: Kokkos::deep_copy(data5_, 0.0); break;
        case 6: Kokkos::deep_copy(data6_, 0.0); break;
        case 7: Kokkos::deep_copy(data7_, 0.0); break;
        default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
      }
    }
    
    //! Copies from the provided DynRankView into the underlying Kokkos::View container storing the unique data.
    void copyDataFromDynRankViewMatchingUnderlying(const ScalarView<DataScalar,DeviceType> &dynRankView) const
    {
//      std::cout << "Entered copyDataFromDynRankViewMatchingUnderlying().\n";
      switch (dataRank_)
      {
        case 1: copyContainer(data1_,dynRankView); break;
        case 2: copyContainer(data2_,dynRankView); break;
        case 3: copyContainer(data3_,dynRankView); break;
        case 4: copyContainer(data4_,dynRankView); break;
        case 5: copyContainer(data5_,dynRankView); break;
        case 6: copyContainer(data6_,dynRankView); break;
        case 7: copyContainer(data7_,dynRankView); break;
        default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
      }
    }
    
    //! returns the true extent of the data corresponding to the logical dimension provided; if the data does not vary in that dimension, returns 1
    KOKKOS_INLINE_FUNCTION int getDataExtent(const ordinal_type &d) const
    {
      for (int i=0; i<numActiveDims_; i++)
      {
        if (activeDims_[i] == d)
        {
          return getUnderlyingViewExtent(i);
        }
        else if (activeDims_[i] > d)
        {
          return 1; // data does not vary in the specified dimension
        }
      }
      return 1; // data does not vary in the specified dimension
    }
    
    /** \brief  Variation modulus accessor.
       \param [in] d - the logical dimension whose variation modulus is requested.
       \return the variation modulus.
     
     The variation modulus is defined as the number of unique entries in the specified dimension.
     This is defined as follows:
     - for CONSTANT variation, the variation modulus is 1
     - for MODULAR variation, the variation modulus is exactly the modulus by which the data repeats in the specified dimension
     - for GENERAL variation, the variation modulus is the extent in the specified dimension
     - for BLOCK_PLUS_DIAGONAL, the variation modulus in the first logical dimension of the matrix is the number of nonzeros in the matrix; in the second logical dimension the variation modulus is 1.
    */
    KOKKOS_INLINE_FUNCTION
    int getVariationModulus(const int &d) const
    {
      return variationModulus_[d];
    }
    
    //! Returns an array with the variation types in each logical dimension.
    KOKKOS_INLINE_FUNCTION
    const Kokkos::Array<DataVariationType,7> & getVariationTypes() const
    {
      return variationType_;
    }
    
    //! Returns a (read-only) value corresponding to the specified logical data location.  If passThroughBlockDiagonalArgs is TRUE, the corresponding arguments are interpreted as entries in the 1D packed matrix rather than as logical 2D matrix row and column.
    template<class ...IntArgs>
    KOKKOS_INLINE_FUNCTION
    return_type getEntryWithPassThroughOption(const bool &passThroughBlockDiagonalArgs, const IntArgs&... intArgs) const
    {
      return getWritableEntryWithPassThroughOption(passThroughBlockDiagonalArgs, intArgs...);
    }
    
    //! Returns a (read-only) value corresponding to the specified logical data location.
    template<class ...IntArgs>
    KOKKOS_INLINE_FUNCTION
    return_type getEntry(const IntArgs&... intArgs) const
    {
      return getEntryWithPassThroughOption(false, intArgs...);
    }
    
    template <bool...> struct bool_pack;

    template <bool... v>
    using all_true = std::is_same<bool_pack<true, v...>, bool_pack<v..., true>>;
    
    template <class ...IntArgs>
    using valid_args = all_true<std::is_integral<IntArgs>{}...>;
    
    static_assert(valid_args<int,long,unsigned>::value, "valid args works");

    //! Returns a value corresponding to the specified logical data location.
    template <class ...IntArgs>
    KOKKOS_INLINE_FUNCTION
#ifndef __INTEL_COMPILER
    // icc has a bug that prevents compilation with this enable_if_t
    // (possibly the same as https://community.intel.com/t5/Intel-C-Compiler/Intel-Compiler-bug-while-deducing-template-arguments-inside/m-p/1164358)
    // so with icc we'll just skip the argument type/count check
    enable_if_t<valid_args<IntArgs...>::value && (sizeof...(IntArgs) <= 7),return_type>
#else
    return_type
#endif
    operator()(const IntArgs&... intArgs) const {
      return getEntry(intArgs...);
    }

    //! Returns the logical extent in the specified dimension.
    KOKKOS_INLINE_FUNCTION
    int extent_int(const int& r) const
    {
      return extents_[r];
    }
    
    template <typename iType>
    KOKKOS_INLINE_FUNCTION constexpr
    typename std::enable_if<std::is_integral<iType>::value, size_t>::type
    extent(const iType& r) const {
      return extents_[r];
    }
    
    //! returns true for containers that have two dimensions marked as BLOCK_PLUS_DIAGONAL for which the non-diagonal block is empty or size 1.
    KOKKOS_INLINE_FUNCTION bool isDiagonal() const
    {
      if (blockPlusDiagonalLastNonDiagonal_ >= 1) return false;
      int numBlockPlusDiagonalTypes = 0;
      for (unsigned r = 0; r<variationType_.size(); r++)
      {
        const auto &entryType = variationType_[r];
        if (entryType == BLOCK_PLUS_DIAGONAL) numBlockPlusDiagonalTypes++;
      }
      // 2 BLOCK_PLUS_DIAGONAL entries, combined with blockPlusDiagonalLastNonDiagonal being -1 or 0 indicates diagonal
      if      (numBlockPlusDiagonalTypes == 2) return true;
      else if (numBlockPlusDiagonalTypes == 0) return false; // no BLOCK_PLUS_DIAGONAL --> not a diagonal matrix
      else INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "Unexpected number of ranks marked as BLOCK_PLUS_DIAGONAL (should be 0 or 2)");
      return false; // statement should be unreachable; included because compilers don't necessarily recognize that fact...
    }
    
    //! Constructs a container with extents matching this, with a single-value underlying View, CONSTANT in each dimension.
    //! \param value  [in] - the constant value.
    //! \return A container with the same logical shape as this, with a single-value underlying View.
    Data<DataScalar,DeviceType> allocateConstantData( const DataScalar &value )
    {
      return Data<DataScalar,DeviceType>(value, this->getExtents());
    }
    
    //! Constructs a container suitable for storing the result of an in-place combination of the two provided data containers.  The two containers must have the same logical shape.
    //! \see storeInPlaceCombination()
    //! \param A  [in] - the first data container.
    //! \param B  [in] - the second data container.  Must have the same logical shape as A.
    //! \return A container with the same logical shape as A and B, with underlying View storage sufficient to store the result of A + B (or any other in-place combination).
    static Data<DataScalar,DeviceType> allocateInPlaceCombinationResult( const Data<DataScalar,DeviceType> &A, const Data<DataScalar,DeviceType> &B )
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(A.rank() != B.rank(), std::invalid_argument, "A and B must have the same logical shape");
      const int rank = A.rank();
      std::vector<DimensionInfo> dimInfo(rank);
      for (int d=0; d<rank; d++)
      {
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(A.extent_int(d) != B.extent_int(d), std::invalid_argument, "A and B must have the same logical shape");
        dimInfo[d] = A.combinedDataDimensionInfo(B, d);
      }
      Data<DataScalar,DeviceType> result(dimInfo);
      return result;
    }
    
    //! Constructs a container suitable for storing the result of a matrix-vector multiply corresponding to the two provided containers.
    //! \see storeMatMat()
    //! \param A_MatData                                            [in] - logically (...,D1,D2)-dimensioned container, where D1,D2 correspond to matrix dimensions.
    //! \param transposeA                                          [in] - if true, A will be transposed prior to being multiplied by B (or B's transpose).
    //! \param B_MatData                                            [in] - logically (...,D3,D4)-dimensioned container, where D3,D4 correspond to matrix dimensions.
    //! \param transposeB                                          [in] - if true, B will be transposed prior to the multiplication by A (or A's transpose).
    static Data<DataScalar,DeviceType> allocateMatMatResult( const bool transposeA, const Data<DataScalar,DeviceType> &A_MatData, const bool transposeB, const Data<DataScalar,DeviceType> &B_MatData )
    {
      // we treat last two logical dimensions of matData as the matrix; last dimension of vecData as the vector
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(A_MatData.rank() != B_MatData.rank(), std::invalid_argument, "AmatData and BmatData have incompatible ranks");
      
      const int D1_DIM = A_MatData.rank() - 2;
      const int D2_DIM = A_MatData.rank() - 1;
      
      const int A_rows = A_MatData.extent_int(D1_DIM);
      const int A_cols = A_MatData.extent_int(D2_DIM);
      const int B_rows = B_MatData.extent_int(D1_DIM);
      const int B_cols = B_MatData.extent_int(D2_DIM);
      
      const int leftRows  = transposeA ? A_cols : A_rows;
      const int leftCols  = transposeA ? A_rows : A_cols;
      const int rightRows = transposeB ? B_cols : B_rows;
      const int rightCols = transposeB ? B_rows : B_cols;
      
      INTREPID2_TEST_FOR_EXCEPTION(leftCols != rightRows, std::invalid_argument, "incompatible matrix dimensions");
      
      Kokkos::Array<int,7> resultExtents;                      // logical extents
      Kokkos::Array<DataVariationType,7> resultVariationTypes; // for each dimension, whether the data varies in that dimension
      
      resultExtents[D1_DIM] = leftRows;
      resultExtents[D2_DIM] = rightCols;
      int resultBlockPlusDiagonalLastNonDiagonal = -1;
      if ( (A_MatData.getVariationTypes()[D1_DIM] == BLOCK_PLUS_DIAGONAL) && (B_MatData.getVariationTypes()[D1_DIM] == BLOCK_PLUS_DIAGONAL) )
      {
        // diagonal times diagonal is diagonal; the result will have the maximum of A and B's non-diagonal block size
        resultVariationTypes[D1_DIM] = BLOCK_PLUS_DIAGONAL;
        resultVariationTypes[D2_DIM] = BLOCK_PLUS_DIAGONAL;
        resultBlockPlusDiagonalLastNonDiagonal = std::max(A_MatData.blockPlusDiagonalLastNonDiagonal(), B_MatData.blockPlusDiagonalLastNonDiagonal());
      }
      
      const int resultRank = A_MatData.rank();
      
      auto A_VariationTypes = A_MatData.getVariationTypes();
      auto B_VariationTypes = B_MatData.getVariationTypes();
      
      Kokkos::Array<int,7> resultActiveDims;
      Kokkos::Array<int,7> resultDataDims;
      int resultNumActiveDims = 0; // how many of the 7 entries are actually filled in
      // the following loop is over the dimensions *prior* to matrix dimensions
      for (int i=0; i<resultRank-2; i++)
      {
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(A_MatData.extent_int(i) != B_MatData.extent_int(i), std::invalid_argument, "A and B extents must match in each non-matrix dimension");
        
        resultExtents[i] = A_MatData.extent_int(i);
        
        const DataVariationType &A_VariationType = A_VariationTypes[i];
        const DataVariationType &B_VariationType = B_VariationTypes[i];
        
        // BLOCK_PLUS_DIAGONAL should only occur in matData, and only in the matrix (final) dimensions
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(A_VariationType == BLOCK_PLUS_DIAGONAL, std::invalid_argument, "unsupported variationType");
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(B_VariationType == BLOCK_PLUS_DIAGONAL, std::invalid_argument, "unsupported variationType");
        
        int dataSize = 0;
        DataVariationType resultVariationType;
        if ((A_VariationType == GENERAL) || (B_VariationType == GENERAL))
        {
          resultVariationType = GENERAL;
          dataSize = resultExtents[i];
        }
        else if ((B_VariationType == CONSTANT) && (A_VariationType == CONSTANT))
        {
          resultVariationType = CONSTANT;
          dataSize = 1;
        }
        else if ((B_VariationType == MODULAR) && (A_VariationType == CONSTANT))
        {
          resultVariationType = MODULAR;
          dataSize = B_MatData.getVariationModulus(i);
        }
        else if ((B_VariationType == CONSTANT) && (A_VariationType == MODULAR))
        {
          resultVariationType = MODULAR;
          dataSize = A_MatData.getVariationModulus(i);
        }
        else
        {
          // both are modular.  We allow this if they agree on the modulus
          auto A_Modulus = A_MatData.getVariationModulus(i);
          auto B_Modulus = B_MatData.getVariationModulus(i);
          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(A_Modulus != B_Modulus, std::invalid_argument, "If both matrices have variation type MODULAR, they must agree on the modulus");
          resultVariationType = MODULAR;
          dataSize = A_Modulus;
        }
        resultVariationTypes[i] = resultVariationType;
        
        if (resultVariationType != CONSTANT)
        {
          resultActiveDims[resultNumActiveDims] = i;
          resultDataDims[resultNumActiveDims]   = dataSize;
          resultNumActiveDims++;
        }
      }
      
      // set things for final dimensions:
      resultExtents[D1_DIM] = leftRows;
      resultExtents[D2_DIM] = rightCols;
      
      if ( (A_MatData.getVariationTypes()[D1_DIM] == BLOCK_PLUS_DIAGONAL) && (B_MatData.getVariationTypes()[D1_DIM] == BLOCK_PLUS_DIAGONAL) )
      {
        // diagonal times diagonal is diagonal; the result will have the maximum of A and B's non-diagonal block size
        resultVariationTypes[D1_DIM] = BLOCK_PLUS_DIAGONAL;
        resultVariationTypes[D2_DIM] = BLOCK_PLUS_DIAGONAL;
        resultBlockPlusDiagonalLastNonDiagonal = std::max(A_MatData.blockPlusDiagonalLastNonDiagonal(), B_MatData.blockPlusDiagonalLastNonDiagonal());
        
        resultActiveDims[resultNumActiveDims]   = resultRank - 2;
        
        const int numDiagonalEntries    = leftRows - (resultBlockPlusDiagonalLastNonDiagonal + 1);
        const int numNondiagonalEntries = (resultBlockPlusDiagonalLastNonDiagonal + 1) * (resultBlockPlusDiagonalLastNonDiagonal + 1);
        
        resultDataDims[resultNumActiveDims] = numDiagonalEntries + numNondiagonalEntries;
        resultNumActiveDims++;
      }
      else
      {
        // pretty much the only variation types that make sense for matrix dims are GENERAL and BLOCK_PLUS_DIAGONAL
        resultVariationTypes[D1_DIM] = GENERAL;
        resultVariationTypes[D2_DIM] = GENERAL;
        
        resultActiveDims[resultNumActiveDims]   = resultRank - 2;
        resultActiveDims[resultNumActiveDims+1] = resultRank - 1;
        
        resultDataDims[resultNumActiveDims]     = leftRows;
        resultDataDims[resultNumActiveDims+1]   = rightCols;
        resultNumActiveDims += 2;
      }
      
      for (int i=resultRank; i<7; i++)
      {
        resultVariationTypes[i] = CONSTANT;
        resultExtents[i]        = 1;
      }
      
      ScalarView<DataScalar,DeviceType> data; // new view will match this one in layout and fad dimension, if any
      auto viewToMatch = A_MatData.getUnderlyingView();
      if (resultNumActiveDims == 1)
      {
        data = getMatchingViewWithLabel(viewToMatch, "Data mat-mat result", resultDataDims[0]);
      }
      else if (resultNumActiveDims == 2)
      {
        data = getMatchingViewWithLabel(viewToMatch, "Data mat-mat result", resultDataDims[0], resultDataDims[1]);
      }
      else if (resultNumActiveDims == 3)
      {
        data = getMatchingViewWithLabel(viewToMatch, "Data mat-mat result", resultDataDims[0], resultDataDims[1], resultDataDims[2]);
      }
      else if (resultNumActiveDims == 4)
      {
        data = getMatchingViewWithLabel(viewToMatch, "Data mat-mat result", resultDataDims[0], resultDataDims[1], resultDataDims[2],
                                        resultDataDims[3]);
      }
      else if (resultNumActiveDims == 5)
      {
        data = getMatchingViewWithLabel(viewToMatch, "Data mat-mat result", resultDataDims[0], resultDataDims[1], resultDataDims[2],
                                        resultDataDims[3], resultDataDims[4]);
      }
      else if (resultNumActiveDims == 6)
      {
        data = getMatchingViewWithLabel(viewToMatch, "Data mat-mat result", resultDataDims[0], resultDataDims[1], resultDataDims[2],
                                        resultDataDims[3], resultDataDims[4], resultDataDims[5]);
      }
      else // resultNumActiveDims == 7
      {
        data = getMatchingViewWithLabel(viewToMatch, "Data mat-mat result", resultDataDims[0], resultDataDims[1], resultDataDims[2],
                                        resultDataDims[3], resultDataDims[4], resultDataDims[5], resultDataDims[6]);
      }
      
      return Data<DataScalar,DeviceType>(data,resultRank,resultExtents,resultVariationTypes,resultBlockPlusDiagonalLastNonDiagonal);
    }
    
    //! Constructs a container suitable for storing the result of a contraction over the final dimensions of the two provided containers.  The two containers must have the same logical shape.
    //! \see storeInPlaceCombination()
    //! \param A  [in] - the first data container.
    //! \param B  [in] - the second data container.  Must have the same logical shape as A.
    //! \param numContractionDims [in] - the number of dimensions over which the contraction should take place.
    //! \return A numContractionDims-rank-lower container with the same logical shape as A and B in all but the last dimensions.
    static Data<DataScalar,DeviceType> allocateContractionResult( const Data<DataScalar,DeviceType> &A, const Data<DataScalar,DeviceType> &B, const int &numContractionDims )
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(A.rank() != B.rank(), std::invalid_argument, "A and B must have the same logical shape");
      const int rank = A.rank();
      const int resultRank = rank - numContractionDims;
      std::vector<DimensionInfo> dimInfo(resultRank);
      for (int d=0; d<resultRank; d++)
      {
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(A.extent_int(d) != B.extent_int(d), std::invalid_argument, "A and B must have the same logical shape");
        dimInfo[d] = A.combinedDataDimensionInfo(B, d);
      }
      Data<DataScalar,DeviceType> result(dimInfo);
      return result;
    }
    
    //! Constructs a container suitable for storing the result of a contraction over the final dimension of the two provided containers.  The two containers must have the same logical shape.
    //! \see storeInPlaceCombination()
    //! \param A  [in] - the first data container.
    //! \param B  [in] - the second data container.  Must have the same logical shape as A.
    //! \return A 1-rank-lower container with the same logical shape as A and B in all but the last dimension.
    static Data<DataScalar,DeviceType> allocateDotProductResult( const Data<DataScalar,DeviceType> &A, const Data<DataScalar,DeviceType> &B )
    {
      return allocateContractionResult(A, B, 1);
    }
    
    //! Constructs a container suitable for storing the result of a matrix-vector multiply corresponding to the two provided containers.
    //! \see storeMatVec()
    static Data<DataScalar,DeviceType> allocateMatVecResult( const Data<DataScalar,DeviceType> &matData, const Data<DataScalar,DeviceType> &vecData, const bool transposeMatrix = false )
    {
      // we treat last two logical dimensions of matData as the matrix; last dimension of vecData as the vector
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(matData.rank() != vecData.rank() + 1, std::invalid_argument, "matData and vecData have incompatible ranks");
      const int vecDim  = vecData.extent_int(vecData.rank() - 1);
      
      const int D1_DIM = matData.rank() - 2;
      const int D2_DIM = matData.rank() - 1;
      
      const int matRows = matData.extent_int(D1_DIM);
      const int matCols = matData.extent_int(D2_DIM);
      
      const int rows  = transposeMatrix ? matCols : matRows;
      const int cols  = transposeMatrix ? matRows : matCols;
      
      const int resultRank = vecData.rank();
      
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(cols != vecDim, std::invalid_argument, "matData column count != vecData dimension");
      
      Kokkos::Array<int,7> resultExtents;                      // logical extents
      Kokkos::Array<DataVariationType,7> resultVariationTypes; // for each dimension, whether the data varies in that dimension
      auto vecVariationTypes = vecData.getVariationTypes();
      auto matVariationTypes = matData.getVariationTypes();
      
      Kokkos::Array<int,7> resultActiveDims;
      Kokkos::Array<int,7> resultDataDims;
      int resultNumActiveDims = 0; // how many of the 7 entries are actually filled in
      // the following loop is over the dimensions *prior* to matrix/vector dimensions
      for (int i=0; i<resultRank-1; i++)
      {
        resultExtents[i] = vecData.extent_int(i);
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vecData.extent_int(i) != matData.extent_int(i), std::invalid_argument, "matData and vecData extents must match in each non-matrix/vector dimension");
        
        const DataVariationType &vecVariationType = vecVariationTypes[i];
        const DataVariationType &matVariationType = matVariationTypes[i];
        
        // BLOCK_PLUS_DIAGONAL should only occur in matData, and only in the matrix (final) dimensions
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vecVariationType == BLOCK_PLUS_DIAGONAL, std::invalid_argument, "unsupported variationType");
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(matVariationType == BLOCK_PLUS_DIAGONAL, std::invalid_argument, "unsupported variationType");
        
        int dataSize = 0;
        DataVariationType resultVariationType;
        if ((vecVariationType == GENERAL) || (matVariationType == GENERAL))
        {
          resultVariationType = GENERAL;
          dataSize = resultExtents[i];
        }
        else if ((matVariationType == CONSTANT) && (vecVariationType == CONSTANT))
        {
          resultVariationType = CONSTANT;
          dataSize = 1;
        }
        else if ((matVariationType == MODULAR) && (vecVariationType == CONSTANT))
        {
          resultVariationType = MODULAR;
          dataSize = matData.getVariationModulus(i);
        }
        else if ((matVariationType == CONSTANT) && (vecVariationType == MODULAR))
        {
          resultVariationType = MODULAR;
          dataSize = matData.getVariationModulus(i);
        }
        else
        {
          // both are modular.  We allow this if they agree on the modulus
          auto matModulus = matData.getVariationModulus(i);
          auto vecModulus = vecData.getVariationModulus(i);
          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(matModulus != vecModulus, std::invalid_argument, "If matrix and vector both have variation type MODULAR, they must agree on the modulus");
          resultVariationType = MODULAR;
          dataSize = matModulus;
        }
        resultVariationTypes[i] = resultVariationType;
        
        if (resultVariationType != CONSTANT)
        {
          resultActiveDims[resultNumActiveDims] = i;
          resultDataDims[resultNumActiveDims]   = dataSize;
          resultNumActiveDims++;
        }
      }
      // for the final dimension, the variation type is always GENERAL
      // (Some combinations, e.g. CONSTANT/CONSTANT *would* generate a CONSTANT result, but constant matrices don't make a lot of sense beyond 1x1 matrices…)
      resultActiveDims[resultNumActiveDims]     = resultRank - 1;
      resultDataDims[resultNumActiveDims]       = rows;
      resultNumActiveDims++;
      
      for (int i=resultRank; i<7; i++)
      {
        resultVariationTypes[i] = CONSTANT;
        resultExtents[i]        = 1;
      }
      resultVariationTypes[resultRank-1] = GENERAL;
      resultExtents[resultRank-1]        = rows;
      
      ScalarView<DataScalar,DeviceType> data;
      if (resultNumActiveDims == 1)
      {
        data = matData.allocateDynRankViewMatchingUnderlying(resultDataDims[0]);
      }
      else if (resultNumActiveDims == 2)
      {
        data = matData.allocateDynRankViewMatchingUnderlying(resultDataDims[0], resultDataDims[1]);
      }
      else if (resultNumActiveDims == 3)
      {
        data = matData.allocateDynRankViewMatchingUnderlying(resultDataDims[0], resultDataDims[1], resultDataDims[2]);
      }
      else if (resultNumActiveDims == 4)
      {
        data = matData.allocateDynRankViewMatchingUnderlying(resultDataDims[0], resultDataDims[1], resultDataDims[2],
                                                             resultDataDims[3]);
      }
      else if (resultNumActiveDims == 5)
      {
        data = matData.allocateDynRankViewMatchingUnderlying(resultDataDims[0], resultDataDims[1], resultDataDims[2],
                                                             resultDataDims[3], resultDataDims[4]);
      }
      else if (resultNumActiveDims == 6)
      {
        data = matData.allocateDynRankViewMatchingUnderlying(resultDataDims[0], resultDataDims[1], resultDataDims[2],
                                                             resultDataDims[3], resultDataDims[4], resultDataDims[5]);
      }
      else // resultNumActiveDims == 7
      {
        data = matData.allocateDynRankViewMatchingUnderlying(resultDataDims[0], resultDataDims[1], resultDataDims[2],
                                                             resultDataDims[3], resultDataDims[4], resultDataDims[5], resultDataDims[6]);
      }
      
      return Data<DataScalar,DeviceType>(data,resultRank,resultExtents,resultVariationTypes);
    }
    
    //! returns an MDRangePolicy over the underlying data extents (but with the logical shape).
    template<int rank>
    enable_if_t<(rank!=1) && (rank!=7), Kokkos::MDRangePolicy<typename DeviceType::execution_space,Kokkos::Rank<rank>> >
    dataExtentRangePolicy()
    {
      using ExecutionSpace = typename DeviceType::execution_space;
      Kokkos::Array<int,rank> startingOrdinals;
      Kokkos::Array<int,rank> extents;
      
      for (int d=0; d<rank; d++)
      {
        startingOrdinals[d] = 0;
        extents[d] = getDataExtent(d);
      }
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<rank>>(startingOrdinals,extents);
      return policy;
    }
    
    //! returns an MDRangePolicy over the first six underlying data extents (but with the logical shape).
    template<int rank>
    enable_if_t<rank==7, Kokkos::MDRangePolicy<typename DeviceType::execution_space,Kokkos::Rank<6>> >
    dataExtentRangePolicy()
    {
      using ExecutionSpace = typename DeviceType::execution_space;
      Kokkos::Array<int,6> startingOrdinals;
      Kokkos::Array<int,6> extents;
      
      for (int d=0; d<6; d++)
      {
        startingOrdinals[d] = 0;
        extents[d] = getDataExtent(d);
      }
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<6>>(startingOrdinals,extents);
      return policy;
    }
    
    template<int rank>
    inline
    enable_if_t<rank==1, Kokkos::RangePolicy<typename DeviceType::execution_space> >
    dataExtentRangePolicy()
    {
      using ExecutionSpace = typename DeviceType::execution_space;
      Kokkos::RangePolicy<ExecutionSpace> policy(ExecutionSpace(),0,getDataExtent(0));
      return policy;
    }
    
    //! Creates a new Data object with the same underlying view, but with the specified logical rank, extents, and variation types.
    Data shallowCopy(const int rank, const Kokkos::Array<int,7> &extents, const Kokkos::Array<DataVariationType,7> &variationTypes) const
    {
      switch (dataRank_)
      {
        case 1: return Data(rank, data1_, extents, variationTypes);
        case 2: return Data(rank, data2_, extents, variationTypes);
        case 3: return Data(rank, data3_, extents, variationTypes);
        case 4: return Data(rank, data4_, extents, variationTypes);
        case 5: return Data(rank, data5_, extents, variationTypes);
        case 6: return Data(rank, data6_, extents, variationTypes);
        case 7: return Data(rank, data7_, extents, variationTypes);
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unhandled dataRank_");
      }
    }
    
    //! Places the result of a contraction along the final dimension of A and B into this data container.
    void storeDotProduct(const Data<DataScalar,DeviceType> &A, const Data<DataScalar,DeviceType> &B)
    {
      const int D_DIM = A.rank() - 1;
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(A.extent_int(D_DIM) != B.extent_int(D_DIM), std::invalid_argument, "A and B have different extents");
      const int vectorComponents = A.extent_int(D_DIM);
      
      // shallow copy of this to avoid implicit references to this in call to getWritableEntry() below
      Data<DataScalar,DeviceType> thisData = *this;
      
      using ExecutionSpace = typename DeviceType::execution_space;
      // note the use of getDataExtent() below: we only range over the possibly-distinct entries
      if (rank_ == 1) // contraction result rank; e.g., (P)
      {
        Kokkos::parallel_for("compute dot product", getDataExtent(0),
        KOKKOS_LAMBDA (const int &pointOrdinal) {
          auto & val = thisData.getWritableEntry(pointOrdinal);
          val = 0;
          for (int i=0; i<vectorComponents; i++)
          {
            val += A(pointOrdinal,i) * B(pointOrdinal,i);
          }
        });
      }
      else if (rank_ == 2) // contraction result rank; e.g., (C,P)
      {
        // typical case for e.g. gradient data: (C,P,D)
        auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{getDataExtent(0),getDataExtent(1)});
        Kokkos::parallel_for("compute dot product", policy,
        KOKKOS_LAMBDA (const int &cellOrdinal, const int &pointOrdinal) {
          auto & val = thisData.getWritableEntry(cellOrdinal, pointOrdinal);
          val = 0;
          for (int i=0; i<vectorComponents; i++)
          {
            val += A(cellOrdinal,pointOrdinal,i) * B(cellOrdinal,pointOrdinal,i);
          }
        });
      }
      else if (rank_ == 3)
      {
        auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{getDataExtent(0),getDataExtent(1),getDataExtent(2)});
        Kokkos::parallel_for("compute dot product", policy,
        KOKKOS_LAMBDA (const int &cellOrdinal, const int &pointOrdinal, const int &d) {
          auto & val = thisData.getWritableEntry(cellOrdinal, pointOrdinal,d);
          val = 0;
          for (int i=0; i<vectorComponents; i++)
          {
            val += A(cellOrdinal,pointOrdinal,d,i) * B(cellOrdinal,pointOrdinal,d,i);
          }
        });
      }
      else
      {
        // TODO: handle other cases
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::logic_error, "rank not yet supported");
      }
    }
    
    //! Places the result of an in-place combination (e.g., entrywise sum) into this data container.
    template<class BinaryOperator>
    void storeInPlaceCombination(const Data<DataScalar,DeviceType> &A, const Data<DataScalar,DeviceType> &B, BinaryOperator binaryOperator);
    
    //! stores the in-place (entrywise) sum, A .+ B, into this container.
    void storeInPlaceSum(const Data<DataScalar,DeviceType> &A, const Data<DataScalar,DeviceType> &B)
    {
      ScalarSumFunctor<DataScalar> sum;
      storeInPlaceCombination(A, B, sum);
    }
    
    //! stores the in-place (entrywise) product, A .* B, into this container.
    void storeInPlaceProduct(const Data<DataScalar,DeviceType> &A, const Data<DataScalar,DeviceType> &B)
    {
      ScalarProductFunctor<DataScalar> product;
      storeInPlaceCombination(A, B, product);
    }
    
    //! stores the in-place (entrywise) difference, A .- B, into this container.
    void storeInPlaceDifference(const Data<DataScalar,DeviceType> &A, const Data<DataScalar,DeviceType> &B)
    {
      ScalarDifferenceFunctor<DataScalar> difference;
      storeInPlaceCombination(A, B, difference);
    }
    
    //! stores the in-place (entrywise) quotient, A ./ B, into this container.
    void storeInPlaceQuotient(const Data<DataScalar,DeviceType> &A, const Data<DataScalar,DeviceType> &B)
    {
      ScalarQuotientFunctor<DataScalar> quotient;
      storeInPlaceCombination(A, B, quotient);
    }
    
    //! Places the result of a matrix-vector multiply corresponding to the two provided containers into this Data container.  This Data container should have been constructed by a call to allocateMatVecResult(), or should match such a container in underlying data extent and variation types.
    void storeMatVec( const Data<DataScalar,DeviceType> &matData, const Data<DataScalar,DeviceType> &vecData, const bool transposeMatrix = false )
    {
      // TODO: add a compile-time (SFINAE-type) guard against DataScalar types that do not support arithmetic operations.  (We support Orientation as a DataScalar type; it might suffice just to compare DataScalar to Orientation, and eliminate this method for that case.)
      // TODO: check for invalidly shaped containers.
      
      const int D1_DIM = matData.rank() - 2;
      const int D2_DIM = matData.rank() - 1;
      
      const int matRows = matData.extent_int(D1_DIM);
      const int matCols = matData.extent_int(D2_DIM);
      
      const int rows  = transposeMatrix ? matCols : matRows;
      const int cols  = transposeMatrix ? matRows : matCols;
      
      // shallow copy of this to avoid implicit references to this in call to getWritableEntry() below
      Data<DataScalar,DeviceType> thisData = *this;
      
      using ExecutionSpace = typename DeviceType::execution_space;
      // note the use of getDataExtent() below: we only range over the possibly-distinct entries
      if (rank_ == 3)
      {
        // typical case for e.g. gradient data: (C,P,D)
        auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{getDataExtent(0),getDataExtent(1),rows});
        Kokkos::parallel_for("compute mat-vec", policy,
        KOKKOS_LAMBDA (const int &cellOrdinal, const int &pointOrdinal, const int &i) {
          auto & val_i = thisData.getWritableEntry(cellOrdinal, pointOrdinal, i);
          val_i = 0;
          for (int j=0; j<cols; j++)
          {
            const auto & mat_ij  = transposeMatrix ? matData(cellOrdinal,pointOrdinal,j,i) : matData(cellOrdinal,pointOrdinal,i,j);
            val_i += mat_ij * vecData(cellOrdinal,pointOrdinal,j);
          }
        });
      }
      else if (rank_ == 2)
      {
        //
        auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{getDataExtent(0),rows});
        Kokkos::parallel_for("compute mat-vec", policy,
        KOKKOS_LAMBDA (const int &vectorOrdinal, const int &i) {
          auto & val_i = thisData.getWritableEntry(vectorOrdinal, i);
          val_i = 0;
          for (int j=0; j<cols; j++)
          {
            const auto & mat_ij  = transposeMatrix ? matData(vectorOrdinal,j,i) : matData(vectorOrdinal,i,j);
            val_i += mat_ij * vecData(vectorOrdinal,j);
          }
        });
      }
      else if (rank_ == 1)
      {
        // single-vector case
        Kokkos::RangePolicy<ExecutionSpace> policy(0,rows);
        Kokkos::parallel_for("compute mat-vec", policy,
        KOKKOS_LAMBDA (const int &i) {
          auto & val_i = thisData.getWritableEntry(i);
          val_i = 0;
          for (int j=0; j<cols; j++)
          {
            const auto & mat_ij  = transposeMatrix ? matData(j,i) : matData(i,j);
            val_i += mat_ij * vecData(j);
          }
        });
      }
      else
      {
        // TODO: handle other cases
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::logic_error, "rank not yet supported");
      }
    }
    
    //! Places the result of a matrix-matrix multiply corresponding to the two provided containers into this Data container.  This Data container should have been constructed by a call to allocateMatMatResult(), or should match such a container in underlying data extent and variation types.
    //! \see allocateMatMat()
    //! \param A_MatData                                            [in] - logically (...,D1,D2)-dimensioned container, where D1,D2 correspond to matrix dimensions.
    //! \param transposeA                                          [in] - if true, A will be transposed prior to being multiplied by B (or B's transpose).
    //! \param B_MatData                                            [in] - logically (...,D3,D4)-dimensioned container, where D3,D4 correspond to matrix dimensions.
    //! \param transposeB                                          [in] - if true, B will be transposed prior to the multiplication by A (or A's transpose).
    void storeMatMat( const bool transposeA, const Data<DataScalar,DeviceType> &A_MatData, const bool transposeB, const Data<DataScalar,DeviceType> &B_MatData )
    {
      // TODO: add a compile-time (SFINAE-type) guard against DataScalar types that do not support arithmetic operations.  (We support Orientation as a DataScalar type; it might suffice just to compare DataScalar to Orientation, and eliminate this method for that case.)
      // TODO: check for invalidly shaped containers.
      
      // we treat last two logical dimensions of matData as the matrix; last dimension of vecData as the vector
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(A_MatData.rank() != B_MatData.rank(), std::invalid_argument, "AmatData and BmatData have incompatible ranks");
      
      const int D1_DIM = A_MatData.rank() - 2;
      const int D2_DIM = A_MatData.rank() - 1;
      
      const int A_rows = A_MatData.extent_int(D1_DIM);
      const int A_cols = A_MatData.extent_int(D2_DIM);
      const int B_rows = B_MatData.extent_int(D1_DIM);
      const int B_cols = B_MatData.extent_int(D2_DIM);
     
      const int leftRows  = transposeA ? A_cols : A_rows;
      const int leftCols  = transposeA ? A_rows : A_cols;
      const int rightCols = transposeB ? B_rows : B_cols;
      
#ifdef INTREPID2_HAVE_DEBUG
      const int rightRows = transposeB ? B_cols : B_rows;
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(leftCols != rightRows, std::invalid_argument, "inner dimensions do not match");
#endif
      
      // shallow copy of this to avoid implicit references to this in call to getWritableEntry() below
      Data<DataScalar,DeviceType> thisData = *this;
      
      using ExecutionSpace = typename DeviceType::execution_space;
      
      const int diagonalStart = (variationType_[D1_DIM] == BLOCK_PLUS_DIAGONAL) ? blockPlusDiagonalLastNonDiagonal_ + 1 : leftRows;
      // note the use of getDataExtent() below: we only range over the possibly-distinct entries
      if (rank_ == 3)
      {
        // (C,D,D), say
        auto policy = Kokkos::RangePolicy<ExecutionSpace>(0,getDataExtent(0));
        Kokkos::parallel_for("compute mat-mat", policy,
        KOKKOS_LAMBDA (const int &matrixOrdinal) {
          for (int i=0; i<diagonalStart; i++)
          {
            for (int j=0; j<rightCols; j++)
            {
              auto & val_ij = thisData.getWritableEntry(matrixOrdinal, i, j);
              val_ij = 0;
              for (int k=0; k<leftCols; k++)
              {
                const auto & left  = transposeA ? A_MatData(matrixOrdinal,k,i) : A_MatData(matrixOrdinal,i,k);
                const auto & right = transposeB ? B_MatData(matrixOrdinal,j,k) : B_MatData(matrixOrdinal,k,j);
                val_ij += left * right;
              }
            }
          }
          for (int i=diagonalStart; i<leftRows; i++)
          {
            auto & val_ii = thisData.getWritableEntry(matrixOrdinal, i, i);
            const auto & left  = A_MatData(matrixOrdinal,i,i);
            const auto & right = B_MatData(matrixOrdinal,i,i);
            val_ii = left * right;
          }
        });
      }
      else if (rank_ == 4)
      {
        // (C,P,D,D), perhaps
        auto policy = Kokkos::MDRangePolicy<ExecutionSpace, Kokkos::Rank<2> >({0,0},{getDataExtent(0),getDataExtent(1)});
        if (underlyingMatchesLogical_) // receiving data object is completely expanded
        {
          Kokkos::parallel_for("compute mat-mat", policy,
          KOKKOS_LAMBDA (const int &cellOrdinal, const int &pointOrdinal) {
            for (int i=0; i<leftRows; i++)
            {
              for (int j=0; j<rightCols; j++)
              {
                auto & val_ij = thisData.getUnderlyingView4()(cellOrdinal,pointOrdinal, i, j);
                val_ij = 0;
                for (int k=0; k<leftCols; k++)
                {
                  const auto & left  = transposeA ? A_MatData(cellOrdinal,pointOrdinal,k,i) : A_MatData(cellOrdinal,pointOrdinal,i,k);
                  const auto & right = transposeB ? B_MatData(cellOrdinal,pointOrdinal,j,k) : B_MatData(cellOrdinal,pointOrdinal,k,j);
                  val_ij += left * right;
                }
              }
            }
          });
        }
        else
        {
          Kokkos::parallel_for("compute mat-mat", policy,
          KOKKOS_LAMBDA (const int &cellOrdinal, const int &pointOrdinal) {
            for (int i=0; i<diagonalStart; i++)
            {
              for (int j=0; j<rightCols; j++)
              {
                auto & val_ij = thisData.getWritableEntry(cellOrdinal,pointOrdinal, i, j);
                val_ij = 0;
                for (int k=0; k<leftCols; k++)
                {
                  const auto & left  = transposeA ? A_MatData(cellOrdinal,pointOrdinal,k,i) : A_MatData(cellOrdinal,pointOrdinal,i,k);
                  const auto & right = transposeB ? B_MatData(cellOrdinal,pointOrdinal,j,k) : B_MatData(cellOrdinal,pointOrdinal,k,j);
                  val_ij += left * right;
                }
              }
            }
            for (int i=diagonalStart; i<leftRows; i++)
            {
              auto & val_ii = thisData.getWritableEntry(cellOrdinal,pointOrdinal, i, i);
              const auto & left  = A_MatData(cellOrdinal,pointOrdinal,i,i);
              const auto & right = B_MatData(cellOrdinal,pointOrdinal,i,i);
              val_ii = left * right;
            }
          });
        }
      }
      else
      {
        // TODO: handle other cases
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::logic_error, "rank not yet supported");
      }
    }
    
    //! returns true for containers that have data; false for those that don't (namely, those that have been constructed by the default constructor).
    KOKKOS_INLINE_FUNCTION constexpr bool isValid() const
    {
      return extents_[0] > 0;
    }
    
    //! Returns the logical rank of the Data container.
    KOKKOS_INLINE_FUNCTION
    unsigned rank() const
    {
      return rank_;
    }
    
 /** \brief sets the logical extent in the specified dimension.  If needed, the underlying data container is resized.
     \param [in] d - the logical dimension in which the extent is to be changed
     \param [in] newExtent - the new extent
     \note Not supported for dimensions in which the variation type is BLOCK_PLUS_DIAGONAL.
     \note If the variation type is MODULAR, the existing modulus must evenly divide the new extent; the underlying data structure will not be resized in this case.
     */
    void setExtent(const ordinal_type &d, const ordinal_type &newExtent)
    {
      INTREPID2_TEST_FOR_EXCEPTION(variationType_[d] == BLOCK_PLUS_DIAGONAL, std::invalid_argument, "setExtent is not supported for BLOCK_PLUS_DIAGONAL dimensions");
      
      if (variationType_[d] == MODULAR)
      {
        bool dividesEvenly = ((newExtent / variationModulus_[d]) * variationModulus_[d] == newExtent);
        INTREPID2_TEST_FOR_EXCEPTION(!dividesEvenly, std::invalid_argument, "when setExtent is called on dimenisions with MODULAR variation, the modulus must divide the new extent evenly");
      }
      
      if ((newExtent != extents_[d]) && (variationType_[d] == GENERAL))
      {
        // then we need to resize; let's determine the full set of new extents
        std::vector<ordinal_type> newExtents(dataRank_,-1);
        for (int r=0; r<dataRank_; r++)
        {
          if (activeDims_[r] == d)
          {
            // this is the changed dimension
            newExtents[r] = newExtent;
          }
          else
          {
            // unchanged; copy from existing
            newExtents[r] = getUnderlyingViewExtent(r);
          }
        }
        
        switch (dataRank_)
        {
          case 1: Kokkos::resize(data1_,newExtents[0]);
            break;
          case 2: Kokkos::resize(data2_,newExtents[0],newExtents[1]);
            break;
          case 3: Kokkos::resize(data3_,newExtents[0],newExtents[1],newExtents[2]);
            break;
          case 4: Kokkos::resize(data4_,newExtents[0],newExtents[1],newExtents[2],newExtents[3]);
            break;
          case 5: Kokkos::resize(data5_,newExtents[0],newExtents[1],newExtents[2],newExtents[3],newExtents[4]);
            break;
          case 6: Kokkos::resize(data6_,newExtents[0],newExtents[1],newExtents[2],newExtents[3],newExtents[4],newExtents[5]);
            break;
          case 7: Kokkos::resize(data7_,newExtents[0],newExtents[1],newExtents[2],newExtents[3],newExtents[4],newExtents[5],newExtents[6]);
            break;
          default: INTREPID2_TEST_FOR_EXCEPTION(true, std::logic_error, "Unexpected dataRank_ value");
        }
      }
      
      extents_[d] = newExtent;
    }
    
    //! Returns true if the underlying container has exactly the same rank and extents as the logical container.
    KOKKOS_INLINE_FUNCTION
    bool underlyingMatchesLogical() const
    {
      return underlyingMatchesLogical_;
    }
  };

  template<class DataScalar, typename DeviceType>
  KOKKOS_INLINE_FUNCTION constexpr unsigned rank(const Data<DataScalar, DeviceType>& D) {
    return D.rank();
  }
}

// we do ETI for doubles and default ExecutionSpace's device_type
extern template class Intrepid2::Data<double,Kokkos::DefaultExecutionSpace::device_type>;

#include "Intrepid2_DataDef.hpp"

#endif /* Intrepid2_Data_h */
