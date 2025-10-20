// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Tests against structured integration facilities - "synthetic" test cases (i.e., no geometry specified).
    \author Nathan V. Roberts
*/

#include "Kokkos_Core.hpp"

#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_UnitTestRepository.hpp>

#include <Intrepid2_TestUtils.hpp>

#include "Intrepid2_ArgExtractor.hpp"
#include "Intrepid2_DataDimensionInfo.hpp"
#include "Intrepid2_DataFunctors.hpp"
#include "Intrepid2_DataVariationType.hpp"

#include "Intrepid2_ScalarView.hpp"

template<class DataScalar,typename DeviceType>
using ScalarView = Intrepid2::ScalarView<DataScalar, DeviceType>;

using ordinal_type = Intrepid2::ordinal_type;

using DimensionInfo = Intrepid2::DimensionInfo;
using DataVariationType = Intrepid2::DataVariationType;

template<class DataScalar, int rank>
using RankExpander = Intrepid2::RankExpander<DataScalar, rank>;

using std::enable_if_t;

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
    
//    ScalarView<DataScalar,DeviceType> zeroView_; // one-entry (zero); used to allow getEntry() to return 0 for off-diagonal entries in BLOCK_PLUS_DIAGONAL
    
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
    {}

    public:

    //! Returns an l-value reference to the specified logical entry in the underlying view.  Note that for variation types other than GENERAL, multiple valid argument sets will refer to the same memory location.  Intended for Intrepid2 developers and expert users only.  If passThroughBlockDiagonalArgs is TRUE, the corresponding arguments are interpreted as entries in the 1D packed matrix rather than as logical 2D matrix row and column.
    template<class ...IntArgs>
    KOKKOS_INLINE_FUNCTION
    reference_type getWritableEntryWithPassThroughOption(const bool &passThroughBlockDiagonalArgs, const IntArgs... intArgs) const
    {}
    
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
    dataRank_(0), extents_({0,0,0,0,0,0,0}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(-1), rank_(dimInfoVector.size())
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
          const bool isBlockPlusDiagonal = (variationType_[d] == Intrepid2::BLOCK_PLUS_DIAGONAL);
          const bool isSecondBlockPlusDiagonal = isBlockPlusDiagonal && blockPlusDiagonalEncountered;
          if (isBlockPlusDiagonal)
          {
            blockPlusDiagonalEncountered = true;
            blockPlusDiagonalLastNonDiagonal_ = dimInfo.blockPlusDiagonalLastNonDiagonal;
          }
          if ((variationType_[d] != Intrepid2::CONSTANT) && (!isSecondBlockPlusDiagonal))
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
         Kokkos::Array<DataVariationType,7> {Intrepid2::GENERAL,Intrepid2::GENERAL,Intrepid2::GENERAL,Intrepid2::GENERAL,Intrepid2::GENERAL,Intrepid2::GENERAL,Intrepid2::GENERAL}, -1)
    {}
    
    //! Constructor that accepts a DynRankView as an argument.  The data belonging to the DynRankView will be copied into a new View of matching dimensions.
    template<size_t rank, class ...DynRankViewProperties>
    Data(const Kokkos::DynRankView<DataScalar,DeviceType, DynRankViewProperties...> &data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
    :
    dataRank_(data.rank()), extents_({1,1,1,1,1,1,1}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
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
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
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
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
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
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
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
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
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
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
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
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
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
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
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
    dataRank_(data.rank), extents_({1,1,1,1,1,1,1}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
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
    dataRank_(1), extents_({1,1,1,1,1,1,1}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(-1), rank_(rank)
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
    dataRank_(0), extents_({0,0,0,0,0,0,0}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(-1), rank_(0)
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
    
    template<int underlyingRank, class ViewScalar>
    KOKKOS_INLINE_FUNCTION
    void setUnderlyingView(const Kokkos::View<ViewScalar, DeviceType> & view)
    {}
    
    //! Returns a DynRankView constructed atop the same underlying data as the fixed-rank Kokkos::View used internally.
    ScalarView<DataScalar,DeviceType> getUnderlyingView() const
    {
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
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
    }
    
    //! Copies from the provided DynRankView into the underlying Kokkos::View container storing the unique data.
    void copyDataFromDynRankViewMatchingUnderlying(const ScalarView<DataScalar,DeviceType> &dynRankView) const
    {}
    
    //! returns the true extent of the data corresponding to the logical dimension provided; if the data does not vary in that dimension, returns 1
    KOKKOS_INLINE_FUNCTION int getDataExtent(const ordinal_type &d) const
    {
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
    
  };

namespace
{



TEUCHOS_UNIT_TEST( PR14546, AllocationIssue )
{
  using DeviceType = Intrepid2::DefaultTestDeviceType;
  using BigDataArray = Kokkos::Array< Intrepid2::Data<double,DeviceType>, 20000 >;
  BigDataArray vectorComponents;
}

} // anonymous namespace
