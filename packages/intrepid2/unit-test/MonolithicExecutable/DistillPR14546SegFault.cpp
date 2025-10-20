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

  template<class DataScalar,typename DeviceType>
  class A {
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
    bool underlyingMatchesLogical_;   // if true, this A object has the same rank and extent as the underlying view
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
    //! Generic data copying method to allow construction of A object from DynRankViews for which deep_copy() to the underlying view would be disallowed.  This method made public to allow CUDA compilation (because it contains a Kokkos lambda).
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
        case 1: data1_ = Kokkos::View<DataScalar*,       DeviceType>("Intrepid2 A", data.extent_int(0)); break;
        case 2: data2_ = Kokkos::View<DataScalar**,      DeviceType>("Intrepid2 A", data.extent_int(0), data.extent_int(1)); break;
        case 3: data3_ = Kokkos::View<DataScalar***,     DeviceType>("Intrepid2 A", data.extent_int(0), data.extent_int(1), data.extent_int(2)); break;
        case 4: data4_ = Kokkos::View<DataScalar****,    DeviceType>("Intrepid2 A", data.extent_int(0), data.extent_int(1), data.extent_int(2), data.extent_int(3)); break;
        case 5: data5_ = Kokkos::View<DataScalar*****,   DeviceType>("Intrepid2 A", data.extent_int(0), data.extent_int(1), data.extent_int(2), data.extent_int(3), data.extent_int(4)); break;
        case 6: data6_ = Kokkos::View<DataScalar******,  DeviceType>("Intrepid2 A", data.extent_int(0), data.extent_int(1), data.extent_int(2), data.extent_int(3), data.extent_int(4), data.extent_int(5)); break;
        case 7: data7_ = Kokkos::View<DataScalar*******, DeviceType>("Intrepid2 A", data.extent_int(0), data.extent_int(1), data.extent_int(2), data.extent_int(3), data.extent_int(4), data.extent_int(5), data.extent_int(6)); break;
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
    A(std::vector<DimensionInfo> dimInfoVector)
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
          case 1: data1_ = Kokkos::View<DataScalar*,       DeviceType>("Intrepid2 A", dataExtents[0]); break;
          case 2: data2_ = Kokkos::View<DataScalar**,      DeviceType>("Intrepid2 A", dataExtents[0], dataExtents[1]); break;
          case 3: data3_ = Kokkos::View<DataScalar***,     DeviceType>("Intrepid2 A", dataExtents[0], dataExtents[1], dataExtents[2]); break;
          case 4: data4_ = Kokkos::View<DataScalar****,    DeviceType>("Intrepid2 A", dataExtents[0], dataExtents[1], dataExtents[2], dataExtents[3]); break;
          case 5: data5_ = Kokkos::View<DataScalar*****,   DeviceType>("Intrepid2 A", dataExtents[0], dataExtents[1], dataExtents[2], dataExtents[3], dataExtents[4]); break;
          case 6: data6_ = Kokkos::View<DataScalar******,  DeviceType>("Intrepid2 A", dataExtents[0], dataExtents[1], dataExtents[2], dataExtents[3], dataExtents[4], dataExtents[5]); break;
          case 7: data7_ = Kokkos::View<DataScalar*******, DeviceType>("Intrepid2 A", dataExtents[0], dataExtents[1], dataExtents[2], dataExtents[3], dataExtents[4], dataExtents[5], dataExtents[6]); break;
          default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid data rank");
        }
      }
      setActiveDims();
    }
    
    //! DynRankView constructor.  Will copy to a View of appropriate rank.
    A(const ScalarView<DataScalar,DeviceType> &data, int rank, Kokkos::Array<int,7> extents, Kokkos::Array<DataVariationType,7> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
    :
    dataRank_(data.rank()), extents_(extents), variationType_(variationType), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
    {
      allocateAndCopyFromDynRankView(data);
      setActiveDims();
    }
    
    //! copy-like constructor for differing device type, but same memory space.  This does a shallow copy of the underlying view.
    template<typename OtherDeviceType, class = typename std::enable_if< std::is_same<typename DeviceType::memory_space, typename OtherDeviceType::memory_space>::value>::type,
                                       class = typename std::enable_if<!std::is_same<DeviceType,OtherDeviceType>::value>::type>
    A(const A<DataScalar,OtherDeviceType> &data)
    :
    dataRank_(data.getUnderlyingViewRank()), extents_(data.getExtents()), variationType_(data.getVariationTypes()), blockPlusDiagonalLastNonDiagonal_(data.blockPlusDiagonalLastNonDiagonal()), rank_(data.rank())
    {
//      std::cout << "Entered copy-like A constructor.\n";
      if (dataRank_ != 0) // dataRank_ == 0 indicates an invalid A object (a placeholder, can indicate zero value)
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
    A(const A<DataScalar,OtherDeviceType> &data)
    :
    dataRank_(data.getUnderlyingViewRank()), extents_(data.getExtents()), variationType_(data.getVariationTypes()), blockPlusDiagonalLastNonDiagonal_(data.blockPlusDiagonalLastNonDiagonal()), rank_(data.rank())
    {
//      std::cout << "Entered copy-like A constructor.\n";
      if (dataRank_ != 0) // dataRank_ == 0 indicates an invalid A object (a placeholder, can indicate zero value)
      {
        const auto view = data.getUnderlyingView();
        switch (dataRank_)
        {
          case 1: data1_ = Kokkos::View<DataScalar*,       DeviceType>("Intrepid2 A", view.extent_int(0)); break;
          case 2: data2_ = Kokkos::View<DataScalar**,      DeviceType>("Intrepid2 A", view.extent_int(0), view.extent_int(1)); break;
          case 3: data3_ = Kokkos::View<DataScalar***,     DeviceType>("Intrepid2 A", view.extent_int(0), view.extent_int(1), view.extent_int(2)); break;
          case 4: data4_ = Kokkos::View<DataScalar****,    DeviceType>("Intrepid2 A", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3)); break;
          case 5: data5_ = Kokkos::View<DataScalar*****,   DeviceType>("Intrepid2 A", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3), view.extent_int(4)); break;
          case 6: data6_ = Kokkos::View<DataScalar******,  DeviceType>("Intrepid2 A", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3), view.extent_int(4), view.extent_int(5)); break;
          case 7: data7_ = Kokkos::View<DataScalar*******, DeviceType>("Intrepid2 A", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3), view.extent_int(4), view.extent_int(5), view.extent_int(6)); break;
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
//    A(const A<DataScalar,ExecSpaceType> &data)
//    :
//    dataRank_(data.getUnderlyingViewRank()), extents_(data.getExtents()), variationType_(data.getVariationTypes()), blockPlusDiagonalLastNonDiagonal_(data.blockPlusDiagonalLastNonDiagonal()), rank_(data.rank())
//    {
//      std::cout << "Entered A copy constructor.\n";
//      if (dataRank_ != 0) // dataRank_ == 0 indicates an invalid A object (a placeholder, can indicate zero value)
//      {
//        const auto view = data.getUnderlyingView();
//        switch (dataRank_)
//        {
//          case 1: data1_ = Kokkos::View<DataScalar*,       DeviceType>("Intrepid2 A - explicit copy constructor(for debugging)", view.extent_int(0)); break;
//          case 2: data2_ = Kokkos::View<DataScalar**,      DeviceType>("Intrepid2 A - explicit copy constructor(for debugging)", view.extent_int(0), view.extent_int(1)); break;
//          case 3: data3_ = Kokkos::View<DataScalar***,     DeviceType>("Intrepid2 A - explicit copy constructor(for debugging)", view.extent_int(0), view.extent_int(1), view.extent_int(2)); break;
//          case 4: data4_ = Kokkos::View<DataScalar****,    DeviceType>("Intrepid2 A - explicit copy constructor(for debugging)", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3)); break;
//          case 5: data5_ = Kokkos::View<DataScalar*****,   DeviceType>("Intrepid2 A - explicit copy constructor(for debugging)", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3), view.extent_int(4)); break;
//          case 6: data6_ = Kokkos::View<DataScalar******,  DeviceType>("Intrepid2 A - explicit copy constructor(for debugging)", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3), view.extent_int(4), view.extent_int(5)); break;
//          case 7: data7_ = Kokkos::View<DataScalar*******, DeviceType>("Intrepid2 A - explicit copy constructor(for debugging)", view.extent_int(0), view.extent_int(1), view.extent_int(2), view.extent_int(3), view.extent_int(4), view.extent_int(5), view.extent_int(6)); break;
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
    A(ScalarView<DataScalar,DeviceType> data)
    :
    A(data,
         data.rank(),
         Kokkos::Array<int,7> {data.extent_int(0),data.extent_int(1),data.extent_int(2),data.extent_int(3),data.extent_int(4),data.extent_int(5),data.extent_int(6)},
         Kokkos::Array<DataVariationType,7> {Intrepid2::GENERAL,Intrepid2::GENERAL,Intrepid2::GENERAL,Intrepid2::GENERAL,Intrepid2::GENERAL,Intrepid2::GENERAL,Intrepid2::GENERAL}, -1)
    {}
    
    //! Constructor that accepts a DynRankView as an argument.  The data belonging to the DynRankView will be copied into a new View of matching dimensions.
    template<size_t rank, class ...DynRankViewProperties>
    A(const Kokkos::DynRankView<DataScalar,DeviceType, DynRankViewProperties...> &data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
    :
    dataRank_(data.rank()), extents_({1,1,1,1,1,1,1}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(blockPlusDiagonalLastNonDiagonal), rank_(rank)
    {
//      std::cout << "Entered a DynRankView A() constructor.\n";
      allocateAndCopyFromDynRankView(data);
      for (unsigned d=0; d<rank; d++)
      {
        extents_[d]       = extents[d];
        variationType_[d] = variationType[d];
      }
      setActiveDims();
    }
        
    template<size_t rank, class ...ViewProperties>
    A(Kokkos::View<DataScalar*,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
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
    A(Kokkos::View<DataScalar**,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
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
    A(Kokkos::View<DataScalar***,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
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
    A(Kokkos::View<DataScalar****,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
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
    A(Kokkos::View<DataScalar*****,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
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
    A(Kokkos::View<DataScalar******,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
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
    A(Kokkos::View<DataScalar*******,DeviceType, ViewProperties...> data, Kokkos::Array<int,rank> extents, Kokkos::Array<DataVariationType,rank> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
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
    A(const unsigned rank, Kokkos::View<ViewScalar,DeviceType, ViewProperties...> data, Kokkos::Array<int,7> extents, Kokkos::Array<DataVariationType,7> variationType, const int blockPlusDiagonalLastNonDiagonal = -1)
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
    A(DataScalar constantValue, Kokkos::Array<int,rank> extents)
    :
    dataRank_(1), extents_({1,1,1,1,1,1,1}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(-1), rank_(rank)
    {
      data1_ = Kokkos::View<DataScalar*,DeviceType>("Constant A",1);
      Kokkos::deep_copy(data1_, constantValue);
      for (unsigned d=0; d<rank; d++)
      {
        extents_[d] = extents[d];
      }
      setActiveDims();
    }
            
    //! default constructor (empty data)
    A()
    :
    dataRank_(0), extents_({0,0,0,0,0,0,0}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT}), blockPlusDiagonalLastNonDiagonal_(-1), rank_(0)
    {
      setActiveDims();
    }
    
  };

namespace
{



TEUCHOS_UNIT_TEST( PR14546, AllocationIssue )
{
  using DeviceType = Intrepid2::DefaultTestDeviceType;
  using BigDataArray = Kokkos::Array< A<double,DeviceType>, 20000 >;
  BigDataArray vectorComponents;
}

} // anonymous namespace
