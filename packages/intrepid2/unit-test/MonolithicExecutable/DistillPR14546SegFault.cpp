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
    
    Kokkos::Array<ordinal_type,7> activeDims_;
    int numActiveDims_; // how many of the 7 entries are actually filled in
    
    //! class initialization method.  Called by constructors.
    void setActiveDims()
    {}

    public:
  public:
            
    //! default constructor (empty data)
    Data()
    :
    dataRank_(0), extents_({0,0,0,0,0,0,0}), variationType_({Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT,Intrepid2::CONSTANT})
    {
      setActiveDims();
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
