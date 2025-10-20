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
    
    bool underlyingMatchesLogical_;   // if true, this A object has the same rank and extent as the underlying view
    Kokkos::Array<ordinal_type,7> activeDims_;
    int numActiveDims_; // how many of the 7 entries are actually filled in
    
    ordinal_type rank_;
  public:
            
    //! default constructor (empty data)
    A()
    :
    dataRank_(0), rank_(0)
    {}
    
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
