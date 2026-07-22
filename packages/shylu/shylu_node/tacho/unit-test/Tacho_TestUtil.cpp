// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include "Tacho.hpp"
#include "Tacho_Util.hpp"
#include "Tacho_MatrixMarket.hpp"

using host_device_type = typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type;
using device_type = typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::type;

TEST( Util, is_complex_type ) {
  EXPECT_FALSE(int(Tacho::ArithTraits<double>::is_complex));
  EXPECT_TRUE(int(Tacho::ArithTraits<std::complex<double> >::is_complex));
  EXPECT_TRUE(int(Tacho::ArithTraits<Kokkos::complex<double> >::is_complex));
}


TEST( Util, tag ) {
  using Tacho::NullTag;
  using Tacho::Partition;
  EXPECT_TRUE(int(Tacho::is_valid_partition_tag<Partition::Top>::value));
  EXPECT_TRUE(int(Tacho::is_valid_partition_tag<Partition::Bottom>::value));
  EXPECT_TRUE(int(Tacho::is_valid_partition_tag<Partition::Left>::value));
  EXPECT_TRUE(int(Tacho::is_valid_partition_tag<Partition::Right>::value));
  EXPECT_TRUE(int(Tacho::is_valid_partition_tag<Partition::TopLeft>::value));
  EXPECT_TRUE(int(Tacho::is_valid_partition_tag<Partition::TopRight>::value));
  EXPECT_TRUE(int(Tacho::is_valid_partition_tag<Partition::BottomLeft>::value));
  EXPECT_TRUE(int(Tacho::is_valid_partition_tag<Partition::BottomRight>::value));

  using Tacho::Uplo;
  EXPECT_TRUE(int(Tacho::is_valid_uplo_tag<Uplo::Upper>::value));
  EXPECT_TRUE(int(Tacho::is_valid_uplo_tag<Uplo::Lower>::value));
  EXPECT_FALSE(int(Tacho::is_valid_uplo_tag<NullTag>::value));

  using Tacho::Side;
  EXPECT_TRUE(int(Tacho::is_valid_side_tag<Tacho::Side::Left>::value));
  EXPECT_TRUE(int(Tacho::is_valid_side_tag<Tacho::Side::Right>::value));
  EXPECT_FALSE(int(Tacho::is_valid_side_tag<NullTag>::value));

  using Tacho::Diag;
  EXPECT_TRUE(int(Tacho::is_valid_diag_tag<Tacho::Diag::Unit>::value));
  EXPECT_TRUE(int(Tacho::is_valid_diag_tag<Tacho::Diag::NonUnit>::value));
  EXPECT_FALSE(int(Tacho::is_valid_diag_tag<NullTag>::value));

  using Tacho::Trans;
  EXPECT_TRUE(int(Tacho::is_valid_trans_tag<Tacho::Trans::Transpose>::value));
  EXPECT_TRUE(int(Tacho::is_valid_trans_tag<Tacho::Trans::ConjTranspose>::value));
  EXPECT_TRUE(int(Tacho::is_valid_trans_tag<Tacho::Trans::NoTranspose>::value));
  EXPECT_FALSE(int(Tacho::is_valid_trans_tag<NullTag>::value));
}


int main (int argc, char *argv[]) {

  int r_val(0);
  Kokkos::initialize(argc, argv);
  {
     const bool detail = false;
     Tacho::printExecSpaceConfiguration<typename device_type::execution_space>("DeviceSpace", detail);
     Tacho::printExecSpaceConfiguration<typename host_device_type::execution_space>("HostSpace", detail);
  
     ::testing::InitGoogleTest(&argc, argv);

     r_val = RUN_ALL_TESTS();
  }
  Kokkos::finalize();
  return r_val;
}
