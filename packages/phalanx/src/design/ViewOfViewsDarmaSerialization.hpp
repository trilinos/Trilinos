// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __EMPIRE_VIEW_OF_VIEWS_HPP__
#define __EMPIRE_VIEW_OF_VIEWS_HPP__

#include <checkpoint/checkpoint.h>
#include <Phalanx_KokkosViewOfViews.hpp>

namespace checkpoint {
  template <typename SerializerT,typename InnerViewType,typename... OuterProps>
  void serialize(SerializerT& s, PHX::ViewOfViews3<2,InnerViewType,OuterProps...>& vov) {
    std::string label = "";
    bool is_initialized = false;
    bool is_synced = false;
    bool safety_check = true;
    size_t extent_0 = 0;
    size_t extent_1 = 0;

    if (s.isSizing() || s.isPacking() || s.isFootprinting()) {
      auto host_view = vov.getViewHost();
      label = std::string(host_view.label());
      is_initialized = vov.isInitialized();
      is_synced = vov.deviceViewIsSynced();
      safety_check = vov.safetyCheck();
      extent_0 = host_view.extent(0);
      extent_1 = host_view.extent(1);
    }

    s | label
      | is_initialized
      | is_synced
      | safety_check
      | extent_0
      | extent_1;

    if (!safety_check)
      vov.disableSafetyCheck();

    // If not initialized, then the inner views don't exist and we can
    // completely ignore inner view data.
    if (is_initialized) {

      if (s.isSizing() || s.isPacking() || s.isFootprinting()) {
        auto host_view = vov.getViewHost();
        for (size_t i=0; i < extent_0; ++i) {
          for (size_t j=0; j < extent_1; ++j) {
            InnerViewType tmp = host_view(i,j);
            TEUCHOS_ASSERT(tmp.extent(0)==10);
            s | tmp;
          }
        }
      }
      else if (s.isUnpacking()) {
        vov.initialize(label,extent_0,extent_1);
        InnerViewType tmp;
        for (size_t i=0; i < extent_0; ++i) {
          for (size_t j=0; j < extent_1; ++j) {
            s | tmp;
            vov.addView(tmp,i,j);
          }
        }
        if (is_synced)
          vov.syncHostToDevice();
      }
    }

  }
}

#endif

TEUCHOS_UNIT_TEST( Serialization, PHX_ViewOfViews3 ) {

  using InnerViewType = Kokkos::View<double*>;
  using VoV = PHX::ViewOfViews3<2,InnerViewType,empire::Device::memory_space>;
  const int inner_size = 10;

  VoV vov("vov serialization test",2,2);
  {
    InnerViewType a("a",inner_size);
    InnerViewType b("b",inner_size);
    InnerViewType c("c",inner_size);
    InnerViewType d("d",inner_size);
    Kokkos::deep_copy(a,2.0);
    Kokkos::deep_copy(b,3.0);
    Kokkos::deep_copy(c,4.0);

    vov.addView(a,0,0);
    vov.addView(b,1,0);
    vov.addView(c,0,1);
    vov.addView(d,1,1);
    vov.syncHostToDevice();
  }
  auto dev = vov.getViewDevice();
  Kokkos::parallel_for("vov serialization",inner_size,KOKKOS_LAMBDA(const int i){
      dev(1,1)(i) = dev(0,0)(i) * dev(1,0)(i) + dev(0,1)(i);
    });

  auto packed_source = checkpoint::serialize(vov);

  VoV vov2;
  ::checkpoint::deserializeInPlace<VoV>( packed_source->getBuffer(), &vov2 );

  auto vov_host = PHX::createHostHostViewOfViews(vov2.getViewHost());

  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
  for (int i=0; i < inner_size; ++i) {
    double val = vov_host(1,1)(i);
    TEST_FLOATING_EQUALITY(val,10.0,tol);
  }
}
