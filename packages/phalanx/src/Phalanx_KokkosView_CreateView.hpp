// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_KOKKOS_VIEW_CREATE_VIEW_HPP
#define PHALANX_KOKKOS_VIEW_CREATE_VIEW_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_KokkosView_HiddenDimensionForSFINAE.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Sacado.hpp"
#include <any>
#include <vector>

namespace PHX {

  enum class ViewCreationMode {AllocateMemory,UseTracker};
  
  /** \brief Creates a view by allocating new memory or binding to an existing tracker.
   *
   *  NOTE: This function operates in two modes set by the mode input
   *  parameter. In AllocateMemory mode, it allocates new memory while
   *  in UseTracker mode it reuses the memory out of the tracker. In
   *  AllocateMemory mode, it also returns the newly allocated tracker
   *  in the tracker parameter. In UseTracker mode it expects the
   *  tracker parameter to contain a valid tracker that is large
   *  enough to hold the view.
   *
   *  \param[in] mode If set to AllocateMemory, then new memory is allocated. If set to UseTracker, then it uses the memory supplied in the tracker. 
   *  \param[in] t FieldTag with size information.
   *  \param[in or out depending on mode] Contains teh tracked associated with the created field. In AllocateMemory, the tracker is in output parameter. In UseTracker, the tracker is an input parameter.
   *  \param[in] derivative_dimensions Hidden dimension sizes required by some scalar types (e.g. in AD scalar types, the derivative length.
   * 
   *  NOTE: To get the map impl with correct padded sizes for kokkos
   *  views, we actually have to allocate the view and then call
   *  span(). It's inefficient, but this only happens during
   *  setup. For more details, see issue:
   *  https://github.com/kokkos/kokkos/issues/2182
   */
  template<typename ScalarT,typename Layout,typename Device>
  typename std::enable_if<PHX::requires_dynamic_hidden_dimension<ScalarT>::value,std::any>::type
  createView(const PHX::ViewCreationMode& mode,
             const PHX::FieldTag& t,
             const std::vector<PHX::index_size_type>& derivative_dimensions,
             Kokkos::Impl::SharedAllocationTracker& tracker) {
    std::any a;

    const PHX::DataLayout& dl = t.dataLayout();
    // DFad type contains a hidden dimension of the size of the number
    // of derivatives.  We add one to this for tracking the value of
    // the actual function.
    const PHX::index_size_type hDim = derivative_dimensions[0] + 1;

    using MemoryType = typename Sacado::ValueType<ScalarT>::type;

    if (dl.rank() == 1) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT*,Layout,Device>(t.identifier(),
                                                      dl.dimension(0),
                                                      hDim);
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT*,Layout,Device>(ptr,
                                                 dl.dimension(0),
                                                 hDim);
      }
    }
    else if (dl.rank() == 2) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT**,Layout,Device>(t.identifier(),
                                                       dl.dimension(0),
                                                       dl.dimension(1),
                                                       hDim);
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT**,Layout,Device>(ptr,
                                                  dl.dimension(0),
                                                  dl.dimension(1),
                                                  hDim);
      }
    }
    else if (dl.rank() == 3) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT***,Layout,Device>(t.identifier(),
                                                        dl.dimension(0),
                                                        dl.dimension(1),
                                                        dl.dimension(2),
                                                        hDim);
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT***,Layout,Device>(ptr,
                                                   dl.dimension(0),
                                                   dl.dimension(1),
                                                   dl.dimension(2),
                                                   hDim);
      }
    }
    else if (dl.rank() == 4) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT****,Layout,Device>(t.identifier(),
                                                         dl.dimension(0),
                                                         dl.dimension(1),
                                                         dl.dimension(2),
                                                         dl.dimension(3),
                                                         hDim);
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT****,Layout,Device>(ptr,
                                                    dl.dimension(0),
                                                    dl.dimension(1),
                                                    dl.dimension(2),
                                                    dl.dimension(3),
                                                    hDim);
      }
    }
    else if (dl.rank() == 5) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT*****,Layout,Device>(t.identifier(),
                                                          dl.dimension(0),
                                                          dl.dimension(1),
                                                          dl.dimension(2),
                                                          dl.dimension(3),
                                                          dl.dimension(4),
                                                          hDim);
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT*****,Layout,Device>(ptr,
                                                     dl.dimension(0),
                                                     dl.dimension(1),
                                                     dl.dimension(2),
                                                     dl.dimension(3),
                                                     dl.dimension(4),
                                                     hDim);
      }
    }
    else if (dl.rank() == 6) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT******,Layout,Device>(t.identifier(),
                                                           dl.dimension(0),
                                                           dl.dimension(1),
                                                           dl.dimension(2),
                                                           dl.dimension(3),
                                                           dl.dimension(4),
                                                           dl.dimension(5),
                                                           hDim);
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT******,Layout,Device>(ptr,
                                                      dl.dimension(0),
                                                      dl.dimension(1),
                                                      dl.dimension(2),
                                                      dl.dimension(3),
                                                      dl.dimension(4),
                                                      dl.dimension(5),
                                                      hDim);
      }
    }
    else if (dl.rank() == 7) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT*******,Layout,Device>(t.identifier(),
                                                            dl.dimension(0),
                                                            dl.dimension(1),
                                                            dl.dimension(2),
                                                            dl.dimension(3),
                                                            dl.dimension(4),
                                                            dl.dimension(5),
                                                            dl.dimension(6),
                                                            hDim);
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT*******,Layout,Device>(ptr,
                                                       dl.dimension(0),
                                                       dl.dimension(1),
                                                       dl.dimension(2),
                                                       dl.dimension(3),
                                                       dl.dimension(4),
                                                       dl.dimension(5),
                                                       dl.dimension(6),
                                                       hDim);
      }
    }
    
    return a;
  }

  /** \brief Creates a view by allocating new memory or binding to an existing tracker.
   *
   *  NOTE: This function operates in two modes set by the mode input
   *  parameter. In AllocateMemory mode, it allocates new memory while
   *  in UseTracker mode it reuses the memory out of the tracker. In
   *  AllocateMemory mode, it also returns the newly allocated tracker
   *  in the tracker parameter. In UseTracker mode it expects the
   *  tracker parameter to contain a valid tracker that is large
   *  enough to hold the view.
   *
   *  \param[in] mode If set to AllocateMemory, then new memory is allocated. If set to UseTracker, then it uses the memory supplied in the tracker. 
   *  \param[in] t FieldTag with size information.
   *  \param[in or out depending on mode] Contains teh tracked associated with the created field. In AllocateMemory, the tracker is in output parameter. In UseTracker, the tracker is an input parameter.
   *  \param[in] derivative_dimensions Hidden dimension sizes required by some scalar types (e.g. in AD scalar types, the derivative length.
   * 
   *  NOTE: To get the map impl with correct padded sizes for kokkos
   *  views, we actually have to allocate the view and then call
   *  span(). It's inefficient, but this only happens during
   *  setup. For more details, see issue:
   *  https://github.com/kokkos/kokkos/issues/2182
   */
  template<typename ScalarT,typename Layout,typename Device>
  typename std::enable_if<!PHX::requires_dynamic_hidden_dimension<ScalarT>::value,std::any>::type
  createView(const PHX::ViewCreationMode& mode,
             const PHX::FieldTag& t,
             const std::vector<PHX::index_size_type>& ,
             Kokkos::Impl::SharedAllocationTracker& tracker) {
    std::any a;
    const PHX::DataLayout& dl = t.dataLayout();

    using MemoryType = typename Sacado::ValueType<ScalarT>::type;

    if (dl.rank() == 1) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT*,Layout,Device>(t.identifier(),
                                                      dl.dimension(0));
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT*,Layout,Device>(ptr,
                                                 dl.dimension(0));
      }
    }
    else if (dl.rank() == 2) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT**,Layout,Device>(t.identifier(),
                                                       dl.dimension(0),
                                                       dl.dimension(1));
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT**,Layout,Device>(ptr,
                                                  dl.dimension(0),
                                                  dl.dimension(1));
      }
    }
    else if (dl.rank() == 3) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT***,Layout,Device>(t.identifier(),
                                                        dl.dimension(0),
                                                        dl.dimension(1),
                                                        dl.dimension(2));
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT***,Layout,Device>(ptr,
                                                   dl.dimension(0),
                                                   dl.dimension(1),
                                                   dl.dimension(2));
      }
    }
    else if (dl.rank() == 4) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT****,Layout,Device>(t.identifier(),
                                                         dl.dimension(0),
                                                         dl.dimension(1),
                                                         dl.dimension(2),
                                                         dl.dimension(3));
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT****,Layout,Device>(ptr,
                                                    dl.dimension(0),
                                                    dl.dimension(1),
                                                    dl.dimension(2),
                                                    dl.dimension(3));
      }
    }
    else if (dl.rank() == 5) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT*****,Layout,Device>(t.identifier(),
                                                          dl.dimension(0),
                                                          dl.dimension(1),
                                                          dl.dimension(2),
                                                          dl.dimension(3),
                                                          dl.dimension(4));
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT*****,Layout,Device>(ptr,
                                                     dl.dimension(0),
                                                     dl.dimension(1),
                                                     dl.dimension(2),
                                                     dl.dimension(3),
                                                     dl.dimension(4));
      }
    }
    else if (dl.rank() == 6) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT******,Layout,Device>(t.identifier(),
                                                           dl.dimension(0),
                                                           dl.dimension(1),
                                                           dl.dimension(2),
                                                           dl.dimension(3),
                                                           dl.dimension(4),
                                                           dl.dimension(5));
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT******,Layout,Device>(ptr,
                                                      dl.dimension(0),
                                                      dl.dimension(1),
                                                      dl.dimension(2),
                                                      dl.dimension(3),
                                                      dl.dimension(4),
                                                      dl.dimension(5));
      }
    }
    else if (dl.rank() == 7) {
      if (mode == ViewCreationMode::AllocateMemory) {
        auto v = Kokkos::View<ScalarT*******,Layout,Device>(t.identifier(),
                                                            dl.dimension(0),
                                                            dl.dimension(1),
                                                            dl.dimension(2),
                                                            dl.dimension(3),
                                                            dl.dimension(4),
                                                            dl.dimension(5),
                                                            dl.dimension(6));
        tracker = v.impl_track();
        a = v;
      }
      else {
        MemoryType* ptr = reinterpret_cast<MemoryType*>(tracker.template get_record<Device>()->data());
        a = Kokkos::View<ScalarT*******,Layout,Device>(ptr,
                                                       dl.dimension(0),
                                                       dl.dimension(1),
                                                       dl.dimension(2),
                                                       dl.dimension(3),
                                                       dl.dimension(4),
                                                       dl.dimension(5),
                                                       dl.dimension(6));
      }
    }

    return a;
  }
}

#endif
