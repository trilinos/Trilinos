// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_KOKKOS_VIEW_ALLOCATION_SIZE_HPP
#define PHALANX_KOKKOS_VIEW_ALLOCATION_SIZE_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_KokkosView_HiddenDimensionForSFINAE.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Sacado.hpp"
#include <any>
#include <vector>

namespace PHX {

  /** \brief Returns the allocation size in bytes for a particular kokkos view.
   *
   *  NOTE: To get the padded sizes for kokkos views, we actually have
   *  to allocate the view and then call span(). It's inefficient, but
   *  this only happens during setup. For more details, see issue:
   *  https://github.com/kokkos/kokkos/issues/2182
   */
  template<typename ScalarT,typename Layout,typename Device>
  typename std::enable_if<PHX::requires_dynamic_hidden_dimension<ScalarT>::value,std::size_t>::type
  getAllocationSize(const PHX::FieldTag& t,const std::vector<PHX::index_size_type>& derivative_dimensions) {
    std::size_t s = 0;

    const PHX::DataLayout& dl = t.dataLayout();
    // DFad type contains a hidden dimension of the size of the number
    // of derivatives.  We add one to this for tracking the value of
    // the actual function.
    const PHX::index_size_type hDim = derivative_dimensions[0] + 1;

    if (dl.rank() == 1)
      s = Kokkos::View<ScalarT*,Layout,Device>(t.identifier(),
                                               dl.dimension(0),
                                               hDim).impl_track().template get_record<Device>()->size();
    else if (dl.rank() == 2)
      s = Kokkos::View<ScalarT**,Layout,Device>(t.identifier(),
                                                dl.dimension(0),
                                                dl.dimension(1),
                                                hDim).impl_track().template get_record<Device>()->size();
    else if (dl.rank() == 3)
      s = Kokkos::View<ScalarT***,Layout,Device>(t.identifier(),
                                                 dl.dimension(0),
                                                 dl.dimension(1),
                                                 dl.dimension(2),
                                                 hDim).impl_track().template get_record<Device>()->size();
    else if (dl.rank() == 4)
      s = Kokkos::View<ScalarT****,Layout,Device>(t.identifier(),
                                                  dl.dimension(0),
                                                  dl.dimension(1),
                                                  dl.dimension(2),
                                                  dl.dimension(3),
                                                  hDim).impl_track().template get_record<Device>()->size();
    else if (dl.rank() == 5)
      s = Kokkos::View<ScalarT*****,Layout,Device>(t.identifier(),
                                                   dl.dimension(0),
                                                   dl.dimension(1),
                                                   dl.dimension(2),
                                                   dl.dimension(3),
                                                   dl.dimension(4),
                                                   hDim).impl_track().template get_record<Device>()->size();
    else if (dl.rank() == 6)
      s = Kokkos::View<ScalarT******,Layout,Device>(t.identifier(),
                                                    dl.dimension(0),
                                                    dl.dimension(1),
                                                    dl.dimension(2),
                                                    dl.dimension(3),
                                                    dl.dimension(4),
                                                    dl.dimension(5),
                                                    hDim).impl_track().template get_record<Device>()->size();
    else if (dl.rank() == 7)
      s = Kokkos::View<ScalarT*******,Layout,Device>(t.identifier(),
                                                     dl.dimension(0),
                                                     dl.dimension(1),
                                                     dl.dimension(2),
                                                     dl.dimension(3),
                                                     dl.dimension(4),
                                                     dl.dimension(5),
                                                     dl.dimension(6),
                                                     hDim).impl_track().template get_record<Device>()->size();
    
    return s;
  }

  template<typename ScalarT,typename Layout,typename Device>
  typename std::enable_if<!PHX::requires_dynamic_hidden_dimension<ScalarT>::value,std::size_t>::type
  getAllocationSize(const PHX::FieldTag& t,const std::vector<PHX::index_size_type>& ) {
    std::size_t s = 0;
    const PHX::DataLayout& dl = t.dataLayout();

    if (dl.rank() == 1)
      s = Kokkos::View<ScalarT*,Layout,Device>(t.identifier(),
                                               dl.dimension(0)).impl_track().template get_record<Device>()->size();
    else if (dl.rank() == 2)
      s = Kokkos::View<ScalarT**,Layout,Device>(t.identifier(),
                                                dl.dimension(0),
                                                dl.dimension(1)).impl_track().template get_record<Device>()->size();
    else if (dl.rank() == 3)
      s = Kokkos::View<ScalarT***,Layout,Device>(t.identifier(),
                                                 dl.dimension(0),
                                                 dl.dimension(1),
                                                 dl.dimension(2)).impl_track().template get_record<Device>()->size();
    else if (dl.rank() == 4)
      s = Kokkos::View<ScalarT****,Layout,Device>(t.identifier(),
                                                  dl.dimension(0),
                                                  dl.dimension(1),
                                                  dl.dimension(2),
                                                  dl.dimension(3)).impl_track().template get_record<Device>()->size();
    else if (dl.rank() == 5)
      s = Kokkos::View<ScalarT*****,Layout,Device>(t.identifier(),
                                                   dl.dimension(0),
                                                   dl.dimension(1),
                                                   dl.dimension(2),
                                                   dl.dimension(3),
                                                   dl.dimension(4)).impl_track().template get_record<Device>()->size();
    else if (dl.rank() == 6)
      s = Kokkos::View<ScalarT******,Layout,Device>(t.identifier(),
                                                    dl.dimension(0),
                                                    dl.dimension(1),
                                                    dl.dimension(2),
                                                    dl.dimension(3),
                                                    dl.dimension(4),
                                                    dl.dimension(5)).impl_track().template get_record<Device>()->size();
    else if (dl.rank() == 7)
      s = Kokkos::View<ScalarT*******,Layout,Device>(t.identifier(),
                                                     dl.dimension(0),
                                                     dl.dimension(1),
                                                     dl.dimension(2),
                                                     dl.dimension(3),
                                                     dl.dimension(4),
                                                     dl.dimension(5),
                                                     dl.dimension(6)).impl_track().template get_record<Device>()->size();
    return s;
  }
}

#endif
