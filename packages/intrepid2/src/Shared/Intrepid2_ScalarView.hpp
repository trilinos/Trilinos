// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  Intrepid2_ScalarView.hpp
//  QuadraturePerformance
//
//  Created by Roberts, Nathan V on 8/24/20.
//

#ifndef Intrepid2_ScalarView_h
#define Intrepid2_ScalarView_h

#include <Kokkos_DynRankView.hpp>
#include "Intrepid2_Utils.hpp"

namespace Intrepid2
{
  //! Intrepid2::ScalarView template requires execution space argument, in contrast to Kokkos::DynRankView, which defaults to the default execution space.  (We allow users to set the execution space for many of our classes; by using this template we can avoid accidentally using the default in members of those classes.
  template<typename Scalar, typename ExecSpaceType>
  using ScalarView = Kokkos::DynRankView<Scalar, ExecSpaceType>;
} // end namespace Intrepid2

#endif /* Intrepid2_ScalarView_h */
