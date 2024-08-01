// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   ETI_FunctionSpaceTools_declare.hpp
    \brief  Explicit Template Instantiation declarations for FunctionSpaceTools.  Each declaration here should be paired with a definition in ETI_FunctionSpaceTools_define.cpp.
    \author Created by N.V. Roberts.
 */

#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Intrepid2_ConfigDefs.hpp"

#include <Kokkos_Core.hpp>

namespace Intrepid2
{
//  extern template class FunctionSpaceTools<Kokkos::DefaultExecutionSpace>;
//
//  template<> extern Data<double,Kokkos::DefaultExecutionSpace>
//  FunctionSpaceTools<Kokkos::DefaultExecutionSpace>::allocateIntegralData<double>(const TransformedVectorData<double,Kokkos::DefaultExecutionSpace> vectorDataLeft,
//                                                                                  const TensorData<double,Kokkos::DefaultExecutionSpace> cellMeasures,
//                                                                                  const TransformedVectorData<double,Kokkos::DefaultExecutionSpace> vectorDataRight);
//
//  template<> extern void
//  FunctionSpaceTools<Kokkos::DefaultExecutionSpace>::integrate<double>(Data<double,Kokkos::DefaultExecutionSpace> integrals,
//                                                                       const TransformedVectorData<double,Kokkos::DefaultExecutionSpace> vectorDataLeft,
//                                                                       const TensorData<double,Kokkos::DefaultExecutionSpace> cellMeasures,
//                                                                       const TransformedVectorData<double,Kokkos::DefaultExecutionSpace> vectorDataRight);
}
