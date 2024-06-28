// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   ETI_FunctionSpaceTools_define.cpp
    \brief  Explicit Template Instantiation for FunctionSpaceTools.  Each definition here should be paired with a declaration in ETI_FunctionSpaceTools_declare.hpp.
    \author Created by N.V. Roberts.
 */

#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Intrepid2_ConfigDefs.hpp"

#include <Kokkos_Core.hpp>

namespace Intrepid2
{
//  using ExecSpaceType = Kokkos::DefaultExecutionSpace;
//
//  template class FunctionSpaceTools<ExecSpaceType>;
//
//  template<> Data<double,ExecSpaceType>
//  FunctionSpaceTools<ExecSpaceType>::allocateIntegralData<double>(const TransformedVectorData<double,ExecSpaceType> vectorDataLeft,
//                                                                  const TensorData<double,ExecSpaceType> cellMeasures,
//                                                                  const TransformedVectorData<double,ExecSpaceType> vectorDataRight);
//
//  template<> void
//  FunctionSpaceTools<ExecSpaceType>::integrate<double>(Data<double,ExecSpaceType> integrals,
//                                                       const TransformedVectorData<double,ExecSpaceType> vectorDataLeft,
//                                                       const TensorData<double,ExecSpaceType> cellMeasures,
//                                                       const TransformedVectorData<double,ExecSpaceType> vectorDataRight);
//
//#ifdef HAVE_INTREPID2_SACADO
//
//#endif
}
