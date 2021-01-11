// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
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
