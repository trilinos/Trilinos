// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 *  Copyright 2008-2009 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include <cusp/array1d.h>
#include <cusp/detail/device/MVmultiply.h>

namespace cusp
{
namespace detail
{
namespace dispatch
{
//////////////////
// Device Paths //
//////////////////
template <typename LinearOperator,
          typename MatrixOrVector1,
          typename MatrixOrVector2>
void MVmultiply(const LinearOperator&  A,
              const MatrixOrVector1& B,
                    MatrixOrVector2& C,
              cusp::device_memory,
              cusp::device_memory,
              cusp::device_memory)
{
    cusp::detail::device::MVmultiply(A, B, C);
}

template <typename LinearOperator,
          typename MatrixOrVector1,
          typename MatrixOrVector2>
void OVmultiply(const LinearOperator&  A,
              const MatrixOrVector1& B,
                    MatrixOrVector2& C,
              cusp::device_memory,
              cusp::device_memory,
              cusp::device_memory)
{
    cusp::detail::device::OVmultiply(A, B, C);
}

template <typename MatrixOrVector,
          typename MatrixOrVector1,
          typename MatrixOrVector2>
void MVdot(const MatrixOrVector&  A,
            const  MatrixOrVector1& B,
                    MatrixOrVector2& C,
              cusp::device_memory,
              cusp::device_memory,
              cusp::device_memory)
{
    cusp::detail::device::MVdot(A, B, C);
}

template <typename ValueType,
          typename MatrixOrVector1,
          typename MatrixOrVector2>
void axpby_array(const ValueType&  A,
            const  MatrixOrVector1& X,
            const ValueType&  B,
            const  MatrixOrVector1& Y,
                  MatrixOrVector2& Z,
              cusp::device_memory,
              cusp::device_memory)
{
    cusp::detail::device::axpby(A, X, B, Y, Z);
}


} // end namespace dispatch
} // end namespace detail
} // end namespace cusp
