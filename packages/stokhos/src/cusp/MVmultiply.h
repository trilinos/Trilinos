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

/*! \file multiply.h
 *  \brief Matrix multiplication
 */

#pragma once

#include <cusp/detail/config.h>

namespace cusp
{

template <typename LinearOperator,
          typename MatrixOrVector1,
          typename MatrixOrVector2>
void MVmultiply(LinearOperator&  A,
                MatrixOrVector1& B,
                MatrixOrVector2& C);

template <typename LinearOperator,
          typename MatrixOrVector1,
          typename MatrixOrVector2>
void OVmultiply(LinearOperator&  A,
                MatrixOrVector1& B,
                MatrixOrVector2& C);


template <typename MatrixOrVector,
          typename MatrixOrVector1,
          typename MatrixOrVector2>
void MVdot(const MatrixOrVector&  A,
              const MatrixOrVector1& B,
              MatrixOrVector2& C);

template <typename ValueType,
          typename MatrixOrVector1,
          typename MatrixOrVector2>
void axpby_array(const ValueType&  a,
                 const MatrixOrVector1& x,
                 const ValueType&  b,
                 const MatrixOrVector1& y,
                 MatrixOrVector2& z);



/*! \}
 */

} // end namespace cusp

#include <cusp/detail/MVmultiply.inl>
