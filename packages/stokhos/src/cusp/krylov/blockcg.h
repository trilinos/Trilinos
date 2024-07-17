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

/*! \file cg.h
 *  \brief Conjugate Gradient (CG) method
 */

#pragma once

#include <cusp/detail/config.h>

namespace cusp
{
namespace krylov
{


/*! \p cg : Conjugate Gradient method
 *
 * Solves the symmetric, positive-definite linear system A x = b with multiple right hand sides
 * using the default convergence criteria.
 */
template <class LinearOperator,
          class Vector>
void blockcg(LinearOperator& A,
        Vector& x,
        Vector& b);

/*! \p cg : Conjugate Gradient method
 *
 * Solves the symmetric, positive-definite linear system A x = b without preconditioning.
 */
template <class LinearOperator,
          class Vector,
          class Monitor>
void blockcg(LinearOperator& A,
        Vector& x,
        Vector& b,
        Monitor& monitor);



/*! \p cg : Conjugate Gradient method
 *
 * Solves the symmetric, positive-definite linear system A x = b
 * with preconditioner \p M.
 *
 *
 */
template <class LinearOperator,
          class Vector,
          class Monitor,
          class Preconditioner>
void blockcg(LinearOperator& A,
        Vector& x,
        Vector& b,
        Monitor& monitor,
        Preconditioner& M);
/*! \}
 */

} // end namespace krylov
} // end namespace cusp

#include <cusp/krylov/blockcg.inl>
