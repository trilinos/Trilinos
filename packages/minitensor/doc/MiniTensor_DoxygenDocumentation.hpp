// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINITENSOR_DOXYGEN_DOCUMENTATION_HPP
#define MINITENSOR_DOXYGEN_DOCUMENTATION_HPP

/*!

\mainpage

\section minitensor_index Table of Contents

- \ref minitensor_overview
- \ref minitensor_modules
- \ref minitensor_usage
- \ref minitensor_copyright
- \ref minitensor_questions

\section minitensor_overview MiniTensor Overview

MiniTensor is a library for the use, manipulation, algebra and optimization
of small vectors and tensors and problems that depend on them. Its purpose
is to provide a compact representation of vector and tensor expressions.
Its emphasis is on ease of use, and accurate algorithms, specifically those
used for the development of constitutive models in finite deformation solid
mechanics. It is fully templated, compatible with Sacado automatic
differentiation types, and usable inside Kokkos device kernels.

MiniTensor is part of the
<a href="https://trilinos.github.io">Trilinos Project</a>.

\section minitensor_modules Modules

MiniTensor is organized as a collection of module headers. Including
MiniTensor.h provides the containers plus linear algebra, geometry and
mechanics; each module header can also be included on its own.

| Header | Contents |
|---|---|
| MiniTensor.h | Umbrella: containers, linear algebra, geometry, mechanics |
| MiniTensor_Vector.h, MiniTensor_Tensor.h, MiniTensor_Tensor3.h, MiniTensor_Tensor4.h, MiniTensor_Matrix.h | \ref minitensor_containers |
| MiniTensor_Traits.h, MiniTensor_Scalar.h | \ref minitensor_traits |
| MiniTensor_Norms.h | \ref minitensor_norms |
| MiniTensor_Inverse.h | \ref minitensor_inverse |
| MiniTensor_Factorizations.h | \ref minitensor_factorizations |
| MiniTensor_MatrixFunctions.h | \ref minitensor_matrix_functions |
| MiniTensor_Rotations.h | \ref minitensor_rotations |
| MiniTensor_LinearAlgebra.h | Umbrella over the five linear-algebra modules above |
| MiniTensor_Quaternion.h | Quaternions and rotation-representation conversions: \ref minitensor_rotations |
| MiniTensor_Geometry.h | \ref minitensor_geometry |
| MiniTensor_Mechanics.h | \ref minitensor_mechanics |
| MiniTensor_Solvers.h, MiniTensor_TestFunctions.h | \ref minitensor_solvers |

\section minitensor_usage Basic Usage

\code{.cpp}
#include <MiniTensor.h>

using T = double;

// 3D identity tensor, static storage.
minitensor::Tensor<T, 3> const I = minitensor::identity<T, 3>(3);

// A deformation-gradient-like tensor; fillers provide zeros, ones,
// random values, ...
minitensor::Tensor<T, 3> F(3, minitensor::Filler::RANDOM_UNIFORM);
F = 0.5 * (F + minitensor::transpose(F)) + I;

// Linear algebra and mechanics operations.
T const j = minitensor::det(F);
minitensor::Tensor<T, 3> const C = minitensor::transpose(F) * F;
minitensor::Tensor<T, 3> const E = 0.5 * (C - I);
minitensor::Tensor<T, 3> const logC = minitensor::log_sym(C);
\endcode

Replace \c double with a Sacado forward-AD type to differentiate through
any of the above.

\section minitensor_copyright Copyright and License

MiniTensor is distributed under the BSD 3-clause license; see the
COPYRIGHT and LICENSE files in the package root, and the Trilinos
<a href="https://trilinos.github.io/about.html#license-and-copyright">
License and Copyright</a> page.

\section minitensor_questions For All Questions and Comments...

   Please contact the authors listed in the License above,
   or open an issue in the Trilinos GitHub repository
   (https://github.com/trilinos/Trilinos).

*/

/*!
\defgroup minitensor_containers Containers
\brief Vectors, second- to fourth-order tensors, and general matrices,
with static (compile-time) or dynamic storage.
*/

/*!
\defgroup minitensor_traits Type Traits and Scalar Utilities
\brief Sacado-aware type promotion machinery, index types, and
scalar-level helpers such as machine_epsilon and not_a_number.
*/

/*!
\defgroup minitensor_norms Norms and Invariants
\brief Scalar measures of tensors: norms, determinant, trace, principal
invariants, and argmax queries.
*/

/*!
\defgroup minitensor_inverse Inverse and Linear Solves
\brief Tensor inversion, row/column manipulation, rank-one updates,
preconditioners and small dense linear solves.
*/

/*!
\defgroup minitensor_factorizations Factorizations
\brief Givens rotations, Cholesky, symmetric eigendecomposition, SVD,
polar decompositions and condition numbers.
*/

/*!
\defgroup minitensor_matrix_functions Matrix Functions
\brief Tensor exponential, logarithm and square root families, integer
powers, and the Baker-Campbell-Hausdorff series.
*/

/*!
\defgroup minitensor_rotations Rotations
\brief Logarithmic and exponential maps on SO(N), quaternions, and
conversions among rotation representations: rotation matrix, unit
quaternion and principal rotation vector.
*/

/*!
\defgroup minitensor_geometry Geometry
\brief Lengths, areas, volumes and centroids of finite-element shapes,
and related parametrization utilities.
*/

/*!
\defgroup minitensor_mechanics Mechanics
\brief Continuum-mechanics operations: deformation and stress measures,
push-forward and pull-back operations.
*/

/*!
\defgroup minitensor_solvers Solvers and Optimization
\brief Nonlinear system solvers, unconstrained and constrained
optimization methods, and benchmark test functions.
*/

#endif //ifndef MINITENSOR_DOXYGEN_DOCUMENTATION_HPP
