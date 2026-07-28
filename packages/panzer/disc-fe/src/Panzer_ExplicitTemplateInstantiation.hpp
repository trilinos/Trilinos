// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION_HPP
#define PANZER_EXPLICIT_TEMPLATE_INSTANTIATION_HPP

#include "Panzer_Traits.hpp"

/**
 * \file Panzer_ExplicitTemplateInstantiation.hpp
 * \brief Convenience macros for explicitly instantiating panzer template classes over all enabled evaluation types.
 *
 * Panzer evaluators, responses, and similar classes are templated on
 * an evaluation type `EvalT` (see panzer::Traits::EvalTypes: Residual,
 * Jacobian, Tangent, and optionally Hessian). Rather than letting every
 * translation unit that uses such a class implicitly instantiate and
 * compile it, each class is explicitly instantiated exactly once, in
 * its own `.cpp` file, using the macros below. This avoids redundant
 * template instantiation across the library and significantly reduces
 * build times.
 *
 * The macros are grouped by how many template arguments the class
 * takes in addition to (or instead of) the evaluation type itself:
 *  - `_ONE_T`      : `name<EvalT>`
 *  - `_TWO_T`      : `name<EvalT, panzer::Traits>`
 *  - `_THREE_T`    : `name<EvalT, panzer::Traits, ExtraT>`
 *  - `_THREE_2U_T` : `name<EvalT, FirstExtraT, SecondExtraT>` (no `panzer::Traits`)
 *  - `_FOUR_T`     : `name<EvalT, panzer::Traits, FirstExtraT, SecondExtraT>`
 *
 * Within each group, the top-level `PANZER_INSTANTIATE_TEMPLATE_CLASS_*`
 * macro is the one client `.cpp` files invoke; it expands to one
 * `template class ...;` instantiation per evaluation type, built out of
 * the per-evaluation-type `_RESIDUAL_`, `_TANGENT_`, `_JACOBIAN_`, and
 * `_HESSIAN_` macros. The `_HESSIAN_` variants expand to nothing unless
 * `Panzer_BUILD_HESSIAN_SUPPORT` is defined, so Hessian instantiations
 * are silently skipped in builds without Hessian support.
 *
 * Example: `PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(panzer::MyEvaluator)`
 * in `Panzer_MyEvaluator.cpp` explicitly instantiates
 * `panzer::MyEvaluator<panzer::Traits::Residual, panzer::Traits>`,
 * `<..::Jacobian, ..>`, `<..::Tangent, ..>`, and (if enabled)
 * `<..::Hessian, ..>`.
 */

// ONE template argument
/// \brief Explicitly instantiates `name<panzer::Traits::Residual>`. \param name the template class to instantiate.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_ONE_T(name) \
  template class name<panzer::Traits::Residual>;

/// \brief Explicitly instantiates `name<panzer::Traits::Tangent>`. \param name the template class to instantiate.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_ONE_T(name) \
  template class name<panzer::Traits::Tangent>;

/// \brief Explicitly instantiates `name<panzer::Traits::Jacobian>`. \param name the template class to instantiate.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_ONE_T(name) \
  template class name<panzer::Traits::Jacobian>;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  /// \brief Explicitly instantiates `name<panzer::Traits::Hessian>`. Only defined when built with Hessian support. \param name the template class to instantiate.
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_ONE_T(name) \
    template class name<panzer::Traits::Hessian>;
#else
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_ONE_T(name)
#endif

/**
 * \brief Explicitly instantiates `name<EvalT>` for every enabled evaluation type EvalT.
 * \param name the template class to instantiate; must take a single template parameter (the evaluation type).
 */
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_ONE_T(name)

// TWO template arguments
/// \brief Explicitly instantiates `name<panzer::Traits::Residual, panzer::Traits>`. \param name the template class to instantiate.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_TWO_T(name) \
  template class name<panzer::Traits::Residual, panzer::Traits>;

/// \brief Explicitly instantiates `name<panzer::Traits::Tangent, panzer::Traits>`. \param name the template class to instantiate.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_TWO_T(name) \
  template class name<panzer::Traits::Tangent, panzer::Traits>;

/// \brief Explicitly instantiates `name<panzer::Traits::Jacobian, panzer::Traits>`. \param name the template class to instantiate.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_TWO_T(name) \
  template class name<panzer::Traits::Jacobian, panzer::Traits>;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  /// \brief Explicitly instantiates `name<panzer::Traits::Hessian, panzer::Traits>`. Only defined when built with Hessian support. \param name the template class to instantiate.
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_TWO_T(name) \
    template class name<panzer::Traits::Hessian, panzer::Traits>;
#else
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_TWO_T(name)
#endif

/**
 * \brief Explicitly instantiates `name<EvalT, panzer::Traits>` for every enabled evaluation type EvalT.
 *
 * This is the macro used by most panzer evaluators, which are
 * templated on both the evaluation type and the traits class.
 * \param name the template class to instantiate; must take the evaluation type and `panzer::Traits` as its two template parameters.
 */
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_TWO_T(name)

// THREE (one user defined) template arguments
/// \brief Explicitly instantiates `name<panzer::Traits::Residual, panzer::Traits, ExtraT>`. \param name the template class to instantiate. \param ExtraT the additional, user-supplied template argument.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_THREE_T(name,ExtraT) \
  template class name<panzer::Traits::Residual, panzer::Traits,ExtraT>;

/// \brief Explicitly instantiates `name<panzer::Traits::Tangent, panzer::Traits, ExtraT>`. \param name the template class to instantiate. \param ExtraT the additional, user-supplied template argument.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_THREE_T(name,ExtraT) \
  template class name<panzer::Traits::Tangent, panzer::Traits,ExtraT>;

/// \brief Explicitly instantiates `name<panzer::Traits::Jacobian, panzer::Traits, ExtraT>`. \param name the template class to instantiate. \param ExtraT the additional, user-supplied template argument.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_THREE_T(name,ExtraT) \
  template class name<panzer::Traits::Jacobian, panzer::Traits,ExtraT>;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  /// \brief Explicitly instantiates `name<panzer::Traits::Hessian, panzer::Traits, ExtraT>`. Only defined when built with Hessian support. \param name the template class to instantiate. \param ExtraT the additional, user-supplied template argument.
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_THREE_T(name,ExtraT) \
    template class name<panzer::Traits::Hessian, panzer::Traits,ExtraT>;
#else
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_THREE_T(name,ExtraT)
#endif

/**
 * \brief Explicitly instantiates `name<EvalT, panzer::Traits, ExtraT>` for every enabled evaluation type EvalT.
 * \param name the template class to instantiate.
 * \param ExtraT the additional, user-supplied template argument (e.g. a specific field or algorithm tag).
 */
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_THREE_T(name,ExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_THREE_T(name,ExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_THREE_T(name,ExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_THREE_T(name,ExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_THREE_T(name,ExtraT)

// THREE (two user defined) template arguments
/// \brief Explicitly instantiates `name<panzer::Traits::Residual, FirstExtraT, SecondExtraT>` (no `panzer::Traits` argument). \param name the template class to instantiate. \param FirstExtraT first user-supplied template argument. \param SecondExtraT second user-supplied template argument.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Residual,FirstExtraT,SecondExtraT>;

/// \brief Explicitly instantiates `name<panzer::Traits::Tangent, FirstExtraT, SecondExtraT>` (no `panzer::Traits` argument). \param name the template class to instantiate. \param FirstExtraT first user-supplied template argument. \param SecondExtraT second user-supplied template argument.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Tangent,FirstExtraT,SecondExtraT>;

/// \brief Explicitly instantiates `name<panzer::Traits::Jacobian, FirstExtraT, SecondExtraT>` (no `panzer::Traits` argument). \param name the template class to instantiate. \param FirstExtraT first user-supplied template argument. \param SecondExtraT second user-supplied template argument.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Jacobian,FirstExtraT,SecondExtraT>;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  /// \brief Explicitly instantiates `name<panzer::Traits::Hessian, FirstExtraT, SecondExtraT>` (no `panzer::Traits` argument). Only defined when built with Hessian support. \param name the template class to instantiate. \param FirstExtraT first user-supplied template argument. \param SecondExtraT second user-supplied template argument.
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
    template class name<panzer::Traits::Hessian,FirstExtraT,SecondExtraT>;
#else
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT)
#endif

/**
 * \brief Explicitly instantiates `name<EvalT, FirstExtraT, SecondExtraT>` for every enabled evaluation type EvalT.
 *
 * Unlike PANZER_INSTANTIATE_TEMPLATE_CLASS_THREE_T, `panzer::Traits`
 * is not passed as a template argument here; both extra arguments are
 * user-supplied.
 * \param name the template class to instantiate.
 * \param FirstExtraT first user-supplied template argument.
 * \param SecondExtraT second user-supplied template argument.
 */
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT)

// FOUR (two user defined) template arguments
/// \brief Explicitly instantiates `name<panzer::Traits::Residual, panzer::Traits, FirstExtraT, SecondExtraT>`. \param name the template class to instantiate. \param FirstExtraT first user-supplied template argument. \param SecondExtraT second user-supplied template argument.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_FOUR_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Residual, panzer::Traits,FirstExtraT,SecondExtraT>;

/// \brief Explicitly instantiates `name<panzer::Traits::Tangent, panzer::Traits, FirstExtraT, SecondExtraT>`. \param name the template class to instantiate. \param FirstExtraT first user-supplied template argument. \param SecondExtraT second user-supplied template argument.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_FOUR_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Tangent, panzer::Traits,FirstExtraT,SecondExtraT>;

/// \brief Explicitly instantiates `name<panzer::Traits::Jacobian, panzer::Traits, FirstExtraT, SecondExtraT>`. \param name the template class to instantiate. \param FirstExtraT first user-supplied template argument. \param SecondExtraT second user-supplied template argument.
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_FOUR_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Jacobian, panzer::Traits,FirstExtraT,SecondExtraT>;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  /// \brief Explicitly instantiates `name<panzer::Traits::Hessian, panzer::Traits, FirstExtraT, SecondExtraT>`. Only defined when built with Hessian support. \param name the template class to instantiate. \param FirstExtraT first user-supplied template argument. \param SecondExtraT second user-supplied template argument.
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_FOUR_T(name,FirstExtraT,SecondExtraT) \
    template class name<panzer::Traits::Hessian, panzer::Traits,FirstExtraT,SecondExtraT>;
#else
  #define PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_FOUR_T(name,FirstExtraT,SecondExtraT)
#endif

/**
 * \brief Explicitly instantiates `name<EvalT, panzer::Traits, FirstExtraT, SecondExtraT>` for every enabled evaluation type EvalT.
 * \param name the template class to instantiate.
 * \param FirstExtraT first user-supplied template argument.
 * \param SecondExtraT second user-supplied template argument.
 */
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_FOUR_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_TANGENT_FOUR_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_FOUR_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_HESSIAN_FOUR_T(name,FirstExtraT,SecondExtraT)

#endif
