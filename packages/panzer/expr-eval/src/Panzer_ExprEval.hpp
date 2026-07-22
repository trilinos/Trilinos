// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EXPR_EVAL_HPP
#define PANZER_EXPR_EVAL_HPP

/** \file Panzer_ExprEval.hpp
 *  \brief Declares the panzer::Expr::Eval templated class
 */

#include <functional>
#include <map>
#include <type_traits>

#include <Teuchos_Reader.hpp>

#include <Kokkos_Core.hpp>

namespace panzer
{

/**
 * \brief Contains all symbols which support panzer::Expr::Eval
 */
namespace Expr
{

/**
 * \brief Denotes the native binary operators in the Teuchos::MathExpr
 *        language.
 */
enum class BinaryOpCode {
  OR,
  AND,
  GT,
  LT,
  GEQ,
  LEQ,
  EQ,
  ADD,
  SUB,
  MUL,
  DIV,
  POW,
};

/**
 * \brief Base class for panzer::Expr::Eval, does everything that is independent
 *        of the Kokkos::View template parameter.
 */
class EvalBase : public Teuchos::Reader {
 public:
  /**
   * \brief Constructor.
   */
  EvalBase();
  /**
   * \brief The type of user-defined functions which are callable in the math language
   *
   * The first argument will be the name that was used to call the function,
   * the second argument is the return value, and the third argument is a vector
   * containing the actual arguments given in the math language.
   * Free functions, functors, and lambdas can all be used to construct a std::function.
   */
  using Function = std::function<void(std::string const& name, Teuchos::any&, std::vector<Teuchos::any>& rhs)>;
  /**
   * \brief Registers an EvalBase::Function, binding it to a name and making it callable.
   */
  void set(std::string const& name, Function const& value);
  /**
   * \brief Get the value of a variable in the symbol map.
   */
  template <typename T>
  T const& get(std::string const& name) const {
    auto it = symbol_map.find(name);
    TEUCHOS_TEST_FOR_EXCEPTION(it == symbol_map.end(), std::logic_error,
        "EvalBase::get: \"" << name << "\" not found");
    return Teuchos::any_ref_cast<T>(it->second);
  }
 protected:
  /**
   * \brief Stores all current symbols including variables and functions
   */
  std::map<std::string, Teuchos::any> symbol_map;

  /// Called at every parsed token in the math language
  void at_shift(Teuchos::any& result, int token, std::string& text) override;
  /// Called at every reduced production in the math language
  void at_reduce(Teuchos::any& result, int prod, std::vector<Teuchos::any>& rhs) override;
  /// Executes the ternary operator, e.g. (a > b) ? a : b
  void ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right);
  /**
   * \brief Executes a binary operator
   * \param code Which binary operator to execute
   * \param result Holds the resulting value
   * \param left Holds the left operand
   * \param right Holds the right operand
   *
   * This function figures out which operands are singular versus plural and
   * then dispatches to the appropriate virtual functions (e.g. many_single_binary_op),
   * which are implemented in Expr::Eval.
   */
  void binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right);
  /// Executes the only native unary operator in the math language, numeric negation via a minus sign
  void neg_op(Teuchos::any& result, Teuchos::any& right);
 protected:

  /** \name Dispatch to template-specialized Expr::Eval code. */
  /** @{ */
  virtual void make_constant(Teuchos::any& result, double const& value) = 0;
  virtual void inspect_arg(Teuchos::any const& arg, bool& is_many, bool& is_bool) = 0;
  virtual void single_single_ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) = 0;
  virtual void single_many_ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) = 0;
  virtual void many_single_ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) = 0;
  virtual void many_many_ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) = 0;
  virtual void single_single_binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) = 0;
  virtual void single_many_binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) = 0;
  virtual void many_single_binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) = 0;
  virtual void many_many_binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) = 0;
  virtual void many_neg_op(Teuchos::any& result, Teuchos::any& right) = 0;
  virtual void single_neg_op(Teuchos::any& result, Teuchos::any& right) = 0;
  /** @} */
};

/// Rebinds a Kokkos::View data type to use a new scalar type
template <typename DataType, typename NewScalarType>
struct RebindDataType {
  /// The new data type, suitable as the first template argument to Kokkos::View
  using type = NewScalarType;
};

template <typename NestedDataType, typename NewScalarType>
struct RebindDataType<NestedDataType*, NewScalarType> {
  using type = typename RebindDataType<NestedDataType, NewScalarType>::type *;
};

template <typename NestedDataType, typename NewScalarType>
struct RebindDataType<NestedDataType[], NewScalarType> {
  using type = typename RebindDataType<NestedDataType, NewScalarType>::type [];
};

template <typename NestedDataType, typename NewScalarType, size_t N>
struct RebindDataType<NestedDataType[N], NewScalarType> {
  using type = typename RebindDataType<NestedDataType, NewScalarType>::type [N];
};

/// Builds on RebindDataType, but acts directly on a Kokkos::View type
template <typename ViewType, typename NewScalarType>
struct RebindViewType;

template <typename DT, typename NewScalarType, typename ... VP>
struct RebindViewType<Kokkos::View<DT, VP ...>, NewScalarType> {
  /// The new Kokkos::View type, whose scalar type is now NewScalarType
  using type = Kokkos::View<typename RebindDataType<DT, NewScalarType>::type, VP ...>;
};

/**
 * \brief Interprets mathematical expressions in a string and evaluates them
 *        using Kokkos::View objects as values and Kokkos::parallel_for
 *        for the operators.
 * This class is mean to support Kokkos-parallel execution of user-provided
 * mathematical expressions.
 * Example uses include evaluating analytic boundary and initial conditions for
 * a PDE problem.
 * The API of this class (namely Eval::read_string or another read_* function inherited
 * from Teuchos::Reader) is only meant to be called once for a given set of
 * evaluation points.
 * Variables should be predefined with point-dependent values all stored in a
 * single Kokkos::View.
 * For example, X coordinates for all points could be stored in a Kokkos::View<double*>,
 * whose extent is the number of points.
 * Values which are the same for all points are still supported, for example the current time.
 * Then if the expression is for example "x^t", this class launches a single parallel_for
 * to raise the coordinate of all points at once to the "t" power.
 *
 * \note This class currently supports rank-1 and rank-2 Kokkos::Views.
 *       It should support all scalar types which support the native binary operators
 *       ('+', '-', '*', '/'), and also pow(x, y).
 *
 * \note Expression results are allocated when an expression is evaluated, and deallocated
 *       after they have been used in all higher-level expressions.
 *       For example, "(a * b) + c" does the following:
 *        - allocate space for "(a * b)"
 *        - compute and store "(a * b)"
 *        - allocate space for "(a * b) + c"
 *        - compute and store "(a * b) + c"
 *        - deallocate space for "(a * b)"
 */
template <typename DT, typename ... VP>
class Eval : public EvalBase {
 public:
  /// The corresponding Kokkos::View type, using the same template arguments are were given to Eval
  using original_view_type = Kokkos::View<DT, VP ...>;
  /// The data type, including dimension information
  using view_data_type = DT;
  /**
   * \brief The scalar type
   *
   * This type should be constructible from a value of type "double" and should be usable
   * with the four basic binary operator ('+', '-', '*', '/') and also the pow() function.
   */
  using scalar_type = typename original_view_type::non_const_value_type;
  /**
   *  \brief One scalar for each evaluation point, read-only
   *
   *  This is the original view type where the scalar type is constant.
   *  This means it is a read-only view.
   *  and will most likely be the type that is returned from read_string().
   */
  using const_view_type = Kokkos::View<typename RebindDataType<view_data_type, scalar_type const>::type, VP ...>;
  /**
   *  \brief One boolean for each evaluation point, read-only
   */
  using const_bool_view_type = Kokkos::View<typename RebindDataType<view_data_type, bool const>::type, VP ...>;
  /**
   *  \brief One scalar (same for all evaluation points)
   */
  using single_view_type = Kokkos::View<scalar_type, VP ...>;
  /**
   *  \brief One scalar (same for all evaluation points), read-only
   */
  using const_single_view_type = Kokkos::View<scalar_type const, VP ...>;
  /**
   *  \brief One boolean (same for all evaluation points)
   */
  using single_bool_view_type = Kokkos::View<bool, VP ...>;
  /**
   *  \brief One boolean (same for all evaluation points), read-only
   */
  using const_single_bool_view_type = Kokkos::View<bool const, VP ...>;

  Eval();

  /**
   *  \brief Assign a boolean value to a variable symbol
   */
  void set(std::string const& name, bool value);
  /**
   *  \brief Assign a scalar value to a variable symbol
   */
  void set(std::string const& name, scalar_type const& value);
  /**
   *  \brief Assign scalar values (one for each evaluation point) to a variable symbol
   */
  void set(std::string const& name, const_view_type const& value);
 protected:
  void make_constant(Teuchos::any& result, double const& value) override;
  void inspect_arg(Teuchos::any const& arg, bool& is_many, bool& is_bool) override;
  void single_single_ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) override;
  void single_many_ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) override;
  void many_single_ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) override;
  void many_many_ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) override;
  void single_single_binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) override;
  void single_many_binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) override;
  void many_single_binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) override;
  void many_many_binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) override;
  void many_neg_op(Teuchos::any& result, Teuchos::any& right) override;
  void single_neg_op(Teuchos::any& result, Teuchos::any& right) override;
};

/**
 *  \brief Add support for functions such as sqrt(), sin(), and cos()
 *
 *  After this call, (eval) is able to interpret strings which include
 *  calls to these functions:
 *    - abs(x)
 *    - exp(x)
 *    - log(x)
 *    - sqrt(x)
 *    - sin(x)
 *    - cos(x)
 *    - tan(x)
 *
 * \note This function is kept separate in case users want to use Expr::Eval
 * with a scalar type that only supports basic operators but doesn't support
 * these higher-level math functions (an integer for example).
 */
template <typename DT, typename ... VP>
void set_cmath_functions(Eval<DT, VP ...>& eval);

}} // end namespace panzer::Expr

#endif // PANZER_EXPR_EVAL_HPP
