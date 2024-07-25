// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EXPR_EVAL_IMPL_HPP
#define PANZER_EXPR_EVAL_IMPL_HPP

#include <Panzer_ExprEval.hpp>

#include <algorithm>
#include <cmath>

namespace panzer
{
namespace Expr
{

struct ScalarTernary {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(bool cond, T const& left, T const& right) {
    return cond ? left : right;
  }
};

struct ScalarOr {
  static KOKKOS_FORCEINLINE_FUNCTION
  bool apply(bool left, bool right) {
    return left || right;
  }
};

struct ScalarAnd {
  static KOKKOS_FORCEINLINE_FUNCTION
  bool apply(bool left, bool right) {
    return left && right;
  }
};

struct ScalarGT {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  bool apply(T const& left, T const& right) {
    return left > right;
  }
};

struct ScalarLT {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  bool apply(T const& left, T const& right) {
    return left < right;
  }
};

struct ScalarGEQ {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  bool apply(T const& left, T const& right) {
    return left >= right;
  }
};

struct ScalarLEQ {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  bool apply(T const& left, T const& right) {
    return left <= right;
  }
};

struct ScalarEQ {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  bool apply(T const& left, T const& right) {
    return left == right;
  }
};

struct ScalarAdd {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(T const& left, T const& right) {
    return left + right;
  }
};

struct ScalarSub {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(T const& left, T const& right) {
    return left - right;
  }
};

struct ScalarMul {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(T const& left, T const& right) {
    return left * right;
  }
};

struct ScalarDiv {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(T const& left, T const& right) {
    return left / right;
  }
};

struct ScalarPow {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(T const& left, T const& right) {
    using std::pow;
    return pow(left, right);
  }
};

struct ScalarNeg {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(T const& right) {
    return -right;
  }
};

// TODO: replace this with .access() after next Kokkos release
template <typename Indexed, size_t IterationRank, size_t IndexedRank = Indexed::rank>
struct Indexer;

template <typename ViewType>
struct Indexer<ViewType, 1, 0> {
  template <typename Integral>
  static KOKKOS_FORCEINLINE_FUNCTION
  typename ViewType::reference_type index(ViewType const& x, Integral) { return x(); }
};

template <typename ViewType>
struct Indexer<ViewType, 1, 1> {
  template <typename Integral>
  static KOKKOS_FORCEINLINE_FUNCTION
  typename ViewType::reference_type index(ViewType const& x, Integral i) { return x(i); }
};

template <typename ViewType>
struct Indexer<ViewType, 2, 0> {
  template <typename Integral>
  static KOKKOS_FORCEINLINE_FUNCTION
  typename ViewType::reference_type index(ViewType const& x, Integral, Integral) { return x(); }
};

template <typename ViewType>
struct Indexer<ViewType, 2, 1> {
  template <typename Integral>
  static KOKKOS_FORCEINLINE_FUNCTION
  typename ViewType::reference_type index(ViewType const& x, Integral i, Integral) { return x(i); }
};

template <typename ViewType>
struct Indexer<ViewType, 2, 2> {
  template <typename Integral>
  static KOKKOS_FORCEINLINE_FUNCTION
  typename ViewType::reference_type index(ViewType const& x, Integral i, Integral j) { return x(i, j); }
};

//TODO: just use std::max once C++14 is the Trilinos standard (which makes std::max constexpr)
template <typename T, typename ... TS>
struct MaxRank;

template <typename T>
struct MaxRank<T> {
  static constexpr size_t value = T::rank;
};

template <typename T, typename ... TS>
struct MaxRank {
  static constexpr size_t left_value = T::rank;
  static constexpr size_t right_value = MaxRank<TS ...>::value;
  static constexpr size_t value = left_value > right_value ? left_value : right_value;
};

template <typename A, typename B>
struct ResultType {
  static constexpr size_t a_rank = A::rank;
  static constexpr size_t b_rank = B::rank;
  using biggest_type = typename std::conditional<(b_rank > a_rank), B, A>::type;
  using const_type = typename RebindViewType<biggest_type, typename biggest_type::const_value_type>::type;
  using type = typename RebindViewType<biggest_type, typename biggest_type::non_const_value_type>::type;
};

template <typename C, typename A, typename B>
struct TernaryResultType {
  static constexpr size_t a_rank = A::rank;
  static constexpr size_t b_rank = B::rank;
  static constexpr size_t c_rank = C::rank;
  using biggest_ab_type = typename std::conditional<(b_rank > a_rank), B, A>::type;
  using biggest_type = typename std::conditional<(c_rank > biggest_ab_type::rank), C, biggest_ab_type>::type;
  using type = typename RebindViewType<biggest_type, typename A::const_value_type>::type;
  using non_const_type = typename RebindViewType<biggest_type, typename A::non_const_value_type>::type;
};

template <typename Op, typename Result, typename Left, typename Right, size_t Rank = Result::rank>
struct BinaryFunctor;

template <typename Op, typename Result, typename Left, typename Right>
struct BinaryFunctor<Op, Result, Left, Right, 0> {
  using NonConstResult = typename RebindViewType<Result, typename Result::non_const_value_type>::type;
  using execution_space = typename Result::execution_space;
  NonConstResult result_;
  Left   left_;
  Right  right_;
  KOKKOS_INLINE_FUNCTION
  void operator()(typename execution_space::size_type) const {
    result_() = Op::apply(left_(), right_());
  }
  BinaryFunctor(std::string const& name, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) {
    left_ = Teuchos::any_cast<Left>(left);
    right_ = Teuchos::any_cast<Right>(right);
    result_ = NonConstResult(Kokkos::ViewAllocateWithoutInitializing(name));
    Kokkos::parallel_for(name, Kokkos::RangePolicy<execution_space>(0, 1), *this);
    result = Result(result_);
  }
};

template <typename Op, typename Result, typename Left, typename Right>
struct BinaryFunctor<Op, Result, Left, Right, 1> {
  using NonConstResult = typename RebindViewType<Result, typename Result::non_const_value_type>::type;
  using execution_space = typename Result::execution_space;
  NonConstResult result_;
  Left   left_;
  Right  right_;
  KOKKOS_INLINE_FUNCTION
  void operator()(typename execution_space::size_type i) const {
    result_(i) =
      Op::apply(
          Indexer<Left, 1>::index(left_, i),
          Indexer<Right, 1>::index(right_, i));
  }
  BinaryFunctor(std::string const& name, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) {
    left_ = Teuchos::any_cast<Left>(left);
    right_ = Teuchos::any_cast<Right>(right);
    auto extent_0 = std::max(left_.extent(0), right_.extent(0));
    result_ = NonConstResult(Kokkos::ViewAllocateWithoutInitializing(name), extent_0);
    Kokkos::parallel_for(name, Kokkos::RangePolicy<execution_space>(0, extent_0), *this);
    result = Result{result_};
  }
};

template <typename Op, typename Result, typename Left, typename Right>
struct BinaryFunctor<Op, Result, Left, Right, 2> {
  using NonConstResult = typename RebindViewType<Result, typename Result::non_const_value_type>::type;
  using execution_space = typename Result::execution_space;
  NonConstResult result_;
  Left   left_;
  Right  right_;
  KOKKOS_INLINE_FUNCTION
  void operator()(typename execution_space::size_type i, typename execution_space::size_type j) const {
    result_(i, j) =
      Op::apply(
          Indexer<Left, 2>::index(left_, i, j),
          Indexer<Right, 2>::index(right_, i, j));
  }
  BinaryFunctor(std::string const& name, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) {
    left_ = Teuchos::any_cast<Left>(left);
    right_ = Teuchos::any_cast<Right>(right);
    auto extent_0 = std::max(left_.extent(0), right_.extent(0));
    auto extent_1 = std::max(left_.extent(1), right_.extent(1));
    result_ = NonConstResult(Kokkos::ViewAllocateWithoutInitializing(name), extent_0, extent_1);
    using policy_type = Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<2>>;
    Kokkos::parallel_for(name, policy_type({0, 0}, {extent_0, extent_1}), *this);
    result = Result{result_};
  }
};

template <typename Cond, typename Left, typename Right, size_t Rank = MaxRank<Cond, Left, Right>::value>
struct TernaryFunctor;

template <typename Cond, typename Left, typename Right>
struct TernaryFunctor<Cond, Left, Right, 1> {
  using NonConstResult = typename TernaryResultType<Cond, Left, Right>::non_const_type;
  using Result = typename TernaryResultType<Cond, Left, Right>::type;
  using execution_space = typename Result::execution_space;
  NonConstResult result_;
  Cond   cond_;
  Left   left_;
  Right  right_;
  KOKKOS_INLINE_FUNCTION
  void operator()(typename execution_space::size_type i) const {
    result_(i) =
          Indexer<Cond, 1>::index(cond_, i) ?
          Indexer<Left, 1>::index(left_, i) :
          Indexer<Right, 1>::index(right_, i);
  }
  TernaryFunctor(std::string const& name, Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) {
    cond_ = Teuchos::any_cast<Cond>(cond);
    left_ = Teuchos::any_cast<Left>(left);
    right_ = Teuchos::any_cast<Right>(right);
    auto extent_0 =
      std::max(cond_.extent(0),
        std::max(left_.extent(0), right_.extent(0)));
    result_ = NonConstResult(Kokkos::ViewAllocateWithoutInitializing(name), extent_0);
    Kokkos::parallel_for(name, Kokkos::RangePolicy<execution_space>(0, extent_0), *this);
    result = Result{result_};
  }
};

template <typename Cond, typename Left, typename Right>
struct TernaryFunctor<Cond, Left, Right, 2> {
  using NonConstResult = typename TernaryResultType<Cond, Left, Right>::non_const_type;
  using Result = typename TernaryResultType<Cond, Left, Right>::type;
  using execution_space = typename Result::execution_space;
  NonConstResult result_;
  Cond   cond_;
  Left   left_;
  Right  right_;
  KOKKOS_INLINE_FUNCTION
  void operator()(typename execution_space::size_type i, typename execution_space::size_type j) const {
    result_(i, j) =
          Indexer<Cond, 2>::index(cond_, i, j) ?
          Indexer<Left, 2>::index(left_, i, j) :
          Indexer<Right, 2>::index(right_, i, j);
  }
  TernaryFunctor(std::string const& name, Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) {
    cond_ = Teuchos::any_cast<Cond>(cond);
    left_ = Teuchos::any_cast<Left>(left);
    right_ = Teuchos::any_cast<Right>(right);
    auto extent_0 =
      std::max(cond_.extent(0),
        std::max(left_.extent(0), right_.extent(0)));
    auto extent_1 =
      std::max(cond_.extent(1),
        std::max(left_.extent(1), right_.extent(1)));
    result_ = NonConstResult(Kokkos::ViewAllocateWithoutInitializing(name), extent_0, extent_1);
    using policy_type = Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<2>>;
    Kokkos::parallel_for(name, policy_type({0, 0}, {extent_0, extent_1}), *this);
    result = Result{result_};
  }
};

template <typename DT, typename ... VP>
Eval<DT, VP ...>::Eval()
  : EvalBase()
{
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::set(std::string const& name, bool value) {
  single_bool_view_type view{Kokkos::ViewAllocateWithoutInitializing{name}};
  auto host_view = Kokkos::create_mirror_view(view);
  host_view() = value;
  Kokkos::deep_copy(view, host_view);
  symbol_map[name] = const_single_bool_view_type{view};
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::set(std::string const& name, scalar_type const& value) {
  single_view_type view{Kokkos::ViewAllocateWithoutInitializing{name}};
  auto host_view = Kokkos::create_mirror_view(view);
  host_view() = value;
  Kokkos::deep_copy(view, host_view);
  symbol_map[name] = const_single_view_type{view};
  bool a, b;
  this->inspect_arg(symbol_map[name], a, b);
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::set(std::string const& name, const_view_type const& value) {
  symbol_map[name] = value;
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::make_constant(Teuchos::any& result, double const& value) {
  single_view_type view{Kokkos::ViewAllocateWithoutInitializing{"constant"}};
  auto host_view = Kokkos::create_mirror_view(view);
  host_view() = value;
  Kokkos::deep_copy(view, host_view);
  result = const_single_view_type{view};
  bool a, b;
  this->inspect_arg(result, a, b);
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::inspect_arg(Teuchos::any const& arg, bool& is_many, bool& is_bool) {
  if (arg.type() == typeid(const_single_bool_view_type)) {
    is_many = false;
    is_bool = true;
  } else if (arg.type() == typeid(const_single_view_type)) {
    is_many = false;
    is_bool = false;
  } else if (arg.type() == typeid(const_view_type)) {
    is_many = true;
    is_bool = false;
  } else if (arg.type() == typeid(const_bool_view_type)) {
    is_many = true;
    is_bool = true;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::ParserFail,
        "value is of illegal type " << arg.typeName() << ", view type is "
        << typeid(const_view_type).name());
  }
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::single_single_ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) {
  using ManyBool = const_bool_view_type;
  using Single = const_single_view_type;
  TernaryFunctor<ManyBool, Single, Single>("many ? single : single", result, cond, left, right);
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::single_many_ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) {
  using ManyBool = const_bool_view_type;
  using Single = const_single_view_type;
  using Many = const_view_type;
  TernaryFunctor<ManyBool, Single, Many>("many ? single : many", result, cond, left, right);
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::many_single_ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) {
  using ManyBool = const_bool_view_type;
  using Single = const_single_view_type;
  using Many = const_view_type;
  TernaryFunctor<ManyBool, Many, Single>("many ? many : single", result, cond, left, right);
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::many_many_ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) {
  using ManyBool = const_bool_view_type;
  using Many = const_view_type;
  TernaryFunctor<ManyBool, Many, Many>("many ? many : many", result, cond, left, right);
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::single_single_binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) {
  using SingleBool = const_single_bool_view_type;
  using Single = const_single_view_type;
  switch (code) {
    case BinaryOpCode::OR:  BinaryFunctor<ScalarOr , SingleBool, SingleBool, SingleBool>("single||single", result, left, right); break;
    case BinaryOpCode::AND: BinaryFunctor<ScalarAnd, SingleBool, SingleBool, SingleBool>("single&&single", result, left, right); break;
    case BinaryOpCode::LT:  BinaryFunctor<ScalarLT , SingleBool, Single, Single>("single< single", result, left, right); break;
    case BinaryOpCode::GT:  BinaryFunctor<ScalarGT , SingleBool, Single, Single>("single> single", result, left, right); break;
    case BinaryOpCode::GEQ: BinaryFunctor<ScalarGEQ, SingleBool, Single, Single>("single>=single", result, left, right); break;
    case BinaryOpCode::LEQ: BinaryFunctor<ScalarLEQ, SingleBool, Single, Single>("single<=single", result, left, right); break;
    case BinaryOpCode::EQ:  BinaryFunctor<ScalarEQ , SingleBool, Single, Single>("single==single", result, left, right); break;
    case BinaryOpCode::MUL: BinaryFunctor<ScalarMul, Single, Single, Single>("single* single", result, left, right); break;
    case BinaryOpCode::DIV: BinaryFunctor<ScalarDiv, Single, Single, Single>("single/ single", result, left, right); break;
    case BinaryOpCode::ADD: BinaryFunctor<ScalarAdd, Single, Single, Single>("single+ single", result, left, right); break;
    case BinaryOpCode::SUB: BinaryFunctor<ScalarSub, Single, Single, Single>("single- single", result, left, right); break;
    case BinaryOpCode::POW: BinaryFunctor<ScalarPow, Single, Single, Single>("single^ single", result, left, right); break;
  }
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::single_many_binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) {
  using Single = const_single_view_type;
  using SingleBool = const_single_bool_view_type;
  using Many = const_view_type;
  using ManyBool = const_bool_view_type;
  switch (code) {
    case BinaryOpCode::OR:  BinaryFunctor<ScalarOr , ManyBool, SingleBool, ManyBool>("single||many", result, left, right); break;
    case BinaryOpCode::AND: BinaryFunctor<ScalarAnd, ManyBool, SingleBool, ManyBool>("single&&many", result, left, right); break;
    case BinaryOpCode::LT:  BinaryFunctor<ScalarLT , ManyBool, Single, Many>("single< many", result, left, right); break;
    case BinaryOpCode::GT:  BinaryFunctor<ScalarGT , ManyBool, Single, Many>("single> many", result, left, right); break;
    case BinaryOpCode::GEQ: BinaryFunctor<ScalarGEQ, ManyBool, Single, Many>("single>=many", result, left, right); break;
    case BinaryOpCode::LEQ: BinaryFunctor<ScalarLEQ, ManyBool, Single, Many>("single<=many", result, left, right); break;
    case BinaryOpCode::EQ:  BinaryFunctor<ScalarEQ , ManyBool, Single, Many>("single==many", result, left, right); break;
    case BinaryOpCode::MUL: BinaryFunctor<ScalarMul, Many, Single, Many>("single* many", result, left, right); break;
    case BinaryOpCode::DIV: BinaryFunctor<ScalarDiv, Many, Single, Many>("single/ many", result, left, right); break;
    case BinaryOpCode::ADD: BinaryFunctor<ScalarAdd, Many, Single, Many>("single+ many", result, left, right); break;
    case BinaryOpCode::SUB: BinaryFunctor<ScalarSub, Many, Single, Many>("single- many", result, left, right); break;
    case BinaryOpCode::POW: BinaryFunctor<ScalarPow, Many, Single, Many>("single^ many", result, left, right); break;
  }
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::many_single_binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) {
  using Single = const_single_view_type;
  using SingleBool = const_single_bool_view_type;
  using Many = const_view_type;
  using ManyBool = const_bool_view_type;
  switch (code) {
    case BinaryOpCode::OR:  BinaryFunctor<ScalarOr , ManyBool, ManyBool, SingleBool>("many||single", result, left, right); break;
    case BinaryOpCode::AND: BinaryFunctor<ScalarAnd, ManyBool, ManyBool, SingleBool>("many&&single", result, left, right); break;
    case BinaryOpCode::LT:  BinaryFunctor<ScalarLT , ManyBool, Many, Single>("many< single", result, left, right); break;
    case BinaryOpCode::GT:  BinaryFunctor<ScalarGT , ManyBool, Many, Single>("many> single", result, left, right); break;
    case BinaryOpCode::GEQ: BinaryFunctor<ScalarGEQ, ManyBool, Many, Single>("many>=single", result, left, right); break;
    case BinaryOpCode::LEQ: BinaryFunctor<ScalarLEQ, ManyBool, Many, Single>("many<=single", result, left, right); break;
    case BinaryOpCode::EQ:  BinaryFunctor<ScalarEQ , ManyBool, Many, Single>("many==single", result, left, right); break;
    case BinaryOpCode::MUL: BinaryFunctor<ScalarMul, Many, Many, Single>("many* single", result, left, right); break;
    case BinaryOpCode::DIV: BinaryFunctor<ScalarDiv, Many, Many, Single>("many/ single", result, left, right); break;
    case BinaryOpCode::ADD: BinaryFunctor<ScalarAdd, Many, Many, Single>("many+ single", result, left, right); break;
    case BinaryOpCode::SUB: BinaryFunctor<ScalarSub, Many, Many, Single>("many- single", result, left, right); break;
    case BinaryOpCode::POW: BinaryFunctor<ScalarPow, Many, Many, Single>("many^ single", result, left, right); break;
  }
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::many_many_binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) {
  using Many = const_view_type;
  using ManyBool = const_bool_view_type;
  switch (code) {
    case BinaryOpCode::OR:  BinaryFunctor<ScalarOr , ManyBool, ManyBool, ManyBool>("many||many", result, left, right); break;
    case BinaryOpCode::AND: BinaryFunctor<ScalarAnd, ManyBool, ManyBool, ManyBool>("many&&many", result, left, right); break;
    case BinaryOpCode::LT:  BinaryFunctor<ScalarLT , ManyBool, Many, Many>("many< many", result, left, right); break;
    case BinaryOpCode::GT:  BinaryFunctor<ScalarGT , ManyBool, Many, Many>("many> many", result, left, right); break;
    case BinaryOpCode::GEQ: BinaryFunctor<ScalarGEQ, ManyBool, Many, Many>("many>=many", result, left, right); break;
    case BinaryOpCode::LEQ: BinaryFunctor<ScalarLEQ, ManyBool, Many, Many>("many<=many", result, left, right); break;
    case BinaryOpCode::EQ:  BinaryFunctor<ScalarEQ , ManyBool, Many, Many>("many==many", result, left, right); break;
    case BinaryOpCode::MUL: BinaryFunctor<ScalarMul, Many, Many, Many>("many* many", result, left, right); break;
    case BinaryOpCode::DIV: BinaryFunctor<ScalarDiv, Many, Many, Many>("many/ many", result, left, right); break;
    case BinaryOpCode::ADD: BinaryFunctor<ScalarAdd, Many, Many, Many>("many+ many", result, left, right); break;
    case BinaryOpCode::SUB: BinaryFunctor<ScalarSub, Many, Many, Many>("many- many", result, left, right); break;
    case BinaryOpCode::POW: BinaryFunctor<ScalarPow, Many, Many, Many>("many^ many", result, left, right); break;
  }
}

template <typename Op, typename Result, size_t Rank = Result::rank>
struct UnaryFunctor;

template <typename Op, typename Result>
struct UnaryFunctor<Op, Result, 0> {
  using NonConstResult = typename RebindViewType<Result, typename Result::non_const_value_type>::type;
  using execution_space = typename Result::execution_space;
  NonConstResult result_;
  Result right_;
  KOKKOS_INLINE_FUNCTION
  void operator()(typename execution_space::size_type i) const {
    result_() = Op::apply(right_());
  }
  UnaryFunctor(std::string const& name, Teuchos::any& result, Teuchos::any& right) {
    right_ = Teuchos::any_cast<Result>(right);
    result_ = NonConstResult(Kokkos::ViewAllocateWithoutInitializing(name));
    Kokkos::parallel_for(name, Kokkos::RangePolicy<execution_space>(0, 1), *this);
    result = Result(result_);
  }
};

template <typename Op, typename Result>
struct UnaryFunctor<Op, Result, 1> {
  using NonConstResult = typename RebindViewType<Result, typename Result::non_const_value_type>::type;
  using execution_space = typename Result::execution_space;
  NonConstResult result_;
  Result right_;
  KOKKOS_INLINE_FUNCTION
  void operator()(typename execution_space::size_type i) const {
    result_(i) = Op::apply(right_(i));
  }
  UnaryFunctor(std::string const& name, Teuchos::any& result, Teuchos::any& right) {
    right_ = Teuchos::any_cast<Result>(right);
    auto extent_0 = right_.extent(0);
    result_ = NonConstResult(Kokkos::ViewAllocateWithoutInitializing(name), extent_0);
    Kokkos::parallel_for(name, Kokkos::RangePolicy<execution_space>(0, extent_0), *this);
    result = Result(result_);
  }
};

template <typename Op, typename Result>
struct UnaryFunctor<Op, Result, 2> {
  using NonConstResult = typename RebindViewType<Result, typename Result::non_const_value_type>::type;
  using execution_space = typename Result::execution_space;
  NonConstResult result_;
  Result right_;
  KOKKOS_INLINE_FUNCTION
  void operator()(typename execution_space::size_type i, typename execution_space::size_type j) const {
    result_(i, j) = Op::apply(right_(i, j));
  }
  UnaryFunctor(std::string const& name, Teuchos::any& result, Teuchos::any& right) {
    right_ = Teuchos::any_cast<Result>(right);
    auto extent_0 = right_.extent(0);
    auto extent_1 = right_.extent(1);
    result_ = NonConstResult(Kokkos::ViewAllocateWithoutInitializing(name), extent_0, extent_1);
    using policy_type = Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<2>>;
    Kokkos::parallel_for(name, policy_type({0, 0}, {extent_0, extent_1}), *this);
    result = Result(result_);
  }
};

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::many_neg_op(Teuchos::any& result, Teuchos::any& right) {
  UnaryFunctor<ScalarNeg, const_view_type>("-many", result, right);
}

template <typename DT, typename ... VP>
void Eval<DT, VP ...>::single_neg_op(Teuchos::any& result, Teuchos::any& right) {
  UnaryFunctor<ScalarNeg, const_single_view_type>("-single", result, right);
}

struct ScalarAbs {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(T const& right) {
    using std::abs;
    return abs(right);
  }
};

struct ScalarExp {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(T const& right) {
    using std::exp;
    return exp(right);
  }
};

struct ScalarLog {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(T const& right) {
    using std::log;
    return log(right);
  }
};

struct ScalarSqrt {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(T const& right) {
    using std::sqrt;
    return sqrt(right);
  }
};

struct ScalarSin {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(T const& right) {
    using std::sin;
    return sin(right);
  }
};

struct ScalarCos {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(T const& right) {
    using std::cos;
    return cos(right);
  }
};

struct ScalarTan {
  template <typename T>
  static KOKKOS_FORCEINLINE_FUNCTION
  T apply(T const& right) {
    using std::tan;
    return tan(right);
  }
};

template <typename Op, typename EvalType>
struct UnaryFunction {
  void operator()(std::string const& name, Teuchos::any& result, std::vector<Teuchos::any>& rhs) const {
    auto& right = rhs.at(0);
    using single_type = typename EvalType::const_single_view_type;
    using many_type = typename EvalType::const_view_type;
    if (right.type() == typeid(single_type)) {
      UnaryFunctor<Op, single_type>(name, result, right);
    } else if (right.type() == typeid(many_type)) {
      UnaryFunctor<Op, many_type>(name, result, right);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::ParserFail,
          "Unexpected type " << right.typeName() << " passed to UnaryFunction \"" << name << "\"\n");
    }
  }
};

template <typename DT, typename ... VP>
void set_cmath_functions(Eval<DT, VP ...>& eval) {
  using eval_type = Eval<DT, VP ...>;
  EvalBase& eval_base = eval;
  eval_base.set("abs", EvalBase::Function(UnaryFunction<ScalarAbs, eval_type>{}));
  eval_base.set("exp", EvalBase::Function(UnaryFunction<ScalarExp, eval_type>{}));
  eval_base.set("log", EvalBase::Function(UnaryFunction<ScalarLog, eval_type>{}));
  eval_base.set("sqrt", EvalBase::Function(UnaryFunction<ScalarSqrt, eval_type>{}));
  eval_base.set("sin", EvalBase::Function(UnaryFunction<ScalarSin, eval_type>{}));
  eval_base.set("cos", EvalBase::Function(UnaryFunction<ScalarCos, eval_type>{}));
  eval_base.set("tan", EvalBase::Function(UnaryFunction<ScalarTan, eval_type>{}));
}

}} // end namespace panzer::Expr

#endif // PANZER_EXPR_EVAL_IMPL_HPP
