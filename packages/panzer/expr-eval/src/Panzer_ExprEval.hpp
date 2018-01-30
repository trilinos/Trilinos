// @HEADER                                                                                                                                    
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_EXPR_EVAL_HPP
#define PANZER_EXPR_EVAL_HPP

#include <functional>
#include <map>
#include <type_traits>

#include <Teuchos_Reader.hpp>

#include <Kokkos_Core.hpp>

namespace panzer
{
namespace Expr
{

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

class EvalBase : public Teuchos::Reader {
 public:
  EvalBase();
  using Function = std::function<void(std::string const& name, Teuchos::any&, std::vector<Teuchos::any>& rhs)>;
  void set(std::string const& name, Function const& value);
  template <typename T>
  T const& get(std::string const& name) const {
    auto it = symbol_map.find(name);
    TEUCHOS_TEST_FOR_EXCEPTION(it == symbol_map.end(), std::logic_error,
        "EvalBase::get: \"" << name << "\" not found");
    return Teuchos::any_ref_cast<T>(it->second);
  }
 protected:
  std::map<std::string, Teuchos::any> symbol_map;

  void at_shift(Teuchos::any& result, int token, std::string& text) override;
  void at_reduce(Teuchos::any& result, int prod, std::vector<Teuchos::any>& rhs) override;
  void ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right);
  void binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right);
  void neg_op(Teuchos::any& result, Teuchos::any& right);
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
};

template <typename DataType, typename NewScalarType>
struct RebindDataType {
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

template <typename ViewType, typename NewScalarType>
struct RebindViewType;

template <typename DT, typename NewScalarType, typename ... VP>
struct RebindViewType<Kokkos::View<DT, VP ...>, NewScalarType> {
  using type = Kokkos::View<typename RebindDataType<DT, NewScalarType>::type, VP ...>;
};

template <typename DT, typename ... VP>
class Eval : public EvalBase {
 public:
  using original_view_type = Kokkos::View<DT, VP ...>;
  using view_data_type = typename original_view_type::data_type;
  using scalar_type = typename original_view_type::non_const_value_type;
  using view_type = Kokkos::View<typename RebindDataType<view_data_type, scalar_type>::type, VP ...>;
  using const_view_type = Kokkos::View<typename RebindDataType<view_data_type, scalar_type const>::type, VP ...>;
  using bool_view_type = Kokkos::View<typename RebindDataType<view_data_type, bool>::type, VP ...>;
  using const_bool_view_type = Kokkos::View<typename RebindDataType<view_data_type, bool const>::type, VP ...>;
  using single_view_type = Kokkos::View<scalar_type, VP ...>;
  using const_single_view_type = Kokkos::View<scalar_type const, VP ...>;
  using single_bool_view_type = Kokkos::View<bool, VP ...>;
  using const_single_bool_view_type = Kokkos::View<bool const, VP ...>;

  Eval();

  void set(std::string const& name, bool value);
  void set(std::string const& name, scalar_type const& value);
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

template <typename DT, typename ... VP>
void set_cmath_functions(Eval<DT, VP ...>& eval);

}} // end namespace panzer::Expr

#endif // PANZER_EXPR_EVAL_HPP
