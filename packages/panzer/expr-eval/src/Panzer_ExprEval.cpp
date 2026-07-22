// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Panzer_ExprEval_impl.hpp>

#include <cstdlib>

#include <Teuchos_MathExpr.hpp>

namespace panzer
{
namespace Expr
{

EvalBase::EvalBase()
  : Teuchos::Reader(Teuchos::MathExpr::ask_reader_tables()) {
}

void EvalBase::set(std::string const& name, Function const& value) {
  symbol_map[name] = value;
}

void EvalBase::at_shift(Teuchos::any& result_any, int token, std::string& text) {
  using std::swap;
  switch (token) {
    case Teuchos::MathExpr::TOK_NAME: {
      std::string& result = Teuchos::make_any_ref<std::string>(result_any);
      swap(result, text);
      return;
    }
    case Teuchos::MathExpr::TOK_CONST: {
      this->make_constant(result_any, std::atof(text.c_str()));
      return;
    }
  }
}

void EvalBase::at_reduce(Teuchos::any& result, int prod, std::vector<Teuchos::any>& rhs) {
  using std::swap;
  switch (prod) {
    case Teuchos::MathExpr::PROD_PROGRAM: {
      swap(result, rhs.at(1));
      break;
    }
    case Teuchos::MathExpr::PROD_NO_STATEMENTS:
    case Teuchos::MathExpr::PROD_NO_EXPR:
    case Teuchos::MathExpr::PROD_NEXT_STATEMENT: {
      break;
    }
    case Teuchos::MathExpr::PROD_ASSIGN: {
      std::string const& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
      swap(symbol_map[name], rhs.at(4));
      break;
    }
    case Teuchos::MathExpr::PROD_YES_EXPR:
    case Teuchos::MathExpr::PROD_EXPR:
    case Teuchos::MathExpr::PROD_TERNARY_DECAY:
    case Teuchos::MathExpr::PROD_OR_DECAY:
    case Teuchos::MathExpr::PROD_AND_DECAY:
    case Teuchos::MathExpr::PROD_ADD_SUB_DECAY:
    case Teuchos::MathExpr::PROD_MUL_DIV_DECAY:
    case Teuchos::MathExpr::PROD_POW_DECAY:
    case Teuchos::MathExpr::PROD_NEG_DECAY:
    case Teuchos::MathExpr::PROD_SOME_ARGS:
      swap(result, rhs.at(0));
      break;
    case Teuchos::MathExpr::PROD_TERNARY:
      this->ternary_op(result, rhs.at(0), rhs.at(3), rhs.at(6));
      break;
    case Teuchos::MathExpr::PROD_OR:
      this->binary_op(BinaryOpCode::OR, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_AND:
      this->binary_op(BinaryOpCode::AND, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_GT:
      this->binary_op(BinaryOpCode::GT, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_LT:
      this->binary_op(BinaryOpCode::LT, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_GEQ:
      this->binary_op(BinaryOpCode::GEQ, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_LEQ:
      this->binary_op(BinaryOpCode::LEQ, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_EQ:
      this->binary_op(BinaryOpCode::EQ, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_BOOL_PARENS:
      swap(result, rhs.at(2));
      break;
    case Teuchos::MathExpr::PROD_ADD:
      this->binary_op(BinaryOpCode::ADD, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_SUB:
      this->binary_op(BinaryOpCode::SUB, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_MUL:
      this->binary_op(BinaryOpCode::MUL, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_DIV:
      this->binary_op(BinaryOpCode::DIV, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_POW:
      this->binary_op(BinaryOpCode::POW, result, rhs.at(0), rhs.at(3));
      break;
    case Teuchos::MathExpr::PROD_CALL: {
      std::string const& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
      auto it = symbol_map.find(name);
      TEUCHOS_TEST_FOR_EXCEPTION(it == symbol_map.end(), Teuchos::ParserFail,
          "symbol \"" << name << "\" being called doesn't exist!");
      Function& func = Teuchos::any_ref_cast<Function>(it->second);
      std::vector<Teuchos::any>& args = Teuchos::any_ref_cast<std::vector<Teuchos::any>>(rhs.at(4));
      func(name, result, args);
      break;
    }
    case Teuchos::MathExpr::PROD_NO_ARGS: {
      result = std::vector<Teuchos::any>{};
      break;
    }
    case Teuchos::MathExpr::PROD_FIRST_ARG: {
      std::vector<Teuchos::any>& args = Teuchos::make_any_ref<std::vector<Teuchos::any>>(result);
      args.push_back(Teuchos::any{});
      swap(args.back(), rhs.at(0));
      break;
    }
    case Teuchos::MathExpr::PROD_NEXT_ARG: {
      swap(result, rhs.at(0));
      std::vector<Teuchos::any>& args = Teuchos::any_ref_cast<std::vector<Teuchos::any>>(result);
      args.push_back(Teuchos::any{});
      swap(args.back(), rhs.at(3));
      break;
    }
    case Teuchos::MathExpr::PROD_NEG:
      this->neg_op(result, rhs.at(2));
      break;
    case Teuchos::MathExpr::PROD_VAL_PARENS:
      swap(result, rhs.at(2));
      break;
    case Teuchos::MathExpr::PROD_CONST:
      swap(result, rhs.at(0));
      break;
    case Teuchos::MathExpr::PROD_VAR:
      std::string const& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
      auto it = symbol_map.find(name);
      TEUCHOS_TEST_FOR_EXCEPTION(it == symbol_map.end(), Teuchos::ParserFail,
          "symbol " << name << " being referenced doesn't exist!");
      result = it->second;
      break;
  }
}

void EvalBase::ternary_op(Teuchos::any& result, Teuchos::any& cond, Teuchos::any& left, Teuchos::any& right) {
  bool cond_is_many;
  bool cond_is_bool;
  this->inspect_arg(cond, cond_is_many, cond_is_bool);
  TEUCHOS_TEST_FOR_EXCEPTION(!cond_is_bool, Teuchos::ParserFail,
      "Ternary condition is not of boolean type!");
  bool is_many[2];
  bool is_bool[2];
  this->inspect_arg(left, is_many[0], is_bool[0]);
  this->inspect_arg(right, is_many[1], is_bool[1]);
  TEUCHOS_TEST_FOR_EXCEPTION(is_bool[0], Teuchos::ParserFail,
      "Boolean values in ternary operator not yet supported");
  if (!cond_is_many) {
    auto cond_value = Teuchos::any_cast<Kokkos::View<bool const>>(cond);
    auto host_cond_value = Kokkos::create_mirror_view(cond_value);
    Kokkos::deep_copy(host_cond_value, cond_value);
    if (host_cond_value()) {
      swap(result, left);
    } else {
      swap(result, right);
    }
  } else {
    if (!is_many[0] && !is_many[1]) {
      this->single_single_ternary_op(result, cond, left, right);
    } else if (!is_many[0] && is_many[1]) {
      this->single_many_ternary_op(result, cond, left, right);
    } else if (is_many[0] && !is_many[1]) {
      this->many_single_ternary_op(result, cond, left, right);
    } else if (is_many[0] && is_many[1]) {
      this->many_many_ternary_op(result, cond, left, right);
    }
  }
}

static const char* get_op_syntax(BinaryOpCode code) {
  switch (code) {
    case BinaryOpCode::OR: return "||";
    case BinaryOpCode::AND: return "&&";
    case BinaryOpCode::GT: return ">";
    case BinaryOpCode::LT: return "<";
    case BinaryOpCode::GEQ: return ">=";
    case BinaryOpCode::LEQ: return "<=";
    case BinaryOpCode::EQ: return "==";
    case BinaryOpCode::ADD: return "+";
    case BinaryOpCode::SUB: return "-";
    case BinaryOpCode::MUL: return "*";
    case BinaryOpCode::DIV: return "/";
    case BinaryOpCode::POW: return "^";
  }
  return "";
}

void EvalBase::binary_op(BinaryOpCode code, Teuchos::any& result, Teuchos::any& left, Teuchos::any& right) {
  bool is_many[2];
  bool is_bool[2];
  this->inspect_arg(left, is_many[0], is_bool[0]);
  this->inspect_arg(right, is_many[1], is_bool[1]);
  bool expect_booleans = (code == BinaryOpCode::AND || code == BinaryOpCode::OR);
  TEUCHOS_TEST_FOR_EXCEPTION(is_bool[0] != expect_booleans, Teuchos::ParserFail,
      "Left argument to '" << get_op_syntax(code) << "' is " << (is_bool[0] ? "" : "not") << " boolean!");
  TEUCHOS_TEST_FOR_EXCEPTION(is_bool[1] != expect_booleans, Teuchos::ParserFail,
      "Right argument to '" << get_op_syntax(code) << "' is " << (is_bool[0] ? "" : "not") << " boolean!");
  if (!is_many[0] && !is_many[1]) {
    this->single_single_binary_op(code, result, left, right);
  } else if (!is_many[0] && is_many[1]) {
    this->single_many_binary_op(code, result, left, right);
  } else if (is_many[0] && !is_many[1]) {
    this->many_single_binary_op(code, result, left, right);
  } else if (is_many[0] && is_many[1]) {
    this->many_many_binary_op(code, result, left, right);
  }
}

void EvalBase::neg_op(Teuchos::any& result, Teuchos::any& right) {
  bool is_many;
  bool is_bool;
  this->inspect_arg(right, is_many, is_bool);
  TEUCHOS_TEST_FOR_EXCEPTION(is_bool, Teuchos::ParserFail,
      "Can't negate a boolean");
  if (is_many) {
    this->many_neg_op(result, right);
  } else {
    this->single_neg_op(result, right);
  }
}

}} // end namespace panzer::Expr
