// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_MathExpr.hpp>

namespace Teuchos {

namespace MathExpr {

Language make_language() {
  Language out;
  Language::Productions& prods = out.productions;
  prods.resize(NPRODS);
  prods[PROD_PROGRAM]("program") >> "statements", "expr?";
  prods[PROD_NO_STATEMENTS]("statements");
  prods[PROD_NEXT_STATEMENT]("statements") >> "statements", "statement", ";", "S?";
  prods[PROD_ASSIGN]("statement") >> "name", "S?", "=", "S?", "expr";
  prods[PROD_NO_EXPR]("expr?");
  prods[PROD_YES_EXPR]("expr?") >> "expr";
  prods[PROD_EXPR]("expr") >> "ternary";
  prods[PROD_TERNARY_DECAY]("ternary") >> "add_sub";
  prods[PROD_OR_DECAY]("or") >> "and";
  prods[PROD_AND_DECAY]("and") >> "comp";
  prods[PROD_ADD_SUB_DECAY]("add_sub") >> "mul_div";
  prods[PROD_MUL_DIV_DECAY]("mul_div") >> "neg";
  prods[PROD_NEG_DECAY]("neg") >> "pow";
  prods[PROD_POW_DECAY]("pow") >> "scalar";
  prods[PROD_TERNARY]("ternary")
    >> "or", "?", "S?", "add_sub", ":", "S?", "add_sub";
  prods[PROD_OR]("or") >> "or", "||", "S?", "and";
  prods[PROD_AND]("and") >> "and", "&&", "S?", "comp";
  prods[PROD_GT]("comp") >> "add_sub", ">", "S?", "add_sub";
  prods[PROD_LT]("comp") >> "add_sub", "<", "S?", "add_sub";
  prods[PROD_GEQ]("comp") >> "add_sub", ">=", "S?", "add_sub";
  prods[PROD_LEQ]("comp") >> "add_sub", "<=", "S?", "add_sub";
  prods[PROD_EQ]("comp") >> "add_sub", "==", "S?", "add_sub";
  prods[PROD_BOOL_PARENS]("comp") >> "(", "S?", "or", ")", "S?";
  prods[PROD_ADD]("add_sub") >> "add_sub", "+", "S?", "mul_div";
  prods[PROD_SUB]("add_sub") >> "add_sub", "-", "S?", "mul_div";
  prods[PROD_MUL]("mul_div") >> "mul_div", "*", "S?", "pow";
  prods[PROD_DIV]("mul_div") >> "mul_div", "/", "S?", "pow";
  prods[PROD_POW]("pow") >> "scalar", "^", "S?", "pow";
  prods[PROD_CALL]("scalar")
    >> "name", "S?", "(", "S?", "args?", ")", "S?";
  prods[PROD_NO_ARGS]("args?");
  prods[PROD_SOME_ARGS]("args?") >> "args";
  prods[PROD_FIRST_ARG]("args") >> "ternary";
  prods[PROD_NEXT_ARG]("args") >> "args", ",", "S?", "ternary";
  prods[PROD_NEG]("neg") >> "-", "S?", "neg";
  prods[PROD_VAL_PARENS]("scalar") >> "(", "S?", "ternary", ")", "S?";
  prods[PROD_CONST]("scalar") >> "constant", "S?";
  prods[PROD_VAR]("scalar") >> "name", "S?";
  prods[PROD_NO_SPACES]("S?");
  prods[PROD_SPACES]("S?") >> "spaces";
  out.tokens.resize(NTOKS);
  out.tokens[TOK_SPACE]("spaces", "[ \t\n\r]+");
  out.tokens[TOK_NAME]("name", "[_a-zA-Z][_a-zA-Z0-9]*");
  out.tokens[TOK_ADD]("+", "\\+");
  out.tokens[TOK_SUB]("-", "\\-");
  out.tokens[TOK_MUL]("*", "\\*");
  out.tokens[TOK_DIV]("/", "\\/");
  out.tokens[TOK_POW]("^", "\\^");
  out.tokens[TOK_LPAREN]("(", "\\(");
  out.tokens[TOK_RPAREN](")", "\\)");
  out.tokens[TOK_COMMA](",", ",");
  out.tokens[TOK_CHECK]("?", "\\?");
  out.tokens[TOK_CHOOSE](":", ":");
  out.tokens[TOK_GT](">", ">");
  out.tokens[TOK_LT]("<", "<");
  out.tokens[TOK_GEQ](">=", ">=");
  out.tokens[TOK_LEQ]("<=", "<=");
  out.tokens[TOK_EQ]("==", "==");
  out.tokens[TOK_AND]("&&", "&&");
  out.tokens[TOK_OR]("||", "\\|\\|");
  out.tokens[TOK_CONST]("constant",
      "(0|([1-9][0-9]*))(\\.[0-9]*)?([eE]\\-?[1-9][0-9]*)?");
  out.tokens[TOK_SEMICOLON](";", ";");
  out.tokens[TOK_ASSIGN]("=", "=");
  return out;
}

LanguagePtr ask_language() {
  static LanguagePtr ptr;
  if (ptr.strong_count() == 0) {
    ptr.reset(new Language(make_language()));
  }
  return ptr;
}

ReaderTablesPtr ask_reader_tables() {
  static ReaderTablesPtr ptr;
  if (ptr.strong_count() == 0) {
    LanguagePtr lang = ask_language();
    ptr = make_reader_tables(*lang);
  }
  return ptr;
}

SymbolSetReader::SymbolSetReader():
  Reader(ask_reader_tables())
{
}

SymbolSetReader::~SymbolSetReader()
{
}

void SymbolSetReader::at_shift(any& result, int token, std::string& text) {
  if (token == TOK_NAME) result = text;
}

void SymbolSetReader::at_reduce(any& /* result */, int prod, std::vector<any>& rhs) {
  if (prod == PROD_VAR) {
    std::string& name = any_ref_cast<std::string>(rhs.at(0));
    variable_names.insert(name);
  } else if (prod == PROD_CALL) {
    std::string& name = any_ref_cast<std::string>(rhs.at(0));
    function_names.insert(name);
  }
}

std::set<std::string> get_variables_used(std::string const& expr) {
  SymbolSetReader reader;
  any result;
  reader.read_string(result, expr, "get_variables_used");
  return reader.variable_names;
}

std::set<std::string> get_symbols_used(std::string const& expr) {
  SymbolSetReader reader;
  any result;
  reader.read_string(result, expr, "get_symbols_used");
  auto set = std::move(reader.variable_names);
  set.insert(reader.function_names.begin(), reader.function_names.end());
  return set;
}

class CalcReader : public Reader {
 public:
  CalcReader():Reader(MathExpr::ask_reader_tables()) {
    unary_function_map["sqrt"] = &std::sqrt;
    unary_function_map["sin"] = &std::sin;
    unary_function_map["cos"] = &std::cos;
    unary_function_map["tan"] = &std::tan;
    unary_function_map["asin"] = &std::asin;
    unary_function_map["acos"] = &std::acos;
    unary_function_map["atan"] = &std::atan;
    unary_function_map["exp"] = &std::exp;
    unary_function_map["log"] = &std::log;
    unary_function_map["log10"] = &std::log10;
    binary_function_map["atan2"] = &std::atan2;
  }
  virtual ~CalcReader() = default;
 protected:
  struct CallArgs {
    double a0;
    double a1;
    int n;
  };
  virtual void at_shift(any& result_any, int token, std::string& text) {
    using std::swap;
    switch (token) {
      case MathExpr::TOK_NAME: {
        std::string& result = make_any_ref<std::string>(result_any);
        swap(result, text);
        return;
      }
      case MathExpr::TOK_CONST: {
        result_any = std::atof(text.c_str());
        return;
      }
    }
  }
  virtual void at_reduce(any& result, int prod, std::vector<any>& rhs) {
    using std::swap;
    switch (prod) {
      case MathExpr::PROD_PROGRAM: {
        TEUCHOS_TEST_FOR_EXCEPTION(!rhs.at(1).has_value(), ParserFail,
          "Calculator needs an expression to evaluate!");
        swap(result, rhs.at(1));
        break;
      }
      case MathExpr::PROD_NO_STATEMENTS:
      case MathExpr::PROD_NO_EXPR:
      case MathExpr::PROD_NEXT_STATEMENT: {
        break;
      }
      case MathExpr::PROD_ASSIGN: {
        std::string const& name = any_ref_cast<std::string>(rhs.at(0));
        double value = any_cast<double>(rhs.at(4));
        variable_map[name] = value;
        break;
      }
      case MathExpr::PROD_YES_EXPR:
      case MathExpr::PROD_EXPR:
      case MathExpr::PROD_TERNARY_DECAY:
      case MathExpr::PROD_OR_DECAY:
      case MathExpr::PROD_AND_DECAY:
      case MathExpr::PROD_ADD_SUB_DECAY:
      case MathExpr::PROD_MUL_DIV_DECAY:
      case MathExpr::PROD_POW_DECAY:
      case MathExpr::PROD_NEG_DECAY:
      case MathExpr::PROD_SOME_ARGS:
        swap(result, rhs.at(0));
        break;
      case MathExpr::PROD_TERNARY:
        result = any_cast<bool>(rhs.at(0)) ?
                 any_cast<double>(rhs.at(3)) :
                 any_cast<double>(rhs.at(6));
        break;
      case MathExpr::PROD_OR:
        result = any_cast<bool>(rhs.at(0)) || any_cast<bool>(rhs.at(3));
        break;
      case MathExpr::PROD_AND:
        result = any_cast<bool>(rhs.at(0)) && any_cast<bool>(rhs.at(3));
        break;
      case MathExpr::PROD_GT:
        result = any_cast<double>(rhs.at(0)) > any_cast<double>(rhs.at(3));
        break;
      case MathExpr::PROD_LT:
        result = any_cast<double>(rhs.at(0)) < any_cast<double>(rhs.at(3));
        break;
      case MathExpr::PROD_GEQ:
        result = any_cast<double>(rhs.at(0)) >= any_cast<double>(rhs.at(3));
        break;
      case MathExpr::PROD_LEQ:
        result = any_cast<double>(rhs.at(0)) <= any_cast<double>(rhs.at(3));
        break;
      case MathExpr::PROD_EQ:
        result = any_cast<double>(rhs.at(0)) == any_cast<double>(rhs.at(3));
        break;
      case MathExpr::PROD_BOOL_PARENS:
        result = any_cast<bool>(rhs.at(2));
        break;
      case MathExpr::PROD_ADD:
        result = any_cast<double>(rhs.at(0)) + any_cast<double>(rhs.at(3));
        break;
      case MathExpr::PROD_SUB:
        result = any_cast<double>(rhs.at(0)) - any_cast<double>(rhs.at(3));
        break;
      case MathExpr::PROD_MUL:
        result = any_cast<double>(rhs.at(0)) * any_cast<double>(rhs.at(3));
        break;
      case MathExpr::PROD_DIV:
        result = any_cast<double>(rhs.at(0)) / any_cast<double>(rhs.at(3));
        break;
      case MathExpr::PROD_POW:
        result = std::pow(any_cast<double>(rhs.at(0)), any_cast<double>(rhs.at(3)));
        break;
      case MathExpr::PROD_CALL: {
        std::string& name = any_ref_cast<std::string>(rhs.at(0));
        CallArgs& args = any_ref_cast<CallArgs>(rhs.at(4));
        TEUCHOS_TEST_FOR_EXCEPTION(args.n < 1 || args.n > 2, ParserFail,
            "Only unary and binary functions supported!\n");
        if (args.n == 1) {
          TEUCHOS_TEST_FOR_EXCEPTION(!unary_function_map.count(name), ParserFail,
              "Unknown unary function name \"" << name << "\"\n");
          Unary fptr = unary_function_map[name];
          result = (*fptr)(args.a0);
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(!binary_function_map.count(name), ParserFail,
              "Unknown binary function name \"" << name << "\"\n");
          Binary fptr = binary_function_map[name];
          result = (*fptr)(args.a0, args.a1);
        }
        break;
      }
      case MathExpr::PROD_NO_ARGS: {
        CallArgs& args = make_any_ref<CallArgs>(result);
        args.n = 0;
        break;
      }
      case MathExpr::PROD_FIRST_ARG: {
        CallArgs& args = make_any_ref<CallArgs>(result);
        args.a0 = any_cast<double>(rhs.at(0));
        args.n = 1;
        break;
      }
      case MathExpr::PROD_NEXT_ARG: {
        CallArgs& args = any_ref_cast<CallArgs>(rhs.at(0));
        args.a1 = any_cast<double>(rhs.at(3));
        args.n = 2;
        swap(result, rhs.at(0));
        break;
      }
      case MathExpr::PROD_NEG:
        result = - any_cast<double>(rhs.at(2));
        break;
      case MathExpr::PROD_VAL_PARENS:
        result = any_cast<double>(rhs.at(2));
        break;
      case MathExpr::PROD_CONST:
        result = any_cast<double>(rhs.at(0));
        break;
      case MathExpr::PROD_VAR:
        std::string const& name = any_ref_cast<std::string>(rhs.at(0));
        auto it = variable_map.find(name);
        TEUCHOS_TEST_FOR_EXCEPTION(it == variable_map.end(), ParserFail,
            "variable " << name << " not defined!");
        double value = it->second;
        result = value;
        break;
    }
  }
 private:
  typedef double (*Unary)(double);
  typedef double (*Binary)(double, double);
  std::map<std::string, Unary> unary_function_map;
  std::map<std::string, Binary> binary_function_map;
  std::map<std::string, double> variable_map;
};

Reader* new_calc_reader() {
  return new CalcReader();
}

}  // end namespace MathExpr

}  // end namespace Teuchos
