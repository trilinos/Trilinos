#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <Teuchos_Language.hpp>
#include <Teuchos_Reader.hpp>
#include <Teuchos_MathExpr.hpp>

namespace {

using Teuchos::any;
using Teuchos::any_cast;

struct CallArgs {
  double a0;
  double a1;
  int n;
};

class Reader : public Teuchos::Reader {
 public:
  Reader():Teuchos::Reader(Teuchos::MathExpr::ask_reader_tables()) {
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
  virtual ~Reader() {}
 protected:
  virtual void at_shift(Teuchos::any& result_any, int token, std::string& text) {
    using std::swap;
    switch (token) {
      case Teuchos::MathExpr::TOK_NAME: {
        std::string& result = Teuchos::make_any_ref<std::string>(result_any);
        swap(result, text);
        return;
      }
      case Teuchos::MathExpr::TOK_CONST: {
        result_any = std::atof(text.c_str());
        return;
      }
    }
  }
  virtual void at_reduce(Teuchos::any& result, int prod, std::vector<Teuchos::any>& rhs) {
    using std::swap;
    switch (prod) {
      case Teuchos::MathExpr::PROD_PROGRAM: {
        TEUCHOS_TEST_FOR_EXCEPTION(rhs.at(1).empty(), Teuchos::ParserFail,
          "Calculator needs an expression to evaluate!");
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
        double value = any_cast<double>(rhs.at(4));
        variable_map[name] = value;
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
        result = any_cast<bool>(rhs.at(0)) ?
                 any_cast<double>(rhs.at(3)) :
                 any_cast<double>(rhs.at(6));
        break;
      case Teuchos::MathExpr::PROD_OR:
        result = any_cast<bool>(rhs.at(0)) || any_cast<bool>(rhs.at(3));
        break;
      case Teuchos::MathExpr::PROD_AND:
        result = any_cast<bool>(rhs.at(0)) && any_cast<bool>(rhs.at(3));
        break;
      case Teuchos::MathExpr::PROD_GT:
        result = any_cast<double>(rhs.at(0)) > any_cast<double>(rhs.at(3));
        break;
      case Teuchos::MathExpr::PROD_LT:
        result = any_cast<double>(rhs.at(0)) < any_cast<double>(rhs.at(3));
        break;
      case Teuchos::MathExpr::PROD_GEQ:
        result = any_cast<double>(rhs.at(0)) >= any_cast<double>(rhs.at(3));
        break;
      case Teuchos::MathExpr::PROD_LEQ:
        result = any_cast<double>(rhs.at(0)) <= any_cast<double>(rhs.at(3));
        break;
      case Teuchos::MathExpr::PROD_EQ:
        result = any_cast<double>(rhs.at(0)) == any_cast<double>(rhs.at(3));
        break;
      case Teuchos::MathExpr::PROD_BOOL_PARENS:
        result = any_cast<bool>(rhs.at(2));
        break;
      case Teuchos::MathExpr::PROD_ADD:
        result = any_cast<double>(rhs.at(0)) + any_cast<double>(rhs.at(3));
        break;
      case Teuchos::MathExpr::PROD_SUB:
        result = any_cast<double>(rhs.at(0)) - any_cast<double>(rhs.at(3));
        break;
      case Teuchos::MathExpr::PROD_MUL:
        result = any_cast<double>(rhs.at(0)) * any_cast<double>(rhs.at(3));
        break;
      case Teuchos::MathExpr::PROD_DIV:
        result = any_cast<double>(rhs.at(0)) / any_cast<double>(rhs.at(3));
        break;
      case Teuchos::MathExpr::PROD_POW:
        result = std::pow(any_cast<double>(rhs.at(0)), any_cast<double>(rhs.at(3)));
        break;
      case Teuchos::MathExpr::PROD_CALL: {
        std::string& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
        CallArgs& args = Teuchos::any_ref_cast<CallArgs>(rhs.at(4));
        TEUCHOS_TEST_FOR_EXCEPTION(args.n < 1 || args.n > 2, Teuchos::ParserFail,
            "Only unary and binary functions supported!\n");
        if (args.n == 1) {
          TEUCHOS_TEST_FOR_EXCEPTION(!unary_function_map.count(name), Teuchos::ParserFail,
              "Unknown unary function name \"" << name << "\"\n");
          Unary fptr = unary_function_map[name];
          result = (*fptr)(args.a0);
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(!binary_function_map.count(name), Teuchos::ParserFail,
              "Unknown binary function name \"" << name << "\"\n");
          Binary fptr = binary_function_map[name];
          result = (*fptr)(args.a0, args.a1);
        }
        break;
      }
      case Teuchos::MathExpr::PROD_NO_ARGS: {
        CallArgs& args = Teuchos::make_any_ref<CallArgs>(result);
        args.n = 0;
        break;
      }
      case Teuchos::MathExpr::PROD_FIRST_ARG: {
        CallArgs& args = Teuchos::make_any_ref<CallArgs>(result);
        args.a0 = any_cast<double>(rhs.at(0));
        args.n = 1;
        break;
      }
      case Teuchos::MathExpr::PROD_NEXT_ARG: {
        CallArgs& args = Teuchos::any_ref_cast<CallArgs>(rhs.at(0));
        args.a1 = any_cast<double>(rhs.at(3));
        args.n = 2;
        swap(result, rhs.at(0));
        break;
      }
      case Teuchos::MathExpr::PROD_NEG:
        result = - any_cast<double>(rhs.at(2));
        break;
      case Teuchos::MathExpr::PROD_VAL_PARENS:
        result = any_cast<double>(rhs.at(2));
        break;
      case Teuchos::MathExpr::PROD_CONST:
        result = any_cast<double>(rhs.at(0));
        break;
      case Teuchos::MathExpr::PROD_VAR:
        std::string const& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
        auto it = variable_map.find(name);
        TEUCHOS_TEST_FOR_EXCEPTION(it == variable_map.end(), Teuchos::ParserFail,
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

}  // end anonymous namespace

int main() {
  Reader reader;
  for (std::string line; std::getline(std::cin, line);) {
    Teuchos::any result_any;
    try {
      reader.read_string(result_any, line, "input");
      double value = Teuchos::any_cast<double>(result_any);
      std::cout << value << '\n';
    } catch (const Teuchos::ParserFail& e) {
      std::cerr << e.what() << '\n';
    }
  }
}
