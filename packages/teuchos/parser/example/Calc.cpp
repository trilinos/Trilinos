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

bool operator==(CallArgs const&, CallArgs const&) {
  throw std::logic_error("CallArgs operator== is just to satisfy Teuchos");
  return false;
}

std::ostream& operator<<(std::ostream& os, CallArgs const&) {
  throw std::logic_error("CallArgs operator<< is just to satisfy Teuchos");
  return os;
}

class Reader : public Teuchos::Reader {
 public:
  Reader():Teuchos::Reader(Teuchos::MathExpr::ask_reader_tables()) {
    unary_map["sqrt"] = &std::sqrt;
    unary_map["sin"] = &std::sin;
    unary_map["cos"] = &std::cos;
    unary_map["tan"] = &std::tan;
    unary_map["asin"] = &std::asin;
    unary_map["acos"] = &std::acos;
    unary_map["atan"] = &std::atan;
    unary_map["exp"] = &std::exp;
    unary_map["log"] = &std::log;
    unary_map["log10"] = &std::log10;
    binary_map["atan2"] = &std::atan2;
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
          TEUCHOS_TEST_FOR_EXCEPTION(!unary_map.count(name), Teuchos::ParserFail,
              "Unknown unary function name \"" << name << "\"\n");
          Unary fptr = unary_map[name];
          result = (*fptr)(args.a0);
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(!binary_map.count(name), Teuchos::ParserFail,
              "Unknown binary function name \"" << name << "\"\n");
          Binary fptr = binary_map[name];
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
        throw Teuchos::ParserFail("Variables not supported!\n");
        break;
    }
  }
 private:
  typedef double (*Unary)(double);
  typedef double (*Binary)(double, double);
  std::map<std::string, Unary> unary_map;
  std::map<std::string, Binary> binary_map;
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
