#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <Teuchos_Language.hpp>
#include <Teuchos_Reader.hpp>

namespace {

using Teuchos::any;
using Teuchos::any_cast;

enum {
  PROD_EXPR,
  PROD_ADD_SUB_DECAY,
  PROD_MUL_DIV_DECAY,
  PROD_POW_DECAY,
  PROD_NEG_DECAY,
  PROD_ADD,
  PROD_SUB,
  PROD_MUL,
  PROD_DIV,
  PROD_POW,
  PROD_UNARY_CALL,
  PROD_BINARY_CALL,
  PROD_PARENS,
  PROD_CONST,
  PROD_NEG,
  PROD_NO_SPACES,
  PROD_SPACES
};

enum { NPRODS = PROD_SPACES + 1 };

enum {
  TOK_SPACE,
  TOK_NAME,
  TOK_ADD,
  TOK_SUB,
  TOK_MUL,
  TOK_DIV,
  TOK_POW,
  TOK_LPAREN,
  TOK_RPAREN,
  TOK_COMMA,
  TOK_CONST
};

enum { NTOKS = TOK_CONST + 1 };

Teuchos::Language build_language() {
  Teuchos::Language out;
  Teuchos::Language::Productions& prods = out.productions;
  prods.resize(NPRODS);
  prods[PROD_EXPR]("expr") >> "expr(+-)";
  prods[PROD_ADD_SUB_DECAY]("expr(+-)") >> "expr(*/)";
  prods[PROD_MUL_DIV_DECAY]("expr(*/)") >> "expr(^)";
  prods[PROD_POW_DECAY]("expr(^)") >> "neg-expr";
  prods[PROD_NEG_DECAY]("neg-expr") >> "scalar-expr";
  prods[PROD_ADD]("expr(+-)") >> "expr(+-)", "+", "S?", "expr(*/)";
  prods[PROD_SUB]("expr(+-)") >> "expr(+-)", "-", "S?", "expr(*/)";
  prods[PROD_MUL]("expr(*/)") >> "expr(*/)", "*", "S?", "expr(^)";
  prods[PROD_DIV]("expr(*/)") >> "expr(*/)", "/", "S?", "expr(^)";
  prods[PROD_POW]("expr(^)") >> "expr(^)", "^", "S?", "neg-expr";
  prods[PROD_NEG]("neg-expr") >> "-", "scalar-expr";
  prods[PROD_UNARY_CALL]("scalar-expr") >> "name", "S?", "(", "S?", "expr(+-)", ")", "S?";
  prods[PROD_BINARY_CALL]("scalar-expr") >> "name", "S?", "(", "S?", "expr(+-)", ",", "S?", "expr(+-)", ")", "S?";
  prods[PROD_PARENS]("scalar-expr") >> "(", "S?", "expr(+-)", ")", "S?";
  prods[PROD_CONST]("scalar-expr") >> "constant", "S?";
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
  out.tokens[TOK_CONST]("constant", "(0|([1-9][0-9]*))(\\.[0-9]*)?([eE]\\-?[1-9][0-9]*)?");
  return out;
}

Teuchos::ReaderTablesPtr ask_reader_tables() {
  static Teuchos::ReaderTablesPtr ptr;
  if (ptr.strong_count() == 0) {
    ptr = Teuchos::build_reader_tables(build_language());
  }
  return ptr;
}

class Reader : public Teuchos::Reader {
 public:
  Reader():Teuchos::Reader(ask_reader_tables()) {
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
      case TOK_NAME: {
        std::string& result = Teuchos::make_any_ref<std::string>(result_any);
        swap(result, text);
        return;
      }
      case TOK_CONST: {
        result_any = std::atof(text.c_str());
        return;
      }
    }
  }
  virtual void at_reduce(Teuchos::any& result, int prod, std::vector<Teuchos::any>& rhs) {
    using std::swap;
    switch (prod) {
      case PROD_EXPR:
      case PROD_ADD_SUB_DECAY:
      case PROD_MUL_DIV_DECAY:
      case PROD_POW_DECAY:
      case PROD_NEG_DECAY:
        swap(result, rhs.at(0));
        break;
      case PROD_ADD:
        result = any_cast<double>(rhs.at(0)) + any_cast<double>(rhs.at(3));
        break;
      case PROD_SUB:
        result = any_cast<double>(rhs.at(0)) - any_cast<double>(rhs.at(3));
        break;
      case PROD_MUL:
        result = any_cast<double>(rhs.at(0)) * any_cast<double>(rhs.at(3));
        break;
      case PROD_DIV:
        result = any_cast<double>(rhs.at(0)) / any_cast<double>(rhs.at(3));
        break;
      case PROD_POW:
        result = std::pow(any_cast<double>(rhs.at(0)), any_cast<double>(rhs.at(3)));
        break;
      case PROD_NEG:
        result = - any_cast<double>(rhs.at(1));
        break;
      case PROD_UNARY_CALL: {
        std::string& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
        double arg = any_cast<double>(rhs.at(4));
        TEUCHOS_TEST_FOR_EXCEPTION(!unary_map.count(name), Teuchos::ParserFail,
            "Unknown unary function name \"" << name << "\"\n");
        Unary fptr = unary_map[name];
        result = (*fptr)(arg);
        break;
      }
      case PROD_BINARY_CALL: {
        std::string& name = Teuchos::any_ref_cast<std::string>(rhs.at(0));
        double arg1 = any_cast<double>(rhs.at(4));
        double arg2 = any_cast<double>(rhs.at(7));
        TEUCHOS_TEST_FOR_EXCEPTION(!binary_map.count(name), Teuchos::ParserFail,
            "Unknown binary function name \"" << name << "\"\n");
        Binary fptr = binary_map[name];
        result = (*fptr)(arg1, arg2);
        break;
      }
      case PROD_PARENS:
        result = any_cast<double>(rhs.at(2));
        break;
      case PROD_CONST:
        result = any_cast<double>(rhs.at(0));
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
