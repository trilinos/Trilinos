#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <Teuchos_language.hpp>
#include <Teuchos_reader.hpp>

namespace {

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
  PROD_SPACES,
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
  TOK_CONST,
};

enum { NTOKS = TOK_CONST + 1 };

Teuchos::Language build_language() {
  Teuchos::Language out;
  auto& prods = out.productions;
  prods.resize(NPRODS);
  prods[PROD_EXPR] = {"expr", {"expr(+-)"}};
  prods[PROD_ADD_SUB_DECAY] = {"expr(+-)", {"expr(*/)"}};
  prods[PROD_MUL_DIV_DECAY] = {"expr(*/)", {"expr(^)"}};
  prods[PROD_POW_DECAY] = {"expr(^)", {"neg-expr"}};
  prods[PROD_NEG_DECAY] = {"neg-expr", {"scalar-expr"}};
  prods[PROD_ADD] = {"expr(+-)", {"expr(+-)", "+", "S?", "expr(*/)"}};
  prods[PROD_SUB] = {"expr(+-)", {"expr(+-)", "-", "S?", "expr(*/)"}};
  prods[PROD_MUL] = {"expr(*/)", {"expr(*/)", "*", "S?", "expr(^)"}};
  prods[PROD_DIV] = {"expr(*/)", {"expr(*/)", "/", "S?", "expr(^)"}};
  prods[PROD_POW] = {"expr(^)", {"expr(^)", "^", "S?", "neg-expr"}};
  prods[PROD_NEG] = {"neg-expr", {"-", "scalar-expr"}};
  prods[PROD_UNARY_CALL] = {"scalar-expr", {"name", "S?", "(", "S?", "expr(+-)", ")", "S?"}};
  prods[PROD_BINARY_CALL] = {"scalar-expr", {"name", "S?", "(", "S?", "expr(+-)", ",", "S?", "expr(+-)", ")", "S?"}};
  prods[PROD_PARENS] = {"scalar-expr", {"(", "S?", "expr(+-)", ")", "S?"}};
  prods[PROD_CONST] = {"scalar-expr", {"constant", "S?"}};
  prods[PROD_NO_SPACES] = {"S?", {}};
  prods[PROD_SPACES] = {"S?", {"spaces"}};
  out.tokens.resize(NTOKS);
  out.tokens[TOK_SPACE] = {"spaces", "[ \t\n\r]+"};
  out.tokens[TOK_NAME] = {"name", "[_a-zA-Z][_a-zA-Z0-9]*"};
  out.tokens[TOK_ADD] = {"+", "\\+"};
  out.tokens[TOK_SUB] = {"-", "\\-"};
  out.tokens[TOK_MUL] = {"*", "\\*"};
  out.tokens[TOK_DIV] = {"/", "\\/"};
  out.tokens[TOK_POW] = {"^", "\\^"};
  out.tokens[TOK_LPAREN] = {"(", "\\("};
  out.tokens[TOK_RPAREN] = {")", "\\)"};
  out.tokens[TOK_COMMA] = {",", ","};
  out.tokens[TOK_CONST] = {"constant", "(0|([1-9][0-9]*))(\\.[0-9]*)?([eE]\\-?[1-9][0-9]*)?"};
  return out;
}

Teuchos::ReaderTablesPtr ask_reader_tables() {
  static Teuchos::ReaderTablesPtr ptr;
  if (ptr.use_count() == 0) {
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
    unary_map["log2"] = &std::log2;
    unary_map["cbrt"] = &std::cbrt;
    binary_map["hypot"] = &std::hypot;
    binary_map["atan2"] = &std::atan2;
  }
  Reader(Reader const& other) = default;
  virtual ~Reader() override = default;
 protected:
  virtual Teuchos::any at_shift(int token, std::string& text) override {
    switch (token) {
      case TOK_NAME: return Teuchos::any(std::move(text));
      case TOK_CONST: {
        return Teuchos::any(std::atof(text.c_str()));
      }
    }
    return Teuchos::any();
  }
  virtual bool at_reduce(int prod, std::vector<Teuchos::any>& rhs, Teuchos::any& result) override {
    switch (prod) {
      case PROD_EXPR:
      case PROD_ADD_SUB_DECAY:
      case PROD_MUL_DIV_DECAY:
      case PROD_POW_DECAY:
      case PROD_NEG_DECAY:
        result = std::move(rhs.at(0));
        break;
      case PROD_ADD:
        result = Teuchos::move_value<double>(rhs.at(0)) + Teuchos::move_value<double>(rhs.at(3));
        break;
      case PROD_SUB:
        result = Teuchos::move_value<double>(rhs.at(0)) - Teuchos::move_value<double>(rhs.at(3));
        break;
      case PROD_MUL:
        result = Teuchos::move_value<double>(rhs.at(0)) * Teuchos::move_value<double>(rhs.at(3));
        break;
      case PROD_DIV:
        result = Teuchos::move_value<double>(rhs.at(0)) / Teuchos::move_value<double>(rhs.at(3));
        break;
      case PROD_POW:
        result = std::pow(Teuchos::move_value<double>(rhs.at(0)), Teuchos::move_value<double>(rhs.at(3)));
        break;
      case PROD_NEG:
        result = - Teuchos::move_value<double>(rhs.at(1));
        break;
      case PROD_UNARY_CALL: {
        auto name = Teuchos::move_value<std::string>(rhs.at(0));
        auto arg = Teuchos::move_value<double>(rhs.at(4));
        if (!unary_map.count(name)) {
          std::cerr << "Unknown unary function name \"" << name << "\"\n";
          return false;
        }
        auto fptr = unary_map[name];
        result = (*fptr)(arg);
        break;
      }
      case PROD_BINARY_CALL: {
        auto name = Teuchos::move_value<std::string>(rhs.at(0));
        auto arg1 = Teuchos::move_value<double>(rhs.at(4));
        auto arg2 = Teuchos::move_value<double>(rhs.at(7));
        if (!binary_map.count(name)) {
          std::cerr << "Unknown binary function name \"" << name << "\"\n";
          return false;
        }
        auto fptr = binary_map[name];
        result = (*fptr)(arg1, arg2);
        break;
      }
      case PROD_PARENS:
        result = std::move(rhs.at(2));
        break;
      case PROD_CONST:
        result = std::move(rhs.at(0));
        break;
    }
    return true;
  }
 private:
  typedef double (*Unary)(double);
  typedef double (*Binary)(double, double);
  std::map<std::string, Unary> unary_map;
  std::map<std::string, Binary> binary_map;
};

}  // end anonymous namespace

int main() {
  auto reader = Reader();
  for (std::string line; std::getline(std::cin, line);) {
    bool ok = reader.read_string("input", line);
    if (ok) {
      auto value = reader.move_result<double>();
      std::cout << value << '\n';
    }
  }
}
