#include <Teuchos_MathExpr.hpp>

namespace Teuchos {

namespace MathExpr {

Teuchos::Language make_language() {
  Teuchos::Language out;
  Teuchos::Language::Productions& prods = out.productions;
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
  prods[PROD_MUL_DIV_DECAY]("mul_div") >> "pow";
  prods[PROD_POW_DECAY]("pow") >> "neg";
  prods[PROD_NEG_DECAY]("neg") >> "scalar";
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
  prods[PROD_POW]("pow") >> "pow", "^", "S?", "neg";
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

Teuchos::ReaderTablesPtr ask_reader_tables() {
  static Teuchos::ReaderTablesPtr ptr;
  if (ptr.strong_count() == 0) {
    LanguagePtr lang = ask_language();
    ptr = Teuchos::make_reader_tables(*lang);
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

void SymbolSetReader::at_reduce(any& result, int prod, std::vector<any>& rhs) {
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

}  // end namespace MathExpr

}  // end namespace Teuchos
