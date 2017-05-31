#ifndef TEUCHOS_REGEX_HPP
#define TEUCHOS_REGEX_HPP

#include <Teuchos_language.hpp>
#include <Teuchos_finite_automaton.hpp>
#include <Teuchos_reader_tables.hpp>
#include <Teuchos_reader.hpp>

namespace Teuchos {
namespace regex {

enum {
  PROD_REGEX,
  PROD_UNION_DECAY,
  PROD_UNION,
  PROD_CONCAT_DECAY,
  PROD_CONCAT,
  PROD_QUAL_DECAY,
  PROD_STAR,
  PROD_PLUS,
  PROD_MAYBE,
  PROD_SINGLE_CHAR,
  PROD_ANY,
  PROD_SINGLE_SET,
  PROD_PARENS_UNION,
  PROD_SET_POSITIVE,
  PROD_SET_NEGATIVE,
  PROD_POSITIVE_SET,
  PROD_NEGATIVE_SET,
  PROD_SET_ITEMS_DECAY,
  PROD_SET_ITEMS_ADD,
  PROD_SET_ITEM_CHAR,
  PROD_SET_ITEM_RANGE,
  PROD_RANGE,
};

enum { NPRODS = PROD_RANGE + 1 };

enum {
  TOK_CHAR,
  TOK_DOT,
  TOK_LRANGE,
  TOK_RRANGE,
  TOK_LPAREN,
  TOK_RPAREN,
  TOK_UNION,
  TOK_RANGE,
  TOK_NEGATE,
  TOK_STAR,
  TOK_PLUS,
  TOK_MAYBE,
};

enum { NTOKS = TOK_MAYBE + 1 };

Language build_language();
LanguagePtr ask_language();

FiniteAutomaton build_lexer();

ReaderTablesPtr ask_reader_tables();

FiniteAutomaton build_dfa(std::string const& name, std::string const& regex, int token);

any at_shift_internal(int token, std::string& text);
any at_reduce_internal(int production, std::vector<any>& rhs, int result_token);

class Reader : public Teuchos::Reader {
 public:
  Reader(int result_token_in);
  Reader(Reader const& other) = default;
  virtual ~Reader() override = default;
 protected:
  virtual any at_shift(int token, std::string& text) override;
  virtual any at_reduce(int token, std::vector<any>& rhs) override;
 private:
  int result_token;
};

}  // end namespace regex
}  // end namespace Teuchos

#endif
