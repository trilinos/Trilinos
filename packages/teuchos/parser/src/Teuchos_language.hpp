#ifndef TEUCHOS_LANGUAGE_HPP
#define TEUCHOS_LANGUAGE_HPP

#include <map>
#include <string>
#include <vector>
#include <iosfwd>

#include <Teuchos_grammar.hpp>
#include <Teuchos_finite_automaton.hpp>
#include <Teuchos_reader_tables.hpp>

namespace Teuchos {

struct Language {
  struct Token {
    std::string name;
    std::string regex;
  };
  std::vector<Token> tokens;
  struct Production {
    std::string lhs;
    std::vector<std::string> rhs;
  };
  std::vector<Production> productions;
};

using LanguagePtr = std::shared_ptr<Language>;

GrammarPtr build_grammar(Language const& language);

FiniteAutomaton build_lexer(Language const& language);

ReaderTablesPtr build_reader_tables(Language const& language);

std::ostream& operator<<(std::ostream& os, Language const& lang);

}

#endif
