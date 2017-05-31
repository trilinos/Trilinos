#ifndef TEUCHOS_LANGUAGE_HPP
#define TEUCHOS_LANGUAGE_HPP

#include <string>
#include <vector>
#include <iosfwd>

#include <Teuchos_Grammar.hpp>
#include <Teuchos_FiniteAutomaton.hpp>
#include <Teuchos_ReaderTables.hpp>
#include <Teuchos_RCP.hpp>

namespace Teuchos {

struct Language {
  struct Token {
    std::string name;
    std::string regex;
    void operator()(std::string const& name_in, std::string const& regex_in);
  };
  typedef std::vector<Token> Tokens;
  Tokens tokens;
  typedef std::vector<std::string> RHS;
  /* because C++98 doesn't allow brace initialization, here are some
     convoluted operator overloads to allow us to replace:
       prods[PROD_UNION] = {"union", {"left", "|", "right"}};
     with 
       prods[PROD_UNION]("union") >> "left", "|", "right";
   */
  struct Production;
  struct RHSBuilder {
    Production& prod;
    RHSBuilder(Production& prod_in);
    RHSBuilder& operator,(std::string const& rhs_item);
    RHSBuilder& operator>>(std::string const& rhs_item);
  };
  struct Production {
    std::string lhs;
    RHS rhs;
    RHSBuilder operator()(std::string const& lhs_in);
  };
  typedef std::vector<Production> Productions;
  Productions productions;
};

typedef RCP<const Language> LanguagePtr;

GrammarPtr build_grammar(Language const& language);

void build_lexer(FiniteAutomaton& result, Language const& language);

ReaderTablesPtr build_reader_tables(Language const& language);

std::ostream& operator<<(std::ostream& os, Language const& lang);

}

#endif
