#include <cctype>
#include <stdexcept>
#include <sstream>
#include <typeinfo>

#include <stk_expreval/Lexer.hpp>

namespace stk {
namespace expreval {

template <class T>
T convert_cast(const std::string &s)
{
  /* %TRACE% */  /* %TRACE% */
  std::istringstream is(s.c_str());
  T t = 0;

  is >> t;

  if (!is) {
    std::ostringstream msg;
    msg << "Unable to convert \"" << s << "\" to type " << typeid(T).name();
    throw std::runtime_error(msg.str().c_str());
  }
  

  return t;
}

template double convert_cast<double>(const std::string &);
template float convert_cast<float>(const std::string &);
template int convert_cast<int>(const std::string &);
template unsigned convert_cast<unsigned>(const std::string &);
template long convert_cast<long>(const std::string &);
template unsigned long convert_cast<unsigned long>(const std::string &);

namespace {

bool
is_operator(Token token) {
  return token == TOKEN_PLUS
    || token == TOKEN_MINUS
    || token == TOKEN_MULTIPLY
    || token == TOKEN_DIVIDE
    || token == TOKEN_PERCENT
    || token == TOKEN_EXPONENTIATION
    || token == TOKEN_QUESTION
    || token == TOKEN_COLON
    || token == TOKEN_SEMI
    || token == TOKEN_ASSIGN
    || token == TOKEN_LESS
    || token == TOKEN_GREATER
    || token == TOKEN_ARITHMETIC_OR
    || token == TOKEN_ARITHMETIC_AND
    || token == TOKEN_NOT
    || token == TOKEN_EQUAL
    || token == TOKEN_NOT_EQUAL
    || token == TOKEN_GREATER_EQUAL
    || token == TOKEN_LESS_EQUAL
    || token == TOKEN_LOGICAL_OR
    || token == TOKEN_LOGICAL_AND;
}


} // namespace <unnamed>


LexemVector
tokenize(
  const std::string &		expression)
{
  struct Graph
  {
    char	ch;
    Token	token;
  };

  static Graph graph[] = {
    {'+', TOKEN_PLUS},
    {'-', TOKEN_MINUS},
    {'*', TOKEN_MULTIPLY},
    {'/', TOKEN_DIVIDE},
    {'%', TOKEN_PERCENT},
    {'^', TOKEN_EXPONENTIATION},
    {'?', TOKEN_QUESTION},
    {',', TOKEN_COMMA},
    {':', TOKEN_COLON},
    {';', TOKEN_SEMI},
    {'(', TOKEN_LPAREN},
    {')', TOKEN_RPAREN},
    {'[', TOKEN_LBRACK},
    {']', TOKEN_RBRACK},
    {'=', TOKEN_ASSIGN},
    {'<', TOKEN_LESS},
    {'>', TOKEN_GREATER},
    {'|', TOKEN_ARITHMETIC_OR},
    {'&', TOKEN_ARITHMETIC_AND},
    {'!', TOKEN_NOT}
  };

  struct Digraph
  {
    char	ch1;
    char	ch2;
    Token	token;
  };

  static Digraph digraph[] = {
    {'=', '=', TOKEN_EQUAL},
    {'!', '=', TOKEN_NOT_EQUAL},
    {'>', '=', TOKEN_GREATER_EQUAL},
    {'<', '=', TOKEN_LESS_EQUAL},
    {'|', '|', TOKEN_LOGICAL_OR},
    {'&', '&', TOKEN_LOGICAL_AND}
  };

  LexemVector	lex_vector;

  const char *it = expression.c_str();

  while (*it != '\0') {
    if (std::isspace(*it) || ::iscntrl(*it))
      ++it;

    // Parse constant [0-9]*(\.[0-9*])?(E([+-]?[0-9]*))?
    //   take unary plus and minus into account
    else if (std::isdigit(*it) || *it == '.'
	     || ((*it == '-' || *it == '+')
		 && (std::isdigit(*(it + 1)) || *(it + 1) == '.')
		 && (lex_vector.empty()
		     || lex_vector.back().getToken() == TOKEN_COMMA
		     || lex_vector.back().getToken() == TOKEN_LPAREN
		     || is_operator(lex_vector.back().getToken()))))
    {
      bool is_real = false;
      const char *from = it;
      if (*it == '-' || *it == '+')
	++it;
      while (std::isdigit(*it))
	++it;
      if (*it == '.') {
	is_real = true;
	++it;
	while (std::isdigit(*it))
	  ++it;
      }
      if (*from == '.' && it == from + 1)
	throw std::runtime_error(std::string("'.' is not a valid real number ") + *it);
      if (*it == '.')
	throw std::runtime_error(std::string("'.' is not a valid real number ") + *it);
      if (std::toupper(*it) == 'E') {
	is_real = true;
	++it;
	if (*it == '+' || *it == '-') {
	  ++it;
	}
	while (std::isdigit(*it)) {
	  ++it;
	}
      }
      if (is_real)
	lex_vector.push_back(Lexem(TOKEN_REAL_CONSTANT, from, it));
      else
	lex_vector.push_back(Lexem(TOKEN_INTEGER_CONSTANT, from, it));
    }

    // Parse literal
    else if (*it == '"') {
      std::string s;
      ++it;
      for (; *it && *it != '"'; ++it) {
	if (*it == '\\') {
	  ++it;
	  if (*it)
	    switch (*it) {
	    case '"':
	      s += '"';
	      break;
	    case '\\':
	      s += '\\';
	      break;
	    case 'n':
	      s += '\n';
	      break;
	    default:
	      s += *it;
	    }
	}
	else
	  s += *it;
      }
      ++it;
      lex_vector.push_back(Lexem(TOKEN_LITERAL, s.c_str()));
    }

    // Parse identifier [a-zA-Z][a-zA-Z0-9_.]*
    else if (std::isalpha(*it)) {
      const char *from = it;
      while (std::isalpha(*it) || std::isdigit(*it) || *it == '_' || *it == '.')
	++it;
      lex_vector.push_back(Lexem(TOKEN_IDENTIFIER, from, it));
    }

    // Parse graphs and digraphs
    else if (ispunct(*it)) {
      const char *from = it;
      if (*(it + 1) != '\0') {
	for (size_t i = 0; i < sizeof(digraph)/sizeof(digraph[0]); ++i) {
	  if (*it == digraph[i].ch1 && *(it + 1) == digraph[i].ch2) {
	    ++it; ++it;
	    lex_vector.push_back(Lexem(digraph[i].token, from, it));
	    goto next_token;
	  }
	}
      }

      for (size_t i = 0; i < sizeof(graph)/sizeof(graph[0]); ++i)
	if (*it == graph[i].ch) {
	  ++it;
	  lex_vector.push_back(Lexem(graph[i].token, from, it));
	  goto next_token;
	}

      throw std::runtime_error(std::string("std::expreval::tokenize: Invalid graphic character '") + *it + "'");
    }

    else
      throw std::runtime_error("Impossible expression parse error");

  next_token:
    continue;
  }

  lex_vector.push_back(Lexem(TOKEN_END, ""));

  return lex_vector;
}

} // namespace expreval
} // namespace stk
