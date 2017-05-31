#include <Teuchos_finite_automaton.hpp>
#include <Teuchos_string.hpp>
#include <Teuchos_language.hpp>
#include <Teuchos_parser.hpp>
#include <Teuchos_regex.hpp>
#include <Teuchos_build_parser.hpp>
#include <Teuchos_reader.hpp>
#include <Teuchos_xml.hpp>
#include <Teuchos_yaml.hpp>

#include <iostream>
#include <sstream>

using namespace Teuchos;

static bool accepts(FiniteAutomaton const& fa, std::string const& s, int token = 0) {
  int state = 0;
  for (int i = 0; i < size(s); ++i) {
    auto c = at(s, i);
    if (!is_symbol(c)) {
      return false;
    }
    auto symbol = get_symbol(c);
    state = step(fa, state, symbol);
    if (state == -1) return false;
  }
  return accepts(fa, state) == token;
}

static void test_finite_automaton() {
  auto lower = make_char_range_nfa('a', 'z');
  assert(accepts(lower, "a"));
  assert(accepts(lower, "q"));
  assert(accepts(lower, "z"));
  assert(!accepts(lower, "X"));
  assert(!accepts(lower, "246"));
  assert(!accepts(lower, "abc"));
  assert(!accepts(lower, "\xff"));
  auto single = make_char_single_nfa('q');
  assert(!accepts(single, "a"));
  assert(accepts(single, "q"));
  assert(!accepts(single, "r"));
  assert(!accepts(single, "abc"));
  auto upper = make_char_range_nfa('A', 'Z');
  auto alpha_nfa = unite(lower, upper);
  auto alpha_dfa = make_deterministic(alpha_nfa);
  assert(get_nstates(alpha_dfa) > 2);
  auto alpha = simplify(alpha_dfa);
  assert(get_nstates(alpha) == 2);
  {
    auto num = make_char_range_nfa('0', '9');
    auto under = make_char_single_nfa('_');
    auto under_alpha = unite(under, alpha);
    auto under_alnum = unite(under_alpha, num);
    auto under_alnum_star = star(under_alnum);
    auto ident_nfa = concat(under_alpha, under_alnum_star);
    auto ident_dfa = make_deterministic(ident_nfa);
    assert(get_nstates(ident_dfa) > 2);
    auto identifier = simplify(ident_dfa);
    assert(get_nstates(identifier) == 2);
    assert(accepts(identifier, "_soup"));
    assert(!accepts(identifier, "007"));
    assert(accepts(identifier, "All_The_Case"));
    assert(!accepts(identifier, "fire|hose"));
  }
}

static void test_lr0_language() {
  /* The LR(0) grammar $G_1$ on page 20 of Pager's paper */ 
  Language lang;
  lang.tokens.push_back({"a", "a"});
  lang.tokens.push_back({"(", "\\("});
  lang.tokens.push_back({")", "\\)"});
  lang.tokens.push_back({"+", "\\+"});
  lang.productions.push_back({"G", {"(", "E", ")"}});
  lang.productions.push_back({"E", {"E", "+", "T"}});
  lang.productions.push_back({"E", {"T"}});
  lang.productions.push_back({"T", {"(", "E", ")"}});
  lang.productions.push_back({"T", {"a"}});
  auto grammar = build_grammar(lang);
  auto parser = accept_parser(build_lalr1_parser(grammar));
}

static void test_lalr1_language() {
  /* this is a pretty simple language I found in some
     university's homework assignment which is LALR(1)
     but not SLR(1), so it is a good first test of our
     context set derivation */
  Language lang;
  lang.tokens.push_back({"a", "a"});
  lang.tokens.push_back({"b", "b"});
  lang.tokens.push_back({"c", "c"});
  lang.tokens.push_back({"d", "d"});
  lang.productions.push_back({"S", {"A", "a"}});
  lang.productions.push_back({"S", {"b", "A", "c"}});
  lang.productions.push_back({"S", {"d", "c"}});
  lang.productions.push_back({"S", {"b", "d", "a"}});
  lang.productions.push_back({"A", {"d"}});
  auto grammar = build_grammar(lang);
  auto parser = accept_parser(build_lalr1_parser(grammar));
}

static void test_regex_lexer() {
  auto lexer = regex::build_lexer();
  assert(accepts(lexer, "a", regex::TOK_CHAR));
  assert(accepts(lexer, ".", regex::TOK_DOT));
  assert(accepts(lexer, "?", regex::TOK_MAYBE));
  assert(accepts(lexer, "\\.", regex::TOK_CHAR));
  assert(accepts(lexer, "\\?", regex::TOK_CHAR));
}

static void test_regex_language() {
  auto lang = regex::ask_language();
  auto grammar = build_grammar(*lang);
  auto parser = accept_parser(build_lalr1_parser(grammar));
}

static void test_regex_reader(std::string const& regex,
    std::vector<std::string> const& expect_matches,
    std::vector<std::string> const& expect_non_matches) {
  auto reader = regex::Reader(42);
  bool ok = reader.read_string("test_regex_reader", regex);
  assert(ok);
  auto result = reader.move_result<FiniteAutomaton>();
  for (auto& expect_match : expect_matches) {
    assert(accepts(result, expect_match, 42));
  }
  for (auto& expect_non_match : expect_non_matches) {
    assert(!accepts(result, expect_non_match, 42));
  }
}

static void test_regex_reader() {
  test_regex_reader("a", {"a"}, {"b", "aa"});
  test_regex_reader("ab", {"ab"}, {"a", "b", "abc"});
  test_regex_reader("abc", {"abc"}, {"a", "ab", "abd", "abcd"});
  test_regex_reader("a?", {"a", ""}, {"b", "aa"});
  test_regex_reader("[a]", {"a"}, {"b", "aa"});
  test_regex_reader("[^a]", {"q", "z"}, {"a", "qq"});
  test_regex_reader("[A-Z]", {"A", "Q", "Z"}, {"a", "z", "AA", ">"});
  test_regex_reader("\\?", {"?"}, {"b", "??"});
  test_regex_reader("a+", {"a", "aa", "aaa"}, {"", "??"});
  test_regex_reader("a*", {"", "a", "aa", "aaa"}, {"??", "b"});
  test_regex_reader("foo|bar", {"foo", "bar"}, {"??", "foobar"});
  test_regex_reader("[a-zA-Z0-9]",
      {"a", "q", "z", "A", "Q", "Z", "0", "5", "9"},
      {"?", "_", "-", "+"});
  test_regex_reader("[_a-zA-Z][_a-zA-Z0-9]*",
      {"_soup", "All_The_Case", "l33t"},
      {"007", "fire|hose", "+1"});
  test_regex_reader("\t \n\r", {"\t \n\r"}, {"abc", "  \n\r"});
  test_regex_reader("ba(na)+", {"bana", "banana", "bananana"}, {"banan", "banane"});
  /* detect a C-style comment like this one */
  test_regex_reader("/\\*[^\\*]*\\*+([^/\\*][^\\*]*\\*+)*/",
      {"/**/", "/*\n *\n */", "/* /* /***/"},
      {"/*/"});
  /* <!-- detect an XML-style comment like this one -->
   * note: after reading the XML spec, it seems they don't allow
   *       comments to be this general (contain --), and furthermore it doesn't
   *       seem like the right approach to detect comments in the lexer
   *       because XML's character meanings change so much between comments,
   *       attribute values, and elsewhere.
   *       as such, the following is more of an academic exercise in how to generalize
   *       the C comment rule above, but don't expect to see it in Teuchos::xml */
  test_regex_reader("<!\\-\\-([^\\-]*\\-([^\\-]+\\-)*)\\-([^\\->]([^\\-]*\\-([^\\-]+\\-)*)\\-)*>",
      {"<!-- foo -->", "<!-- <!--\nfoo bar\n>-- -->"},
      {"<!-- --> -->"});
  test_regex_reader("(0|([1-9][0-9]*))(\\.[0-9]*)?([eE]\\-?[1-9][0-9]*)?",
      {"1", "1.0", "1e6", "3.14159", "2.2e3", "0.0", "0.0001"}, {"a", "-1"});
}

static void test_xml_language() {
  auto lang = xml::ask_language();
  auto grammar = build_grammar(*lang);
  auto pip = build_lalr1_parser(grammar);
}

static void test_xml_reader(std::string const& str) {
  auto reader = Reader(xml::ask_reader_tables());
  bool ok = reader.read_string("test_xml_reader", str);
  assert(ok);
}

static void test_xml_reader() {
  test_xml_reader("<the_tag100/>");
  test_xml_reader("<Parameter name=\"force\"/>");
  test_xml_reader("<Parameter name=\"force\"\ttype=\"double\" value=\"1.9\"/>");
  test_xml_reader("<ParameterList name=\"Physics Vars\">\n  <Parameter name=\"force\"\ttype=\"double\" value=\"1.9\"/>\n</ParameterList>");
  test_xml_reader("<P name=\"foo&quot;&#72;bar\"/>");
}

static void test_yaml_language() {
  auto lang = yaml::ask_language();
  auto grammar = build_grammar(*lang);
  auto pip = build_lalr1_parser(grammar);
}

static void test_yaml_reader(std::string const& str) {
  auto reader = Reader(yaml::ask_reader_tables());
  bool ok = reader.read_string("test_yaml_reader", str);
  assert(ok);
}

static void test_yaml_reader() {
  test_yaml_reader("---\nfoo:bar\n...\n");
  test_yaml_reader("%YAML 1.2\n---\nfoo:bar\n...\n");
  test_yaml_reader("---\nfoo:bar\nfar:boo\n...\n");
  test_yaml_reader("---\nfoo:\n  bar:42\n  baz:  100\n...\n");
  test_yaml_reader("---\nfoo:\n  bar  :42\n  baz:  100\n...\n");
  test_yaml_reader("---\n foo:bar\n...\n");
  test_yaml_reader("---\n\"Don Knuth\":bar\n...\n");
  test_yaml_reader("---\n\"Don \\\"Maverick\\\" Knuth\":bar\n...\n");
  test_yaml_reader("---\n'never say never':true\n...\n");
  test_yaml_reader("---\n'never say ''never''':true\n...\n");
  test_yaml_reader("---\n - foo\n - bar\n...\n");
  test_yaml_reader("---\n1:\n - a\n - b\n...\n");
  test_yaml_reader("---\na: {1: 1, 2: 4, 3: 9}\n...\n");
  test_yaml_reader("---\na: [1, 2, 3]\n...\n");
  test_yaml_reader("---\na: {1: [0, 1], 2: [0, 1, 2]}\n...\n");
  test_yaml_reader("---\nassocs: [[bottom, 1], [top, 42]]\n...\n");
  test_yaml_reader("---\npressure: -1.9e-6\nvolume: 0.7e+10\n...\n");
  auto reader = DebugReader(yaml::ask_reader_tables());
  bool ok = reader.read_string("debug_yaml_reader", "---\npressure: -1.9e-6\nvolume: 0.7e+10\n...\n");
  assert(ok);
}

int main() {
  std::string a("  ");
  std::string b("");
  assert(0 == a.compare(0, 0, b));
  test_finite_automaton();
  test_lr0_language();
  test_lalr1_language();
  test_regex_lexer();
  test_regex_language();
  test_regex_reader();
  test_xml_language();
  test_xml_reader();
  test_yaml_language();
  test_yaml_reader();
}
