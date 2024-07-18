// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <iomanip>
#include <ios>
#include <sstream>
#include <cctype>
#include <fstream>

#include "Teuchos_YamlParser_decl.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_YamlParameterListCoreHelpers.hpp"
#include "Teuchos_TwoDArray.hpp"

#include "Teuchos_Reader.hpp"

#ifdef HAVE_TEUCHOSPARAMETERLIST_YAMLCPP
#include "yaml-cpp/yaml.h"
#endif // HAVE_TEUCHOSPARAMETERLIST_YAMLCPP
#include "Teuchos_YAML.hpp"


namespace Teuchos {

#ifdef HAVE_TEUCHOSPARAMETERLIST_YAMLCPP

/* see https://github.com/jbeder/yaml-cpp/issues/261
   there are times when we want to insist that a parameter
   value be interpreted as a string despite it being parseable
   as a number.
   the standard way to do this in YAML is to put the number in quotes,
   i.e. '1e-3' instead of 1e-3.
   however, the usual YAML::Node::as<T> system doesn't respect quoting
   when trying to cast to numbers.
   so, this is our own version of as<T>, called quoted_as<T>, using
   the Tag workaround suggested in the issue linked above. */

template <typename T>
struct QuotedAs {
  static T eval(::YAML::Node const& node) {
    // this "!" tag apparently denotes that the value was quoted
    if (node.Tag() == "!") {
      throw std::runtime_error("quoted_as from quoted string to number");
    }
    return node.as<T>();
  }
};

template <>
struct QuotedAs<std::string> {
  // only a cast to string will succeed if quoted
  static std::string eval(::YAML::Node const& node) { return node.as<std::string>(); }
};

template <typename T>
static T quoted_as(::YAML::Node const& node) { return QuotedAs<T>::eval(node); }

template<typename T>
Teuchos::Array<T> getYamlArray(const ::YAML::Node& node)
{
  Teuchos::Array<T> arr;
  for(::YAML::const_iterator it = node.begin(); it != node.end(); it++)
  {
    arr.push_back(quoted_as<T>(*it));
  }
  return arr;
}

bool checkYamlTwoDArrayIsRagged(const ::YAML::Node& node)
{
  bool ragged = false;
  for (::YAML::const_iterator it = node.begin(); it != node.end(); ++it)
  {
    if (it->size() != node.begin()->size())
    {
      ragged=true;
    }
  }
  return ragged;
}

template<typename T> Teuchos::TwoDArray<T> getYamlTwoDArray(const ::YAML::Node& node)
{
  Teuchos::TwoDArray<T> arr;
  typename Teuchos::TwoDArray<T>::size_type i, j;
  arr.resizeRows(node.size());
  arr.resizeCols(node.begin()->size());
  i = 0;
  for (::YAML::const_iterator rit = node.begin(); rit != node.end(); ++rit)
  {
    j = 0;
    for (::YAML::const_iterator cit = rit->begin(); cit != rit->end(); ++cit)
    {
      arr(i, j) = quoted_as<T>(*cit);
      ++j;
    }
   ++i;
  }
  return arr;
}

int getYamlArrayDim(const ::YAML::Node& node)
{
  int ndim = 0;
  if (node.Type() == ::YAML::NodeType::Sequence)
  {
    ++ndim;
    if (node.begin()->Type() == ::YAML::NodeType::Sequence)
    {
      ++ndim;
      if (node.begin()->begin()->Type() == ::YAML::NodeType::Sequence)
      {
        ++ndim;
      }
    }
  }
  return ndim;
}

template <typename tarray_t, typename T>
tarray_t getYaml2DRaggedArray(::YAML::Node node, int ndim, std::string key)
{
  tarray_t base_arr;
  if (ndim == 2) {
    Teuchos::Array<T> sub_arr;
    for (::YAML::const_iterator it1 = node.begin(); it1 != node.end(); ++it1) {
      for (::YAML::const_iterator it2 = it1->begin(); it2 != it1->end(); ++it2) {
        sub_arr.push_back(quoted_as<T>(*it2));
      } base_arr.push_back(sub_arr);
    sub_arr.clear();
    }
  }
  else
  {
    throw YamlSequenceError(std::string("MDArray \"" + key + "\" must have dim 2."));
  }
  return base_arr;
}

// This handles the requested use case of a list of 2D arrays; further nesting would require a getYaml4DArray() function,
// which could be straightforwardly implemented along the lines of the below function.

template <typename tarray_t, typename T>
tarray_t getYaml3DArray(::YAML::Node node, int ndim, std::string key)
{
  tarray_t base_arr;
  if (ndim == 3) {
    Teuchos::Array<Teuchos::Array<T>> sub_arr;
    Teuchos::Array<T> sub_sub_arr;
    for (::YAML::const_iterator it1 = node.begin(); it1 != node.end(); ++it1) {
      for (::YAML::const_iterator it2 = it1->begin(); it2 != it1->end(); ++it2) {
        for (::YAML::const_iterator it3 = it2->begin(); it3 != it2->end(); ++it3) {
          sub_sub_arr.push_back(quoted_as<T>(*it3));
        } sub_arr.push_back(sub_sub_arr);
      sub_sub_arr.clear();
      } base_arr.push_back(sub_arr);
    sub_arr.clear();
    }
  }
  else
  {
    throw YamlSequenceError(std::string("MDArray \"" + key + "\" must have dim 3."));
  }
  return base_arr;
}

template <typename T>
void safe_set_entry(ParameterList& list, std::string const& name_in, T const& entry_in) {
  TEUCHOS_TEST_FOR_EXCEPTION(list.isParameter(name_in), ParserFail,
      "Parameter \"" << name_in << "\" already exists in list \"" << list.name() << "\"\n");
  list.set(name_in, entry_in);
}
#endif // HAVE_TEUCHOSPARAMETERLIST_YAMLCPP

std::string remove_trailing_whitespace(std::string const& in) {
  std::size_t new_end = 0;
  for (std::size_t ri = 0; ri < in.size(); ++ri) {
    std::size_t i = in.size() - 1 - ri;
    if (in[i] != ' ' && in[i] != '\t') {
      new_end = i + 1;
      break;
    }
  }
  return in.substr(0, new_end);
}

std::string remove_trailing_whitespace_and_newlines(std::string const& in) {
  std::size_t new_end = 0;
  for (std::size_t ri = 0; ri < in.size(); ++ri) {
    std::size_t i = in.size() - 1 - ri;
    if (in[i] != ' ' && in[i] != '\t' && in[i] != '\n' && in[i] != '\r') {
      new_end = i + 1;
      break;
    }
  }
  return in.substr(0, new_end);
}

template <typename T>
bool is_parseable_as(std::string const& text) {
  std::istringstream ss(text);
  T val;
  ss >> std::noskipws >> val;
  return ss.eof() && !ss.fail();
}

template <>
bool is_parseable_as<int>(std::string const& text) {
  std::istringstream ss(text);
  using LL = long long;
  LL val;
  ss >> std::noskipws >> val;
  return ss.eof() && !ss.fail() &&
    (val >= LL(std::numeric_limits<int>::min())) &&
    (val <= LL(std::numeric_limits<int>::max()));
}

template <typename T>
T parse_as(std::string const& text) {
  std::istringstream ss(text);
  T value;
  ss >> value;
  return value;
}

// http://en.cppreference.com/w/cpp/string/byte/tolower
static char my_tolower(char ch)
{
  return std::tolower(static_cast<unsigned char>(ch));
}

// http://en.cppreference.com/w/cpp/string/byte/isdigit
static bool my_isdigit(char ch)
{
  return std::isdigit(static_cast<unsigned char>(ch));
}

template <>
bool is_parseable_as<bool>(std::string const& text) {
  std::string lower;
  for (std::size_t i = 0; i < text.size(); ++i) {
    lower.push_back(my_tolower(text[i]));
  }
  return lower == "true" || lower == "yes" ||
         lower == "false" || lower == "no";
}

template <>
bool parse_as<bool>(std::string const& text) {
  std::string lower;
  for (std::size_t i = 0; i < text.size(); ++i) {
    lower.push_back(my_tolower(text[i]));
  }
  return !(lower == "false" || lower == "no");
}

struct PLPair {
  std::string key;
  ParameterEntry value;
};

struct Scalar {
  enum Source {
    RAW,
    DQUOTED,
    SQUOTED,
    BLOCK
  };
  /* order matters, a higher type should be convertible to a lower type */
  enum Type {
    STRING = 0,
    DOUBLE = 1,
    LONG_LONG = 2,
    INT    = 3,
    BOOL   = 4
  };
  int source;
  int tag_type;
  std::string text;
  int infer_type() const {
    if (tag_type != -1) {
      return tag_type;
    }
    if (source != RAW) {
      return STRING;
    }
    if (is_parseable_as<bool>(text)) {
      return BOOL;
    }
    if (is_parseable_as<int>(text)) {
      return INT;
    }
    if (is_parseable_as<long long>(text)) {
      return LONG_LONG;
    }
    if (is_parseable_as<double>(text)) {
      return DOUBLE;
    }
    return STRING;
  }
};

bool operator==(Scalar const&, Scalar const&) { return false; }
std::ostream& operator<<(std::ostream& os, Scalar const&) { return os; }

void safe_set_entry(ParameterList& list, std::string const& name_in, ParameterEntry const& entry_in) {
  TEUCHOS_TEST_FOR_EXCEPTION(list.isParameter(name_in), ParserFail,
      "Parameter \"" << name_in << "\" already exists in list \"" << list.name() << "\"\n");
  list.setEntry(name_in, entry_in);
}

namespace YAMLParameterList {

class Reader : public Teuchos::Reader {
 public:
  Reader():Teuchos::Reader(Teuchos::YAML::ask_reader_tables()) {}
  virtual ~Reader() {}
 protected:
  enum {
    TRIM_NORMAL,
    TRIM_DASH
  };
  virtual void at_shift(any& result_any, int token, std::string& text) {
    using std::swap;
    switch (token) {
      case Teuchos::YAML::TOK_NEWLINE: {
        std::string& result = make_any_ref<std::string>(result_any);
        swap(result, text);
        break;
      }
      case Teuchos::YAML::TOK_SPACE:
      case Teuchos::YAML::TOK_OTHER: {
        result_any = text.at(0);
        break;
      }
    }
  }
  virtual void at_reduce(any& result_any, int prod, std::vector<any>& rhs) {
    using std::swap;
    switch (prod) {
      case Teuchos::YAML::PROD_DOC:
      case Teuchos::YAML::PROD_DOC2: {
        std::size_t offset = prod == Teuchos::YAML::PROD_DOC2 ? 1 : 0;
        TEUCHOS_ASSERT(rhs.at(offset).has_value());
        swap(result_any, rhs.at(offset));
        TEUCHOS_ASSERT(result_any.type() == typeid(ParameterList));
        break;
      }
      case Teuchos::YAML::PROD_TOP_BMAP: {
        TEUCHOS_ASSERT(rhs.at(0).has_value());
        TEUCHOS_ASSERT(rhs.at(0).type() == typeid(PLPair));
        PLPair& pair = any_ref_cast<PLPair>(rhs.at(0));
        any& pair_rhs_any = pair.value.getAny(/* set isUsed = */ false);
        result_any = pair_rhs_any;
        break;
      }
      case Teuchos::YAML::PROD_TOP_FIRST: {
        if (rhs.at(0).type() == typeid(ParameterList)) {
          swap(result_any, rhs.at(0));
        }
        break;
      }
      case Teuchos::YAML::PROD_TOP_NEXT: {
        if (rhs.at(1).type() == typeid(ParameterList)) {
          TEUCHOS_TEST_FOR_EXCEPTION(rhs.at(0).has_value(), ParserFail,
              "Can't specify multiple top-level ParameterLists in one YAML file!\n");
          swap(result_any, rhs.at(1));
        } else {
          swap(result_any, rhs.at(0));
        }
        break;
      }
      case Teuchos::YAML::PROD_BMAP_FIRST:
      case Teuchos::YAML::PROD_FMAP_FIRST: {
        TEUCHOS_ASSERT(rhs.at(0).type() == typeid(PLPair));
        map_first_item(result_any, rhs.at(0));
        TEUCHOS_ASSERT(result_any.type() == typeid(ParameterList));
        break;
      }
      case Teuchos::YAML::PROD_BMAP_NEXT: {
        map_next_item(result_any, rhs.at(0), rhs.at(1));
        break;
      }
      case Teuchos::YAML::PROD_FMAP_NEXT: {
        map_next_item(result_any, rhs.at(0), rhs.at(3));
        break;
      }
      case Teuchos::YAML::PROD_BMAP_SCALAR:
      case Teuchos::YAML::PROD_FMAP_SCALAR:
      case Teuchos::YAML::PROD_FMAP_FMAP:
      case Teuchos::YAML::PROD_FMAP_FSEQ: {
        int scalar_type = interpret_tag(rhs.at(3));
        map_item(result_any, rhs.at(0), rhs.at(4), scalar_type);
        break;
      }
      case Teuchos::YAML::PROD_BMAP_BSCALAR: {
        map_item(result_any, rhs.at(0), rhs.at(3), Scalar::STRING);
        break;
      }
      case Teuchos::YAML::PROD_BMAP_BVALUE: {
        map_item(result_any, rhs.at(0), rhs.at(4));
        break;
      }
      case Teuchos::YAML::PROD_BVALUE_EMPTY: {
        result_any = ParameterList();
        break;
      }
      case Teuchos::YAML::PROD_BVALUE_BMAP:
      case Teuchos::YAML::PROD_BVALUE_BSEQ: {
        swap(result_any, rhs.at(1));
        break;
      }
      case Teuchos::YAML::PROD_BMAP_FMAP: {
        map_item(result_any, rhs.at(0), rhs.at(4));
        break;
      }
      case Teuchos::YAML::PROD_BMAP_FSEQ: {
        TEUCHOS_ASSERT(rhs.at(4).type() == typeid(Array<Scalar>) ||
            rhs.at(4).type() == typeid(Array<Array<Scalar>>));
        int scalar_type = interpret_tag(rhs.at(3));
        map_item(result_any, rhs.at(0), rhs.at(4), scalar_type);
        break;
      }
      case Teuchos::YAML::PROD_BSEQ_FIRST: {
        seq_first_item(result_any, rhs.at(0));
        break;
      }
      case Teuchos::YAML::PROD_BSEQ_NEXT: {
        seq_next_item(result_any, rhs.at(0), rhs.at(1));
        break;
      }
      case Teuchos::YAML::PROD_BSEQ_SCALAR: {
        swap(result_any, rhs.at(3));
        Scalar& scalar = any_ref_cast<Scalar>(result_any);
        scalar.tag_type = interpret_tag(rhs.at(2));
        break;
      }
      case Teuchos::YAML::PROD_BSEQ_BSCALAR: {
        swap(result_any, rhs.at(2));
        break;
      }
      case Teuchos::YAML::PROD_BSEQ_BMAP:
      case Teuchos::YAML::PROD_BSEQ_BMAP_TRAIL:
      case Teuchos::YAML::PROD_BSEQ_FMAP: {
        throw ParserFail("Can't interpret a map inside a sequence as a Teuchos Parameter");
      }
      case Teuchos::YAML::PROD_BSEQ_BSEQ: {
        swap(result_any, rhs.at(3));
        break;
      }
      case Teuchos::YAML::PROD_BSEQ_BSEQ_TRAIL: {
        swap(result_any, rhs.at(4));
        break;
      }
      case Teuchos::YAML::PROD_BSEQ_FSEQ: {
        swap(result_any, rhs.at(3));
        break;
      }
      case Teuchos::YAML::PROD_FMAP: {
        swap(result_any, rhs.at(2));
        break;
      }
      case Teuchos::YAML::PROD_FMAP_EMPTY: {
        result_any = ParameterList();
        break;
      }
      case Teuchos::YAML::PROD_FSEQ: {
        swap(result_any, rhs.at(2));
        TEUCHOS_ASSERT(result_any.type() == typeid(Array<Scalar>) ||
            result_any.type() == typeid(Array<Array<Scalar>>));
        break;
      }
      case Teuchos::YAML::PROD_FSEQ_EMPTY: {
        result_any = Array<Scalar>();
        break;
      }
      case Teuchos::YAML::PROD_FSEQ_FIRST: {
        seq_first_item(result_any, rhs.at(0));
        break;
      }
      case Teuchos::YAML::PROD_FSEQ_NEXT: {
        seq_next_item(result_any, rhs.at(0), rhs.at(3));
        break;
      }
      case Teuchos::YAML::PROD_FSEQ_SCALAR: {
        swap(result_any, rhs.at(1));
        Scalar& scalar = any_ref_cast<Scalar>(result_any);
        scalar.tag_type = interpret_tag(rhs.at(0));
        break;
      }
      case Teuchos::YAML::PROD_FSEQ_FSEQ:
      case Teuchos::YAML::PROD_FSEQ_FMAP: {
        swap(result_any, rhs.at(1));
        break;
      }
      case Teuchos::YAML::PROD_SCALAR_QUOTED:
      case Teuchos::YAML::PROD_MAP_SCALAR_QUOTED: {
        swap(result_any, rhs.at(0));
        break;
      }
      case Teuchos::YAML::PROD_SCALAR_RAW:
      case Teuchos::YAML::PROD_MAP_SCALAR_RAW: {
        Scalar& scalar = make_any_ref<Scalar>(result_any);
        TEUCHOS_ASSERT(rhs.at(0).has_value());
        scalar.text = any_ref_cast<std::string>(rhs.at(0));
        scalar.text += any_ref_cast<std::string>(rhs.at(1));
        if (prod == Teuchos::YAML::PROD_MAP_SCALAR_RAW) {
          scalar.text += any_ref_cast<std::string>(rhs.at(2));
        }
        scalar.text = remove_trailing_whitespace(scalar.text);
        scalar.source = Scalar::RAW;
        scalar.tag_type = -1;
        break;
      }
      case Teuchos::YAML::PROD_SCALAR_HEAD_OTHER:
      case Teuchos::YAML::PROD_SCALAR_HEAD_DOT:
      case Teuchos::YAML::PROD_SCALAR_HEAD_DASH:
      case Teuchos::YAML::PROD_SCALAR_HEAD_DOT_DOT: {
        std::size_t offset;
        if (prod == Teuchos::YAML::PROD_SCALAR_HEAD_OTHER) offset = 0;
        else if (prod == Teuchos::YAML::PROD_SCALAR_HEAD_DOT_DOT) offset = 2;
        else offset = 1;
        char second = any_cast<char>(rhs.at(offset));
        std::string& result = make_any_ref<std::string>(result_any);
        if (prod == Teuchos::YAML::PROD_SCALAR_HEAD_DOT) result += '.';
        else if (prod == Teuchos::YAML::PROD_SCALAR_HEAD_DASH) result += '-';
        else if (prod == Teuchos::YAML::PROD_SCALAR_HEAD_DOT_DOT) result += "..";
        result += second;
        break;
      }
      case Teuchos::YAML::PROD_SCALAR_DQUOTED:
      case Teuchos::YAML::PROD_SCALAR_SQUOTED: {
        std::string& first = any_ref_cast<std::string>(rhs.at(1));
        std::string& rest = any_ref_cast<std::string>(rhs.at(2));
        Scalar& scalar = make_any_ref<Scalar>(result_any);
        scalar.text += first;
        scalar.text += rest;
        if (prod == Teuchos::YAML::PROD_SCALAR_DQUOTED) {
          scalar.source = Scalar::DQUOTED;
        } else if (prod == Teuchos::YAML::PROD_SCALAR_SQUOTED) {
          scalar.source = Scalar::SQUOTED;
        }
        scalar.tag_type = -1;
        break;
      }
      case Teuchos::YAML::PROD_MAP_SCALAR_ESCAPED_EMPTY: {
        result_any = std::string();
        break;
      }
      case Teuchos::YAML::PROD_MAP_SCALAR_ESCAPED_NEXT: {
        swap(result_any, rhs.at(0));
        std::string& str = any_ref_cast<std::string>(result_any);
        str += ',';
        str += any_ref_cast<std::string>(rhs.at(2));
        break;
      }
      case Teuchos::YAML::PROD_TAG: {
        swap(result_any, rhs.at(2));
        break;
      }
      case Teuchos::YAML::PROD_BSCALAR: {
        std::size_t parent_indent_level =
          this->symbol_indentation_stack.at(
              this->symbol_indentation_stack.size() - 5);
        std::string& header = any_ref_cast<std::string>(rhs.at(0));
        std::string& leading_empties_or_comments =
          any_ref_cast<std::string>(rhs.at(2));
        std::string& rest = any_ref_cast<std::string>(rhs.at(4));
        std::string& content = make_any_ref<std::string>(result_any);
        std::string comment;
        handle_block_scalar(
            parent_indent_level,
            header, leading_empties_or_comments, rest,
            content, comment);
        break;
      }
      case Teuchos::YAML::PROD_BSCALAR_FIRST: {
        swap(result_any, rhs.at(0));
        break;
      }
      // all these cases reduce to concatenating two strings
      case Teuchos::YAML::PROD_BSCALAR_NEXT:
      case Teuchos::YAML::PROD_BSCALAR_LINE:
      case Teuchos::YAML::PROD_DESCAPE_NEXT:
      case Teuchos::YAML::PROD_SESCAPE_NEXT: {
        swap(result_any, rhs.at(0));
        std::string& str = any_ref_cast<std::string>(result_any);
        str += any_ref_cast<std::string>(rhs.at(1));
        break;
      }
      case Teuchos::YAML::PROD_BSCALAR_INDENT: {
        swap(result_any, rhs.at(1));
        break;
      }
      case Teuchos::YAML::PROD_BSCALAR_HEADER_LITERAL:
      case Teuchos::YAML::PROD_BSCALAR_HEADER_FOLDED: {
        std::string& result = make_any_ref<std::string>(result_any);
        if (prod == Teuchos::YAML::PROD_BSCALAR_HEADER_LITERAL) {
          result += "|";
        } else {
          result += ">";
        }
        std::string& rest = any_ref_cast<std::string>(rhs.at(1));
        result += rest;
        break;
      }
      case Teuchos::YAML::PROD_DESCAPE: {
        std::string& str = make_any_ref<std::string>(result_any);
        std::string& rest = any_ref_cast<std::string>(rhs.at(2));
        str += any_cast<char>(rhs.at(1));
        str += rest;
        break;
      }
      case Teuchos::YAML::PROD_SESCAPE: {
        std::string& str = make_any_ref<std::string>(result_any);
        std::string& rest = any_ref_cast<std::string>(rhs.at(2));
        str += '\'';
        str += rest;
        break;
      }
      case Teuchos::YAML::PROD_OTHER_FIRST:
      case Teuchos::YAML::PROD_SPACE_PLUS_FIRST: {
        std::string& str = make_any_ref<std::string>(result_any);
        str.push_back(any_cast<char>(rhs.at(0)));
        break;
      }
      case Teuchos::YAML::PROD_SCALAR_TAIL_SPACE:
      case Teuchos::YAML::PROD_SCALAR_TAIL_OTHER:
      case Teuchos::YAML::PROD_DESCAPED_DQUOTED:
      case Teuchos::YAML::PROD_DQUOTED_COMMON:
      case Teuchos::YAML::PROD_SQUOTED_COMMON:
      case Teuchos::YAML::PROD_ANY_COMMON:
      case Teuchos::YAML::PROD_COMMON_SPACE:
      case Teuchos::YAML::PROD_COMMON_OTHER:
      case Teuchos::YAML::PROD_BSCALAR_HEAD_OTHER: {
        swap(result_any, rhs.at(0));
        break;
      }
      // all these cases reduce to appending a character
      case Teuchos::YAML::PROD_DQUOTED_NEXT:
      case Teuchos::YAML::PROD_SQUOTED_NEXT:
      case Teuchos::YAML::PROD_ANY_NEXT:
      case Teuchos::YAML::PROD_SCALAR_TAIL_NEXT:
      case Teuchos::YAML::PROD_SPACE_STAR_NEXT:
      case Teuchos::YAML::PROD_SPACE_PLUS_NEXT:
      case Teuchos::YAML::PROD_BSCALAR_HEAD_NEXT: {
        TEUCHOS_TEST_FOR_EXCEPTION(!rhs.at(0).has_value(), ParserFail,
            "leading characters in " << prod << ": any was empty\n");
        swap(result_any, rhs.at(0));
        std::string& str = any_ref_cast<std::string>(result_any);
        str += any_cast<char>(rhs.at(1));
        break;
      }
      case Teuchos::YAML::PROD_DQUOTED_EMPTY:
      case Teuchos::YAML::PROD_SQUOTED_EMPTY:
      case Teuchos::YAML::PROD_ANY_EMPTY:
      case Teuchos::YAML::PROD_DESCAPE_EMPTY:
      case Teuchos::YAML::PROD_SESCAPE_EMPTY:
      case Teuchos::YAML::PROD_SCALAR_TAIL_EMPTY:
      case Teuchos::YAML::PROD_SPACE_STAR_EMPTY:
      case Teuchos::YAML::PROD_BSCALAR_HEAD_EMPTY: {
        result_any = std::string();
        break;
      }
      case Teuchos::YAML::PROD_DESCAPED_DQUOT:
      case Teuchos::YAML::PROD_SQUOTED_DQUOT:
      case Teuchos::YAML::PROD_ANY_DQUOT: {
        result_any = '"';
        break;
      }
      case Teuchos::YAML::PROD_DESCAPED_SLASH:
      case Teuchos::YAML::PROD_SQUOTED_SLASH:
      case Teuchos::YAML::PROD_ANY_SLASH: {
        result_any = '\\';
        break;
      }
      case Teuchos::YAML::PROD_SCALAR_TAIL_SQUOT:
      case Teuchos::YAML::PROD_DQUOTED_SQUOT:
      case Teuchos::YAML::PROD_ANY_SQUOT: {
        result_any = '\'';
        break;
      }
      case Teuchos::YAML::PROD_COMMON_COLON: {
        result_any = ':';
        break;
      }
      case Teuchos::YAML::PROD_SCALAR_TAIL_DOT:
      case Teuchos::YAML::PROD_COMMON_DOT: {
        result_any = '.';
        break;
      }
      case Teuchos::YAML::PROD_SCALAR_TAIL_DASH:
      case Teuchos::YAML::PROD_COMMON_DASH:
      case Teuchos::YAML::PROD_BSCALAR_HEAD_DASH: {
        result_any = '-';
        break;
      }
      case Teuchos::YAML::PROD_COMMON_PIPE: {
        result_any = '|';
        break;
      }
      case Teuchos::YAML::PROD_COMMON_LSQUARE: {
        result_any = '[';
        break;
      }
      case Teuchos::YAML::PROD_COMMON_RSQUARE: {
        result_any = ']';
        break;
      }
      case Teuchos::YAML::PROD_COMMON_LCURLY: {
        result_any = '{';
        break;
      }
      case Teuchos::YAML::PROD_COMMON_RCURLY: {
        result_any = '}';
        break;
      }
      case Teuchos::YAML::PROD_COMMON_RANGLE: {
        result_any = '>';
        break;
      }
      case Teuchos::YAML::PROD_COMMON_COMMA: {
        result_any = ',';
        break;
      }
      case Teuchos::YAML::PROD_COMMON_PERCENT: {
        result_any = '%';
        break;
      }
      case Teuchos::YAML::PROD_COMMON_EXCL: {
        result_any = '!';
        break;
      }
    }
  }
  void map_first_item(any& result_any, any& first_item) {
    ParameterList& list = make_any_ref<ParameterList>(result_any);
    TEUCHOS_ASSERT(first_item.has_value());
    PLPair& pair = any_ref_cast<PLPair>(first_item);
    safe_set_entry(list, pair.key, pair.value);
  }
  void map_next_item(any& result_any, any& items, any& next_item) {
    using std::swap;
    swap(result_any, items);
    ParameterList& list = any_ref_cast<ParameterList>(result_any);
    PLPair& pair = any_ref_cast<PLPair>(next_item);
    safe_set_entry(list, pair.key, pair.value);
  }
  void map_item(any& result_any, any& key_any, any& value_any, int scalar_type = -1) {
    using std::swap;
    PLPair& result = make_any_ref<PLPair>(result_any);
    {
      std::string& key = any_ref_cast<Scalar>(key_any).text;
      swap(result.key, key);
    }
    resolve_map_value(value_any, scalar_type);
    if (value_any.type() == typeid(bool)) {
      bool value = any_cast<bool>(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(int)) {
      int value = any_cast<int>(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(long long)) {
      long long value = any_cast<long long>(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(double)) {
      double value = any_cast<double>(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(std::string)) {
      std::string& value = any_ref_cast<std::string >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(Array<int>)) {
      Array<int>& value = any_ref_cast<Array<int> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(Array<long long>)) {
      Array<long long>& value = any_ref_cast<Array<long long> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(Array<double>)) {
      Array<double>& value = any_ref_cast<Array<double> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(Array<std::string>)) {
      Array<std::string>& value = any_ref_cast<Array<std::string> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(TwoDArray<int>)) {
      TwoDArray<int>& value = any_ref_cast<TwoDArray<int> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(TwoDArray<long long>)) {
      TwoDArray<long long>& value = any_ref_cast<TwoDArray<long long> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(TwoDArray<double>)) {
      TwoDArray<double>& value = any_ref_cast<TwoDArray<double> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(TwoDArray<std::string>)) {
      TwoDArray<std::string>& value = any_ref_cast<TwoDArray<std::string> >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(ParameterList)) {
      ParameterList& value = any_ref_cast<ParameterList>(value_any);
      ParameterList& result_pl = result.value.setList();
      swap(result_pl, value);
      result_pl.setName(result.key);
    } else {
      std::string msg = "unexpected YAML map value type ";
      msg += value_any.type().name();
      msg += " for key \"";
      msg += result.key;
      msg += "\"\n";
      throw ParserFail(msg);
    }
  }
  void resolve_map_value(any& value_any, int scalar_type = -1) const {
    if (value_any.type() == typeid(Scalar)) {
      Scalar& scalar_value = any_ref_cast<Scalar>(value_any);
      if (scalar_type == -1) {
        scalar_type = scalar_value.infer_type();
      }
      if (scalar_type == Scalar::BOOL) {
        value_any = parse_as<bool>(scalar_value.text);
      } else if (scalar_type == Scalar::INT) {
        value_any = parse_as<int>(scalar_value.text);
      } else if (scalar_type == Scalar::LONG_LONG) {
        value_any = parse_as<long long>(scalar_value.text);
      } else if (scalar_type == Scalar::DOUBLE) {
        value_any = parse_as<double>(scalar_value.text);
      } else {
        value_any = scalar_value.text;
      }
    } else if (value_any.type() == typeid(Array<Scalar>)) {
      Array<Scalar>& scalars = any_ref_cast<Array<Scalar> >(value_any);
      if (scalar_type == -1) {
        if (scalars.size() == 0) {
          throw ParserFail("implicitly typed arrays can't be empty\n"
                           "(need to determine their element type)\n");
        }
        /* Teuchos::Array uses std::vector but doesn't account for std::vector<bool>,
           so it can't store bools */
        scalar_type = Scalar::INT;
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          scalar_type = std::min(scalar_type, scalars[i].infer_type());
        }
      }
      if (scalar_type == Scalar::INT) {
        Array<int> result(scalars.size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          result[i] = parse_as<int>(scalars[i].text);
        }
        value_any = result;
      } else if (scalar_type == Scalar::LONG_LONG) {
        Array<long long> result(scalars.size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          result[i] = parse_as<long long>(scalars[i].text);
        }
        value_any = result;
      } else if (scalar_type == Scalar::DOUBLE) {
        Array<double> result(scalars.size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          result[i] = parse_as<double>(scalars[i].text);
        }
        value_any = result;
      } else if (scalar_type == Scalar::STRING) {
        Array<std::string> result(scalars.size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          result[i] = scalars[i].text;
        }
        value_any = result;
      }
    } else if (value_any.type() == typeid(Array<Array<Scalar>>)) {
      Array<Array<Scalar>>& scalars = any_ref_cast<Array<Array<Scalar>> >(value_any);
      if (scalar_type == -1) {
        if (scalars.size() == 0) {
          throw ParserFail("implicitly typed 2D arrays can't be empty\n"
                           "(need to determine their element type)\n");
        }
        /* Teuchos::Array uses std::vector but doesn't account for std::vector<bool>,
           so it can't store bools */
        scalar_type = Scalar::INT;
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          if (scalars[0].size() == 0) {
            throw ParserFail("implicitly typed 2D arrays can't have empty rows\n"
                             "(need to determine their element type)\n");
          }
          if (scalars[i].size() != scalars[0].size()) {
            throw ParserFail("2D array: sub-arrays are different sizes");
          }
          for (Teuchos_Ordinal j = 0; j < scalars[i].size(); ++j) {
            int item_type = scalars[i][j].infer_type();
            scalar_type = std::min(scalar_type, item_type);
          }
        }
      }
      if (scalar_type == Scalar::INT) {
        TwoDArray<int> result(scalars.size(), scalars[0].size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          for (Teuchos_Ordinal j = 0; j < scalars[0].size(); ++j) {
            result(i, j) = parse_as<int>(scalars[i][j].text);
          }
        }
        value_any = result;
      } else if (scalar_type == Scalar::LONG_LONG) {
        TwoDArray<long long> result(scalars.size(), scalars[0].size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          for (Teuchos_Ordinal j = 0; j < scalars[0].size(); ++j) {
            result(i, j) = parse_as<long long>(scalars[i][j].text);
          }
        }
        value_any = result;
      } else if (scalar_type == Scalar::DOUBLE) {
        TwoDArray<double> result(scalars.size(), scalars[0].size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          for (Teuchos_Ordinal j = 0; j < scalars[0].size(); ++j) {
            result(i, j) = parse_as<double>(scalars[i][j].text);
          }
        }
        value_any = result;
      } else if (scalar_type == Scalar::STRING) {
        TwoDArray<std::string> result(scalars.size(), scalars[0].size());
        for (Teuchos_Ordinal i = 0; i < scalars.size(); ++i) {
          for (Teuchos_Ordinal j = 0; j < scalars[0].size(); ++j) {
            result(i, j) = scalars[i][j].text;
          }
        }
        value_any = result;
      }
    }
  }
  int interpret_tag(any& tag_any) {
    if (tag_any.type() != typeid(std::string)) return -1;
    std::string& text = any_ref_cast<std::string>(tag_any);
    if (text.find("bool") != std::string::npos) return Scalar::BOOL;
    else if (text.find("int") != std::string::npos) return Scalar::INT;
    else if (text.find("double") != std::string::npos) return Scalar::DOUBLE;
    else if (text.find("string") != std::string::npos) return Scalar::STRING;
    else {
      std::string msg = "Unable to parse type tag \"";
      msg += text;
      msg += "\"\n";
      throw ParserFail(msg);
    }
  }
  void seq_first_item(any& result_any, any& first_any) {
    using std::swap;
    if (first_any.type() == typeid(Scalar)) {
      Array<Scalar>& a = make_any_ref<Array<Scalar> >(result_any);
      Scalar& v = any_ref_cast<Scalar>(first_any);
      a.push_back(Scalar());
      swap(a.back(), v);
    } else if (first_any.type() == typeid(Array<Scalar>)) {
      Array<Array<Scalar>>& a = make_any_ref<Array<Array<Scalar>> >(result_any);
      Array<Scalar>& v = any_ref_cast<Array<Scalar> >(first_any);
      a.push_back(Array<Scalar>());
      swap(a.back(), v);
    } else {
      throw Teuchos::ParserFail(
          "bug in YAMLParameterList::Reader: unexpected type for first sequence item");
    }
  }
  void seq_next_item(any& result_any, any& items, any& next_item) {
    using std::swap;
    swap(result_any, items);
    if (result_any.type() == typeid(Array<Scalar>)) {
      Array<Scalar>& a = any_ref_cast<Array<Scalar> >(result_any);
      Scalar& val = any_ref_cast<Scalar>(next_item);
      a.push_back(Scalar());
      swap(a.back(), val);
    } else if (result_any.type() == typeid(Array<Array<Scalar>>)) {
      Array<Array<Scalar>>& a = any_ref_cast<Array<Array<Scalar>> >(result_any);
      Array<Scalar>& v = any_ref_cast<Array<Scalar> >(next_item);
      a.push_back(Array<Scalar>());
      swap(a.back(), v);
    } else {
      throw Teuchos::ParserFail(
          "bug in YAMLParameterList::Reader: unexpected type for next sequence item");
    }
  }
  /* block scalars are a super complicated mess, this function handles that mess */
  void handle_block_scalar(
      std::size_t parent_indent_level,
      std::string const& header,
      std::string const& leading_empties_or_comments,
      std::string const& rest,
      std::string& content,
      std::string& comment) {
    /* read the header, resulting in: block style, chomping indicator, and indentation indicator */
    char style;
    char chomping_indicator;
    std::size_t indentation_indicator = 0;
    style = header[0];
    std::stringstream ss(header.substr(1,std::string::npos));
    if (header.size() > 1 && my_isdigit(header[1])) {
      ss >> indentation_indicator;
      // indentation indicator is given as a relative number, but we need it in absolute terms
      indentation_indicator += parent_indent_level;
    }
    if (!(ss >> chomping_indicator)) chomping_indicator = '\0';
    /* get information about newlines, indentation level, and comment from
       the leading_empties_or_comments string */
    std::size_t first_newline = leading_empties_or_comments.find_first_of("\r\n");
    std::string newline;
    if (first_newline > 0 && leading_empties_or_comments[first_newline - 1] == '\r') {
      newline = "\r\n";
    } else {
      newline = "\n";
    }
    std::size_t keep_beg = first_newline + 1 - newline.size();
    if (leading_empties_or_comments[0] == '#') {
      comment = leading_empties_or_comments.substr(1, keep_beg);
    }
    // according to the YAML spec, a tab is content, not indentation
    std::size_t content_beg = leading_empties_or_comments.find_first_not_of("\r\n ");
    if (content_beg == std::string::npos) content_beg = leading_empties_or_comments.size();
    std::size_t newline_before_content = leading_empties_or_comments.rfind("\n", content_beg);
    std::size_t num_indent_spaces = (content_beg - newline_before_content) - 1;
    /* indentation indicator overrides the derived level of indentation, in case the
       user wants to keep some of that indentation as content */
    if (indentation_indicator > 0) {
      TEUCHOS_TEST_FOR_EXCEPTION(num_indent_spaces < indentation_indicator,
          Teuchos::ParserFail,
          "Indentation indicator " << indentation_indicator << " > leading spaces " << num_indent_spaces);
      num_indent_spaces = indentation_indicator;
    }
    /* prepend the content from the leading_empties_or_comments to the rest */
    content = leading_empties_or_comments.substr(keep_beg, std::string::npos);
    content += rest;
    /* per Trilinos issue #2090, there can be trailing comments after the block
       scalar which are less indented than it, but they will be included in the
       final NEWLINE token.
       this code removes all contiguous trailing lines which are less indented
       than the content.
     */
    while (true) {
      auto last_newline = content.find_last_of("\n", content.size() - 2);
      if (last_newline == std::string::npos) break;
      std::size_t num_spaces = 0;
      for (auto ispace = last_newline + 1;
           ispace < content.size() && content[ispace] == ' ';
           ++ispace) {
        ++num_spaces;
      }
      if (num_spaces >= num_indent_spaces) break;
      content.erase(content.begin() + last_newline + 1, content.end());
    }
    /* remove both indentation and newlines as dictated by header information */
    std::size_t unindent_pos = 0;
    while (true) {
      std::size_t next_newline = content.find_first_of("\n", unindent_pos);
      if (next_newline == std::string::npos) break;
      std::size_t start_cut = next_newline + 1;
      /* folding block scalars remove newlines */
      if (style == '>') start_cut -= newline.size();
      std::size_t end_cut = next_newline + 1;
      /* the actual amount of indentation in the content varies, start by
         marking it all for removal */
      while (end_cut < content.size() && content[end_cut] == ' ') {
        ++end_cut;
      }
      /* but don't remove more than the actual indent number */
      end_cut = std::min(next_newline + 1 + num_indent_spaces, end_cut);
      /* cut this (newline?)+indentation out of the content */
      content = content.substr(0, start_cut) +
        content.substr(end_cut, std::string::npos);
      unindent_pos = start_cut;
    }
    if (chomping_indicator != '+') {
      content = remove_trailing_whitespace_and_newlines(content);
      if (chomping_indicator != '-') content += newline;
    }
    if (style == '|') {
      // if not already, remove the leading newline
      content = content.substr(newline.size(), std::string::npos);
    }
  }
};

} // end namespace YAMLParameterList

/* Helper functions */

void updateParametersFromYamlFile(const std::string& yamlFileName,
                                  const Teuchos::Ptr<Teuchos::ParameterList>& paramList)
{
  //load the YAML file in as a new param list
  Teuchos::RCP<Teuchos::ParameterList> updated = YAMLParameterList::parseYamlFile(yamlFileName);
  if (paramList->name() == "ANONYMOUS") {
    paramList->setName(updated->name());
  }
  //now update the original list (overwriting values with same key)
  paramList->setParameters(*updated);
}

void updateParametersFromYamlCString(const char* const data,
                                     const Teuchos::Ptr<Teuchos::ParameterList>& paramList,
                                     bool overwrite)
{
  Teuchos::RCP<Teuchos::ParameterList> updated = YAMLParameterList::parseYamlText(data, "CString");
  if(overwrite)
  {
    if (paramList->name() == "ANONYMOUS") {
      paramList->setName(updated->name());
    }
    paramList->setParameters(*updated);
  }
  else
  {
    paramList->setParametersNotAlreadySet(*updated);
  }
}

void updateParametersFromYamlString(const std::string& yamlData,
                                  const Teuchos::Ptr<Teuchos::ParameterList>& paramList,
                                  bool overwrite,
                                  const std::string& name)
{
  Teuchos::RCP<Teuchos::ParameterList> updated = YAMLParameterList::parseYamlText(yamlData, name);
  if(overwrite)
  {
    if (paramList->name() == "ANONYMOUS") {
      paramList->setName(updated->name());
    }
    paramList->setParameters(*updated);
  }
  else
  {
    paramList->setParametersNotAlreadySet(*updated);
  }
}

Teuchos::RCP<Teuchos::ParameterList> getParametersFromYamlFile(const std::string& yamlFileName)
{
  return YAMLParameterList::parseYamlFile(yamlFileName);
}

Teuchos::RCP<Teuchos::ParameterList> getParametersFromYamlString(const std::string& yamlStr)
{
  std::stringstream ss(yamlStr);
  return YAMLParameterList::parseYamlStream(ss);
}

void writeParameterListToYamlOStream(
  const ParameterList &paramList,
  std::ostream &yamlOut
  )
{
  YAMLParameterList::writeYamlStream(yamlOut, paramList);
}

void writeParameterListToYamlFile(
  const ParameterList &paramList,
  const std::string &yamlFileName
  )
{
  YAMLParameterList::writeYamlFile(yamlFileName, paramList);
}

std::string convertXmlToYaml(const std::string& xmlFileName)
{
  //load the parameter list from xml
  Teuchos::RCP<Teuchos::ParameterList> toConvert = Teuchos::getParametersFromXmlFile(xmlFileName);
  //replace the file extension ".xml" with ".yaml", or append it if there was no extension
  std::string yamlFileName;
  if(xmlFileName.find(".xml") == std::string::npos)
  {
    yamlFileName = xmlFileName + ".yaml";
  }
  else
  {
    yamlFileName = xmlFileName.substr(0, xmlFileName.length() - 4) + ".yaml";
  }
  YAMLParameterList::writeYamlFile(yamlFileName, *toConvert);
  return yamlFileName;
}

void convertXmlToYaml(const std::string& xmlFileName, const std::string& yamlFileName)
{
  Teuchos::RCP<Teuchos::ParameterList> toConvert = Teuchos::getParametersFromXmlFile(xmlFileName);
  YAMLParameterList::writeYamlFile(yamlFileName, *toConvert);
}

void convertXmlToYaml(std::istream& xmlStream, std::ostream& yamlStream)
{
  //read xmlStream into a string until EOF
  std::istreambuf_iterator<char> begin(xmlStream);
  std::istreambuf_iterator<char> end;
  std::string xmlString(begin, end);
  Teuchos::RCP<Teuchos::ParameterList> toConvert = Teuchos::getParametersFromXmlString(xmlString);
  YAMLParameterList::writeYamlStream(yamlStream, *toConvert);
}

namespace YAMLParameterList
{

Teuchos::RCP<Teuchos::ParameterList> parseYamlText(const std::string& text, const std::string& name)
{
#ifdef HAVE_TEUCHOSPARAMETERLIST_YAMLCPP
  auto yaml_input = ::YAML::LoadAll(text); // std::vector<::YAML::Node>
  return readParams(yaml_input);
#else
  any result;
  Teuchos::YAMLParameterList::Reader reader;
  reader.read_string(result, text, name);
  ParameterList& pl = any_ref_cast<ParameterList>(result);
  return Teuchos::rcp(new ParameterList(pl));
#endif // HAVE_TEUCHOSPARAMETERLIST_YAMLCPP
}

Teuchos::RCP<Teuchos::ParameterList> parseYamlFile(const std::string& yamlFile)
{
#ifdef HAVE_TEUCHOSPARAMETERLIST_YAMLCPP
  auto yaml_input = ::YAML::LoadAllFromFile(yamlFile);
  return readParams(yaml_input);
#else
  any result;
  Teuchos::YAMLParameterList::Reader reader;
  reader.read_file(result, yamlFile);
  ParameterList& pl = any_ref_cast<ParameterList>(result);
  return Teuchos::rcp(new ParameterList(pl));
#endif // HAVE_TEUCHOSPARAMETERLIST_YAMLCPP
}

Teuchos::RCP<Teuchos::ParameterList> parseYamlStream(std::istream& yaml)
{
#ifdef HAVE_TEUCHOSPARAMETERLIST_YAMLCPP
  auto yaml_input = ::YAML::LoadAll(yaml);
  return readParams(yaml_input);
#else
  any result;
  Teuchos::YAMLParameterList::Reader reader;
  reader.read_stream(result, yaml, "parseYamlStream");
  ParameterList& pl = any_ref_cast<ParameterList>(result);
  return Teuchos::rcp(new ParameterList(pl));
#endif // HAVE_TEUCHOSPARAMETERLIST_YAMLCPP
}

// The following three functions (readParams, processMapNode, and processKeyValueNode)
// were previously removed from Trilinos in PR 1779 (Teuchos: use Parser, not yaml-cpp, to read YAML PL).

#ifdef HAVE_TEUCHOSPARAMETERLIST_YAMLCPP

Teuchos::RCP<Teuchos::ParameterList> readParams(std::vector<::YAML::Node>& lists)
{
  Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList); //pl is the root ParameterList to be returned
  // If there is exactly one element in "lists", assume it is the anonymous top-level parameter list
  // If there are more than one, place them all in the anonymous top-level list
  for(size_t i = 0; i < lists.size(); i++)
  {
    processMapNode(lists[i], *pl, true);
  }
  return pl;
}

void processMapNode(const ::YAML::Node& node, Teuchos::ParameterList& parent, bool topLevel)
{
  if (node.Type() != ::YAML::NodeType::Map)
  {
    throw YamlStructureError("All top-level elements of the YAML file must be maps.");
  }
  if (topLevel)
  {
    parent.setName(node.begin()->first.as<std::string>());
    processMapNode(node.begin()->second, parent);
  }
  else
  {
    for (::YAML::const_iterator i = node.begin(); i != node.end(); i++)
    {
      // make sure the key type is a string
      if(i->first.Type() != ::YAML::NodeType::Scalar)
      {
        throw YamlKeyError("Keys must be YAML scalars (int, double, or string)");
      }
      // if this conversion fails and throws for any reason (shouldn't), let the caller handle it
      const std::string key = quoted_as<std::string>(i->first);
      processKeyValueNode(key, i->second, parent, topLevel);
    }
  }
}

void processKeyValueNode(const std::string& key, const ::YAML::Node& node, Teuchos::ParameterList& parent, bool topLevel)
{
  // node (value) type can be a map (for nested param lists),
  // a scalar (int, double, string), or a sequence of doubles (vector<double>)
  if(node.Type() == ::YAML::NodeType::Scalar)
  {
    try
    {
      safe_set_entry<int>(parent, key, quoted_as<int>(node));
    }
    catch(...)
    {
      try
      {
        safe_set_entry<long long>(parent, key, quoted_as<long long>(node));
      }
      catch(...)
      {
        try
        {
          safe_set_entry<double>(parent, key, quoted_as<double>(node));
        }
        catch(...)
        {
          try
          {
            bool raw_bool = quoted_as<bool>(node);

            /* yaml-cpp parses ON/OFF as a bool, but the in-house parser does not.
               To preserve backwards compatibility, make sure the string passes
               the in-house parser's is_parseable_as<bool> function (which protects
               against the ON/OFF case).
               Otherwise, a failure is observed in YAML_ConvertFromXML unit test.*/

            std::string raw_string = quoted_as<std::string>(node);
            if (is_parseable_as<bool>(raw_string))
            {
              safe_set_entry<bool>(parent, key, raw_bool);
            }
            else
            {
              safe_set_entry<std::string>(parent, key, raw_string);
            }
          }
          catch(...)
          {
            safe_set_entry<std::string>(parent, key, quoted_as<std::string>(node));
          }
        }
      }
    }
  }
  else if(node.Type() == ::YAML::NodeType::Map)
  {
    if(topLevel)
    {
      processMapNode(node, parent);
    }
    else
    {
      Teuchos::ParameterList& sublist = parent.sublist(key);
      processMapNode(node, sublist);
    }
  }
  else if(node.Type() == ::YAML::NodeType::Sequence)
  {
    int ndim = getYamlArrayDim(node);
    if (ndim == 1)
    {
      ::YAML::Node const& first_value = *(node.begin());
      try
      {
        quoted_as<int>(first_value);
        safe_set_entry<Teuchos::Array<int>>(parent, key, getYamlArray<int>(node));
      }
      catch(...)
      {
        try
        {
          quoted_as<double>(first_value);
          safe_set_entry<Teuchos::Array<double>>(parent, key, getYamlArray<double>(node));
        }
        catch(...)
        {
          try
          {
            quoted_as<std::string>(first_value);
            safe_set_entry<Teuchos::Array<std::string>>(parent, key, getYamlArray<std::string>(node));
          }
          catch(...)
          {
            throw YamlSequenceError(std::string("Array \"") + key + "\" must contain int, double, bool or string");
          }
        }
      }
    }
    else if (ndim == 2)
    {
      bool is_ragged = checkYamlTwoDArrayIsRagged(node);
      ::YAML::Node const& first_value = *(node.begin()->begin());
      try
      {
        quoted_as<int>(first_value);
        using arr_t = Teuchos::Array<Teuchos::Array<int>>;
        if (is_ragged) {
          safe_set_entry<arr_t>(parent, key, getYaml2DRaggedArray<arr_t, int>(node, ndim, key));
        } else {
          safe_set_entry<Teuchos::TwoDArray<int>>(parent, key, getYamlTwoDArray<int>(node));
        }
      }
      catch(...)
      {
        try
        {
          quoted_as<double>(first_value);
          using arr_t = Teuchos::Array<Teuchos::Array<double>>;
          if (is_ragged) {
            safe_set_entry<arr_t>(parent, key, getYaml2DRaggedArray<arr_t, double>(node, ndim, key));
          } else {
            safe_set_entry<Teuchos::TwoDArray<double>>(parent, key, getYamlTwoDArray<double>(node));
          }
        }
        catch(...)
        {
          try
          {
            quoted_as<std::string>(first_value);
            using arr_t = Teuchos::Array<Teuchos::Array<std::string>>;
            if (is_ragged) {
              safe_set_entry<arr_t>(parent, key, getYaml2DRaggedArray<arr_t, std::string>(node, ndim, key));
            } else {
              safe_set_entry<Teuchos::TwoDArray<std::string>>(parent, key, getYamlTwoDArray<std::string>(node));
            }
          }
          catch(...)
          {
            throw YamlSequenceError(std::string("TwoDArray \"") + key + "\" must contain int, double, bool or string");
          }
        }
      }
    }
    else if (ndim == 3)
    {
      ::YAML::Node const& first_value = *(node.begin()->begin()->begin());
      try
      {
        quoted_as<int>(first_value);
        using arr_t = Teuchos::Array<Teuchos::Array<Teuchos::Array<int>>>;
        safe_set_entry<arr_t>(parent, key, getYaml3DArray<arr_t, int>(node, ndim, key));
      }
      catch(...)
      {
        try
        {
          quoted_as<double>(first_value);
          using arr_t = Teuchos::Array<Teuchos::Array<Teuchos::Array<double>>>;
          safe_set_entry<arr_t>(parent, key, getYaml3DArray<arr_t, double>(node, ndim, key));

        }
        catch(...)
        {
          try
          {
            quoted_as<std::string>(first_value);
            using arr_t = Teuchos::Array<Teuchos::Array<Teuchos::Array<std::string>>>;
            safe_set_entry<arr_t>(parent, key, getYaml3DArray<arr_t, std::string>(node, ndim, key));

          }
          catch(...)
          {
            throw YamlSequenceError(std::string("3DArray \"") + key + "\" must contain int, double, bool or string");
          }
        }
      }
    }
  }
  else if(node.Type() == ::YAML::NodeType::Null)
  {
    // treat NULL as empty sublist (not an error)
    parent.sublist(key);
  }
  else
  {
    // Undefined
    throw YamlUndefinedNodeError("Value type in a key-value pair must be one of: int, double, string, array, sublist.");
  }
}

#endif // HAVE_TEUCHOSPARAMETERLIST_YAMLCPP

void writeYamlStream(std::ostream& yaml, const Teuchos::ParameterList& pl)
{
  // warn the user if floats/doubles with integer values will be printed incorrectly
  std::ios_base::fmtflags flags = yaml.flags();
  // make temporary stringstream to test flags
  std::ostringstream testStream;
  testStream.flags(flags);
  double testVal = 1;
  testStream << testVal;
  bool popFlags = false;
  if(testStream.str() == "1")
  {
    // must add showpoint to flags while writing yaml
    // this will always disambiguate (double) n and n, even if stream precision is 0
    // prints as "n.0000" where the number of trailing zeros is the stream precision
    // note: in YAML, "5." is a double but not an int
    std::cout << "Warning: yaml stream format flags would confuse double with integer value with int.\n";
    std::cout << "Setting std::ios::showpoint on the stream to fix this (will restore flags when done)\n";
    std::ios_base::fmtflags flagsCopy = flags;
    flagsCopy |= std::ios::showpoint;
    popFlags = true;
  }
  yaml << "%YAML 1.1\n---\n";
  yaml << pl.name() << ':';
  if(pl.numParams() == 0)
  {
    yaml << " { }\n";
  }
  else
  {
    writeParameterList(pl, yaml, 2);
  }
  yaml << "...\n";
  // restore flags
  if(popFlags)
  {
    yaml.flags(flags);
  }
}

void writeYamlFile(const std::string& yamlFile, const Teuchos::ParameterList& pl)
{
  std::ofstream yaml(yamlFile.c_str());
  /* set default floating-point style:
     1. 17 decimal places to ensure the value remains the same
     2. scientific: this prevents floating-point values that happen
        to be integers from being printed as integers, because YAML
        will then read that value back typed as an integer.
   */
  yaml << std::scientific << std::setprecision(17);
  writeYamlStream(yaml, pl);
}

void writeParameterList(const Teuchos::ParameterList& pl, std::ostream& yaml, int indentLevel)
{
  if(pl.begin() == pl.end())
  {
    yaml << "{ }\n";
  }
  else
  {
    yaml << '\n';
    for(PLIter it = pl.begin(); it != pl.end(); it++)
    {
      writeParameter(pl.name(it), pl.entry(it), yaml, indentLevel);
    }
  }
}

template <typename T>
struct YamlWrite {
  static void write(T const& x, std::ostream& stream) {
    stream << x;
  }
};

template <>
struct YamlWrite<double> {
  static void write(double const& x, std::ostream& stream) {
    generalWriteDouble(x, stream);
  }
};

template <>
struct YamlWrite<std::string> {
  static void write(std::string const& x, std::ostream& stream) {
    generalWriteString(x, stream);
  }
};

template <typename T>
void writeYamlTwoDArray(Teuchos::TwoDArray<T> const& arr, std::ostream& stream)
{
  typename Teuchos::TwoDArray<T>::size_type i, j;
  stream << '[';
  for (i = 0; i < arr.getNumRows(); ++i)
  {
    if (i) stream << ", ";
    stream << '[';
    for (j = 0; j < arr.getNumCols(); ++j)
    {
      if (j) stream << ", ";
      YamlWrite<T>::write(arr(i, j), stream);
    }
    stream << ']';
  }
  stream << ']';
}

void writeParameter(const std::string& paramName, const Teuchos::ParameterEntry& entry, std::ostream& yaml, int indentLevel)
{
  for(int i = 0; i < indentLevel; i++)
  {
    yaml << ' ';
  }
  generalWriteString(paramName, yaml);
  yaml << ": ";
  if(entry.isList())
  {
    writeParameterList(Teuchos::getValue<Teuchos::ParameterList>(entry), yaml, indentLevel + 2);
    return;
  }
  else if(entry.isArray())
  {
    yaml << '[';
    if(entry.isType<Teuchos::Array<int> >())
    {
      Teuchos::Array<int>& arr = Teuchos::getValue<Teuchos::Array<int> >(entry);
      for(int i = 0; i < arr.size(); i++)
      {
        yaml << arr[i];
        if(i != arr.size() - 1)
          yaml << ", ";
      }
    }
    if(entry.isType<Teuchos::Array<long long> >())
    {
      Teuchos::Array<long long>& arr = Teuchos::getValue<Teuchos::Array<long long> >(entry);
      for(int i = 0; i < arr.size(); i++)
      {
        yaml << arr[i];
        if(i != arr.size() - 1)
          yaml << ", ";
      }
    }
    else if(entry.isType<Teuchos::Array<double> >())
    {
      Teuchos::Array<double>& arr = Teuchos::getValue<Teuchos::Array<double> >(entry);
      for(int i = 0; i < arr.size(); i++)
      {
        generalWriteDouble(arr[i], yaml);
        if(i != arr.size() - 1)
          yaml << ", ";
      }
    }
    else if(entry.isType<Teuchos::Array<std::string> >())
    {
      Teuchos::Array<std::string>& arr = Teuchos::getValue<Teuchos::Array<std::string> >(entry);
      for(int i = 0; i < arr.size(); i++)
      {
        generalWriteString(arr[i], yaml);
        if(i != arr.size() - 1)
          yaml << ", ";
      }
    }
    yaml << ']';
  }
  else if(entry.isTwoDArray())
  {
    if(entry.isType<Teuchos::TwoDArray<int> >())
    {
      writeYamlTwoDArray<int>(
          Teuchos::getValue<Teuchos::TwoDArray<int> >(entry), yaml);
    }
    if(entry.isType<Teuchos::TwoDArray<long long> >())
    {
      writeYamlTwoDArray<long long>(
          Teuchos::getValue<Teuchos::TwoDArray<long long> >(entry), yaml);
    }
    else if(entry.isType<Teuchos::TwoDArray<double> >())
    {
      writeYamlTwoDArray<double>(
          Teuchos::getValue<Teuchos::TwoDArray<double> >(entry), yaml);
    }
    else if(entry.isType<Teuchos::TwoDArray<std::string> >())
    {
      writeYamlTwoDArray<std::string>(
          Teuchos::getValue<Teuchos::TwoDArray<std::string> >(entry), yaml);
    }
  }
  else if(entry.isType<int>())
  {
    yaml << Teuchos::getValue<int>(entry);
  }
  else if(entry.isType<long long>())
  {
    yaml << Teuchos::getValue<long long>(entry);
  }
  else if(entry.isType<double>())
  {
    generalWriteDouble(Teuchos::getValue<double>(entry), yaml);
  }
  else if(entry.isType<std::string>())
  {
    std::string& str = Teuchos::getValue<std::string>(entry);
    if(strchr(str.c_str(), '\n'))
    {
      yaml << "|";
      // if the content has leading spaces, automatic indentation
      // detection would fail, in which case we must emit
      // an indentation indicator
      std::size_t first_non_newline_pos = str.find_first_not_of("\r\n");
      if (first_non_newline_pos != std::string::npos &&
          str[first_non_newline_pos] == ' ') {
        yaml << "2";
      }
      if (str[str.size() - 1] != '\n') yaml << "-";
      yaml << "\n";
      //for each line, apply indent then print the line verbatim
      size_t index = 0;
      while(true)
      {
        size_t next = str.find('\n', index);
        for(int i = 0; i < indentLevel + 2; i++)
        {
          yaml << ' ';
        }
        if(next == std::string::npos)
        {
          yaml << str.substr(index, std::string::npos);
          break;
        }
        else
        {
          yaml << str.substr(index, next - index) << '\n';
        }
        index = next + 1;
      }
    }
    else
    {
      generalWriteString(str, yaml);
    }
  }
  else if(entry.isType<bool>())
  {
    yaml << (Teuchos::getValue<bool>(entry) ? "true" : "false");
  }
  yaml << '\n';
}

void generalWriteString(const std::string& str, std::ostream& yaml)
{
  // default to single quoting
  if(stringNeedsQuotes(str))
  {
    yaml << '\'';
    for (std::size_t i = 0; i < str.size(); ++i) {
      if (str[i] == '\'') yaml << "''";
      else yaml << str[i];
    }
    yaml << '\'';
  }
  else
  {
    yaml << str;
  }
}

void generalWriteDouble(double d, std::ostream& yaml)
{
  yaml << d;
}

static bool containsSpecialCharacters(std::string const& s) {
  char const* const control_chars = ":'{}[],&*#?|<>=!%@\\";
  return s.find_first_of(control_chars) != std::string::npos;
}

bool stringNeedsQuotes(const std::string& s)
{
  return s.empty() ||
         containsSpecialCharacters(s) ||
         is_parseable_as<bool>(s) ||
         is_parseable_as<int>(s) ||
         is_parseable_as<long long>(s) ||
         is_parseable_as<double>(s);
}

} //namespace YAMLParameterList

} //namespace Teuchos
