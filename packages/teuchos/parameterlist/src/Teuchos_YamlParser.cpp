// @HEADER
//
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Dan Ibanez        (daibane@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef TEUCHOS_YAMLPARSER_DEF_H_
#define TEUCHOS_YAMLPARSER_DEF_H_

#include <iostream>
#include <iomanip>
#include <ios>
#include <sstream>
#include <cctype>

#include "Teuchos_YamlParser_decl.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_YamlParameterListCoreHelpers.hpp"
#include "Teuchos_TwoDArray.hpp"

#include "Teuchos_Reader.hpp"
#include "Teuchos_YAML.hpp"

namespace Teuchos {

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

bool is_parseable_as_bool(std::string const& text) {
  std::string lower;
  for (std::size_t i = 0; i < text.size(); ++i) {
    lower.push_back(my_tolower(text[i]));
  }
  return lower == "true" || lower == "yes" ||
         lower == "false" || lower == "no";
}

bool parse_as_bool(std::string const& text) {
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
    INT    = 2,
    BOOL   = 3
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
    if (is_parseable_as_bool(text)) {
      return BOOL;
    }
    if (is_parseable_as<int>(text)) {
      return INT;
    }
    if (is_parseable_as<double>(text)) {
      return DOUBLE;
    }
    return STRING;
  }
};

bool operator==(PLPair const&, PLPair const&) { return false; }
std::ostream& operator<<(std::ostream& os, PLPair const&) { return os; }

bool operator==(Scalar const&, Scalar const&) { return false; }
std::ostream& operator<<(std::ostream& os, Scalar const&) { return os; }

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
        TEUCHOS_ASSERT(!rhs.at(offset).empty());
        swap(result_any, rhs.at(offset));
        TEUCHOS_ASSERT(result_any.type() == typeid(ParameterList));
        break;
      }
      case Teuchos::YAML::PROD_TOP_BMAP: {
        TEUCHOS_ASSERT(!rhs.at(0).empty());
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
          TEUCHOS_TEST_FOR_EXCEPTION(!rhs.at(0).empty(), ParserFail,
              "Can't specify multiple top-level ParameterLists in one YAML file!\n");
          swap(result_any, rhs.at(1));
        } else {
          swap(result_any, rhs.at(0));
        }
        break;
      }
      case Teuchos::YAML::PROD_BMAP_FIRST: {
        TEUCHOS_ASSERT(rhs.at(0).type() == typeid(PLPair));
        map_first_item(result_any, rhs.at(0));
        TEUCHOS_ASSERT(result_any.type() == typeid(ParameterList));
        break;
      }
      case Teuchos::YAML::PROD_BMAP_NEXT: {
        map_next_item(result_any, rhs.at(0), rhs.at(1));
        break;
      }
      case Teuchos::YAML::PROD_BMAP_SCALAR: {
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
            rhs.at(4).type() == typeid(Array<Array<Scalar> >));
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
            result_any.type() == typeid(Array<Array<Scalar> >));
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
        TEUCHOS_ASSERT(!rhs.at(0).empty());
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
        TEUCHOS_TEST_FOR_EXCEPTION(rhs.at(0).empty(), ParserFail,
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
    TEUCHOS_ASSERT(!first_item.empty());
    PLPair& pair = any_ref_cast<PLPair>(first_item);
    list.setEntry(pair.key, pair.value);
  }
  void map_next_item(any& result_any, any& items, any& next_item) {
    using std::swap;
    swap(result_any, items);
    ParameterList& list = any_ref_cast<ParameterList>(result_any);
    PLPair& pair = any_ref_cast<PLPair>(next_item);
    list.setEntry(pair.key, pair.value);
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
    } else if (value_any.type() == typeid(double)) {
      double value = any_cast<double>(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(std::string)) {
      std::string& value = any_ref_cast<std::string >(value_any);
      result.value = ParameterEntry(value);
    } else if (value_any.type() == typeid(Array<int>)) {
      Array<int>& value = any_ref_cast<Array<int> >(value_any);
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
        value_any = parse_as_bool(scalar_value.text);
      } else if (scalar_type == Scalar::INT) {
        value_any = parse_as<int>(scalar_value.text);
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
    } else if (value_any.type() == typeid(Array<Array<Scalar> >)) {
      Array<Array<Scalar> >& scalars = any_ref_cast<Array<Array<Scalar> > >(value_any);
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
      Array<Array<Scalar> >& a = make_any_ref<Array<Array<Scalar> > >(result_any);
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
    } else if (result_any.type() == typeid(Array<Array<Scalar> >)) {
      Array<Array<Scalar> >& a = any_ref_cast<Array<Array<Scalar> > >(result_any);
      Array<Scalar>& v = any_ref_cast<Array<Scalar> >(next_item);
      a.push_back(Array<Scalar>());
      swap(a.back(), v);
    } else {
      throw Teuchos::ParserFail(
          "bug in YAMLParameterList::Reader: unexpected type for next sequence item");
    }
  }
  void handle_block_scalar(
      std::size_t parent_indent_level,
      std::string const& header,
      std::string const& leading_empties_or_comments,
      std::string const& rest,
      std::string& content,
      std::string& comment) {
    char style;
    char chomping_indicator;
    std::size_t indentation_indicator = 0;
    style = header[0];
    std::stringstream ss(header.substr(1,std::string::npos));
    if (header.size() > 1 && my_isdigit(header[1])) {
      ss >> indentation_indicator;
      indentation_indicator += parent_indent_level;
    }
    if (!(ss >> chomping_indicator)) chomping_indicator = '\0';
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
    if (indentation_indicator > 0) {
      TEUCHOS_TEST_FOR_EXCEPTION(num_indent_spaces < indentation_indicator,
          Teuchos::ParserFail,
          "Indentation indicator " << indentation_indicator << " > leading spaces " << num_indent_spaces);
      num_indent_spaces = indentation_indicator;
    }
    content = leading_empties_or_comments.substr(keep_beg, std::string::npos);
    content += rest;
    std::size_t unindent_pos = 0;
    while (true) {
      std::size_t next_newline = content.find_first_of("\n", unindent_pos);
      if (next_newline == std::string::npos) break;
      std::size_t start_cut = next_newline + 1;
      if (style == '>') start_cut -= newline.size();
      std::size_t end_cut = next_newline + 1;
      // this scanning and min are needed for trailing lines that are less indented
      while (end_cut < content.size() && content[end_cut] == ' ')
        ++end_cut;
      end_cut = std::min(end_cut, next_newline + 1 + num_indent_spaces);
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
  Teuchos::YAMLParameterList::Reader reader;
  any result;
  reader.read_string(result, text, name);
  ParameterList& pl = any_ref_cast<ParameterList>(result);
  return Teuchos::rcp(new ParameterList(pl));
}

Teuchos::RCP<Teuchos::ParameterList> parseYamlFile(const std::string& yamlFile)
{
  Teuchos::YAMLParameterList::Reader reader;
  any result;
  reader.read_file(result, yamlFile);
  ParameterList& pl = any_ref_cast<ParameterList>(result);
  return Teuchos::rcp(new ParameterList(pl));
}

Teuchos::RCP<Teuchos::ParameterList> parseYamlStream(std::istream& yaml)
{
  Teuchos::YAMLParameterList::Reader reader;
  any result;
  reader.read_stream(result, yaml, "parseYamlStream");
  ParameterList& pl = any_ref_cast<ParameterList>(result);
  return Teuchos::rcp(new ParameterList(pl));
}

void writeYamlStream(std::ostream& yaml, const Teuchos::ParameterList& pl)
{
  //warn the user if floats/doubles with integer values will be printed incorrectly
  std::ios_base::fmtflags flags = yaml.flags();
  //make temporary stringstream to test flags
  std::ostringstream testStream;
  testStream.flags(flags);
  double testVal = 1;
  testStream << testVal;
  bool popFlags = false;
  if(testStream.str() == "1")
  {
    //must add showpoint to flags while writing yaml
    //this will always disambiguate (double) n and n, even if stream precision is 0
    //prints as "n.0000" where the number of trailing zeros is the stream precision
    //note: in YAML, "5." is a double but not an int
    std::cout << "Warning: yaml stream format flags would confuse double with integer value with int.\n";
    std::cout << "Setting std::ios::showpoint on the stream to fix this (will restore flags when done)\n";
    std::ios_base::fmtflags flagsCopy = flags;
    flagsCopy |= std::ios::showpoint;
    popFlags = true;
  }
  yaml << "%YAML 1.1\n---\n";
  yaml << "ANONYMOUS:";         //original top-level list name is not stored by ParameterList
  if(pl.numParams() == 0)
  {
    yaml << " { }\n";
  }
  else
  {
    writeParameterList(pl, yaml, 2);
  }
  yaml << "...\n";
  //restore flags
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
         is_parseable_as_bool(s) ||
         is_parseable_as<int>(s) ||
         is_parseable_as<double>(s);
}

} //namespace YAMLParameterList

} //namespace Teuchos

#endif
