// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_YAML.hpp"

#include <iostream>

namespace Teuchos {
namespace YAML {

namespace {

Language::Productions make_productions() {
  Language::Productions prods;
  prods.resize(NPRODS);
  prods[PROD_DOC]("doc") >> "top_items";
  prods[PROD_DOC2]("doc") >> "NEWLINE", "top_items";
  prods[PROD_TOP_FIRST]("top_items") >> "top_item";
  prods[PROD_TOP_NEXT]("top_items") >> "top_items", "top_item";
  prods[PROD_TOP_DIRECT]("top_item") >> "%", "any*", "NEWLINE";
  prods[PROD_TOP_BEGIN]("top_item") >> "-", "-", "-", "NEWLINE";
  prods[PROD_TOP_END]("top_item") >> ".", ".", ".", "NEWLINE";
  prods[PROD_TOP_BMAP]("top_item") >> "bmap_item";
  prods[PROD_BMAP_FIRST]("bmap_items") >> "bmap_item";
  prods[PROD_BMAP_NEXT]("bmap_items") >> "bmap_items", "bmap_item";
  prods[PROD_BMAP_SCALAR]("bmap_item") >> "scalar", ":", "WS*", "tag?", "map_scalar", "NEWLINE";
  prods[PROD_BMAP_BSCALAR]("bmap_item") >> "scalar", ":", "WS*", "bscalar";
  prods[PROD_BMAP_BVALUE]("bmap_item") >> "scalar", ":", "WS*", "NEWLINE", "bvalue";
  prods[PROD_BVALUE_EMPTY]("bvalue");
  prods[PROD_BVALUE_BMAP]("bvalue") >> "INDENT", "bmap_items", "DEDENT";
  /* TODO: allow a tag in this */
  prods[PROD_BVALUE_BSEQ]("bvalue") >> "INDENT", "bseq_items", "DEDENT";
  prods[PROD_BMAP_FMAP]("bmap_item") >> "scalar", ":", "WS*", "tag?", "fmap", "NEWLINE";
  prods[PROD_BMAP_FSEQ]("bmap_item") >> "scalar", ":", "WS*", "tag?", "fseq", "NEWLINE";
  prods[PROD_BSEQ_FIRST]("bseq_items") >> "bseq_item";
  prods[PROD_BSEQ_NEXT]("bseq_items") >> "bseq_items", "bseq_item";
  prods[PROD_BSEQ_SCALAR]("bseq_item") >> "-", "WS+", "tag?", "scalar", "NEWLINE";
  prods[PROD_BSEQ_BSCALAR]("bseq_item") >> "-", "WS+", "bscalar";
  prods[PROD_BSEQ_BMAP]("bseq_item") >> "-", "NEWLINE", "INDENT", "bmap_items", "DEDENT";
  prods[PROD_BSEQ_BMAP_TRAIL]("bseq_item") >> "-", "WS+", "NEWLINE", "INDENT", "bmap_items", "DEDENT";
  prods[PROD_BSEQ_BSEQ]("bseq_item") >> "-", "NEWLINE", "INDENT", "bseq_items", "DEDENT";
  prods[PROD_BSEQ_BSEQ_TRAIL]("bseq_item") >> "-", "WS+", "NEWLINE", "INDENT", "bseq_items", "DEDENT";
  prods[PROD_BSEQ_FMAP]("bseq_item") >> "-", "WS+", "tag?", "fmap", "NEWLINE";
  prods[PROD_BSEQ_FSEQ]("bseq_item") >> "-", "WS+", "tag?", "fseq", "NEWLINE";
  prods[PROD_FMAP]("fmap") >> "{", "WS*", "fmap_items", "}", "WS*";
  prods[PROD_FMAP_EMPTY]("fmap") >> "{", "WS*", "}", "WS*";
  prods[PROD_FMAP_FIRST]("fmap_items") >> "fmap_item";
  prods[PROD_FMAP_NEXT]("fmap_items") >> "fmap_items", ",", "WS*", "fmap_item";
  prods[PROD_FMAP_SCALAR]("fmap_item") >> "scalar", ":", "WS*", "tag?", "scalar";
  prods[PROD_FMAP_FMAP]("fmap_item") >> "scalar", ":", "WS*", "tag?", "fmap";
  prods[PROD_FMAP_FSEQ]("fmap_item") >> "scalar", ":", "WS*", "tag?", "fseq";
  prods[PROD_FSEQ]("fseq") >> "[", "WS*", "fseq_items", "]", "WS*";
  prods[PROD_FSEQ_EMPTY]("fseq") >> "[", "WS*", "]", "WS*";
  prods[PROD_FSEQ_FIRST]("fseq_items") >> "fseq_item";
  prods[PROD_FSEQ_NEXT]("fseq_items") >> "fseq_items", ",", "WS*", "fseq_item";
  prods[PROD_FSEQ_SCALAR]("fseq_item") >> "tag?", "scalar";
  prods[PROD_FSEQ_FMAP]("fseq_item") >> "tag?", "fmap";
  prods[PROD_FSEQ_FSEQ]("fseq_item") >> "tag?", "fseq";
  prods[PROD_SCALAR_RAW]("scalar") >> "scalar_head", "scalar_tail*";
  prods[PROD_SCALAR_QUOTED]("scalar") >> "scalar_quoted";
  prods[PROD_MAP_SCALAR_RAW]("map_scalar") >> "scalar_head", "scalar_tail*", "map_scalar_escaped*";
  prods[PROD_MAP_SCALAR_QUOTED]("map_scalar") >> "scalar_quoted";
  prods[PROD_SCALAR_DQUOTED]("scalar_quoted") >> "\"", "dquoted*", "descape*", "\"", "WS*";
  prods[PROD_SCALAR_SQUOTED]("scalar_quoted") >> "'", "squoted*", "sescape*", "'", "WS*";
  prods[PROD_SCALAR_HEAD_OTHER]("scalar_head") >> "OTHERCHAR";
  prods[PROD_SCALAR_HEAD_DOT]("scalar_head") >> ".", "OTHERCHAR";
  prods[PROD_SCALAR_HEAD_DASH]("scalar_head") >> "-", "OTHERCHAR";
  prods[PROD_SCALAR_HEAD_DOT_DOT]("scalar_head") >> ".", ".", "OTHERCHAR";
  prods[PROD_MAP_SCALAR_ESCAPED_EMPTY]("map_scalar_escaped*");
  prods[PROD_MAP_SCALAR_ESCAPED_NEXT]("map_scalar_escaped*") >> "map_scalar_escaped*", ",", "scalar_tail*";
  prods[PROD_TAG_EMPTY]("tag?");
  prods[PROD_TAG]("tag?") >> "!", "!", "OTHERCHAR+", "WS+";
  prods[PROD_BSCALAR]("bscalar") >> "bscalar_header", "WS*", "NEWLINE", "INDENT", "bscalar_items", "DEDENT";
  prods[PROD_BSCALAR_FIRST]("bscalar_items") >> "bscalar_item";
  prods[PROD_BSCALAR_NEXT]("bscalar_items") >> "bscalar_items", "bscalar_item";
  prods[PROD_BSCALAR_LINE]("bscalar_item") >> "any*", "NEWLINE";
  prods[PROD_BSCALAR_INDENT]("bscalar_item") >> "INDENT", "bscalar_items", "DEDENT";
  prods[PROD_BSCALAR_HEADER_LITERAL]("bscalar_header") >> "|", "bscalar_head_c*";
  prods[PROD_BSCALAR_HEADER_FOLDED]("bscalar_header") >> ">", "bscalar_head_c*";
  prods[PROD_BSCALAR_HEAD_EMPTY]("bscalar_head_c*");
  prods[PROD_BSCALAR_HEAD_NEXT]("bscalar_head_c*") >> "bscalar_head_c*", "bscalar_head_c";
  prods[PROD_BSCALAR_HEAD_OTHER]("bscalar_head_c") >> "OTHERCHAR";
  prods[PROD_BSCALAR_HEAD_DASH]("bscalar_head_c") >> "-";
  prods[PROD_DQUOTED_EMPTY]("dquoted*");
  prods[PROD_DQUOTED_NEXT]("dquoted*") >> "dquoted*", "dquoted";
  prods[PROD_SQUOTED_EMPTY]("squoted*");
  prods[PROD_SQUOTED_NEXT]("squoted*") >> "squoted*", "squoted";
  prods[PROD_ANY_EMPTY]("any*");
  prods[PROD_ANY_NEXT]("any*") >> "any*", "any";
  prods[PROD_DESCAPE_EMPTY]("descape*");
  prods[PROD_DESCAPE_NEXT]("descape*") >> "descape*", "descape";
  prods[PROD_DESCAPE]("descape") >> "\\", "descaped", "dquoted*";
  prods[PROD_SESCAPE_EMPTY]("sescape*");
  prods[PROD_SESCAPE_NEXT]("sescape*") >> "sescape*", "sescape";
  prods[PROD_SESCAPE]("sescape") >> "'", "'", "squoted*";
  prods[PROD_SCALAR_TAIL_EMPTY]("scalar_tail*");
  prods[PROD_SCALAR_TAIL_NEXT]("scalar_tail*") >> "scalar_tail*", "scalar_tail";
  prods[PROD_OTHER_FIRST]("OTHERCHAR+") >> "OTHERCHAR";
  prods[PROD_OTHER_NEXT]("OTHERCHAR+") >> "OTHERCHAR+", "OTHERCHAR";
  prods[PROD_SCALAR_TAIL_SPACE]("scalar_tail") >> "WS";
  prods[PROD_SCALAR_TAIL_DOT]("scalar_tail") >> ".";
  prods[PROD_SCALAR_TAIL_DASH]("scalar_tail") >> "-";
  prods[PROD_SCALAR_TAIL_SQUOT]("scalar_tail") >> "'";
  prods[PROD_SCALAR_TAIL_OTHER]("scalar_tail") >> "OTHERCHAR";
  prods[PROD_DESCAPED_DQUOT]("descaped") >> "\"";
  prods[PROD_DESCAPED_SLASH]("descaped") >> "\\";
  prods[PROD_DESCAPED_DQUOTED]("descaped") >> "dquoted";
  prods[PROD_DQUOTED_COMMON]("dquoted") >> "common";
  prods[PROD_DQUOTED_SQUOT]("dquoted") >> "'";
  prods[PROD_SQUOTED_COMMON]("squoted") >> "common";
  prods[PROD_SQUOTED_DQUOT]("squoted") >> "\"";
  prods[PROD_SQUOTED_SLASH]("squoted") >> "\\";
  prods[PROD_ANY_COMMON]("any") >> "common";
  prods[PROD_ANY_DQUOT]("any") >> "\"";
  prods[PROD_ANY_SQUOT]("any") >> "'";
  prods[PROD_ANY_SLASH]("any") >> "\\";
  prods[PROD_COMMON_SPACE]("common") >> "WS";
  prods[PROD_COMMON_COLON]("common") >> ":";
  prods[PROD_COMMON_DOT]("common") >> ".";
  prods[PROD_COMMON_DASH]("common") >> "-";
  prods[PROD_COMMON_PIPE]("common") >> "|";
  prods[PROD_COMMON_LSQUARE]("common") >> "[";
  prods[PROD_COMMON_RSQUARE]("common") >> "]";
  prods[PROD_COMMON_LCURLY]("common") >> "{";
  prods[PROD_COMMON_RCURLY]("common") >> "}";
  prods[PROD_COMMON_RANGLE]("common") >> ">";
  prods[PROD_COMMON_COMMA]("common") >> ",";
  prods[PROD_COMMON_PERCENT]("common") >> "%";
  prods[PROD_COMMON_EXCL]("common") >> "!";
  prods[PROD_COMMON_OTHER]("common") >> "OTHERCHAR";
  prods[PROD_SPACE_STAR_EMPTY]("WS*");
  prods[PROD_SPACE_STAR_NEXT]("WS*") >> "WS*", "WS";
  prods[PROD_SPACE_PLUS_FIRST]("WS+") >> "WS";
  prods[PROD_SPACE_PLUS_NEXT]("WS+") >> "WS+", "WS";
  return prods;
}

} // end anonymous namespace

Language make_language() {
  Language out;
  Language::Productions& prods = out.productions;
  prods = make_productions();
  Language::Tokens& toks = out.tokens;
  toks.resize(NTOKS);
  toks[TOK_NEWLINE]("NEWLINE", "((#[^\r\n]*)?\r?\n[ \t]*)+");
  toks[TOK_INDENT]("INDENT", "((#[^\r\n]*)?\r?\n[ \t]*)+");
  toks[TOK_DEDENT]("DEDENT", "((#[^\r\n]*)?\r?\n[ \t]*)+");
  toks[TOK_SPACE]("WS", "[ \t]");
  toks[TOK_COLON](":", ":");
  toks[TOK_DOT](".", "\\.");
  toks[TOK_DASH]("-", "\\-");
  toks[TOK_DQUOT]("\"", "\"");
  toks[TOK_SQUOT]("'", "'");
  toks[TOK_SLASH]("\\", "\\\\");
  toks[TOK_PIPE]("|", "\\|");
  toks[TOK_LSQUARE]("[", "\\[");
  toks[TOK_RSQUARE]("]", "\\]");
  toks[TOK_LCURLY]("{", "{");
  toks[TOK_RCURLY]("}", "}");
  toks[TOK_RANGLE](">", ">");
  toks[TOK_COMMA](",", ",");
  toks[TOK_PERCENT]("%", "%");
  toks[TOK_EXCL]("!", "!");
  toks[TOK_OTHER]("OTHERCHAR", "[^ \t:\\.\\-\"'\\\\\\|\\[\\]{}>,%#!\n\r]");
  return out;
}

LanguagePtr ask_language() {
  static LanguagePtr ptr;
  if (ptr.strong_count() == 0) {
    ptr.reset(new Language(make_language()));
  }
  return ptr;
}

ReaderTablesPtr ask_reader_tables() {
  static ReaderTablesPtr ptr;
  if (ptr.strong_count() == 0) {
    ptr = make_reader_tables(*(YAML::ask_language()));
  }
  return ptr;
}

}  // end namespace YAML
}  // end namespace Teuchos
