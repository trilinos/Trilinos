#include "Teuchos_YAML.hpp"

#include <iostream>

namespace Teuchos {
namespace YAML {

Language make_language() {
  Language out;
  Language::Productions& prods = out.productions;
  Language::Tokens& toks = out.tokens;
  prods.resize(NPRODS);
  toks.resize(NTOKS);
  toks[TOK_NEWLINE]("NEWLINE", "]NEWLINE[");
  toks[TOK_INDENT]("INDENT", "]INDENT[");
  toks[TOK_DEDENT]("DEDENT", "]DEDENT[");
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
  toks[TOK_COMMA](",", ",");
  toks[TOK_PERCENT]("%", "%");
  toks[TOK_POUND]("#", "#");
  toks[TOK_EXCL]("!", "!");
  toks[TOK_OTHER]("OTHERCHAR", "[^ \t:\\.\\-\"'\\\\\\|\\[\\]{},%#!\n\r]");
  prods[PROD_DOC]("doc") >> "top_items";
  prods[PROD_TOP_FIRST]("top_items") >> "top_item";
  prods[PROD_TOP_NEXT]("top_items") >> "top_items", "top_item";
  prods[PROD_TOP_DIRECT]("top_item") >> "%", "any*", "NEWLINE";
  prods[PROD_TOP_BEGIN]("top_item") >> "-", "-", "-", "NEWLINE", "comment*";
  prods[PROD_TOP_END]("top_item") >> ".", ".", ".", "NEWLINE";
  prods[PROD_TOP_BMAP]("top_item") >> "bmap_item";
  prods[PROD_BMAP_FIRST]("bmap_items") >> "bmap_item", "comment*";
  prods[PROD_BMAP_NEXT]("bmap_items") >> "bmap_items", "bmap_item", "comment*";
  prods[PROD_BMAP_SCALAR]("bmap_item") >> "scalar", ":", "WS*", "tag?", "scalar", "NEWLINE";
  prods[PROD_BMAP_BSCALAR]("bmap_item") >> "scalar", ":", "WS*", "bscalar";
  prods[PROD_BMAP_BMAP]("bmap_item") >> "scalar", ":", "WS*", "NEWLINE", "INDENT", "comment*", "bmap_items", "DEDENT";
  prods[PROD_BMAP_BSEQ]("bmap_item") >> "scalar", ":", "WS*", "NEWLINE", "INDENT", "comment*", "bseq_items", "DEDENT";
  prods[PROD_BMAP_FMAP]("bmap_item") >> "scalar", ":", "WS*", "tag?", "fmap", "NEWLINE";
  prods[PROD_BMAP_FSEQ]("bmap_item") >> "scalar", ":", "WS*", "tag?", "fseq", "NEWLINE";
  prods[PROD_BSEQ_FIRST]("bseq_items") >> "bseq_item", "comment*";
  prods[PROD_BSEQ_NEXT]("bseq_items") >> "bseq_items", "bseq_item", "comment*";
  prods[PROD_BSEQ_SCALAR]("bseq_item") >> "-", "WS+", "tag?", "scalar", "NEWLINE";
  prods[PROD_BSEQ_BSCALAR]("bseq_item") >> "-", "WS+", "bscalar";
  prods[PROD_BSEQ_BMAP]("bseq_item") >> "-", "NEWLINE", "INDENT", "comment*", "bmap_items", "DEDENT";
  prods[PROD_BSEQ_BMAP_TRAIL]("bseq_item") >> "-", "WS+", "NEWLINE", "INDENT", "comment*", "bmap_items", "DEDENT";
  prods[PROD_BSEQ_BSEQ]("bseq_item") >> "-", "NEWLINE", "INDENT", "comment*", "bseq_items", "DEDENT";
  prods[PROD_BSEQ_BSEQ_TRAIL]("bseq_item") >> "-", "WS+", "NEWLINE", "INDENT", "comment*", "bseq_items", "DEDENT";
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
  prods[PROD_SCALAR_NORMAL]("scalar") >> "OTHERCHAR", "rest*";
  prods[PROD_SCALAR_DOT]("scalar") >> ".", "OTHERCHAR", "rest*";
  prods[PROD_SCALAR_DASH]("scalar") >> "-", "OTHERCHAR", "rest*";
  prods[PROD_SCALAR_DQUOTED]("scalar") >> "\"", "dquoted*", "descape*", "\"";
  prods[PROD_SCALAR_SQUOTED]("scalar") >> "'", "squoted*", "sescape*", "'";
  prods[PROD_COMMENT_EMPTY]("comment*");
  prods[PROD_COMMENT_NEXT]("comment*") >> "comment*", "#", "any*", "NEWLINE";
  prods[PROD_TAG_EMPTY]("tag?");
  prods[PROD_TAG]("tag?") >> "!", "!", "OTHERCHAR+", "WS+";
  prods[PROD_BSCALAR]("bscalar") >> "|", "WS*", "NEWLINE", "INDENT", "bscalar_items", "DEDENT";
  prods[PROD_BSCALAR_FIRST]("bscalar_items") >> "bscalar_item";
  prods[PROD_BSCALAR_NEXT]("bscalar_items") >> "bscalar_items", "bscalar_item";
  prods[PROD_BSCALAR_LINE]("bscalar_item") >> "any*", "NEWLINE";
  prods[PROD_BSCALAR_INDENT]("bscalar_item") >> "INDENT", "bscalar_items", "DEDENT";
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
  prods[PROD_REST_EMPTY]("rest*");
  prods[PROD_REST_NEXT]("rest*") >> "rest*", "rest";
  prods[PROD_OTHER_FIRST]("OTHERCHAR+") >> "OTHERCHAR";
  prods[PROD_OTHER_NEXT]("OTHERCHAR+") >> "OTHERCHAR+", "OTHERCHAR";
  prods[PROD_REST_SPACE]("rest") >> "WS";
  prods[PROD_REST_DOT]("rest") >> ".";
  prods[PROD_REST_DASH]("rest") >> "-";
  prods[PROD_REST_OTHER]("rest") >> "OTHERCHAR";
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
  prods[PROD_COMMON_COMMA]("common") >> ",";
  prods[PROD_COMMON_PERCENT]("common") >> "%";
  prods[PROD_COMMON_POUND]("common") >> "#";
  prods[PROD_COMMON_EXCL]("common") >> "!";
  prods[PROD_COMMON_OTHER]("common") >> "OTHERCHAR";
  prods[PROD_SPACE_STAR_EMPTY]("WS*");
  prods[PROD_SPACE_START_NEXT]("WS*") >> "WS*", "WS";
  prods[PROD_SPACE_PLUS_EMPTY]("WS+") >> "WS";
  prods[PROD_SPACE_PLUS_NEXT]("WS+") >> "WS+", "WS";
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
