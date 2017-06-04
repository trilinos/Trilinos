#include "Teuchos_YAML.hpp"

namespace Teuchos {
namespace yaml {

Language make_language() {
  Language out;
  Language::Productions& prods = out.productions;
  Language::Tokens& toks = out.tokens;
  prods.resize(NPRODS);
  toks.resize(NTOKS);
  prods[PROD_DOC]("doc") >> "EQDENT", "prolog", "toplevel", "epilog";
  prods[PROD_PROLOG]("prolog") >> "directive?", "DOC_START";
  prods[PROD_NO_DIRECT]("directive?");
  prods[PROD_ONE_DIRECT]("directive?") >> "directive";
  prods[PROD_DIRECT]("directive") >> "DIRECTIVE", "EQDENT";
  prods[PROD_EPILOG]("epilog") >> "doc_end";
  prods[PROD_DOC_END]("doc_end") >> "DOC_END", "EQDENT";
  prods[PROD_TOP_BMAP]("toplevel") >> "EQDENT", "block_map_items", "EQDENT?";
  prods[PROD_TOP_BSEQ]("toplevel") >> "EQDENT", "block_sequence_items", "EQDENT?";
  prods[PROD_TOP_BLOCK]("toplevel") >> "block_collective";
  prods[PROD_BMAP]("block_collective") >> "INDENT", "block_map_items", "DEDENT";
  prods[PROD_BSEQ]("block_collective") >> "INDENT", "block_sequence_items", "DEDENT";
  prods[PROD_BMAP_ONE_ITEM]("block_map_items") >> "block_map_item";
  prods[PROD_BMAP_ITEMS]("block_map_items") >> "block_map_items", "EQDENT", "block_map_item";
  prods[PROD_BSEQ_ITEM]("block_sequence_items") >> "block_sequence_item";
  prods[PROD_BSEQ_ITEMS]("block_sequence_items") >> "block_sequence_items", "EQDENT", "block_sequence_item";
  prods[PROD_BSEQ_SCALAR]("block_sequence_item") >> "comments", "BLOCK_SEQ", "scalar", "S?";
  prods[PROD_BMAP_ITEM]("block_map_item") >> "comments", "scalar", "S?", ":", "S?", "block_map_value";
  prods[PROD_BMAP_SCALAR]("block_map_value") >> "scalar", "S?";
  prods[PROD_BMAP_BLOCK]("block_map_value") >> "block_collective";
  prods[PROD_BMAP_FLOW]("block_map_value") >> "flow_collective";
  prods[PROD_FSEQ_EMPTY]("flow_collective") >> "[", "S?", "]";
  prods[PROD_FMAP_EMPTY]("flow_collective") >> "{", "S?", "}";
  prods[PROD_FSEQ]("flow_collective") >> "[", "S?", "flow_sequence_items", "]";
  prods[PROD_FMAP]("flow_collective") >> "{", "S?", "flow_map_items", "}";
  prods[PROD_FSEQ_ITEM]("flow_sequence_items") >> "flow_sequence_item";
  prods[PROD_FSEQ_ITEMS]("flow_sequence_items") >> "flow_sequence_items", "FLOW_SEP", "S?", "flow_sequence_item";
  prods[PROD_FSEQ_SCALAR]("flow_sequence_item") >> "scalar", "S?";
  prods[PROD_FSEQ_FLOW]("flow_sequence_item") >> "flow_collective", "S?";
  prods[PROD_FMAP_ONE_ITEM]("flow_map_items") >> "flow_map_item";
  prods[PROD_FMAP_ITEMS]("flow_map_items") >> "flow_map_items", "FLOW_SEP", "S?", "flow_map_item";
  prods[PROD_FMAP_ITEM]("flow_map_item") >> "scalar", "S?", ":", "S?", "flow_map_value";
  prods[PROD_FMAP_SCALAR]("flow_map_value") >> "scalar", "S?";
  prods[PROD_FMAP_FLOW]("flow_map_value") >> "flow_collective", "S?";
  prods[PROD_NO_SPACE]("S?");
  prods[PROD_SPACE]("S?") >> "S";
  prods[PROD_NO_COMMENTS]("comments");
  prods[PROD_COMMENTS]("comments") >> "comments", "COMMENT", "EQDENT";
  prods[PROD_NO_EQDENT]("EQDENT?");
  prods[PROD_EQDENT]("EQDENT?") >> "EQDENT";
  prods[PROD_RAW]("scalar") >> "RAW_SCALAR";
  prods[PROD_DQUOTED]("scalar") >> "DOUBLE_QUOTED";
  prods[PROD_SQUOTED]("scalar") >> "SINGLE_QUOTED";
  toks[TOK_NODENT]("NODENT", "]NODENT[");
  toks[TOK_INDENT]("INDENT", "]INDENT[");
  toks[TOK_DEDENT]("DEDENT", "]DEDENT[");
  toks[TOK_EQDENT]("EQDENT", "]EQDENT[");
  toks[TOK_SPACES]("S", "[ \t]+");
  toks[TOK_COMMENT]("COMMENT", "#[^\n\r]*");
  toks[TOK_COLON](":", ":");
  toks[TOK_DOC_START]("DOC_START", "\\-\\-\\-");
  toks[TOK_DOC_END]("DOC_END", "\\.\\.\\.");
  toks[TOK_DIRECT]("DIRECTIVE", "%[^\n\r]*");
  toks[TOK_RAW]("RAW_SCALAR", "(\\-[^ \t:\n\r%'\",{}\\[\\]]|[^ \t:\n\r\\-%'\",\\[\\]{}])([^:\n\r%'\",{}\\[\\]]*[^ \t:\n\r%'\",\\[\\]{}])?");
  toks[TOK_DQUOTED]("DOUBLE_QUOTED", "\"(\\\\.|[^\\\\\"])*\"");
  toks[TOK_SQUOTED]("SINGLE_QUOTED", "'[^']*(''[^']*)*'");
  toks[TOK_BSEQ]("BLOCK_SEQ", "\\-[ \t]+");
  toks[TOK_FSEP]("FLOW_SEP", ",");
  toks[TOK_LSQUARE]("[", "\\[");
  toks[TOK_RSQUARE]("]", "\\]");
  toks[TOK_LCURLY]("{", "{");
  toks[TOK_RCURLY]("}", "}");
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
    ptr = make_reader_tables(*(yaml::ask_language()));
  }
  return ptr;
}

}  // end namespace yaml
}  // end namespace Teuchos
