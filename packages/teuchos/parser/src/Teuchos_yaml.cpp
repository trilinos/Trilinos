#include "Teuchos_yaml.hpp"

namespace Teuchos {
namespace yaml {

Language build_language() {
  Language out;
  auto& prods = out.productions;
  auto& toks = out.tokens;
  prods.push_back({"doc", {"EQDENT", "prolog", "toplevel", "epilog"}});
  prods.push_back({"prolog", {"directive?", "DOC_START"}});
  prods.push_back({"directive?", {}});
  prods.push_back({"directive?", {"directive"}});
  prods.push_back({"directive", {"DIRECTIVE", "EQDENT"}});
  prods.push_back({"epilog", {"doc_end"}});
  prods.push_back({"doc_end", {"DOC_END", "EQDENT"}});
  prods.push_back({"toplevel", {"EQDENT", "block_map_items", "EQDENT?"}});
  prods.push_back({"toplevel", {"EQDENT", "block_sequence_items", "EQDENT?"}});
  prods.push_back({"toplevel", {"block_collective"}});
  prods.push_back({"block_collective", {"INDENT", "block_map_items", "DEDENT"}});
  prods.push_back({"block_collective", {"INDENT", "block_sequence_items", "DEDENT"}});
  prods.push_back({"block_map_items", {"block_map_item"}});
  prods.push_back({"block_map_items", {"block_map_items", "EQDENT", "block_map_item"}});
  prods.push_back({"block_sequence_items", {"block_sequence_item"}});
  prods.push_back({"block_sequence_items", {"block_sequence_items", "EQDENT", "block_sequence_item"}});
  prods.push_back({"block_sequence_item", {"BLOCK_SEQ", "scalar", "S?"}});
  prods.push_back({"block_map_item", {"scalar", "S?", ":", "S?", "block_map_value"}});
  prods.push_back({"block_map_value", {"scalar", "S?"}});
  prods.push_back({"block_map_value", {"block_collective"}});
  prods.push_back({"block_map_value", {"flow_collective"}});
  prods.push_back({"flow_collective", {"[", "flow_sequence_items", "]"}});
  prods.push_back({"flow_collective", {"{", "flow_map_items", "}"}});
  prods.push_back({"flow_sequence_items", {"flow_sequence_item"}});
  prods.push_back({"flow_sequence_items", {"flow_sequence_items", "FLOW_SEP", "S?", "flow_sequence_item"}});
  prods.push_back({"flow_sequence_item", {"scalar", "S?"}});
  prods.push_back({"flow_sequence_item", {"flow_collective", "S?"}});
  prods.push_back({"flow_map_items", {"flow_map_item"}});
  prods.push_back({"flow_map_items", {"flow_map_items", "FLOW_SEP", "S?", "flow_map_item"}});
  prods.push_back({"flow_map_item", {"scalar", "S?", ":", "S?", "flow_map_value"}});
  prods.push_back({"flow_map_value", {"scalar", "S?"}});
  prods.push_back({"flow_map_value", {"flow_collective", "S?"}});
  prods.push_back({"S?", {}});
  prods.push_back({"S?", {"S"}});
  prods.push_back({"EQDENT?", {}});
  prods.push_back({"EQDENT?", {"EQDENT"}});
  prods.push_back({"scalar", {"RAW_SCALAR"}});
  prods.push_back({"scalar", {"DOUBLE_QUOTED"}});
  prods.push_back({"scalar", {"SINGLE_QUOTED"}});
  toks.push_back({"NODENT", "]NODENT["});
  toks.push_back({"INDENT", "]INDENT["});
  toks.push_back({"DEDENT", "]DEDENT["});
  toks.push_back({"EQDENT", "]EQDENT["});
  toks.push_back({"S", "[ \t]+"});
  toks.push_back({":", ":"});
  toks.push_back({"DOC_START", "\\-\\-\\-"});
  toks.push_back({"DOC_END", "\\.\\.\\."});
  toks.push_back({"DIRECTIVE", "%[^\n\r]*"});
  toks.push_back({"RAW_SCALAR", "(\\-[^ \t:\n\r%'\",{}\\[\\]]|[^ \t:\n\r\\-%'\",\\[\\]{}])([^:\n\r%'\",{}\\[\\]]*[^ \t:\n\r%'\",\\[\\]{}])?"});
  toks.push_back({"DOUBLE_QUOTED", "\"(\\\\.|[^\\\\\"])*\""});
  toks.push_back({"SINGLE_QUOTED", "'[^']*(''[^']*)*'"});
  toks.push_back({"BLOCK_SEQ", "\\-[ \t]+"});
  toks.push_back({"FLOW_SEP", ","});
  toks.push_back({"[", "\\["});
  toks.push_back({"]", "\\]"});
  toks.push_back({"{", "{"});
  toks.push_back({"}", "}"});
  return out;
}

LanguagePtr ask_language() {
  static LanguagePtr ptr;
  if (ptr.use_count() == 0) {
    ptr.reset(new Language(build_language()));
  }
  return ptr;
}

ReaderTablesPtr ask_reader_tables() {
  static ReaderTablesPtr ptr;
  if (ptr.use_count() == 0) {
    ptr = build_reader_tables(*(yaml::ask_language()));
  }
  return ptr;
}

}  // end namespace yaml
}  // end namespace Teuchos
