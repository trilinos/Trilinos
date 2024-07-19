// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_XML.hpp"

namespace Teuchos {
namespace XML {

Language make_language() {
  Language out;
  Language::Productions& prods = out.productions;
  prods.resize(NPRODS);
  prods[PROD_DOC]("document") >> "toplevels";
  prods[PROD_TOPLEVEL]("toplevels") >> "toplevel";
  prods[PROD_TOPLEVELS]("toplevels") >> "toplevel", "toplevels";
  prods[PROD_TOPLEVELS_MISC]("toplevels") >> "toplevel", "Misc", "toplevels";
  prods[PROD_TOPLEVEL_ELEMENT]("toplevel") >> "element";
  prods[PROD_TOPLEVEL_XMLDECL]("toplevel") >> "XMLDecl";
  prods[PROD_ELEMENT_EMPTY]("element") >> "EmptyElemTag";
  prods[PROD_ELEMENT]("element") >> "STag", "content";
  prods[PROD_XMLDECL]("XMLDecl") >> "<", "?", "Name", "TagFill", "?", ">";
  prods[PROD_STAG]("STag") >> "<", "Name", "TagFill", ">";
  prods[PROD_ETAG]("ETag") >> "<", "/", "Name", "S?", ">";
  prods[PROD_EMPTY_TAG]("EmptyElemTag") >> "<", "Name", "TagFill", "/", ">";
  prods[PROD_CONTENT]("content") >> "CharData?", "ContentItem*", "ETag";
  prods[PROD_NO_CONTENT_ITEMS]("ContentItem*");
  prods[PROD_CONTENT_ITEMS]("ContentItem*") >> "ContentItem*", "ContentItem", "CharData?";
  prods[PROD_CONTENT_ELEMENT]("ContentItem") >> "element";
  prods[PROD_CONTENT_REF]("ContentItem") >> "Reference";
  prods[PROD_CONTENT_COMMENT]("ContentItem") >> "Comment";
  prods[PROD_NO_CHARDATA]("CharData?");
  prods[PROD_CHARDATA]("CharData?") >> "CharData?", "DataChar";
  prods[PROD_TAGFILL]("TagFill") >> "Attributes", "S?";
  prods[PROD_NO_ATTS]("Attributes");
  prods[PROD_ATTS]("Attributes") >> "Attributes", "S", "Attribute";
  prods[PROD_ATT]("Attribute") >> "Name", "Eq", "AttValue";
  prods[PROD_EQ]("Eq") >> "S?", "=", "S?";
  prods[PROD_ATTVALUE_D]("AttValue") >> "\"", "DQuoteds", "\"";
  prods[PROD_ATTVALUE_S]("AttValue") >> "'", "SQuoteds", "'";
  prods[PROD_NO_DQUOTS]("DQuoteds");
  prods[PROD_DQUOTS]("DQuoteds") >> "DQuoteds", "DQuoted";
  prods[PROD_DQUOT_CHAR]("DQuoted") >> "DQuotedChar";
  prods[PROD_DQUOT_REF]("DQuoted") >> "Reference";
  prods[PROD_NO_SQUOTS]("SQuoteds");
  prods[PROD_SQUOTS]("SQuoteds") >> "SQuoteds", "SQuoted";
  prods[PROD_SQUOT_CHAR]("SQuoted") >> "SQuotedChar";
  prods[PROD_SQUOT_REF]("SQuoted") >> "Reference";
  prods[PROD_NAME]("Name") >> "NameFirstChar", "NameChars";
  prods[PROD_NAME_FIRST_LETTER]("NameFirstChar") >> "Letter";
  prods[PROD_NAME_FIRST_UNDER]("NameFirstChar") >> "_";
  prods[PROD_NAME_FIRST_COLON]("NameFirstChar") >> ":";
  prods[PROD_NO_NAME_CHARS]("NameChars");
  prods[PROD_NAME_CHARS]("NameChars") >> "NameChars", "NameChar";
  prods[PROD_NAME_LETTER]("NameChar") >> "Letter";
  prods[PROD_NAME_DIGIT]("NameChar") >> "Digit";
  prods[PROD_NAME_DOT]("NameChar") >> ".";
  prods[PROD_NAME_DASH]("NameChar") >> "-";
  prods[PROD_NAME_UNDER]("NameChar") >> "_";
  prods[PROD_NAME_COLON]("NameChar") >> ":";
  prods[PROD_NO_MISCS]("Miscs");
  prods[PROD_MISCS]("Miscs") >> "Miscs", "Misc";
  prods[PROD_MISC_COMMENT]("Misc") >> "Comment";
  prods[PROD_MISC_SPACE]("Misc") >> "S";
  prods[PROD_COMMENT]("Comment") >> "<", "!", "-", "-", "Commenteds", "-", "-", ">";
  prods[PROD_NO_COMMENTED]("Commenteds");
  prods[PROD_COMMENTED]("Commenteds") >> "Commenteds", "Commented";
  prods[PROD_COMMENT_CHAR]("Commented") >> "CommentChar";
  prods[PROD_COMMENT_DASH]("Commented") >> "-", "CommentChar";
  prods[PROD_ENT_REF]("Reference") >> "&", "Name", ";";
  prods[PROD_CHAR_REF]("Reference") >> "&", "#", "Digits", ";";
  prods[PROD_ONE_DIGIT]("Digits") >> "Digit";
  prods[PROD_DIGITS]("Digits") >> "Digits", "Digit";
  prods[PROD_NO_SPACES]("S?");
  prods[PROD_YES_SPACES]("S?") >> "S";
  prods[PROD_ONE_SPACE]("S") >> "Space";
  prods[PROD_SPACES]("S") >> "S", "Space";
  prods[PROD_DQUOTED_COMMON]("DQuotedChar") >> "CommonChar";
  prods[PROD_DQUOTED_SQUOT]("DQuotedChar") >> "'";
  prods[PROD_DQUOTED_RSQUARE]("DQuotedChar") >> "]";
  prods[PROD_DQUOTED_DASH]("DQuotedChar") >> "-";
  prods[PROD_SQUOTED_CHAR]("SQuotedChar") >> "CommonChar";
  prods[PROD_SQUOTED_DQUOT]("SQuotedChar") >> "\"";
  prods[PROD_SQUOTED_RSQUARE]("SQuotedChar") >> "]";
  prods[PROD_SQUOTED_DASH]("SQuotedChar") >> "-";
  prods[PROD_DATA_COMMON]("DataChar") >> "CommonChar";
  prods[PROD_DATA_SQUOT]("DataChar") >> "'";
  prods[PROD_DATA_DQUOT]("DataChar") >> "\"";
  prods[PROD_DATA_DASH]("DataChar") >> "-";
  prods[PROD_COMMENT_COMMON]("CommentChar") >> "CommonChar"; 
  prods[PROD_COMMENT_LANGLE]("CommentChar") >> "<"; 
  prods[PROD_COMMENT_AMP]("CommentChar") >> "&"; 
  prods[PROD_COMMENT_SQUOT]("CommentChar") >> "'"; 
  prods[PROD_COMMENT_DQUOT]("CommentChar") >> "\""; 
  prods[PROD_COMMENT_RSQUARE]("CommentChar") >> "]"; 
  prods[PROD_COMMON_SPACE]("CommonChar") >> "Space";
  prods[PROD_COMMON_LETTER]("CommonChar") >> "Letter";
  prods[PROD_COMMON_DIGIT]("CommonChar") >> "Digit";
  prods[PROD_COMMON_EXCL]("CommonChar") >> "!";
  prods[PROD_COMMON_POUND]("CommonChar") >> "#";
  prods[PROD_COMMON_DOT]("CommonChar") >> ".";
  prods[PROD_COMMON_SLASH]("CommonChar") >> "/";
  prods[PROD_COMMON_COLON]("CommonChar") >> ":";
  prods[PROD_COMMON_SEMICOLON]("CommonChar") >> ";";
  prods[PROD_COMMON_RANGLE]("CommonChar") >> ">";
  prods[PROD_COMMON_QUESTION]("CommonChar") >> "?";
  prods[PROD_COMMON_EQUAL]("CommonChar") >> "=";
  prods[PROD_COMMON_LSQUARE]("CommonChar") >> "[";
  prods[PROD_COMMON_UNDER]("CommonChar") >> "_";
  prods[PROD_COMMON_OTHER]("CommonChar") >> "OtherChar";
  Language::Tokens& toks = out.tokens;
  toks.resize(NTOKS);
  toks[TOK_SPACE]("Space", "[ \t\r\n]");
  toks[TOK_LETTER]("Letter", "[a-zA-Z]");
  toks[TOK_DIGIT]("Digit", "[0-9]");
  toks[TOK_EXCL]("!", "!");
  toks[TOK_DQUOTE]("\"", "\"");
  toks[TOK_SQUOTE]("'", "'");
  toks[TOK_POUND]("#", "#");
  toks[TOK_AMP]("&", "&");
  toks[TOK_DASH]("-", "\\-");
  toks[TOK_DOT](".", "\\.");
  toks[TOK_SLASH]("/", "/");
  toks[TOK_COLON](":", ":");
  toks[TOK_SEMICOLON](";", ";");
  toks[TOK_LANGLE]("<", "<");
  toks[TOK_RANGLE](">", ">");
  toks[TOK_QUESTION]("?", "\\?");
  toks[TOK_EQUAL]("=", "=");
  toks[TOK_LSQUARE]("[", "\\[");
  toks[TOK_RSQUARE]("]", "\\]");
  toks[TOK_UNDER]("_", "_");
  toks[TOK_OTHER]("OtherChar", "[$%\\(\\)\\*\\+,@\\\\\\^`{}\\|~]");
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
    LanguagePtr lang = ask_language();
    ptr = make_reader_tables(*lang);
  }
  return ptr;
}

}  // end namespace XML
}  // end namespace Teuchos
