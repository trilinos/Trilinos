// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_YAML_HPP
#define TEUCHOS_YAML_HPP

/*! \file Teuchos_YAML.hpp
    \brief A TeuchosParser Language for a subset of YAML

This is a grammar for a subset of the YAML language.
Since YAML is indentation-sensitive, it is not context-free.
An extension has been made to Teuchos::Reader to emit
INDENT and DEDENT tokens in order to make indentation detectable
by a context-free grammar.

Known limitations of this grammar compared to the full YAML language,
in particular those which are relevant to Teuchos::ParameterList:

<ol>
<li> INDENT and DEDENT tokens must be matched, meaning the indentation of
the YAML file itself must be nicely nested.
Examples of valid YAML that this grammar cannot handle, and workarounds:

\code{.yaml}
# nested block sequences compressed onto one line:
- - - one
    - two
    - three
  - ten
  - twenty
- one hundred
# do not compress the first line:
-
  -
    - one
    - two
    - three
  - ten
  - twenty
- one hundred
\endcode
\code{.yaml}
# comments not indented the same as subsequent lines:
all:
   one: 1
   two: 2
#bad comment
   three: 3
# indent the comments:
all:
   one: 1
   two: 2
   #good comment
   three: 3
\endcode
\code{.yaml}
# flow sequences and maps that span multiple lines
timesteps: [1, 2, 3,
   4, 5]
# keep it on one line, or use block sequences
timesteps: [1, 2, 3, 4, 5]
timesteps:
  - 1
  - 2
  - 3
  - 4
  - 5
\endcode
\code{.yaml}
# block scalars with non-nested indentation
inline file: |
   if (a == 5) {
      strange();
    indentation();
  }
# ensure block scalars have nested indentation
inline file: |
  if (a == 5) {
    strange();
    indentation();
  }
\endcode

<li> Scalars can start with a '.' or '-', but not two such symbols in a row

\code{.yaml}
# ".." at the beginning of a scalar
filepath: ../cube.exo
# quote the scalar
filepath: "../cube.exo"
\endcode

</ol>
*/

#include <Teuchos_Language.hpp>
#include <Teuchos_ReaderTables.hpp>

namespace Teuchos {
namespace YAML {

enum {
  PROD_DOC,
  PROD_DOC2,
  PROD_TOP_FIRST,
  PROD_TOP_NEXT,
  PROD_TOP_DIRECT,
  PROD_TOP_BEGIN,
  PROD_TOP_END,
  PROD_TOP_BMAP,
  PROD_BMAP_FIRST,
  PROD_BMAP_NEXT,
  PROD_BMAP_SCALAR,
  PROD_BMAP_BSCALAR,
  PROD_BMAP_BVALUE,
  PROD_BVALUE_EMPTY,
  PROD_BVALUE_BMAP,
  PROD_BVALUE_BSEQ,
  PROD_BMAP_FMAP,
  PROD_BMAP_FSEQ,
  PROD_BSEQ_FIRST,
  PROD_BSEQ_NEXT,
  PROD_BSEQ_SCALAR,
  PROD_BSEQ_BSCALAR,
  PROD_BSEQ_BMAP,
  PROD_BSEQ_BMAP_TRAIL,
  PROD_BSEQ_BSEQ,
  PROD_BSEQ_BSEQ_TRAIL,
  PROD_BSEQ_FMAP,
  PROD_BSEQ_FSEQ,
  PROD_FMAP,
  PROD_FMAP_EMPTY,
  PROD_FMAP_FIRST,
  PROD_FMAP_NEXT,
  PROD_FMAP_SCALAR,
  PROD_FMAP_FMAP,
  PROD_FMAP_FSEQ,
  PROD_FSEQ,
  PROD_FSEQ_EMPTY,
  PROD_FSEQ_FIRST,
  PROD_FSEQ_NEXT,
  PROD_FSEQ_SCALAR,
  PROD_FSEQ_FMAP,
  PROD_FSEQ_FSEQ,
  PROD_SCALAR_RAW,
  PROD_SCALAR_QUOTED,
  PROD_MAP_SCALAR_RAW,
  PROD_MAP_SCALAR_QUOTED,
  PROD_SCALAR_DQUOTED,
  PROD_SCALAR_SQUOTED,
  PROD_SCALAR_HEAD_OTHER,
  PROD_SCALAR_HEAD_DOT,
  PROD_SCALAR_HEAD_DASH,
  PROD_SCALAR_HEAD_DOT_DOT,
  PROD_MAP_SCALAR_ESCAPED_EMPTY,
  PROD_MAP_SCALAR_ESCAPED_NEXT,
  PROD_TAG_EMPTY,
  PROD_TAG,
  PROD_BSCALAR,
  PROD_BSCALAR_FIRST,
  PROD_BSCALAR_NEXT,
  PROD_BSCALAR_LINE,
  PROD_BSCALAR_INDENT,
  PROD_BSCALAR_HEADER_LITERAL,
  PROD_BSCALAR_HEADER_FOLDED,
  PROD_BSCALAR_HEAD_EMPTY,
  PROD_BSCALAR_HEAD_NEXT,
  PROD_BSCALAR_HEAD_OTHER,
  PROD_BSCALAR_HEAD_DASH,
  PROD_DQUOTED_EMPTY,
  PROD_DQUOTED_NEXT,
  PROD_SQUOTED_EMPTY,
  PROD_SQUOTED_NEXT,
  PROD_ANY_EMPTY,
  PROD_ANY_NEXT,
  PROD_DESCAPE_EMPTY,
  PROD_DESCAPE_NEXT,
  PROD_DESCAPE,
  PROD_SESCAPE_EMPTY,
  PROD_SESCAPE_NEXT,
  PROD_SESCAPE,
  PROD_SCALAR_TAIL_EMPTY,
  PROD_SCALAR_TAIL_NEXT,
  PROD_OTHER_FIRST,
  PROD_OTHER_NEXT,
  PROD_SCALAR_TAIL_SPACE,
  PROD_SCALAR_TAIL_DOT,
  PROD_SCALAR_TAIL_DASH,
  PROD_SCALAR_TAIL_SQUOT,
  PROD_SCALAR_TAIL_OTHER,
  PROD_DESCAPED_DQUOT,
  PROD_DESCAPED_SLASH,
  PROD_DESCAPED_DQUOTED,
  PROD_DQUOTED_COMMON,
  PROD_DQUOTED_SQUOT,
  PROD_SQUOTED_COMMON,
  PROD_SQUOTED_DQUOT,
  PROD_SQUOTED_SLASH,
  PROD_ANY_COMMON,
  PROD_ANY_DQUOT,
  PROD_ANY_SQUOT,
  PROD_ANY_SLASH,
  PROD_COMMON_SPACE,
  PROD_COMMON_COLON,
  PROD_COMMON_DOT,
  PROD_COMMON_DASH,
  PROD_COMMON_PIPE,
  PROD_COMMON_LSQUARE,
  PROD_COMMON_RSQUARE,
  PROD_COMMON_LCURLY,
  PROD_COMMON_RCURLY,
  PROD_COMMON_RANGLE,
  PROD_COMMON_COMMA,
  PROD_COMMON_PERCENT,
  PROD_COMMON_EXCL,
  PROD_COMMON_OTHER,
  PROD_SPACE_STAR_EMPTY,
  PROD_SPACE_STAR_NEXT,
  PROD_SPACE_PLUS_FIRST,
  PROD_SPACE_PLUS_NEXT
};

enum { NPRODS = PROD_SPACE_PLUS_NEXT + 1 };

enum {
  TOK_NEWLINE,
  TOK_INDENT,
  TOK_DEDENT,
  TOK_SPACE,
  TOK_COLON,
  TOK_DOT,
  TOK_DASH,
  TOK_DQUOT,
  TOK_SQUOT,
  TOK_SLASH,
  TOK_PIPE,
  TOK_LSQUARE,
  TOK_RSQUARE,
  TOK_LCURLY,
  TOK_RCURLY,
  TOK_RANGLE,
  TOK_COMMA,
  TOK_PERCENT,
  TOK_EXCL,
  TOK_OTHER
};

enum { NTOKS = TOK_OTHER + 1 };

Language make_language();
LanguagePtr ask_language();
ReaderTablesPtr ask_reader_tables();

}  // end namespace yaml
}  // end namespace Teuchos

#endif
