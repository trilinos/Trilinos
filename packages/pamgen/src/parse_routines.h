// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef parse_routinesH
namespace PAMGEN_NEVADA {

Token Parse_User_Defined_Element_Density(Token_Stream *token_stream, int);
Token Parse_User_Defined_Geometry_Transformation(Token_Stream *token_stream, int);
Token Parse_Decomposition_Strategy(Token_Stream *token_stream, int);
Token Parse_Topology_Modification(Token_Stream *token_stream, int);
Token Parse_Inline_Mesh_Tok(Token_Stream *token_stream, int value);
Token Parse_Inline_Mesh_3D_Tok(Token_Stream *token_stream, int value);
Token Parse_Inline_Mesh_2D_Tok(Token_Stream *token_stream, int value);

void Allow_New_Mesh_Specification();

}//end of namespace PAMGEN_NEVADA
#define parse_routinesH
#endif
