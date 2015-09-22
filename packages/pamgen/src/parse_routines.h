#ifndef parse_routinesH
namespace PAMGEN_NEVADA {

Token Parse_User_Defined_Element_Density(Token_Stream *token_stream, int);
Token Parse_User_Defined_Geometry_Transformation(Token_Stream *token_stream, int);
Token Parse_Decomposition_Strategy(Token_Stream *token_stream, int);
Token Parse_Topology_Modification(Token_Stream *token_stream, int);
Token Parse_Inline_Mesh_Tok(Token_Stream *token_stream, int value);
Token Parse_Inline_Mesh_3D_Tok(Token_Stream *token_stream, int value);
Token Parse_Inline_Mesh_2D_Tok(Token_Stream *token_stream, int value);
}//end of namespace PAMGEN_NEVADA
#define parse_routinesH
#endif
