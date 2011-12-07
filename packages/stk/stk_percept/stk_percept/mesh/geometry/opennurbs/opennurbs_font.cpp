/* $NoKeywords: $ */
/*
//
// Copyright (c) 1993-2007 Robert McNeel & Associates. All rights reserved.
// Rhinoceros is a registered trademark of Robert McNeel & Assoicates.
//
// THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY.
// ALL IMPLIED WARRANTIES OF FITNESS FOR ANY PARTICULAR PURPOSE AND OF
// MERCHANTABILITY ARE HEREBY DISCLAIMED.
//				
// For complete openNURBS copyright information see <http://www.opennurbs.org>.
//
////////////////////////////////////////////////////////////////
*/

#include "opennurbs.h"

//srkenno@sandia.gov: split this into two files (see opennurbs_font_1.cpp) to workaround a potential Mac/Darwin gcc4.4 compiler/linker bug

ON_OBJECT_IMPLEMENT( ON_Font, ON_Object, "4F0F51FB-35D0-4865-9998-6D2C6A99721D" );

ON_Font::ON_Font()
{
  Defaults();
}

ON_Font::~ON_Font()
{
}

void ON_Font::Defaults()
{
  m_font_name.Empty();
  m_font_weight = 0;
  m_font_italic = false;
  m_font_underlined = false;
  m_linefeed_ratio = m_default_linefeed_ratio;
  m_font_index = -1;
  memset(&m_font_id,0,sizeof(m_font_id));
  memset( &m_facename, 0, sizeof( m_facename));
  m_I_height = 0;
#if defined(ON_OS_WINDOWS_GDI)
  memset(&m_logfont,0,sizeof(m_logfont));
  m_logfont.lfHeight = normal_font_height;
  m_logfont.lfCharSet = default_charset;
#endif
  SetFontFaceName(L"Arial");
  SetFontWeight(normal_weight);
}

//////////////////////////////////////////////////////////////////////
//
// ON_Object overrides

ON_BOOL32 ON_Font::IsValid( ON_TextLog* text_log ) const
{
  return ( m_font_name.Length() > 0 
           && m_font_index >= 0 
           && m_facename[0] > 32 
           && m_facename[64] == 0
           );
}

void ON_Font::Dump( ON_TextLog& dump ) const
{
  const wchar_t* name = FontName();
  if ( !name )
    name = L"";
  dump.Print("font index = %d\n",m_font_index);
  dump.Print("font name = \"%S\"\n",name);
  dump.Print("font face name = \"%S\"\n",m_facename);
  dump.Print("font weight = \"%d\"\n",m_font_weight);
  dump.Print("font is italic = \"%d\"\n",m_font_italic);
  dump.Print("font is underlined = \"%d\"\n",m_font_underlined);
  dump.Print("font linefeed ratio = \"%g\"\n", m_linefeed_ratio);
}

ON_BOOL32 ON_Font::Write(
       ON_BinaryArchive& file // serialize definition to binary archive
     ) const
{
  bool rc = file.Write3dmChunkVersion(1,2);
  while(rc)
  {
    rc = file.WriteInt(m_font_index);
    if  (!rc) break;
    rc = file.WriteString(m_font_name);
    if  (!rc) break;
    {
      // 18 October 2002 Dale Lear:
      //   Lowell, wchar_t has different sizes on different OSs.
      //   When writing a wchar_t string, you should use one
      //   of the WriteString functions.  This function must continue
      //   to use WriteShort(64,...) so old files will remain valid.
      unsigned short sh[64];
      memset(sh,0,sizeof(sh));
      int i;
      for ( i = 0; i < 64 && i < face_name_size-1; i++ )
        sh[i] = m_facename[i];
      rc = file.WriteShort(64, sh);
      if  (!rc) break;
    }

    // 1.1 additions
    rc = file.WriteInt( m_font_weight);
    if  (!rc) break;
    rc = file.WriteInt( m_font_italic);
    if  (!rc) break;
    rc = file.WriteDouble( m_linefeed_ratio);
    if  (!rc) break;

    // 1.2 addition
    rc = file.WriteUuid( m_font_id);
    if (!rc) break;

    // 1.3 addition
//    rc = file.WriteInt( m_font_underlined);
//    if  (!rc) break;

    break;
  }

  return rc;
}

ON_BOOL32 ON_Font::Read(
       ON_BinaryArchive& file // restore definition from binary archive
     )
{
  Defaults();
  m_font_index = -1;
  int major_version = 0;
  int minor_version = 0;
  bool rc = file.Read3dmChunkVersion(&major_version,&minor_version);
  if ( rc && major_version == 1 ) 
  {
    int i;
    for(;;)
    {
      rc = file.ReadInt( &m_font_index );
      if  (!rc) break;
      rc = file.ReadString( m_font_name );
      if  (!rc) break;

      {
        // 18 October 2002 Dale Lear:
        //   Lowell, wchar_t has different sizes on different OSs.
        //   When writing a wchar_t string, you should use one
        //   of the WriteString functions.  This function must continue
        //   to use ReadShort(64,...) so old files will remain valid.
        unsigned short sh[64];
        rc = file.ReadShort(64, sh);
        if (!rc) break;

        wchar_t facename[65];
        for ( i = 0; i < 64; i++ )
        {
          facename[i] = sh[i];
        }
        facename[64] = 0;
        SetFontFaceName(facename);
      }

      if( minor_version >= 1 )
      {
        rc = file.ReadInt( &i );
        if (!rc) break;
        SetFontWeight(i);

        rc = file.ReadInt( &i);
        if (!rc) break;
        SetIsItalic(i?true:false);

        rc = file.ReadDouble( &m_linefeed_ratio );
        if (!rc) break;

        if ( minor_version >= 2 )
        {
          rc = file.ReadUuid( m_font_id );
          if (!rc) break;
        }
        //if ( minor_version >= 3 )
        //{
        //  rc = file.ReadInt( &i);
        //  if (!rc) break;
        //  SetUnderlined(i?true:false);
        //}
      }

      break;
    }
  }
  else
  {
    ON_ERROR("ON_Font::Read - get newer version of opennurbs");
    rc = false;
  }

  return rc;
}

// Ratio of linefeed to character height
const double ON_Font::m_default_linefeed_ratio = 1.6;

// This must be an 'I' or 'H', but we have not tested 'H'.
// There are problems with any other upper case character.
// In particular, the standard 'M' does not work.
const int ON_Font::m_metrics_char = 'I';




