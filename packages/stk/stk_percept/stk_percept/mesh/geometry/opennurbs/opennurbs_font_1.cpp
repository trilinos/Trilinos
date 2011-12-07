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

//srkenno@sandia.gov: split this from opennurbs_font.cpp to workaround a potential Mac/Darwin gcc4.4 compiler/linker bug

//////////////////////////////////////////////////////////////////////
//
// Interface
void ON_Font::SetFontName( const wchar_t* s )
{
  m_font_name = s;
}

void ON_Font::SetFontName( const char* s )
{
  m_font_name = s;
}

void ON_Font::GetFontName( ON_wString& s ) const
{
  s = m_font_name;
}

const wchar_t* ON_Font::FontName() const
{
  const wchar_t* s = m_font_name;
  return s;
}

#if defined(ON_OS_WINDOWS_GDI)
static 
int CALLBACK ON__IsSymbolFontFaceNameHelper( ENUMLOGFONTEX*, NEWTEXTMETRICEX*, DWORD, LPARAM)
{
  // If the fontname in the logfont structure has
  // a corresponding symbol font on the system, 
  // set the  lfCharSet member to SYMBOL_CHARSET, 
  // otherwise DEFAULT_CHARSET
  // The input logfont structure may be modified.
  return 7;
}
#endif

bool ON_Font::IsSymbolFontFaceName( const wchar_t* s)
{
  bool rc = false;

#if defined(ON_OS_WINDOWS_GDI)
  if( s && s[0])
  {
    HDC hdc = ::GetDC( NULL);    
    if( hdc)
    {      
      LOGFONT logfont;
      memset( &logfont, 0, sizeof( logfont));
      int i;
      for ( i = 0; i < LF_FACESIZE && s[i]; i++ )
      {
        logfont.lfFaceName[i] = s[i];
      }
      logfont.lfCharSet = ON_Font::symbol_charset;
      if( 7 == ::EnumFontFamiliesEx( hdc, &logfont, (FONTENUMPROC)ON__IsSymbolFontFaceNameHelper, 0, 0))
      {
        rc = true;
      }    
      ::ReleaseDC( NULL, hdc);
    }
  }
#endif

  return rc;
}

bool ON_Font::SetFontFaceName( const wchar_t* s )
{
  int i;
  memset( &m_facename, 0, sizeof(m_facename) );
  if ( s)
  {
    for ( i = 0; i < face_name_size-1 && s[i]; i++ )
    {
      m_facename[i] = s[i];
    }
  }

#if defined(ON_OS_WINDOWS_GDI)
  memset( &m_logfont.lfFaceName, 0, sizeof(m_logfont.lfFaceName) );
  if ( s )
  {
    for ( i = 0; i < LF_FACESIZE && s[i]; i++ )
    {
      m_logfont.lfFaceName[i] = s[i];
    }
  }

  m_logfont.lfCharSet = ON_Font::IsSymbolFontFaceName( s)
                      ? ((unsigned char)ON_Font::symbol_charset)
                      : ((unsigned char)ON_Font::default_charset);
#endif

  m_I_height = 0;

  return( m_facename[0] ? true : false);
}

bool ON_Font::SetFontFaceName( const char* s )
{
  ON_wString wstr(s);
  const wchar_t* w = wstr;
  return SetFontFaceName(w);
}

double ON_Font::AscentRatio() const
{
  return ((double)normal_font_height) / ((double)HeightOfI());
}

void ON_Font::GetFontFaceName( ON_wString& s ) const
{
  s = m_facename;
}

const wchar_t* ON_Font::FontFaceName() const
{
  const wchar_t* s = m_facename;
  return s;
}

void ON_Font::SetFontIndex(int i)
{
  m_font_index = i;
}

int ON_Font::FontIndex() const
{
  return m_font_index;
}

double ON_Font::LinefeedRatio() const
{
  return m_linefeed_ratio;
}

void ON_Font::SetLinefeedRatio( double d)
{
  m_linefeed_ratio = d;
}

int ON_Font::FontWeight() const
{
  return m_font_weight;
}

void ON_Font::SetFontWeight( int w)
{
  if ( w != m_font_weight )
  {
    if ( w < 0 )
    {
      w = 0;
    }

    m_font_weight = w;

#if defined(ON_OS_WINDOWS_GDI)
    m_logfont.lfWeight = m_font_weight;
#endif

    m_I_height = 0;
  }
}

bool ON_Font::IsItalic() const
{
  return m_font_italic;
}

void ON_Font::SetIsItalic( bool b)
{
  SetItalic( b);
}


void ON_Font::SetItalic( bool b)
{
  if ( m_font_italic != b )
  {
    m_font_italic = b?true:false;
#if defined(ON_OS_WINDOWS_GDI)
    m_logfont.lfItalic = m_font_italic;
#endif
    m_I_height = 0;
  }
}


bool ON_Font::IsBold() const
{
  return (FontWeight() >= bold_weight);
}


void ON_Font::SetBold( bool bBold )
{
  SetFontWeight( bBold ? bold_weight : normal_weight);
}


bool ON_Font::IsUnderlined() const
{
  return m_font_underlined;
}


void ON_Font::SetUnderlined( bool b)
{
  if ( m_font_underlined != b )
  {
    m_font_underlined = b?true:false;
#if defined(ON_OS_WINDOWS_GDI)
    m_logfont.lfUnderline = m_font_underlined;
#endif
  }
}



/*
Returns:
  Height of the 'I' character when the font is drawn 
  with m_logfont.lfHeight = 256.
*/
int ON_Font::HeightOfI() const
{
  if ( m_I_height  <= 0 )
  {
    // Default is height of Arial 'I'.  If we are running
    // on Windows, then we calculate the actual height of
    // an 'I' in the font.
    //   The ..ON_Font::normal_font_height/256 is here 
    //   so this code will continue to work correctly 
    //   if somebody changes ON_Font::normal_font_height.
    int I_height = (166*ON_Font::normal_font_height)/256;

#if defined(ON_OS_WINDOWS_GDI)
    if ( m_logfont.lfFaceName[0] )
    {
      // Get the height of an 'I'
      HDC hdc = ::GetDC( NULL);
      if (hdc)
      {
        LOGFONT logfont = m_logfont;
        logfont.lfHeight = normal_font_height;
        HFONT font = ::CreateFontIndirect( &logfont);
        if ( font )
        {
          wchar_t str[2];
          str[0] = ON_Font::m_metrics_char;
          str[1] = 0;
          HFONT oldfont = (HFONT)::SelectObject( hdc, font);
          ::SetBkMode( hdc, TRANSPARENT);
          ::BeginPath(hdc);
          ::ExtTextOut( hdc, 0, 0, 0, NULL, (LPCSTR)str, 1, NULL);
          ::EndPath( hdc);
          int numPoints = ::GetPath( hdc, NULL, NULL, 0);

          if( numPoints > 2)
          {
            // Allocate room for the points & point types
            LPPOINT pPoints = (LPPOINT)onmalloc( numPoints * sizeof(*pPoints) );
            LPBYTE pTypes = (LPBYTE)onmalloc( numPoints * sizeof(*pTypes) );
            if ( pTypes && pPoints)
            {
              // Get the points and types from the current path
              numPoints = ::GetPath( hdc, pPoints, pTypes, numPoints);
              if( numPoints > 2)
              {
                int ymin = pPoints[0].y;
                int ymax = ymin;
                int k;
                for( k = 1; k < numPoints; k++)
                {
                  if( pPoints[k].y < ymin)
                    ymin = pPoints[k].y;
                  else if( pPoints[k].y > ymax)
                    ymax = pPoints[k].y;
                }
                I_height = ymax - ymin + 1;
              }
            }
            onfree( pPoints);
            onfree( pTypes);
          }
          ::SelectObject( hdc, oldfont);
          ::DeleteObject( font);
        }
      }
      ::ReleaseDC( NULL, hdc);
    }
#endif
    const_cast<ON_Font*>(this)->m_I_height = I_height;
  }
  return m_I_height;
}


int ON_Font::HeightOfLinefeed() const
{
  return ( (int)( ceil(m_linefeed_ratio*HeightOfI()) ) );
}

#if defined(ON_OS_WINDOWS_GDI)

#pragma message( " --- OpenNURBS including Windows LOGFONT support in ON_Font" )

bool ON_Font::SetLogFont( const LOGFONT& logfont )
{
  if ( &m_logfont != &logfont )
  {
    memcpy(&m_logfont,&logfont,sizeof(m_logfont));
  }
    
  // synch persistent fields
  m_font_weight = m_logfont.lfWeight;
  m_font_italic = (m_logfont.lfItalic?true:false);
  memset(&m_facename[0],0,sizeof(m_facename));
  int i;
  for ( i = 0; i < face_name_size && i < LF_FACESIZE; i++ )
  {
    m_facename[i] = (wchar_t)m_logfont.lfFaceName[i]; 
  }
  m_facename[face_name_size-1] = 0;

  m_I_height = 0;

  return true;
}

const LOGFONT& ON_Font::LogFont() const
{
  return m_logfont;
}

ON_Font::ON_Font( const LOGFONT& logfont )
{
  Defaults();
  SetLogFont(logfont);
}

ON_Font& ON_Font::operator=( const LOGFONT& logfont )
{
  if ( &m_logfont == &logfont )
  {
    LOGFONT lf = logfont;
    SetLogFont(lf);
  }
  else
  {
    SetLogFont(logfont);
  }
  return *this;
}

bool ON_Font::CompareFontCharacteristics( ON_Font& other_font, bool bCompareName) const
{
  if( bCompareName && m_font_name.CompareNoCase( other_font.m_font_name))
    return false;

  if( m_font_weight != other_font.m_font_weight)
    return false;
  
  if( m_font_italic != other_font.m_font_italic)
    return false;
  
  if( m_font_underlined != other_font.m_font_underlined)
    return false;
  
  if( m_linefeed_ratio != other_font.m_linefeed_ratio)
    return false;
  
  if( _wcsicmp( m_facename, other_font.m_facename))
    return false;

  return true;
}


#endif
