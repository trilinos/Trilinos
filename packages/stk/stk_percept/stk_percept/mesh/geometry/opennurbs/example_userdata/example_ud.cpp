// uncomment the "ON_DLL_IMPORTS" define to use opennurbs as a Windows DLL
//#define ON_DLL_IMPORTS
#include "../opennurbs.h"
#include "example_ud.h"

CExampleWriteUserData::CExampleWriteUserData()
{
  m_sn = ++m__sn;
  m_userdata_uuid = Id();
  m_userdata_copycount = 1;
}

// {2532FB4A-DED9-4600-B6A4-1568504B66A5}
static const ON_UUID ExampleWriteUserData_Id = 
{ 0x2532fb4a, 0xded9, 0x4600, { 0xb6, 0xa4, 0x15, 0x68, 0x50, 0x4b, 0x66, 0xa5 } };

CExampleWriteUserData::CExampleWriteUserData( const char* s)
{
  m_sn = ++m__sn;
  m_userdata_uuid = Id();
  m_application_uuid = ExampleWriteUserData_Id;
  m_userdata_copycount = 1;
  m_str = s;
}

CExampleWriteUserData::CExampleWriteUserData(const CExampleWriteUserData& src) : ON_UserData(src), m_str(src.m_str)
{
  m_sn = ++m__sn;
}

CExampleWriteUserData& CExampleWriteUserData::operator=(const CExampleWriteUserData& src)
{
  if ( this != &src )
  {
    ON_UserData::operator=(src);
    m_str = src.m_str;
  }
  return *this;
}

CExampleWriteUserData::~CExampleWriteUserData()
{
  m_sn = -abs(m_sn);
}

void CExampleWriteUserData::Dump( ON_TextLog& text_log ) const
{
  ON_UserData::Dump(text_log);
  text_log.PushIndent();
  const wchar_t* s = m_str;
  if ( 0 == s ) 
    s = L"";
  text_log.Print("m_str: %S\n",s);
  text_log.Print("m_sn: %d\n",m_sn);
  text_log.PopIndent();
}

ON_BOOL32 CExampleWriteUserData::GetDescription( ON_wString& description )
{
  description = L"example_write.exe user data";
  return true;
}


ON_BOOL32 CExampleWriteUserData::Archive() const
{
  return true;
}

ON_BOOL32 CExampleWriteUserData::Write(ON_BinaryArchive& file) const
{
  return file.WriteString(m_str);
}

ON_BOOL32 CExampleWriteUserData::Read(ON_BinaryArchive& file)
{
  return file.ReadString(m_str);
}

int CExampleWriteUserData::m__sn = 0;

ON_OBJECT_IMPLEMENT(CExampleWriteUserData,ON_UserData,"DADD17C5-706D-44ea-9B13-7D9D2C56D085");

ON_UUID CExampleWriteUserData::Id()
{
  // {6FC7CDF1-751E-4fa0-9D86-73E84D416DD7}
  static const ON_UUID id = 
  { 0x6fc7cdf1, 0x751e, 0x4fa0, { 0x9d, 0x86, 0x73, 0xe8, 0x4d, 0x41, 0x6d, 0xd7 } };
  return id;
}

