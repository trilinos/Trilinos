#if !defined(OPENNURBS_EXAMPLE_UD_INC_)
#define OPENNURBS_EXAMPLE_UD_INC_

class CExampleWriteUserData : public ON_UserData
{
  static int m__sn;
  ON_OBJECT_DECLARE(CExampleWriteUserData);

public:
  static ON_UUID Id();

  CExampleWriteUserData();
  CExampleWriteUserData( const char* s);
  CExampleWriteUserData(const CExampleWriteUserData& src);
  CExampleWriteUserData& operator=(const CExampleWriteUserData& src);
  ~CExampleWriteUserData();
  void Dump( ON_TextLog& text_log ) const;
  ON_BOOL32 GetDescription( ON_wString& description );
  ON_BOOL32 Archive() const;
  ON_BOOL32 Write(ON_BinaryArchive& file) const;
  ON_BOOL32 Read(ON_BinaryArchive& file);

  ON_wString m_str;
  int m_sn;
};

#endif
