/* #undef BUILD_SHARED_LIBS */
#if defined (_WIN32) && defined (BUILD_SHARED_LIBS)
#  if defined(ANASAZI_LIB_EXPORTS_MODE)
#    define ANASAZI_LIB_DLL_EXPORT __declspec(dllexport)
#  else
#    define ANASAZI_LIB_DLL_EXPORT __declspec(dllimport)
#  endif
#else
#  define ANASAZI_LIB_DLL_EXPORT
#endif
