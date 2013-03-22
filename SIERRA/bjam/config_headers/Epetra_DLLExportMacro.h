/* #undef BUILD_SHARED_LIBS */
#if defined (_WIN32) && defined (BUILD_SHARED_LIBS)
#  if defined(EPETRA_LIB_EXPORTS_MODE)
#    define EPETRA_LIB_DLL_EXPORT __declspec(dllexport)
#  else
#    define EPETRA_LIB_DLL_EXPORT __declspec(dllimport)
#  endif
#else
#  define EPETRA_LIB_DLL_EXPORT
#endif
