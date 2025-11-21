
#if defined (_WIN32) && defined (BUILD_SHARED_LIBS)
#  if defined(TEUCHOSCORE_LIB_EXPORTS_MODE)
#    define TEUCHOSCORE_LIB_DLL_EXPORT __declspec(dllexport)
#  else
#    define TEUCHOSCORE_LIB_DLL_EXPORT __declspec(dllimport)
#  endif
#else
#  define TEUCHOSCORE_LIB_DLL_EXPORT
#endif

#if defined (_WIN32) && defined (BUILD_SHARED_LIBS)
#  if defined(TEUCHOSCOMM_LIB_EXPORTS_MODE)
#    define TEUCHOSCOMM_LIB_DLL_EXPORT __declspec(dllexport)
#  else
#    define TEUCHOSCOMM_LIB_DLL_EXPORT __declspec(dllimport)
#  endif
#else
#  define TEUCHOSCOMM_LIB_DLL_EXPORT
#endif

#if defined (_WIN32) && defined (BUILD_SHARED_LIBS)
#  if defined(TEUCHOSPARAMETERLIST_LIB_EXPORTS_MODE)
#    define TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT __declspec(dllexport)
#  else
#    define TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT __declspec(dllimport)
#  endif
#else
#  define TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
#endif

#if defined (_WIN32) && defined (BUILD_SHARED_LIBS)
#  if defined(TEUCHOSNUMERICS_LIB_EXPORTS_MODE)
#    define TEUCHOSNUMERICS_LIB_DLL_EXPORT __declspec(dllexport)
#  else
#    define TEUCHOSNUMERICS_LIB_DLL_EXPORT __declspec(dllimport)
#  endif
#else
#  define TEUCHOSNUMERICS_LIB_DLL_EXPORT
#endif

/* Defining this to allow the backwards compatibility testing to work. This is not needed for library code */
#define TEUCHOS_LIB_DLL_EXPORT

/* There is not export stuff used in the remainder subpackage yet. */
