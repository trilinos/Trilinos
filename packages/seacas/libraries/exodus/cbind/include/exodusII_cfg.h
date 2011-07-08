#ifndef __exodusII_cfg_h
#define __exodusII_cfg_h

#define EXODUSII_VERSION_MAJOR 5
#define EXODUSII_VERSION_MINOR 05
#define EXODUSII_VERSION_PATCH 0
#define EXODUSII_VERSION "5.05.0"

#if defined(_WIN32) && !defined(WIN32)
# define WIN32
#endif

/* #undef EXODUSII_BUILD_SHARED_LIBS */
#if defined(WIN32) && defined(EXODUSII_BUILD_SHARED_LIBS)
#  if defined(exoIIc_EXPORTS)
#    define EXODUS_EXPORT __declspec( dllexport ) extern
#  else
#    define EXODUS_EXPORT __declspec( dllimport ) extern
#  endif
#else
#  define EXODUS_EXPORT extern
#endif

#endif

