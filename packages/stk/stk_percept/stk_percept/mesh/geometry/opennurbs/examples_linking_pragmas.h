#if defined(_MSC_VER)

// This file is specific to Micrsoft's compiler.
// It contains linking pragmas for building the opennurbs examples.

#pragma once

#if defined(ON_DLL_EXPORTS)
// If you get the following error, your compiler settings
// indicate you are building opennurbs as a DLL. This file
// is used for linking with opennurbs.
#error This file contains linking pragmas for using opennurbs.
#endif


#if defined(WIN64)

// x64 (64 bit) static libraries

#if defined(NDEBUG)

// Release x64 (64 bit) libs
#pragma message( " --- Opennurbs examples Release x64 (64 bit) build." )
#if defined(ON_DLL_IMPORTS)
#pragma comment(lib, "../x64/Release/opennurbs.lib")
#else
#pragma comment(lib, "../zlib/x64/Release/zlib.lib")
#pragma comment(lib, "../x64/Release/opennurbs_staticlib.lib")
#endif

#else // _DEBUG

// Debug x64 (64 bit) libs
#pragma message( " --- Opennurbs examples Debug x64 (64 bit) build." )
#if defined(ON_DLL_IMPORTS)
#pragma comment(lib, "../x64/Debug/opennurbs.lib")
#else
#pragma comment(lib, "../zlib/x64/Debug/zlib.lib")
#pragma comment(lib, "../x64/Debug/opennurbs_staticlib.lib")
#endif

#endif // NDEBUG else _DEBUG

#else // WIN32

// x86 (32 bit) static libraries

#if defined(NDEBUG)

// Release x86 (32 bit) libs
#pragma message( " --- Opennurbs examples Release x86 (32 bit) build." )
#if defined(ON_DLL_IMPORTS)
#pragma comment(lib, "../Release/opennurbs.lib")
#else
#pragma comment(lib, "../zlib/Release/zlib.lib")
#pragma comment(lib, "../Release/opennurbs_staticlib.lib")
#endif

#else // _DEBUG

// Debug x86 (32 bit) libs
#pragma message( " --- Opennurbs examples Debug x86 (32 bit) build." )
#if defined(ON_DLL_IMPORTS)
#pragma comment(lib, "../Debug/opennurbs.lib")
#else
#pragma comment(lib, "../zlib/Debug/zlib.lib")
#pragma comment(lib, "../Debug/opennurbs_staticlib.lib")
#endif

#endif // NDEBUG else _DEBUG

#endif // WIN64 else WIN32

#endif
