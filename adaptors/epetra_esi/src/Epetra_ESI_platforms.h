#ifndef _Epetra_ESI_platforms_h_
#define _Epetra_ESI_platforms_h_

//
//This header does any platform-specific checks and #defines as necessary.
//

#if defined(__PGI) && defined (__i386)
//Compiling on ASCI Red (janus, a.k.a. TFLOP)
#define ESI_NO_COMPLEX
#define TFLOP
#endif

#if (defined(__SUNPRO_CC) && __SUNPRO_CC < 0x500)
#define ESI_SIMULATE_BOOL
#endif

#endif

