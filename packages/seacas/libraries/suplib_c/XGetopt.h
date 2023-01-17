// XGetopt.h  Version 1.2
//
// Author:  Hans Dietrich
//          hdietrich2@hotmail.com
//
// This software is released into the public domain.
// You are free to use it in any way you like.
//
// This software is provided "as is" with no expressed
// or implied warranty.  I accept no liability for any
// damage or loss of business that this software may cause.
//
///////////////////////////////////////////////////////////////////////////////
#pragma once

#ifdef __cplusplus
extern "C" {
#endif
extern int   optind;
extern int   optopt;
extern int   opterr;
extern char *optarg;

int getopt(int argc, char *const argv[], const char *optstring);
#ifdef __cplusplus
} /* close brackets on extern "C" declaration */
#endif
